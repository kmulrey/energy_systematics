import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser
import cPickle
import re
from scipy.signal import hilbert
from scipy.signal import resample
import scipy.fftpack as fftp
import os
import pycrtools as cr
import test_process_func as prf # NB. Testing / rewriting process_func

def getJonesMatrix(azimuth, zenith, frequencies):
    # Produce Jones matrix from antenna model, for position (azimuth, zenith) on the sky
    # as a function of 'frequencies' (freq. axis)
    
    #read tables antenna model and initialize values
    vt=np.loadtxt(os.environ["LOFARSOFT"] + "/data/lofar/antenna_response_model/LBA_Vout_theta.txt", skiprows=1)
    vp=np.loadtxt(os.environ["LOFARSOFT"] + "/data/lofar/antenna_response_model/LBA_Vout_phi.txt", skiprows=1)
    cvt = cr.hArray(vt[:, 3] + 1j * vt[:, 4])
    cvp = cr.hArray(vp[:, 3] + 1j * vp[:, 4])
    fstart = 10.0 * 1.e6 # Default numbers used for computing Jones matrix from antenna model
    fstep = 1.0 * 1.e6
    fn = 101
    ttstart = 0.0
    ttstep = 5.0
    ttn = 19
    pstart = 0.0
    pstep = 10.0
    pn = 37

    jones_matrix = cr.hArray(complex, dimensions=(len(frequencies), 2, 2))
    for k,f in enumerate(frequencies):
        if (f>1e7 and f<1e8):
            cr.hGetJonesMatrix(jones_matrix[k], f, 180-(azimuth/np.pi*180),90-(zenith/np.pi*180), cvt, cvp, fstart, fstep, fn, ttstart, ttstep, ttn, pstart, pstep, pn)
    jm=jones_matrix.toNumpy()

    return jm

def measurePulseEnergy(timeseries, half_width, peak_time):
    # Input: timeseries, as 2D numpy-array of shape (npol, nsamples) where npol = # polarizations
    # Pulse power is integrated from the peak position, including 'half_width' samples to each side
    # peak_time = the time of the maximum in the strongest polarization
    
    # rotate the time series to have the peak position at sample 'half_width'
    shifted_timeseries = np.roll(timeseries, int(half_width-peak_time), axis=1)
    # so we can integrate from sample 0 up to (incl) 2*half_width, i.e. 2*half_width+1 samples
    power = np.sum(np.square(shifted_timeseries[:, 0:(2*half_width+1)]), axis=-1) * 5.0e-9
    # Power is from voltage, i.e. after applying antenna model. Therefore, no Z0 factor

    return power

def ProcessData(datadir,fileno, lorafile_suffix='', debug_testpulse=False, debug_lofarpulse=None):
    # set lorafile_suffix if a different run of LORAsimulation is desired
    lowco=30
    hico=80
    nantennas=160
    B_inc = 1.1837
    
    debug_figure = None
    onskypower=np.zeros([nantennas,2])
    antenna_position=np.zeros([nantennas,3])
    filteredpower=np.zeros([nantennas,2])
    power=np.zeros([nantennas,2])
    power11=np.zeros([nantennas,2])
    power21=np.zeros([nantennas,2])
    power41=np.zeros([nantennas,2])
    peak_time=np.zeros([nantennas,2])
    peak_bin=np.zeros([nantennas,2])
    peak_amplitude=np.zeros([nantennas,2])
    pol_angle=np.zeros([nantennas])
    pol_angle_filt=np.zeros([nantennas])

    longfile = '{0}/DAT{1}.long'.format(datadir,str(fileno).zfill(6))
    steerfile = '{0}/steering/RUN{1}.inp'.format(datadir,str(fileno).zfill(6))
    listfile = open('{0}/steering/SIM{1}.list'.format(datadir,str(fileno).zfill(6)))
    lorafile = '{0}/DAT{1}{2}.lora'.format(datadir,str(fileno).zfill(6), lorafile_suffix)
    
    longdata=np.genfromtxt(longfile, skip_header=2, skip_footer=5, usecols=(0,2,3))
    xlength=np.argmax(np.isnan(longdata[:,0])) # Array length of longitudinal profile data
    Xground=xlength*1.0 # Profile is in steps of 10 g/cm2 (can be set in steering file), so Xground = xlength * 10
    # NB. Setting to 1.0 g/cm2 now (Sept 2017).
    profile = longdata[0:xlength,:]
    longprofile=np.zeros([4000,3]) #The profiles have different lengths, unlikely to exceed 400... (maybe inclined events??) NB. When running with 1 g/cm2 resolution instead of 10, make it 4000.
    longprofile[0:xlength,:]=np.array([profile[:,0],profile[:,1]+profile[:,2],profile[:,2]-profile[:,1]]).T
    # Tables in .long file: X, # electrons, # positrons
    # So: longprofile contains (X, # electrons+positrons, net charge) for each step.
    
    hillas = np.genfromtxt(re.findall("PARAMETERS.*",open(longfile,'r').read()))[2:]
    zenith=(np.genfromtxt(re.findall("THETAP.*",open(steerfile,'r').read()))[1])*np.pi/180. #rad; CORSIKA coordinates
    azimuth=np.mod(np.genfromtxt(re.findall("PHIP.*",open(steerfile,'r').read()))[1],360)*np.pi/180.  #rad; CORSIKA coordinates
    energy=np.genfromtxt(re.findall("ERANGE.*",open(steerfile,'r').read()))[1] #GeV

    az_rot=3*np.pi/2+azimuth    #conversion from CORSIKA coordinates to standard axis coordinates: 0=east, pi/2=north (Corsika has 0=north, pi/2=west)
    zen_rot=zenith

    lines = listfile.readlines()
    debug_XYZ_power = [] # cannot avoid returning this parameter when not in debug mode...?
    for j in np.arange(nantennas):
        antenna_position[j] = (lines[j].split(" ")[2:5]) #read antenna position...
        antenna_file = lines[j].split(" ")[5]   #... and output filename from the antenna list file
        coreasfile = '{0}/SIM{1}_coreas/raw_{2}.dat'.format(datadir,str(fileno).zfill(6),antenna_file[:-1]) #drop the \n from the string!
        data=np.genfromtxt(coreasfile)
        data[:,1:]*=2.99792458e4 # convert Ex, Ey and Ez (not time!) to Volt/meter
        dlength=data.shape[0]
        poldata=np.ndarray([dlength,2])
        XYZ=np.zeros([dlength,3])
        XYZ[:,0]=-data[:,2] #conversion from CORSIKA coordinates to 0=east, pi/2=north
        XYZ[:,1]=data[:,1]
        XYZ[:,2]=data[:,3]
        
        # In debug mode (debug_testpulse==True), replace XYZ data by either a test pulse (0-100 MHz)
        # or by LOFAR calibrated XYZ data for the given event (to be passed along when calling ProcessData)
        
        if debug_testpulse and (debug_lofarpulse is not None):
            # LOFAR pulse is time series in X,Y,Z for one antenna (currently). Upsample to 0.1 ns sampling rate
            # and cut out pulse to a length corresponding to dlength
            #if j == 0:
            # hack direction to (pseudo-)zenith
            broad_window = 1024
            zen_rot = 45.0 * np.pi/180
            zenith = 45.0 * np.pi/180
            azimuth = (180.0 - 0.0) * np.pi/180
            az_rot = 3*np.pi/2 + azimuth
            # Get pulse from LOFAR with index j modulo nof LOFAR antennas...
            pulse_index = j % (len(debug_lofarpulse) / 3)
            thispulse = debug_lofarpulse[pulse_index*3:(pulse_index+1)*3]
            start = 32768 - broad_window
            end = 32768 + broad_window
            thispulse = thispulse[:, start:end]
            lofarpulse_upsampled = resample(thispulse, len(thispulse[0])*50, axis=1)
            maxpos = np.median(np.argmax(lofarpulse_upsampled, axis=1)) # median over X,Y,Z
            start = int(maxpos - dlength/2)
            lofarpulse_cutout = lofarpulse_upsampled[:, start:start+dlength]
        
            XYZ = np.copy(lofarpulse_cutout).T # this is the test data to be injected
            
            # Get power in XYZ polarizations, for comparison
            XYZ_filtered = resample(XYZ, len(lofarpulse_cutout[0]) / 50, axis=-2)

            XYZ_power_timeseries = XYZ_filtered**2
            sum_XYZ_power_timeseries = np.sum(XYZ_power_timeseries, axis=0)
            # Get strongest polarization, then use that to get position (sample) of maximum power
            strongest_pol_XYZ = np.argmax(sum_XYZ_power_timeseries)

            maxpos = np.argmax(XYZ_power_timeseries[:, strongest_pol_XYZ]) # The time of maximum power in the strongest polarization, required for measurePulseEnergy
            #print maxpos
            Z0 = 120.0 * np.pi # almost
            XYZ_power = measurePulseEnergy(XYZ_filtered.T, 5, maxpos) # transpose to get shape (npol, nsamples)
            XYZ_power /= Z0
            debug_XYZ_power.append(XYZ_power)
            if (j==0): print 'XYZ power = (%1.5e, %1.5e, %1.5e), total = %1.5e' % (XYZ_power[0], XYZ_power[1], XYZ_power[2], np.sum(XYZ_power))
        
            data[:, 2] = - XYZ[:, 0] # overwrite Coreas data with test pulse
            data[:, 1] = XYZ[:, 1]
            data[:, 3] = XYZ[:, 2]
        
        elif debug_testpulse:
            testpulse = np.zeros(1 + dlength / 50)
            testpulse[dlength/2] += 100.0
            testpulse_upsampled = resample(testpulse, len(testpulse)*50)
            testpulse_upsampled = testpulse_upsampled[0:dlength]
            XYZ[:, 0] = np.copy(testpulse_upsampled)
            XYZ[:, 1] = np.copy(testpulse_upsampled)
            XYZ[:, 2] = np.copy(testpulse_upsampled)
        
        
        
        # Convert to, v, vxB, vxvxB coordinates to compute Stokes parameters and polarization angle
        UVW = prf.GetUVW(XYZ, 0, 0, 0, zen_rot, az_rot, B_inc)
        Stokes=prf.stokes_parameters(UVW[:,0],UVW[:,1],fftp.hilbert(UVW[:,0]),fftp.hilbert(UVW[:,1]))
        pol_angle[j]=prf.polarization_angle(Stokes)
        UVWfilt=prf.FreqFilter(UVW, lowco, hico, data[1,0]-data[0,0])
        Stokesfilt=prf.stokes_parameters(UVWfilt[:,0],UVWfilt[:,1],fftp.hilbert(UVWfilt[:,0]),fftp.hilbert(UVWfilt[:,1]))
        pol_angle_filt[j]=prf.polarization_angle(Stokesfilt)
        # Convert to on-sky coordinates (n, theta, phi) to prepare for application of antenna model
        # poldata[:,0] = -1.0/np.sin(zen_rot)*data[:,3] # -1/sin(theta) *z
        # NB. Not guaranteed to hold, only if v-component is really zero... otherwise use the full formula here:
        poldata[:, 0] = -data[:, 2] * np.cos(zen_rot)*np.cos(az_rot) + data[:, 1] * np.cos(zen_rot)*np.sin(az_rot) - np.sin(zen_rot)*data[:, 3]
        poldata[:,1] = np.sin(az_rot)*data[:,2] + np.cos(az_rot)*data[:,1] # -sin(phi) *x + cos(phi)*y in coREAS 0=positive y, 1=negative x
        spec=np.fft.rfft(poldata, axis=-2)
        # Apply antenna model
        tstep = data[1,0]-data[0,0]
        onskypower[j]=np.array([np.sum(poldata[:,0]*poldata[:,0]),np.sum(poldata[:,1]*poldata[:,1])])*tstep
        # On sky power: this is E-field, so power = tstep * sum |E|^2 / Z0
        
        #freqhi = 0.5/tstep/1e6 # MHz
        #freqstep = freqhi/(dlength/2+1) # MHz
        #frequencies = np.arange(0,freqhi,freqstep)*1e6 # Hz WRONG
        #frequencies = np.arange(0,dlength/2+1)*freqstep*1e6 WRONG
        frequencies = np.fft.rfftfreq(dlength, tstep) # Ends at 5000 MHz as it should for tstep=0.1 ns
        freqstep = (frequencies[1] - frequencies[0]) / 1.0e6 # MHz
        
        
        #katie--> changing zenith
        if j == 0: # Jones matrix will be the same for every antenna
            jm = getJonesMatrix(azimuth, zenith, frequencies)
            #print '\n\n  getting new theta value!!! {0} -> {1}'.format(zenith*180/np.pi,(zenith-(5*np.pi/180.0))*180/np.pi)
        
        instr_spec=np.ndarray([dlength/2+1,2],dtype=complex)
        instr_spec[:,0] = jm[:,0,0] * spec[:,0] + jm[:,0,1] * spec[:,1]
        instr_spec[:,1] = jm[:,1,0] * spec[:,0] + jm[:,1,1] * spec[:,1]
        #Apply window and reduce maximum frequency to acquire downsampled signal
        fb = int(np.floor(lowco/freqstep))
        lb = int(np.floor(hico/freqstep)+1)
        window = np.zeros([1,dlength/2+1,1])
        window[0,fb:lb+1,0]=1
        pow0=np.abs(instr_spec[:,0])*np.abs(instr_spec[:,0])
        pow1=np.abs(instr_spec[:,1])*np.abs(instr_spec[:,1])
        ospow0=np.abs(spec[:,0])*np.abs(spec[:,0])
        ospow1=np.abs(spec[:,1])*np.abs(spec[:,1])
        power[j]=np.array([np.sum(pow0[fb:lb+1]),np.sum(pow1[fb:lb+1])])/(dlength/2.)*tstep
        filteredpower[j]=np.array([np.sum(ospow0[fb:lb+1]),np.sum(ospow1[fb:lb+1])])/(dlength/2.)*tstep
        # assume that simulated time resolution is higher than LOFAR time resolution (t_step=5 ns)

        maxfreqbin= int(np.floor(tstep/5e-9 * dlength/2.)+1) # Apply frequency bandpass only up to 100 MHz i.e. LOFAR maximum
        shortspec=np.array([instr_spec[0:maxfreqbin,0]*window[0,0:maxfreqbin,0],instr_spec[0:maxfreqbin,1]*window[0,0:maxfreqbin,0]])
        filt=np.fft.irfft(shortspec, axis=-1)
        # after downsampling, renormalize the signal!
        # nodig??????
        # test power[j] vs power11 (of volledige som)
        dlength_new=filt.shape[1]
        filt *= 1.0*dlength_new/dlength
        # to calculate the time of arrival upsample with a factor 5
        # zelfde factor als in pipeline!!!!
        
        filt_upsampled=resample(filt,16*dlength_new,axis=-1)
        # compute hilbert enevelope
        hilbenv=np.abs(hilbert(filt,axis=-1))
        hilbenv_upsampled=np.abs(hilbert(filt_upsampled,axis=-1))
        # peak_time is the bin where the maximum is located; NOT the actual time of the peak!
        peak_bin[j]=np.argmax(hilbenv,axis=-1)
        peak_time[j]=np.argmax(hilbenv_upsampled,axis=-1)*1e-9 #in seconds
        peak_amplitude[j]=np.max(hilbenv_upsampled,axis=-1)
        if (peak_amplitude[j,0]>peak_amplitude[j,1]):
            pt=peak_bin[j,0]
        else:
            pt=peak_bin[j,1]
        # for 3 different window size, the total power is calculated. The window is allowed to `wrap around', so some voodoo is needed to determine the range:
        
        d=filt.shape[1]
        rng=5
        a=int(np.round(np.max([0,pt-rng])))
        b=int(np.round(pt+rng+1))
        c=int(np.round(np.min([d,pt+d-rng])))
        power11[j]=(np.sum(np.square(filt[:,a:b]),axis=-1)+np.sum(np.square(filt[:,c:d]),axis=-1))*5e-9
        
        test_power11 = measurePulseEnergy(filt, rng, pt) # This is the one to use (AC)
        power11[j] = test_power11 # so put it into power11 array
        if j == 0:
            print 'Power 0/1 polarization: '
            print power11[0]
        #        print 'Antenna %d' % j
        #if j == 98:
        #    import pdb; pdb.set_trace()
        
        #assert abs(test_power11[0] - power11[j][0]) / power11[j][0] < 1.0e-10
        #assert abs(test_power11[1] - power11[j][1]) / power11[j][1] < 1.0e-10

        rng=10
        a=int(np.round(np.max([0,pt-rng])))
        b=int(np.round(pt+rng+1))
        c=int(np.round(np.min([d,pt+d-rng])))
        power21[j]=(np.sum(np.square(filt[:,a:b]),axis=-1)+np.sum(np.square(filt[:,c:d]),axis=-1))*5e-9

        test_power21 = measurePulseEnergy(filt, rng, pt)
        power21[j] = test_power21
        #assert abs(test_power21[0] - power21[j][0]) / power21[j][0] < 1.0e-10
        #assert abs(test_power21[1] - power21[j][1]) / power21[j][1] < 1.0e-10

        rng=20
        a=int(np.round(np.max([0,pt-rng])))
        b=int(np.round(pt+rng+1)) # No guarantee that this will not go over 80 (= length of time trace) !
        c=int(np.round(np.min([d,pt+d-rng]))) # See at antenna 150 of event 164228080
        # where a = 43, b=84, c=80, d=80 hence the c:d sum is empty, and the length is 37 instead of 41...
        power41[j]=(np.sum(np.square(filt[:,a:b]),axis=-1)+np.sum(np.square(filt[:,c:d]),axis=-1))*5e-9
        #        if j == 150:
        #    import pdb; pdb.set_trace()
        test_power41 = measurePulseEnergy(filt, rng, pt)
        power41[j] = test_power41
        #        print 'ant %d: absdiff 0 = %e, absdiff 1 = %e, reldiff 0 = %e, reldiff 1 = %e' % (j, test_power41[0] - power41[j][0], test_power41[1] - power41[j][1], abs(test_power41[0] - power41[j][0]) / power41[j][0], abs(test_power41[1] - power41[j][1]) / power41[j][1])
        #assert abs(test_power41[0] - power41[j][0]) / power41[j][0] < 1.0e-10
        #assert abs(test_power41[1] - power41[j][1]) / power41[j][1] < 1.0e-10

        #        shifted_timeseries = np.roll(filt, int(rng-pt), axis=-1)
        #        test_power41 = np.sum(np.square(shifted_timeseries[:, 0:(2*rng+1)]), axis=-1) * 5.0e-9
        
        #        assert max(abs(test_power41 - power41[j])) < 1.0e-30
        #if c < 80:
        #import pdb; pdb.set_trace()

        # In debug-pulse mode, save time series for one antenna for inspection, plus max/min levels and power in 11 samples
        if debug_lofarpulse is not None and (j==0):
            # make figure from 'filt' and print power levels and max/min levels
            plt.figure()
            plt.plot(filt[0], lw=2, label='pol 0')
            plt.plot(filt[1], lw=2, label='pol 1')
            maxlevels = np.max(filt, axis=1)
            minlevels = np.min(filt, axis=1)
            plt.text(0.1, 0.12, 'Power-11 = (%1.3e, %1.3e)' % (power11[0][0], power11[0][1]), transform=plt.gca().transAxes)
            plt.text(0.1, 0.07, 'Max      = (%1.3e, %1.3e)' % (maxlevels[0], maxlevels[1]), transform=plt.gca().transAxes)
            plt.text(0.1, 0.02, 'Min      = (%1.3e, %1.3e)' % (minlevels[0], minlevels[1]), transform=plt.gca().transAxes)
            plt.title('Time series for LOFAR pulse (XYZ-calibrated) converted to 0/1 pol\nTo be compared to LOFAR (test) pulse in cr_physics pipeline')
            #plt.savefig('/vol/astro7/lofar/sim/pipeline/test_analysis/DEBUG_%d_filtered_01_timeseries.pdf' % )
            debug_figure = plt.gcf()

    # read particle data and energy deposit from Geant4 simulation of LORA
    pdata = np.genfromtxt(lorafile)
    particle_radius=5*pdata[:,0]+2.5
    energy_deposit=pdata[:,1]
   
    ### convert CORSIKA to AUGER coordinates (AUGER y = CORSIKA x, AUGER x = - CORSIKA y
    temp=np.copy(antenna_position)
    antenna_position[:,0], antenna_position[:,1], antenna_position[:,2] = -temp[:,1]/100., temp[:,0]/100., temp[:,2]/100.
   
    #temp=np.copy(peak_amplitude)
    #peak_amplitude[:,:,0], peak_amplitude[:,:,1] = -temp[:,:,1], temp[:,:,0]
   
    #temp=np.copy(peak_time)
    #peak_time[:,:,0], peak_time[:,:,1] = temp[:,:,1], temp[:,:,0]
   
    azimuth=3*np.pi/2+azimuth #  convert to +x = east (phi=0), +y = north (phi=90)
    ###
    #import pdb; pdb.set_trace()
    debug_XYZ_power = np.array(debug_XYZ_power)

    return (zenith, azimuth, energy, hillas, longprofile, Xground, antenna_position, onskypower, filteredpower, power, power11, power21, power41, peak_time, peak_amplitude, particle_radius, energy_deposit, pol_angle, pol_angle_filt, debug_XYZ_power, debug_figure)

   
