import matplotlib
matplotlib.use('Agg')
import numpy as np
from optparse import OptionParser
import matplotlib.pyplot as plt
from matplotlib import cm
import cPickle
import scipy.interpolate as intp
import scipy.optimize as opt
from scipy.special import gamma
import random
import os
import glob
import ldf_db_rev_updated as dbldf # To re-fetch LOFAR data from CR Database

def setCommandlineOptions(cmdline_options):
    global options

    options = cmdline_options
    return

def setOptions(options): # either from command-line options, or through another script importing fit_analysis_updated
    global eventid, inputdir, iteration, outputdir, randomseed, doFetch, doRewrite, radio_only_fit
    # The use of 'global' here is wrong... need to change to Class instead.
    if type(options) == type( (1, 2) ): # if it is tuple, return it directly
        # 'options' must in this case be (eventid, inputdir, iteration, outputdir, randomseed, doFetch, doRewrite, radio_only_fit)
        (eventid, inputdir, iteration, outputdir, randomseed, doFetch, doRewrite, radio_only_fit) = options
        return

    eventid = int(options.event)
    
    inputdir = options.inputdir
    iteration = options.iteration
    if iteration == 'latest':
        inputfiles_alliterations = glob.glob(os.path.join(inputdir, 'filtered/SIM{0}*.filt'.format(eventid))) # Gets iterations 1, 2, ...
        inputfile_iterationzero = glob.glob(os.path.join(inputdir, 'filtered/SIM{0}.filt'.format(eventid))) # See if iteration zero is also there
        all_iterations = []
        if len(inputfile_iterationzero) == 1:
            all_iterations.append(0)
        if len(inputfiles_alliterations) == 0:
            raise ValueError("No file SIM<eventid>.filt found!")
        all_iterations.extend([int(file[1+file.find("_"):file.find(".")]) for file in inputfiles_alliterations if file.find("_") > 0]) # Get numbers from file names...

        iteration = max(all_iterations) # The latest one
        #for file in inputfiles_alliterations:
        #        if file.find("_") > 0:
        #            all_iterations.append(int(file[1+file.find("_"):file.find(".")]))
    else:
        iteration = int(iteration)
        
        outputdir = options.outputdir
        randomseed = int(options.randomseed)
        
        doFetch = options.fetch_lofardata
        doRewrite = options.rewrite_lofardata
        if doRewrite:
            doFetch = True # need to fetch in order to rewrite the file with LOFAR data

    radio_only_fit = options.radio_only_fit

    return (eventid, inputdir, iteration, outputdir, randomseed, doFetch, doRewrite, radio_only_fit)

if __name__ == "__main__": # if this is executed as script, not as module: run the analysis

    parser = OptionParser()
    parser.add_option("-n", "--event", default = "0", help = "filename of database")
    parser.add_option("-t", "--iteration", default="latest", help="Iteration number of the simulations to read results from, or 'latest' to get the latest one")
    parser.add_option("-i", "--inputdir", default="/vol/astro3/lofar/sim/pipeline/run", help="Base directory where input files are located, i.e. in <inputdir>/data and <inputdir>/filtered")
    parser.add_option("-o", "--outputdir", default="/vol/astro3/lofar/sim/pipeline/run/analysis")
    parser.add_option("--randomseed", default=2017, help="Set random seed for initial choice of core positions")

    parser.add_option("-f", "--fetch-lofardata", default=False, action="store_true", help="Re-fetch LOFAR measured data from CR Database (instead of from previously saved file)")
    parser.add_option("-w", "--rewrite-lofardata", default=False, action="store_true", help="Rewrite file containing LOFAR measured data (from CR Database)")
    parser.add_option("--radio-only-fit", default=False, action="store_true", help="Perform fit using only radio data; default is to use both radio and particle data")
    #parser.add_option("--test", default=False, action="store_true", help="Run debugging test")
    parser.add_option("--pickle-fancy-footprint", default=False, action="store_true", help="Produce pickle file to make fancy hi-res footprint")
    (options, args) = parser.parse_args()

    (eventid, inputdir, iteration, outputdir, randomseed, doFetch, doRewrite, radio_only_fit) = setOptions(options)


def setFetchRewrite(fetch, rewrite):
    global doFetch, doRewrite

    doFetch = fetch
    doRewrite = rewrite

    return

def GetPower_UVW(amp, power, cx,cy,zen,az):
    #first transform to vector
    vec=np.sign(amp)*np.sqrt(power)
    #no get UVW coordinates
    vec2=GetUVW(vec,cx,cy,zen,az)
    #transform to power
    return np.square(vec2)

def GetXYZ(pos, zen, az):
    inc=1.1837
    B = np.array([0,np.cos(inc),-np.sin(inc)])
    v = np.array([-np.cos(az)*np.sin(zen),-np.sin(az)*np.sin(zen),-np.cos(zen)])
    #print v
    vxB = np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]])
    vxB = vxB/np.linalg.norm(vxB)
    vxvxB = np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])
    return pos[0]*vxB+pos[1]*vxvxB+pos[2]*v

def GetUVW(pos, cx, cy, zen, az):
    relpos = pos-np.array([cx,cy,7.6])
    inc=1.1837
    B = np.array([0,np.cos(inc),-np.sin(inc)])
    v = np.array([-np.cos(az)*np.sin(zen),-np.sin(az)*np.sin(zen),-np.cos(zen)])
    #print v
    vxB = np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]])
    vxB = vxB/np.linalg.norm(vxB)
    vxvxB = np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])
    return np.array([np.inner(vxB,relpos),np.inner(vxvxB,relpos),np.inner(v,relpos)]).T

def readLOFARdata(datafile):
    global doFetch, doRewrite # they get updated if info from a data file is requested but doesn't exist
    # Uses parameters from the command line: eventid, doFetch, doRewrite
    if not doFetch and not os.path.exists(datafile): # Want to read from file but it is not there: fetch from database instead, and write to file.
        doFetch = True; doRewrite = True
        print 'LOFAR data file not found; fetching LOFAR data from CR Database and writing to file: %s' % datafile
    
    if doFetch: # Get info from CR Database using 'dbldf' module
        print 'Fetching LOFAR data from CR Database'
        nofAttempts = 10 # CR Database sometimes doesn't connect, so retry if this happens
        thisAttempt = 0
        readout_OK = False
        while not readout_OK and (thisAttempt < nofAttempts):
            try:
                core_x, core_y , stname, antenna_ids, positions, dist, x_err, signals, power11, power21, power41, rms, noisepower, pulse_delay_fit_residual, time, lora_x, lora_y, lora_dens, az, elev, elev_lora = dbldf.GetLDF(options.event)
                readout_OK = True
            except Exception as e:
                thisAttempt += 1
                print 'Database connection failed at attempt %d, error message is %s' % (thisAttempt, e)
                import time
                time.sleep(20 + 100*thisAttempt*np.random.rand() ) # sleep for 20 to 120 seconds before retrying
        if not readout_OK:
            raise ValueError("Could not get a database connection after %d attempts!" % nofAttempts)

        # Convert az/el to az/zenith in radians and with different convention (which?)
        lofarData = (core_x, core_y , stname, antenna_ids, positions, dist, x_err, signals, power11, power21, power41, rms, noisepower, pulse_delay_fit_residual, time, lora_x, lora_y, lora_dens, (450-az)/180.*np.pi, (90-elev)/180.*np.pi,(90-elev_lora)/180.*np.pi)
        if doRewrite: # Rewrite 'datafile' with values just read in from CR Database
            print 'Rewriting file containing LOFAR data'
            f = open(datafile, 'w')
            cPickle.dump((core_x, core_y , stname, antenna_ids, positions, dist, x_err, signals, power11, power21, power41, rms, noisepower, pulse_delay_fit_residual, time, lora_x, lora_y, lora_dens, (450-az)/180.*np.pi, (90-elev)/180.*np.pi,(90-elev_lora)/180.*np.pi),f)
            f.close()

    else: # Just read in from pre-produced file
        print 'Reading LOFAR data from pre-produced file'
        f = open(datafile,'r')
        lofarData = cPickle.load(f)
        f.close()
    
    return lofarData


def reverseAnalysis(eventno, iteration, inputdir, outputdir, randomseed, flagging=True, plots=True, verbose=True, outfile="reverse", simmode=False, simevent=0, saveplt=False, showplt=True, simcorex=0, simcorey=0):
   
    print "start"

    realxmax=0
    random.seed(randomseed) # Take a fixed random seed to be reproducable, can be updated in command line parameter
    np.random.seed(randomseed)
    
    twentyone=False

    datafile = os.path.join(inputdir, 'data/dbrev{0}.dat'.format(eventno))
    print 'using datafile {0}'.format(datafile)
    iterationSuffix = '_{0}'.format(iteration) if iteration > 0 else ''
    # read simulation results with this iteration number; no iteration number specified if it is zero.
    simfile = os.path.join(inputdir, 'filtered/SIM{0}{1}.filt'.format(eventno, iterationSuffix))
    fitoutfile = os.path.join(outputdir, 'reco{0}{1}.dat'.format(eventno, iterationSuffix))

    plotfile = os.path.join(outputdir, '{0}_{1}a{2}.png'.format(outfile, eventno, iterationSuffix))
    plotfile2 = os.path.join(outputdir, '{0}_{1}b{2}.png'.format(outfile, eventno, iterationSuffix))
    plotfile3 = os.path.join(outputdir, '{0}_{1}c{2}.png'.format(outfile, eventno, iterationSuffix))
    plotfile_catalog = os.path.join(outputdir, '{0}_{1}catalog{2}.pdf'.format(outfile, eventno, iterationSuffix))
    (core_x, core_y, station_name, antenna_ids, positions, dist, x_err, signal, dpower11, dpower21, dpower41, rms, noisepower, pulse_delay_fit_residual, data_time, lora_x, lora_y, lora_dens, data_azimuth, data_zenith, lora_zenith) = readLOFARdata(datafile)
    print '\n\n LORA DATA: {0}  \n\n'.format(lora_dens)

#    f=open(datafile,'r')
#    core_x, core_y , station_name, positions, dist, x_err, signal, dpower11, dpower21, dpower41, rms, noisepower, data_time, lora_x, lora_y, lora_dens, data_azimuth, data_zenith, lora_zenith = cPickle.load(f)

    if (simmode):
        core_x = simcorex
        core_y = simcorey

    dpower=dpower11
    samples=11

    if (twentyone): 
        spower=dpower21
        samples=21

    antenna_ids = [id for id in antenna_ids if int(id[len(id)-2::]) % 2 == 0] # Identify antennas by their 0-dipole.
    antenna_ids = np.array(antenna_ids)

    print 'There are %d antenna ids' % len(antenna_ids)

    nant=len(dist)
    print 'There are %d antennas' % nant
    assert nant == len(antenna_ids)
    
    #import pdb; pdb.set_trace()
    strongest_dipole=0
    if (np.max(dpower[:,1])>np.max(dpower[:,0])): strongest_dipole=1
    # Take only this dipole for pulse delay fit residual
    pulse_delay_fit_residual = pulse_delay_fit_residual[strongest_dipole::2]

    dtotpower = (dpower[:,0]+dpower[:,1])
    max_dtotpower=np.max(dtotpower)
    dpower_nonzero = (np.abs(dpower)+dpower)/2.   
    dsigma = 1.9 * np.sqrt(4*dpower_nonzero*noisepower+ 2*samples*noisepower*noisepower) # dpower = data power for 1 polarization, noisepower = power of noise per bin in oone polarization
    dsigma = np.sqrt(dsigma*dsigma+np.square(0.05*dpower)) #account for 5% fluctuations in gain from antenna to antenna

    dsigmatot = np.sqrt(dsigma[:,0]*dsigma[:,0]+dsigma[:,1]*dsigma[:,1]) 

    #dsigmatot = np.sqrt(dsigmatot*dsigmatot + 0.025*0.025*np.max(dtotpower)*np.max(dtotpower)) #account for 2.5% max interpolation error
    dmeanpower = np.sum(dtotpower)/len(dtotpower)
    dmeannoise = np.array([np.sum(noisepower[:,0])/len(noisepower),np.sum(noisepower[:,1])/len(noisepower)])
    
    ant_sel=np.ones([nant],dtype=bool)
    if (flagging):
        if verbose: print 'Flagging antennas according to a 3 sigma criterion, a minimum power and a pulse timing within 8 ns of the plane-wave fit (per station)...'
        for i in np.unique(station_name):
            #print station_name
            selection=station_name==i
            stpower=dtotpower[selection]
            stmed=np.median(stpower)
            uprange=np.percentile(stpower,80)-stmed
            downrange=stmed-np.percentile(stpower,20)
            
            ant_sel[selection]=(stpower>(stmed-3*downrange))*(stpower<(stmed+3*uprange))*(dtotpower[selection]>0.02*max_dtotpower)
        if verbose: print "Antennas flagged: ", nant-np.sum(ant_sel), " out of ", nant
        nof_flagged = nant - np.sum(ant_sel)
        # Of the remaining antennas, flag all of them that have pulse timings not within 8 ns
        ant_sel_timedelays = (np.abs(pulse_delay_fit_residual) < 8.0e-9)
        ant_sel *= ant_sel_timedelays

        nof_flagged_timedelays = nant - np.sum(ant_sel) - nof_flagged
        if verbose: print "Antennas flagged additionally, due to outlying pulse timing: %d" % nof_flagged_timedelays

    nsel_ant = np.sum(ant_sel)

    nstations=len(lora_x)

    lora_positions=np.zeros([nstations,3])
    lora_positions[:,0]=lora_x
    lora_positions[:,1]=lora_y

    lora_dens=np.array(lora_dens)

    eff_area=0.9*np.cos(lora_zenith) # 0.9 m2 area
    lora_err=np.sqrt(lora_dens*eff_area)/eff_area
    for i in np.arange(nstations):
        if (lora_err[i]>lora_dens[i]): lora_err[i]=lora_dens[i]
      
    lora_err=lora_err[np.nonzero(lora_dens)]
    lora_positions=lora_positions[np.nonzero(lora_dens)]
    lora_dens=lora_dens[np.nonzero(lora_dens)]

    ConversionFactor = eff_area * 6.7 # Density => #particles => energy deposit (MeV)
    lora_dens=lora_dens*ConversionFactor 
    lora_err=lora_err*ConversionFactor

    nstations=len(lora_dens)

    if verbose: print 'Reading in LORA info'
    g=open(simfile,'r')
    siminfo = cPickle.load(g)

    zenith=siminfo['zenith']
    azimuth=siminfo['azimuth']
    energy=siminfo['energy']
    hillas=siminfo['hillas'].T
    longprofile=siminfo['longprofile']
    Xground=siminfo['Xground']
    sim_antenna_position=siminfo['antenna_position']
    sim_power=siminfo['totpower']
    sim_power11=siminfo['power11bins']
    sim_power21=siminfo['power21bins']
    sim_power41=siminfo['power41bins']
    sim_time=siminfo['peak_time']
    sim_ampl=siminfo['peak_amplitude']
    cors_r=siminfo['particle_radius']
    cors_dens=siminfo['energy_deposit']
    primary_type=siminfo['primary'][:,0]

    nsim_prot=np.sum(primary_type==14)

    sim_power=sim_power11
    if (twentyone): sim_power=sim_power21

    sim_tot_power=np.sum(sim_power,axis=2)
    nsim=sim_tot_power.shape[0]
    nsimant=sim_tot_power.shape[1]

    data_zenith=zenith[0]
    data_azimuth=azimuth[0]

    pos_sim_UVW = np.zeros([nsim,nsimant,3])
      
    rbf= np.ndarray([nsim],dtype=object)
    rbf_0= np.ndarray([nsim],dtype=object)
    rbf_1= np.ndarray([nsim],dtype=object)

    if verbose: print 'Interpolate simulations'
    for i in np.arange(nsim):
        pos_sim_UVW[i,:,:] = GetUVW(sim_antenna_position[i,:,:],0,0, zenith[0], azimuth[0])
        selection=np.array(np.isfinite(sim_power[i,:,0]))*np.array(np.isfinite(sim_power[i,:,1]))
        rbf[i] = intp.Rbf(pos_sim_UVW[i,selection,0], pos_sim_UVW[i,selection,1], sim_tot_power[i,selection],smooth =0,function='quintic')
        rbf_0[i] = intp.Rbf(pos_sim_UVW[i,selection,0], pos_sim_UVW[i,selection,1], sim_power[i,selection,0],smooth =0,function='quintic')
        rbf_1[i] = intp.Rbf(pos_sim_UVW[i,selection,0], pos_sim_UVW[i,selection,1], sim_power[i,selection,1],smooth =0,function='quintic')

    part_intp=np.ndarray([nsim],dtype=object)

    if (verbose): print "#    N_ch    R_m    s"
    for i in np.arange(nsim):
        r_range=np.concatenate([[0],cors_r[i]])
        dens_range=np.concatenate([[cors_dens[i,0]],cors_dens[i]])
        part_intp[i]=intp.interp1d(r_range, dens_range)

    # NKG fit function
    def lora_fitfunc(ratio,lpos,cx,cy,az,zen,n): 
        pos_lora_UVW = GetUVW(lpos, cx, cy, zen, az)
        r = np.sqrt(pos_lora_UVW[:,0]*pos_lora_UVW[:,0]+pos_lora_UVW[:,1]*pos_lora_UVW[:,1])
        if (np.max(r)<995): return ratio*part_intp[n](r)
        else: return 0

    lora_errfunc_restrained = lambda p,lofar_pow,lofar_err, lofar_pow01,lofar_err01, lofar_pos,lora_data, lora_err, lora_pos, cx, cy, az, zen, n: (lora_fitfunc(p[0], lora_pos,cx+p[1],cy+p[2],az,zen,n) - lora_data)/lora_err

    lora_errfunc_fixedcore_restrained = lambda p, cx_offset, cy_offset, lofar_pow,lofar_err, lofar_pow01,lofar_err01, lofar_pos,lora_data, lora_err, lora_pos, cx, cy, az, zen, n: (lora_fitfunc(p, lora_pos,cx+cx_offset,cy+cy_offset,az,zen,n) - lora_data)/lora_err # cx_offset and cy_offset are the offsets to core position from radio (or combined) fit


    # 2D radio fit function
    def radio_fitfunc(f,dpos,n,cx,cy,az,zen):
        pos_ant_UVW = GetUVW(dpos, cx, cy, zen, az)
        interp_power = rbf[n](pos_ant_UVW[:,0],pos_ant_UVW[:,1])
        return f*interp_power

    # 2D radio fit function   
    def radio_fitfunc01(f0, f1, dpos,n,cx,cy,az,zen):
        pos_ant_UVW = GetUVW(dpos, cx, cy, zen, az)
        interp_power0 = rbf_0[n](pos_ant_UVW[:,0],pos_ant_UVW[:,1])
        interp_power1 = rbf_1[n](pos_ant_UVW[:,0],pos_ant_UVW[:,1])
        return f0*interp_power0+f1*interp_power1

    # 2D radio fit function   
    def radio_fitfunc0(f,dpos,n,cx,cy,az,zen):
        pos_ant_UVW = GetUVW(dpos, cx, cy, zen, az)
        interp_power = rbf_0[n](pos_ant_UVW[:,0],pos_ant_UVW[:,1])
        return f*interp_power

    # 2D radio fit function   
    def radio_fitfunc1(f,dpos,n,cx,cy,az,zen):
        pos_ant_UVW = GetUVW(dpos, cx, cy, zen, az)
        interp_power = rbf_1[n](pos_ant_UVW[:,0],pos_ant_UVW[:,1])
        return f*interp_power

    # radio fit p=[pratio,xoff,yoff]   
    radio_errfunc = lambda p,lofar_pow,lofar_err, lofar_pow01,lofar_err01, lofar_pos,lora_data, lora_err, lora_pos, cx, cy, az, zen, n: (radio_fitfunc(p[0],lofar_pos,n,cx+p[1],cy+p[2],az,zen) - lofar_pow)/lofar_err

    # combined error function p = [pratio,dratio,xoff,yoff] 
    def combined_errfunc(p,*args):
        a = radio_errfunc([p[0],p[2],p[3]],*args)
        b = lora_errfunc_restrained([p[1],p[2],p[3]],*args)
        c = np.concatenate((a,b))
        return c

    #repeat original lora fit
    #lora_orig_fitparam, covar_lora_orig = opt.leastsq(lora_errfunc_free_cf,[1e8,50], args=(lora_dens,lora_err,lora_positions,core_x,core_y,data_azimuth,data_zenith,0))

    combchi2=np.zeros([nsim])
    radiochi2=np.zeros([nsim])
    radiochi2_d=np.zeros([nsim])
    lorachi2=np.zeros([nsim])
    p_ratio=np.zeros([nsim])
    p_ratio0=np.zeros([nsim])
    p_ratio1=np.zeros([nsim])
    d_ratio=np.zeros([nsim])
    xoffset=np.zeros([nsim])
    yoffset=np.zeros([nsim])
    simpower=np.zeros([nant,nsim])
    simpower0=np.zeros([nant,nsim])
    simpower1=np.zeros([nant,nsim])


    if (simmode):
        realxmax=hillas[2,simevent]
        pos_ant_UVW = GetUVW(positions, core_x, core_y, data_zenith, data_azimuth) #HACK: get pos... change core?
        relerror=dsigmatot/dtotpower
        dtotpower=rbf[simevent](pos_ant_UVW[:,0],pos_ant_UVW[:,1]) # HACK: set power to simulation #simevent ... add noise!
        dsigmatot=dtotpower*relerror
        for i in np.arange(nant):
            dtotpower[i]=dtotpower[i]+random.gauss(0,dsigmatot[i])
        lora_dens = lora_fitfunc(1.0,lora_positions,core_x, core_y,data_azimuth,data_zenith,simevent)
        lora_err = np.sqrt(lora_dens)
        for i in np.arange(nstations):
            lora_dens[i]=lora_dens[i]+random.gauss(0,lora_err[i])


    def FitRoutine_radio_only(fit_args,iterations=1, dosim=False, debug=False):
        paramset=np.ndarray([iterations,3]) # 3 free parameters in fit
        itchi2=np.zeros(iterations)
        for j in np.arange(iterations):
            if verbose: print 'Fit routine, iteration %d...' % j
            initguess=[3e9, 200*np.random.rand()-100, 200*np.random.rand()-100] # Radio power scale factor, core X, core Y
            if (dosim): initguess=[1,200*np.random.rand()-100,200*np.random.rand()-100] #CHANGE VALUES ACCRODING TO MODE!!!
            paramset[j], covar = opt.leastsq(radio_errfunc, initguess, args=fit_args)
            itchi2[j] = np.sum(np.square(radio_errfunc(paramset[j],*fit_args)))
            if (verbose): print itchi2[j], paramset[j]

        if debug:
            import pdb; pdb.set_trace()
        bestiter=np.argmin(itchi2)
        fitparam=paramset[bestiter]
        # For the optimal radio fit, obtain the optimal d_ratio (# particles scale factor)
        #print 'Doing d_Ratio fit'
        initguess_d_ratio = 1.0
        (xoffset, yoffset) = (fitparam[1], fitparam[2])
        particle_fit_args = (xoffset, yoffset) + fit_args # concatenates tuples
        optimal_d_ratio, covar = opt.leastsq(lora_errfunc_fixedcore_restrained, initguess_d_ratio, args=particle_fit_args)
        optimal_d_ratio = optimal_d_ratio[0]
        combined_fit_param = [fitparam[0], optimal_d_ratio, fitparam[1], fitparam[2]]
        
        chi2_comb=np.sum(np.square(combined_errfunc(combined_fit_param, *fit_args))) # Get combined chi2 from radio fit + optimal d_ratio after radio fit
        chi2_rad=np.sum(np.square(radio_errfunc(fitparam, *fit_args)))
        chi2_part=np.sum(np.square(lora_errfunc_restrained([optimal_d_ratio, fitparam[1], fitparam[2]],*fit_args)))
        
        return chi2_comb, chi2_rad, chi2_part, fitparam[0], optimal_d_ratio, fitparam[1], fitparam[2]



    def FitRoutine(fit_args,iterations=1, dosim=False, debug=False):
        paramset=np.ndarray([iterations,4])
        itchi2=np.zeros(iterations)
        for j in np.arange(iterations):
            if verbose: print 'Fit routine, iteration %d...' % j
            initguess=[3e9,1,200*np.random.rand()-100,200*np.random.rand()-100] # Radio power scale factor, # particles scale factor, core X, core Y
            if (dosim): initguess=[1,1,200*np.random.rand()-100,200*np.random.rand()-100] #CHANGE VALUES ACCRODING TO MODE!!!
            paramset[j], covar = opt.leastsq(combined_errfunc,initguess, args=fit_args)
            itchi2[j]=np.sum(np.square(combined_errfunc(paramset[j],*fit_args)))
            if (verbose): print itchi2[j], paramset[j]

        if debug:
            import pdb; pdb.set_trace()
        bestiter=np.argmin(itchi2)
        fitparam=paramset[bestiter]
        chi2_comb=np.sum(np.square(combined_errfunc(fitparam,*fit_args)))
        chi2_rad=np.sum(np.square(radio_errfunc([fitparam[0],fitparam[2],fitparam[3]],*fit_args)))
        chi2_part=np.sum(np.square(lora_errfunc_restrained([fitparam[1],fitparam[2],fitparam[3]],*fit_args)))
        return chi2_comb, chi2_rad, chi2_part, fitparam[0], fitparam[1], fitparam[2], fitparam[3]

    if verbose: print 'Run fit routine...'
    Desired_FitRoutine = FitRoutine_radio_only if radio_only_fit else FitRoutine

    niterations=3
    for i in np.arange(nsim): # Loop over nsim, i.e. number of simulated showers. Do fit with core as free parameters, to interpolated LDF footprint of this shower. Get chi^2, repeat for all showers.
        if verbose: print 'Simulated shower %d out of %d' % (i, nsim)
        fit_args=(dtotpower[ant_sel],dsigmatot[ant_sel],dpower[ant_sel],dsigma[ant_sel],positions[ant_sel],lora_dens,lora_err,lora_positions,core_x,core_y,data_azimuth,data_zenith,i)
        
        #debug = False if i != 27 else True
        
        # Do radio-only fit if option given on command line
        #print 'Radio only: %s' % radio_only_fit
        combchi2[i], radiochi2[i], lorachi2[i], p_ratio[i], d_ratio[i], xoffset[i], yoffset[i] = Desired_FitRoutine(fit_args,niterations, simmode)

    Desired_Chi2 = radiochi2 if radio_only_fit else combchi2

    if radio_only_fit:
        print "Selecting radio chi2 for parabola fitting"
    else:
        print "Selecting combined chi2 for parabola fitting"


    if (simmode):
        combchi2[simevent]=1000
        radiochi2[simevent]=1000
    bestsim=np.argmin(Desired_Chi2)
    xoff=xoffset[bestsim]
    yoff=yoffset[bestsim]
    #import pdb; pdb.set_trace()

    #ADDITIONAL FLAGGING
    if verbose: print "Additional flagging for outliers from the fit, redo the fit afterwards."
    fitparam=np.array([p_ratio[bestsim],d_ratio[bestsim],xoffset[bestsim],yoffset[bestsim]])
    fit_args=(dtotpower,dsigmatot,dpower,dsigma,positions,lora_dens,lora_err,lora_positions,core_x,core_y,data_azimuth,data_zenith,bestsim)
    rad_rel_err= np.square(radio_errfunc([p_ratio[bestsim],xoffset[bestsim],yoffset[bestsim]],*fit_args))
    ant_sel=ant_sel*(rad_rel_err<16)

    nsel_ant = np.sum(ant_sel)
    if verbose: print "Antennas flagged: ", nant-np.sum(ant_sel), " out of ", nant

    ndf_lora=nstations-4
    ndf_comb = nsel_ant+nstations-5
    ndf_radio = nsel_ant-4
    ndf_radio_d = 2*nsel_ant-4    
    if verbose: print 'Re-doing fits without flagged antennas'
    
    for i in np.arange(nsim):
        if verbose: print 'Simulated shower %d out of %d' % (i, nsim)
        fit_args=(dtotpower[ant_sel],dsigmatot[ant_sel],dpower[ant_sel],dsigma[ant_sel],positions[ant_sel],lora_dens,lora_err,lora_positions,core_x,core_y,data_azimuth,data_zenith,i)
        combchi2[i], radiochi2[i], lorachi2[i], p_ratio[i], d_ratio[i], xoffset[i], yoffset[i] = Desired_FitRoutine(fit_args,niterations, simmode)

        if (verbose): print i, p_ratio[i], hillas[2,i], xoffset[i], yoffset[i], combchi2[i], radiochi2[i]

    p_ratio0=p_ratio
    p_ratio1=p_ratio

    if (simmode):
        combchi2[simevent]=1000
        radiochi2[simevent]=1000
    bestsim=np.argmin(Desired_Chi2)
    xoff=xoffset[bestsim]
    yoff=yoffset[bestsim]

    #import pdb; pdb.set_trace()

    # Fixed: take core position of best-fitting simulation for the evaluation of simulated power at each LOFAR antenna.
    # To be used for radio anti-bias cut and making plots.
    for i in np.arange(nsim):
        simpower[:,i]=radio_fitfunc(p_ratio[i],positions,i,core_x+xoff,core_y+yoff,azimuth[0],zenith[0])
        simpower0[:,i]=radio_fitfunc0(p_ratio[i],positions,i,core_x+xoff,core_y+yoff,azimuth[0],zenith[0])
        simpower1[:,i]=radio_fitfunc1(p_ratio[i],positions,i,core_x+xoff,core_y+yoff,azimuth[0],zenith[0])

    # AC: store radio error (LOFAR_data - fit) in output for later analysis per antenna (syst. deviations).
    fit_residuals = radio_errfunc([p_ratio[bestsim],xoffset[bestsim],yoffset[bestsim]],*fit_args)



    pos_ant_UVW = GetUVW(positions, core_x+xoff, core_y+yoff, data_zenith, data_azimuth)
    axdist_ant = np.sqrt(pos_ant_UVW[:,0]*pos_ant_UVW[:,0]+pos_ant_UVW[:,1]*pos_ant_UVW[:,1])
    orig_pos_lora_UVW = GetUVW(lora_positions, core_x, core_y, data_zenith, data_azimuth)
    orig_axdist_lora = np.sqrt(orig_pos_lora_UVW[:,0]*orig_pos_lora_UVW[:,0]+orig_pos_lora_UVW[:,1]*orig_pos_lora_UVW[:,1])
    pos_lora_UVW = GetUVW(lora_positions, core_x+xoff, core_y+yoff, data_zenith, data_azimuth)
    axdist_lora = np.sqrt(pos_lora_UVW[:,0]*pos_lora_UVW[:,0]+pos_lora_UVW[:,1]*pos_lora_UVW[:,1])

    no150 = np.argmin(np.abs(axdist_ant - 150))

    if verbose: print "starting parabola fit procedure, selection of points and actual fit"
    Desired_ndf = ndf_radio if radio_only_fit else ndf_comb

    urange=hillas[2,bestsim]+100
    drange=hillas[2,bestsim]-100
    chirange = Desired_Chi2[bestsim] + 0.5 * Desired_ndf

    # Two methods to select points to participate in fit:

    fit_selection=np.zeros(nsim) # all points that have lower chi2 values on on side only
    for i in np.arange(nsim):
        if (np.sum(Desired_Chi2[(hillas[2,:]>hillas[2,i])] < Desired_Chi2[i])==0): fit_selection[i]=fit_selection[i]+1
        if (np.sum(Desired_Chi2[(hillas[2,:]<hillas[2,i])] < Desired_Chi2[i])==0): fit_selection[i]=fit_selection[i]+1

    fit_selection=fit_selection*(hillas[2,:]>drange)*(hillas[2,:]<urange)*(Desired_Chi2<chirange)
    if (simmode): fit_selection[simevent]=0

    fo=2
    chi2fitparam_1=np.zeros([3])
    fitconverged_1=False
    if (np.sum(fit_selection)>fo):
        chi2fitparam_1, res_1, rank_1, singval_1, rcond_1 = np.polyfit(hillas[2,(fit_selection>0)],Desired_Chi2[(fit_selection>0)],fo,full=True)
        pfit_1 = np.poly1d(chi2fitparam_1)
        r1=pfit_1.deriv().r
        if (r1[(r1>drange)*(r1<urange)].size<3): 
            if (r1[(r1>drange)*(r1<urange)].size>0): fitconverged_1=True

    xmaxreco=0

    if (fitconverged_1):
        xmaxreco = r1[(r1>drange)*(r1<urange)][0] # ?

    print 'Reconstructed Xmax = %4.3f g/cm2' % xmaxreco

    def CScale(a,i):
        b = (a-np.min(a))/(np.max(a)-np.min(a))
        return (b[i],1-b[i],0)                  

    if (plots):
        if verbose: print 'Now doing the first panel plot...'
        ###
        if options.pickle_fancy_footprint:
            print 'Creating bypass pickle file with input params for footprint-on-ground plot...'
            bypass_file = '/vol/astro7/lofar/sim/pipeline/test_fancy_footprint/%d_input_params.pickle' % eventno
            
            outfile = open(bypass_file, 'w')
            plot_input_params = (positions, core_x, core_y, xoff, yoff, sim_tot_power, bestsim, nsimant, sim_antenna_position, p_ratio, dtotpower, lora_positions, lora_dens)
            
            cPickle.dump(plot_input_params, outfile)
            outfile.close()
        ###
                     
        f, ((ax1, ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2,3,figsize=(15,10))

        dist_scale=1.2*np.max(np.sqrt((positions[:,0]-core_x-xoff)*(positions[:,0]-core_x-xoff)+(positions[:,1]-core_y-yoff)*(positions[:,1]-core_y-yoff)))
        ti = np.linspace(-dist_scale, dist_scale, 150)
        XI, YI = np.meshgrid(ti, ti)
        selection = np.unique(np.array(np.isfinite(sim_tot_power[bestsim,:]))*np.arange(nsimant))
        rbf_ground = intp.Rbf(sim_antenna_position[bestsim,selection,0], sim_antenna_position[bestsim,selection,1], sim_tot_power[bestsim,selection], smooth =0, function='quintic')
        ZI = rbf_ground(XI, YI)*p_ratio[bestsim]
        maxp = np.max([np.max(dtotpower),np.max(ZI)])
        ax1.pcolor(XI+core_x+xoff, YI+core_y+yoff, ZI,vmax=maxp, vmin=0,cmap=cm.jet)
        dist_scale=1.2*np.max(np.sqrt((positions[:,0]-core_x-xoff)*(positions[:,0]-core_x-xoff)+(positions[:,1]-core_y-yoff)*(positions[:,1]-core_y-yoff)))

        ax1.scatter(positions[:,0],positions[:,1],50,c=dtotpower,vmax=maxp, vmin=0,cmap=cm.jet)
        ax1.scatter(sim_antenna_position[bestsim,:,0]+core_x+xoff,sim_antenna_position[bestsim,:,1]+core_y+yoff,10, c=sim_tot_power[bestsim,:]*p_ratio[bestsim],vmax=maxp, vmin=0,cmap=cm.jet)
        #ax1.scatter(newp[:,0]+core_x+xoff,newp[:,1]+core_y+yoff,10, c=totpower*p_ratio[bestsim],vmax=maxp, vmin=0,cmap=cm.jet)
        ax1.scatter(lora_positions[:,0],lora_positions[:,1],250*lora_dens/np.max(lora_dens),marker='p',c='white')
        ax1.scatter(core_x,core_y,marker='o')
        ax1.scatter(core_x+xoff,core_y+yoff,marker='+')
        ax1.set_xlim((-dist_scale+core_x+xoff,dist_scale+core_x+xoff))
        ax1.set_ylim((-dist_scale+core_y+yoff,dist_scale+core_y+yoff))
        ax1.set_title("Footprint on ground")
        ax1.set_xlabel("West-East (m)")
        ax1.set_ylabel("South-North (m)")


        ax2.scatter(hillas[2,:],Desired_Chi2[:],50,color='b')
        if (nsim_prot!=nsim): ax2.scatter(hillas[2,nsim_prot:], Desired_Chi2[nsim_prot:],50,color='m')
        ax2.scatter(hillas[2,bestsim], Desired_Chi2[bestsim],50,color='r')
        #ax2.scatter(hillas[2,:],radiochi2[:],50,marker='^',color='g')
        #ax2.scatter(hillas[2,:],lorachi2[:],50,marker='p', color='black')
        #ax2.scatter(hillas[2,:],lorachi2[:],50,marker='p', color='black')
        #ax2.axhline(combchi2[bestsim]+1,linestyle='--',color='r')
        #ax2.set_yscale('log')
        ax2.set_title("Fit Quality after core optimization")
        ax2.set_ylabel("reduced Chi^2")
        ax2.set_xlabel(r"X_max ($g/cm^2$)")

        dist_scale = 1.2*np.max(axdist_ant)
        ti = np.linspace(-dist_scale, dist_scale, 150)
        XI, YI = np.meshgrid(ti, ti)

        ZI_power_showerplane = rbf[bestsim](XI, YI) * p_ratio[bestsim]
        maxp = np.max([np.max(dtotpower),np.max(ZI)])
        
        def plotRadiationProfileInShowerPlane(ax=None, eventIDtext=None):
            # A function to be able to do a plot more than once without duplicate code...
            ax.pcolor(XI, YI, ZI_power_showerplane, vmax=maxp, vmin=0,cmap=cm.jet, linewidth=0, rasterized=True)
            ax.scatter(pos_sim_UVW[bestsim,:,0],pos_sim_UVW[bestsim,:,1],10,sim_tot_power[bestsim,:]*p_ratio[bestsim], vmax=maxp, vmin=0,cmap=cm.jet)
            #ax.scatter(newp_UVW[:,0],newp_UVW[:,1],10,totpower[bestsim,:]*p_ratio[bestsim], vmax=maxp, vmin=0,cmap=cm.jet)
            ax.scatter(pos_ant_UVW[:,0],pos_ant_UVW[:,1],50,dtotpower, vmax=maxp, vmin=0,cmap=cm.jet)
            ax.scatter(0,0,marker='+')
            ax.set_xlim((-dist_scale,dist_scale))
            ax.set_ylim((-dist_scale,dist_scale))
            #ax.set_title("Radiation profile in shower plane")
            ax.set_xlabel("vxB (m)")
            ax.set_ylabel("vx(vxB) (m)")
            if eventIDtext is not None:
                ax.text(-0.9*dist_scale, 0.9*dist_scale, eventIDtext, color='white')
                ax.text(-0.9*dist_scale, 0.8*dist_scale, 'Sim E = %1.2e eV' % (energy[0]*1.0e9), color='white', fontsize=10)
                ax.text(-0.9*dist_scale, 0.7*dist_scale,  'd_ratio         = %1.3f' % d_ratio[bestsim], color='white', fontsize=10)
                ax.text(-0.9*dist_scale, 0.6*dist_scale,  'sqrt(p_ratio) = %1.3f' % np.sqrt(p_ratio[bestsim]), color='white', fontsize=10)

        plotRadiationProfileInShowerPlane(ax3)

        """
        # HACK: save this plot separately as pdf for demonstration purposes
        #        extent = ax3.get_window_extent().transformed(plt.gcf().dpi_scale_trans.inverted())
        #fig.savefig('ax2_figure.png', bbox_inches=extent)
        this_figure = plt.gcf()
        new_fig = plt.figure()
        plt.pcolor(XI, YI, ZI,vmax=maxp, vmin=0,cmap=cm.jet)
        plt.scatter(pos_sim_UVW[bestsim,:,0],pos_sim_UVW[bestsim,:,1],10,sim_tot_power[bestsim,:]*p_ratio[bestsim], vmax=maxp, vmin=0,cmap=cm.jet)
        #ax3.scatter(newp_UVW[:,0],newp_UVW[:,1],10,totpower[bestsim,:]*p_ratio[bestsim], vmax=maxp, vmin=0,cmap=cm.jet)
        plt.scatter(pos_ant_UVW[:,0],pos_ant_UVW[:,1],50,dtotpower, vmax=maxp, vmin=0,cmap=cm.jet)
        plt.scatter(0,0,marker='+')
        plt.xlim((-dist_scale,dist_scale))
        plt.ylim((-dist_scale,dist_scale))
        #ax3.set_title("Radiation profile in shower plane")
        plt.xlabel("Position along vxB (m)")
        plt.ylabel("Position along vx(vxB) (m)")
        plt.axes().set_aspect('equal', 'box')
        mm = cm.ScalarMappable(cmap=cm.jet)
        mm.set_array(ZI / 1e-15)
        cbar = plt.colorbar(mm)
        cbar.set_label('Intensity (a.u.)', rotation=90)
        # Save figure
        plt.savefig('demo_plot_LDF_sim_lofar.png', dpi=300)#, bbox_inches=extent.expanded(1.4, 1.3))
        
        plt.figure(this_figure)
        """
        def plotLDFradio(ax=None, plot_title=False):
            ax.errorbar(axdist_ant[ant_sel],dtotpower[ant_sel],dsigmatot[ant_sel],linestyle='',marker='o',color='r')
            if (nant-nsel_ant>0): ax4.errorbar(axdist_ant[ant_sel<1],dtotpower[ant_sel<1],dsigmatot[ant_sel<1],linestyle='',marker='x',color='black')
            ax.plot(axdist_ant[ant_sel],simpower[ant_sel,bestsim],linestyle='',marker='o',color='b')
            if plot_title:
                ax.set_title("Lateral distribution radio signal")
            ax.set_ylabel("Total intensity")
            ax.set_xlabel("Distance (m)")

        plotLDFradio(ax=ax4, plot_title=True)

        ax5.errorbar(orig_axdist_lora,lora_dens, lora_err, linestyle='',marker='o',color='r')
        ax5.errorbar(axdist_lora,lora_dens, lora_err, linestyle='',marker='o',color='b')
        hival = np.max([np.max(orig_axdist_lora),np.max(axdist_lora)])
        dval=np.linspace(10,hival,100)
        fitrange=np.min([hival,995])
        fitval=np.linspace(0,fitrange,100)
        #ax5.plot(dval,particle_fitfunc([lora_orig_fitparam[0],lora_orig_fitparam[1],1.7],dval), color='r')
        #ax5.plot(dval,particle_fitfunc([d_ratio[bestsim],mol_radius[bestsim],sh_age[bestsim]],dval), color='b')
        ax5.plot(fitval,d_ratio[bestsim]*part_intp[bestsim](fitval),color='b')
        ax5.set_yscale('log')
        ax5.set_title("LORA ldf after core optimization")
        ax5.set_ylabel(r"particle density ($\mathrm{m}^-2$)")
        ax5.set_xlabel("distance (m)")
        #ax5.text(0.2,0.9,'orig. $N_e$={0:.2e}, $R_m$={1:.2f}'.format(lora_orig_fitparam[0],lora_orig_fitparam[1]),transform=ax5.transAxes,color='r')
        #ax5.text(0.8,0.8,'$s$=1.7',transform=ax5.transAxes,color='r')
        #ax5.text(0.2,0.7,'new  $N_e$={0:.2e}, $R_m$={1:.2f}'.format(d_ratio[bestsim],mol_radius[bestsim]),transform=ax5.transAxes,color='b')
        #ax5.text(0.75,0.6,'$s$={0:.2f}'.format(sh_age[bestsim]),transform=ax5.transAxes,color='b')

        ax6.text(0.1,0.8,'event no. {0}'.format(eventno))
        ax6.text(0.1,0.7,'zenith: {0:.2f}, azimuth: {1:.2f}'.format(data_zenith/np.pi*180, np.fmod(data_azimuth/np.pi*180,360)))
        #ax6.text(0.1,0.6,'sim zenith: {0:.2f}, azimuth: {1:.2f}'.format(zenith[0]/np.pi*180, azimuth[0]/np.pi*180))
        ax6.text(0.1,0.6,'energy: {0:.2e} GeV (x {1:.1e})'.format(energy[0],d_ratio[bestsim]))
        ax6.text(0.1,0.5,'original core: ({0:.2f}, {1:.2f}) m'.format(core_x,core_y))
        ax6.text(0.1,0.4,'new core: ({0:.2f}, {1:.2f}) m'.format(xoff+core_x, yoff+core_y))
        ax6.text(0.1,0.3,'radio red. chi2: {0:.3f} ({1:.3f})'.format(radiochi2[bestsim]/(ndf_radio+1e-25),radiochi2_d[bestsim]/(ndf_radio_d+1e-25)))
        ax6.text(0.1,0.2,'comb red. chi2: {0:.3f}'.format(combchi2[bestsim]/(ndf_comb+1e-25)))
        ax6.text(0.1,0.1,'Xmax: {0:.0f} g/cm^2'.format(hillas[2,bestsim]))
        ax6.set_xticks([])
        ax6.set_yticks([])
        if (saveplt): plt.savefig(plotfile, dpi=300)
        plt.close()

        if verbose: print 'Doing second plot...'
        f, ((ax10, ax11,ax12),(ax13,ax14,ax15)) = plt.subplots(2,3,figsize=(15,10))

        dist_scale = 1.2*np.max(axdist_ant)
        ti = np.linspace(-dist_scale, dist_scale, 150)
        XI, YI = np.meshgrid(ti, ti)
        ZI = rbf[bestsim](XI, YI)*p_ratio[bestsim]
        maxp = np.max([np.max(dtotpower),np.max(ZI)]) 
        ax10.pcolor(XI, YI, ZI,vmax=maxp, vmin=0,cmap=cm.jet)
        ax10.scatter(pos_sim_UVW[bestsim,:,0],pos_sim_UVW[bestsim,:,1],10,sim_tot_power[bestsim,:]*p_ratio[bestsim], vmax=maxp, vmin=0,cmap=cm.jet)
        #ax3.scatter(newp_UVW[:,0],newp_UVW[:,1],10,totpower[bestsim,:]*p_ratio[bestsim], vmax=maxp, vmin=0,cmap=cm.jet)
        ax10.scatter(pos_ant_UVW[:,0],pos_ant_UVW[:,1],50,dtotpower, vmax=maxp, vmin=0,cmap=cm.jet)
        ax10.scatter(0,0,marker='+')
        ax10.set_xlim((-dist_scale,dist_scale))
        ax10.set_ylim((-dist_scale,dist_scale))
        ax10.set_title("Radiation profile in shower plane")
        ax10.set_xlabel("vxB (m)")
        ax10.set_ylabel("vx(vxB) (m)")

        dist_scale = 1.2*np.max(axdist_ant)
        ti = np.linspace(-dist_scale, dist_scale, 150)
        XI, YI = np.meshgrid(ti, ti)
        ZI = rbf_0[bestsim](XI, YI)*p_ratio[bestsim]  
        maxp = np.max([np.max(dpower[:,0]),np.max(ZI)]) 
        ax11.pcolor(XI, YI, ZI,vmax=maxp, vmin=0,cmap=cm.jet)
        ax11.scatter(pos_sim_UVW[bestsim,:,0],pos_sim_UVW[bestsim,:,1],10,sim_power[bestsim,:,0]*p_ratio0[bestsim], vmax=maxp, vmin=0,cmap=cm.jet)
        #ax3.scatter(newp_UVW[:,0],newp_UVW[:,1],10,totpower[bestsim,:]*p_ratio[bestsim], vmax=maxp, vmin=0,cmap=cm.jet)
        ax11.scatter(pos_ant_UVW[:,0],pos_ant_UVW[:,1],50,dpower[:,0], vmax=maxp, vmin=0,cmap=cm.jet)
        ax11.scatter(0,0,marker='+')
        ax11.set_xlim((-dist_scale,dist_scale))
        ax11.set_ylim((-dist_scale,dist_scale))
        ax11.set_title("Radiation profile in shower plane")
        ax11.set_xlabel("vxB (m)")
        ax11.set_ylabel("vx(vxB) (m)")

        dist_scale = 1.2*np.max(axdist_ant)
        ti = np.linspace(-dist_scale, dist_scale, 150)
        XI, YI = np.meshgrid(ti, ti)
        ZI = rbf_1[bestsim](XI, YI)*p_ratio[bestsim]
        maxp = np.max([np.max(dpower[:,1]),np.max(ZI)]) 
        ax12.pcolor(XI, YI, ZI,vmax=maxp, vmin=0,cmap=cm.jet)
        ax12.scatter(pos_sim_UVW[bestsim,:,0],pos_sim_UVW[bestsim,:,1],10,sim_power[bestsim,:,1]*p_ratio1[bestsim], vmax=maxp, vmin=0,cmap=cm.jet)
        #ax3.scatter(newp_UVW[:,0],newp_UVW[:,1],10,totpower[bestsim,:]*p_ratio[bestsim], vmax=maxp, vmin=0,cmap=cm.jet)
        ax12.scatter(pos_ant_UVW[:,0],pos_ant_UVW[:,1],50,dpower[:,1], vmax=maxp, vmin=0,cmap=cm.jet)
        ax12.scatter(0,0,marker='+')
        ax12.set_xlim((-dist_scale,dist_scale))
        ax12.set_ylim((-dist_scale,dist_scale))
        ax12.set_title("Radiation profile in shower plane")
        ax12.set_xlabel("vxB (m)")
        ax12.set_ylabel("vx(vxB) (m)")

        maxp = np.max([np.max(dtotpower),np.max(simpower[:,bestsim])])
        ax13.errorbar(axdist_ant[ant_sel],dtotpower[ant_sel]/maxp,dsigmatot[ant_sel]/maxp,linestyle='',marker='o',color='r')
        if (nant-nsel_ant>0): ax13.errorbar(axdist_ant[ant_sel<1],dtotpower[ant_sel<1]/maxp,dsigmatot[ant_sel<1]/maxp,linestyle='',marker='x',color='black')
        ax13.plot(axdist_ant[ant_sel],simpower[ant_sel,bestsim]/maxp,linestyle='',marker='o',color='b')
        ax13.set_title("Lateral distribution radio signal TOTAL")
        ax13.set_ylabel("total power (a.u.)")
        ax13.set_xlabel("distance (m)")

        ax14.errorbar(axdist_ant[ant_sel],dpower[ant_sel,0]/maxp,dsigma[ant_sel,0]/maxp,linestyle='',marker='o',color='r')
        if (nant-nsel_ant>0): ax14.errorbar(axdist_ant[ant_sel<1],dpower[ant_sel<1,0]/maxp,dsigma[ant_sel<1,0]/maxp,linestyle='',marker='x',color='black')
        ax14.plot(axdist_ant[ant_sel],simpower0[ant_sel,bestsim]/maxp,linestyle='',marker='o',color='b')
        ax14.set_title("Lateral distribution radio signal 0")
        ax14.set_ylabel("total power (a.u.)")
        ax14.set_xlabel("distance (m)")

        ax15.errorbar(axdist_ant[ant_sel],dpower[ant_sel,1]/maxp,dsigma[ant_sel,1]/maxp,linestyle='',marker='o',color='r')
        if (nant-nsel_ant>0): ax15.errorbar(axdist_ant[ant_sel<1],dpower[ant_sel<1,1]/maxp,dsigma[ant_sel<1,1]/maxp,linestyle='',marker='x',color='black')
        ax15.plot(axdist_ant[ant_sel],simpower1[ant_sel,bestsim]/maxp,linestyle='',marker='o',color='b')
        ax15.set_title("Lateral distribution radio signal 1")
        ax15.set_ylabel("total power (a.u.)")
        ax15.set_xlabel("distance (m)")
        if (saveplt): plt.savefig(plotfile3, dpi=300)
        plt.close()

        if verbose: print 'Doing third plot...'
        f, ((ax16, ax17,ax18),(ax19,ax20,ax21)) = plt.subplots(2,3,figsize=(15,10))
        #ax16.scatter(sh_age,hillas[2,:])
        #ax16.scatter(sh_age,hillas[2,:])
        # Diagnostic plot to review the effect of core shifts on radio chi2 and Xmax ??
        ax16.scatter(xoffset[:],yoffset[:],200/radiochi2[:]*ndf_radio,vmin=np.min(hillas[2,:]), vmax=np.max(hillas[2,:]),c=hillas[2,:])
        ax16.scatter(xoffset[:],yoffset[:],200/radiochi2[:]*ndf_radio,vmin=np.min(hillas[2,:]), vmax=np.max(hillas[2,:]),c=hillas[2,:])
        if (nsim_prot!=nsim): ax16.scatter(xoffset[nsim_prot:],yoffset[nsim_prot:],200/Desired_Chi2[nsim_prot:],vmin=np.min(hillas[2,:]), vmax=np.max(hillas[2,:]),c=hillas[2,nsim_prot:],marker='^')

        for i in np.arange(nsim):
            ax17.plot(longprofile[i,:,1]/hillas[0,i],c=CScale(Desired_Chi2,i))

        t = np.linspace(drange, urange, 100)

        radius=np.ndarray([nsim])
        for i in np.arange(nsim):
            #selection=np.array(np.isfinite(sim_power[i,:,0]))*np.array(np.isfinite(sim_power[i,:,1]))
            hiSNR = ( (sim_tot_power[i,:]*p_ratio[bestsim]<9*dmeannoise[1]) * np.isfinite(sim_tot_power[i,:]) )
            n=(np.argmax(hiSNR))
            #print sim_tot_power[i,n]
            radius[i]=np.sqrt(pos_sim_UVW[i,n,0]*pos_sim_UVW[i,n,0]+pos_sim_UVW[i,n,1]*pos_sim_UVW[i,n,1])
        minrad=np.min(radius)
        #print radius, minrad
        hiSNR = ( (sim_tot_power[bestsim,:]*p_ratio[bestsim]<9*dmeannoise[1]) * np.isfinite(sim_tot_power[bestsim,:]) )
        n=(np.argmax(hiSNR))
        bestradius=np.sqrt(pos_sim_UVW[i,n,0]*pos_sim_UVW[i,n,0]+pos_sim_UVW[i,n,1]*pos_sim_UVW[i,n,1])

        dist_scale = 1.2*bestradius
        ti = np.linspace(-dist_scale, dist_scale, 150)
        XI, YI = np.meshgrid(ti, ti)
        ZI = rbf[bestsim](XI, YI)*p_ratio[bestsim]
        maxp = np.max([np.max(dtotpower),np.max(ZI)]) 


        ax18.pcolor(XI, YI, ZI,vmax=maxp, vmin=0,cmap=cm.jet)
        #ax18.scatter(pos_sim_UVW[bestsim,hiSNR,0],pos_sim_UVW[bestsim,hiSNR,1],30,sim_tot_power[bestsim,hiSNR]*p_ratio[bestsim], vmax=maxp, vmin=0,cmap=cm.jet)
        #ax18.scatter(pos_sim_UVW[bestsim,:,0],pos_sim_UVW[bestsim,:,1],10,sim_tot_power[bestsim,:]*p_ratio[bestsim], vmax=maxp, vmin=0,cmap=cm.jet)
        circle=plt.Circle((0,0),minrad,color='r',fill=False)
        circle2=plt.Circle((0,0),bestradius,color='g',fill=False, ls='dashed')
        ax18.add_patch(circle)
        ax18.add_patch(circle2)
        #ax3.scatter(newp_UVW[:,0],newp_UVW[:,1],10,totpower[bestsim,:]*p_ratio[bestsim], vmax=maxp, vmin=0,cmap=cm.jet)
        ax18.scatter(pos_ant_UVW[:,0],pos_ant_UVW[:,1],50,dtotpower, vmax=maxp, vmin=0,cmap=cm.jet)
        ax18.scatter(0,0,marker='+')
        ax18.set_xlim((-dist_scale,dist_scale))
        ax18.set_ylim((-dist_scale,dist_scale))
        ax18.set_title("Radiation profile in shower plane")
        ax18.set_xlabel("vxB (m)")
        ax18.set_ylabel("vx(vxB) (m)")


        if (fitconverged_1): # Plot the reduced chi^2 vs Xmax and the parabola fit
            ax19.scatter(hillas[2,:],Desired_Chi2[:]/(Desired_ndf + 1e-25),50,color='b')
            ax19.scatter(hillas[2,(fit_selection>0)],Desired_Chi2[(fit_selection>0)]/(Desired_ndf + 1e-25),50,color='m')
            ax19.plot(t,pfit_1(t)/Desired_ndf,'-',color='m')
            ax19.set_ylim((0.8*np.min(Desired_Chi2[(fit_selection>0)]/(Desired_ndf+1e-25)),1.2*np.max(Desired_Chi2[(fit_selection>0)]/(Desired_ndf+1e-25))))
            ax19.set_xlim(0.8*np.min(hillas[2,(fit_selection>0)]),1.2*np.max(hillas[2,(fit_selection>0)]))

        dist_scale = 1.2*np.max(axdist_ant)
        ti = np.linspace(-dist_scale, dist_scale, 150)
        XI, YI = np.meshgrid(ti, ti)
        bestmap = np.ones([150,150])
        totalmap = np.ones([150,150])
        SNRant=np.zeros([nsim,nant])
        hiSNRsel=np.ones([nant],dtype=bool)
        lora_dens_all=np.zeros([nsim,nstations])

        #import pdb; pdb.set_trace()
        
        for i in np.arange(nsim):
            ZI = np.maximum(rbf_0[i](XI, YI),rbf_1[i](XI, YI))*p_ratio[bestsim]
            map=(ZI>6*np.max(dmeannoise))
            if (i==bestsim): bestmap=map
            totalmap*=map
            SNRant[i]=np.maximum(simpower0[:,i],simpower1[:,i])/np.max(dmeannoise) # (AC) p_ratio factor is taken into account through 'radio_fitfunc', used to produce simpower0 / 1 (see line 445)
            hiSNRsel*=SNRant[i]>6
            lora_dens_all[i] = lora_fitfunc(d_ratio[bestsim],lora_positions,core_x+xoff, core_y+yoff, data_azimuth, data_zenith, i)

        pulses_per_station=[]
        for i in np.unique(station_name):
            selection=station_name==i
            if verbose: print "Station: %s has %d pulses that pass radio-bias cut" % (i, np.sum(hiSNRsel[selection]))
            pulses_per_station.append(np.sum(hiSNRsel[selection]))


        lora_dens_best = lora_fitfunc(d_ratio[bestsim],lora_positions,core_x+xoff, core_y+yoff, data_azimuth, data_zenith, bestsim)
        lora_dens_worst = np.min(lora_dens_all,axis=0)

        ax20.pcolor(XI, YI, bestmap,vmax=1, vmin=0,cmap=cm.jet)
        ax20.scatter(pos_ant_UVW[:,0],pos_ant_UVW[:,1],20,c='black')
        ax20.scatter(pos_ant_UVW[ant_sel,0],pos_ant_UVW[ant_sel,1],20,c='red')
        ax20.scatter(pos_ant_UVW[SNRant[bestsim]>9,0],pos_ant_UVW[SNRant[bestsim]>9,1],20,c='green')
        ax20.scatter(pos_lora_UVW[:,0],pos_lora_UVW[:,1],40,c='red')
        ax20.scatter(pos_lora_UVW[lora_dens_best>7,0],pos_lora_UVW[lora_dens_best>7,1],40,c='green')

        ax21.pcolor(XI, YI, totalmap,vmax=1, vmin=0,cmap=cm.jet)
        ax21.scatter(pos_ant_UVW[:,0],pos_ant_UVW[:,1],20,c='black')
        ax21.scatter(pos_ant_UVW[ant_sel,0],pos_ant_UVW[ant_sel,1],20,c='red')
        ax21.scatter(pos_ant_UVW[hiSNRsel,0],pos_ant_UVW[hiSNRsel,1],20,c='green')
        ax21.scatter(pos_lora_UVW[:,0],pos_lora_UVW[:,1],40,c='red')
        ax21.scatter(pos_lora_UVW[lora_dens_worst>7,0],pos_lora_UVW[lora_dens_worst>7,1],40,c='green')

        if (saveplt): plt.savefig(plotfile2, dpi=300)

        if verbose: print 'Doing catalog plot...'
        
        f, (ax31, ax32, ax33) = plt.subplots(1,3,figsize=(15,5))
        plotRadiationProfileInShowerPlane(ax31, eventIDtext='ID %d' % eventid)
        # plot reduced chi2 with parabola on larger scale than ax19
        Xmax_axis = hillas[2, :]
        Xmax_outer_envelope_points = hillas[2, (fit_selection>0)]
        reduced_chi2 = Desired_Chi2[:] / (Desired_ndf + 1e-25)
        reduced_chi2_outer_envelope_points = Desired_Chi2[(fit_selection>0)]/(Desired_ndf + 1e-25)
        
        ax32.scatter(Xmax_axis, reduced_chi2, 50, color='b')
            
        ax32.scatter(Xmax_outer_envelope_points, reduced_chi2_outer_envelope_points, 50, color='m')
        if fitconverged_1:
            ax32.plot(t,pfit_1(t)/Desired_ndf,'-',color='m', lw=2)

        best_reduced_chi2 = np.min(reduced_chi2)
        
        ax32.set_ylim((0.6*best_reduced_chi2), 3.5*best_reduced_chi2)
        #ax32.set_xlim(0.8*np.min(Xmax_outer_envelope_points)),1.2*np.max(Xmax_outer_envelope_points)
        if len(Xmax_outer_envelope_points) > 0:
            lower_xmaxlim = min(550.0, min(Xmax_outer_envelope_points))
            upper_xmaxlim = max(850.0, max(Xmax_outer_envelope_points))
        else:
            lower_xmaxlim = 550.0
            upper_xmaxlim = 850.0
        ax32.set_xlim(lower_xmaxlim, upper_xmaxlim)
        ax32.set_ylabel("Reduced $\chi^2$")
        ax32.set_xlabel(r"$X_{\mathrm{max}}$ (g/$\mathrm{cm}^2$)")
        ax32.text(1.03*lower_xmaxlim, 0.95*(3.5*best_reduced_chi2), "Xmax reco = %4.1f" % xmaxreco)
    
        plotLDFradio(ax=ax33)
        if (saveplt): plt.savefig(plotfile_catalog) #, dpi=300)

        if (showplt): plt.show()
        plt.close()
        if verbose: print 'Plotting completed.'

    if verbose: print 'Writing reconstruction info...'
    fitconverged=np.array([fitconverged_1])
    analysisinfo={'antenna_ids': antenna_ids[ant_sel], 'lofar_power': dtotpower[ant_sel], 'lofar_sigma': dsigmatot[ant_sel], 'fit_residuals': fit_residuals, 'fitconverged':fitconverged, 'xmaxreco':xmaxreco, 'xmaxbest':hillas[2,bestsim], 'xmax': hillas[2,:], 'realxmax':realxmax, 'combchi2':combchi2[bestsim]/(ndf_comb+1e-25), 'radiochi2': radiochi2[bestsim]/(ndf_radio+1e-25), 'all_radiochi2': radiochi2 / (ndf_radio+1e-25), 'fitparam':chi2fitparam_1, 'p_ratio':p_ratio[bestsim], 'p_ratio0':p_ratio0[bestsim], 'pratio1': p_ratio1[bestsim], 'd_ratio':d_ratio[bestsim], 'd_ratios':d_ratio, 'simtotpower': sim_tot_power, 'rbf':rbf[bestsim](xoff+core_x, yoff+core_y), 'dtotpower':dtotpower[no150], 'xoff':xoff, 'yoff':yoff, 'core_x':core_x, 'core_y':core_y, 'energy':energy[0], 'zenith':zenith[0], 'azimuth':azimuth[0], 'dmeannoise':dmeannoise, 'ant_fraction':1.0*nsel_ant/nant, 'nof_flagged_timedelays': nof_flagged_timedelays, 'nsim':nsim, 'nsim_prot':nsim_prot, 'SNRant': SNRant, 'hiSNRsel': hiSNRsel, 'pulses_per_station': pulses_per_station, 'lora_dens_best': lora_dens_best, 'lora_dens_worst': lora_dens_worst, 'lora_dens': lora_dens, 'station_name': station_name}
    fout=open(fitoutfile,"w")
    cPickle.dump(analysisinfo,fout)

    if verbose: print 'Analysis completed.'
    return fitconverged, xmaxreco, hillas[2,bestsim], hillas[2,:], realxmax, combchi2[bestsim]/(ndf_comb+1e-25), radiochi2[bestsim]/(ndf_radio+1e-25), chi2fitparam_1, p_ratio[bestsim], p_ratio0[bestsim], p_ratio1[bestsim], d_ratio[bestsim], d_ratio, sim_tot_power,rbf[bestsim](xoff+core_x, yoff+core_y),dtotpower[no150], xoff, yoff, core_x, core_y, energy[0], zenith[0], azimuth[0], dmeannoise, 1.0*nsel_ant/nant, nsim, nsim_prot        

if __name__ == "__main__": # if this is executed as script, not as module: run the analysis
    reverseAnalysis(eventid, iteration, inputdir, outputdir, randomseed, outfile="reco", saveplt=True, verbose=True, showplt=False, plots=True)
 
