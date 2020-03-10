import os
import glob
from optparse import OptionParser
#import process as radio
import test_process as test_radio
import cPickle
import numpy as np
import matplotlib.pyplot as plt
import subprocess
plt.ion()

parser = OptionParser()
parser.add_option("-n", "--eventid", default = "0", help = "Event ID to process")
parser.add_option("-d", "--datadir", default="/vol/astro3/lofar/sim/pipeline/events", help="Directory where simulations are located")
parser.add_option("--lorafile-suffix", default="", help="Optional suffix for lora file, e.g. DAT000001_GeVfix.lora, to distinguish different runs of LORAsimulation (for testing)")
parser.add_option("--force-reprocess", default=False, action="store_true", help="Force reprocessing of simulated data")
parser.add_option("--test", default=False, action="store_true", help="Run debugging test")
parser.add_option("--debug-test-pulse", default=False, action="store_true", help="Replace simulated CoREAS data by a test pulse")
parser.add_option("--debug-lofar-pulse", default=False, action="store_true", help="Replace simulated CoREAS data by LOFAR data, or the processed test pulse from LOFAR-CR pipeline. Does not include LOFAR antenna positions (yet)")
parser.add_option("-w", "--writedir", default="/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/events", help="Directory where simulations are located")

(options, args) = parser.parse_args()
eventid = int(options.eventid)
datadir = options.datadir
doTestProcessFunction = options.test
force_reprocess = options.force_reprocess
lorafile_suffix = options.lorafile_suffix
writedir = options.writedir
lorafile_suff='lora'
# read in LOFAR calibrated pulse block, if desired with debug-lofar-pulse option
lofar_pulse = None
if options.debug_lofar_pulse:
    options.debug_test_pulse = True # if not set, set to True
    # Read database info, for crp_plotfiles key which contains the most actual results directory (i.e. the one used in  the database)
    import psycopg2 # for changing status in database
    from pycrtools import crdatabase as crdb
    dbManager = crdb.CRDatabase("crdb", host="coma00.science.ru.nl", user="crdb", password="crdb", dbname="crdb")
    db = dbManager.db
    print '### Replacing simulated data by LOFAR data processed by pipeline (calibrated-xyz) ###'
    print 'Reading event data, id = %d ...' % eventid
    event_data = crdb.Event(db=db, id=eventid)
    #event_data.stations[0].polarization['0']["crp_integrated_pulse_power"]
    station = [st for st in event_data.stations if st.stationname=='CS002'][0]
    crp_integrated_pulse_power_01 = station.polarization['0']['crp_integrated_pulse_power']
    crp_integrated_pulse_power_xyz = station.polarization['xyz']['crp_integrated_pulse_power']
    crp_plotfiles = event_data["crp_plotfiles"]
    results_dir = os.path.split(crp_plotfiles[0])[0]
    print 'Reading in timeseries data from directory %s' % results_dir
    print '(for now, reading in one antenna from CS002 and copying...)'
    xyz_timeseries = np.load(os.path.join(results_dir, 'xyz_calibrated_pulse_block-%d-CS002.npy' % eventid))

    lofar_pulse = xyz_timeseries # just take antenna 0, XYZ


#print 'done'


#eventdirs=glob.glob("/vol/astro3/lofar/sim/pipeline/events/*")

#for eventdir in eventdirs:
#     print eventdir
#     eventno=int(eventdir.split("/")[-1])
print '\n\nin filter jobs____________'
print datadir
print eventid
dirs=glob.glob(datadir + "/{0}/*/coreas/*".format(eventid))
print dirs
print dirs[0]
print dirs[0].split('events/')
print dirs[0].split('events')[1]
print writedir + dirs[0].split('events')[1]
new_dir= writedir + dirs[0].split('events')[1]

subprocess.check_output(['mkdir', '-p', new_dir])

print 'new_dir: ',new_dir

if len(dirs) == 0:
    raise ValueError("No directories with simulations found for event %d" % eventid)

for d in dirs:
    print 'starting loop'
    print d
    showerfiles=glob.glob(d+"/DAT??????")
    print showerfiles
    for showerfile in showerfiles:
        showerno=int(showerfile[-6:])
        #check if all files are ready
        
        
        
        print '\n\n'
        print os.path.isdir(d+"/SIM{0}_coreas".format(str(showerno).zfill(6)))
        print os.path.isfile(d+"/DAT{0}{1}.lora".format(str(showerno).zfill(6), lorafile_suffix))
        print os.stat(d+"/DAT{0}.long".format(str(showerno).zfill(6))).st_size>0
        print os.path.isfile(d+"/DAT{0}.long".format(str(showerno).zfill(6)))
        print os.path.isfile(d+"/steering/RUN{0}.inp".format(str(showerno).zfill(6)))
        print os.path.isfile(d+"/steering/SIM{0}.list".format(str(showerno).zfill(6)))
        print '\n\n'

        if (os.path.isdir(d+"/SIM{0}_coreas".format(str(showerno).zfill(6))) and
            os.path.isfile(d+"/DAT{0}{1}.lora".format(str(showerno).zfill(6), lorafile_suffix)) and
            (os.stat(d+"/DAT{0}.long".format(str(showerno).zfill(6))).st_size>0) and
            os.path.isfile(d+"/DAT{0}.long".format(str(showerno).zfill(6))) and
            os.path.isfile(d+"/steering/RUN{0}.inp".format(str(showerno).zfill(6))) and
            os.path.isfile(d+"/steering/SIM{0}.list".format(str(showerno).zfill(6))) ):

            outfile=new_dir+"/DAT{0}.filt".format(str(showerno).zfill(6))
            if (not os.path.isfile(outfile)) or (os.path.getsize(outfile) == 0) or force_reprocess:
                print 'Processing simulated data for event %d, shower %d, output dir %s, outfile %s' % (eventid, showerno, d, outfile)
                    #                try:
                (zenith, azimuth, energy, hillas, longprofile, Xground, antenna_position, onskypower, filteredpower, power, power11, power21, power41, peak_time, peak_amplitude, particle_radius, energy_deposit, pol_angle, pol_angle_filt, debug_XYZ_power, debug_lofarpulse_figure) = test_radio.ProcessData(d, showerno, lorafile_suffix=lorafile_suffix, debug_testpulse=options.debug_test_pulse, debug_lofarpulse=lofar_pulse)
                pickfile = open(outfile, 'w')
                cPickle.dump((zenith, azimuth, energy, hillas, longprofile, Xground, antenna_position, onskypower, filteredpower, power, power11, power21, power41, peak_time, peak_amplitude, particle_radius, energy_deposit, pol_angle, pol_angle_filt), pickfile)
                pickfile.close()
                
                if debug_lofarpulse_figure is not None:
                    plt.figure(debug_lofarpulse_figure.number)
                    plt.savefig('/vol/astro7/lofar/sim/pipeline/test_analysis_radio_only/DEBUG_lofarpulse_filtered_01_timeseries_%d.pdf' % eventid)
                    plt.figure()
                    nofantennas = len(crp_integrated_pulse_power_01.ravel())
                    x_axis = np.arange(len(crp_integrated_pulse_power_01.ravel()))
                    plt.scatter(x_axis, crp_integrated_pulse_power_01.ravel(), c='b', marker='s', label='Power from cr_physics')
                    plt.scatter(x_axis, power11.ravel()[0:nofantennas], c='r', marker='o', label='Power from sim pipeline')
                    # Make plot of power11 versus LOFAR crp_integrated_pulse_power for 0/1
                    plt.xlabel('Antenna number')
                    plt.ylabel('Signal energy [ J/m2 ]')
                    max_energy = np.max(crp_integrated_pulse_power_01.ravel())
                    plt.ylim(-0.2 * max_energy, 1.5 * max_energy)
                    plt.legend()
                    plt.tight_layout()
                    
                    plt.figure()
                    plt.scatter(x_axis, (power11.ravel()[0:nofantennas] - crp_integrated_pulse_power_01.ravel()) / crp_integrated_pulse_power_01.ravel(), c='r', label='Rel. difference cr_physics and sim pipeline 0/1 power')
                    plt.legend()
                    plt.xlabel('Antenna number (XYZ)')
                    plt.ylabel('Relative energy difference')
                    plt.tight_layout()
                                
                    # Make plot of XYZ power (sim pipeline) versus LOFAR crp_integrated_pulse_power for XYZ pol.
                    plt.figure()
                    nofantennas = len(crp_integrated_pulse_power_xyz.ravel())
                    x_axis = np.arange(len(crp_integrated_pulse_power_xyz.ravel()))
    
                    plt.scatter(x_axis, crp_integrated_pulse_power_xyz.ravel(), c='b', marker='s', label='Power from cr_physics')
                    plt.scatter(x_axis, debug_XYZ_power.ravel()[0:nofantennas], c='r', marker='o', label='Power from sim pipeline')
                    plt.xlabel('Antenna number')
                    plt.ylabel('Signal energy [ J/m2 ]')
                    max_energy = np.max(crp_integrated_pulse_power_xyz.ravel())
                    plt.ylim(-0.2 * max_energy, 1.5 * max_energy)
                    plt.legend()
                    plt.tight_layout()

                    plt.figure()
                    plt.scatter(x_axis, (debug_XYZ_power.ravel()[0:nofantennas] - crp_integrated_pulse_power_xyz.ravel()) / crp_integrated_pulse_power_xyz.ravel(), c='r', label='Rel. difference sim pipeline minus cr_physics 0/1 power')
                    plt.legend()
                    plt.xlabel('Antenna number (XYZ)')
                    plt.ylabel('Relative energy difference')
                    plt.tight_layout()
                    import pdb; pdb.set_trace()
                                
                                
                                              
                doFootprintTest = False
                if doFootprintTest:
                    # Make interpolated footprint plot of testProfile
                    import matplotlib.pyplot as plt
                    plt.ion()
                    plt.figure()
                    power_total = power11[:, 0] + power11[:, 1]
                    x_array = antenna_position[:, 0]
                    y_array = antenna_position[:, 1]
                    
                    plt.scatter(x_array, y_array, c=power_total)
                    
                    import pdb; pdb.set_trace()
                    """                    import scipy.interpolate as intp
                    plot_profile = 1.00000001 * testProfile / np.max(testProfile) # normalize to 1.0 max for plotting (a.u.)
                    Interpolation = intp.Rbf(x_array[distanceUpTo250Indices], y_array[distanceUpTo250Indices], plot_profile, smooth=0, function='quintic')
                    smooth_x_array = np.arange(-250.0, 250.0, 1.0)
                    smooth_y_array = np.arange(-250.0, 250.0, 1.0)
                    #            imarray_xy = np.array( (smooth_x_array, smooth_y_array) )
                    imarray = Interpolation(smooth_x_array, smooth_y_array)
                    imarray = np.zeros( (500, 500) )
                    for ii, x in enumerate(smooth_x_array):
                        for jj, y in enumerate(smooth_y_array):
                            imarray[ii, jj] = Interpolation(x, y)
                            plt.imshow(imarray.T, cmap=plt.cm.cubehelix_r, extent=[-250.0, 250.0, -250.0, 250.0], origin='lower')
                    cbar = plt.colorbar()
                    cbar.solids.set_rasterized(True)
                    plt.scatter(x_array, y_array, c='w', s=40)
                    plt.grid()
                    plt.xlabel('Meters in v x B direction')
                    plt.ylabel('Meters in v x (v x B) direction')
                    plt.savefig('/vol/astro2/users/acorstanje/test/footprint_lowerXmax.pdf')
                    """
                
                
                if doTestProcessFunction: # Test updates of the ProcessData function
                    standard_power11 = power11
                    standard_pol_angle = pol_angle
                    standard_pol_angle_filt = pol_angle_filt
                    
                    (zenith, azimuth, energy, hillas, longprofile, Xground, antenna_position, onskypower, filteredpower, power, power11, power21, power41, peak_time, peak_amplitude, particle_radius, energy_deposit, pol_angle, pol_angle_filt) = test_radio.ProcessData(d, showerno, outfile)
                    test_power11 = power11
                    test_pol_angle = pol_angle
                    test_pol_angle_filt = pol_angle_filt

                    import pdb; pdb.set_trace()

#                except:
#                    print "error!"

     

