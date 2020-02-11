import os
import glob
import cPickle
import sys
from optparse import OptionParser
import numpy as np

parser = OptionParser()
parser.add_option("-n", "--eventid", default = "0", help = "Event ID to process")
parser.add_option("-i", "--iteration", default="latest", help="Iteration number of simulation to process, or 'all' to process all of them, or 'latest' to only process the latest")
parser.add_option("-d", "--outputdir", default="/vol/astro3/lofar/sim/pipeline/run/filtered", help="Output directory of combined shower data")
parser.add_option("--datadir", default="/vol/astro3/lofar/sim/pipeline/events", help="Base dir where the simulated events are stored")
#parser.add_option("--force-reprocess", default=False, action="store_true", help="Force reprocessing of simulated data")
#parser.add_option("--test", default=False, action="store_true", help="Run debugging test")

(options, args) = parser.parse_args()
eventid = int(options.eventid)
iteration = options.iteration
if iteration != "all" and iteration != "latest":
    iteration = int(iteration)
print 'Iteration setting: %s' % iteration
outputdir = options.outputdir
datadir = options.datadir

print 'output dir: {0}'.format(outputdir)
print 'data dir: {0}'.format(datadir)


def RepresentsInt(s): # check if string is convertible to 'int'...
    try:
        int(s)
        return True
    except ValueError:
        return False

eventdir = os.path.join(datadir, "{0}".format(eventid))

print 'event dir: {0}'.format(eventdir)

# Get all available iterations
all_iterations = glob.glob(os.path.join(eventdir, '*'))
all_iterations = [int(os.path.split(x)[1]) for x in all_iterations if RepresentsInt(os.path.split(x)[1])] # make integers
if iteration == 'all':
    iterations_to_process = all_iterations
elif iteration == 'latest':
    iterations_to_process = [max(all_iterations)]
elif iteration in all_iterations:
    iterations_to_process = [iteration]
else:
    raise ValueError('Given iteration not found!')

print 'Going to process iteration numbers: '
print iterations_to_process

#eventdirs=glob.glob("/vol/astro3/lofar/sim/pipeline/events/*")
#for eventdir in eventdirs:
#try:
iterationCompleted = False # check if at least one iteration completes
for thisIteration in iterations_to_process:
    print 'Event directory: %s, iteration %d' % (eventdir, thisIteration)
    #eventno=int(eventdir.split("/")[-1])
    filt_files_p=glob.glob(eventdir+"/{0}/coreas/proton/DAT??????.filt".format(thisIteration))
    filt_files_Fe=glob.glob(eventdir+"/{0}/coreas/iron/DAT??????.filt".format(thisIteration))
    nshow_p=len(filt_files_p)
    nshow=nshow_p+len(filt_files_Fe)
    print 'Number of showers found: %d' % nshow
    nantennas=160
    if nshow==0:
        if len(iterations_to_process) > 1:
            print 'No showers found! Continuing with next iteration'
            continue
        else:
            raise ValueError('No showers found in iteration %d! Quitting.' % thisIteration)

        #print 'No showers found! Exiting.'
        #sys.exit(1)
    #    if (nshow>0):
    antenna_position=np.zeros([nshow,nantennas,3])
    onskypower=np.zeros([nshow,nantennas,2])
    filteredpower=np.zeros([nshow,nantennas,2])
    power=np.zeros([nshow,nantennas,2])
    power11=np.zeros([nshow,nantennas,2])
    power21=np.zeros([nshow,nantennas,2])
    power41=np.zeros([nshow,nantennas,2])
    peak_time=np.zeros([nshow,nantennas,2]) 
    peak_amplitude=np.zeros([nshow,nantennas,2]) 
    zenith=np.zeros([nshow])
    azimuth=np.zeros([nshow])
    energy=np.zeros([nshow])
    hillas=np.zeros([nshow,6])
    longprofile=np.zeros([nshow,4000,3]) # NB. Changed from 400 to 4000 for 1 g/cm2 longitudinal resolution!
    Xground=np.zeros([nshow])
    particle_radius=np.zeros([nshow,200])
    energy_deposit=np.zeros([nshow,200])
    pol_angle=np.zeros([nshow,nantennas])
    pol_angle_filt=np.zeros([nshow,nantennas])
    ptype=np.zeros([nshow,nantennas],dtype=int)

    i=0
    for filt_file in filt_files_p:
       f=open(filt_file,"r")
       zenith[i], azimuth[i], energy[i], hillas[i], longprofile[i], Xground[i], antenna_position[i], onskypower[i], filteredpower[i], power[i], power11[i], power21[i], power41[i], peak_time[i], peak_amplitude[i], particle_radius[i], energy_deposit[i], pol_angle[i], pol_angle_filt[i] = cPickle.load(f)
       ptype[i]=14
       i+=1
    for filt_file in filt_files_Fe:
       f=open(filt_file,"r")
       zenith[i], azimuth[i], energy[i], hillas[i], longprofile[i], Xground[i], antenna_position[i], onskypower[i],filteredpower[i], power[i], power11[i], power21[i], power41[i], peak_time[i], peak_amplitude[i], particle_radius[i], energy_deposit[i], pol_angle[i], pol_angle_filt[i] = cPickle.load(f)
       ptype[i]=5626
       i+=1

    siminfo={'zenith':zenith,'azimuth':azimuth,'energy':energy,'hillas':hillas,'longprofile':longprofile,'Xground':Xground,
             'antenna_position':antenna_position, 'onskypower':onskypower, 'totpower':power, 'filteredpower':filteredpower, 'power11bins':power11, 'power21bins':power21, 'power41bins':power41,
             'peak_time': peak_time, 'peak_amplitude': peak_amplitude, 'particle_radius': particle_radius, 'energy_deposit': energy_deposit,
             'pol_angle': pol_angle, 'pol_angle_filt': pol_angle_filt, 'primary': ptype}

    iterationSuffix = '_{0}'.format(thisIteration) if thisIteration > 0 else ''
    # Suffix _1 for iteration 1, etc. No suffix in file name for iteration 0
    outfile = os.path.join(outputdir, "SIM{0}{1}.filt".format(eventid, iterationSuffix) )
    f = open(outfile, "wb")
    cPickle.dump(siminfo, f)
    f.close()
    iterationCompleted = True


if not iterationCompleted:
    raise ValueError('No valid iterations to process!')

