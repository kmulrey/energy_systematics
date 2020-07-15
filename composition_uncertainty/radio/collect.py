import os
import glob
import cPickle
import sys
from optparse import OptionParser
import numpy as np

BASE_PATH='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/t2b_events/'
RESULTS_PATH='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/radio/radio_results/'
NEWSIM_PATH='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/t2b_events/'
DATA_DIR=BASE_PATH+'events/'
SIMULATION_DIR=BASE_PATH+'run/'
OUTPUT_DIR_RADIO_ONLY=RESULTS_PATH+'/production_analysis_radio/'
WRITE_FILT='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/radio/events/'
COLLECT_DIR='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/radio/filtered/'
LORA_DIR='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/data_cal_final/'

eventid = 148663780
type='proton'
outputdir=COLLECT_DIR
datadir=BASE_PATH


print 'output dir: {0}'.format(outputdir)
print 'data dir: {0}'.format(datadir)


eventdir = os.path.join(datadir, "{0}".format(eventid))

print 'event dir: {0}'.format(eventdir)
for l in np.arange(1):
    print 'Event directory: %s' % (eventdir)
    #eventno=int(eventdir.split("/")[-1])
    filt_files=glob.glob(eventdir+"/coreas/{0}/DAT??????.filt".format(type))
    #filt_files_p=glob.glob(eventdir+"/coreas/proton/DAT??????.filt")
    #filt_files_He=glob.glob(eventdir+"/coreas/helium/DAT??????.filt")
    #filt_files_O=glob.glob(eventdir+"/coreas/oxygen/DAT??????.filt")
    #filt_files_Fe=glob.glob(eventdir+"/coreas/iron/DAT??????.filt")

    #nshow_p=len(filt_files_p)
    #nshow_Fe=len(filt_files_Fe)
    #nshow_He=len(filt_files_He)
    #nshow_O=len(filt_files_O)

    nshow=len(filt_files)
    print 'Number of showers found: %d' % nshow
    nantennas=160

    '''
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

    '''
