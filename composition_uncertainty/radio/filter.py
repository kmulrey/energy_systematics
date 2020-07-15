import numpy as np
import glob, os, sys
from multiprocessing import Pool
from optparse import OptionParser
import sys

import process as process

BASE_PATH='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/t2b_events/'
RESULTS_PATH='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/radio/radio_results/'
NEWSIM_PATH='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/t2b_events/'
DATA_DIR=BASE_PATH+'events/'
SIMULATION_DIR=BASE_PATH+'run/'
OUTPUT_DIR_RADIO_ONLY=RESULTS_PATH+'/production_analysis_radio/'
WRITE_FILT='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/radio/events/'
COLLECT_DIR='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/radio/filtered/'
LORA_DIR='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/data_cal_final/'

eventid = 122146757
datadir = NEWSIM_PATH
writedir = BASE_PATH
lorafile_suffix=''

print('\n\nin filter jobs____________')
print(datadir)
print(eventid)
dirs=glob.glob(datadir + "/{0}/coreas/*".format(eventid))
print(datadir + "/{0}/*/coreas/*".format(eventid))

for i in np.arange(1):
  d=dirs[i]
  print('starting loop')
  print(d)
  showerfiles=glob.glob(d+"/DAT??????")
  new_dir= writedir + d.split('events')[1]


  for showerfile in showerfiles:
      showerno=int(showerfile[-6:])
      #check if all files are ready
      if (os.path.isdir(d+"/SIM{0}_coreas".format(str(showerno).zfill(6))) and
          os.path.isfile(d+"/DAT{0}{1}.lora".format(str(showerno).zfill(6), lorafile_suffix)) and
          os.path.isfile(d+"/DAT{0}.long".format(str(showerno).zfill(6))) and
          os.path.isfile(d+"/steering/RUN{0}.inp".format(str(showerno).zfill(6))) and
          os.path.isfile(d+"/steering/SIM{0}.list".format(str(showerno).zfill(6))) ):

          outfile=new_dir+"/DAT{0}.filt".format(str(showerno).zfill(6))
          #if (not os.path.isfile(outfile)) or (os.path.getsize(outfile) == 0):
          #    print('Processing simulated data for event %d, shower %d, output dir %s, outfile %s' % (eventid, showerno, d, outfile))
          (zenith, azimuth, energy, hillas, longprofile, Xground, antenna_position, onskypower, filteredpower, power, power11, power21, power41, peak_time, peak_amplitude, particle_radius, energy_deposit, pol_angle, pol_angle_filt, debug_XYZ_power, debug_lofarpulse_figure) = process.ProcessData(d, showerno, lorafile_suffix=lorafile_suffix)
          pickfile = open(outfile, 'w')
          cPickle.dump((zenith, azimuth, energy, hillas, longprofile, Xground, antenna_position, onskypower, filteredpower, power, power11, power21, power41, peak_time, peak_amplitude, particle_radius, energy_deposit, pol_angle, pol_angle_filt), pickfile)
          pickfile.close()
      else:
          print('missing info')
