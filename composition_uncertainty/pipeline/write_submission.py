import numpy as np
from optparse import OptionParser
import os
import os.path
from os import path
import pickle
import re
import sys
from os import listdir
from os.path import isfile, join
import glob
import __future__



parser = OptionParser()
#parser.add_option("-i", "--eventindex", type="int", help="event number", default="0")
parser.add_option("-e", "--event", type="int", help="event number", default="110395704")

options, arguments = parser.parse_args()
event=options.event


proton_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/pipeline/jobs_proton/'
helium_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/pipeline/jobs_helium/'
oxygen_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/pipeline/jobs_oxygen/'
iron_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/pipeline/jobs_iron/'
max_sim=15
write_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/pipeline/run_directory/corsika_runs/'


proton_list=np.genfromtxt('xmax_lists/runs_'+str(int(event))+'_proton.txt')
helium_list=np.genfromtxt('xmax_lists/runs_'+str(int(event))+'_helium.txt')
oxygen_list=np.genfromtxt('xmax_lists/runs_'+str(int(event))+'_oxygen.txt')
iron_list=np.genfromtxt('xmax_lists/runs_'+str(int(event))+'_iron.txt')

use_P_list=proton_list[0:max_sim]
use_He_list=helium_list[0:max_sim]
use_O_list=oxygen_list[0:max_sim]
use_Fe_list=iron_list[0:max_sim]

file=open(write_dir+'run_{0}.q'.format(event),'w')
file.write('#! /bin/bash\n')
file.write('#SBATCH --time=1-00:00:00\n')
file.write('cd {0}\n'.format(proton_dir))

file.write('sbatch --array ')
for i in np.arange(len(use_P_list)-1):
        file.write('{0},'.format(int(proton_list[i][0])))
file.write('{1} {0}_corsika_proton.q\n'.format(event,int(proton_list[i+1][0])))

file.write('sbatch --array ')
for i in np.arange(len(use_He_list)-1):
        file.write('{0},'.format(int(helium_list[i][0])))
file.write('{1} {0}_corsika_helium.q\n'.format(event,int(helium_list[i+1][0])))

file.write('sbatch --array ')
for i in np.arange(len(use_O_list)-1):
        file.write('{0},'.format(int(oxygen_list[i][0])))
file.write('{1} {0}_corsika_oxygen.q\n'.format(event,int(oxygen_list[i+1][0])))

file.write('sbatch --array ')
for i in np.arange(len(use_Fe_list)-1):
        file.write('{0},'.format(int(iron_list[i][0])))
file.write('{1} {0}_corsika_iron.q\n'.format(event,int(iron_list[i+1][0])))

file.close()
