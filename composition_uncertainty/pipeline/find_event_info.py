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
parser.add_option("-i", "--eventindex", type="int", help="event number", default="0")

options, arguments = parser.parse_args()
eventindex=options.eventindex




reco_dir='/vol/astro7/lofar/sim/pipeline/production_analysis_Dec2019/'
event_list=np.genfromtxt('/home/kmulrey/cr_physics_new')
dat_files=glob.glob(reco_dir+'*dat')
event=event_list[eventindex]
for file in glob.glob(reco_dir+'*'+str(int(event))+'*dat'):
    reco_file_name=file


reco_file=open(reco_file_name,"rb")
reco_info=pickle.load(reco_file)#, encoding='latin1')
reco_file.close()


core_x=reco_info['core_x']+reco_info['xoff']
core_y=reco_info['core_y']+reco_info['yoff']
sim_energy=reco_info['energy']
xmax_fit=reco_info['xmax']


dir3='/vol/astro7/lofar/sim/pipeline/events/'+str(int(event))+'/3/coreas/proton/steering/'
dir2='/vol/astro7/lofar/sim/pipeline/events/'+str(int(event))+'/2/coreas/proton/steering/'
dir1='/vol/astro7/lofar/sim/pipeline/events/'+str(int(event))+'/1/coreas/proton/steering/'
dir0='/vol/astro7/lofar/sim/pipeline/events/'+str(int(event))+'/0/coreas/proton/steering/'

dir_list=[]
if path.exists(dir3):
    dir_list.append(dir3)
if path.exists(dir2):
    dir_list.append(dir2)
if path.exists(dir1):
    dir_list.append(dir1)
if path.exists(dir0):
    dir_list.append(dir0)
    
found_run=0
cnt=0
inp_use=''
for d in np.arange(len(dir_list)):
    while found_run==0 and cnt<4:
        if str(path.exists(dir_list[cnt])):
            RUN_files=glob.glob(dir_list[cnt]+'RUN*.inp')
            with open(RUN_files[0]) as fp:
                line = fp.readline()
                while line:
                    line = fp.readline()
                    if 'ERANGE' in line:
                        e_RUN=float(line.strip().split()[1])
                    if 'SEED' in line:
                        seed_RUN=float(line.strip().split()[1])
                    if 'THETAP' in line:
                        theta_RUN=float(line.strip().split()[1])
                    if 'PHIP' in line:
                        azimuth_RUN=float(line.strip().split()[1])+270.0
        if int(e_RUN)==int(sim_energy):
            found_run=1
            inp_use=RUN_files[0]
        cnt=cnt+1




print(inp_use)

def write_file(event, azimuth, zenith, energy, seed, type):



    part_id=''
    if type=='proton':
        part_id='14'
        outfile=open(proton_dir+event+'_conex_'+type+'.q','w')
    
    if type=='helium':
        part_id='402'
        outfile=open(helium_dir+event+'_conex_'+type+'.q','w')

    if type=='oxygen':
        part_id='1608'
        outfile=open(oxygen_dir+event+'_conex_'+type+'.q','w')

    if type=='iron':
        part_id='5626'
        outfile=open(iron_dir+event+'_conex_'+type+'.q','w')





    outfile.write('#! /bin/bash\n')
    outfile.write('#SBATCH --time=1-00:00:00\n')
    outfile.write('export RUNNR=`printf "%06d" $SLURM_ARRAY_TASK_ID`\n')
    outfile.write('export FLUPRO=/vol/optcoma/cr-simulations/fluka64\n')
    outfile.write('mkdir -p {0}/events/{1}/conex/{2}/steering/\n'.format(base_dir,event,type))
    outfile.write('rm -rf /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
    outfile.write('mkdir -p /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
    outfile.write('python /vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/geninp.py --atmosphere --atmfile=/vol/astro7/lofar/sim/pipeline/atmosphere_files/ATMOSPHERE_{4}.DAT -r $RUNNR -s {0} -u {1} -a {2} -z {3} -t {5} -c True -d /scratch/kmulrey/{4}/{5}/$RUNNR/ > /scratch/kmulrey/{4}/{5}/$RUNNR/RUN$RUNNR.inp\n'.format(seed,energy,azimuth,zenith,event,part_id))

    outfile.write('cd /vol/optcoma/cr-simulations/corsika_production/run/\n')
    outfile.write('./corsika77100Linux_QGSII_fluka_thin_conex < //scratch/kmulrey/{0}/{1}/$RUNNR/RUN$RUNNR.inp\n'.format(event,part_id))
    outfile.write('cd /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
    outfile.write('mv RUN$RUNNR.inp {0}/events/{1}/conex/{2}/steering/RUN$RUNNR.inp\n'.format(base_dir,event,type))
    outfile.write('mv *.long {0}events/{1}/conex/{2}/\n'.format(base_dir,event,type))
    outfile.write('rm -rf /scratch/kmulrey/{0}/{1}/$RUNNR/*\n'.format(event,part_id))

    outfile.close()



write_file(str(int(event)),azimuth_RUN, theta_RUN, e_RUN, seed_RUN,'proton')
write_file(str(int(event)),azimuth_RUN, theta_RUN, e_RUN, seed_RUN,'helium')
write_file(str(int(event)),azimuth_RUN, theta_RUN, e_RUN, seed_RUN,'oxygen')
write_file(str(int(event)),azimuth_RUN, theta_RUN, e_RUN, seed_RUN,'iron')









