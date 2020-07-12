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

delta_xmax=5

def return_xmax(file):
    longfile=open(file,'r')
    hold=''
    #print(longfile)
    for line in longfile:
        if "PARAMETERS" in line:
            hold=line
            break
    xmax=float(hold.split()[2:][2])
    #print(file.split('/')[11].split('.')[0].split('DAT')[1])
    RUNNR=file.split('/')[11].split('.')[0].split('DAT')[1]
    longfile.close()
    return RUNNR,xmax


reco_dir='/vol/astro7/lofar/sim/pipeline/production_analysis_Dec2019/'
proton_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/jobs_proton/'
helium_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/jobs_helium/'
oxygen_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/jobs_oxygen/'
iron_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/jobs_iron/'

base_dir='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/'


for file in glob.glob(reco_dir+'*'+str(int(event))+'*dat'):
    reco_file_name=file
print(reco_file_name)

reco_file=open(reco_file_name,"rb")
reco_info=pickle.load(reco_file)
reco_file.close()
print(reco_info.keys())

core_x=reco_info['core_x']+reco_info['xoff']
core_y=reco_info['core_y']+reco_info['yoff']
target_xmax=reco_info['xmaxreco']

sim_dir_proton='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/'+str(int(event))+'/conex/proton/'
sim_dir_helium='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/'+str(int(event))+'/conex/helium/'
sim_dir_oxygen='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/'+str(int(event))+'/conex/oxygen/'
sim_dir_iron='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/'+str(int(event))+'/conex/iron/'

proton_inp_dir='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/'+str(int(event))+'/corsika/proton/steering/'
helium_inp_dir='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/'+str(int(event))+'/corsika/helium/steering/'
oxygen_inp_dir='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/'+str(int(event))+'/corsika/oxygen/steering/'
iron_inp_dir='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/'+str(int(event))+'/corsika/iron/steering/'



longfiles_proton = [f for f in listdir(sim_dir_proton) if isfile(join(sim_dir_proton, f))]
longfiles_helium = [f for f in listdir(sim_dir_helium) if isfile(join(sim_dir_helium, f))]
longfiles_oxygen = [f for f in listdir(sim_dir_oxygen) if isfile(join(sim_dir_oxygen, f))]
longfiles_iron = [f for f in listdir(sim_dir_iron) if isfile(join(sim_dir_iron, f))]

proton_list_RUNNR=[]
proton_list_xmax=[]

helium_list_RUNNR=[]
helium_list_xmax=[]

oxygen_list_RUNNR=[]
oxygen_list_xmax=[]

iron_list_RUNNR=[]
iron_list_xmax=[]

for i in np.arange(len(longfiles_proton)):
    RUNNR,xmax=return_xmax(sim_dir_proton+longfiles_proton[i])
    if np.abs(xmax-target_xmax)<delta_xmax:
        proton_list_RUNNR.append(RUNNR)
        proton_list_xmax.append(xmax)
        
for i in np.arange(len(longfiles_helium)):
    RUNNR,xmax=return_xmax(sim_dir_helium+longfiles_helium[i])
    if np.abs(xmax-target_xmax)<delta_xmax:
        helium_list_RUNNR.append(RUNNR)
        helium_list_xmax.append(xmax)
        
for i in np.arange(len(longfiles_oxygen)):
    RUNNR,xmax=return_xmax(sim_dir_oxygen+longfiles_oxygen[i])
    if np.abs(xmax-target_xmax)<delta_xmax:
        oxygen_list_RUNNR.append(RUNNR)
        oxygen_list_xmax.append(xmax)
        
for i in np.arange(len(longfiles_iron)):
    RUNNR,xmax=return_xmax(sim_dir_iron+longfiles_iron[i])
    if np.abs(xmax-target_xmax)<delta_xmax:
        iron_list_RUNNR.append(RUNNR)
        iron_list_xmax.append(xmax)
        
outfile_proton=open('xmax_lists/runs_'+str(int(event))+'_proton.txt','w')
outfile_helium=open('xmax_lists/runs_'+str(int(event))+'_helium.txt','w')
outfile_oxygen=open('xmax_lists/runs_'+str(int(event))+'_oxygen.txt','w')
outfile_iron=open('xmax_lists/runs_'+str(int(event))+'_iron.txt','w')

for i in np.arange(len(proton_list_xmax)):
    outfile_proton.write('{0}   {1}\n'.format(proton_list_RUNNR[i], proton_list_xmax[i]))
for i in np.arange(len(helium_list_xmax)):
    outfile_helium.write('{0}   {1}\n'.format(helium_list_RUNNR[i], helium_list_xmax[i]))
for i in np.arange(len(oxygen_list_xmax)):
    outfile_oxygen.write('{0}   {1}\n'.format(oxygen_list_RUNNR[i], oxygen_list_xmax[i]))
for i in np.arange(len(iron_list_xmax)):
    outfile_iron.write('{0}   {1}\n'.format(iron_list_RUNNR[i], iron_list_xmax[i]))
    
outfile_proton.close()
outfile_helium.close()
outfile_oxygen.close()
outfile_iron.close()


# write proton files

for e in np.arange(len(proton_list_RUNNR)):
    RUNNR=str(int(proton_list_RUNNR[e])).zfill(6)
    proton_filepath=sim_dir_proton+'steering/'+'RUN'+RUNNR+'.inp'
    outputfile=open(proton_inp_dir+'RUN'+RUNNR+'.inp','w')
    with open(proton_filepath) as fp:
        line = fp.readline()
        cnt = 1
        while line:
            if 'PAROUT' in line:
                outputfile.write('PAROUT  T F\n')
            elif 'CASCADE' in line:
                outputfile.write('CASCADE F F F\n')
            else:
                outputfile.write(line)
     
            line = fp.readline()
            cnt += 1

    outputfile.close()
    
# write helium files

for e in np.arange(len(helium_list_RUNNR)):
    RUNNR=str(int(helium_list_RUNNR[e])).zfill(6)
    helium_filepath=sim_dir_helium+'steering/'+'RUN'+RUNNR+'.inp'
    outputfile=open(helium_inp_dir+'RUN'+RUNNR+'.inp','w')
    with open(helium_filepath) as fp:
        line = fp.readline()
        cnt = 1
        while line:
            if 'PAROUT' in line:
                outputfile.write('PAROUT  T F\n')
            elif 'CASCADE' in line:
                outputfile.write('CASCADE F F F\n')
            else:
                outputfile.write(line)
     
            line = fp.readline()
            cnt += 1

    outputfile.close()
    
# write oxygen files

for e in np.arange(len(oxygen_list_RUNNR)):
    RUNNR=str(int(oxygen_list_RUNNR[e])).zfill(6)
    oxygen_filepath=sim_dir_oxygen+'steering/'+'RUN'+RUNNR+'.inp'
    outputfile=open(oxygen_inp_dir+'RUN'+RUNNR+'.inp','w')
    with open(oxygen_filepath) as fp:
        line = fp.readline()
        cnt = 1
        while line:
            if 'PAROUT' in line:
                outputfile.write('PAROUT  T F\n')
            elif 'CASCADE' in line:
                outputfile.write('CASCADE F F F\n')
            else:
                outputfile.write(line)
     
            line = fp.readline()
            cnt += 1

    outputfile.close()
    
# write iron files

for e in np.arange(len(iron_list_RUNNR)):
    RUNNR=str(int(iron_list_RUNNR[e])).zfill(6)
    iron_filepath=sim_dir_iron+'steering/'+'RUN'+RUNNR+'.inp'
    outputfile=open(iron_inp_dir+'RUN'+RUNNR+'.inp','w')
    with open(iron_filepath) as fp:
        line = fp.readline()
        cnt = 1
        while line:
            if 'PAROUT' in line:
                outputfile.write('PAROUT  T F\n')
            elif 'CASCADE' in line:
                outputfile.write('CASCADE F F F\n')
            else:
                outputfile.write(line)
     
            line = fp.readline()
            cnt += 1

    outputfile.close()
    
    
def write_file(event, type):


    part_id=''
    if type=='proton':
        part_id='14'
        outfile=open(proton_dir+event+'_corsika_'+type+'.q','w')
        inp_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/jobs_proton/'
     
    if type=='helium':
        part_id='402'
        outfile=open(helium_dir+event+'_corsika_'+type+'.q','w')
        inp_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/jobs_helium/'

    if type=='oxygen':
        part_id='1608'
        outfile=open(oxygen_dir+event+'_corsika_'+type+'.q','w')
        inp_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/jobs_oxygen/'

    if type=='iron':
        part_id='5626'
        outfile=open(iron_dir+event+'_corsika_'+type+'.q','w')
        inp_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/jobs_iron/'

    outfile.write('#! /bin/bash\n')
    outfile.write('#SBATCH --time=2-00:00:00\n')
 
    outfile.write('umask 002\n')
    outfile.write('use geant\n')
    outfile.write('export LOFARSOFT=/vol/optcoma/pycrtools\n')
    outfile.write('G4WORKDIR=$LOFARSOFT/LORA_simulation\n')


    outfile.write('cd /vol/optcoma/geant4_9.6_install/share/Geant4-9.6.4/geant4make/\n')
    outfile.write('./geant4make.sh\n')
    outfile.write('cd {0}/events/\n'.format(base_dir))

    outfile.write('export RUNNR=`printf "%06d" $SLURM_ARRAY_TASK_ID`\n')
    outfile.write('export FLUPRO=/vol/optcoma/cr-simulations/fluka64\n')
    outfile.write('mkdir -p {0}/events/{1}/corsika/{2}/steering/\n'.format(base_dir,event,type))
    outfile.write('rm -rf /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
    outfile.write('mkdir -p /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))

    outfile.write('cp {2}/RUN$RUNNR.inp /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id,proton_inp_dir))
    outfile.write('cd /vol/optcoma/cr-simulations/corsika_production/run/\n')
    outfile.write('./corsika77100Linux_QGSII_fluka_thin_conex < //scratch/kmulrey/{0}/{1}/$RUNNR/RUN$RUNNR.inp\n'.format(event,part_id))
    outfile.write('cd /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
    outfile.write('mv RUN$RUNNR.inp {0}/events/{1}/corsika/{2}/steering/RUN$RUNNR.inp\n'.format(base_dir,event,type))
    outfile.write('mv *.long {0}events/{1}/corsika/{2}/\n'.format(base_dir,event,type))
    outfile.write('/vol/optcoma/pycrtools/LORA_simulation/DAT2txt DAT$RUNNR DAT$RUNNR.tmp\n')
    outfile.write('/vol/optcoma/pycrtools/LORA_simulation/LORA_simulation DAT$RUNNR.tmp DAT$RUNNR.lora\n')
    outfile.write('rm DAT$RUNNR.tmp\n')
    #outfile.write('mv SIM$RUNNR.reas {0}events/{1}/coreas/{2}/steering/SIM$RUNNR.reas\n'.format(base_dir,event,type))
    #outfile.write('mv SIM$RUNNR.list {0}events/{1}/coreas/{2}/steering/SIM$RUNNR.list\n'.format(base_dir,event,type))
    outfile.write('cp -r * {0}events/{1}/corsika/{2}/\n'.format(base_dir,event,type))
    outfile.write('rm -rf /scratch/kmulrey/{0}/{1}/$RUNNR/*\n'.format(event,part_id))

    outfile.close()


write_file(str(int(event)),'proton')
write_file(str(int(event)),'helium')
write_file(str(int(event)),'oxygen')
write_file(str(int(event)),'iron')
