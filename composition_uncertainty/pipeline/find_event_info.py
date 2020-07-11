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
        if int(e_RUN)==int(sim_energy):
            found_run=1
            inp_use=RUN_files[0]
        cnt=cnt+1




print(inp_use)
