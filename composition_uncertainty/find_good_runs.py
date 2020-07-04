import numpy as np
from optparse import OptionParser
import os
import pickle
import os.path
import re
import sys
from os import listdir
from os.path import isfile, join

delta_xmax=5


parser.add_option("-e", "--event", type="string", help="event number", default="95166806")
parser.add_option("-x", "--xmax", type="int", help="target xmax", default=650)

options, arguments = parser.parse_args()
event=options.event
target_xmax=options.xmax

#event='95166806'
#target_xmax=636


def return_xmax(file):
    longfile=open(file,'r')
    hold=''
    for line in longfile:
        if "PARAMETERS" in line:
            hold=line
            break
    xmax=float(hold.split()[2:][2])
    RUNNR=file.split('/')[12].split('.')[0].split('DAT')[1]
    longfile.close()
    return RUNNR,xmax



sim_dir_proton='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/events/'+event+'/conex/proton/'
sim_dir_helium='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/events/'+event+'/conex/helium/'
sim_dir_oxygen='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/events/'+event+'/conex/oxygen/'
sim_dir_iron='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/events/'+event+'/conex/iron/'

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


print(len(proton_list_xmax))
print(len(helium_list_xmax))
print(len(oxygen_list_xmax))
print(len(iron_list_xmax))


outfile_proton=open('xmax_lists/runs_'+event+'_proton.txt','w')
outfile_helium=open('xmax_lists/runs_'+event+'_helium.txt','w')
outfile_oxygen=open('xmax_lists/runs_'+event+'_oxygen.txt','w')
outfile_iron=open('xmax_lists/runs_'+event+'_iron.txt','w')

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
