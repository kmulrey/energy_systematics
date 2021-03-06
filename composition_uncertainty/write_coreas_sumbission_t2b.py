import numpy as np
from optparse import OptionParser



parser = OptionParser()
parser.add_option("-e", "--event", type="string", help="event number", default="95166806")

options, arguments = parser.parse_args()
event=options.event

proton_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/jobs_proton/'
helium_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/jobs_helium/'
oxygen_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/jobs_oxygen/'
iron_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/jobs_iron/'

base_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/'



proton_list=np.genfromtxt('xmax_lists/runs_'+event+'_proton.txt')
helium_list=np.genfromtxt('xmax_lists/runs_'+event+'_helium.txt')
oxygen_list=np.genfromtxt('xmax_lists/runs_'+event+'_oxygen.txt')
iron_list=np.genfromtxt('xmax_lists/runs_'+event+'_iron.txt')
 
outfile=open('submit_coreas_{0}.sh'.format(event),'w')

outfile.write('#! /bin/bash\n')

max_runs=15
countP=0
countHe=0
countO=0
countFe=0

outfile.write('cd /user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/jobs_proton\n')
outfile.write('qsub -t ')
for e in np.arange(len(proton_list)-1):
    if countP<(max_runs-1) and countP<(len(proton_list)-1):
        outfile.write('{0},'.format(int(proton_list[countP][0])))
        countP=countP+1
outfile.write('{0}  {1}_coreas_proton.q\n'.format(int(proton_list[countP][0]),event))


outfile.write('cd /user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/jobs_helium\n')
outfile.write('qsub -t ')
for e in np.arange(len(helium_list)-1):
    if countHe<(max_runs-1) and countHe<(len(helium_list)-1):
        outfile.write('{0},'.format(int(helium_list[countHe][0])))
        countHe=countHe+1
outfile.write('{0}  {1}_coreas_helium.q\n'.format(int(helium_list[countHe][0]),event))


outfile.write('cd /user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/jobs_oxygen\n')
outfile.write('qsub -t ')
for e in np.arange(len(oxygen_list)-1):
    if countO<(max_runs-1) and countO<(len(oxygen_list)-1):
        outfile.write('{0},'.format(int(oxygen_list[countO][0])))
        countO=countO+1
outfile.write('{0}  {1}_coreas_oxygen.q\n'.format(int(oxygen_list[countO][0]),event))



outfile.write('cd /user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/jobs_iron\n')
outfile.write('qsub -t ')
for e in np.arange(len(iron_list)-1):
    if countFe<(max_runs-1) and  countFe<(len(iron_list)-1):
        outfile.write('{0},'.format(int(iron_list[countFe][0])))
        countFe=countFe+1
outfile.write('{0}  {1}_coreas_iron.q\n'.format(int(iron_list[countFe][0]),event))








outfile.close()

 #qsub -t 1-500 $event\_conex_proton.q
