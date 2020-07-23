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
outfile.write('qsub -t ')

max_runs=15
countP=0
countHe=0
countO=0
countFe=0

for e in np.arange(len(proton_list)-1):
    if countP<(max_runs-1):
        outfile.write('{0},'.format(proton_list[countP][0]))
        countP=countP+1
outfile.write('{0}  {1}_coreas_proton.q'.format(proton_list[countP+1][0],event))

outfile.close()

 #qsub -t 1-500 $event\_conex_proton.q
