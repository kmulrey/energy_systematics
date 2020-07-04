import numpy as np
from optparse import OptionParser


parser = OptionParser()
parser.add_option("-e", "--event", type="string", help="event number", default="95166806")

options, arguments = parser.parse_args()
event=options.event


proton_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/events/'+event+'/conex/proton/steering/'
helium_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/events/'+event+'/conex/helium/steering/'
oxygen_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/events/'+event+'/conex/oxygen/steering/'
iron_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/events/'+event+'/conex/iron/steering/'

proton_inp_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/jobs_proton/'+event+'/'
helium_inp_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/jobs_helium/'+event+'/'
oxygen_inp_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/jobs_oxygen/'+event+'/'
iron_inp_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/jobs_iron/'+event+'/'


proton_list=np.genfromtxt('xmax_lists/runs_'+event+'_proton.txt')
helium_list=np.genfromtxt('xmax_lists/runs_'+event+'_helium.txt')
oxygen_list=np.genfromtxt('xmax_lists/runs_'+event+'_oxygen.txt')
iron_list=np.genfromtxt('xmax_lists/runs_'+event+'_iron.txt')


for e in np.arange(len(proton_list)):
    RUNNR=str(int(proton_list[e][0])).zfill(6)
    proton_filepath=proton_dir+'RUN'+RUNNR+'.inp'
    outputfile=open(proton_inp_dir+'RUN'+RUNNR+'.inp','w')
    with open(proton_filepath) as fp:
        line = fp.readline()
        cnt = 1
        while line:
            if 'PAROUT' in line:
                outputfile.write('PAROUT  T F\n')
            elif 'CASCADE' in line:
                #outputfile.write('CASCADE F F F\n')
                print("0")
            else:
                outputfile.write(line)
     
            line = fp.readline()
            cnt += 1

    outputfile.close()
    
for e in np.arange(len(helium_list)):
    RUNNR=str(int(helium_list[e][0])).zfill(6)
    helium_filepath=helium_dir+'RUN'+RUNNR+'.inp'
    outputfile=open(helium_inp_dir+'RUN'+RUNNR+'.inp','w')
    with open(helium_filepath) as fp:
        line = fp.readline()
        cnt = 1
        while line:
            if 'PAROUT' in line:
                outputfile.write('PAROUT  T F\n')
            elif 'CASCADE' in line:
                #outputfile.write('CASCADE F F F\n')
                print("0")
            else:
                outputfile.write(line)
     
            line = fp.readline()
            cnt += 1

    outputfile.close()


for e in np.arange(len(oxygen_list)):
    RUNNR=str(int(oxygen_list[e][0])).zfill(6)
    oxygen_filepath=oxygen_dir+'RUN'+RUNNR+'.inp'
    outputfile=open(oxygen_inp_dir+'RUN'+RUNNR+'.inp','w')
    with open(oxygen_filepath) as fp:
        line = fp.readline()
        cnt = 1
        while line:
            if 'PAROUT' in line:
                outputfile.write('PAROUT  T F\n')
            elif 'CASCADE' in line:
                #outputfile.write('CASCADE F F F\n')
                print("0")
            else:
                outputfile.write(line)
 
            line = fp.readline()
            cnt += 1

    outputfile.close()
    
for e in np.arange(len(iron_list)):
    RUNNR=str(int(iron_list[e][0])).zfill(6)
    iron_filepath=iron_dir+'RUN'+RUNNR+'.inp'
    outputfile=open(iron_inp_dir+'RUN'+RUNNR+'.inp','w')
    with open(iron_filepath) as fp:
        line = fp.readline()
        cnt = 1
        while line:
            if 'PAROUT' in line:
                outputfile.write('PAROUT  T F\n')
            elif 'CASCADE' in line:
                #outputfile.write('CASCADE F F F\n')
                print("0")
            else:
                outputfile.write(line)
     
            line = fp.readline()
            cnt += 1

    outputfile.close()
