import numpy as np

proton_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/qgsjet_runs/jobs_proton_part2/'
iron_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/qgsjet_runs/jobs_iron_part2/'

base_dir='/vol/astro7/lofar/sim/pipeline_qgsjet/'

def write_file(event, type):


    
    part_id=''
    if type=='proton':
        part_id='14'
        outfile=open(proton_dir+event+'_geant_'+type+'.q','w')

    if type=='iron':
        part_id='5626'
        outfile=open(iron_dir+event+'_geant_'+type+'.q','w')





    outfile.write('#! /bin/bash\n')
    outfile.write('#SBATCH --time=1-00:00:00\n')
    outfile.write('use geant\n')
    outfile.write('export RUNNR=`printf "%06d" $SLURM_ARRAY_TASK_ID`\n')

    
    outfile.write('cd {0}/events/{1}/0/coreas/{2}/\n'.format(base_dir,event,type))
    outfile.write('rm -rf /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
    outfile.write('mkdir -p /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
    
    outfile.write('cp DAT$RUNNR /scratch/kmulrey/{0}/{1}/$RUNNR \n'.format(event,part_id))
    outfile.write('cd /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))

    outfile.write('export LOFARSOFT=/vol/optcoma/pycrtools\n')
    outfile.write('G4WORKDIR=$LOFARSOFT/LORA_simulation\n')
    outfile.write('. /vol/optcoma/geant4_9.6_install/share/Geant4-9.6.4/geant4make/geant4make.sh\n')
    
    outfile.write('/vol/optcoma/pycrtools/LORA_simulation/DAT2txt DAT$RUNNR DAT$RUNNR.tmp\n')
    outfile.write('/vol/optcoma/pycrtools/LORA_simulation/LORA_simulation DAT$RUNNR.tmp DAT$RUNNR.lora\n')
    
    outfile.write('rm DAT$RUNNR.tmp\n')
    outfile.write('cp -r * {0}events/{1}/0/coreas/{2}/\n'.format(base_dir,event,type))
    outfile.write('rm -rf /scratch/kmulrey/{0}/{1}/$RUNNR/*\n'.format(event,part_id))

    outfile.close()



event=[60409606,61271909,63246671,65490891,65214973,64960703,70904490,80495081,85083852,84432712,83467151]


for i in np.arange(len(event)):
    write_file(str(int(event[i])),'proton')
    write_file(str(int(event[i])),'iron')
