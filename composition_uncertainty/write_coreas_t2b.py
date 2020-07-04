import numpy as np

proton_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/jobs_proton/'
helium_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/jobs_helium/'
oxygen_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/jobs_oxygen/'
iron_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/jobs_iron/'

base_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/'
run_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/run_files/'

def write_file(event, type):


    part_id=''
    if type=='proton':
        part_id='14'
        outfile=open(proton_dir+event+'_coreas_'+type+'.q','w')
        
    if type=='helium':
        part_id='402'
        outfile=open(helium_dir+event+'_coreas_'+type+'.q','w')

    if type=='oxygen':
        part_id='1608'
        outfile=open(oxygen_dir+event+'_coreas_'+type+'.q','w')

    if type=='iron':
        part_id='5626'
        outfile=open(iron_dir+event+'_coreas_'+type+'.q','w')




    inp_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/jobs_'+type+'/'+event+'/'

    outfile.write('#! /bin/bash\n')

    outfile.write('export RUNNR=`printf "%06d" $PBS_ARRAYID`\n')
    outfile.write('mkdir -p {0}/events/{1}/coreas/{2}/steering/\n'.format(base_dir,event,type))
    outfile.write('rm -rf /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
    outfile.write('mkdir -p /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
   
    outfile.write('cp {2}/RUN$RUNNR.inp /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id,inp_dir))
    outfile.write('cp {2}/SIM.reas /scratch/kmulrey/{0}/{1}/$RUNNR/SIM$RUNNR.reas\n'.format(event,part_id,run_dir))
    outfile.write('cp {2}/SIM{0}.list /scratch/kmulrey/{0}/{1}/$RUNNR/SIM$RUNNR.reas\n'.format(event,part_id,run_dir))   
   
    outfile.write('cd /user/kmulrey/software/corsika-77100/run/\n')
    outfile.write('./corsika77100Linux_QGSII_urqmd_thin_coreas < //scratch/kmulrey/{0}/{1}/$RUNNR/RUN$RUNNR.inp\n'.format(event,part_id))
    outfile.write('cd /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
    outfile.write('mv RUN$RUNNR.inp {0}/events/{1}/corsika/{2}/steering/RUN$RUNNR.inp\n'.format(base_dir,event,type))
    outfile.write('mv *.long {0}events/{1}/corsika/{2}/\n'.format(base_dir,event,type))
    
    outfile.write('mv SIM$RUNNR.reas {0}events/{1}/coreas/{2}/steering/SIM$RUNNR.reas\n'.format(base_dir,event,type))
    outfile.write('mv SIM$RUNNR.list {0}events/{1}/coreas/{2}/steering/SIM$RUNNR.list\n'.format(base_dir,event,type))
    outfile.write('cp -r * {0}events/{1}/corsika/{2}/\n'.format(base_dir,event,type))
    outfile.write('rm -rf /scratch/kmulrey/{0}/{1}/$RUNNR/*\n'.format(event,part_id))

    outfile.close()
    
    
   
   



#event=[61254596,95166806,126484310,177295087,]
#azimuth=[-434.547641748+270,-191.888958177+270,-131.364106846+270,-411.974034426+270]
#zenith=[21.659887910,35.9479424193,36.4747650732,22.475475965]
#energy=[91947841.7402,197551979.48,785829352.566,214482624.249]
#seed=[64704,35764,28210,59522]

event=[95166806,126484310]
for i in np.arange(len(event)):
    write_file(str(int(event[i])),'proton')
    write_file(str(int(event[i])),'helium')
    write_file(str(int(event[i])),'oxygen')
    write_file(str(int(event[i])),'iron')

