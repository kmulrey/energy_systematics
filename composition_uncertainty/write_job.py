import numpy as np

proton_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/jobs_proton/'
helium_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/jobs_helium/'
oxygen_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/jobs_oxygen/'
iron_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/jobs_iron/'

base_dir='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/'

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
    #outfile.write('#SBATCH --time=1-00:00:00\n')
    #outfile.write('#SBATCH --output {0}/run/output/{1}_coreas_{2}-%j\n'.format(base_dir,event,part_id))
    #outfile.write('#SBATCH --error {0}/run/output/{1}_coreas_{2}-ERROR-%j\n'.format(base_dir,event,part_id))

    #outfile.write('umask 002\n')
    #outfile.write('use geant\n')
    #outfile.write('cd /vol/optcoma/geant4_9.6_install/share/Geant4-9.6.4/geant4make/\n')
    #outfile.write('./geant4make.sh\n')
    #outfile.write('cd {0}/events/\n'.format(base_dir))

    outfile.write('export RUNNR=`printf "%06d" $SLURM_ARRAY_TASK_ID`\n')
    #outfile.write('export FLUPRO=/vol/optcoma/cr-simulations/fluka64\n')
    outfile.write('mkdir -p {0}/events/{1}/conex/{2}/steering/\n'.format(base_dir,event,type))
    outfile.write('rm -rf /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
    outfile.write('mkdir -p /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
    outfile.write('python /vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/composition_uncertainty/geninp.py --atmosphere --atmfile=/vol/astro7/lofar/sim/pipeline/atmosphere_files/ATMOSPHERE_{4}.DAT -r $RUNNR -s {0} -u {1} -a {2} -z {3} -t {5} -c True -d /scratch/kmulrey/{4}/{5}/$RUNNR/ > /scratch/kmulrey/{4}/{5}/$RUNNR/RUN$RUNNR.inp\n'.format(seed,energy,azimuth,zenith,event,part_id))

    #outfile.write('cp /vol/astro7/lofar/sim/pipeline/run/SIM.reas /scratch/kmulrey/{0}/{1}/$RUNNR/SIM$RUNNR.reas\n'.format(event,part_id))
    #outfile.write('cp /vol/astro7/lofar/sim/pipeline/run/SIM{0}.list /scratch/kmulrey/{0}/{1}/$RUNNR/SIM$RUNNR.list\n'.format(event,part_id))
    outfile.write('cd /vol/optcoma/cr-simulations/corsika_production/run/\n')
    outfile.write('./corsika77100Linux_QGSII_fluka_thin_conex < //scratch/kmulrey/{0}/{1}/$RUNNR/RUN$RUNNR.inp\n'.format(event,part_id))
    outfile.write('cd /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
    outfile.write('mv RUN$RUNNR.inp {0}/events/{1}/conex/{2}/steering/RUN$RUNNR.inp\n'.format(base_dir,event,type))
    outfile.write('mv *.long {0}events/{1}/conex/{2}/\n'.format(base_dir,event,type))
    #outfile.write('/vol/optcoma/pycrtools/LORA_simulation/DAT2txt DAT$RUNNR DAT$RUNNR.tmp\n')
    #outfile.write('/vol/optcoma/pycrtools/LORA_simulation/LORA_simulation DAT$RUNNR.tmp DAT$RUNNR.lora\n')
    #outfile.write('rm DAT$RUNNR.tmp\n')
    #outfile.write('mv SIM$RUNNR.reas {0}events/{1}/coreas/{2}/steering/SIM$RUNNR.reas\n'.format(base_dir,event,type))
    #outfile.write('mv SIM$RUNNR.list {0}events/{1}/coreas/{2}/steering/SIM$RUNNR.list\n'.format(base_dir,event,type))
    #outfile.write('cp -r * {0}events/{1}/conex/{2}/\n'.format(base_dir,event,type))
    outfile.write('rm -rf /scratch/kmulrey/{0}/{1}/$RUNNR/*\n'.format(event,part_id))

    outfile.close()



event=[61254596,95166806,126484310,177295087,]
azimuth=[-434.547641748+270,-191.888958177+270,-131.364106846+270,-411.974034426+270]
zenith=[21.659887910,35.9479424193,36.4747650732,22.475475965]
energy=[91947841.7402,197551979.48,785829352.566,214482624.249]
seed=[64704,35764,28210,59522]



#azimuth= -149.856203773
#zenith=25.1606454355
#energy=166070995.766
#seed=25056
for i in np.arange(len(event)):
    write_file(str(int(event[i])),azimuth[i], zenith[i], energy[i], seed[i],'proton')
    write_file(str(int(event[i])),azimuth[i], zenith[i], energy[i], seed[i],'helium')
    write_file(str(int(event[i])),azimuth[i], zenith[i], energy[i], seed[i],'oxygen')
    write_file(str(int(event[i])),azimuth[i], zenith[i], energy[i], seed[i],'iron')
