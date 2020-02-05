proton_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/epos_runs/jobs_proton/'
iron_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/epos_runs/jobs_iron/'

base_dir='/vol/astro7/lofar/sim/pipeline_epos/'

def write_file(event, azimuth, zenith, energy, seed,type):


    
    part_id=''
    if type=='proton':
        part_id='14'
    if type=='iron':
        part_id='5626'




    outfile=open(proton_dir+event+'_coreas_'+type+'.q','w')

    outfile.write('#! /bin/bash\n')
    outfile.write('#SBATCH --time=7-00:00:00\n')
    outfile.write('#SBATCH --output {0}/run/output/'+event+'_coreas_'+part_id+'-%j\n'.format(base_dir))
    outfile.write('#SBATCH --error {0}/run/output/'+event+'_coreas_'+part_id+'-ERROR-%j\n'.format(base_dir))



    outfile.write('. /vol/optcoma/geant4_9.6_install/share/Geant4-9.6.4/geant4make/geant4make.sh\n')
    outfile.write('export RUNNR=`printf "%06d" $SLURM_ARRAY_TASK_ID`\n')
    outfile.write('export FLUPRO=/vol/optcoma/cr-simulations/fluka64\n')
    outfile.write('cd {0}/run/'.format(base_dir))
    outfile.write('mkdir -p {0}/events/{1}/coreas/{2}/steering/\n'.format(base_dir,event,type))
    outfile.write('rm -rf /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
    outfile.write('mkdir -p /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
    outfile.write('python /vol/optcoma/pycrtools/src/PyCRTools/extras/geninp.py --atmosphere --atmfile=/vol/astro7/lofar/sim/pipeline/atmosphere_files/ATMOSPHERE_{4}.DAT -r $RUNNR -s {0} -u {1} -a {2} -z {3} -t {5} -d /scratch/kmulrey/{4}/{5}/$RUNNR/ > /scratch/kmulrey/{4}/{5}/$RUNNR/RUN$RUNNR.inp\n'.format(seed,energy,azimuth,zenith,event,part_id))

    outfile.write('cp /vol/astro7/lofar/sim/pipeline/run/SIM.reas /scratch/kmulrey/{0}/{1}/$RUNNR/SIM$RUNNR.reas\n'.format(event,part_id))
    outfile.write('cp /vol/astro7/lofar/sim/pipeline/run/SIM{0}.list /scratch/kmulrey/{0}/{1}/$RUNNR/SIM$RUNNR.list\n'.format(event,part_id))
    outfile.write('cd /vol/optcoma/cr-simulations/corsika_production_epos/run\n')
    outfile.write('./corsika77100Linux_EPOS_urqmd_thin_conex_coreas < //scratch/kmulrey/{0}/{1}/$RUNNR/RUN$RUNNR.inp\n'.format(event,part_id))
    outfile.write('cd /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
    outfile.write('mv RUN$RUNNR.inp {0}/events/{1}/coreas/{2}/steering/RUN$RUNNR.inp\n'.format(base_dir,event,type))
    outfile.write('mv *.long {0}events/{1}/coreas/{2}/\n'.format(base_dir,event,type))
    outfile.write('/vol/optcoma/pycrtools/LORA_simulation/DAT2txt DAT$RUNNR DAT$RUNNR.tmp\n')
    outfile.write('/vol/optcoma/pycrtools/LORA_simulation/LORA_simulation DAT$RUNNR.tmp DAT$RUNNR.lora\n')
    outfile.write('rm DAT$RUNNR.tmp\n')
    outfile.write('mv SIM$RUNNR.reas {0}events/{1}/coreas/{2}/steering/SIM$RUNNR.reas\n'.format(base_dir,event,type))
    outfile.write('mv SIM$RUNNR.list {0}events/{1}/coreas/{2}/steering/SIM$RUNNR.list\n'.format(base_dir,event,type))
    outfile.write('cp -r * {0}events/{1}/coreas/{2}/\n'.format(base_dir,event,type))
    outfile.write('rm -rf /scratch/kmulrey/{0}/{1}/$RUNNR/*\n'.format(event,part_id))

    outfile.close()
            


azimuth= -149.856203773
zenith=25.1606454355
energy=166070995.766
seed=25056

write_file(str(60309606),azimuth, zenith, energy, seed,'proton')
