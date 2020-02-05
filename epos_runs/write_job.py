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
    outfile.write('mkdir -p {0}/events/'+event+'/coreas/'+type+'/steering/\n'.format(base_dir))
    outfile.write('rm -rf /scratch/kmulrey/'+event+'/'+part_id+'/$RUNNR\n')
    outfile.write('mkdir -p /scratch/kmulrey/'+event+'/'+part_id+'/$RUNNR\n')
    outfile.write('python /vol/optcoma/pycrtools/src/PyCRTools/extras/geninp.py --atmosphere --atmfile=/vol/astro7/lofar/sim/pipeline/atmosphere_files/ATMOSPHERE_'+event+'.DAT -r $RUNNR -s {0} -u {1} -a {2} -z {3} -t '.format(seed,energy,azimuth,zenith)+part_id+' -d /scratch/kmulrey/'+event+'/'+part_id+'/$RUNNR/ > /scratch/kmulrey/'+event+'/'+part_id+'/$RUNNR/RUN$RUNNR.inp\n')







    outfile.close()


azimuth= -149.856203773
zenith=25.1606454355
energy=166070995.766
seed=25056

write_file(str(60309606),azimuth, zenith, energy, seed,'proton')

'''

python /vol/optcoma/pycrtools/src/PyCRTools/extras/geninp.py --atmosphere --atmfile=/vol/astro7/lofar/sim/pipeline/atmosphere_files/ATMOSPHERE_188735056.DAT -r $RUNNR -s 25056 -u 166070995.766 -a -149.856203773 -z 25.1606454355 -t 5626 -d /scratch/crsim_pipeline/188735056/5626/$RUNNR/ > /scratch/crsim_pipeline/188735056/5626/$RUNNR/RUN$RUNNR.inp
cp /vol/astro7/lofar/sim/pipeline/run/SIM.reas /scratch/crsim_pipeline/188735056/5626/$RUNNR/SIM$RUNNR.reas
cp /vol/astro7/lofar/sim/pipeline/run/SIM188735056.list /scratch/crsim_pipeline/188735056/5626/$RUNNR/SIM$RUNNR.list
cd /vol/optcoma/cr-simulations/corsika_production/run
./corsika77100Linux_QGSII_fluka_thin_coreas < //scratch/crsim_pipeline/188735056/5626/$RUNNR/RUN$RUNNR.inp
cd /scratch/crsim_pipeline/188735056/5626/$RUNNR
mv RUN$RUNNR.inp /vol/astro7/lofar/sim/pipeline/events/188735056/2/coreas//iron/steering/RUN$RUNNR.inp
mv *.long /vol/astro7/lofar/sim/pipeline/events/188735056/2/coreas//iron/
# /vol/optcoma/cr-simulations/LORAtools/DAT2txt DAT$RUNNR DAT$RUNNR.tmp
# /vol/optcoma/cr-simulations/LORAtools/LORA_simulation DAT$RUNNR.tmp DAT$RUNNR.lora
# rm DAT$RUNNR.tmp
mv SIM$RUNNR.reas /vol/astro7/lofar/sim/pipeline/events/188735056/2/coreas//iron/steering/SIM$RUNNR.reas
mv SIM$RUNNR.list /vol/astro7/lofar/sim/pipeline/events/188735056/2/coreas//iron/steering/SIM$RUNNR.list
cp -r * /vol/astro7/lofar/sim/pipeline/events/188735056/2/coreas//iron/
rm -rf /scratch/crsim_pipeline/188735056/5626/$RUNNR/*
'''
