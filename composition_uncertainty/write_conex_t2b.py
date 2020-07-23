import numpy as np

proton_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/jobs_proton/'
helium_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/jobs_helium/'
oxygen_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/jobs_oxygen/'
iron_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/jobs_iron/'

base_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/'
run_dir='/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/run_files/'

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
    outfile.write('export RUNNR=`printf "%06d" $PBS_ARRAYID`\n')

    outfile.write('mkdir -p {0}/events/{1}/conex/{2}/steering/\n'.format(base_dir,event,type))
    outfile.write('rm -rf /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
    outfile.write('mkdir -p /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
    
    
    
    outfile.write('python /user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/geninp.py --atmosphere --atmfile=/user/kmulrey/energy_systematics/energy_systematics/composition_uncertainty/run_files/ATMOSPHERE_{4}.DAT -r $RUNNR -s {0} -u {1} -a {2} -z {3} -t {5} -c True -d /scratch/kmulrey/{4}/{5}/$RUNNR/ > /scratch/kmulrey/{4}/{5}/$RUNNR/RUN$RUNNR.inp\n'.format(seed,energy,azimuth,zenith,event,part_id))

    outfile.write('cp {2}/SIM.reas /scratch/kmulrey/{0}/{1}/$RUNNR/SIM$RUNNR.reas\n'.format(event,part_id,run_dir))
    outfile.write('cp {2}/SIM{0}.list /scratch/kmulrey/{0}/{1}/$RUNNR/SIM$RUNNR.list\n'.format(event,part_id,run_dir))
    outfile.write('cd /user/kmulrey/software/corsika-77100/run/\n')
    outfile.write('./corsika77100Linux_QGSII_urqmd_thin_conex_coreas < //scratch/kmulrey/{0}/{1}/$RUNNR/RUN$RUNNR.inp\n'.format(event,part_id))
    outfile.write('cd /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
    outfile.write('mv RUN$RUNNR.inp {0}/events/{1}/conex/{2}/steering/RUN$RUNNR.inp\n'.format(base_dir,event,type))
    outfile.write('mv *.long {0}events/{1}/conex/{2}/\n'.format(base_dir,event,type))
    outfile.write('rm -rf /scratch/kmulrey/{0}/{1}/$RUNNR/*\n'.format(event,part_id))

    outfile.close()


event=[216660305,192006137,157888051,81147431,48480819,132649890,102912996,196164498,110169343]
seed=[44370,97808,37864,64284,31506,98560,82054,12052,70862]
energy=[112738808.826,306519397.535,159167985.034,1185809599.6,63264163.3962,434581545.376,647051721.228,304057589.037,166614702.194]
zenith=[26.52,24.06,26.8,54.73,23.43,39.73,30.92,25.75,13.63]
azimuth=[-98.38+270,-202.67+270,-138.64+270,-367.29+270,-208.10+270,-293.59+270,-102.75+270,-97.73+270,-189.01+270]



#azimuth= -149.856203773
#zenith=25.1606454355
#energy=166070995.766
#seed=25056
for i in np.arange(len(event)):
    write_file(str(int(event[i])),azimuth[i], zenith[i], energy[i], seed[i],'proton')
    write_file(str(int(event[i])),azimuth[i], zenith[i], energy[i], seed[i],'helium')
    write_file(str(int(event[i])),azimuth[i], zenith[i], energy[i], seed[i],'oxygen')
    write_file(str(int(event[i])),azimuth[i], zenith[i], energy[i], seed[i],'iron')
