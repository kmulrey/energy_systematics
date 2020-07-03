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
    
    
    
    outfile.write('python /user/kmulrey/energy_systematics/composition_uncertainty/geninp.py --atmosphere --atmfile=/vol/astro7/lofar/sim/pipeline/atmosphere_files/ATMOSPHERE_{4}.DAT -r $RUNNR -s {0} -u {1} -a {2} -z {3} -t {5} -c True -d /scratch/kmulrey/{4}/{5}/$RUNNR/ > /scratch/kmulrey/{4}/{5}/$RUNNR/RUN$RUNNR.inp\n'.format(seed,energy,azimuth,zenith,event,part_id))

    
    outfile.write('cd /user/kmulrey/software/corsika-77100/run/\n')
    outfile.write('./corsika77100Linux_QGSII_urqmd_thin_conex < //scratch/kmulrey/{0}/{1}/$RUNNR/RUN$RUNNR.inp\n'.format(event,part_id))
    outfile.write('cd /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
    outfile.write('mv RUN$RUNNR.inp {0}/events/{1}/conex/{2}/steering/RUN$RUNNR.inp\n'.format(base_dir,event,type))
    outfile.write('mv *.long {0}events/{1}/conex/{2}/\n'.format(base_dir,event,type))
    outfile.write('rm -rf /scratch/kmulrey/{0}/{1}/$RUNNR/*\n'.format(event,part_id))

    outfile.close()



event=[70988116,86129434,120768260,122146757,148663780,48361669,156964925,175485680,174699876,158978461]

seed=[33914,12386,55760,54718,97540,50344,68200,35910,9898,21856]

energy=[767818416.752,523262814.697,275284442.453,182305582.389,132722570.742,713537600.0,202181825.107,220948238.494,196392329.376,542183573.418]

zenith=[41.08,30.82,29.99,28.97,23.61,37.47,22.15,24.73,28.89,34.85]

azimuth=[-308.45+270,-142.30+270,-100.32+270,-145.26+270,-103.83+270,-138.05+270,-156.61+270,-237.53+270,-434.78+270,-208.25+270]



#azimuth= -149.856203773
#zenith=25.1606454355
#energy=166070995.766
#seed=25056
for i in np.arange(len(event)):
    write_file(str(int(event[i])),azimuth[i], zenith[i], energy[i], seed[i],'proton')
    write_file(str(int(event[i])),azimuth[i], zenith[i], energy[i], seed[i],'helium')
    write_file(str(int(event[i])),azimuth[i], zenith[i], energy[i], seed[i],'oxygen')
    write_file(str(int(event[i])),azimuth[i], zenith[i], energy[i], seed[i],'iron')
