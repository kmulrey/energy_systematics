import numpy as np

proton_dir='/user/kmulrey/energy_systematics/energy_systematics/sibyll_runs/jobs_proton/'
iron_dir='/user/kmulrey/energy_systematics/energy_systematics/sibyll_runs/jobs_iron/'

base_dir='/user/kmulrey/energy_systematics/pipeline_sibyll/'

run_dir='/user/kmulrey/energy_systematics/energy_systematics/run_files/'

def write_file(event, azimuth, zenith, energy, seed, type):


    
    part_id=''
    if type=='proton':
        part_id='14'
        outfile=open(proton_dir+event+'_coreas_'+type+'.sh','w')

    if type=='iron':
        part_id='5626'
        outfile=open(iron_dir+event+'_coreas_'+type+'.sh','w')





    outfile.write('#! /bin/bash\n')
    #outfile.write('#SBATCH --time=7-00:00:00\n')
    #outfile.write('#SBATCH --output {0}/run/output/{1}_coreas_{2}-%j\n'.format(base_dir,event,part_id))
    #outfile.write('#SBATCH --error {0}/run/output/{1}_coreas_{2}-ERROR-%j\n'.format(base_dir,event,part_id))

    #outfile.write('umask 002\n')
    #outfile.write('use geant\n')
    #outfile.write('cd /vol/optcoma/geant4_9.6_install/share/Geant4-9.6.4/geant4make/\n')
    #outfile.write('./geant4make.sh\n')
    #outfile.write('cd {0}/run/\n'.format(base_dir))

    outfile.write('export RUNNR=`printf "%06d" $PBS_ARRAYID`\n')
    #outfile.write('export FLUPRO=/vol/optcoma/cr-simulations/fluka64\n')
    outfile.write('mkdir -p {0}/events/{1}/coreas/{2}/steering/\n'.format(base_dir,event,type))
    outfile.write('rm -rf /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
    outfile.write('mkdir -p /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
    outfile.write('python /user/kmulrey/energy_systematics/energy_systematics/sibyll_runs/geninp.py --atmosphere --atmfile={6}/ATMOSPHERE_{4}.DAT -r $RUNNR -s {0} -u {1} -a {2} -z {3} -t {5} -d /scratch/kmulrey/{4}/{5}/$RUNNR/ > /scratch/kmulrey/{4}/{5}/$RUNNR/RUN$RUNNR.inp\n'.format(seed,energy,azimuth,zenith,event,part_id,run_dir))

    outfile.write('cp {2}/SIM.reas /scratch/kmulrey/{0}/{1}/$RUNNR/SIM$RUNNR.reas\n'.format(event,part_id,run_dir))
    outfile.write('cp {2}/SIM{0}.list /scratch/kmulrey/{0}/{1}/$RUNNR/SIM$RUNNR.list\n'.format(event,part_id,run_dir))
    outfile.write('cd /user/kmulrey/software/corsika-77100/run/\n')
    outfile.write('./corsika77100Linux_SIBYLL_urqmd_thin_coreas < //scratch/kmulrey/{0}/{1}/$RUNNR/RUN$RUNNR.inp\n'.format(event,part_id))
    outfile.write('cd /scratch/kmulrey/{0}/{1}/$RUNNR\n'.format(event,part_id))
    outfile.write('mv RUN$RUNNR.inp {0}/events/{1}/coreas/{2}/steering/RUN$RUNNR.inp\n'.format(base_dir,event,type))
    outfile.write('mv *.long {0}events/{1}/coreas/{2}/\n'.format(base_dir,event,type))
    #outfile.write('/vol/optcoma/pycrtools/LORA_simulation/DAT2txt DAT$RUNNR DAT$RUNNR.tmp\n')
    #outfile.write('/vol/optcoma/pycrtools/LORA_simulation/LORA_simulation DAT$RUNNR.tmp DAT$RUNNR.lora\n')
    #outfile.write('rm DAT$RUNNR.tmp\n')
    outfile.write('mv SIM$RUNNR.reas {0}events/{1}/coreas/{2}/steering/SIM$RUNNR.reas\n'.format(base_dir,event,type))
    outfile.write('mv SIM$RUNNR.list {0}events/{1}/coreas/{2}/steering/SIM$RUNNR.list\n'.format(base_dir,event,type))
    outfile.write('cp -r * {0}events/{1}/coreas/{2}/\n'.format(base_dir,event,type))
    outfile.write('rm -rf /scratch/kmulrey/{0}/{1}/$RUNNR/*\n'.format(event,part_id))

    outfile.close()



#event=[60409606,61271909,63246671,65490891,65214973,64960703,70904490,80495081,85083852,84432712,83467151]
#azimuth=[57.8777037378,5.62278563,24.5338416226,85.323975652,-45.2948503091,-11.61947731,25.876497776,-170.239430922,130.84941074,14.554336601,-32.732258602]
#zenith=[28.7984121752,36.3926767814,17.6159726846,28.5091835026,38.2803974379,35.5384552615,14.4052883375,40.6923782667,20.8766959694,39.5359400205,47.2151295936]
#energy=[159421449.975,656804646.049,470693726.644,102292495.459,325078049.841,151927036.517,127390466.343,319912958.115,101648514.077,273404353.855,201546777.962]
#seed=[43356,6726,53394,72474,38098,15078,60740,25706,51352,37712,91526]


event=[150115145,188215122,100201518,207028458,113916058,153763491,127818680,149538803,193057293,123439073,207593143,177113844,203043496,123314680,134522987,193286518,202265278,206408479]

azimuth=[243.358,348.360,184.474,23.287,132.083,109.010 ,149.408,133.552,95.833,186.568,65.901,1.317,56.750,116.460,100.776,342.553,193.960,28.006]

zenith=[26.8184875341,28.1175462672,42.3372983422,52.8230718319,40.5384025412,33.3745994365,46.5667172663,42.1285255204,26.3185841019,21.7540816693,15.1280779128,22.443394534,30.2746094802,57.3195671292,53.3000500667,26.2660271676,27.388715137,41.7229692358]

energy=[208937740.163,163095093.566,130226748.638,346015958.684,208136732.35,53876280.0,436526843.101,184555120.311,396401102.006,109545194.746,390078470.267,79887320.3463,61392638.8162,720613353.187,320084652.997,145116318.801,117853548.095,66451906.3009]

seed=[1663,1620,2152,9844,2507,8891,6154,4792,1120,4713,1308,9386,896,555,2376,1129,2302,2021]






#azimuth= -149.856203773
#zenith=25.1606454355
#energy=166070995.766
#seed=25056
for i in np.arange(len(event)):
    write_file(str(int(event[i])),azimuth[i], zenith[i], energy[i], seed[i],'proton')
    write_file(str(int(event[i])),azimuth[i], zenith[i], energy[i], seed[i],'iron')
