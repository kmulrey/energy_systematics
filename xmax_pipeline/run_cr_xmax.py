import numpy as np
import glob, os, sys
from multiprocessing import Pool

file=open('/home/kmulrey/cr_physics_new','r')
events=np.genfromtxt(file)
file.close()

def waitAndHandleErrors(process, name):
    print 'Now running script %s...' % name
    print '----------------------------------------------------'
    #for line in iter(process.stdout.readline,""):
    #sys.stdout.write(line)
    #sys.stdout.flush()
    #output, error = process.communicate()
    process.wait()
    if process.returncode != 0:
        print "Error running {0}:".format(name)
        #        print output
        #print error
        failedlog = open(scripts_directory+'/../log/failed.txt', 'a')
        failedlog.write(str(eventid)+'\n')
        failedlog.close()
        #errorlog = open(scripts_directory+'/../log/error-{0}.txt'.format(eventid), 'w')
        #errorlog.write(error)
        #errorlog.close()
        sys.exit()
    else:
        #print output
        print '{0} finished normally.'.format(name)


BASE_PATH='/vol/astro7/lofar/sim/pipeline/'
RESULTS_PATH='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/results/'
NEWSIM_PATH='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/'
DATA_DIR=BASE_PATH+'events/'
SIMULATION_DIR=BASE_PATH+'run/'
OUTPUT_DIR_RADIO_ONLY=RESULTS_PATH+'production_analysis_radio_only_cal_FINAL/'
LOG_DIR=RESULTS_PATH+'log/'
WRITE_FILT='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/events/'#_thetaMINUS1/
COLLECT_DIR='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/filtered/'#_thetaMINUS1/
LORA_DIR='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/data/'


def run_event(event):

    try:
        runCommand =    'python cr_xmaxfit.py --event={0} --lorafile-suffix=_GeVfix --datadir={1} --simulationdir={2} --iteration=0  --outputdir-radio-only={3} --logdir={4} --filtdir={5} --writedir={6} --collectdir={7} --loradir={8}'.format(int(event),DATA_DIR,SIMULATION_DIR,OUTPUT_DIR_RADIO_ONLY,LOG_DIR,NEWSIM_PATH,WRITE_FILT,COLLECT_DIR,LORA_DIR)
    
        print 'Running command: %s' % runCommand
        #process = subprocess.Popen([runCommand], shell=True)#, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        retcode = os.system(runCommand)
        if retcode != 0:
            print 'Error running fit_analysis_updated.py (radio-only)!'
            sys.exit()
        #waitAndHandleErrors(process, 'fit_analysis_updated.py')
        print 'cr_xmaxfit.py completed.'
    
    
    except:
        print 'issue running {0}'.format(int(event))
    
    
    
    
    
###############

#event=int(events[11])
#event=196796518

#print '--------> event {0}'.format(event)

#use=[int(events[13]),int(events[14]),int(events[15]),int(events[16])]
p = Pool(12)
p.map(run_event,events)
#run_event(event)







'''

p = Pool(14)
p.map(process_and_save, muon_calib_files)






runCommand = '/usr/bin/python -u '+scripts_directory+'/fit_analysis_updated.py --event={0} --iteration={1} --inputdir={2} --outputdir={3} --randomseed={4} --loradir={5} --radio-only-fit {6} {7}'.format(eventid, iteration, collect_outputdir, outputdir_radio_only, randomseed, simulationdir, doFetchLofarData, doRewriteLofarData)

print 'Running command: %s' % runCommand
#process = subprocess.Popen([runCommand], shell=True)#, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
retcode = os.system(runCommand)
if retcode != 0:
    print 'Error running fit_analysis_updated.py (radio-only)!'
    sys.exit()
    #waitAndHandleErrors(process, 'fit_analysis_updated.py')

print 'cr_xmaxfit.py completed.'
'''
