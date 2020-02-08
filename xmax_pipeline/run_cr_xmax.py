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


def run_event(event):


    runCommand =    'python cr_xmaxfit.py --event={0} --lorafile-suffix=_GeVfix --datadir=$DATA_DIR --simulationdir=$SIMULATION_DIR --iteration=0  --outputdir-radio-only=$OUTPUT_DIR_RADIO_ONLY --mcvsmcdir=$MCVSMC_DIR --logdir=$LOG_DIR --filtdir=$NEWSIM_PATH --writedir=$WRITE_FILT --collectdir=$COLLECT_DIR'.format(event)
    
    print 'Running command: %s' % runCommand
    #process = subprocess.Popen([runCommand], shell=True)#, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    retcode = os.system(runCommand)
    if retcode != 0:
        print 'Error running fit_analysis_updated.py (radio-only)!'
        sys.exit()
        #waitAndHandleErrors(process, 'fit_analysis_updated.py')

    print 'cr_xmaxfit.py completed.'

event=int(events[0])

print '--------> event {0}'.format(event)

run_event(event)







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
