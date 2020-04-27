import numpy as np
import glob, os, sys
from multiprocessing import Pool
from optparse import OptionParser


file=open('/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/energyScale/radio/lofar_events/energy_events.txt','r')
#file=open('/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/energyScale/radio/lofar_events/hadronic_events.txt','r')

events=np.genfromtxt(file)

file.close()
#events=np.asarray([60409606,63246671,65214973,80495081,84432712])

'''
parser = OptionParser()
parser.add_option("-i", "--eventindex", default = "0", help = "Event ID to process")
(options, args) = parser.parse_args()
eventindex = int(options.eventindex)
event=events[eventindex]
'''





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
OUTPUT_DIR_RADIO_ONLY=RESULTS_PATH+'/production_analysis_radio_only_LORA_new_baseline/'
LOG_DIR=RESULTS_PATH+'log/'
WRITE_FILT='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/events/'
COLLECT_DIR='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/filtered/'
LORA_DIR='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/data_baseline/'


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
    
    
    
    
    
########################################################

event=80495081
run_event(event)


#p = Pool(12)#
#p.map(run_event,events)








