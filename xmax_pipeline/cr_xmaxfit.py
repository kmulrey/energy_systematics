# Basic setup for CR Xmax-fit pipeline using updated scripts from Stijn
# Author: A. Corstanje, Jan 2017 (a.corstanje@astro.ru.nl)
#
# TODO: set fixed directories for the scripts called below. Now assumed to be in same directory.

import os
import glob
import sys
from optparse import OptionParser
import subprocess
import pycrtools as cr
from pycrtools import crdatabase as crdb
from pycrtools.tasks import Task

# Directory (fixed) where the scripts are located
#scripts_directory = '/vol/astro3/lofar/sim/pipeline/scripts' # run here also for data on /vol/astro7
scripts_directory = '/vol/optcoma/pycrtools/src/PyCRTools/xmax_pipeline' # run here also for data on /vol/astro7

parser = OptionParser()
parser.add_option("-n", "--event", default = "0", help = "filename of database")
parser.add_option("-t", "--iteration", default="latest", help="Iteration number of the simulations to read results from in fit_analysis. Specify a number, or 'latest'.")
parser.add_option("-d", "--datadir", default="/vol/astro3/lofar/sim/pipeline/events", help="Base dir where the simulated events are stored")
parser.add_option("-i", "--simulationdir", default="/vol/astro3/lofar/sim/pipeline/run", help="Base directory where pre-processed LOFAR and simulations data are stored, i.e. in <simulationdir>/data and <simulationdir>/filtered")
parser.add_option("-o", "--outputdir", default="/vol/astro3/lofar/sim/pipeline/test_analysis", help="Output dir for analysis results") # set to pipeline/run/analysis for production run (or create another dir)
parser.add_option("--outputdir-radio-only", default="/vol/astro7/lofar/sim/pipeline/production_analysis_radio_only", help="Output dir for analysis results with RADIO ONLY fit") # set to pipeline/run/analysis for production run (or create another dir)
parser.add_option("-m", "--mcvsmcdir", default="/vol/astro3/lofar/sim/pipeline/test_mcvsmc", help="Output dir for MC-vs-MC analysis")
parser.add_option("-l", "--logdir", default="/vol/astro3/lofar/sim/pipeline/log", help="Output dir for log-files")
parser.add_option("-r", "--randomseed", default=2017, help="Set random seed for initial choice of core positions")
parser.add_option("-f", "--fetch-lofardata", default=False, action="store_true", help="Re-fetch LOFAR measured data from CR Database (instead of from previously saved file)")
parser.add_option("-w", "--rewrite-lofardata", default=False, action="store_true", help="Rewrite file containing LOFAR measured data (from CR Database)")
parser.add_option("--test", default=False, action="store_true", help="Run debugging test")
parser.add_option("--mcvsmc", default=False, action="store_true", help="Run the MC vs MC analysis (slow)")
parser.add_option("--radio-only-mcvsmc", default=False, action="store_true", help="Use radio-only fit for MC-vs-MC")
parser.add_option("--lorafile-suffix", default="", help="Optional suffix for LORA simulation files e.g. DAT000001_GeVfix.lora (testing)")

parser.add_option("--force-reprocess", default=False, action="store_true", help="Force reprocessing of simulation files to produce .filt files (filterjobs_perevent). Required if setting a different lorafile-suffix")
parser.add_option("--debug-lofar-pulse", default=False, action="store_true", help="Replace Coreas simulation data by the LOFAR pulse (calibrated XYZ from cr_physics pipeline) from one antenna")
parser.add_option("--set-status-done", default=False, action="store_true", help="Set event status to XMAXFIT_DONE even when no mcvsmc is (re)done (not recommended)")

# Get command-line parameters
(options, args) = parser.parse_args()
eventid = int(options.event)
iteration = options.iteration
if iteration != 'latest': # which iteration is the latest, will be handled per event in the scripts below
    iteration = int(iteration)
datadir = options.datadir
simulationdir = options.simulationdir
outputdir = options.outputdir
outputdir_radio_only = options.outputdir_radio_only
mcvsmcdir = options.mcvsmcdir
logdir = options.logdir
randomseed = int(options.randomseed)
doFetch = options.fetch_lofardata
doRewrite = options.rewrite_lofardata
doMCvsMC = options.mcvsmc
radio_only_mcvsmc = options.radio_only_mcvsmc

logfile = os.path.join(logdir, "cr_xmaxfit-{0}-pipeline.txt".format(eventid))

def getLatestIteration():
    def RepresentsInt(s): # check if string is convertible to 'int'...
        try:
            int(s)
            return True
        except ValueError:
            return False

    eventdir = os.path.join(datadir, "{0}".format(eventid))
    # Get all available iterations
    all_iterations = glob.glob(os.path.join(eventdir, '*'))
    all_iterations = [int(os.path.split(x)[1]) for x in all_iterations if RepresentsInt(os.path.split(x)[1])] # make integers
    iteration_to_process = max(all_iterations)

    return iteration_to_process

if iteration == 'latest':
    iteration = getLatestIteration()

iterationSuffix = '_{0}'.format(iteration) if iteration > 0 else ''


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

additional_flags = ""
if options.lorafile_suffix != "":
    additional_flags += "--lorafile-suffix="+options.lorafile_suffix
if options.force_reprocess:
    additional_flags += " --force-reprocess"
if options.debug_lofar_pulse:
    additional_flags += " --debug-lofar-pulse"

# Run filterjobs_perevent in a subprocess and wait for it to finish
runCommand = 'python -u '+scripts_directory+'/filterjobs_perevent.py --eventid={0} --datadir={1} {3} > {2}'.format(eventid, datadir, logfile, additional_flags)
print 'Running command: %s' % runCommand
process = subprocess.Popen([runCommand], shell=True)#, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
waitAndHandleErrors(process, 'filterjobs_perevent.py')

# Run collectfiles_perevent.py
collect_outputdir = os.path.join(options.simulationdir, 'filtered') # subdirectory 'filtered' for combined simulation results

runCommand = 'python -u '+scripts_directory+'/collectfiles_perevent.py --event={0} --iteration={1} --outputdir={2} --datadir={3} >> {4}'.format(eventid, "all", collect_outputdir, datadir, logfile)
print 'Running command: %s' % runCommand
process = subprocess.Popen([runCommand], shell=True)#), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
waitAndHandleErrors(process, 'collectfiles_perevent.py')

doFetchLofarData = '--fetch-lofardata' if doFetch else ''
doRewriteLofarData = '--rewrite-lofardata' if doRewrite else ''

# Run the Xmax fit analysis with RADIO ONLY fit procedure
if os.path.exists(os.path.join(outputdir_radio_only, 'reco{0}{1}.dat'.format(eventid, iterationSuffix))):
    print 'Fit analysis (radio-only) already done for event %d iteration %d, skipping...' % (eventid, iteration)
else:
    runCommand = '/usr/bin/python -u '+scripts_directory+'/fit_analysis_updated.py --event={0} --iteration={1} --inputdir={2} --outputdir={3} --randomseed={4} --radio-only-fit {5} {6} >> {7}'.format(eventid, iteration, simulationdir, outputdir_radio_only, randomseed, doFetchLofarData, doRewriteLofarData, logfile)
    print 'Running command: %s' % runCommand
    #process = subprocess.Popen([runCommand], shell=True)#, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    retcode = os.system(runCommand)
    if retcode != 0:
        print 'Error running fit_analysis_updated.py (radio-only)!'
        sys.exit()
    #waitAndHandleErrors(process, 'fit_analysis_updated.py')

# Run the Xmax fit analysis with COMBINED radio & particles fit procedure
if os.path.exists(os.path.join(outputdir, 'reco{0}{1}.dat'.format(eventid, iterationSuffix))):
    print 'Fit analysis (combined fit) already done for event %d iteration %d, skipping...' % (eventid, iteration)
else:
    runCommand = '/usr/bin/python -u '+scripts_directory+'/fit_analysis_updated.py --event={0} --iteration={1} --inputdir={2} --outputdir={3} --randomseed={4} {5} {6} >> {7}'.format(eventid, iteration, simulationdir, outputdir, randomseed, doFetchLofarData, doRewriteLofarData, logfile)
    print 'Running command: %s' % runCommand
    #process = subprocess.Popen([runCommand], shell=True)#, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    retcode = os.system(runCommand)
    if retcode != 0:
        print 'Error running fit_analysis_updated.py (combined fit)!'
        sys.exit()

if doMCvsMC:
    # Run the MC-vs-MC analysis if not already done for this event
    if os.path.exists(os.path.join(mcvsmcdir, 'meth{0}{1}.dat'.format(eventid, iterationSuffix))):
        print 'MC vs MC already done for event %d iteration %d, skipping...' % (eventid, iteration)
    else: # do mc vs mc analysis
        if doRewriteLofarData:
            doFetchLofarData = ''
        doRadioOnly = '--radio-only-fit' if radio_only_mcvsmc else ''
        runCommand = '/usr/bin/python -u '+scripts_directory+'/mcvsmc_updated.py --event={0} --iteration={1} {2} {3} --inputdir={4} --recodir={5} --outputdir={6} >> {7}'.format(eventid, iteration, doFetchLofarData, doRadioOnly, simulationdir, outputdir, mcvsmcdir, logfile)
        print 'Running command: %s' % runCommand
        #process = subprocess.Popen([runCommand], shell=True) #, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        retcode = os.system(runCommand)
        if retcode != 0:
            print 'Error running mcvsmc_updated.py!'
            sys.exit()
elif not options.set_status_done:
    sys.exit()
# Only set status to XMAXFIT_DONE if mc-vs-mc has been done (unless overridden by option)


#waitAndHandleErrors(process, 'mcvsmc_updated.py')

dbManager = crdb.CRDatabase("crdb", host="coma00.science.ru.nl", user="crdb", password="crdb", dbname="crdb")
db = dbManager.db

print 'Reading database event data for eventid %d ...' % eventid
event = crdb.Event(db=db, id=eventid)
print 'done'

event.simulation_status = "XMAXFIT_DONE"
event.write()

#eventid = 212496991
"""
print 'Reading database event data for eventid %d ...' % eventid
event = crdb.Event(db=db, id=eventid)
print 'done'
    
thisIteration = event['simulation_current_iteration']
thisSimulationStatus = event.simulation_status

print 'Setting simulation_status from %s to XMAXFIT_DONE' % thisSimulationStatus
event.simulation_status = "XMAXFIT_DONE"
event.write()
"""
print 'cr_xmaxfit.py completed.'
