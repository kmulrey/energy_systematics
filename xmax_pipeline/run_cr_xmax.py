import numpy as np
import ROOT, glob, os, sys
from multiprocessing import Pool

file=np.open('/home/kmulrey/cr_physics_new','r')
events=np.genfromtxt(file)
file.close()


print events










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
