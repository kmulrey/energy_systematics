import numpy as np
import glob
import fit_analysis_updated as analysis
import cPickle
from optparse import OptionParser
import os


parser = OptionParser()
parser.add_option("-n", "--event", default = "0", help = "filename of database")
parser.add_option("-t", "--iteration", default=0, help="Iteration number of the simulations to read results from")
parser.add_option("--randomseed", default=2017, help="Set random seed for initial choice of core positions")

parser.add_option("-b", "--inputdir", default="/vol/astro3/lofar/sim/pipeline/run", help="Base directory where input files are located, i.e. in <inputdir>/data and <inputdir>/filtered")

parser.add_option("-r", "--recodir", default="/vol/astro3/lofar/sim/pipeline/test_analysis", help="Directory where reconstruction data (output of fit_analysis) are located")
parser.add_option("-o", "--outputdir", default="/vol/astro3/lofar/sim/pipeline/test_mcvsmc")
parser.add_option("-f", "--fetch-lofardata", default=False, action="store_true", help="Re-fetch LOFAR measured data from CR Database (instead of from previously saved file)")
parser.add_option("-w", "--rewrite-lofardata", default=False, action="store_true", help="Rewrite file containing LOFAR measured data (from CR Database)")
parser.add_option("--radio-only-fit", default=False, action="store_true", help="Perform fit using only radio data; default is to use both radio and particle data")

#parser.add_option("-r", "--randomseed", default=2017, help="Set random seed for initial choice of core positions")

#parser.add_option("-i", "--iteration", default = "0", help = "iteration number")
#parser.add_option("-n", "--event", default = "0", help = "event number")
(options, args) = parser.parse_args()
eventid = int(options.event)
iteration = int(options.iteration)
inputdir = options.inputdir
recodir = options.recodir
outputdir = options.outputdir

doFetch = options.fetch_lofardata # should be set to the same value as when doing fit_analysis on data.
doRewrite = options.rewrite_lofardata
if doRewrite:
    doFetch = True # need to fetch in order to rewrite the file with LOFAR data
radio_only_fit = options.radio_only_fit

randomseed = 2017
analysis.setOptions( (eventid, inputdir, iteration, outputdir, randomseed, doFetch, doRewrite, radio_only_fit) )

#analysis.setFetchRewrite(doFetch, doRewrite) # Set values in fit_analysis

iterationSuffix = '_{0}'.format(iteration) if iteration > 0 else ''
analysis_file = os.path.join(recodir, 'reco{0}{1}.dat'.format(eventid, iterationSuffix))

#if (nit=0): analysis_file='/vol/astro3/lofar/sim/pipeline/test_analysis/reco{0}.dat'.format(ev)
#else: analysis_file='/vol/astro3/lofar/sim/pipeline/test_analysis/reco{0}_{1}.dat'.format(ev,nit)

anfile=open(analysis_file)
aninfo=cPickle.load(anfile)

nsample=3

corex=aninfo['core_x']
corey=aninfo['core_y']

outputfilename = os.path.join(outputdir, 'meth{0}{1}.dat'.format(eventid, iterationSuffix))
#outputfilename='/vol/astro3/lofar/sim/pipeline/run/test_analysis/meth{0}_{1}.dat'.format(ev,nit)

#find out nsim!

#basedir = "/vol/astro3/lofar/sim/pipeline/run"

simfile = os.path.join(inputdir, 'filtered/SIM{0}{1}.filt'.format(eventid, iterationSuffix))
#if (nit==0): simfile = os.path.join(basedir, 'filtered/SIM{0}.filt'.format(ev))
#else: simfile = os.path.join(basedir, 'filtered/SIM{0}_{1}.filt'.format(ev,nit))

g=open(simfile,'r')
siminfo= cPickle.load(g)
primary_type=siminfo['primary'][:,0]
print primary_type
nsimprot=np.sum(primary_type==14)
nsimiron=np.sum(primary_type>14)
l=nsimprot+nsimiron

print nsimprot, nsimiron, l

xreco=np.zeros([l,nsample])
xbest=np.zeros([l,nsample])
xreal=np.zeros([l,nsample])
cchi2=np.zeros([l,nsample])
rchi2=np.zeros([l,nsample])
p_ratio=np.zeros([l,nsample])
p_ratio0=np.zeros([l,nsample])
p_ratio1=np.zeros([l,nsample])
d_ratio=np.zeros([l,nsample])
dratio=np.zeros([l,nsample,40])
sim_tot_power=np.zeros([l,nsample,40,160])
rbf=np.zeros([l,nsample])
d150=np.zeros([l,nsample])
nch=np.zeros([l,nsample])
mr=np.zeros([l,nsample])
sa=np.zeros([l,nsample])
xoff=np.zeros([l,nsample])
yoff=np.zeros([l,nsample])
xcore=np.zeros([l,nsample])
ycore=np.zeros([l,nsample])
xmaxsim=np.zeros([l,nsample,40])
energy=np.zeros([l,nsample])
zenith=np.zeros([l,nsample])
azimuth=np.zeros([l,nsample])
noisepower=np.zeros([l,nsample,2])
antratio=np.zeros([l,nsample])
nsim=np.zeros([l,nsample])
nsim_prot=np.zeros([l,nsample])

#inputdir = "/vol/astro3/lofar/sim/pipeline/run"
#outputdir = "/vol/astro3/lofar/sim/pipeline/run/analysis"
#import pdb; pdb.set_trace()

for i in np.arange(l):
    #i = 12 # one test shower in the middle somewhere
    for j in np.arange(nsample):
        randomseed = j
        #randomseed = eventid * (100*i + j) # Too large, < 2^32 needed
        print 'Doing shower %d/%d, sample %d/%d' % (i, l, j, nsample)
        conv, xreco[i,j], xbest[i,j], dum3, xreal[i,j], cchi2[i,j], rchi2[i,j], fitparam, p_ratio[i,j], p_ratio0[i,j], p_ratio1[i,j], d_ratio[i,j], dum1, dum2,rbf[i,j], d150[i,j], xoff[i,j], yoff[i,j], xcore[i,j], ycore[i,j], energy[i,j], zenith[i,j], azimuth[i,j], noisepower[i,j], antratio[i,j], nsim[i,j], nsim_prot[i,j] = analysis.reverseAnalysis(eventid, iteration, inputdir, outputdir, randomseed, simmode=True, simevent=i, outfile="tmp", saveplt=True, showplt=False, plots=True, simcorex=corex, simcorey=corey, verbose=True)

   
pickfile = open(outputfilename,'w')
cPickle.dump((nsimprot, nsimiron, eventid, xreco, xreal, xbest, cchi2, rchi2, p_ratio, p_ratio0, p_ratio1, d_ratio, dratio, sim_tot_power, rbf, d150, nch, mr, sa, xoff, yoff, xcore, ycore, energy, zenith, azimuth, noisepower, antratio, nsim, nsim_prot),pickfile)


# correct core & correct energy: P_radio = f * P_sim
# correct core & wrong energy: P_radio = f * P_sim * (E_real/E_sim)

# f = P_rad (Er) / P_sim (Es)
# P_sim (E)/E^2 = constant
# P_sim (Er) = P_sim (Es) / Es^2 * Er^2

# const(?) = P_rad(Er) / P_sim (Er) = Prad (Er) / PSim (Es) * Es^2 / Er^2
# = f * Esim^2/ Ereal^2

# f = (P_radio / P_sim) * (E_sim / E_real)
# plot (P_radio / P_sim) * (E_sim / E_new_core) => f should be close to unity

# At high energy, f_corr gets smaller! => F ~ E^-1 = P_data/P_sim ===> E / E^2 ??? 
