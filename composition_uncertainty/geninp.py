from optparse import OptionParser
import os

# Parse commandline options
parser = OptionParser()
parser.add_option("-r", "--runnumber", default = "0", help = "number of run")
parser.add_option("-s", "--seed", default = "1", help = "seed 1")
parser.add_option("-u", "--energy", default = "1e8", help = "cr energy (GeV)")
parser.add_option("-t", "--type", default = "14", help = "particle type ")
parser.add_option("-a", "--azimuth", default = "0", help = "azimuth (degrees; AUGER definition)")
parser.add_option("-z", "--zenith", default = "0", help = "zenith (degrees; AUGER definition) ")
parser.add_option("-d", "--dir", default = "/vol/astro/lofar/sim/test/", help = "output dir")
parser.add_option("-c", "--conex", default = False, action="store_true", help="run CONEX")
parser.add_option("-p", "--parallel", default=False, action="store_true", help="Run in parallel (MPI) so set PARALLEL flag")
parser.add_option("--parallel-lowcut", default=1000.0, help="Low cutoff for parallellization of subshowers (GeV)")
parser.add_option("--parallel-highcut", default=1000000.0, help="High cutoff for parallellization of subshowers (GeV)")
parser.add_option("--atmosphere", default=False, action="store_true", help="Use atmospheric data file with refractive index profile, for new version of CoREAS")
parser.add_option("--atmfile", default="ATMOSPHERE_USSTD.DAT", help="File containing the 5 layer parameters (a, b, c) for the atmosphere, and the full refractive-index profile")
parser.add_option("--atmlayerfile", default=None, help="File containing atmosphere layer definitions (1 line)")
parser.add_option("--first-interaction", default=None, help="Set first interaction height (cm), optional")

(options, args) = parser.parse_args()

print "RUNNR", int(options.runnumber),"                               run number" 
print "EVTNR   1                              number of first shower event"
print "SEED",int(options.seed)+10*int(options.runnumber)," 0 0"
print "SEED",int(options.seed)+10*int(options.runnumber)+1," 0 0"
print "SEED",int(options.seed)+10*int(options.runnumber)+2," 0 0"
print "NSHOW   1                              number of showers to generate"
print "PRMPAR",options.type,"                             primary particle"
print "ERANGE",float(options.energy),float(options.energy),"                    range of energy" 
print "THETAP",float(options.zenith),float(options.zenith),"                     range of zenith angle (degree)"
print "PHIP", -270+float(options.azimuth), -270+float(options.azimuth),"                      range of azimuth angle (degree)"
print "ECUTS   3.000E-01 3.000E-01 4.010E-04 4.010E-04"
print "ELMFLG  T   T                          em. interaction flags (NKG,EGS)"
print "THIN    1.000E-06",1e-6*float(options.energy)," 0.000E+00"
print "THINH   1.000E+02 1.000E+02"
print "OBSLEV  760                       observation level (in cm)"
print "ECTMAP  1.E5                           cut on gamma factor for printout"
print "STEPFC  1.0                            mult. scattering step length fact."
print "MUADDI  T                              additional info for muons"
print "MUMULT  T                              muon multiple scattering angle"
print "MAXPRT  1                              max. number of printed events"
print "MAGNET  18.6   45.6                   magnetic field auger"
print "LONGI   T  1.  T  T                    longit.distr. & step size & fit & out"
print "RADNKG  5.e5                           outer radius for NKG lat.dens.distr."
print "DIRECT", options.dir,"                             output directory"
print "DATBAS  F                              write .dbase file"
print "USER    kmulrey                     user"


if options.atmosphere:
    if options.atmlayerfile is not None:
        infile = open(options.atmlayerfile, 'r')
        atmlay_line = infile.readline()
        infile.close()
        print atmlay_line.split('\n')[0]
    
    #print "ATMLAY 400000   1000000   4000000   10000000                           height of bottom of layers 2 to 5 (cm)"
    print "ATMFILE " + options.atmfile

print "PAROUT  F F"
print "CASCADE T T T"
print "EXIT                                   terminates input"
