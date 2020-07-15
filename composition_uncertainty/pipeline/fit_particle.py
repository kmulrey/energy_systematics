import numpy as np
from optparse import OptionParser
import os
import pickle
import os.path
import re
import sys
from os import listdir
from os.path import isfile, join
from os import walk
import glob
import scipy.interpolate as intp
import matplotlib.pyplot as plt
import scipy.optimize as opt

parser = OptionParser()
#parser.add_option("-i", "--eventindex", type="int", help="event number", default="0")
parser.add_option("-e", "--event", type="int", help="event number", default="127944584")

options, arguments = parser.parse_args()
event=options.event


proton_dir='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/'+event+'/conex/proton/steering/'
helium_dir='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/'+event+'/conex/helium/steering/'
oxygen_dir='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/'+event+'/conex/oxygen/steering/'
iron_dir='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/'+event+'/conex/iron/steering/'

proton_runnr=[]
proton_files=glob.glob('/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/'+event+'/corsika/proton/*.lora')
for i in np.arange(len(proton_files)):
    proton_runnr.append(proton_files[i].split('/')[11].split('.')[0].split('DAT')[1])
helium_runnr=[]
helium_files=glob.glob('/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/'+event+'/corsika/helium/*.lora')
for i in np.arange(len(helium_files)):
    helium_runnr.append(helium_files[i].split('/')[11].split('.')[0].split('DAT')[1])
oxygen_runnr=[]
oxygen_files=glob.glob('/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/'+event+'/corsika/oxygen/*.lora')
for i in np.arange(len(oxygen_files)):
    oxygen_runnr.append(oxygen_files[i].split('/')[11].split('.')[0].split('DAT')[1])
iron_runnr=[]
iron_files=glob.glob('/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/'+event+'/corsika/iron/*.lora')
for i in np.arange(len(iron_files)):
    iron_runnr.append(iron_files[i].split('/')[11].split('.')[0].split('DAT')[1])


print(len(proton_runnr),len(helium_runnr),len(oxygen_runnr),len(iron_runnr))












reco_dir='/vol/astro7/lofar/sim/pipeline/production_analysis_Dec2019/'
dat_files=glob.glob(reco_dir+'*dat')
for file in glob.glob(reco_dir+'*'+str(int(event))+'*dat'):
    reco_file_name=file
print(reco_file_name)

reco_file=open(reco_file_name,"rb")
reco_info=pickle.load(reco_file)
reco_file.close()


def GetUVW(pos, cx, cy, zen, az):
    relpos = pos-np.array([cx,cy,7.6])
    inc=1.1837
    B = np.array([0,np.cos(inc),-np.sin(inc)])
    v = np.array([-np.cos(az)*np.sin(zen),-np.sin(az)*np.sin(zen),-np.cos(zen)])
    #print v
    vxB = np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]])
    vxB = vxB/np.linalg.norm(vxB)
    vxvxB = np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])
    return np.array([np.inner(vxB,relpos),np.inner(vxvxB,relpos),np.inner(v,relpos)]).T


datafile='/vol/astro3/lofar/sim/pipeline/run/data/dbrev'+event+'.dat'
#simfile='/vol/astro7/lofar/sim/pipeline/run/filtered/SIM95166806_2.filt'
simfile='/vol/astro7/lofar/sim/pipeline/run/filtered/SIM'+event+'.filt'

with open(datafile, 'rb') as f:
    u = pickle._Unpickler(f)
    u.encoding = 'latin1'
    p = u.load()
(core_x, core_y, station_name, antenna_ids, positions, dist, x_err, signal, dpower11, dpower21, dpower41, rms, noisepower, pulse_delay_fit_residual, data_time, lora_x, lora_y, lora_dens, data_azimuth, data_zenith, lora_zenith) = p


with open(simfile, 'rb') as f:
    u = pickle._Unpickler(f)
    u.encoding = 'latin1'
    siminfo = u.load()
    
    
nstations=len(lora_x)
lora_positions=np.zeros([nstations,3])
lora_positions[:,0]=lora_x
lora_positions[:,1]=lora_y

lora_dens=np.array(lora_dens)

eff_area=0.9*np.cos(lora_zenith) # 0.9 m2 area
lora_err=np.sqrt(lora_dens*eff_area)/eff_area
for i in np.arange(nstations):
    if (lora_err[i]>lora_dens[i]):
        lora_err[i]=lora_dens[i]
      
lora_err=lora_err[np.nonzero(lora_dens)]
lora_positions=lora_positions[np.nonzero(lora_dens)]
lora_dens=lora_dens[np.nonzero(lora_dens)]

ConversionFactor = eff_area * 6.7 # Density => #particles => energy deposit (MeV)
lora_dens=lora_dens*ConversionFactor
lora_err=lora_err*ConversionFactor

nstations=len(lora_dens)

pos_lora_UVW = GetUVW(lora_positions, coreX, coreY, data_zenith, data_azimuth)
lora_r=np.sqrt(pos_lora_UVW[:,0]*pos_lora_UVW[:,0]+pos_lora_UVW[:,1]*pos_lora_UVW[:,1])


proton_xmax_dir='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/'+event+'/corsika/proton/'
helium_xmax_dir='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/'+event+'/corsika/helium/'
oxygen_xmax_dir='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/'+event+'/corsika/oxygen/'
iron_xmax_dir='/vol/astro7/lofar/kmulrey/sim/composition_uncertainty/events/'+event+'/corsika/iron/'

def returnLORA(RUNNR,dir0,v):
    if v==1:
        lorafile=dir0+'DAT'+str(int(RUNNR)).zfill(6)+'_GeVfix.lora'
    if v==2:
        lorafile=dir0+'DAT'+str(int(RUNNR)).zfill(6)+'.lora'
    longfile=dir0+'DAT'+str(int(RUNNR)).zfill(6)+'.long'
    longfile=open(longfile,'r')
    hold=''
    for line in longfile:
        if "PARAMETERS" in line:
            hold=line
            break
    xmax=float(hold.split()[2:][2])
    #print((hold.split()[2:8]))
    hillas=np.asarray([float(hold.split()[2]),float(hold.split()[3]),float(hold.split()[4]),float(hold.split()[5]),float(hold.split()[6]),float(hold.split()[7])])
    #RUNNR=file.split('/')[12].split('.')[0].split('DAT')[1]
    longfile.close()
    #print(hillas)
    pdata = np.genfromtxt(lorafile)
    cors_r=5*pdata[:,0]+2.5
    cors_dens=pdata[:,1]
    #hillas=0
    r_range=cors_r
    dens_range=cors_dens
    part_intp=intp.interp1d(r_range, dens_range)
    initguess_d_ratio=1.0
    # NKG fit function
    def lora_fitfunc(ratio,lpos,cx,cy,az,zen,part_intp):
        pos_lora_UVW = GetUVW(lpos, cx, cy, zen, az)
        r = np.sqrt(pos_lora_UVW[:,0]*pos_lora_UVW[:,0]+pos_lora_UVW[:,1]*pos_lora_UVW[:,1])
        if (np.max(r)<995): return ratio*part_intp(r)
        else: return 0

    #lora_errfunc_restrained = lambda p,lofar_pow,lofar_err, lofar_pow01,lofar_err01, lofar_pos,lora_data, lora_err, lora_pos, cx, cy, az, zen, n: (lora_fitfunc(p[0], lora_pos,cx+p[1],cy+p[2],az,zen,n) - lora_data)/lora_err

    lora_errfunc_fixedcore_restrained = lambda p,lora_data, lora_err, lora_pos, cx, cy, az, zen: (lora_fitfunc(p, lora_pos,cx,cy,az,zen,part_intp) - lora_dens)/lora_err # cx_offset and cy_offset are the offsets to core position from radio (or combined) fit

    particle_fit_args=(lora_dens,lora_err,lora_positions,core_x,core_y,data_azimuth,data_zenith)
    optimal_d_ratio, covar = opt.leastsq(lora_errfunc_fixedcore_restrained, initguess_d_ratio, args=particle_fit_args)
    return part_intp,hillas,optimal_d_ratio
