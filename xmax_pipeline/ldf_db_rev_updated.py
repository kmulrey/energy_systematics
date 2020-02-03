import pycrtools as cr
from optparse import OptionParser
import psycopg2
import cPickle as pickle
import re
import numpy as np
import time
#import matplotlib.ticker as Ticker
from pycrtools import crdatabase as crdb
from pycrtools.tasks import Task
#from ROOT import *

Task.task_write_parfiles = False

# -----------
# Script to draw LDF from 2nd stage pipeline, total intensity as well as footprint, needs database specified"

def unpickle_parameter(db_parameter):

    if db_parameter is None:
        return "-1"
    else:
        return pickle.loads(re.sub('"', "'", str(db_parameter)))
        

#----------------------------------------------------------------------------
# Using cuts on Good quality of LORA reconstruction, will yield true when reconstruction is reliable

def LoraQualityPassed(moliere, corex, corey, elevation):
    
    quality = False 
    
    if np.sqrt(abs(corex)**2 + abs(corey)**2) < 150:
        if moliere < 100 and moliere > 20:
            if elevation > 55:
                quality = True
                
    return quality         

def GetDistance(core,direction,positions):


        """
        Calculating the distance of a station to the shower core in the shower plane.
        """
        
        positions = cr.hArray(positions)
        core = cr.hArray(core)
        
        distances = cr.hArray(copy=positions)
        finaldistance = cr.hArray(copy=positions)

        theta = cr.radians(direction[1])           # spherical coordinates: 0 = horizon, 90 = vertical
        phi = cr.radians(450 - direction[0])        # recalculate to 0 = east counter-clockwise


        distances = positions - core
        axis = cr.hArray([cr.cos(theta)*cr.cos(phi),cr.cos(theta)*cr.sin(phi),cr.sin(theta)])  # shower axis

        finaldistance[...].crossproduct(distances[...],axis)
        dist = cr.hVectorLength(finaldistance[...])

        return dist

def GetTotalDistanceUncertainty(core,coreuncertainties,positions,direction,directionuncertainties,distances):

        """
        Propagating the uncertatinties of the core and the directions into the uncertainties of the distance to the shower axis.

        """

        # Assuming differences in z-coordinates = 0

        dist = distances.vec()
        positions = cr.hArray(positions)
        
        pos1 = cr.hArray_transpose(positions)[0].vec()
        pos2 = cr.hArray_transpose(positions)[1].vec()

        core1 = float(core[0])
        core2 = float(core[1])

        phi = cr.radians(450 -direction[0])        # recalculate to 0 = east counter-clockwise
        theta = cr.radians(direction[1])           # spherical coordinates: 0 = horizon, 90 = vertical


        ecore1 = float(coreuncertainties[0])
        ecore2 = float(coreuncertainties[1])
        covcore = float(coreuncertainties[2])

        ephi = cr.radians(float(directionuncertainties[0]))
        etheta = cr.radians(float(directionuncertainties[1]))
        covangles = float(directionuncertainties[2])

        # short cuts

        cost = cr.cos(theta)
        sint = cr.sin(theta)
        cos2t = cr.cos(theta)**2
        sin2t = cr.sin(theta)**2
        cosp = cr.cos(phi)
        sinp = cr.sin(phi)
        cos2p = cr.cos(phi)**2
        sin2p = cr.sin(phi)**2

        m1 = pos1 - core1
        m2 = pos2 - core2

        sqrfac = 0.5/dist

        # partical derivatives

        diffc1 = -2.*(m1*sin2t + m1*cos2t*sin2p - m2*cos2t*sinp*cosp) * sqrfac
        diffc2 = -2.*(m2*sin2t + m2*cos2t*cos2p - m1*cos2t*sinp*cosp) * sqrfac
        diffp = 2.*(m1*m1*cos2t*sinp*cosp - m2*m2*cos2t*cosp*sinp- m1*m2*cos2t*(cos2p-sin2p)) *sqrfac

        difft = 2.*((m2*m2)*sint*cost + (m1*m1)*sint*cost - (m1*m1)*sin2p*cost*sint - (m2*m2)*cos2p*cost*sint + m1*m2*sint*cost*sinp*cosp) * sqrfac

        # adding contributions

        err =        diffc1*diffc1 * ecore1*ecore1
        err = err +  diffc2*diffc2 * ecore2*ecore2
        err = err +  diffp*diffp * ephi*ephi
        err = err +  difft*difft * etheta*etheta
        err = err + 2* difft*diffp *covangles
        err = err + 2* diffc1*diffc2 *covcore

        err.sqrt()

        return err


def GetLDF(eventno):     
    # Open CR Database connection to read out data pipeline parameters for this event
    #    nofAttempts = 3
    #thisAttempt = 0
    #database_connection = False
    #while not database_connection and (thisAttempt < nofAttempts):
    #    try:
    #            database_connection = True
    #    except:
    #        thisAttempt += 1
    #        print 'Database connection failed at attempt %d' % thisAttempt
    #        import time
    #        time.sleep(20 + 100*np.random.rand() ) # sleep for 20 to 120 seconds before retrying

#    if not database_connection:
#        raise ValueError("No database connection after {0} attempts!".format(nofAttempts))

    dbManager = crdb.CRDatabase(host='coma00.science.ru.nl', user='crdb', password='crdb', dbname='crdb')

    db = dbManager.db
    event = crdb.Event(db = db, id = eventno)

    """gotEvent = False
        nofAttempts = 3
        thisAttempt = 0
        while not gotEvent and (thisAttempt < nofAttempts):
            try:
                gotEvent = True
            except:
                thisAttempt += 1
                print 'Get event failed at attempt %d' % thisAttempt
                import time
                time.sleep(20 + 100*np.random.rand() )

        if not gotEvent:
            raise ValueError("Could not read event data after {0} attempts!".format(nofAttempts))
    """
    lora_elevation=0
    try:
        energy = event["lora_energy"]
        core_x = event["lora_core_x"]
        core_y = event["lora_core_y"]
        azimuth = event["lora_azimuth"]
        lora_elevation = event["lora_elevation"] 
        moliere = event["lora_moliere"]
        lora_x = event["lora_posx"]
        lora_y = event["lora_posy"]
        lora_dens = event["lora_particle_density__m2"]

        if LoraQualityPassed(moliere, core_x, core_y, lora_elevation):
            print "GOOD quality event (LORA)", energy, azimuth, lora_elevation, core_x, core_y, moliere
        else:
            print "LORA Quality not passed", energy, azimuth, lora_elevation, core_x, core_y, moliere
    except:
        print "Skipping event, no LORA data"
        return -1 # This should be tested for when calling the function... Nothing to return for a file without LORA info.
 
    # Loop over all stations in event
    stations = []
    for f in event.datafiles:
       stations.extend(f.stations)

    positions = []
    selected_dipoles = []
    delays = []
    amplitude = []
    rms = []
    power11 = []
    power21 =[]
    power41 =[]
    noisepower = []
    stationname = []
    pulse_direction = []
    pulse_delay_fit_residual = []
    time = []
    for station in stations:
        if station.status == "GOOD":
            try:
                p = station.polarization["xyz"]
                p0 = station.polarization["0"]
                positions.append(station["local_antenna_positions"])
                selected_dipoles.append(station["crp_selected_dipoles"])
                amplitude.append(p["crp_pulse_peak_amplitude"])
                power11.append(p0["crp_integrated_pulse_power"])
                power21.append(p0["crp_integrated_pulse_power_wide"])
                power41.append(p0["crp_integrated_pulse_power_double_wide"])
                noisepower.append(p0["crp_integrated_noise_power"])
                rms.append(p["crp_rms"])
                stationname.append([station.stationname]*len(p["crp_rms"]))
                pulse_direction.append(station["crp_pulse_direction"])
                pulse_delay_fit_residual.append(station["crp_pulse_delay_fit_residual"])
                time.append(station["crp_pulse_time"])
            except:
                print "Do not have all pulse parameters for station", station.stationname

    #print positions
    positions = np.vstack(positions)
    selected_dipoles = np.hstack(selected_dipoles)
    amplitude = np.vstack(amplitude)
    rms = np.vstack(rms)
    power11 = np.vstack(power11)
    power21 = np.vstack(power21)
    power41 = np.vstack(power41)
    noisepower = np.vstack(noisepower)
    pulse_delay_fit_residual = np.hstack(pulse_delay_fit_residual)
    #stationname = np.array(stationname)
    stationname=np.array(sum(stationname,[]))
    time = np.vstack(time)
    #print time[:,0]
    lof_azimuth=np.mean(pulse_direction,axis=0)[0]
    lof_elevation=np.mean(pulse_direction,axis=0)[1]

    shape = positions.shape
    positions = positions.reshape((shape[0] / 2, 2, shape[1]))[:,0]
    positions = positions.copy()

    #uncer = np.sqrt(rms[:,0]**2+rms[:,1]**2 + rms[:,2]**2)
    #total = np.sqrt(amplitude[:,0]*signals[:,0]+signals[:,1]*signals[:,1]+signals[:,2]*signals[:,2])

    #print "Returning event:", eventno
    dist = GetDistance([core_x,core_y,0],[azimuth,lora_elevation],positions)
    x_err =  GetTotalDistanceUncertainty([core_x,core_y,0],[5.,5.,0],positions,[azimuth,lora_elevation],[2,2,0],dist)

    dist = np.array(dist)
    x_err = np.array(x_err)
    #print positions.shape, amplitude.shape, power.shape, rms.shape
    return core_x, core_y , stationname, selected_dipoles, positions, dist, x_err, amplitude, power11, power21, power41, rms, noisepower, pulse_delay_fit_residual, time[:,0], lora_x, lora_y, lora_dens, lof_azimuth, lof_elevation, float(lora_elevation)

           
             
  
