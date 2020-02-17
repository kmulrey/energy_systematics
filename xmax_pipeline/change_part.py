
def readLOFARdata(datafile):
    global doFetch, doRewrite # they get updated if info from a data file is requested but doesn't exist
    # Uses parameters from the command line: eventid, doFetch, doRewrite
    if not doFetch and not os.path.exists(datafile): # Want to read from file but it is not there: fetch from database instead, and write to file.
        doFetch = True; doRewrite = True
        print 'LOFAR data file not found; fetching LOFAR data from CR Database and writing to file: %s' % datafile

    if doFetch: # Get info from CR Database using 'dbldf' module
        print 'Fetching LOFAR data from CR Database'
        nofAttempts = 10 # CR Database sometimes doesn't connect, so retry if this happens
        thisAttempt = 0
        readout_OK = False
        while not readout_OK and (thisAttempt < nofAttempts):
            try:
                core_x, core_y , stname, antenna_ids, positions, dist, x_err, signals, power11, power21, power41, rms, noisepower, pulse_delay_fit_residual, time, lora_x, lora_y, lora_dens, az, elev, elev_lora = dbldf.GetLDF(options.event)
                readout_OK = True
            except Exception as e:
                thisAttempt += 1
                print 'Database connection failed at attempt %d, error message is %s' % (thisAttempt, e)
                import time
                time.sleep(20 + 100*thisAttempt*np.random.rand() ) # sleep for 20 to 120 seconds before retrying
        if not readout_OK:
            raise ValueError("Could not get a database connection after %d attempts!" % nofAttempts)
        
        # Convert az/el to az/zenith in radians and with different convention (which?)
        lofarData = (core_x, core_y , stname, antenna_ids, positions, dist, x_err, signals, power11, power21, power41, rms, noisepower, pulse_delay_fit_residual, time, lora_x, lora_y, lora_dens, (450-az)/180.*np.pi, (90-elev)/180.*np.pi,(90-elev_lora)/180.*np.pi)
        if doRewrite: # Rewrite 'datafile' with values just read in from CR Database
            print 'Rewriting file containing LOFAR data'
            f = open(datafile, 'w')
            cPickle.dump((core_x, core_y , stname, antenna_ids, positions, dist, x_err, signals, power11, power21, power41, rms, noisepower, pulse_delay_fit_residual, time, lora_x, lora_y, lora_dens, (450-az)/180.*np.pi, (90-elev)/180.*np.pi,(90-elev_lora)/180.*np.pi),f)
            f.close()

else: # Just read in from pre-produced file
    print 'Reading LOFAR data from pre-produced file'
        f = open(datafile,'r')
        lofarData = cPickle.load(f)
        f.close()
    
    return lofarData
    
    
    
    
    
    
    

(core_x, core_y, station_name, antenna_ids, positions, dist, x_err, signal, dpower11, dpower21, dpower41, rms, noisepower, pulse_delay_fit_residual, data_time, lora_x, lora_y, lora_dens, data_azimuth, data_zenith, lora_zenith) = readLOFARdata(datafile)
