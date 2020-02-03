import numpy as np
from scipy.signal import resample

def FreqFilter(data, lowfreq, hifreq, tstep):
    # Bandpass filter on timeseries 'data' with timestep tstep.
    # Data formatted as 2d array, to be processed along the first dimension ('axis 0').
    # Time series is resampled to 5 ns sampling period.
    
    spec = np.fft.rfft(data, axis=0)
    dlength = data.shape[0]
    freqs = np.fft.rfftfreq(dlength, tstep) # Frequency axis, in Hz (for tstep in s).
    freqstep = freqs[1] - freqs[0]
    # Window function = 1 for indices inside our bandpass, 0 otherwise. Apply roundoff by widening the band.
    window_function = (freqs > lowfreq*1e6 - freqstep) & (freqs < hifreq*1e6 + freqstep) # Boolean array
    filtered_spectrum = spec * window_function[:, np.newaxis] # Apply window function on the correct axis
    filtered = np.fft.irfft(filtered_spectrum, axis=0) # Filtered time series
    
    filtered_resampled = resample(filtered, int(np.round(dlength * tstep / 5.0e-9)), axis=-2) # Time series resampled to 5 ns sample spacing

    return filtered_resampled


"""
def FreqFilter(data,lowfreq,hifreq,tstep):
    dlength=data.shape[0]
    spec=np.fft.rfft(data, axis=-2)
    freqhi = 0.5/tstep/1e6 # MHz
    freqstep = freqhi/(dlength/2+1) # MHz
    fb = np.floor(lowfreq/freqstep)
    lb = np.floor(hifreq/freqstep)+1
    window = np.zeros([dlength/2+1,1])
    window[fb:lb+1,:]=1
    import pdb; pdb.set_trace()
    maxfreqbin= np.floor(tstep/5e-9 * dlength/2.)+1
    shortspec=spec[0:maxfreqbin,:]*window[0:maxfreqbin,:]
    filt=np.fft.irfft(shortspec, axis=-2)
    return filt
"""
def GetUVW(pos, cx, cy, cz, zen, az, Binc):
    # Return U, V, W coordinates, along vxB, vx(vxB) and v respectively
    relpos = pos - np.array([cx,cy,cz]) # Position relative to core
    B = np.array([0,np.cos(Binc),-np.sin(Binc)]) # Magnetic field vector
    v = np.array([-np.cos(az)*np.sin(zen),-np.sin(az)*np.sin(zen),-np.cos(zen)]) # Incoming vector, from spherical coordinates (note minus sign)

    #vxB = np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]])
    vxB = np.cross(v, B)
    vxB /= np.linalg.norm(vxB)
    #vxvxB = np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])
    vxvxB = np.cross(v, vxB)
    
    return np.array([np.inner(vxB,relpos),np.inner(vxvxB,relpos),np.inner(v,relpos)]).T

def GetAlpha(zen,az,Binc):
    B = np.array([0,np.cos(Binc),-np.sin(Binc)])
    v = np.array([-np.cos(az)*np.sin(zen),-np.sin(az)*np.sin(zen),-np.cos(zen)])
    #vxB = np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]])
    vxB = np.cross(v, B)
    vxB = vxB/np.linalg.norm(vxB)
    #vxvxB = np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])
    vxvxB = np.cross(v, vxB)
    
    return np.arccos(np.inner(np.asarray(B), np.asarray(v)) / (np.linalg.norm(B) * np.linalg.norm(v)))

def stokes_parameters(x, y, hx, hy):
    """Stokes parameters given timeseries *x*, *y* in two orthogonal
polarisations and their Hilbert transforms *hx* and *hy*. The *x* and
*y* axis are along vxB and vxvxB respectively.
    """
    n = x.shape[0]

    I = (1./n) * np.sum(x*x + hx*hx + y*y + hy*hy)
    Q = (1./n) * np.sum(x*x + hx*hx - y*y - hy*hy)
    U = (2./n) * np.sum(x*y + hx*hy)
    V = (2./n) * np.sum(hx*y - x*hy)

    return np.array([I, Q, U, V])

def polarization_angle(S):

    Q = S[1]
    U = S[2]

    psi = (1./2.)*np.arctan2(U, Q)

    return psi
