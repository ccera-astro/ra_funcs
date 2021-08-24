#
#
# A collection of functions suitable for some radio astronomy DSP
#  applications
#
import math
import numpy
import ephem
import time
import os
import sys

# Given:
#  frequency(Hz), baselines(meters), declination(deg),
#  latitude(decimal degrees as a float)
#
# Return the fringe period, in seconds
#
#
def fperiod(freq,baseline,decln,latitude):
    
    C=299792000.0
    Lambda = C/freq
    #
    # Convert baseline into fringe-spacing in degrees
    #
    fwidth= (math.degrees(Lambda))/baseline
    
    #
    # 240 seconds (4 minutes) per degree on the celestial equator
    #
    fwidth *= (4.0 * 60.0)
    
    #
    # Adjust for declination and local latitude
    # Takes longer for source to transit through 'fwidth' at higher
    # declinations
    #
    fwidth /= math.cos(math.radians(decln))
    return fwidth

#
# Given coner-frequency in Hz, sample-rate in Hz
#
# Return an appropriate 'alpha' value for a single-pole IIR
#  filter
#
def getalpha(corner, srate):
    q = math.pow(math.e,-2.0*3.14159*(corner/srate))
    alpha = 1.0 - q
    return alpha

#
# Given longitude(decimal degrees as a float)
#
# Return the current sidereal time as a string with
#  "," separated tokens
#
def cur_sidereal(longitude):
    longstr = "%02d" % int(longitude)
    longstr = longstr + ":"
    longitude = abs(longitude)
    frac = longitude - int(longitude)
    frac *= 60
    mins = int(frac)
    longstr += "%02d" % mins
    longstr += ":00"
    x = ephem.Observer()
    x.date = ephem.now()
    x.long = longstr
    jdate = ephem.julian_date(x)
    tokens=str(x.sidereal_time()).split(":")
    hours=int(tokens[0])
    minutes=int(tokens[1])
    seconds=int(float(tokens[2]))
    sidt = "%02d,%02d,%02d" % (hours, minutes, seconds)
    return (sidt)

#
# Return mask in FFTW3 order given list of RFI frequencies, and given FFT size
#
# We provide for both RAW FFT (which will be complex), and post-mag**2,
#  which will be floats
#
def rfi_mask(srate,freq,rfilist,fftsize,iscomplex):
    
    #
    # Make up appropriate array
    #
    if (iscomplex == True):
        rv = [complex(1.0,0.0)]*fftsize
    else:
        rv = [1.0]*fftsize
        
    #
    # Bin width
    #
    binw = float(srate)/float(fftsize)
    
    #
    # Frequency limits
    #
    low = freq-(srate/2.0)
    high = freq+(srate/2.0)
    
    #
    # Determine correct "zero" value to stuff into array
    #
    zerov = complex(0.0,0.0) if iscomplex == True else 0.0
    
    #
    # For each entry in the RFI list
    #
    for r in rfilist:
        
        #
        # Within limits?
        #
        if (r >= low and r <= high):
            
            #
            # Compute the index into the mask
            # 
            # Recall that the ordering in FFTW3 is:
            #  [positive-frequencies,negative-frequencies]
            #
            # So we'll have to do a bit of index manipulation to
            #  place the zero in the correct place in the output array
            #
            #
            ndx = r - freq
            ndx /= binw
            ndx = int(ndx)
            
            #
            # If negative, adjust
            #
            if (ndx < 0):
                ndx += int(fftsize/2)
                rv[ndx] = zerov
            else:
                rv[ndx] = zerov
    return (rv)
