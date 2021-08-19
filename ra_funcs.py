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
#  latitude(deg)
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
    q = math.pow(math.e,-2.0*(corner/srate))
    alpha = 1.0 - q
    return alpha

#
# Given longitude(deg)
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
