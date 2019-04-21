#!/usr/env python


###############################################################################
# Some small functions that do conversion between different coordinate systems,
# as well as do projections which can be used for plotting such as galactic map.
#
# To be consistent, all functions inputs and outputs are in RADIANS.
# You can always convert degrees to radians quickly using math.radians() before
# passing them to the functions.
#
# eq2cart(ra, dec, r) Equatorial to Cartesian
# cart2eq(x, y, z)    Cartesian to Equatorial
# eq2gal(ra, dec)     Equatorial to Galactic
# gal2eq(l, b)        Galactic to Equatorial
#
# Here are three different types of projection.
# Notice that all three functions takes traditional longitude and latitude
# definition in radians, which means lon is in [-pi:pi] from the meridian,
# and lat is in [-pi/2:pi/2] from the equator.
# So, for example, if you would like to have galactic map centered on the
# galactic center, before passing galactic longitude l to the function you
# should do:
# l = l if l <= math.pi else l - 2 * math.pi
#
# mollweide(lon, lat) Mollweide Projection
# aitoff(lon, lat)    Aitoff Projection
# hammer(lon, lat)    Hammer Projection
#
# by MQQ @ Vanderbilt   
# last modified: Nov 09 2011
###############################################################################

import os
import sys
import math
from operator import itemgetter, attrgetter

# RA(radians),Dec(radians),distance(kpc) of Galactic center in J2000
Galactic_Center_Equatorial=(math.radians(266.40510), math.radians(-28.936175), 8.33)

# RA(radians),Dec(radians) of Galactic Northpole in J2000
Galactic_Northpole_Equatorial=(math.radians(192.859508), math.radians(27.128336))

# ################################################################
# Convert equatorial coordinates to cartesian coordinates
# ################################################################
def eq2cart(ra, dec, r):
    """
    Convert Equatorial coordinates to Cartesian coordinates.
    Return a tuple (x, y, z) in the same unit of the input distance. 
    
    Keywords arguments:
    ra  -- Right Ascension (in radians)
    dec -- Declination (in radians)
    r   -- Distance
    
    """

    x = r * math.cos(ra) * math.cos(dec)
    y = r * math.sin(ra) * math.cos(dec)
    z = r * math.sin(dec)

    return x, y, z
# #The end of eq2cart ###################################################


# #############################################################
# Convert cartesian to equatorial, return RA and Dec in radians
# ##################################################################
def cart2eq(x, y, z):
    """
    Convert Cartesian coordinates to Equatorial coordinates
    
    Keywords arguments:
    x -- x coordinate
    y -- y coordinate
    z -- z coordinate

    Return a tuple (ra, dec, z):
    ra  -- Right Ascension (in radians)
    dec -- Declination (in radians)
    r   -- Distance in the same unit of input (x, y, z)
    """

    r = math.sqrt(x * x + y * y + z * z)
    ra = math.atan2(y, x)
    ra = ra if ra >= 0 else (2.0 * math.pi + ra)
    dec = math.asin(z / r)

    return ra, dec, r
# #The end of cart2eq #################################################


# #######################################################################
# Convert Equatorial coordinates to Galactic coordinates in the epoch J2000
# #########################################################################
def eq2gal(ra,dec):
    """
    Convert Equatorial coordinates to Galactic Coordinates in the epch J2000.
    
    Keywords arguments:
    ra  -- Right Ascension (in radians)
    dec -- Declination (in radians)

    Return a tuple (l, b):
    l -- Galactic longitude (in radians)
    b -- Galactic latitude (in radians)
    """

    alpha = Galactic_Northpole_Equatorial[0]
    delta = Galactic_Northpole_Equatorial[1]
    la = math.radians(33.0)
    
    b = math.asin(math.sin(dec) * math.sin(delta) +
                  math.cos(dec) * math.cos(delta) * math.cos(ra - alpha))

    l = math.atan2(math.sin(dec) * math.cos(delta) - 
                   math.cos(dec) * math.sin(delta) * math.cos(ra - alpha), 
                   math.cos(dec) * math.sin(ra - alpha)
                   ) + la

    l = l if l >= 0 else (l + math.pi * 2.0)

    l = l % (2.0 * math.pi)

    return l, b
# #The end of eq2gal   ###################################################


# ########################################################################
# Convert Galactic coordinates to Equatorial coordinates in the epoch J2000
# ########################################################################
def gal2eq(l, b):
    """
    Convert Galatic coordinates to Equatorial Coordinates in the epch J2000.
    
    Keywords arguments:
    l -- Galactic longitude (in radians)
    b -- Galactic latitude (in radians)

    Return a tuple (ra, dec):
    ra  -- Right Ascension (in radians)
    dec -- Declination (in radians)
    """

    alpha = Galactic_Northpole_Equatorial[0]
    delta = Galactic_Northpole_Equatorial[1]
    la = math.radians(33.0)

    dec = math.asin(math.sin(b) * math.sin(delta) +
                    math.cos(b) * math.cos(delta) * math.sin(l - la))

    ra = math.atan2(math.cos(b) * math.cos(l - la), 
                    math.sin(b) * math.cos(delta) - 
                    math.cos(b) * math.sin(delta) * math.sin(l - la) 
                    ) + alpha

    ra = ra if ra>=0 else (ra + math.pi * 2.0)

    ra = ra % (2.0 * math.pi)
    
    return ra, dec
# #The end of gal2eq ####################################################


# #####################################################################
# Mollweide Projection
# http://en.wikipedia.org/wiki/Mollweide_projection
# #####################################################################
def mollweide(lon, lat):
    """
    Make Mollweide map projection.
    
    Take traditional longitude and latitude in radians and return a
    tuple (x, y).
    
    Notice that traditionally longitude is in [-pi:pi] from the meridian,
    and latitude is in [-pi/2:pi/2] from the equator. So, for example, if
    you would like to make a galactic map projection centered on the galactic
    center, before passing galactic longitude l to the function you should
    first do:
    l = l if l <= math.pi else l - 2 * math.pi

    Keyword arguments:
    lon -- Traditional longitude in radians, in range [-pi:pi]
    lat -- Traditional latitude in radians, in range [-pi/2:pi/2]
    """

    # check if the input values are in the range
    if lon > math.pi or lon < -math.pi or lat > math.pi / 2 or lat < -math.pi /2 :
        sys.stdout.write('Mollweide: Input longitude and latitude out of range.\n')
        sys.stdout.write('           lon: [-pi,pi]; lat: [-pi/2,pi/2].\n')
        sys.exit()
        
    # bypass lat = +/- pi/2, otherwise division by zero may occur
    if lat == math.pi / 2.0:
        theta = math.pi / 2.0
    elif lat == -math.pi / 2.0:
        theta = -math.pi / 2.0
    else:
        # a simple Newton-Raphson iteration
        # to solve the equation 2x + sin(2x) = pi*sin(lat) 
        x1 = lat

        while 1:

            x2 = x1 - ((2.0 * x1 + math.sin(2.0 * x1) - math.pi * math.sin(lat)) /
                       (2.0 + 2.0 * math.cos(2.0 * x1)))
            
            # break the loop when desired accuracy achieved 
            if math.fabs(x2 - x1) < 1.0e-10: 
                break
            else:
                x1 = x2
        
        theta = x2

    x = 2.0 * math.sqrt(2.0) / math.pi * lon * math.cos(theta)

    y = math.sqrt(2.0) * math.sin(theta)
    
    return x,y
# #The end of mollweide ###############################################


# ######################################################################
# Aitoff Projection
# http://en.wikipedia.org/wiki/Aitoff_projection
# #######################################################################
def aitoff(lon, lat):
    """
    Make Aitoff map projection.
    
    Take traditional longitude and latitude in radians and return a
    tuple (x, y).
    
    Notice that traditionally longitude is in [-pi:pi] from the meridian,
    and latitude is in [-pi/2:pi/2] from the equator. So, for example, if
    you would like to make a galactic map projection centered on the galactic
    center, before passing galactic longitude l to the function you should
    first do:
    l = l if l <= math.pi else l - 2 * math.pi

    Keyword arguments:
    lon -- Traditional longitude in radians, in range [-pi:pi]
    lat -- Traditional latitude in radians, in range [-pi/2:pi/2]
    """

    # check if the input values are in the range
    if lon > math.pi or lon < -math.pi or lat > math.pi / 2 or lat < -math.pi /2 :
        sys.stdout.write('Aitoff: Input longitude and latitude out of range.\n')
        sys.stdout.write('           lon: [-pi,pi]; lat: [-pi/2,pi/2].\n')
        sys.exit()

    # take care of the sigularity at (0, 0), otherwise division by zero may happen
    if lon == 0 and lat ==0:
        return 0.0, 0.0
    
    # a quick inline unnormalized sinc function, with discontinuity removed 
    sinc = lambda x: 0 if x == 0 else math.sin(x) / x

    alpha = math.acos(math.cos(lat) * math.cos(lon / 2.0))

    # the sinc function used here is the unnormalized sinc function
    x = 2.0 * math.cos(lat) * math.sin(lon / 2.0) / sinc(alpha)

    y = math.sin(lat) / sinc(alpha)

    return x, y
# #The end of aitoff ##################################################


# ####################################################################
# Hammer Projection
# http://en.wikipedia.org/wiki/Hammer_projection
# ####################################################################
def hammer(lon, lat):
    """
    Make Hammer map projection.
    
    Take traditional longitude and latitude in radians and return a
    tuple (x, y).
    
    Notice that traditionally longitude is in [-pi:pi] from the meridian,
    and latitude is in [-pi/2:pi/2] from the equator. So, for example, if
    you would like to make a galactic map projection centered on the galactic
    center, before passing galactic longitude l to the function you should
    first do:
    l = l if l <= math.pi else l - 2 * math.pi

    Keyword arguments:
    lon -- Traditional longitude in radians, in range [-pi:pi]
    lat -- Traditional latitude in radians, in range [-pi/2:pi/2]
    """

    # check if the input values are in the range
    if lon > math.pi or lon < -math.pi or lat > math.pi / 2 or lat < -math.pi /2 :
        sys.stdout.write('Hammer: Input longitude and latitude out of range.\n')
        sys.stdout.write('           lon: [-pi,pi]; lat: [-pi/2,pi/2].\n')
        sys.exit()

    s = math.sqrt(1.0 + math.cos(lat) * math.cos(lon / 2.0))

    x = 2.0 * math.sqrt(2.0) * math.cos(lat) * math.sin(lon / 2.0) / s

    y = math.sqrt(2.0) * math.sin(lat) / s
    
    return x, y
# #The end of hammer #################################################

