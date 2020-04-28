# -*- coding: utf-8 -*-
"""
Functions for computing solar orbital and orientation parameters

The ecliptic, the apparent solar orbit (in reality, the orbit of Earth about
the Sun), is tilted approximately 23.5 deg with respect to the celestial
equator. The sun is tilted at an angle of 7.25 deg with respect to the
ecliptic plane. Both the celestial coordinates and the ecliptic coordinates
have their zero meridian (RA = 0 h, lambda = 0 deg) at the point where the
Sun's orbit rises above the celestial plane.

Contents
========

Utilities
---------
Provided for convenience::

  radians(angle)  - Convert degrees to radians
  truncate(angle) - Bounds an angle between 0 and 360.

Functions Related to Time
-------------------------
Converts time to other formats::

  Carrington_Number(jd)        - Carrington Number gives the number of rotations.
  calcJD(y,m,d,t)              - Fractional Julian date
  Julian_centuries_since_2000(jd) - Julian centuries from Jan 1, 2000.
  time_series(jd)                 - List of powers of Julian centuries.

Functions Related to the Ecliptic
---------------------------------
T is a tuple with powers (0,1,2,3) of Julian centuries::

  ecliptic_obliquity(T) - Inclination of ecliptic plane w.r.t. celestial equator
  Eccentricity_of_orbit(T)- Eccentricity of the orbit ellipse.
  sun_geom_mean_long(T)   - Geometric Mean Ecliptic Longitude (deg) of Sun
  sun_mean_anomaly(T)     - Position of Sun relative to perigee
  equation_of_the_center(M,T) - Difference between mean anomaly and true anomaly
  Sun_radius_vector(e, tar)   - Distance between the Sun and the Earth
  semidiameter(distance)      - Semi-diameter in arcsec

Functions concerning Nutation
-----------------------------
Handle nutation of the Earth's orbit::

  lunar_long_asc_node(t) - longitude of the mean ascending node of the lunar
                           orbit on the ecliptic; t in Julian centuries
  longitude_nutation(omega) - Nutation in longitude

Functions for Heliographic Coordinates
--------------------------------------
These provide coordinates in the Sun-based lat/long system::

  sun_coords_long_ascen_node(jd) - Ecliptic longitude at which solar equator
                                   intersects the ecliptic plane
  Sun_central_longitude(anomaly, inclin, Lsun) - Longitude of center of disk
  Sun_central_latitude(longitude, tilt) - Sun's latitude where central meridian
                                          crosses ecliptic plane
  lat_long(radius, polar_position_angle, radial_distance, position_angle)
                                    - latitude and longitude of an active region

Geometrical Functions
---------------------
manage projections to/from plane of the sky ::

  axis_tilt_projection(longitude,tilt) - Projection of Z-axis of a tilted frame
  xy_to_ra_dec(x,y,P,ra0,dec0)    - Convert XY coordinates in a Sun-aligned
                                    frame to RA and dec
  ra_dec_to_xy(ra,dec,P,ra0,dec0) - Convert ra, dec to xy-coordinates aligned
                                    with Sun

Comprehensive Parameter Calculation
-----------------------------------
Handles a lot of calculations at once::

  calc_solar(jd) - Overall function to calculate various solar data.

Notes
=====
Times in this module for slowly varying parameters are in Julian centuries from
Jan. 1, 2000.  This differs from times in module Astronomy, which are based on
Jan. 1, 1900.

Resources
=========
A form for calculating sunspot coordinates::

  http://www.nature1st.net/bogan/astro/sun/sunspots.html

Formulae for the Solar Orientation angles were taken from
Astronomical Algorithms by Meeus.
  
Formulae for Calculating the Sunspot lat-long taken from Smith's Astronomical 
Calculation for Calculators from::

  http://www.idialstars.com/fipl.htm

Supporting definitions can be found at::

  http://www.astro.washington.edu/docs/idl/cgi-bin/getpro/library32.html?GET_SUN
  http://sspg1.bnsc.rl.ac.uk/SEG/Coordinates/angles.htm
  http://farside.ph.utexas.edu/syntaxis/Almagest/node34.html
  
"""
from Astronomy import julian_date
import math

radian = 180./math.pi
sidereal_rotation_period = 25.38 # days

# -----------------------------utilities-------------------------------
def radians(degrees):
  """
  Convert degrees to radians

  @type degrees : float

  @return: float
  """
  return math.pi/180

def truncate(angle):
  """
  Bounds an angle between 0 and 360.

  @type angle : float

  @return: float
  """
  truncated = fmod(angle,360.)
  if truncated < 0.:
    truncated += 360.
  return truncated

# -------------------------functions related to time ----------------------------

def Carrington_Number(jd):
  """
  Carrington Number gives the number of rotations.
  
  Carrington #1 was on November 8, 1853

  @param jd : Julian date including UT
  @type  jd : float
  """
  return floor((jd - 2398140.22710)/27.2752316)

def calcJD(y,m,d,t):
  """
  Julian date with fraction for UT

  @param y : for digit year
  @type  y : int

  @param m : month number
  @type  m : int

  @param d : day number
  @type  d : int

  @param t : UT, hrs
  @type  t : float

  @return: float
  """
  return julian_date(y,day_of_year(y,m,d)) + t/24.

def Julian_centuries_since_2000(jd):
  """
  Julian centuries from Jan 1, 2000.

  @param jd : Julian date
  @type  jd : float

  @return: float
  """
  t = (jd-2451545.)/36525.
  return t

def time_series(jd):
  """
  List of powers of Julian centuries.

  To have the ndex match the power, T[0] = 1

  @param jd : Julian date
  @type  jd : float

  @return: list of float
  """
  t = Julian_centuries_since_2000(jd)
  t2=t*t
  t3=t*t2
  return 1.0, t, t2, t3

#---------------Functions Related to the Ecliptic---------------------------

def ecliptic_obliquity(T):
  """
  Inclination of ecliptic plane w.r.t. celestial equator

  @param T : powers of Julian centuries
  @type  T : list of float

  @return: float
  """
  return 23.4392911-0.0130042*T[1]-0.0000164*T[2]+0.0000504*T[3]

def Eccentricity_of_orbit(T):
  """
  Eccentricity of the orbit ellipse.

  Eccentricity e is the ratio of half the distance between the foci c to
  the semi-major axis a: e=c/a. For example, an orbit with e=0 is circular,
  e=1 is parabolic, and e between 0 and 1 is elliptic.

  The eccentricity of Earth's orbit changes slowly with time.

  @param T : list of powers of Julian centuries from 2000.0
  @type  T : list of float

  @return: degrees
  """
  return 0.01675104 - 0.0000418*T[1] - 0.000000126*T[2]

def sun_geom_mean_long(T):
  """
  Geometric Mean Ecliptic Longitude (deg) of Sun

  @param T : list of powers of Julian centuries
  @type  T : list of float

  @return: float
  """
  return 280.46645 + 36000.76983*T[1] + 0.0003032*T[2]

def sun_mean_anomaly(T):
  """
  Position of Sun relative to perigee

  @param T : list of powers of Julian centuries
  @type  T : list of float

  @return: float
  """
  return 357.52910 + 35999.05030*T[1] - 0.0001559*T[2] - 0.00000048*T[3]

def equation_of_the_center(M,T):
  """
  Difference between mean anomaly and true anomaly

  @param M : Sun's mean anomaly
  @type  M : float

  @param T : list of powers of Julian centuries
  @type  T : list of float

  @return: float
  """
  Mr = M/radian
  return   (1.914600 - 0.004817*T[1] - 0.000014*T[2])*math.sin(Mr)   \
         + (0.019993 - 0.000101*T[1])                *math.sin(2*Mr) \
         +  0.000290                                 *math.sin(3*Mr)

def Sun_radius_vector(e, tar):
  """
  Distance between the Sun and the Earth

  There are a set of higher accuracy terms not included here.

  @param e : eccentricity of the orbit
  @type  e : float

  @param tar : Sun's true anomaly (radians)
  @type  tar : float

  @return: AU
  """
  return 1.0000002*(1. - e**2)/(1. + e*math.cos(tar))

def semidiameter(distance):
  """
  Semi-diameter in arcsec

  @param distance : distance between Sun and Earth, in AU
  @type  distance : float

  @return: float
  """
  return 959.63/distance

#--------------------Functions concerning Nutation-----------------------

def lunar_long_asc_node(t):
  """
  longitude of the mean ascending node of the lunar orbit on the ecliptic

  This is used to calculate nutation.

  @param t : Julian centuries
  @type  t : float

  @return: float
  """
  return 125.04-1934.136*t

def longitude_nutation(omega):
  """
  Nutation in longitude

  @param omega : longitude of Moon's ascending node
  @type  omega : float

  @return: float
  """
  return -0.00569 - 0.00478*math.sin(omega/radian)

#----------------------Functions for Heliographic Coordinates----------------

def sun_coords_long_ascen_node(jd):
  """
  Ecliptic longitude at which solar equator intersects the ecliptic plane

  For higher longitudes the equator is above the ecliptic plane, until
  the ecliptic longitude is 180 deg higher.

  @param jd : Julian date including UT
  @type  jd : float
  
  @return: deg
  """
  return 73.6667 + 1.3958333*(jd-2396758)/36525

def Sun_central_longitude(anomaly,inclin,Lsun):
  """
  Longitude of center of disk

  @param anomaly : position of Sun from perigee along its orbit
  @type  anomaly : float

  @param inclin : inclination of the orbit
  @type  inclin : float

  @param Lsun : true geometric ecliptic longitude (deg)
  @type  Lsun : float

  @return: longitude in degrees
  """
  etay = -math.sin(anomaly)*math.cos(inclin)
  etax = -math.cos(anomaly)
  eta = (math.atan2(etay,etax))*radian; # longitude of the Sun's central meridian
  L0 = (eta - Lsun)
  L0 = fmod(L0,360.)
  if L0 < 0.:
    L0 += 360.
  return L0

def Sun_central_latitude(longitude, tilt):
  """
  Sun's latitude where central meridian crosses ecliptic plane

  @param longitude : ecliptic longitude of the Sun
  @type  longitude : float

  @param tilt : tilt of Sun's axis w.r.t. celestial north pole
  @type  tilt : float

  @return: float
  """
  B0r = math.asin( math.sin(longitude)*math.sin(tilt));   # central latitude
  return B0r*radian;

def lat_long(radius,polar_position_angle,radial_distance,position_angle):
  """
  latitude and longitude of an active region

  @param radius : radius of Sun

  @param polar_position_angle : apparent tilt of polar axis

  @param radial_distance :

  @param position_angle :

  @return: (latitude, longitude)
  """
  Pr = polar_position_angle*radian
  position_angle = position_angle/radian
  r2 = math.asin(radial_distance/radius)
  dp = (Pr - position_angle)
  B = math.asin(math.sin(B0r)*math.cos(r2) + math.cos(B0r)*math.sin(r2)*math.cos(dp))
  L = math.asin(math.sin(r2)*math.sin(dp)/math.cos(B)) + L0r
  L = radian * L
  B = radian * B
  return B, L

# ---------------------- geometrical functions ---------------------------

def axis_tilt_projection(longitude,tilt):
  """
  Project of Z-axis of a tilted frame

  The tilt is in the YZ plane towards -Y

  @param longitude : longitude in XY plane from the X-axis (Y = 0)
  @type  longitude : float

  @param tilt : tilt angle measured CCW about X-axis
  @type  tilt : flat

  @return: radians
  """
  tan_projected_tilt = - math.cos(longitude)*math.tan(tilt)
  projected_tilt = math.atan(tan_projected_tilt)
  return projected_tilt

def xy_to_ra_dec(x,y,P,ra0,dec0):
  """
  Convert XY coordinates in a Sun-aligned frame to RA and dec

  @param x : degrees perpendicular to solar axis
  @type  x : float

  @param y : degrees parallel to the solar axis
  @type  y : float

  @param P : polar angle of solar axis from celestial north
  @type  P : float

  @param ra0 : right ascension of Sun's center in hours
  @type  ra0 : float

  @param dec0 : declination of Sun's center in degrees
  @type  dec0 : float

  @return: [ra,dec] corresponding to x and y
  """
  Pr = P/radian
  dec0r = dec0/radian
  ra  = ra0  + (x*math.cos(Pr) + y*math.sin(Pr))/math.cos(dec0r)/15
  dec = dec0 -  x*math.sin(Pr) + y*math.cos(Pr)
  return ra,dec

def ra_dec_to_xy(ra, dec, P, ra0, dec0):
  """
  Convert ra, dec to xy-coordinates aligned with Sun

  @param ra : right ascension, hrs
  @type  ra : float

  @param dec : declination, degrees
  @type  dec : float

  @param P : polar angle of solar axis from celestial north
  @type  P : float

  @param ra0 : right ascension of Sun's center in hours
  @type  ra0 : float

  @param dec0 : declination of Sun's center in degrees
  @type  dec0 : float

  @return: [x,y] in degrees corresponding to ra and dec
  """
  Pr = P/radian
  dec0r = dec0/radian
  x = -15*(ra-ra0)*math.cos(dec0r)*math.cos(Pr) + (dec-dec0)*math.sin(Pr)
  y =  15*(ra-ra0)*math.cos(dec0r)*math.sin(Pr) + (dec-dec0)*math.cos(Pr)
  return x,y

def calc_solar(jd):
  """
  Overall function to calculate various solar data.

  @param jd : fractional Julian date
  @type  jd : float

  @return: (Long_of_disk_center, Carrington_No, central_latitude, Pos_angle)
  """
  T = time_series(jd)

  # Geometric Mean Ecliptic Longitude (deg) of Sun
  L0 = sun_geom_mean_long(T)

  # Where the Sun is in its orbit determines how its rotation axis appears
  # on Earth
  # Mean anomaly (deg) measured around the orbit with 0 at perigee
  # and 180 at apogee.
  M = sun_mean_anomaly(T)
  # Sun's equation of center (deg) = Sun's ecliptic longitude
  # It is the way the true anomaly varies about the mean anomaly
  C = equation_of_the_center(M,T)
  # Sun's true geometric ecliptic longitude (deg)
  sunL = L0 + C
  # Sun's true anomaly (deg)
  v = M + C
  long_perigee = truncate(sunL - v)

  # Apparent ecliptic longitude (deg) from true ecliptic longitude
  omega = lunar_long_asc_node(T[1])
  longitude = fmod(sunL + longitude_nutation(omega),360.)

  ecc = Eccentricity_of_orbit(T)
  distance = Sun_radius_vector(ecc, v/radian)
  radius = semidiameter(distance)

  # Ecliptic longitude of ascending node of the Sun's equator
  long_asc_node = (sun_coords_long_ascen_node(jd))

  relative_long_asc_node = fmod((longitude-long_asc_node),360.) 

  # True Obliquity of the ecliptic (deg)
  obliq = ecliptic_obliquity(T)

  sun_inclination = 7.25
  inc = sun_inclination/radian # Inclination of the Sun's equator

  x = axis_tilt_projection(longitude/radian,obliq/radian)
  y = axis_tilt_projection(relative_long_asc_node/radian,inc)

  Pr = (x + y)     # Position angle
  P = Pr*radian    # in degrees

  B0 = Sun_central_latitude(relative_long_asc_node/radian, inc)

  jd0 = 2398220.0 # Julian date when Sun's longitude was 0 degrees
  Lsun = (jd-jd0)*360/sidereal_rotation_period
  L0 = Sun_central_longitude(relative_long_asc_node/radian, inc, Lsun)

  CarrNo = Carrington_Number(jd)
  result = {}
  result["CarrNo"]           = CarrNo
  result["L0"]               = L0
  result["B0"]               = B0
  result["ecliptic tilt"]    = x*radian
  result["solar tilt"]       = y*radian
  result["polar angle"]      = P
  result["obliquity"]        = obliq
  result["eclip long"]       = longitude
  result["sun incl"]         = sun_inclination
  result["sun perigee long"] = long_perigee
  result["anomaly"]          = fmod(v,360)
  result["sun asc node"]     = long_asc_node
  result["sun inclin"]       = sun_inclination
  result["eccentricity"]     = ecc
  result["distance"]         = distance
  result["semidiameter"]     = radius/60.
  return result

if __name__ == "__main__":
  DOY = day_of_year(2012,5,21)
  jd = julian_date(2012,DOY)+(1./24)
  print(calc_solar(jd))
  print("Bogan's calculator gives:")
  print("2123, 57.888144736760296, -1.93752, -19.0764")
