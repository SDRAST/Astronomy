"""
Classes and functions for astronomical calculations

This is a complete re-write of the package to base it on astropy instead of 
heritage code.  The following were removed::

 * List Months,       with long month names
 * default_astro_dir, for source catalogs
 * cat_path,          for source catalogs
 * eph_url,           IERS site for UT1-UTC
 * short_month(),     which truncates the month name
 * long_month(),      extracts the month name from Months by number
 * parse_iers_line(), for earth orientation parameters
 * get_cartesian_coordinates(), now in MonitorControl.Configurations
 * show_catalogs()
 * get_sources()
 * get_coords()
 * geodetic_coords,   now in MonitorControl.Configurations
 

Basic Terminology
=================
Time
----
Terrestrial Time (TT) is a modern astronomical time standard defined by the 
IAU, primarily for time-measurements of astronomical observations made from 
the surface of Earth. TT continues Terrestrial Dynamical Time (TDT) which in
turn succeeded ephemeris time (ET). The purpose for which ET was designed is to
be free of the irregularities in the rotation of Earth.  The unit of TT is the 
SI second, the definition of which is currently based on the caesium atomic 
clock.

TT is distinct from UTC, the time scale used as a basis for civil 
purposes, Coordinated Universal Time (UTC). TT indirectly underlies UTC, via 
International Atomic Time (TAI). Each leap second that is introduced into UTC 
causes UTC to diverge a little further from TT.

Epoch
-----
An epoch is a moment in time used as a reference point for some time-varying 
astronomical quantity, such as the celestial coordinates.

Equinox
-------
In its broadest sense this is the moment when the day and night are exactly
equally long.  This defines the equator (and the ecliptic), latitude zero.
This moment does not occur at exactly the same time each year so a mean
equinox time is defined.

Besselian Year
--------------
The beginning of a Besselian year to be the moment at which the mean longitude
of the Sun, including the effect of aberration and measured from the mean 
equinox of the date, is exactly 280 degrees.  This moment falls near the 
beginning of the corresponding Gregorian (modern calendar) year.  A 
"Besselian epoch" can be calculated from the Julian date.  The beginnings of
some commonly used Besselian years are::

  B1900.0 = JD 2415020.3135 = 1900 January 0.8135 TT
  B1950.0 = JD 2433282.4235 = 1950 January 0.9235 TT

FK4
---
The Fourth Fundamental Catalogue is a reference frame published in 1963 based
on 1,535 stars in various equinoxes from 1950.0 to 1975.0.

FK5
---
The Fifth Fundamental Catalog is a reference frame published in 1988 with 
updated new positions for the 1,535 stars.

ICRS
----
The International Celestial Reference System (ICRS) is the current standard 
celestial reference system adopted by the IAU. Its origin is at the 
barycenter of the solar system, with axes that are intended to be fixed 
in space. The ICRS is based on the International Celestial Reference Frame
(ICRF) which consists of the positions of 3414 compact radio sources measured 
using VLBI. 

Celestial Coordinates
=====================

Brandon Rhodes gives an excellent explanation in the documentation for 
PyEphem.  Summarizing:

Astrometric Geocentric
  Mean geocentric position for the epoch of the specified star atlas
Apparent Geocentric
  current geocentric position at the date and time of observation
Apparent Topocentric
  current position for the observatory at the date and time of observation

Ecliptic Coordinates
====================
``astropy`` does not yet support ecliptic coordinates so we use ``pyephem`` for that.
``https://github.com/firelab/met_utils/blob/master/sun.py`` has some extensions for
``astropy`` for that and other things.

"""
import logging
logger = logging.getLogger(__name__)
import math
import datetime
import os.path

import astropy
from astropy.coordinates import AltAz, EarthLocation, SkyCoord  # High-level
from astropy.coordinates import ICRS, FK4, FK5            # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude          # Angles
from astropy.time import Time
import astropy.units as u
from astropy.units import cds
from astropy.coordinates.name_resolve import get_icrs_coordinates, \
                                             NameResolveError
import ephem
from ephem import Ecliptic, Equatorial

from Astronomy.coordconv import coordconv
from DatesTimes import calendar_date, MJD
from Math.geometry import Circular

from math import pi
c = cds.c.in_units('m / s')
pc = cds.pc.in_units('m')
AU = cds.AU.in_units('m')

R_earth = cds.Rgeo.in_units('m') # 6731e3 # m
R_sun = cds.Rsun.in_units('m')   # 695800e3 # m

# for backward compatibility
import DatesTimes as DT
day_of_week = DT.day_of_week
day_of_year = DT.day_of_year
MJD         = DT.MJD
julian_date = DT.julian_date
leap_year   = DT.leap_year

from . import DSN_coordinates as C
get_cartesian_coordinates = C.get_cartesian_coordinates
get_geodetic_coords = C.get_geodetic_coords

# planets recognized by module ephem
Planets = ['Jupiter', 'Mars', 'Mercury', 'Moon', 'Neptune', 'Pluto',
           'Saturn', 'Sun', 'Uranus', 'Venus']


def refresh_ierstab(url=None, force=False):
  """
  Updates the UT1-UTC offset table from 'url'.
  
  Deals with this warning from astropy:
  WARNING: Tried to get polar motions for times after IERS data is valid. 
  Defaulting to polar motion from the 50-yr mean for those.
  If you need enough precision such that this matters (~<10 arcsec), you can 
  download the latest IERS predictions.
  
  @param url : not used
  @param force : not used
  @return: None
  
  """
  from astropy.utils.data import download_file
  from astropy.utils import iers
  iers.IERS.iers_table = iers.IERS_A.open(download_file(iers.IERS_A_URL,
                                                          cache=True))

def B_epoch_to_J(ra50, dec50, format=None):
  """
  Convert B1950 coordinates to J2000
  
  This is a convenience to avoid having to format angles before calling
  the coordinate conversion routine.  The need often comes up when an
  operator asks for coordinates in J2000 or decimal format.

  Examples
  ========  
  Example
  -------
  Basic::
  
   In [1]: import Astronomy as A
   In [2]: A.B_epoch_to_J('''00h02m29.056400s''',
                          '''54d11'43.187000"''', 'decimal')
   Out[2]: (0.084544, 54.4736)
   In [3]: A.B_epoch_to_J('''00h02m29.056400s''',
                          '''54d11'43.187000"''', 'formatted')
   Out[3]: [u'00h05m04.359s', u'+54d28m25.056s']
   In [4]: A.B_epoch_to_J('''00h02m29.056400s''',
                          '''54d11'43.187000"''')
   Out[4]: ([0, 5, 4.359], [54, 28, 25.056])

  Compare to the VLA Calibrator List::
  
    J2000 A 00h05m04.363531s 54d28'24.926230"
    B1950 A 00h02m29.056400s 54d11'43.187000"
  
  Example
  -------
  Crossing the midnight boundary::
  
    In [8]: A.B_epoch_to_J('''23h58m34.865400s''',
                           '''18d57'51.753000"''')
    Out[8]: ([0, 1, 8.6169], [19, 14, 33.9321])
    
  Compare to the VLA Calibrator List::
  
    J2000 A 00h01m08.621563s 19d14'33.801860"
    B1950 A 23h58m34.865400s 18d57'51.753000"
  
  Example
  -------
  Negative declination::
  
    In [10]: A.B_epoch_to_J('''00h00m48.4200s''',
                            '''-17d43'54.000"''', 'formatted')
    Out[10]: [u'00h03m21.9921s', u'-17d27m11.6511s']


  @param ra50 : string
    For example: '00h02m29.056400s'

  @param dec50 : string
    For example : '''54d11'43.187000"'''

  @param format : string
    'decimal', 'formatted', None

  @return: tuple of strings
    See notes for details.
    
  """
  coordstr = ra50+" "+dec50
  logger.debug("B_epoch_to_J: 1950 coordinates: %s", coordstr)
  coords = SkyCoord(coordstr, frame="fk4", unit=(u.hourangle, u.deg))
  if format == None:
    rastr, decstr = coords.fk5.to_string('hmsdms').split()
    h = rastr.split('h')[0]
    m = rastr.split('h')[1].split('m')[0]
    s = rastr.split('h')[1].split('m')[1][:-1]
    ralist = [int(h),int(m),float(s)]
    d = decstr.split('d')[0]
    m = decstr.split('d')[1].split("m")[0]
    s = decstr.split('d')[1].split("m")[1][:-1]
    declist = [int(d),int(m),float(s)]
    return ralist,declist
  if ( format[0] == "d" ):
    rastr, decstr = coords.fk5.to_string().split()
    return float(rastr)/15, float(decstr)
  elif ( format[0] == "f"):
    return coords.fk5.to_string('hmsdms').split()

def J_epoch_to_B(ra2000, dec2000, format):
  """
  Convert J2000 coordinates to B1950.
  
  See B_epoch_to_J for documentation

  Notes
  =====
  A test::
  
    In [1]: import Astronomy as A
    In [2]: A.B_epoch_to_J('''00h02m29.056400s''',
                           '''54d11'43.187000"''', 'formatted')
    Out[2]: [u'00h05m04.359s', u'+54d28m25.056s']
    In [4]: A.J_epoch_to_B('''00h05m04.359s''',
                           '''+54d28\'25.056"''', 'formatted')
    Out[4]: [u'00h02m29.0564s', u'+54d11m43.1869s']

  Seems to be good to about an arcsecond.

  @param ra2000 : string

  @param dec2000 : string

  @param format : string

  @return: tuple of strings
  """
  coordstr = ra2000+" "+dec2000
  coords = SkyCoord(coordstr, frame="fk5")
  if format == None:
    rastr, decstr = coords.fk4.to_string('hmsdms').split()
    h = rastr.split('h')[0]
    m = rastr.split('h')[1].split('m')[0]
    s = rastr.split('h')[1].split('m')[1][:-1]
    ralist = [int(h),int(m),float(s)]
    d = decstr.split('d')[0]
    m = decstr.split('d')[1].split("m")[0]
    s = decstr.split('d')[1].split("m")[1][:-1]
    declist = [int(d),int(m),float(s)]
    return ralist,declist
  if ( format[0] == "d" ):
    rastr, decstr = coords.fk4.to_string().split()
    return float(rastr)/15, float(decstr)
  elif ( format[0] == "f"):
    return coords.fk4.to_string('hmsdms').split()

def YMD_to_MJD(year,month,day):
  """
  Modified Julian date from calendar date

  @param year : four digit year
  @type  year : int

  @param month : int

  @param day : int

  @return: long
  """
  return MJD(year,month,day)

def v_sun(mjd, ra_source, dec_source):
  """
  LSR velocity of the sun

  Notes
  =====
  The velocity of the Sun with respect to the Local Standard of Rest is
  the conventional value of 20.0 km/s toward
  RA, dec = 18h03m50.29s, +30d00'16.8 (J2000)

  @param mjd : int
    mean Julian day

  @param ra_source : float
    Right ascension in decimal hours

  @param dec_source : float
    Declination in decimal degrees

  @return: float
    the projection onto the line of sight of the velocity of the
    Sun with respect to the Local Standard of Rest on the date given by
    mjd at UT = 0h.
  """
  LSR = SkyCoord('''18h03m50.29s +30d00'16.8"''')
  # Precess solar motion from J2000 to epoch
  mjd = Time(mjd, format='mjd')
  LSRmjd = LSR.transform_to(FK5(equinox=mjd))
  ra_epoch = LSRmjd.ra.radian
  dec_epoch = LSRmjd.dec.radian
  # Cartesian coordinates
  v_x = 20.0 * math.cos(ra_epoch) * math.cos(dec_epoch)
  v_y = 20.0 * math.sin(ra_epoch) * math.cos(dec_epoch)
  v_z = 20.0 * math.sin(dec_epoch)
  # Source direction in Cartesian coordinates
  ra = ra_source*pi/12
  dec = dec_source*pi/180
  x = math.cos(dec) * math.cos(ra)
  y = math.cos(dec) * math.sin(ra)
  z = math.sin(dec)
  # Dot product
  v_sun = -v_x * x - v_y * y - v_z * z
  return v_sun

def decimal_day_to_tuple(day):
  """
  Float day to (h,m,s) tuple

  @param day : float

  @return: tuple (int, int, float)
  """
  t = Longitude(day, unit='cycle')
  return (int(t.hms.h), int(t.hms.m), t.hms.s)

def decimal_day_to_HMS(day):
  """
  Formatted time to 0.01 seconds

  @param day : float

  @return: string
    'HH:MM:SS.SS"
  """
  
  return "%02d:%02d:%05.2f" % decimal_day_to_tuple(day)

def greenwich_sidereal_time(year,doy):
  """
  Approximate sidereal time at Greenwich

  @param year : year of date
  @type  year : int

  @param doy : day of year
  @type  doy : float

  @return: Greenwich sidereal time in hours
  """
  year_from_1966 = year-1966
  dt = (year_from_1966*365 + int((year_from_1966 + 1)/4.) + int(doy)-1)/36525.
  dst = 0.278329562 + (8640184.67*dt+0.0929*dt**2)/86400
  gst0 = dst % 1 # GST on Jan. 0 of current year
  return 24*(gst0 + (doy % 1)/0.997269566) % 24

def time_aliases(year, UTdoy, obs_long):
  """
  Time as days since 1900, centuries since 1900 and LST

  Notes
  =====
  This simplifies the output from ManGord.atime().

  @param year : int

  @param UTdoy : float
    Day of year and niversal time

  @param obs_long : float
    Observatory west longitude in degrees

  @return: tuple of floats
    (number of days since 1900.5, the same
    in Julian centuries, and the local sidereal time in decimal days
  """
  doy = int(UTdoy)
  date_tuple = calendar_date(year,doy)
  h,m,s = decimal_day_to_tuple(UTdoy-doy)
  dt_tuple = date_tuple+(h,m,int(s),int((s-int(s))*1e6))
  logger.debug("dt_tuple: %s", dt_tuple)
  time = datetime.datetime( *dt_tuple )
  t = Time(time)
  days_since_1900 = t.mjd - MJD(1900,1,1) + 1
  try:
    lst = t.sidereal_time('mean',longitude=-obs_long*u.deg)
  except IndexError:
    logger.warning(" Times is outside of range covered by IERS table.")
    t.delta_ut1_utc = 0.
    lst = t.sidereal_time('mean', longitude = -obs_long*u.deg)
  julian_centuries_since_1900 = days_since_1900/36525.
  return days_since_1900, julian_centuries_since_1900, lst.cycle

def current_to_ecliptic(julian_centuries_since_1900,\
                        current_ra,\
                        current_dec):
  """
  Apparent ecliptic coordinates from apparent RA and dec

  Notes
  =====
  Checked against http://lambda.gsfc.nasa.gov/toolbox/tb_coordconv.cfm

  @param julian_centuries_since_1900 : float

  @param current_ra : float
    radians

  @param current_dec : float
    radians

  @return: tuple of floats
    (eclip_long, eclip_lat) in radians
  """
  t = Time(julian_centuries_since_1900*36525 + MJD(1900,1,1) - 1, format='mjd')
  sc = SkyCoord(ra=current_ra, dec=current_dec, unit='radian',
                frame='icrs')
  ra, dec = sc.to_string('hmsdms').split()
  ra = ra.replace('h',':').replace('m',':')[:-1]
  dec = dec.replace('d',':').replace('m',':')[:-1]
  eq = Equatorial(str(ra), str(dec), epoch=1950)
  eclip = Ecliptic(eq)
  return eclip.long, eclip.lat

def ecliptic_to_J2000(elong, elat, mjd):
  """
  convert from ecliptic coordinates to RA and dec
  """
  t = Time(mjd, format='mjd')
  # t is probably needed for convert from current to J2000
  eq = Ecliptic(str(elong), str(elat), epoch=2000)
  ra, dec = Equatorial(eq).to_radec()
  #print "Astronomy:",type(ra),type(dec)
  #print "Astronomy:", type(ra.real), type(dec.real)
  return ra.real, dec.real
  
def HaDec_to_AzEl(HourAngle, Declination, Latitude):
  """
  Celestial to horizon coordinates
  
  Notes
  =====
  Sidereal Time verified with http://tycho.usno.navy.mil/sidereal.html

  @param HourAngle : float
    decimal hours

  @param Declination : float
    decimal degrees

  @param Latitude : float
    decimal degrees

  @return: tuple
    Azimuth and elevation in degrees.
  """
  t = Time.now()
  l = EarthLocation(lon=0*u.deg, lat=Latitude*u.deg)
  t.location = l
  t.delta_ut1_utc = 0
  lst = t.sidereal_time('mean')
  c = SkyCoord(ra=(lst.hour-HourAngle)*u.hourangle, dec=Declination*u.deg,
               obstime=t)
  c.location = l
  return c.altaz.az.deg, c.altaz.alt.deg
  
def AzEl_to_HaDec(Azimuth, Elevation, Latitude):
  """
  converts from azimuth and elevation to hour angle and declination
  
  @param Azimuth : degrees, clockwise from north
  @type  Azimuth : float

  @param Elevation : degrees above horizon
  @type  Elevation : float
  

  @param Latitude : east longitude in degrees
  @type  Latitude : float

  @return: tuple, Hour angle in float hours and declination in float degrees
  """
  azr = Azimuth*pi/180
  elr = Elevation*pi/180
  latr = Latitude*pi/180
  har,decr = coordconv(pi, pi/2-latr, 0, latr, azr, elr)
  return har*12/pi, decr*180/pi

def AzEl_to_RaDec(azimuth,elevation,latitude,longitude,date_time):
  """
  Convert azimuth and elevation to right ascension and declination

  @param azimuth : east from north (clockwise) in degrees
  @type  azimuth : float

  @param elevation : above the horizon in degrees
  @type  elevation :

  @param latitude : above the equator, in degrees
  @type  latitude : float

  @param longitude : west from Greenwich
  @type  longitude : float

  @param date_time : (year, DOY) tuple, with fractional day of year
  @type  date_time : (int, float)

  @return: (RA (hrs), dec (degs))
  
  Notes
  =====
  Horizon coordinates to celestial
  
  HA and decl. define a point on the sky with respect to the local meridian,
  which corresponds to a RA equal to the local sidereal time.  It doesn't
  matter what the observer's longitude is. The relationship between the RA of
  the sky point and the LST stays the same. So we can perform the calculation
  for longitude 0.
  
  LST is the ST at Greenwich (long 0 deg) minus the west longitude of the 
  local meridian::
    
    LST = GST - long
                    hr
                    
  HA is positive to the west, so is the LST minus the RA of the point in the
  sky::
  
    HA = LST - RA

  Consider a hypothetical observer on the Earth at longitude zero and the
  latitude of the actual observer. The sidereal time at the hypothetical
  observer's location is the Greenwich sidereal time, which is the actual
  observer's LST plus the west longitude. Then the HA at the actual observer's
  position is the HA at the hypothetical observer's position minus the actual
  observer's longitude in hours::
  
    HA    = HA    - long
      act     hyp       hr
  
  Given the azimuth, elevation and time we compute the RA and dec of the sky
  position with respect to the hypothetical observer, and then the HA::
  
    az, el, lat, time --> RA   ,dec
                            hyp
  
    HA    = GST - RA
      hyp           hyp
  
  The RA of the same az,el w.r.t. the actual observer is then::
  
    RA    = RA    - long
      act     hyp       hr
  
    HA    = GST - RA    - long  
      act           hyp       hr
          = GST - RA
                    act
  
  """
  LST = greenwich_sidereal_time(*date_time)-longitude/15.
  HA,dec  = AzEl_to_HaDec(azimuth, elevation, latitude)
  RA = math.fmod(LST - HA, 24.)
  if RA < 0:
    RA += 24.
  return RA,dec
  
def J2000_to_apparent(MJD, UT, ra2000, dec2000):
  """
  Apparent right ascension and declination

  @param MJD : int
    mean Julian day

  @param UT : float
    decimal hours

  @param ra2000 : float
    radians

  @param dec2000 : float
    radians

  @return: tuple
  """
  t = Time(MJD+UT/24., format='mjd')
  coords = SkyCoord(ra=ra2000*u.rad, dec=dec2000*u.rad, frame='icrs')
  c = coords.transform_to(FK5(equinox=t))
  return c.ra.hourangle, c.dec.deg
  
def get_sky_coords(RA, Dec):
  """
  return SkyCoord from both float and str inputs
  """
  if type(RA) == str and type(Dec) == str:
    skypos = SkyCoord(RA, Dec, unit=(u.hourangle, u.deg))
  elif type(RA) == float and type(Dec) == float:
    skypos = SkyCoord(RA*u.hour, Dec*u.degree)
  else:
    raise RuntimeError(RA, dec, "cannot be parsed")
  return skypos

def get_altaz(RA, Dec, time, location):
  """
  Simplest conversion using astropy
  
  @param RA : right ascension in degrees
  @type  RA : float
  
  @param Dec : declination in degrees
  @type  Dec : float
  
  @param time : time of observation
  @type  time : datetime.datetime object
  
  @param location : location of the observatory
  @type  location : astropy.Location object
  """
  skypos = get_sky_coords(RA, Dec)
  #logger.debug("get_altaz: called for RA,dec: %s", skypos)
  skypos.obstime = Time(time)
  skypos.location = location
  altaz = skypos.altaz.az.deg, skypos.altaz.alt.deg
  #logger.debug("get_altaz: az,el: %s", altaz)
  return altaz

def object_az_el(source, site, year, doy):
  """
  Compute object's position in alt-az for a given site and time
  
  Also returns the source's apparent coordinates
  
  @param source : a source name recognized by Simbad
  @type  source : str
  
  @param site : DSN station number
  @type  site : int
  
  @param year : four digit year
  @type  year : int
  
  @param doy : DOY with UT as a fraction of a day
  @type  doy : float
  """
  try:
    coords = get_icrs_coordinates(source)
  except NameResolveError as details:
    raise NameResolveError(details)
  module_logger.debug("Sky coords: %s", coords)
  
  try:
    dss = C.DSS(site)
    module_logger.debug("DSS-%d: %f, %f", site, dss.long*180/pi, dss.lat*180/pi)
  except KeyError:
    raise KeyError('%d is not a valid DSS station' % site)
  loc = EarthLocation(dss.long*u.rad, dss.lat*u.rad)
  module_logger.debug("Site coords: %s", loc)
  
  if doy:
    mjd = MJD(year,doy)
  else:
    raise RuntimeError("no DOY given")
  tt = Time(mjd, format='mjd')
  module_logger.debug("ISO time = %s", tt.iso)
  tt.delta_ut1_utc = 0
  coords.obstime = tt
  coords.location = loc
  return coords.altaz
  
def ha_rise_set(el_limit, lat, dec):
  """
  Hour angle from transit for rising and setting.

  Returns pi for a source that never sets and 0 for a source always below
  the horizon.

  @param el_limit : the elevation limit in radians
  @type  el_limit : float

  @param lat : the observatory latitude in radians
  @type lat : float

  @param dec : the source declination in radians
  @type  dec : float

  @return: hour angle from transit in radians
  """
  cos_ha = (math.sin(el_limit) - math.sin(lat)*math.sin(dec)) \
           /(math.cos(lat)*math.cos(dec))
  if cos_ha <= -1:
    # never sets
    return pi
  elif cos_ha >= 1:
    # never visible
    return 0
  else:
    return math.acos(cos_ha)

def great_circle_angle(lat1,long1,lat2,long2):
  """
  Great circle angle between two points on a sphere

  Given two positions on a sphere, return the angle between them
  in degrees.
  
  @param lat1 : initial latitude in degrees
  @type  lat1 : float 

  @param long1 : initial longitude in degrees
  @type  long1 : float

  @param lat2 : final latitude in degrees
  @type  lat2 : float

  @param long2 : final longitude in degrees
  @type  long2 : float 

  @return: float
  """
  dlongr = (long1 - long2)*pi/180.
  lat1r = lat1*pi/180.
  long1r = long1*pi/180.
  lat2r = lat2*pi/180.
  long2r = long2*pi/180.
  term1 = math.cos(lat2r)*math.sin(dlongr)
  term2 = math.cos(lat1r)*math.sin(lat2r)
  term3 = math.sin(lat1r)*math.cos(lat2r)*math.cos(dlongr)
  numer = math.sqrt(math.pow(term1,2)+math.pow(term2-term3,2))
  term4 = math.sin(lat1r)*math.sin(lat2r)
  term5 = math.cos(lat1r)*math.cos(lat2r)*math.cos(dlongr)
  denom = term4+term5
  cangler = math.atan2(numer,denom)
  return cangler*180/pi

def great_circle_distance(central_angle):
  """
  Great circle distance in km on Earth from great circle angle

  Notes
  =====
  Assumes a spherical earth of radius 6372.795 km.
  This is accurate to about 0.5%

  @param central_angle : angle subtended at the center of the earth
  @type  central_angle : float

  @return: float
  """
  return 6372.795*pi*central_angle/180.

def horizon_range(planet_radius, observer_height):
  """
  Range to horizon

  Both arguments in the same units and range will be in those units too.
  """
  circle = Circular(planet_radius)
  #cos_theta = float(planet_radius)/(planet_radius + observer_height)
  #theta = math.acos(cos_theta)
  #return planet_radius*math.sin(theta)
  return circle.tangent_distance(observer_height)

def source_start_stop(HAstart, HAend, HAlimit):
  """
  Find when sources are up
  
  This only handles tracks less than 24 hr long.
  
  The result codes are::
  
    A  - always up
    N  - never up
    R  - rises
    RS - rises and sets
    S  - sets
    SR - sets and rises (up twive during a track
    
  @param HAstart - HA at beginning of track
  @type  HAstart - float or int
  
  @param HAend - HA at end of track
  @type  HAend - float or int
  
  @return: str
  """
  if abs(HAstart) < HAlimit and abs(HAend) < HAlimit:
    if HAend > HAstart:
      result = "A"
    else:
      result = "SR"
  elif abs(HAstart) > HAlimit and abs(HAend)> HAlimit:
    if math.copysign(1,HAstart) == math.copysign(1,HAend):
      if (HAend > HAstart):
        result = "N"
      else:
        result = "RS"
    else:
      if (HAend > HAstart):
        result = "RS"
      else:
        result = "N"
  elif ((HAstart < 0 and HAstart < -HAlimit) or
        (HAstart > 0 and HAstart >  HAlimit)):
    result = "R"
  elif (((HAstart > 0 and HAstart <  HAlimit) or
         (HAstart < 0 and HAstart > -HAlimit)) and
        ((HAend   > 0 and HAend   >  HAlimit) or
         (HAend   < 0 and HAend   < -HAlimit))):
    result = "S"
  elif (((HAstart < 0 and HAstart < -HAlimit) or
         (HAstart > 0 and HAstart >  HAlimit)) and
        ((HAend   < 0 and HAend   < -HAlimit) or
         (HAend   > 0 and HAend > HAlimit))):
    result = "SR"
  else:
    result = None
  return result

def galactic_offsets_to_celestial(RA, Dec, glongoff=3, glatoff=0):
  """
  Converts offsets in Galactic coordinates to celestial
  
  The defaults were chosen by Pineda for Galactic plane survey
  
  @param RA : FK5 right ascension in degrees
  @type  RA : float
  
  @param Dec : FK5 declination in degrees
  @type  Dec : float
  
  @param glongoff : Galactic longitude (degrees)
  @type  glongoff : float
  
  @param glatoff : Galactic latitude offset (degrees)
  @type  glatoff : float
  """
  # Initialize SkyCoord Object
  gc = SkyCoord(ra=RA*u.degree, dec=Dec*u.degree, frame='fk5') 
  # Add offset in galactic longitude
  glongOFF=gc.galactic.b.value + glongoff
  # Add offset in galactic latitude
  glat=gc.galactic.l.value + glatoff
  # Convert reference position to RA/DEC
  radecoff=SkyCoord(l=glat*u.degree, b=glongOFF*u.degree, frame='galactic')
  # Replace Raoff and Decoff with new values 
  RAoff = radecoff.fk5.ra.value - RA
  Decoff = radecoff.fk5.dec.value - Dec     
  return RAoff, Decoff
