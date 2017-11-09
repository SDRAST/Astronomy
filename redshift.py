# -*- coding: utf-8 -*-
"""
module Astronomy.redshift - functions for computing to/from redshift

The redshift is the fractional amount by which a wavelength is shifted
longward (to the red) due to the recessional velocity of the emitter::
      lambda - lambda_0
  z = -----------------
          lambda_0

For large velocities special relativity should be taken into account
when computing the redshift from recessional velocity and vice-versa::
      ( c + V_r ) 1/2
  z = ( ------- )     - 1
      ( c - V_r )

             2
  V_r   (z+1)  - 1
  --- = ----------
   c         2
        (z+1)  + 1

The radial velocity which gives rise to redshift can be calculated
according to various conventions.  The most convenient is the radio
astronomy convention which computes redshift as::
      f_0 - f
  z = -------
        f_0
        
The velocity shift in the radio convention is computed as::
               df
  dV      = -c ---
    radio      f_0
  
The optical astronomy convention leads to::
                 df  (f_0)2 
  dV        = -c --- (---)
    optical      f_0 ( f )
  

Reference
=========

Redshift
--------
Lang, Astrophys. Formulae, eqns. 2-227 -- 2-229

Cosmic Time Scale
-----------------
1/Ho; Lang, Astrophys. Formulae, eqn. 5-82

Doppler Shift
-------------
astropy documentation for Spectral Doppler Equivalencies
"""
import astropy.units as u
import logging

from astropy.coordinates import EarthLocation, SkyCoord
from math import pi, sqrt
from novas import compat as novas
from numpy import array

from Astronomy import MJD, v_sun
from MonitorControl.Configurations.coordinates import DSS

c            = 3e5;  # km/s
Ho           = 73.8 # ± 2.40 (km/s)/Mpc
Mpc          = 3.0856e19;  # km
sec_per_year = 365.0*24*60*60
To           = Mpc/Ho/sec_per_year # age of universe,

logger = logging.getLogger(__name__)

# ---------------------------- redshift functions -----------------------------

def v_recess (distance, hubble_constant=Ho):
  """
  velocity of recession for a constant rate of expansion = Hubble constant

  @param distance : Mpc
  @type  distance : float

  @param hubble_constant : optional, (km/s)/Mpc
  @type  hubble_constant : float

  @return: recessional velocity in km/s
  """
  return hubble_constant * distance

def v_redshift (z):
  """
  velocity of recession for a given redshift

  @param z : redshift
  @type  z : float

  @return: recessional velocity in km/s (float)
  """
  return c * ((z+1)**2 - 1) / ((z+1)**2 + 1)

def redshift (v_recess=None, dop_fac=None):
  """
  Redshift corresponding to a recessional velocity

  @param v_recess : recessional velocity in km/s
  @type  v_recess : float

  @param dop_fac : f/f_0 or lambda_0/lambda
  @type  dop_fac : float

  @return: redshift (float)
  """
  if v_recess == None and dop_fac == None:
    raise RuntimeError("No Doppler velocity or Doppler factor given.")
  elif v_recess:
    return sqrt((c + v_recess)/(c - v_recess)) - 1
  else:
    return (1-dop_fac)/dop_fac

def dop_fac (z):
  """
  Redshifted frequency as a ratio: f/f_0

  @param z : redshift
  @type  z : float

  @return: fractional frequency reduction (float)
  """
  return 1./(z+1)

def wave_stretch (z):
  """
  Fractional lengthening of wavelength

  @param z : redshift
  @type  z : float

  @return: fractional wavelength increase (float)
  """
  return z+1
  

def distance (z, hubble_constant=Ho):
  """
  Distance of a cosmically redshifted source, assuming constant expansion
  
  @param z : redshift
  @type  z : float

  @param hubble_constant : optional
  @type  hubble_constant : float

  @return: distance in Mpc (float)
  """
  v_r = v_redshift(z)
  return v_r/hubble_constant

def age (z, hubble_constant=Ho):
  """
  Cosmic time scale = age if acceleration = 0

  @param z : redshift
  @type  z : float

  @param hubble_constant : optional
  @type  hubble_constant : float

  @return: time scale (age) in years (float
  """
  today = Mpc/hubble_constant/sec_per_year
  d = distance( z, hubble_constant)
  return (today - d*Mpc/c/sec_per_year)/1.0e6

# ---------------------------- doppler shift functions ------------------------
    
def doppler_radio(rel_freq, f_ref):
  """
  Returns the Doppler shift for a frequency offset by the radio convention
  
  @param rel_freq : frequency offset from the reference frequency
  @type  rel_freq : float or nparray of float
  
  @param f_ref : reference frequency
  @type  f_ref : float
  
  @return: float
  """
  return -c*rel_freq/f_ref
      
def doppler_optical(rel_freq, f_ref):
  """
  Returns the Doppler shift for a frequency offset by the optical convention
  
  @param rel_freq : frequency offset from the reference frequency
  @type  rel_freq : float or nparray of float
  
  @param f_ref : reference frequency
  @type  f_ref : float
  
  @return: float
  """
  return doppler_radio(rel_freq, f_ref)*(f_ref/(f_ref+rel_freq))**2
    
def doppler_relat(rel_freq, f_ref):
  """
  Returns the relativistic Doppler shift for a frequency offset
  
  @param rel_freq : frequency offset from the reference frequency
  @type  rel_freq : float or nparray of float
  
  @param f_ref : reference frequency
  @type  f_ref : float
  
  @return: float
  """
  r = (f_ref+rel_freq)/f_ref
  return doppler_radio(rel_freq, f_ref)*(4*r)/(1+r**2)**2

# ------------------------------- local standard of rest ----------------------


def V_LSR(RA, dec, dss, timedate):
  """
  Computes the velocity of the local standard of rest w.r.t. the observer
  
  @param ra : J2000 right ascension as a float or as "12h34m56.78s"
  @type  ra : float or str
  
  @param dec : J2000 declination as a float or as "-12d34m56.78s"
  @type  dec : float or str
  
  @param observer : DSN station
  @type  observer : int
  
  @param timedate : date/time of the observation
  @type  timedate : datetime object
  """
  if type(RA) == unicode and type(dec) == unicode:
    skypos = SkyCoord(RA, dec, unit=(u.hourangle, u.deg))
  elif type(RA) == float and type(dec) == float:
    skypos = SkyCoord(RA*u.hour,dec*u.degree)
  else:
    raise RuntimeError(RA, dec, "cannot be parsed")
  logger.debug("V_LSR: sky pos: %s", skypos)
  ra2000, dec2000 = skypos.ra.hour, skypos.dec.deg
  logger.debug("V_LSR: J2000 coordinates are %f, %f", ra2000, dec2000)
  sourcename = "%5.2f%+5.2f" % (ra2000,dec2000)
  cat_entry = novas.make_cat_entry(sourcename,
                                   "",0,
                                   ra2000, dec2000,
                                   0, 0, 0, 0)
  source = novas.make_object(2, 0, sourcename, cat_entry)
  station = DSS(dss)
  logger.debug("V_LSR: station lat=%f", station.lat*180/pi)
  logger.debug("V_LSR: station long=%f", station.lon*180/pi)
  observer = novas.make_observer_on_surface(station.lat*180/pi,
                                         station.lon*180/pi,
                                         station.elev, 0, 0)
  jd = novas.julian_date(timedate.year, timedate.month, timedate.day,
                         timedate.hour+timedate.minute/60.)
  mjd = MJD(timedate.year, timedate.month, timedate.day)
  earth = novas.make_object(0, 3, 'Earth', None)
  urthpos,urthvel = novas.ephemeris((jd,0), earth, origin=0)
  (obspos,obsvel) = novas.geo_posvel(jd,0,observer,0)
  totvel = tuple(array(urthvel)+array(obsvel))
  (srcpos,srcvel) = novas.starvectors(cat_entry)
  V = novas.rad_vel(source, srcpos, srcvel, totvel,0,0,0)
  logger.debug("V_LSR: velocity of observer w.r.t. Sun= %.2f", V)
  return V+v_sun(mjd,ra2000/15.,dec2000)


