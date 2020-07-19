# -*- coding: utf-8 -*-
"""
Module Ephem - extends module pyephem for radio astronomy

Module Ephem extends module ephem (package pyephem).  Module 'ephem' is
imported under its own name, so these two are equivalent::

    import Ephem.ephem
    import ephem

Function ``calibrator()`` returns a planet or radio source suitable for 
calibrating fluxes::

  In [2]: venus = calibrator("Venus")                                             
  In [3]: type(venus)                                                             
  Out[3]: ephem.Venus

  In [8]: c3_273 = calibrator("3C273")                                            
  In [9]: type(c3_273)                                                            
  Out[9]: Astronomy.Ephem.quasar.Quasar
  
  In [10]: from Radio_Astronomy.radio_flux import get_calibrator_flux             
  In [12]: import datetime                                                        
  In [13]: now = datetime.datetime.now()                                          
  In [14]: get_calibrator_flux('Venus', 8.4, now)                                 
  Out[14]: (28.9862540298923, 'Planet')

Class ``DSS``, a subclass of ``Observer`` is imported from ``.DSS_coordinates``
and is useful in this context.
 
  In [15]: from Astronomy.Ephem import DSS                                        
  In [16]: dss43 = DSS(43)                                                        
  In [17]: type(dss43)                                                            
  Out[17]: Astronomy.DSN_coordinates.DSS
  In [18]: print(dss43.long, dss43.lat, dss43.elev)
  -211:01:11.8 -35:24:14.3 688.867

Class ``Pulsar`` provides pulsar physical data::

  In [1]: from Astronomy.Ephem import Pulsar                                      
  In [2]: psr = Pulsar('B0531+21')                                                
  In [3]: type(psr)                                                               
  Out[3]: Astronomy.Ephem.pulsar.Pulsar
  In [4]: psr.properties                                                          
  Out[4]: {'POSEPOCH': '40675',      'DM': '56.791',       'DIST_A': '2.0',
           'S600': '211',            'S925': '45',         'S1400': '14',
           'NGLT': '25',             'RM': '-42.3',        'SPINDX': '-3.1',
           'DIST_DM': '2.49',        'PEPOCH': '40000.00', 'DIST_AMN': '1.5',
           'SURVEY': 'misc,ar4,gb4', 'W50': '3.0',         'DIST_AMX': '2.5',
           'F0': '30.2254370',       'F1': '-3.86228E-10', 'F2': '1.2426E-20',
           'F3': '-0.64E-30',        'EPHEM': 'DE405',     'DIST_DM1': '1.74',
           'TAU_SC': '1.51e-06',     'W10': '4.7',         'S400': '646',
           'ASSOC': 'SNR:Crab[ccl+69],GRS:1FGL_J0534.5+2200[aaa+10g],\
                     GRS:HESS_J0534+220[aab+06b]',
           'TYPE': 'HE[cdt69,fhm+69,hjm+70]'}

Notes
=====

module ephem
------------

http://rhodesmill.org/brandon/projects/pyephem-manual.html

Example of use::
  
  In [1]: from ephem import *
  In [2]: ba = city("Buenos Aires")                                               
  In [3]: venus = Venus(ba)                                                       
  In [4]: venus.compute("2020/6/11 16:00:00") 
                 
  In [5]: print(venus.a_ra, venus.a_dec, venus.a_epoch)                           
  4:28:37.11 20:25:09.7 2000/1/1 12:00:00
  In [6]: print(venus.g_ra, venus.g_dec)                                          
  4:29:46.74 20:27:41.9
  In [7]: print(venus.ra, venus.dec)                                              
  4:29:46.74 20:27:41.9


About Celestial Coordinates
---------------------------

* a_ra, a_dec: Astrometric geocentric position for epoch (e.g. J2000)
* g_ra, g_dec: Apparent geocentric position for date specified in
  the compute(), after correcting for precession, relativistic
  deflection, nutation and aberration
* ra, dec: Apparent topocentric position, after correction for parallax
  and refraction. Set the ``Observer`` attribute pressure to zero if
  you want PyEphem to ignore the effects of atmospheric refraction

"""
import logging
import datetime
from math import pi

import ephem

module_logger = logging.getLogger(__name__)

from Astrophysics.Pulsars import pulsar_data as PD
import Radio_Astronomy as RA
from Radio_Astronomy import michigan, vla_cal

#try:
#    from Astronomy.DSN_coordinates import DSS
#except ImportError as err:
#    module_logger.error(("Couldn't import support version of DSS Observer."
#                         "Falling back to unsupported version in this package."))
#    from .dss import DSS

J2000 = ephem.Date("2000/1/1 00:00:00")

# planets recognized by module ephem
Planets = ['Jupiter', 'Mars', 'Mercury', 'Moon', 'Neptune', 'Pluto',
           'Saturn', 'Sun', 'Uranus', 'Venus']

try:
    Jnames = list(PD.data.keys()) # pulsar Julian epoch names
    Jnames.sort()

    cal_dict = vla_cal.get_cal_dict() # VLA calibrators
    Bname_dict, cat_3C_dict = vla_cal.VLA_name_xref(cal_dict) # name dictionaries
    Bnames = vla_cal.Jnames_to_B(Bname_dict) # B names keyed on J names
except Exception as err:
    module_logger.error("Couldn't process pulsar or vla calibration data: {}".format(err))

class EphemException(Exception):
  """
  Exception class for Ephem module
  """
  def __init__(self, value=None, details=None):
    """
    Creates an instance of EphemException()
    """
    self.value = value
    self.details = details

  def __str__(self):
    """
    EphemException() instance message
    """
    msg = "error"
    if self.value:
      msg += ": " + repr(self.value)
    if self.details:
      msg += ", " + self.details
    return repr(msg)

# these are here so they can see Jnames, Bname_dict, etc.

from .pulsar import Pulsar
from .quasar import Quasar
from .serializable_body import SerializableBody

__all__ = [
        "Quasar",
        "Pulsar",
        "SerializableBody",
        "EphemException",
        "calibrator",
        "Planets",
        "PD"
]

def calibrator(name):
  """
  Creates an Ephem (pyephem) Body() instance for a flux calibrator
  """
  global calsource # needed for exec() below
  try:
    # planet?
    idx = Planets.index(name)
    module_logger.debug("calibrator: %s is planet %d", name, idx)
  except:
    module_logger.info("calibrator: %s is not planet", name)
    try:
      calsource = Quasar(name)
      return calsource
    except Exception as details:
      # not a Quasar either
      module_logger.error("calibrator: not a quasar: {}".format(details))
      return None  
  # It's a planet
  try:
    code_line = 'global calsource; calsource = ephem.'+name+'()'
    module_logger.debug("calibrator: executing '%s'", code_line)
    exec(code_line)
  except Exception as errmsg:
    module_logger.info("calibrator: exec error: {}".format(errmsg))
  if calsource:
    module_logger.debug("calibrator: is %s", calsource)
    return calsource
  else:
    module_logger.debug("calibrator: exec() left 'calsource' as None")
    return None    

