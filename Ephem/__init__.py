# -*- coding: utf-8 -*-
"""
Module Ephem - extends module pyephem for radio astronomy

Module Ephem extends module ephem (package pyephem).  Everything from
module 'ephem' is inherited.  In addition, it provides the classes
``Quasar``, ``Pulsar`` and ``DSS``.  It also provides pulsar physical
data and quasar fluxes.  Use like this:

.. code-block:: python

    >>> from Astronomy import Ephem
    >>> x = Ephem.Pulsar('J1809-0743')

Notes
=====

module ephem
------------

http://rhodesmill.org/brandon/projects/pyephem-manual.html

About Celestial Coordinates
---------------------------

* a_ra, a_dec: Astrometric geocentric position for epoch (e.g. J2000)
* g_ra, g_dec: Apparent geocentric position for date specified in
  the compute(), after correcting for precession, relativistic
  deflection, nutation and aberration
* ra, dec: Apparent topocentric position, after correction for parallax
  and refraction. Set the ``Observer`` attribute pressure to zero if
  you want PyEphem to ignore the effects of atmospheric refraction

Diagnostics
-----------

To see what this module does internally:

.. code-block:: python

    >>> from Astronomy import Ephem
    >>> Ephem.diag = True
    >>> x = Ephem.Pulsar('J1809-0743')

"""
import logging
import datetime
from math import pi

import ephem

module_logger = logging.getLogger(__name__)

from Astrophysics.Pulsars import pulsar_data as PD
from Radio_Astronomy import michigan, vla_cal

try:
    from Astronomy.DSN_coordinates import DSS
except ImportError as err:
    module_logger.error(("Couldn't import support version of DSS Observer."
                         "Falling back to unsupported version in this package."))
    from .dss import DSS

J2000 = ephem.Date("2000/1/1 00:00:00")

diag = False

# planets recognized by module ephem
Planets = ['Jupiter', 'Mars', 'Mercury', 'Moon', 'Neptune', 'Pluto',
           'Saturn', 'Sun', 'Uranus', 'Venus']

try:
    Jnames = PD.data.keys() # pulsar Julian epoch names
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
    if value:
      self.value = value
    if details:
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
        "PD",
        "DSS"
]

def calibrator(name):
  """
  Creates an Ephem (pyephem) Body() instance for a flux calibrator
  """
  try:
    # planet?
    idx = Planets.index(name)
    module_logger.debug("calibrator: %s is planet %d", name, idx)
    code_line = 'calibrator = ephem.'+name+'()'
    module_logger.debug("calibrator: executing '%s'", code_line)
    exec(code_line)
    return calibrator
  except Exception as details:
    module_logger.debug("calibrator: not a planet: {}".format(details))
    # quasar?
    try:
      calibrator = Quasar(name)
      return calibrator
    except Exception as details:
      module_logger.debug("calibrator: not a quasar: {}".format(details))
      return None
