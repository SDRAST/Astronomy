import logging

import ephem

from . import PD, Jnames

module_logger = logging.getLogger(__name__)

class Pulsar(ephem.FixedBody):
  """
  ephem.FixedBody() with pulsar properties.

  The pulsar astrophysical properties can be accessed with:

  .. code-block:: python

      >>> x.properties
      {'CLK': 'UTC(NIST)',
       'DIST_DM': '+17.82',
       'DIST_DM1': '6.33',
       'DM': '240.70',
       'EPHEM': 'DE200',
       'FINISH': '51889.103',
       'P0': '0.313885674748',
       'P1': '1.521E-16',
       'PEPOCH': '51650.000',
       'S1400': '0.29',
       'START': '51396.419',
       'SURVEY': 'swmb,pksmb',
       'W50': '12'}

  To find out what a property is:

  .. code-block:: python

    >>> print Ephem.PD.help['S1400']
    
  (PD stands for 'pulsar data'.)
  """
  def __init__(self, name):
    """
    Create a Pulsar() instance from its name.

    The preferred name form is like 'J2305+4707'.  However, 'B2303+46'
    is also exceptable, as is '2305+4707'.  In the latter case, it must
    be a Julian epoch name.

    Notes
    =====
    Instance attributes
    -------------------
    properties - a dictionary of pulsar physical properties

    Inherited methods and attributes
    --------------------------------
    All methods and attributes of ephem.FixedBody()
    """
    super(Pulsar,self).__init__()
    if name[0] == "J":
      self.name = name
    elif name[0] == "B":
      self.name = PD.names[name]
    else:
      self.name = "J"+name
    try:
      Jnames.index(self.name)
    except:
      raise self.name
    if diag:
      print(self.name)
    pulsar_data = PD.data[self.name]
    [ra,dec]    = PD.equatorial(pulsar_data)
    if diag:
      print("From catalogue:", ra, dec)
      print("In ephem format:",ephem.hours(str(ra)),ephem.degrees(str(dec)))
    self._ra    = ephem.hours(str(ra))     # pulsar_data['RAJ']  or derived
    self._dec   = ephem.degrees(str(dec)) # pulsar_data['DECJ'] or derived
    self._pmra  = float(pulsar_data['PMRA'])
    self._pmdec = float(pulsar_data['PMDEC'])
    self.period = PD.period(pulsar_data)
    self.dpdt   = PD.period_change_rate(pulsar_data)
    if diag:
      print("As attributes:",self._ra, self._dec)
    self._epoch = J2000
    self._class = "L"
    # to initialize ra, dec to something
    self.compute("2000/1/1 00:00:00")
    keys = list(pulsar_data.keys())
    keys = remove_item(keys,'RAJ')
    keys = remove_item(keys,'DECJ')
    keys = remove_item(keys,'ELONG')
    keys = remove_item(keys,'ELAT')
    keys = remove_item(keys,'PMRA')
    keys = remove_item(keys,'PMDEC')
    keys = remove_item(keys,'PSRB')
    keys = remove_item(keys,'P0')
    keys = remove_item(keys,'P1')
    if diag:
      print("Properties:",keys)
    self.properties = {}
    for prop in keys:
      self.properties[prop] = pulsar_data[prop]
    self.pressure = 0
