"""
This is deprecated.

Instead of
```python
from support.Ephem import DSS
```

you should do the following:

```python
from MonitorControl.Configurations.coordinates import DSS
```

"""
import logging

import ephem

import Astronomy
# from MonitorControl.Configurations import coordinates as coord

dsn = Astronomy.get_geodetic_coords()
xyz = Astronomy.get_cartesian_coordinates()

module_logger = logging.getLogger(__name__)

class DSS(ephem.Observer):
  """
  Observer class for DSN stations
  """
  def __init__(self, number):
    """
    Creates an instance of a DSN station Observer()

    Attributes
    ==========
    lon       - east longitude (degrees)
    lat       - north latitude (degrees)
    elevation - altitude (meters)
    timezone  - difference (hours) from UTC
    name      -
    diam      - diameter (m)

    @param number : station number
    """
    super(DSS,self).__init__()
    module_logger.debug("Initiating DSS(%d)", number)
    self.lon = -dsn[number][0]*pi/180.
    self.lat = dsn[number][1]*pi/180.
    self.elevation = dsn[number][2]
    self.timezone = dsn[number][3]
    self.name = dsn[number][4]
    self.diam = dsn[number][5]
    self.xyz = xyz["DSS %2d" % number]
