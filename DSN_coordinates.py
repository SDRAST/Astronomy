"""
Coordinates for DSN and affiliated antennas

The DSN uses west longitudes, whereas most astronomical packages like 
``astropy`` and ``ephem`` use east longitudes.


Contents
========

In this module::

  class DSS
  function DSSLocation
  
"""
import ephem
import math

import astropy.units as u
import astropy.coordinates

Complex = {"Canberra":  [34, 35, 36, 43, 45],
           "Goldstone": [13, 14, 15, 24, 25, 26],
           "Madrid":    [53, 54, 55, 63, 65]}
complexID = {"Canberra":  "CDSCC",
             "Goldstone": "GDSCC",
             "Madrid":    "MDSCC"}
complexCode = {"Canberra": 40,
               "Goldstone": 10,
               "Madrid": 60}

class DSS(ephem.Observer):
  """
  Observer suclass class for DSN stations based on xephem Observer

  For use with xephem.  

  Attributes
  ==========
  
  of this class::
    
    lon       - east longitude (radians)
    lat       - north latitude (radians)
    elevation - altitude (meters)
    timezone  - difference (hours) from UTC
    name      -
    diam      - diameter (m)
    xyz       - geocentric Cartesian coordinates
    
  Notes
  =====
  
  * ``long`` is an alias for ``lon``.
  * Longitudes are measured eastwards.
      
  Examples
  ========

    In [8]: import ephem as E 
   
    In [18]: LA = E.city("Los Angeles")
    In [19]: LA.long, LA.lat
    Out[19]: (-2.0637416211957023, 0.5943236044502173)
    In [20]: print(LA.long, LA.lat)
    -118:14:37.3 34:03:08.0
    
    In [25]: from Astronomy.DSN_coordinates import DSS
    In [26]: dss14 = DSS(14)
    In [27]: dss14.long, dss14.lat
    Out[27]: (-2.0400918530711474, 0.6182990806837911)
    In [28]: print(dss14.long, dss14.lat)
    -116:53:19.2 35:25:33.3

  """
  def __init__(self, number):
    """
    Creates an instance of a DSN station Observer()
    
    Note the longitude sign change from DSN usage to astronomical convention.

    @param number : station number
    """
    super(DSS,self).__init__()
    self.number = number
    dsn = get_geodetic_coords()
    self.lon       = -dsn[number][0]*math.pi/180.
    self.lat       =  dsn[number][1]*math.pi/180.
    self.elevation =  dsn[number][2]
    self.timezone  =  dsn[number][3]
    self.name      =  dsn[number][4]
    self.diam      =  dsn[number][5]
    xyz = get_cartesian_coordinates()
    self.xyz = xyz["DSS %2d" % number]

  def copy(self):
    return self.__class__(self.number)

def DSSLocation(dss):
  """
  Gives DSN station location as astropy EarthLocation
  
  Notes that ``DSS`` longitudes are positive to the east and negative to the 
  west.

  @param dss : DSN station
  @type  dss : DSS class instance or station number

  @return: EarthLocation instance
  """
  if isinstance(dss, DSS):
    station = dss
    number = dss.number
  else:
    station = DSS(dss)
    number = dss
  loc = astropy.coordinates.EarthLocation(lat= (station.lat*180/math.pi)*u.deg,
                                          lon= (station.long*180/math.pi)*u.deg,
                                          height=station.elevation*u.m)
  loc.number = number
  return loc

geodetic_coords = {
12:  (116.804607   ,  35.2999592   , 988.93 ,-8, 34,"Goldstone Echo"),
13:  (116.794009   ,  35.2477189   ,1093.54 ,-8, 26,"Goldstone Venus"),
14:  (116.888653   ,  35.4259278   ,1031.81 ,-8, 70,"Goldstone Mars"),
15:  (116.887193   ,  35.42185386  , 973.945,-8, 34,"Goldstone DSS-15"),
16:  (116.8736477  ,  35.34153998  , 944.711,-8, 26,"Goldstone DSS-16"),
19:  (107.6        ,  34.1         ,2124.   ,-7,100,"VLA Soccoro NM"),
21:  (118.172743   ,  34.199966    , 100.   ,-8,  1,"JPL"),
22:  (118.172743   ,  34.199966    , 100.   ,-8,  1,"JPL"),
23:  (116.87286099 ,  35.33955093  , 946.086,-8, 34,"Goldstone DSS-23"),
24:  (116.8747921  ,  35.3398932   , 952.156,-8, 34,"Goldstone DSS-24"),
25:  (116.8753616  ,  35.33761236  , 960.862,-8, 34,"Goldstone DSS-25"),
26:  (116.873015   ,  35.33568948  , 970.159,-8, 34,"Goldstone DSS-26"),
27:  (116.7766484  ,  35.23827236  ,1053.203,-8, 34,"Goldstone DSS-27"),
28:  (116.778889   ,  35.2382726   ,1065.382,-8, 34,"GAVRT DSS-28"),
32:  (243.80849    , -31.048225    , 252.260, 7, 35,"New Norcia DSA-1"),
33:  (211.0169082  , -35.4004847361, 684.099,10, 11,"Canberra DSS-33"),
34:  (211.018035581, -35.3984788417, 692.020,10, 34,"Canberra DSS-34"),
35:  (211.0185442,   -35.2143052,    694.889,10, 34,"Canberra DSS-35"),
36:  (211.0214558,   -35.2136127,    685.503,10, 34,"Canberra DSS-36" ),
42:  (211.019943   , -35.402232    , 656.08 ,10, 34,"Canberra retired"),
43:  (211.019942862, -35.403983527 , 688.867,10, 70,"Canberra DSS-43"),
45:  (211.022333334, -35.4         , 674.347,10, 34,"Canberra DSS-45"),
46:  (211.016918317, -35.4050106361, 676.812,10, 26,"Canberra DSS-46"),
48:  (211.7376993  , -32.9999629   , 392    ,10, 64,"Parkes"),
53:  (  4.249652169,  40.4273572916, 826.791, 1, 11,"Madrid DSS-53"),
54:  (  4.249654703,  40.427355656 , 837.051, 1, 34,"Madrid DSS-54"),
55:  (  4.2526333  ,  40.4242959028, 819.061, 1, 34,"Madrid DSS-55"),
61:  (  4.247700556,  40.429921388 , 740    , 1, 26,"Madrid DSS-61"),
62:  (  4.36755    ,  40.452688    , 738.87 , 1, 35,"Cebreros DSA-2"),
63:  (  4.247700556,  40.429921388 , 864.816, 1, 70,"Madrid DSS-63"),
65:  (  4.250698897,  40.4272063583, 833.854, 1, 34,"Madrid DSS-65"),
66:  (  4.251417631,  40.4299748722, 849.874, 1, 26,"Madrid DSS-66"),
71:  (  0.,            0.,           0.,      0, 34,"MIL 71"),
74:  (243.80849    , -31.048225    , 252.260, 7, 35,"New Norcia DSA-1"),
83:  (  4.36755    ,  40.452688    , 738.87 , 1, 35,"Cebreros DSA-2"),
84:  ( 69.398      , -35.776       ,1550.   ,-3, 35,"Malargue DSA-3"),
95:  ( 77.368619   ,  12.901631    ,3000.   ,5.5,32,"Indian DSN")}

def get_geodetic_coords(dss=0,observatory=None):
  """
  Observatory's geodetic coordinates

  An observatory may be specified with a DSS number or a name.

  If no observatory is specified, this returns a dictionary in which
  the station IDs are the keys and the values are tuples.

  If the observatory is specified it is used as a key and only that
  observatory's coordinates are returned as a tuple.

  If the observatory is not found, an empty dictionary is returned.

  Notes
  -----
  For an explanation of Earth coordinate systems so
  http://dsnra.jpl.nasa.gov/Antennas/Antennas.html#anchor950381

  Example::

    >>> station_data = get_geodedic_coords()
    >>> long,lat,alt,tz,name, diam = station_data[13]
    >>> print long,lat
    116.794009 35.2477189

  Todo
  ----
  Make both functions more forgiving in the matter of observatory names

  References
  ----------
  https://deepspace.jpl.nasa.gov/dsndocs/810-005/301/301K.pdf Tables 5 and 6

  @param dss : int
    optional DSS numbercomplexCode

  @param observatory : string
    optional observatory name

  @return: dictionary of tuples
    (longitude, latitude, altitude(m), time zone, name, diameter)
    keyed to the station name(s), if any
  """
  if dss == 0 and observatory == None:
    result = {}
    for key in geodetic_coords.keys():
      long,lat,elev,tz,diam,name = geodetic_coords[key]
      result[key] = long,lat,elev,tz,name,diam
    return result
  elif dss == 0 and observatory != None:
    for key in geodetic_coords.keys():
      long,lat,elev,tz,diam,name = geodetic_coords[key]
      if observatory.strip().lower() == name.strip().lower():
        return long,lat,elev,tz,name,diam
  elif type(dss) == int:
    long,lat,elev,tz,diam,name = geodetic_coords[dss]
    return long,lat,elev,tz,name,diam

def get_cartesian_coordinates(station=None):
  """
  Get the Cartesian coordinates of a DSN station, or a dictionary of all

  @param station : 13 or 'Venus' or 'DSS 13'
  @type  station : int or string

  @return: tuple or dict of tuples or None, Cartesian coordinates in meters

  This creates a dictionary with the Cartesian coordinates in meters
  of the DSN stations using the ITRF1993 (assuming subreflector-fixed
  configuration). If a valid station ID is given, it returns the coordinates
  of that station.  None is returned if the ID is invalid.  If no ID is
  given, the entire dictionary is returned.

  The entry for DSS 21 is a rough one for JPL.

  Notes
  -----  
  For an explanation of the coordinate systems see
  http://dsnra.jpl.nasa.gov/Antennas/Antennas.html#anchor950381
  
  Some other stations have been added.
  
  An example to get a baseline length::
  
    >>> x1,y1,z1 = get_cartesian_coordinates('DSS 24')
    >>> x2,y2,z2 = get_cartesian_coordinates('DSS 13')
    >>> print math.sqrt(math.pow(x2-x1,2) + math.pow(y2-y1,2)+math.pow(z2-z1,2))
    12621.4825356

  References
  ----------
  https://deepspace.jpl.nasa.gov/dsndocs/810-005/301/301K.pdf Table 2
  
  """

  coordinates = \
    {"GB OVLBI": ( 884084.2636, -4924578.7481, 3943734.3354),
     "DSS 12":   (-2350443.812, -4651980.837, +3665630.988),
     "Echo":     (-2350443.812, -4651980.837, +3665630.988),
     "DSS 13":   (-2351112.491, -4655530.714, +3660912.787),
     "Venus":    (-2351112.491, -4655530.714, +3660912.787),
     "DSS 14":   (-2353621.251, -4641341.542, +3677052.370),
     "Mars":     (-2353621.251, -4641341.542, +3677052.370),
     "DSS 15": (-2353538.790, -4641649.507, +3676670.043), \
     "DSS 16": (-2354763.158, -4646787.462, +3669387.069), \
     "DSS 17": (-2354730.357, -4646751.776, +3669440.659), \
     "DSS 21": (-2350000.,    -4700000.   , +3700000.   ), \
     "DSS 22": (-2350000.,    -4700000.   , +3700000.   ), \
     "DSS 23": (-2354757.567, -4646934.675, +3669207.824), \
     "DSS 24": (-2354906.495, -4646840.128, +3669242.317), \
     "DSS 25": (-2355022.066, -4646953.636, +3669040.895), \
     "DSS 26": (-2354890.967, -4647166.925, +3668872.212), \
     "DSS 27": (-2349915.260, -4656756.484, +3660096.529), \
     "DSS 28": (-2350101.849, -4656673.447, +3660103.577), \
     "DSS 32": (-1666531.345, +5209373.6709,-3270605.479), \
     "DSS 33": (-4461083.514, +2682281.745, -3674570.392), \
     "DSS 34": (-4461146.756, +2682439.293, -3674393.542),
     "DSS 35": (-4461273.084, +2682568.922, -3674152.089),
     "DSS 36": (-4461168.415, +2682814.657, -3674083.901),
     "DSS 42": (-4460981.016, +2682413.525, -3674582.072), \
     "DSS 43": (-4460894.585, +2682361.554, -3674748.580), \
     "DSS 45": (-4460935.250, +2682765.710, -3674381.402), \
     "DSS 46": (-4460828.619, +2682129.556, -3674975.508), \
     "Parkes": (-4554231.843, +2816758.983, -3454036.065), \
     "DSS 48": (-4554231.843, +2816758.983, -3454036.065), \
     "DSS 53": (+4849330.129,  -360338.092, +4114758.766), \
     "DSS 54": (+4849434.555,  -360724.108, +4114618.643), \
     "DSS 55": (+4849525.318,  -360606.299, +4114494.905), \
     "DSS 61": (+4849245.211,  -360278.166, +4114884.445), \
     "DSS 62": (+4846692.106,  -370171.532, +4116842.926), \
     "DSS 63": (+4849092.647,  -360180.569, +4115109.113), \
     "DSS 65": (+4849336.730,  -360488.859, +4114748.775), \
     "DSS 66": (+4849148.543,  -360474.842, +4114995.021), \
     "MIL 71": (0,0,0),
     "DSS 74": (0,0,0),
     "DSS 83": (0,0,0),
     "DSS 84": (0,0,0),
     "DSS 95": (0,0,0)}
  if type(station) == int:
    try:
      return coordinates["DSS %2d" % station]
    except:
      module_logger.error("DSS %2d is not known", station, exc_info=True)
      return None
  elif type(station) == str:
    try:
      return coordinates[station]
    except:
      module_logger.error("Invalid DSS ID: %s", station, exc_info=True)
      return None
  else:
    return coordinates

def DSN_complex_of(dss):
  """
  Returns the Complex to which a station belongs

  @param dss : DSN station number
  @type  dss : int or str

  @return: str
  """
  if type(dss) == str:
    dss = int(dss[-2:])
#     dss = int(float(dss[-2:]))
  for cmplx in Complex.keys():
    try:
      Complex[cmplx].index(dss)
      return cmplx
    except:
      pass
  return None
