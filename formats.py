"""Converts between conventional angular notation and decimal notation,
either decimal angles or radians as designated."""

import string
import math

def parse_hms_delimited_angle(angle):
    """
    Parses RA strings delimited with h, m and s.
    
    As in VLA Calibrator List: 00h01m08.621563s
    Example::
      >>> parse_hms_delimited_angle('00h01m08.621563s')
      ['00', '01', '08.621563']
    
    @param angle : R.A. formatted like 00h01m08.621563s
    @type  angle : str
    
    @return: list of str
    """
    temp = angle.split('h')
    h = temp[0]
    the_rest = temp[1]
    temp = the_rest.split('m')
    m = temp[0]
    s = temp[1].rstrip('s')
    return [h,m,s]

def parse_dms_delimited_angle(angle):
    """
    Parses decl. strings delimited with d, ' and ".
    
    As in VLA Calibrator List: 19d14'33.801860"
    Example::
      >>> parse_dms_delimited_angle('''19d14'33.801860"''')
      ['19', '14', '33.801860']
    
    @param angle : decl. formatted like 00h01m08.621563s
    @type  angle : str
    
    @return: list of str
    """
    temp = angle.split('d')
    d = temp[0]
    the_rest = temp[1]
    temp = the_rest.split("""'""")
    m = temp[0]
    s = temp[1].rstrip('''"''')
    return [d,m,s]
    
def parse_colon_delimited_angles(rastring, decstring):
    """
    Parses RA and decl. strings delimited with colons
    
    Input data like ('1:38:48.0','41:24:23')
    and it will return [1.64666666667 41.4063888889]
       
    @param rastring : R.A. formatted like '1:38:48.0'
    @type  rastring : str
    
    @param decstring : decl. formatted like '41:24:23'
    @type  decstring : str
    
    @return: list of float
    """
    ralist = rastring.split(":")
    if len(ralist) == 3:
        rah = int(ralist[0])
        ram = int(ralist[1])
        ras = float(ralist[2])
    elif len(ralist) == 2:
        rah = int(ralist[0])
        ram = int(ralist[1])
        ras = 0
    else:
        rad = float(ralist[0])
        ram = 0
        ras = 0
    if ( rastring[0] == "-" ):
        ram = -ram
        ras = -ras
    ra = (rah+(ram+(ras/60.))/60.)
    declist = decstring.split(":")
    if len(declist) == 3:
        decd = int(declist[0])
        decm = int(declist[1])
        decs = float(declist[2])
    elif len(declist) == 2:
        decd = int(declist[0])
        decm = int(declist[1])
        decs = 0
    else:
        decd = float(declist[0])
        decm = 0
        decs = 0
    if ( decstring[0] == "-" ):
        decm = -decm
        decs = -decs
    dec = (decd+(decm+(decs/60.))/60.)
    return [ra, dec]
    
def hms_delimited_angle_to_rads(angle):
    """
    Converts R.A. formatted as 00h01m08.621563s to radians
    
    As in VLA Calibrator List::
      >>> hms_delimited_angle_to_rads('00h01m08.621563s')
      0.0049903008842279899
    
    @param angle : R.A. formatted like 00h01m08.621563s
    @type  angle : str
    
    @return: float
    """
    parsed = parse_hms_delimited_angle(angle)
    sign = 1
    if ( parsed[0].find('-') > -1 ):
        sign = -1
    h = abs(float(parsed[0]))
    m = float(parsed[1])
    s = float(parsed[2])
    return sign*(h + (m + s/60.)/60.)*math.pi/12.
    
def dms_delimited_angle_to_rads(angle):
    """
    converts decl. formatted as 19d14'33.801860 to radians
    As in VLA Calibrator List: 19d14'33.801860"
    Example::
      >>> dms_delimited_angle_to_rads('''19d14'33.801860"''')
      0.33584886884199222
    
    @param angle : decl. formatted like 00h01m08.621563s
    @type  angle : str
    
    @return: float
    """
    parsed = parse_dms_delimited_angle(angle)
    sign = 1
    if ( parsed[0].find('-') > -1 ):
        sign = -1
    d = abs(float(parsed[0]))
    m = float(parsed[1])
    s = float(parsed[2])
    return sign*(d + (m + s/60.)/60.)*math.pi/180.

