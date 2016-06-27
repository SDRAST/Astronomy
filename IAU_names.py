# -*- coding: utf-8 -*-
"""
module IAU_names - IAU source name conversion and matching

IAU source nomenclature rules are here:
http://cdsweb.u-strasbg.fr/Dic/iau-spec.html

Coordinate based IAU source names can have the following forms.  Lowercase
letters indicate a decimal part, except s which means sign.
* Explicit use of a decimal point is recommended.
* A source type designation is also recommended, such as H2O, PSR, QSO, etc.
* No spaces may appear in the designator.
* Leading zeros are required.

For sources in our galaxy::
  G LLL[.ll]sBB[.bb]
  G LLLllsBBbb
    The G is required. The decimal part of the coordinates, to any precision,
    is optional.
    
For sidereal sources::
  J HHMM[SS[.s]]sDDMM[SS[.s]]
  J HHMMsDDMM[.m]
  
"""
diag = False

import re
from math import pow

def split_on_sign(designator):
  """
  Separate the ra/l part from the dec/b part of an IAU designation,

  It's not the responsibility of this function to filter out bad
  designations.  These are examples of valid ones::
    G001.2-00.3
    G001-00
    J1426.8+6950
    J1302-6350
    B004848-4242.8
    004848-4242.8
    0048-427
    0048-42

  @param designator : an IAU source designation
  @type  designator : str

  @return: (x-coordinate, y-coordinate, sign for y)
  """
  if diag:
    print "split_on_sign: received",designator
  if re.search('\+',designator):
    if diag:
      print "Splitting on +"
    parts = designator.split('+')
    parts.append('+')
  elif re.search('-',designator):
    if diag:
      print "Splitting on -"
    parts = designator.split('-')
    parts.append('-')
  else:
    if diag:
      print "No split"
    parts = ["",designator,""]
  # pulsar names can have a "-1, -2, ..." suffix.
  if diag:
    print "Parts found:",parts
  if len(parts) > 3:
    parts = parts[:3]
  elif re.search('-',parts[1]):
    subparts = parts[1].split("-")
    parts = parts[:1] + [subparts[0]]
  # Pulsars may have alphabetic suffixes of one more more characters
  else:
    test = parts[1]
    for i in range(1,len(test)):
      if diag:
        print i,"Testing",test
      if test[-1].isalpha():
        test = test[:-1]
      else:
        parts[1] = test
        break
    if diag:
      print "stripped is",parts[1]
  return parts

def parse_decimal_angle(angle):
  """
  Parse an angle given in degrees, with or without a decimal point.

  A decimal point may be implied by the length of the string.

  @param angle : angle in degrees to be converted
  @type  angle : str

  @return: decimal angle (float)
  """
  if re.search('.',angle) or len(angle) == 2 or len(angle) == 3:
    da = float(angle)
  elif len(angle) > 3:
    exponent = len(angle)-3
    divisor = pow(10,exponent)
    da = float(angle)/divisor
  else:
    # Must be a missing leading zero.
    da = float(angle)
  return da

def parse_sexagesimal_angle(angle):
  """
  Parse a sexagesimal angle, with or without an explicit decimal point.

  The position of a decimal point is inferred from the length of the string

  @param angle : angle as DD, DDMM, DDMM.m, DDMMm, DDMMSS, DDMMSS.s, DDMMSSs
  @type  angle : str

  @return: decimal angle
  """
  if diag:
    print "parse_sexagesimal_angle: >"+angle+"<"
  if re.search('\.',angle):
    DDMMSS,second_fraction = angle.split('.')
    second_fraction = float("0."+second_fraction)
  else:
    DDMMSS = angle
    second_fraction = 0.
  # At the very least, the first four digits should be unambiguous
  if diag:
    print "parse_sexagesimal_angle: >"+DDMMSS+"<"
  DD = int(DDMMSS[:2])
  if len(DDMMSS) == 2:
    # This should never happen for right ascension
    aa = DD
  elif len(DDMMSS) == 3:
    # This should never happen for right ascension
    aa = DD + float(DDMMSS[2])/10
  else:
    MM = int(DDMMSS[2:4])
    aa = DD + MM/60.
  if len(DDMMSS) == 6:
    # It should be this if there was a decimal
    SS = int(DDMMSS[4:6])
    aa += (SS+second_fraction)/3600.
  elif len(DDMMSS) == 5:
    MM += float(DDMMSS[4])/10
    aa += MM/60.
  elif len(DDMMSS) > 6:
    SS = int(DDMMSS[4:6])
    extra_digits = len(DDMMSS) - 6
    divisor = pow(10,extra_digits)
    second_fraction = float(DDMMSS[6:])/divisor
    aa += (SS+second_fraction)/3600.
  return aa
  
def parse_IAU_name(name):
  """
  Parse an IAU source designation.

  @param name : IAU source name
  @type  name : str

  @return: (flag letter, x-coordinate, y-coordinate)
  """
  # First see if there is a source type acronym
  if diag:
    print "parse_IAU_name: received",name
  parts = name.split()
  if len(parts) == 1:
    designation = parts[0]
  elif len(parts) == 2:
    acronym, designation = parts
  else:
    raise("Invalid format: "+name)
  # Now process the designation
  flag = designation[0].upper()
  if flag == "G":
    # Galactic coordinates
    longitude,latitude,sign = split_on_sign(name[1:])
    X = parse_decimal_angle(longitude)
    Y = parse_decimal_angle(latitude)
  elif flag == "J":
    # Julian epoch celestial coordinates
    ra,dec,sign = split_on_sign(name[1:])
    X = parse_sexagesimal_angle(ra)
    Y = parse_sexagesimal_angle(dec)
  elif flag == "B":
    # Besselian epoch celestial coordinates
    ra,dec,sign = split_on_sign(name[1:])
    X = parse_sexagesimal_angle(ra)
    Y = parse_sexagesimal_angle(dec)
  elif designation[0].isdigit():
    # This should be Besselian but who knows?
    # If it is Besselian there should be at least four digits in RA
    # otherwise it could be galactic
    x,y,sign = split_on_sign(name)
    if len(x) > 3:
      X = parse_sexagesimal_angle(x)
      Y = parse_sexagesimal_angle(y)
      flag = "B"
    else:
      X = parse_decimal_angle(x)
      Y = parse_decimal_angle(y)
      flag = "G"
  else:
    return "?",None,None
  if sign == "-":
    Y = -Y
  return flag,X,Y

def match_IAU_names(name1, name2, template):
  """
  Compare two IAU designations.

  Examples of a valid template for one digit of precision is "%4.1f %4.1f"
  """
  if diag:
    print "match_IAU_names:",name1,name2
  flag1,x1,y1 = parse_IAU_name(name1)
  flag2,x2,y2 = parse_IAU_name(name2)
  if flag1 != flag2:
    return False
  string1 = template % (x1,y1)
  string2 = template % (x2,y2)
  if string1 == string2:
    return True
  else:
    return False

def find_IAU_name_match(name,namelist,template):
  """
  Check a name against all the names in a list and return the match
  """
  for nm in namelist:
    if match_IAU_names(name, nm, template):
      return nm
  return None
