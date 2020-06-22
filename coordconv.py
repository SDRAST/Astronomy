"""
Conversion of spherical coordinates

This subroutine converts the longitude-like (A1) and latitude-like (A2) 
coordinates of a point on a sphere into the corresponding coordinates (A2,B2)
in a different coordinate system that is specified by the coordinates of its
origin (long_orig,lat_orig) and its north pole (AP,BP) in the original
coordinate system.

Examples of use
=============== 
1. HOUR ANGLE, DECLINATION --> AZIMUTH, ELEVATION::

       Az, El = coordconv(pi, pi/2-Lat, 0., Lat, HA, Dec)
       Azimuth will be in the range -Pi/2 to Pi/2

2. AZIMUTH, ELEVATION --> HOUR ANGLE, DECLINATION::

       HA, Dec = coordconv(pi, pi-2-Lat, 0., Lat, Az, El)

3. RIGHT ASCENSION, DECLINATION --> GALACTIC LONG & LAT (System I) 

  Refer to Kraus, P., RADIO ASTRONOMY, McGraw Hill, New York, 1966, or 
  Allen, C.W., ASTROPHYSICAL QUANTITIES, Athlone Press, London, (1963) for the
  definition of System I::

       AP=(12. + 40./60.)*PI/12. 
       BP= 28.           *PI/180.
       long_orig=(18. + 40./60.)*PI/12. 
       lat_orig= O.

4. GALACTIC LONG & LAT (System I) --> RIGHT ASCENSION, DECLINATION 

In general, whenever we know the forward transformation (i.e. example 3 above)
we may do the reverse transformation with at most two preliminary call to
coordconv() to calculate the coordinates in system 2 of the pole and origin
which are given in system 1. Often, it is possible to get the needed
coordinates by inspection. In this particular instance::

       APP= 6. *PI/12. 
       BPP=28. *PI/180.
       new_long_orig,new_lat_orig = coordconv(long_orig, lat_orig, AP, BP, 0, 0)
       RA, Dec = coordconv(new_long_orig, new_lat_orig, APP, BPP, lI, bI) 

5. RIGHT ASCENSION, DECLINATION --> GALACTIC LONG & LAT (System II)::

       AP=(12. + 49./60.)*PI/12. 
       BP= 27.4          *PI/180.
       long_orig=(17. + 42.4/60)*PI/12. 
       lat_orig=-(28.0+55./60.)*PI/180.
       lII, bII = coordconv(long_orig, lat_orig, AP, BP, RA, Dec) 

6. GALACTIC LONG & LAT (System II) --> RIGHT ASCENSION, DECLINATION::

       APP, BPP = coordconv(long_orig, lat_orig, AP, BP, 0., Pi/2)
       new_long_orig, new_lat_orig = coordconv(long_orig, lat_orig,
                                               AP, BP, 0., 0.)
       RA, Dec = coordconv(new_long_orig, new_lat_orig, APP, BPP, lII, bII) 

7. RIGHT ASCENSION, DECLINATION --> ECLIPTIC LATITUDE AND LONGITUDE::

       The obliquity of the ecliptic (EPS) depends on the epoch and may
       be obtained from the AMERICAN EPHEMERIS AND NAUTICAL ALMANAC. 
       EPS=23.443*PI/180.
       RA and Dec are for the epoch of observation.
       Ecl.Lat., Ecl.Long = coordconv(0., 0., -Pi/2, Pi/2-EPS, RA, Dec)

Note that the input parameters are partially redundant. For example, if AP, BP,
and long_orig are specified, then there are only two discrete values possible
for lat_orig (except for a few degenerate special cases). The two possible
values of lat_orig may be calculated from::

    sin(lat_orig)=[sin(BP) +/- 2.*cos(BP)**2
                                 *cos(AP-long_orig)
                                 *sqrt{1.+cos(AP-long_orig)**2}]
            / [sin(BP)**2 + {cos(BP)*cos(AP-long_orig)}**2] 
            
If AP, BP, and lat_orig are known, then the two possible values of long_orig
may be calculated from:: 

    cos(AP-long_orig)=[1-sin(lat_orig)*sin(BP)]/[cos(lat_orig)*cos(BP)] 
    
If, instead of long_orig and lat_orig, the longitude of the ascending node is
known in both the old (AN1) and new (AN2) coordinate systems, then long_orig
and lat_orig may be calculated by a preliminary call to COORD::

    long_orig, lat_orig = coordonv(0., 0., AN1-AP, BP, -AN2, 0)
"""
import math

def coordconv(long_orig, lat_orig, AP, BP, A1, B1):
  """
  converts between many lat-long style coordinate systems
  
  @param long_orig : long. of origin of new coordinate system (radians)
  @type  long_orig : float
  
  @param lat_orig : lat. of origin of new coordinate system (radians)
  @type  lat_orig : float
  
  @param AP : long. of pole of  new coordinate system (radians)
  @type  AP : float
  
  @param BP : lat. of pole of  new coordinate system (radians)
  @type  BP : float
  
  @param A1 : long. of point in old coordinate system (radians)
  @type  A1 : float
  
  @param B1 : lat. of point in old coordinate system (radians)
  @type  B1 : float

  @return: (float,float) - long., lat. of point in new coordinate system (rad)
  """
  sin_lat_orig = math.sin(lat_orig)
  cos_lat_orig = math.cos(lat_orig)
  SBP = math.sin(BP)
  CBP = math.cos(BP)
  SB1 = math.sin(B1)
  CB1 = math.cos(B1)

  SB2 = SBP*SB1 + CBP*CB1*math.cos(AP-A1) 
  B2  = math.asin(SB2) 
  CB2 = math.cos(B2)
  SAA = math.sin(AP-A1)*CB1/CB2 
  CAA = (SB1-SB2*SBP)/(CB2*CBP) 
  CBB = sin_lat_orig/CBP 
  SBB = math.sin(AP-long_orig)*cos_lat_orig 
  SA2 = SAA*CBB-CAA*SBB 
  CA2 = CAA*CBB+SAA*SBB 
 
  # There are two formulae for TA2O2 [tan(A2/2)], one usable every- 
  # where except near CA2 = +1, and the other everywhere except near
  # CA2 = -1. 
  if CA2 <= 0:
    TA2O2 = (1.0-CA2)/SA2
  else:
    TA2O2 = SA2/(1.0+CA2)
  A2=2.*math.atan(TA2O2)
  return A2,B2

