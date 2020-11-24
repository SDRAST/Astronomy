# -*- coding: utf-8 -*-
"""
Set of functions using matplotlib to produce 3-dimensional plots of orbits,
meridia, radial vectors, annotations, etc.

Heliographic coordinates are the latitude and longitude of a feature on the 
Sun's surface. Heliographic latitude is an object's angular distance north or
south of the solar equator; heliographic longitude can be measured east or west
of the central solar meridian, or given in terms of the Carrington rotation 
number.
"""
import pylab as p
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D # matplotlib 0.99

def unit_vector(long_asc_node, inclination, longitude):
  """
  Cartesian components of a unit vector.

  The Cartesian components of a unit vector which lies at a given
  longitude along an inclined plane measured from the line of
  intersection of the planes which is defined by the longitude
  of the ascending node measured from the X-axis in the xy plane.

  @param long_asc_node : longitude of ascending node from X-axis in XY-plane
  @type  long_asc_node : float (degrees)

  @param inclination : tilt of the inclined plane's axis from Z-axis
  @type  inclination : float (degrees)

  @param longitude : longitude measured on inclined plane from intersection
  @type  longitude : float (degrees)
  """
  inc = inclination*p.pi/180
  theta = longitude*p.pi/180
  lan = long_asc_node*p.pi/180

  dx = p.cos(theta)
  dy = p.sin(theta)*p.cos(inc)
  dz = p.sin(theta)*p.sin(inc)
  xp =  dx*p.cos(lan) - dy*p.sin(lan)
  yp =  dx*p.sin(lan) + dy*p.cos(lan)
  return p.array([xp,yp,dz])

def make_axis(ax, maxval, direction, npts):
  """
  Show an axis of a Euclidian space

  @param maxval : radius of vector
  @type  maxval : float

  @param direction : x, y or z
  @type  direction : str

  @param npts : number of points
  @type  npts : int
  """
  if direction.lower() == 'x':
    orbit_vector(ax,maxval,0, 0, 0,npts,"k-")
  elif direction.lower() == 'y':
    orbit_vector(ax,maxval,0, 0,90,npts,"k-")
  elif direction.lower() == 'z':
    orbit_vector(ax,maxval,0,90,90,npts,"k-")
  else:
    print("Direction",direction,"is invalid")
    return False
  return True

def orbit_vector(ax, maxval, long_asc_node, inclination, longitude,
                 npts, line_format):
  """
  Draw a radial vector

  The vector is drawn outwards from the center along a plane which
  intersects the xy plane at the longitude of the ascending node with the
  specified inclination.  The longitude of the vector in this plane is
  measured from the intersection.

  @param maxval : radius
  @type  maxval : float

  @param long_asc_node : longitude (deg) where vector crosses the orbit
  @type  long_asc_node : float

  @param inclination : tilt (deg) of the orbital axis
  @type  inclination : float

  @param longitude : longitude measured on inclined plane from intersection
  @type  longitude : float (degrees)

  @param npts : number of points
  @type  npts : int

  @param line_format : standard matplotlib symbol code, like "k-"
  @type  line_format : str
  """
  V = unit_vector(long_asc_node,inclination,longitude)
  x = V[0]*p.linspace(0,maxval,npts)
  y = V[1]*p.linspace(0,maxval,npts)
  z = V[2]*p.linspace(0,maxval,npts)
  space = ax.plot(x,y,z,line_format,lw=1)
  return space

def make_arc(ax,radius,long_asc_node,inclination,start,end,npts,line_format):
  """
  Draw an arc centered on the origin

  The arc lies on a plane which intersects the xy plane at the longitude of
  the ascending node with the specified inclination.  The start and end
  points are angles measured along the inclined plane from this intersection.

  @param radius : radius of the arc
  @type  radius : float

  @param long_asc_node : longitude (deg) where vector crosses the orbit
  @type  long_asc_node : float

  @param inclination : tilt of the arc's plane (deg)
  @type  inclination : float

  @param start : first longitude measured on inclined plane from intersection
  @type  start : float (degrees)

  @param end : last longitude measured on inclined plane from intersection
  @type  end : float (degrees)

  @param npts : number of points
  @type  npts : int

  @param line_format : standard matplotlib symbol code, like "k-"
  @type  line_format : str
  """
  thetas = p.linspace(start,end,npts)
  xp, yp, z = [], [], []
  for theta in thetas:
    V = radius*unit_vector(long_asc_node,inclination,theta)
    xp.append(V[0])
    yp.append(V[1])
    z.append(V[2])
  space = ax.plot(xp,yp,z,line_format,lw=1)
  return True

def make_orbit(ax, radius, long_asc_node, inclination, npts, line_format):
  """
  Draw a full circle around the origin

  @param radius : radius of the arc
  @type  radius : float

  @param inclination : tilt of the arc's plane (deg)
  @type  inclination : float

  @param npts : number of points
  @type  npts : int

  @param line_format : standard matplotlib symbol code, like "k-"
  @type  line_format : str
  """
  space = make_arc(ax,radius,long_asc_node,inclination,0,359,npts,line_format)
  return space

def orbit_text(ax,radius,long_asc_node,inclination,longitude,text):
  """
  Position text using orbital coordinates.
  
  @param radius : radial distance of text
  @type  radius : float (degrees)

  @param long_asc_node : longitude (deg) where vector crosses the orbit
  @type  long_asc_node : float

  @param inclination : tilt of the inclined plane's axis
  @type  inclination : float (degrees)

  @param longitude : longitude along the inclined plane from the intersection
  @type  longitude : float (degrees)

  @param text :
  @type  text : str
  """
  V = radius*unit_vector(long_asc_node,inclination,longitude)
  ax.text3D(V[0],V[1],V[2],text)

def show_celestial(fig, obliquity, long_perigee, anomaly, Title=None):
  """
  Show an orbit in the celestial coordinate system

  @param fig : figure number
  @type  fig : int

  @param obliquity : tilt of the orbit's axis w.r.t. celestial north
  @type  obliquity : float

  @param long_perigee : longitude of the perigee
  @type  long_perigee : float

  @param anomaly : position of object along orbit from perigee
  @type  anomaly : float

  @return: figure Axis instance
  """
  fig = plt.figure(fig,figsize=(8, 6))

  ax = Axes3D(fig)

  # make axes
  space = make_orbit(ax, 10, 0, 0, 360, "k-")
  make_orbit(ax,10,0,obliquity,360,"b-")
  make_axis(ax,11,'x',12)
  make_axis(ax,11,'y',12)
  make_axis(ax,11,'z',12)
  ax.text3D(11,0,0,r"$\alpha=0^h, \lambda=0^{\circ}$")
  ax.text3D(0,11,0,r"$\alpha=6^h$")
  ax.text3D(0,0,11,"N")

  # show the perigee
  orbit_vector(ax,11, 0, obliquity, long_perigee, 12, "b-")
  orbit_text(ax,  11, 0, obliquity, long_perigee, "perigee")

  # show the Sun
  orbit_vector(ax,11,         0, obliquity, long_perigee+anomaly,  12, "y-")
  sun_pos = 10*unit_vector(0, obliquity, long_perigee+anomaly)
  ax.plot3D([sun_pos[0]], [sun_pos[1]], [sun_pos[2]],'yo')
  orbit_text(ax,  11, long_perigee+anomaly, obliquity, 0, "Sun")
  make_arc(ax,10, long_perigee+anomaly,  90, 0, obliquity,20, "y:")
  orbit_vector(ax,10,  0, 0, long_perigee+anomaly, 20, "y:")

  # show the angles
  make_arc(ax,  2,   0, obliquity,   0,         long_perigee, 300, "r-")
  orbit_text(ax,2.2, long_perigee/2, obliquity, 0,            r"$\omega$")

  make_arc(ax,2.5,   0, obliquity,           long_perigee, long_perigee+anomaly,
           300,"g-")
  orbit_text(ax,3,   long_perigee+anomaly/2, obliquity,    0, r"$M$")

  solar_longitude = p.fmod(long_perigee+anomaly,360)
  make_arc(ax,  4,0, obliquity, 0, solar_longitude,   300, "c-")
  orbit_text(ax,4.2, obliquity, 0, solar_longitude/2, r"$\lambda$")

  orbit_vector(ax,11, 0, 90+obliquity, 90, 12,"b-")
  orbit_text(ax,  11, 0, 90+obliquity, 90, r"$\beta=90^{\circ}$")

  make_arc(ax, 6,     90, 90, 90,  90+obliquity, 300, "b-")
  orbit_text(ax,5.3,  90,90-obliquity/2, 90, r"$\epsilon$")

  # show arc from 6 h to 18 hcthrough ecliptic pole in XZ plane
  make_arc(ax,10,     90, 90,  0, 180,            300, "k--")

  make_arc(ax,10, long_perigee, 90, 180-obliquity, 180,           300, 'b:')
  make_arc(ax,10, long_perigee, 90,    -obliquity, 180-obliquity, 180, 'b--')
  ax.clabel(space, fontsize=9, inline=1)
  ax.set_xlabel("x")
  ax.set_ylabel("y")
  ax.set_zlabel("z")
  if Title:
    ax.set_title(Title)
  return ax

def show_heliographic(fig, inclin_sun, long_asc_node, anomaly, sun_tilt,
                      Title=None):
  """
  Show a celestial coordinate system in heliographic coordinates

  @param fig : figure number
  @type  fig : int

  @param inclin_sun : tilt of Sun's polar axis w.r.t. ecliptic north
  @type  inclin_sun : float

  @param long_asc_node : ecliptic long. where Sun's equator rises over ecliptic
  @type  long_asc_node : float

  @param sun_tilt : projection of Sun's polar axis on plane of the sky
  @type  sun_tilt : float

  @return: Axis() instance
  """
  fig = plt.figure(fig)

  ax = Axes3D(fig)
  make_axis(ax,11,'x',12)
  make_axis(ax,11,'y',12)
  make_axis(ax,11,'z',12)
  ax.set_xlabel("x'")
  ax.set_ylabel("y'")
  ax.set_zlabel("z'")
  space = make_orbit(ax,10,0,0,360,"k-")
  ax.text3D(11,0,0,r"$\lambda=0^{\circ}$")
  ax.text3D(0,11,0,r"$\lambda=90^{\circ}$")
  ax.text3D(0,0,11.5,r"$B=90^{\circ}$")

  make_orbit(ax,10, long_asc_node, inclin_sun, 360,"b-")
  orbit_vector(ax,11, 0, 0, long_asc_node, 12, "b-")
  orbit_text(ax,11, long_asc_node, 0, 0, r"$\lambda=75^{\circ}$")

  sun_pos = 10*unit_vector(long_asc_node, inclin_sun, anomaly)
  ax.plot3D([sun_pos[0]], [sun_pos[1]], [sun_pos[2]],'yo')
  orbit_text(ax,  11, long_asc_node, inclin_sun, anomaly, "Sun")
  orbit_vector(ax,11,  0, -inclin_sun, long_asc_node+anomaly,  12, "y-")
  orbit_vector(ax,10,  0, 0, long_asc_node+anomaly, 20, "y:")

  # The pole of the heliographic coordinates lies on a great circle that
  # passes through the ecliptic pole and ecliptic equator at the longitude
  # of the ascending node and has an obliquity of 90 deg
  orbit_vector(ax,11, long_asc_node, 90+inclin_sun, 90, 12,"b-")
  orbit_text(ax,  10.8,            90, 90+inclin_sun, 90, r"$\beta=90^{\circ}$")

  make_arc(ax,10, long_asc_node,     90+inclin_sun, 0,           90,
           200, "b--")
  orbit_vector(ax,11, long_asc_node, inclin_sun, 90, 12, "b")

  make_arc(ax,10, long_asc_node+ 90, 90,            +inclin_sun, 90+inclin_sun,
           200, "b--")
  make_arc(ax,10, long_asc_node+180, 90-inclin_sun, 0,           90,
           200, "b--")
  make_arc(ax,10, long_asc_node+270, 90,            -inclin_sun, 90-inclin_sun,
           200, "b--")

  make_arc(ax,10, long_asc_node+anomaly,    90+sun_tilt, 0, 180, 200, "r-")
  make_arc(ax,10, long_asc_node+anomaly,    90         , 0, 180, 200, "k-")

  make_arc(ax,8,     long_asc_node+90, 90, 90, 90+inclin_sun,10,'m-')
  orbit_text(ax,8.5, long_asc_node+90, 90, 90+inclin_sun/2, r"$\iota$")
  if Title:
    ax.set_title(Title)
  return ax
