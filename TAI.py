"""
Provides delta T = TAI - UTC

Example::
  In [1]: from Astronomy.TAI import dTA
  In [2]: dTA(2450630)
  Out[2]: 31.0

from http://maia.usno.navy.mil/ser7/tai-utc.dat
"""
from datetime import datetime
from matplotlib.dates import datestr2num, num2date
from numpy import array

def init_table():
  """
  """
  global data
  filename = "/usr/local/RATools/Astronomy/TAI-UTC.dat"
  deltaTfile = open(filename,'r')
  lines = deltaTfile.readlines()
  deltaTfile.close()
  delim1 = "=JD"
  delim2 = "TAI-UTC="
  delim3 = "S +"
  datalist = []
  for line in lines:
    date = num2date(datestr2num(line[:line.index(delim1)]))
    JD = float(line[line.index(delim1)+len(delim1):line.index(delim2)])
    deltaT = float(line[line.index(delim2)+len(delim2):line.index(delim3)])
    datalist.append( (date, JD, deltaT) )
  data = array(datalist, dtype=[('date',datetime),('JD',float),('dTA',float)])

def dTA(jd):
  """
  This does not extrapolate into future years.
  """
  for row in range(len(data)):
    if jd > data['JD'][row-1] and jd < data['JD'][row]:
      return data['dTA'][row]
    else:
      continue
  return data['dTA'][-1]
      
init_table()
  
