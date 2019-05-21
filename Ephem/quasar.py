import logging

import ephem

from . import Jnames, Bname_dict, cat_3C_dict, cal_dict, Bnames
#from support.Ephem import EphemException

module_logger = logging.getLogger(__name__)

class Quasar(ephem.FixedBody):
  """
  ephem.FixedBody() with "J" or "B" name
  """
  def __init__(self,name):
    """
    Create an instance of a Quasar()

    The preferred name form is like 'J2305+4707'.  However, 'B2303+46'
    is also exceptable, as is '2305+4707'.  In the latter case, it must
    be a Julian epoch name.

    The name may also be given as "3C84" or "3C 84"

    Notes
    =====
    Instance methods and attributes
    -------------------------------
    flux -     radio flux in janskys at a specified frequency
    flux_ref - source for radio flux data
    freq -     frequency for which flux was calculated
    date -     date for which flux was calculated, if based on Michigan data
    Jname -    Julian epoch name
    Bname -    Besselian epoch name

    Inherited methods and attributes
    --------------------------------
    All methods and attributes of ephem.FixedBody()
    name - Julian epoch name without leading "J"
    """
    super(Quasar,self).__init__()
    if name[0] == "J":
      self.Jname = self.name
      self.name = name[1:]
      try:
        self.Bname = "B"+Bnames[self.name]
      except KeyError:
        raise(EphemException(name,"unknown J-name"))
    elif name[0] == "B":
      self.Bname = name
      try:
        self.name = Bname_dict[name[1:]]
      except KeyError:
        raise(EphemException(name,"unknown B-name"))
      self.Jname = "J"+self.name
    elif name[:2] == "3C":
      name = name.replace(" ","")
      try:
        self.name = cat_3C_dict[name]
      except KeyError:
        raise(EphemException(name,"unknown 3C source"))
      self.Jname = "J"+self.name
      self.Bname = "B"+Bnames[self.name]
    else:
      try:
        Jnames.index("J"+name)
      except:
        raise(EphemException(name,"unknown quasar"))
      self.name = name
      self.Jname = "J"+self.name
      self.Bname = "B"+Bnames[self.name]
    quasar_data = cal_dict[self.name]
    self._ra = hours(str(quasar_data['ra']))
    self._dec = degrees(str(quasar_data['dec']))
    if diag:
      print self._ra, self._dec
    self._epoch = J2000
    self._class = "Q"
    # to initialize ra, dec to something
    self.compute("2000/1/1 00:00:00")
    self.freq = None
    self.flux = None
    self.flux_ref = None
    self.date = None
    self.pressure = 0

  def interpolate_flux(self,freq,date=None):
    """
    Flux of source at given frequency in GHz.

    The Michigan data are interpolated if possible, that is, the source is
    in the data base and the frequency requested is between 4 and 15 GHz.
    Otherwise the VLA data are used.
    If the date is not given, it is assumed to be now.
    """
    self.freq = freq
    if self.freq > 4 and self.freq < 15:
      if date == None:
        self.date = datetime.datetime.now()
      else:
        self.date = date
      try:
        # in Michigan data
        michigan.Bnames.index(self.Bname[1:])
        self.flux = michigan.polate_flux(self.Bname[1:],
                                         date2num(self.date),
                                         freq)
        if diag:
          print "Michigan flux is",self.flux
        self.flux_ref = "Michigan"
        return self.flux
      except ValueError:
        # Not in the Michigan catalog or no date given
        if diag:
          print self.name,"=",self.Bname,"is not in the Michigan catalog"
    else:
      print "Outside Michigan frequency range"
    # try VLA cals
    try:
      cal_data = vla_cal.get_cal_data(self.Jname[1:])
      freqs,fluxes = vla_cal.get_flux_data(cal_data)
      self.flux = vla_cal.interpolate_flux(freqs,fluxes,freq)
      self.flux_ref = "VLA"
      if diag:
        print ref,"flux is",self.flux
      return self.flux
    except ValueError:
      self.flux = None
      self.flux_ref = None
      return None
