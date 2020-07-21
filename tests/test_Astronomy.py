"""
unittest for Astronomy
"""
import unittest
import datetime
import Astronomy

class testDatesTimes(unittest.TestCase):

  def test_B_epoch_to_J(self):
    self.assertEqual(Astronomy.B_epoch_to_J('''00h02m29.056400s''',
                                            '''54d11'43.187000"''',
                                            'decimal'),
                     (0.084544, 54.4736))
                     
    self.assertEqual(Astronomy.B_epoch_to_J('''00h02m29.056400s''',
                                            '''54d11'43.187000"''',
                                            'formatted'),
                     ['00h05m04.359s', '+54d28m25.056s'])
    
    self.assertEqual(Astronomy.B_epoch_to_J('''00h02m29.056400s''',
                                            '''54d11'43.187000"'''),
                     ([0, 5, 4.359], [54, 28, 25.056]))

    self.assertEqual(Astronomy.B_epoch_to_J('''23h58m34.865400s''',
                                  '''18d57'51.753000"'''),
                     ([0, 1, 8.6169], [19, 14, 33.9321]))
    self.assertEqual(Astronomy.B_epoch_to_J('''00h00m48.4200s''',
                                            '''-17d43'54.000"''',
                                            'formatted'),
                     [u'00h03m21.9921s', u'-17d27m11.6511s'])

  def test_J_epoch_to_B(self):
    self.assertEqual(Astronomy.B_epoch_to_J('''00h02m29.056400s''',
                                            '''54d11'43.187000"''',
                                            'formatted'),
                     [u'00h05m04.359s', u'+54d28m25.056s'])
    
    self.assertEqual(Astronomy.J_epoch_to_B('''00h05m04.359s''',
                                            '''+54d28'25.056"''', 
                                            'formatted'),
                     [u'00h02m29.0564s', u'+54d11m43.1869s'])
  
  def test_obs_ra_to_cirs_ra(self):
    """
    This tests that a source with RA = LST and dec = latitude is at the zenith
    
    from
    http://reionization.org/wp-content/uploads/2013/03/HERA_Memo46_lst2ra.html
    """
    from astropy import units
    from astropy.coordinates import Angle, EarthLocation, SkyCoord
    from astropy.time import Time
    # the sidereal time of observation
    mjd = 55780.1
    latitude = Angle('-26d42m11.94986s')
    longitude = Angle('116d40m14.93485s')
    obs_time = Time(mjd, format='mjd', location = (longitude, latitude))
    lst_apparent = obs_time.sidereal_time('apparent')
    print('lst_apparent type is %s' % type(lst_apparent.hour))
    print('lst_apparent is %f' % lst_apparent.hour)
    self.assertEqual(lst_apparent.to_string(), '7h11m46.2716s')
    # convert observed RA (=LST) to CIRS RA
    cirs_ra = Astronomy.obs_ra_to_cirs_ra(lst_apparent.hour,
                                          Time(mjd, format='mjd'),
                                          longitude, latitude)
    print('CIRS RA - observed RA (= LST): %f' % (cirs_ra-lst_apparent.hour))
    self.assertAlmostEqual(cirs_ra, 7.1859741)
    # check that RA=LST, dec=lat is at zenith
    loc_obj = EarthLocation.from_geodetic(lon=longitude, lat=latitude)
    obs_zenith_coord = SkyCoord(ra=cirs_ra, dec=latitude, frame='cirs',
                                unit=("hour","deg"),
                                obstime=Time(mjd, format='mjd'),
                                location = loc_obj)
    obs_zenith_altaz = obs_zenith_coord.transform_to('altaz')
    obs_ZD = (obs_zenith_altaz.alt - Angle('90d')).to_string(unit=units.degree, 
                                                          sep=('deg', 'm', 's'))
    self.assertEqual(obs_ZD, '-0deg00m00.5114s')
  
  def test_HaDec_to_AzEl(self):
    Az,El = Astronomy.HaDec_to_AzEl( 0,0,30)
    self.assertAlmostEqual(Az,-180.0)
    self.assertAlmostEqual(El,  60.0)
    Az,El = Astronomy.HaDec_to_AzEl(-6,0,30)
    self.assertAlmostEqual(Az, 90.0)
    self.assertAlmostEqual(El,  0.0)
    Az,El = Astronomy.HaDec_to_AzEl( 6,0,30)
    self.assertAlmostEqual(Az, -90.0)
    self.assertAlmostEqual(El ,  0.0)

  def test_AzEl_to_HaDec(self):
    Ha,Dec = Astronomy.AzEl_to_HaDec(-180,60,30)
    self.assertAlmostEqual(Ha, 0)
    self.assertAlmostEqual(Dec, 0)
if __name__ == "__main__":
  unittest.main()
