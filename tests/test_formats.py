"""
unittest for Astronomy
"""
import unittest
import Astronomy.formats

class testDatesTimes(unittest.TestCase):

  def test_dms_delimited_angle_to_rads(self):
    self.assertEqual(Astronomy.formats.dms_delimited_angle_to_rads('''19d14'33.801860"'''),
                     0.33584886884199222)
    
  def test_hms_delimited_angle_to_rads(self):
    self.assertEqual(Astronomy.formats.hms_delimited_angle_to_rads('00h01m08.621563s'),
                     0.0049903008842279899)
  
  def test_parse_dms_delimited_angle(self):
    self.assertEqual(Astronomy.formats.parse_dms_delimited_angle('''19d14'33.801860"'''),
                     ['19', '14', '33.801860'])
  
  def test_parse_hms_delimited_angle(self):
    self.assertEqual(Astronomy.formats.parse_hms_delimited_angle('00h01m08.621563s'),
                     ['00', '01', '08.621563'])
  
if __name__ == "__main__":
  unittest.main()
