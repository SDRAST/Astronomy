import unittest

import Pyro4

from Astronomy.Ephem.serializable_body import SerializableBody
from support.pyro.tests.tests_support_pyro4 import test_case_with_server

def create_body():
    b = SerializableBody()
    b._ra = 1.44399937429555
    b._dec = 0.23617770979137057
    b.info["name"] = "0528+134"
    b.info["flux"] = {
        "X": "2.48",
        "C": "3.41",
        "L": "1.67",
        "K": "3.20",
        "S": "2.97",
        "W": "3.00"
    }
    return b

@Pyro4.expose
class Server(object):

    def get_obj(self):
        return create_body()

    def take_obj(self, body):
        return body

class TestSerializableBody(unittest.TestCase):

    def test_to_dict(self):
        b = create_body()
        d = b.to_dict()
        self.assertTrue("az" not in d)
        self.assertTrue("alt" not in d)
        self.assertTrue("ra" not in d)
        self.assertTrue("dec" not in d)

    def test_from_dict(self):
        b = create_body()
        d = b.to_dict()
        b1 = SerializableBody.from_dict(d)
        self.assertTrue(b._ra == b1._ra)
        self.assertTrue(b._dec == b1._dec)
        self.assertTrue(b.info == b1.info)

    def test_observer_info(self):
        from MonitorControl.Configurations import coordinates
        dss43 = coordinates.DSS(43)
        observer_info = {
            "lat":dss43.lat,
            "lon":dss43.lon,
            "elevation":dss43.elevation,
            "epoch":dss43.epoch,
            "date":dss43.date
        }
        b = create_body()
        b1 = create_body()

        b.compute(dss43)
        b1.observer_info = observer_info
        b1.compute(b1.get_observer())

        self.assertTrue(b1.az == b.az)
        self.assertTrue(b1.alt == b.alt)
        self.assertTrue(b1.ra == b.ra)
        self.assertTrue(b1.dec == b.dec)

        b1_dss43 = b1.get_observer()
        self.assertTrue(b1_dss43.lat == dss43.lat)
        self.assertTrue(b1_dss43.lon == dss43.lon)
        self.assertTrue(b1_dss43.elevation == dss43.elevation)
        self.assertTrue(b1_dss43.epoch == dss43.epoch)
        self.assertTrue(b1_dss43.date == dss43.date)

class TestRegisterWithPyro4(test_case_with_server(Server)):

    def setUp(self):
        self.proxy = Pyro4.Proxy(self.uri)

    def test_send_body(self):
        b = create_body()
        b1 = self.proxy.take_obj(b)

    def test_receive_body(self):
        b = self.proxy.get_obj()

if __name__ == "__main__":
    unittest.main()
