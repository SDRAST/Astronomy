"""
Module serializable_class

Support for making pyephem classes serializable for Pyro

Examples:

This creates a serializable FixedBody class
    
.. code-block:: python

    class SerializableBody(ephem.FixedBody, SerializedClass):
        def __init__(self, info=None, name="", observer_info=None):
            ephem.FixedBody.__init__(self) 
            self._class = '\x02'
            
    SerializableBody.register_with_Pyro5()

    body = SerializableBody(
        name="0521-365",
        info={"flux": {"C": "6.76",
                       "K": "3.50",
                       "L": "16",
                       "S": "12.50",
                       "W": "3.55",
                       "X": "5.19"}})
        body._class = ephem.FixedBody._class
        body._ra = 1.4092068106703093
        body._dec = -0.6363220818328413
        body.compute()
        
This creates a serializable Planet class (`e.g.` Sun, Moon, Venus, ...)

.. code-block:: python

    import ephem   
    from Astronomy.Ephem.serializable_class import SerializedClass
    class SerializableVenus(ephem.Venus, SerializedClass): 
        def __init__(self, info=None, name="", observer_info=None): 
            ephem.Venus.__init__(self) 
            self.compute() 
            self.info = info 
            self.observer_info = info 
            self._ra = self.a_ra 
            self._dec = self.a_dec 

    SerializableVenus.register_with_Pyro5()

    venus = SerializableVenus()
    from Astronomy.DSN_coordinates import DSS
    gavrt = DSS(28)
    gavrt.date = ephem.now()
    venus.compute(gavrt) 
    venus.to_dict()
    {'info': None,             'observer_info': None,       'name': 'Venus',
     '__class__': 'SerializableVenus',
     '_ra': 3.4982651747476,   '_dec': -0.116797529229786,   
     'ra':  3.502862013642421, 'dec':  -0.11866434539848866,
     'az': 3.7932510375976562, 'alt': 0.7079898715019226, 'el': 0.7079898715019}
"""
import ephem

#__all__ = ["SerializerBase"]


class SerializedClass(object):
    """
    Superclass to convert ``pyephem`` ``Body`` to and from a Python dictionary.

    Attributes:
    
        info (dict): Additional information about source.
        name (str): Name of source
        observer_info (dict): dictionary containing information about
            an ephem.Observer object.
    """
    def __init__(self, info=None, name="", observer_info=None):
        """
        initialize object as a SerializedClass
        """
        if info is None:
            self.info = {}
        self.name = name
        self.observer_info = observer_info

    def to_dict(self):
        """
        Dump object to a dictionary.
        """
        return_dict = {
            "info": self.info,
            "observer_info":self.observer_info,
            "name":self.name,
            "_ra": float(self._ra),
            "_dec": float(self._dec),
            "__class__":self.__class__.__name__
        }
        if hasattr(self, "ra"):
            return_dict["ra"] = float(self.ra)
            return_dict["dec"] = float(self.dec)

        if hasattr(self, "az"):
            return_dict["az"] = float(self.az)
            return_dict["alt"] = float(self.alt)
            return_dict["el"] = float(self.alt)

        return return_dict

    @classmethod
    def from_dict(cls, src_dict):
        """
        Class method for recreating SerializableBody objects from a dictionary.

        Args:
            src_dict (dict): dictionary returned from ``to_dict`` instance method.
        Returns:
            SerializableBody
        """
        obj = cls()
        obj._ra = src_dict["_ra"]
        obj._dec = src_dict["_dec"]
        obj.info = src_dict["info"]
        obj.observer_info = src_dict["observer_info"]
        obj.name = src_dict["name"]
        return obj

    @classmethod
    def register_with_Pyro5(cls):
        """
        Make sure that we can serialize and deserialize this class when
        sending SerializableBody between Pyro5 clients and servers.
        """
        import Pyro5.serializers as pyro

        def from_dict(name, src_dict):
            return cls.from_dict(src_dict)

        pyro.SerializerBase.register_class_to_dict(cls, cls.to_dict)
        try:
          pyro.SerializerBase.register_dict_to_class(cls.__name__, from_dict)
        except AttributeError:
          pass
          
    def get_observer(self):
        """
        Return an ephem.Observer object from the body's observer_info dictionary.
        """
        if self.observer_info is not None:
            observer = ephem.Observer()
            observer.lon = self.observer_info["lon"]
            observer.lat = self.observer_info["lat"]
            observer.elevation = self.observer_info["elevation"]
            observer.date = self.observer_info["date"]
            observer.epoch = self.observer_info["epoch"]
            return observer

