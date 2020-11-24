"""
provides a class for serializing planets
"""

import ephem

from .serializable_class import SerializedClass
planets = ["MERCURY", "VENUS", "MARS", "JUPITER", "SATURN", "URANUS"," NEPTUNE",
           "PLUTO",   "SUN",   "MOON"]
class SerializablePlanet(ephem.Planet, SerializedClass): 
    """
    serialized subclass of ``Planet``
    """
    def __init__(self, name="", info=None, observer_info=None):
        """
        This initializes a planet implicitly for ``ephem.now()``
        """
        self.__planet__ = planets.index(name.upper())
        if not name:
            raise RuntimeError("A solar system body name is needed")
        if name.capitalize() == "Venus":
            ephem.Venus.__init__(self)
        elif name.capitalize() == "Sun":
            ephem.Sun.__init__(self)
        elif name.capitalize() == "Moon":
            ephem.Moon.__init__(self)
        elif name.capitalize() == "Mars":
            ephem.Mars.__init__(self)
        elif name.capitalize() == "Jupiter":
            ephem.Jupiter.__init__(self)
        elif name.capitalize() == "Mercury":
            ephem.Mercury.__init__(self)
        elif name.capitalize() == "Uranus":
            ephem.Uranus.__init__(self)
        elif name.capitalize() == "Saturn":
            ephem.Saturn.__init__(self)
        elif name.capitalize() == "Pluto":
            ephem.Pluto.__init__(self)
        else:
            raise RuntimeError("%s is not a known Planet class")
        self.compute() 
        self.info = info 
        self.observer_info = info 
        self._ra = self.a_ra 
        self._dec = self.a_dec 

SerializablePlanet.register_with_Pyro5()
