.. Astronomy documentation master file, created by
   sphinx-quickstart on Sat Oct 26 09:26:43 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Single Dish Radio Astronomy Software Tools
==========================================

For an overview of SDRAST and the current status please visit https://sdrast.github.io/.

.. automodapi:: Astronomy

.. raw:: html

  <HR>

.. automodapi:: Astronomy.solar

.. image:: orbit3D.png

+-----------------+-------------------------------------+
|                 |Slowly Varying                       |
+=================+=====================================+
|:math:`\omega`   | ecliptic longitude of the perigee   |
+-----------------+-------------------------------------+
|:math:`\epsilon` | obliquity of the eclipstic          |
+-----------------+-------------------------------------+

+-----------------+-------------------------------------+
|                 | Annually Varying                    |
+=================+=====================================+
|:math:`\alpha`   | right ascension                     |
+-----------------+-------------------------------------+
|:math:`\delta`   | declination                         |
+-----------------+-------------------------------------+
|:math:`M`        | ecliptic longitude of the Sun       |
+-----------------+-------------------------------------+
|:math:`\lambda`  | Sun's anomaly                       |
+-----------------+-------------------------------------+

+------------+-------------------------------------------------+
|            | Projected on Sky toward Sun                     |
+============+=================================================+
|:math:`P`   | tilt of the ecliptic w.r.t. equator             |
+------------+-------------------------------------------------+
|            | tilt of the Sun w.r.t. ecliptic                 |
+------------+-------------------------------------------------+
|:math:`L_0` | longitude of the solar disk center              |
+------------+-------------------------------------------------+
|:math:`B_0` | latitude of the solar disk center               |
+------------+-------------------------------------------------+

.. automodapi:: Astronomy.orbits3D

.. raw:: html

  <HR>
  
.. automodapi:: Astronomy.coordconv

.. raw:: html

  <HR>

.. automodapi:: Astronomy.DSN_coordinates

.. raw:: html

  <HR>

.. automodapi:: Astronomy.Ephem

.. raw:: html

  <HR>

.. automodapi:: Astronomy.formats

.. raw:: html

  <HR>

.. automodapi:: Astronomy.northpolar

.. raw:: html

  <HR>

.. automodapi:: Astronomy.redshift

.. raw:: html

  <HR>

.. automodapi:: Astronomy.ATCA_caldb

.. raw:: html

  <HR>

.. automodapi:: Astronomy.TAI

.. raw:: html

  <HR>

.. automodapi:: Astronomy.SIMBAD

.. toctree::
   :maxdepth: 2



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
