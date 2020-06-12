Astronomy
=========

Classes and functions for astronomical calculations


Astronomy calculations, such as source positions in various coordinate
systems, input and output of astronomical data in standard formats, and also 
source data lookup from the SIMBAD database, are generally handled by 
`astropy <https://www.astropy.org/>`_.

The earliest version of this package was based on the FORTRAN program DOPSET
by Manchester and Gordon [MG1970]_ for Doppler shift with respect to the LSR [#]_.  
For solar system bodies, `PyEphem
<https://rhodesmill.org/pyephem/>`_ support was later added. For convenience, 
a module ``Astronomy.Ephem``
provides extensions to ``ephem``, with
subclass ``DSS`` of ``ephem.Observer`` for DSN and ESA stations,
and subclasses ``Pulsar`` and ``Quasar`` of ``ephem.FixedBody``.

The ``Astronomy.Ephem`` submodule also contains an extension to the PyEphem base 
data type that can be serialized, so as to be sent and received over network 
connections.

``Astronomy`` also has modules for redshift calculations,
solar coordinates, :math:`\Delta T=TAI-UTC`, and more.

**Links**

`Project <https://github.com/SDRAST/Astronomy/>`_ (repository).

`Source code documentation <https://sdrast.github.io/Astronomy/>`_ (API).

.. rubric:: Footnotes

.. [#] For one-time Doppler calculation, the `Online Dopset Tool <http://www.vla.nrao.edu/astro/guides/dopset/>`_.

.. rubric:: Bibliography

.. [MG1970] Manchester, R. N. & Gordon, M. *DOPSET: A Computer Program to Calculate Doppler-Corrected Reference Frequencies for Spectral Line and Pulsar Observations*, Technical Report, National Radio Astronomy Observatory, 1970.

