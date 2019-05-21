.. Opencroplib documentation master file, created by
   sphinx-quickstart
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Opencroplib: easy calculation of evapotranspiration processes with Python
=========================================================================

Opencroplib is a collection of Python libraries for simulating evapostranspiration processes. It includes code for extremely precise ephemeris calculations, and more.


Prerequisites for use
---------------------

Opencroplib requires Python together with some basic libraries.

Python anaconda is recommended, as it includes all the required libraries::

    Python 3.4.2 (default, Oct  8 2014, 14:38:51)
    [GCC 4.9.1] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>>


If you need to, you can download Python from the `Python.org download page <http://python.org/download/>`_.

Examples
========


Estimate of clear-sky radiation
-------------------------------

Once you calculate azimuth and altitude of the sun, you can predict the direct irradiation from the sun using Pysolar. ``get_radiation_direct()`` returns a value in watts per square meter. As of version 0.7, the function is *not* smart enough to return zeros at night. It does account for the scattering of light by the atmosphere, though it uses an atmospheric model based on data taken in the United States.::

   latitude_deg = 42.206 # positive in the northern hemisphere
   longitude_deg = -71.382 # negative reckoning west from prime meridian in Greenwich, England
   date = datetime.datetime(2007, 2, 18, 15, 13, 1, 130320, tzinfo=datetime.timezone.utc)
   altitude_deg = get_altitude(latitude_deg, longitude_deg, date)
   radiation.get_radiation_direct(date, altitude_deg)

Results in

   909.582292149944


References
==========

`Abstract <https://www.sciencedirect.com/science/article/pii/S0034425709002090>`_ | `PDF Download <http://quantalab.ias.csic.es/pdf/paper_berni_RSE_2009.pdf>`_ | Berni, J. A. J., Zarco-Tejada, P. J., Sepulcre-Cantó, G., Fereres, E., and Villalobos, F. (2009). Mapping canopy conductance and CWSI in olive orchards using high resolution thermal remote sensing imagery. Remote Sens. Environ. 113, 2380–2388. doi:10.1016/j.rse.2009.06.018.



Source Code Repository
======================

The source code for Opencrop lib can be found in :
`https://github.com/OpenAgriTech/opencroplib <https://github.com/OpenAgriTech/opencroplib>`_

Contributors
============

Many people have contributed to Opencroplib since its inception, initially as Excel macros

Thanks to Francisco Villalobos, Luca Testi and Elias Fereres for their ideas, contributions of code, bugfixes, documentation, and general encouragement.

.. toctree::
   :maxdepth: 2

