[![DOI](https://zenodo.org/badge/154101466.svg)](https://zenodo.org/badge/latestdoi/154101466)

OpenCropLib
==============

Python Library for simulating physiological processes in plants.
The aim of this project is to implement the energy balance model
developed in Berni et al. 2009 that allows to estimate canopy conductance
from a combination of canopy temperature and basic micro-metheorological
observations.

The different components of the model can also be used in the simulation and
estimation of physical processes such as aerodynamic modeling, calculation
of solar radiation, etc.

# Background

The initial work in developing a simple energy balance model for its application with high-resolution thermal data and infrared thermometers started in my PhD. It was published in [Berni et al., 2009](http://dx.doi.org/10.1016/j.rse.2009.06.018) and was implemented as an Excel VB Macro. In 2017, I started porting the functions from the Excel macro into Python and developing a library that could use `pandas` data frames and efficiently process large weather and sensor data datasets. Some of the functions were adapted from PyETo to deal with `pandas` dataframes.

The library is structured into different modules:
* `radiation` for modeling and estimating various processes such as short/longwave radiation
* `atmosphere` for all the vapor pressure calculations, dew point, air density, etc. 
* `aerodynamic` as the source for different implementations of aerodynamic resistance calculations
* `evapotranspiration` for calculating ETo using FAO56 or ASCE
* `plant` for the calculation of canopy conductance, simulating canopy temperature, or sensible and latent heat calculations. 

# Installation

```
pip install opencroplib
```

# Development

The project is managed with [uv](https://docs.astral.sh/uv/):

```
uv sync          # create .venv with the library and the dev dependencies
uv run pytest    # run the test suite
uv build         # build the sdist and wheel
```

The tests in `tests/test_eto_reference_libraries.py` cross-validate the ETo functions against
[refet](https://github.com/WSWUP/RefET) (ASCE-EWRI 2005) and [pyet](https://github.com/pyet-org/pyet).
`tests/test_radiation.py` checks the solar geometry against the NREL SPA and the clear-sky model against
[pvlib](https://github.com/pvlib/pvlib-python) and GRASS GIS r.sun. These reference libraries are installed as dev
dependencies.

# Examples

Sample notebooks are available in the `sample_notebook` folder. 

# Terms of use and citation
This library is licensed under the MIT license.

You can cite this library as:

Jose A. Jimenez-Berni. (2023). OpenAgriTech/opencroplib: v0.1.0 (v0.1.0). Zenodo. https://doi.org/10.5281/zenodo.8364969


# References
Allen, R. G., Pereira, L. S., Raes, D., and Smith, M. (1998). Crop evapotranspiration: guidelines for computing crop 
water requirements. Roma: Food and Agriculture Organization of the United Nations Available at: 
https://appgeodb.nancy.inra.fr/biljou/pdf/Allen_FAO1998.pdf.

Berni, J. A. J., Zarco-Tejada, P. J., Sepulcre-Cantó, G., Fereres, E., and Villalobos, F. (2009). Mapping canopy 
conductance and CWSI in olive orchards using high resolution thermal remote sensing imagery. Remote Sens. Environ. 113, 
2380–2388. doi:10.1016/j.rse.2009.06.018.

https://github.com/woodcrafty/PyETo

https://github.com/WSWUP/RefET 
