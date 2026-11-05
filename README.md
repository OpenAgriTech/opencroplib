[![DOI](https://zenodo.org/badge/154101466.svg)](https://zenodo.org/badge/latestdoi/154101466)
[![tests](https://github.com/OpenAgriTech/opencroplib/actions/workflows/tests.yml/badge.svg)](https://github.com/OpenAgriTech/opencroplib/actions/workflows/tests.yml)

OpenCropLib
===========

Python library for simulating physiological and micrometeorological processes in crops.

Its core is the canopy energy balance model of [Berni et al. (2009)](https://doi.org/10.1016/j.rse.2009.06.018),
which estimates **canopy conductance** from canopy temperature (thermal imagery or infrared thermometers) and basic
weather observations, and the inverse problem, **simulating canopy temperature** for a given conductance. The
building blocks of the model are exposed as independent functions: vapour pressure and psychrometrics, aerodynamic
resistance, solar geometry, clear-sky and longwave radiation, and reference evapotranspiration (FAO-56 and ASCE).

All functions accept scalars, NumPy arrays and pandas Series, so they can be applied directly to large weather and
sensor datasets.

## Contents

- [Background](#background)
- [Installation](#installation)
- [Modules](#modules)
- [Units](#units)
- [Usage examples](#usage-examples)
- [Validation](#validation)
- [Upgrading from 0.1.x](#upgrading-from-01x)
- [Development](#development)
- [Citation](#citation)
- [References](#references)

## Background

The simple energy balance model for high-resolution thermal data and infrared thermometers was developed during my
PhD and published in [Berni et al. (2009)](https://doi.org/10.1016/j.rse.2009.06.018), where it was implemented as an
Excel VBA macro. In 2017 I started porting it to Python as a library that works with `pandas` and can process large
weather and sensor datasets. Some of the FAO-56 functions were adapted from [PyETo](https://github.com/woodcrafty/PyETo)
and the atmospheric functions from [meteolib](http://python.hydrology-amsterdam.nl/) (Waterloo, 2014).

## Installation

```
pip install opencroplib
```

To install the latest development version from GitHub:

```
pip install "git+https://github.com/OpenAgriTech/opencroplib"
```

Requires Python 3.11 or newer, NumPy, pandas and pytz.

## Modules

| Module | Contents |
|---|---|
| `opencroplib.atmophere` | Saturation and actual vapour pressure, vapour pressure deficit, slope of the saturation curve, psychrometric constant, latent heat of vaporisation, specific heat, air density, volumetric heat capacity, potential temperature, dew point. (The module name is spelled `atmophere`.) |
| `opencroplib.aerodynamic` | Aerodynamic resistance (FAO-56 neutral profile, and Viney (1991) for unstable conditions), roughness length and displacement height (Raupach 1994), wind speed conversion between heights and surfaces. |
| `opencroplib.radiation` | Solar declination, elevation and azimuth; extraterrestrial irradiance; ESRA clear-sky beam, diffuse and global irradiance; downwelling longwave radiation (Brutsaert, Swinbank, Berdahl–Martin, Idso, FAO-56 Brunt-type, with cloud correction); net radiation and available energy; FAO-56 daily and hourly radiation terms; Planck functions for thermal sensors. |
| `opencroplib.plant` | Canopy conductance from canopy temperature (Berni et al. 2009), canopy temperature simulation for a given conductance, sensible and latent heat fluxes. |
| `opencroplib.evapotranspiration` | Reference evapotranspiration: FAO-56 hourly and daily, and ASCE-EWRI hourly (short reference). |
| `opencroplib.time_utils` | Day of the year, equation of time and local solar time. |

## Units

The library uses the following units unless a docstring says otherwise:

| Quantity | Unit |
|---|---|
| Temperature | °C |
| Vapour pressure, VPD, air pressure | Pa |
| Psychrometric constant, slope of the saturation curve | Pa K⁻¹ |
| Relative humidity | % |
| Radiation and energy fluxes (hourly and instantaneous) | W m⁻² |
| Radiation in the FAO-56 *daily* functions | MJ m⁻² day⁻¹ |
| Aerodynamic and canopy resistance | s m⁻¹ |
| Canopy conductance | mm s⁻¹ |
| Wind speed | m s⁻¹ |
| Heights | m |
| Latitude, longitude | degrees (+N, +E) |
| Solar angles returned by the functions | radians |
| Reference evapotranspiration | mm h⁻¹ (hourly), mm day⁻¹ (daily) |

The FAO-56 daily helpers (`net_out_lw_daily`, `daily_potential_toa`, ...) follow FAO-56 and expect vapour pressure in
kPa where indicated.

## Usage examples

### Vapour pressure and psychrometrics

```python
from opencroplib.atmophere import es_calc, vpd_calc, psychrometric_gamma

es_calc(25.0)                               # saturation vapour pressure [Pa]
# 3166.82...
vpd_calc([20.0, 30.0], [60.0, 30.0])        # vapour pressure deficit [Pa]
# array([ 934.83, 2969.91])
psychrometric_gamma(25.0, 50.0, 101325.0)   # psychrometric constant [Pa K-1]
# 67.57...
```

### Daily reference evapotranspiration (FAO-56)

Daily ET₀ from station data with mean temperature, relative humidity, wind speed at 2 m and solar radiation. Net
radiation is estimated with FAO-56 Eqs. 21, 37, 38, 39 and 40.

```python
import pandas as pd
from opencroplib.atmophere import ea_calc
from opencroplib.radiation import (daily_potential_toa, daily_clear_sky_irradiance, net_out_lw_daily,
                                   net_in_sol_rad, net_rad)
from opencroplib.evapotranspiration import et0_FAO56_daily

# Three days from the IAS-CSIC weather station (Cordoba, Spain)
weather = pd.DataFrame(
    {"tmean": [11.7, 12.7, 15.8], "tmax": [20.1, 18.9, 20.8], "tmin": [3.3, 6.6, 10.9],
     "rh": [60.9, 74.1, 64.8], "wind": [1.7, 1.3, 2.0], "rs": [22.4, 18.2, 24.5]},   # rs in MJ m-2 day-1
    index=pd.date_range("2000-03-29", periods=3, freq="D"))
latitude, altitude, pressure = 37.857, 117.0, 100500.0

ra = daily_potential_toa(latitude, weather.index.day_of_year)       # extraterrestrial radiation
rso = daily_clear_sky_irradiance(altitude, ra)                       # clear-sky radiation
ea_kpa = ea_calc(weather.tmean, weather.rh) / 1000.0                 # actual vapour pressure [kPa]
rnl = net_out_lw_daily(weather.tmin, weather.tmax, weather.rs, rso, ea_kpa)
weather["rn"] = net_rad(net_in_sol_rad(weather.rs), rnl)              # MJ m-2 day-1
weather["et0"] = et0_FAO56_daily(weather.tmean, weather.rh, pressure, weather.wind, weather.rs, weather.rn)
weather[["rn", "et0"]]
#                rn   et0
# 2000-03-29  10.96  3.06
# 2000-03-30   9.75  2.53
# 2000-03-31  12.30  3.70
```

> **Note.** `et0_FAO56_daily` currently computes the vapour pressure deficit from the daily *mean* temperature and
> relative humidity. Because the saturation vapour pressure is a convex function of temperature, this underestimates
> the daily VPD: over 26 years of IAS-CSIC data it gives ET₀ about 4 % lower than the station's reported value.
> Supporting the other FAO-56 humidity inputs (Tmax/Tmin with RH, dew point, RHmax/RHmin, measured VPD) is planned.

If the wind is measured at another height over grass, pass it with `z_wind`; it is converted to 2 m with FAO-56 Eq. 47.

### Hourly reference evapotranspiration (FAO-56 and ASCE)

Hourly functions take radiation in W m⁻². When `soil_heat_flux` is not given, G = 0.1 Rn during the day
(`solar_radiation > 0`) and 0.5 Rn at night.

```python
import numpy as np
from opencroplib.evapotranspiration import et0_FAO56, et0_ASCE

t = np.array([18.0, 32.0, 24.0])          # air temperature [C]
rh = np.array([60.0, 30.0, 50.0])         # relative humidity [%]
u10 = np.array([2.0, 4.0, 1.5])           # wind speed at 10 m over grass [m s-1]
rs = np.array([0.0, 850.0, 0.0])          # solar radiation [W m-2]
rn = np.array([-60.0, 560.0, -40.0])      # net radiation [W m-2]

et0_FAO56(t, rh, 101325.0, u10, rs, rn, z_wind=10.0)   # mm h-1
# array([0.02, 0.69, 0.03])
et0_ASCE(t, rh, 101325.0, u10, rs, rn, z_wind=10.0)    # mm h-1 (Cd = 0.24 day, 0.96 night)
# array([0.02, 0.73, 0.03])
```

### Canopy conductance from canopy temperature (Berni et al. 2009)

Canopy conductance is obtained by inverting the canopy energy balance (Berni et al. 2009, Eq. 4) with the measured
canopy temperature, the available energy and the aerodynamic resistance. Here the aerodynamic resistance of an
orchard (4 m tall, reference height 6 m) is computed with the stability correction of Viney (1991), from a wind speed
measured at 2 m over grass.

```python
import numpy as np
from opencroplib.aerodynamic import aerodynamic_resistance_viney
from opencroplib.radiation import longwave_downwelling, available_energy
from opencroplib.plant import canopy_conductance

air_t, rh, pressure = 32.0, 30.0, 100500.0       # C, %, Pa
wind, shortwave = 2.5, 850.0                     # m s-1 at 2 m over grass, W m-2
canopy_t = np.array([31.0, 33.5, 36.0])          # e.g. three trees in a thermal image [C]

ra = aerodynamic_resistance_viney(wind, air_t, canopy_t, canopy_height=4.0, measurement_height=6.0)
rn = available_energy(shortwave, 0.15, longwave_downwelling(air_t, rh), canopy_t, emissivity=0.98)
gc = canopy_conductance(canopy_t, air_t, rh, pressure, rn, ra)
ra    # s m-1
# array([12.77, 10.31,  8.54])
gc    # mm s-1
# array([16.72,  7.54,  0.78])
```

`canopy_conductance` returns NaN where no physical conductance exists (latent heat flux Rn − H ≤ 0, or a negative
canopy resistance).

The aerodynamic resistance functions take a `von_karman` argument (default 0.41, as in FAO-56). At low conductance
the result is sensitive to it: through the sensible heat flux, the 5 % difference in ra between k = 0.40 and 0.41
changes gc by 10–30 % at 1–3 mm s⁻¹. Use the value of the study you want to reproduce. By default the wind is assumed
to be measured at 2 m over reference grass and is transferred to the crop through a 100 m blending height; pass
`wind_over_grass=False` if it was measured at `measurement_height` over the crop. A `ValueError` is raised if the
reference height is inside or too close to the canopy.

### Simulating canopy temperature

`canopy_temperature_sim` solves the canopy energy balance exactly for a given canopy conductance, recomputing the
net radiation with the canopy temperature at each step. The aerodynamic resistance can be a function of canopy
temperature, so that atmospheric stability is included. It is the exact inverse of `canopy_conductance`:

```python
from opencroplib.plant import canopy_temperature_sim

def ra_of_tc(tc):
    return aerodynamic_resistance_viney(wind, air_t, tc, canopy_height=4.0, measurement_height=6.0)

gc_true = np.array([0.0, 2.0, 5.0, 10.0])        # mm s-1; 0 = non-transpiring canopy (upper limit)
tc = canopy_temperature_sim(shortwave, air_t, rh, pressure, gc_true, ra_of_tc, albedo=0.15)
tc
# array([36.34, 35.49, 34.36, 32.76])

rn = available_energy(shortwave, 0.15, longwave_downwelling(air_t, rh), tc, emissivity=0.98)
canopy_conductance(tc, air_t, rh, pressure, rn, ra_of_tc(tc))
# array([nan,  2.,  5., 10.])    (gc = 0 has no latent heat flux to invert)
```

The non-transpiring limit (gc = 0) and a well-watered reference give the bounds used, for example, in crop water stress
indices.

### Solar position and clear-sky irradiance

```python
import numpy as np
import pandas as pd
from opencroplib.time_utils import get_day_of_year, get_solar_time
from opencroplib.radiation import (solar_elevation, solar_azimuth, potential_irradiance, clear_sky_beam,
                                   clear_sky_diffuse)

latitude, longitude, altitude = 37.857, -4.803, 117.0
times = pd.date_range("2024-06-21 06:00", "2024-06-21 21:00", freq="3h", tz="Europe/Madrid")
doy = get_day_of_year(times)
solar_h = get_solar_time(longitude, times)        # local solar time [h]

sun = pd.DataFrame({
    "solar_h": solar_h,
    "elevation": np.degrees(solar_elevation(latitude, doy, solar_h)),
    "azimuth": np.degrees(solar_azimuth(latitude, doy, solar_h)),   # from south, west positive
    "global": potential_irradiance(latitude, doy, solar_h, altitude, t_lk=3.0),
    "beam": clear_sky_beam(latitude, doy, solar_h, altitude, t_lk=3.0),
    "diffuse": clear_sky_diffuse(latitude, doy, solar_h, t_lk=3.0),
}, index=times)
sun.round(1)
#                            solar_h  elevation  azimuth  global   beam  diffuse
# 2024-06-21 06:00:00+02:00      3.7      -10.0   -130.5     0.0    0.0      0.0
# 2024-06-21 09:00:00+02:00      6.7       21.5   -103.6   324.4  252.0     72.5
# 2024-06-21 12:00:00+02:00      9.7       56.7    -74.6   877.1  772.1    105.0
# 2024-06-21 15:00:00+02:00     12.7       73.4     32.9  1023.5  918.4    105.1
# 2024-06-21 18:00:00+02:00     15.7       41.5     89.2   668.5  569.6     98.9
# 2024-06-21 21:00:00+02:00     18.7        7.0    114.4    79.9   44.9     35.0
```

`t_lk` is the Linke turbidity factor: about 2–3 for clear rural skies and 4–6 for turbid or polluted air. The default
of 1 corresponds to a clean, dry atmosphere. The azimuth follows Duffie and Beckman (from south, east negative, west
positive); add 180° for the azimuth clockwise from north.

### Thermal sensors

```python
from opencroplib.radiation import planck, iplanck

radiance = planck(30.0 + 273.15)        # spectral radiance at 10.25 um [W m-2 sr-1 um-1]
round(radiance, 3)
# 10.367
round(iplanck(radiance) - 273.15, 2)    # brightness temperature [C]
# 30.0
```

More examples are in the docstrings of each function. The `sample_notebooks` folder contains a notebook that computes
daily ET₀ for the IAS-CSIC station; it was written for older versions of pandas and IPython.

## Validation

The test suite (`uv run pytest`) compares the library with published examples and independent implementations:

| Component | Reference | Result |
|---|---|---|
| FAO-56 and ASCE ET₀ equations | [refet](https://github.com/WSWUP/RefET) (ASCE-EWRI 2005), fed the same Δ, γ and VPD | Identical to machine precision (hourly and daily) |
| FAO-56 and ASCE ET₀, complete | refet and [pyet](https://github.com/pyet-org/pyet) with their own thermodynamics | Mean difference 0.5 % (daytime max 2 %), from the psychrometric constant (Bringfelt 1986 vs 0.000665 P) |
| Hourly and daily ET₀ | FAO-56 Examples 18 and 19 | Hourly 0.627 mm h⁻¹ (book 0.63) |
| Wind height conversion | FAO-56 Eq. 47 (refet) | Within 0.03 % |
| ESRA clear-sky beam and diffuse irradiance | GRASS GIS r.sun (`rsunlib.c`) | Within 10⁻¹² W m⁻² |
| Clear-sky global irradiance | Ineichen and Perez (2002) in pvlib | Within ±10 % (different model; typically +5 %) |
| Solar elevation and azimuth | NREL SPA (Reda and Andreas 2004) in pvlib, 4 years hourly, 5 latitudes | Elevation 0.22° RMS (max 0.6°); azimuth 0.3° RMS for elevations below 70° |
| Solar declination | Spencer (1971) in pvlib | Identical |
| Equation of time | Spencer (1971) in pvlib | Within 1.2 min |
| Day of year and solar time | `datetime` and `pandas.DatetimeIndex` paths, 5 time zones | Identical |
| Planck functions | Planck's law with `scipy.constants` | Within 10⁻¹² (relative); exact round trip |
| Extraterrestrial radiation, net longwave | FAO-56 Examples 8, 10 and 19 | 32.19 MJ m⁻² day⁻¹ (book 32.2), 3.53 (3.5), 0.137 MJ m⁻² h⁻¹ (0.138) |
| Canopy conductance and canopy temperature | Round trip: simulate Tc for a known gc and invert | gc recovered within 0.1 % |

## Upgrading from 0.1.x

Version 0.2.0 fixes several bugs that change numerical results. Results obtained with the following functions in
0.1.x should be recomputed:

| Function | Change in 0.2.0 |
|---|---|
| `plant.canopy_conductance` | The inversion used the VPD where the actual vapour pressure is needed. The error is zero at 50 % RH and grows with the distance from it: in dry Mediterranean summers gc was overestimated about 3× (median) and was sometimes negative or infinite. |
| `plant.canopy_temperature_sim` | Now an exact solution of the energy balance (it was a 5-iteration linearisation); gc = 0 gives the non-transpiring limit. |
| `atmophere.rho_cp_dry` | Used temperature instead of pressure for the air density (about 5000× too small). |
| `evapotranspiration.et0_FAO56`, `et0_ASCE` (hourly) | Wind coefficient 9.4 instead of 37, and the default soil heat flux was ~0. Annual totals change by +30 to +46 %. |
| `et0_*` with `z_wind != 2` | The wind height conversion went in the wrong direction. |
| `radiation.potential_irradiance` | Miscoded ESRA beam model; now returns global (beam + diffuse) irradiance. |
| `radiation.net_out_lw_hourly` | Vapour pressure in Pa instead of kPa and wrong Stefan–Boltzmann constant; the old values were unusable. |
| `radiation.solar_declination`, `solar_elevation`, `solar_azimuth` | Declination amplitude was 24.1° instead of 23.44°; the azimuth is now signed. |
| `aerodynamic.aerodynamic_resistance_viney` | von Kármán constant 0.41 instead of 0.40 (−4.8 % in ra; use `von_karman=0.40` to reproduce the old value), corrected wind transfer, NaN for missing temperatures. |
| `time_utils.get_day_of_year` | For a single `datetime` it returned the days *after* 1 January (0 for 1 January), one less than for a `DatetimeIndex` and than the radiation functions expect. Now 1-based in both cases. |
| `time_utils.get_solar_time` | For a single `datetime` it could return negative hours or hours above 24; now always in [0, 24), as for a `DatetimeIndex`. |
| `radiation.planck`, `iplanck` | Both use the exact SI (2019) constants; they used different CODATA sets, so the round trip was off by ~1.4 mK. |

Other changes: input arrays are no longer modified in place, list inputs are accepted, aerodynamic functions raise
`ValueError` for reference heights inside the canopy, and `canopy_conductance` returns NaN instead of non-physical
values. The commit messages give the details and magnitudes for each change.

## Development

The project is managed with [uv](https://docs.astral.sh/uv/):

```
uv sync          # create .venv with the library and the dev dependencies
uv run pytest    # run the test suite
uv build         # build the sdist and wheel
```

The reference implementations used by the tests ([refet](https://github.com/WSWUP/RefET),
[pyet](https://github.com/pyet-org/pyet) and [pvlib](https://github.com/pvlib/pvlib-python)) are installed as dev
dependencies.

## Citation

This library is licensed under the MIT license. If you use it, please cite the software:

Jose A. Jimenez-Berni. (2023). OpenAgriTech/opencroplib: v0.1.0 (v0.1.0). Zenodo.
https://doi.org/10.5281/zenodo.8364969

and, when using the canopy conductance model:

Berni, J.A.J., Zarco-Tejada, P.J., Sepulcre-Cantó, G., Fereres, E., Villalobos, F. (2009). Mapping canopy conductance
and CWSI in olive orchards using high resolution thermal remote sensing imagery. *Remote Sensing of Environment* 113,
2380–2388. https://doi.org/10.1016/j.rse.2009.06.018

## References

### Energy balance and canopy conductance

- Berni, J.A.J., Zarco-Tejada, P.J., Sepulcre-Cantó, G., Fereres, E., Villalobos, F. (2009). Mapping canopy conductance
  and CWSI in olive orchards using high resolution thermal remote sensing imagery. *Remote Sensing of Environment* 113,
  2380–2388. https://doi.org/10.1016/j.rse.2009.06.018

### Reference evapotranspiration

- Allen, R.G., Pereira, L.S., Raes, D., Smith, M. (1998). *Crop evapotranspiration: Guidelines for computing crop
  water requirements*. FAO Irrigation and Drainage Paper 56. FAO, Rome. https://www.fao.org/4/x0490e/x0490e00.htm
- ASCE-EWRI (2005). *The ASCE Standardized Reference Evapotranspiration Equation*. Allen, R.G., Walter, I.A.,
  Elliott, R.L., Howell, T.A., Itenfisu, D., Jensen, M.E., Snyder, R.L. (eds.). American Society of Civil Engineers,
  Reston, VA.

### Aerodynamic resistance

- Raupach, M.R. (1994). Simplified expressions for vegetation roughness length and zero-plane displacement as
  functions of canopy height and area index. *Boundary-Layer Meteorology* 71, 211–216.
- Viney, N.R. (1991). An empirical expression for aerodynamic resistance in the unstable boundary layer.
  *Boundary-Layer Meteorology* 56, 381–393.

### Vapour pressure and psychrometrics

- Bringfelt, B. (1986). *Test of a forest evapotranspiration model*. Meteorology and Climatology Reports 52, SMHI,
  Norrköping, Sweden.
- Buck, A.L. (1981). New equations for computing vapor pressure and enhancement factor. *Journal of Applied
  Meteorology* 20, 1527–1532. Coefficients as revised in Buck Research (1996), *CR-1A User's Manual*, Appendix 1.
- Goff, J.A., Gratch, S. (1946). Low-pressure properties of water from −160 to 212 °F. *Transactions of the American
  Society of Heating and Ventilating Engineers* 52, 95–122.
- Goff, J.A. (1957). Saturation pressure of water on the new Kelvin temperature scale. *Transactions of the American
  Society of Heating and Ventilating Engineers* 63, 347–354.
- Waterloo, M.J. (2014). *Meteorology and evaporation modules documentation* (meteolib).
  http://python.hydrology-amsterdam.nl/moduledoc/Meteorologyandevaporationmodules.pdf
- World Meteorological Organization (1984). *Technical Regulations*, WMO No. 49, Appendix A.

### Longwave radiation

- Alados, I., Foyo-Moreno, I., Alados-Arboledas, L. (2012). Estimation of downwelling longwave irradiance under
  all-sky conditions. *International Journal of Climatology* 32, 781–793.
- Berdahl, P., Martin, M. (1984). Emissivity of clear skies. *Solar Energy* 32, 663–664.
- Brunt, D. (1932). Notes on radiation in the atmosphere. *Quarterly Journal of the Royal Meteorological Society* 58,
  389–420.
- Brutsaert, W. (1975). On a derivable formula for long-wave radiation from clear skies. *Water Resources Research* 11,
  742–744.
- Crawford, T.M., Duchon, C.E. (1999). An improved parameterization for estimating effective atmospheric emissivity
  for use in calculating daytime downwelling longwave radiation. *Journal of Applied Meteorology* 38, 474–480.
- Idso, S.B. (1981). A set of equations for full spectrum and 8- to 14-µm and 10.5- to 12.5-µm thermal radiation from
  cloudless skies. *Water Resources Research* 17, 295–304.
- Swinbank, W.C. (1963). Long-wave radiation from clear skies. *Quarterly Journal of the Royal Meteorological Society*
  89, 339–348.

### Solar geometry and clear-sky irradiance

- Duffie, J.A., Beckman, W.A. (2013). *Solar Engineering of Thermal Processes*, 4th ed. Wiley.
- Hofierka, J., Šúri, M. (2002). The solar radiation model for Open source GIS: implementation and applications.
  *Proceedings of the Open Source GIS – GRASS Users Conference*, Trento, Italy.
- Iqbal, M. (1983). *An Introduction to Solar Radiation*. Academic Press.
- Kasten, F., Young, A.T. (1989). Revised optical air mass tables and approximation formula. *Applied Optics* 28,
  4735–4738.
- Rigollier, C., Bauer, O., Wald, L. (2000). On the clear sky model of the ESRA — European Solar Radiation Atlas —
  with respect to the Heliosat method. *Solar Energy* 68, 33–48.
- Spencer, J.W. (1971). Fourier series representation of the position of the sun. *Search* 2(5), 172.

### Reference implementations used for validation

- Anderson, K.S., Hansen, C.W., Holmgren, W.F., Jensen, A.R., Mikofski, M.A., Driesse, A. (2023). pvlib python: 2023
  project update. *Journal of Open Source Software* 8(92), 5994. https://doi.org/10.21105/joss.05994
- GRASS GIS r.sun module. https://github.com/OSGeo/grass/tree/main/raster/r.sun
- Ineichen, P., Perez, R. (2002). A new airmass independent formulation for the Linke turbidity coefficient.
  *Solar Energy* 73, 151–157.
- Reda, I., Andreas, A. (2004). Solar position algorithm for solar radiation applications. *Solar Energy* 76,
  577–589.
- RefET: ASCE-EWRI standardized reference ET, Western States Water Use Program. https://github.com/WSWUP/RefET
- Vremec, M., Collenteur, R.A., Birk, S. (2024). PyEt v1.3.1: a Python package for the estimation of potential
  evapotranspiration. *Geoscientific Model Development* 17, 7083–7103. https://doi.org/10.5194/gmd-17-7083-2024

### Software this library builds on

- PyETo, by Mark Richards. https://github.com/woodcrafty/PyETo
- meteolib, by Maarten J. Waterloo (see Waterloo 2014 above).
