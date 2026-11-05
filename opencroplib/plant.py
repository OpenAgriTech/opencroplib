# -*- coding: utf-8 -*-
#
# The MIT License (MIT)
#
#    Copyright (c) 2017-present Jose A. Jimenez-Berni
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
#
import numpy as np
from .atmophere import psychrometric_gamma, es_calc, ea_calc, rho_cp_dry, rho_cp_moist
from .radiation import longwave_downwelling, available_energy

__author__ = "Jose A. Jimenez-Berni  <berni@ias.csic.es>"
__version__ = "0.1.3"
__date__ = "October 2018"

"""


"""

# Canopy temperature simulation: the root is searched in [Ta - 30, Ta + 30] Celsius
_TC_BRACKET = 30.0
_TC_TOLERANCE = 1e-8
_TC_MAX_ITERATIONS = 100


def _restore_type(values, *templates):
    """
    Return numpy scalars for 0-d results, and a pandas Series (with the index of the first matching 1-d template)
    when one of the inputs was a Series of the same shape. Otherwise return the ndarray.
    """
    values = np.asarray(values)
    if values.ndim == 0:
        return values[()]
    for template in templates:
        if getattr(template, "ndim", None) == 1 and hasattr(template, "index") and np.shape(template) == values.shape:
            return type(template)(values, index=template.index)
    return values


def canopy_conductance(canopy_temperature, air_temperature, relative_humidity, atmospheric_pressure, available_energy,
                       aerodynamic_resistance):
    """
    Function to calculate the canopy conductance from micrometeorological parameters, by inverting the
    energy balance of the canopy (Berni et al. 2009, Eq. 4):

        rc = ra * (e*(Tc) - ea) / (gamma * (ra * Rn / (rho * cp) - (Tc - Ta))) - ra

    which is equivalent to rc = rho * cp * (e*(Tc) - ea) / (gamma * LE) - ra, with LE = Rn - H and
    H = rho * cp * (Tc - Ta) / ra.

    The result is NaN where the latent heat flux LE = Rn - H <= 0 or where rc <= 0, as no physical (finite, positive)
    conductance exists in those cases.

    Reference: Berni, J.A.J., Zarco-Tejada, P.J., Sepulcre-Cantó, G., Fereres, E., Villalobos, F., 2009. Mapping
    canopy conductance and CWSI in olive orchards using high resolution thermal remote sensing imagery. Remote Sensing
    of Environment 113, 2380–2388. https://doi.org/10.1016/j.rse.2009.06.018

    Parameters
    ----------
    canopy_temperature : Union[ndarray, float]
        Temperature of the canopy [Celsius]
    air_temperature : Union[ndarray, float]
        Temperature of the air [Celsius]
    relative_humidity : Union[ndarray, float]
        Relative humidity [%]
    atmospheric_pressure : Union[ndarray, float]
        Atmospheric pressure [Pa]
    available_energy : Union[ndarray, float]
        Net radiation or available energy at the canopy [W/m2]
    aerodynamic_resistance : Union[ndarray, float]
        Aerodynamic resistance of the canopy [s/m]

    Returns
    -------
    Union[ndarray, float]: canopy conductance [mm/s]

    """
    gam = psychrometric_gamma(air_temperature, relative_humidity, atmospheric_pressure)
    rho_cp = rho_cp_moist(air_temperature, relative_humidity, atmospheric_pressure)

    with np.errstate(divide="ignore", invalid="ignore"):
        latent_heat = available_energy - rho_cp * (canopy_temperature - air_temperature) / aerodynamic_resistance
        canopy_resistance = rho_cp * (es_calc(canopy_temperature) - ea_calc(air_temperature, relative_humidity)) / (
                gam * latent_heat) - aerodynamic_resistance
        conductance = 1000.0 / canopy_resistance
        valid = (latent_heat > 0) & (canopy_resistance > 0)

    return _restore_type(np.where(valid, conductance, np.nan), conductance)


def canopy_conductance_2(canopy_temperature, air_temperature, relative_humidity, atmospheric_pressure, available_energy,
                         aerodynamic_resistance):
    """
    Deprecated alias of canopy_conductance, kept for backwards compatibility. See canopy_conductance.

    Returns
    -------
    Union[array, float]: canopy conductance [mm/s]

    """
    return canopy_conductance(canopy_temperature, air_temperature, relative_humidity, atmospheric_pressure,
                              available_energy, aerodynamic_resistance)


def sensible_heat_flux_2(air_temperature, relative_humidity, atmospheric_pressure, canopy_temperature,
                         aerodynamic_resistance):
    """
    Sensible heat flux using the volumetric heat capacity of dry air (see sensible_heat_flux for moist air).
    The aerodynamic resistance is clamped to a minimum of 0.01 s/m (on a copy, the input is not modified).

    Returns
    -------
    Union[ndarray, float]: sensible heat flux [W/m2]
    """
    aerodynamic_resistance = np.maximum(aerodynamic_resistance, 0.01)
    tc_ta = canopy_temperature - air_temperature
    return rho_cp_dry(air_temperature, atmospheric_pressure) * tc_ta / aerodynamic_resistance


def sensible_heat_flux(air_temperature, relative_humidity, atmospheric_pressure, canopy_temperature,
                       aerodynamic_resistance):
    pCpp = rho_cp_moist(air_temperature, relative_humidity, atmospheric_pressure)
    tc_ta = canopy_temperature - air_temperature
    return pCpp * tc_ta / aerodynamic_resistance


def latent_heat_flux(air_temperature, relative_humidity, atmospheric_pressure, canopy_temperature,
                     aerodynamic_resistance, canopy_resistance):
    pCpp = rho_cp_moist(air_temperature, relative_humidity, atmospheric_pressure)
    es_canopy = es_calc(canopy_temperature)
    e_air = ea_calc(air_temperature, relative_humidity)
    gam = psychrometric_gamma(air_temperature, relative_humidity, atmospheric_pressure)
    return pCpp * (es_canopy - e_air) / (gam * (aerodynamic_resistance + canopy_resistance))


def canopy_temperature_sim(shortwave, air_temperature, relative_humidity, atmospheric_pressure, gc,
                           aerodynamic_resistance, albedo=0.23, longwave=None, emissivity=0.98):
    """
    Simulates the canopy temperature for a given canopy conductance and micro-meteorological parameters

    Solves exactly (vectorised bisection) the canopy energy balance

        Rn(Tc) - H(Tc) - LE(Tc) = 0

    with Rn from available_energy (recomputed with Tc on every step), H from sensible_heat_flux and LE from
    latent_heat_flux with rc = 1000 / gc. It is therefore the exact inverse of canopy_conductance.
    The root is searched in [Ta - 30, Ta + 30]; NaN is returned where there is no root in that interval.

    Berni et al. (2009)


    Parameters
    ----------
    shortwave: Union[array, float]
        Downwelling shortwave radiation [W/m2]
    albedo: Union[array, float]
        Canopy albedo, default 0.23 [0.0-1.0]
    air_temperature: Union[array, float]
        Air temperature [Celsius]
    relative_humidity: Union[array, float]
        Relative humidity [%]
    atmospheric_pressure: Union[array, float]
        Atmospheric pressure [Pa]
    gc: Union[array, float]
        Canopy conductance [mm/s]. gc = 0 gives the non-transpiring (upper) limit. Negative values give NaN.
    aerodynamic_resistance: Union[array, float, callable]
        Aerodynamic resistance [s/m], or a function ra(canopy_temperature) returning it, so that stability
        corrections can be included, e.g.
        ``lambda tc: aerodynamic_resistance_viney(u, ta, tc, h, z)``. The function receives an ndarray of canopy
        temperatures [Celsius] with the broadcast shape of the inputs.
    longwave: Union[array, float]
        Downwelling longwave radiation [W/m2]. If not provider will be calculated using air temperature and humidity
        assuming clear skies.
    emissivity: Union[array, float]
        Canopy emissivity. Default 0.98

    Returns
    -------
    Union[array, float]: simulated canopy temperature [Celsius]

    """
    # If longwave radiation is not provided we calculate it from using the model and defaults for clear skies
    if longwave is None:
        longwave = longwave_downwelling(air_temperature, relative_humidity)

    if callable(aerodynamic_resistance):
        ra_function = aerodynamic_resistance
    else:
        def ra_function(_):
            return aerodynamic_resistance

    gc_array = np.asarray(gc, dtype=float)
    with np.errstate(divide="ignore"):
        canopy_resistance = np.where(gc_array > 0, 1000.0 / gc_array, np.inf)
    canopy_resistance = np.where(gc_array < 0, np.nan, canopy_resistance)

    def residual(tc):
        ra = ra_function(tc)
        rn = available_energy(shortwave, albedo, longwave, tc, emissivity)
        h = sensible_heat_flux(air_temperature, relative_humidity, atmospheric_pressure, tc, ra)
        le = latent_heat_flux(air_temperature, relative_humidity, atmospheric_pressure, tc, ra, canopy_resistance)
        return np.asarray(rn - h - le, dtype=float)

    ta = np.asarray(air_temperature, dtype=float)
    f_lo = residual(ta - _TC_BRACKET)
    f_hi = residual(ta + _TC_BRACKET)
    lo = np.broadcast_to(ta - _TC_BRACKET, f_lo.shape).astype(float)
    hi = np.broadcast_to(ta + _TC_BRACKET, f_lo.shape).astype(float)
    bracketed = np.sign(f_lo) * np.sign(f_hi) <= 0

    for _ in range(_TC_MAX_ITERATIONS):
        mid = 0.5 * (lo + hi)
        f_mid = residual(mid)
        move_lo = np.sign(f_mid) == np.sign(f_lo)
        lo = np.where(move_lo, mid, lo)
        f_lo = np.where(move_lo, f_mid, f_lo)
        hi = np.where(move_lo, hi, mid)
        if not np.any(hi - lo > _TC_TOLERANCE):
            break

    canopy_temperature = np.where(bracketed, 0.5 * (lo + hi), np.nan)
    return _restore_type(canopy_temperature, air_temperature, shortwave, gc, relative_humidity)
