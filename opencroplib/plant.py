# -*- coding: utf-8 -*-
#
#    Copyright 2017 Jose A. Jimenez-Berni
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.
#
import numpy as np
from .atmophere import vpd_calc, psychrometric_gamma, specific_heat_cp, es_calc, ea_calc, delta_calc, rho_cp_dry, \
    rho_cp_moist
from .radiation import longwave_downwelling, available_energy

__author__ = "Jose A. Jimenez-Berni  <berni@ias.csic.es>"
__version__ = "0.1.3"
__date__ = "October 2018"

"""


"""


def canopy_conductance(canopy_temperature, air_temperature, relative_humidity, atmospheric_pressure, available_energy,
                       aerodynamic_resistance):
    """
    Function to calculate the canopy conductance from micrometeorological parameters

    Berni et al. (2009)

    Input:
        - canopy_temperature: temperature of the canopy [Celsius]
        - air_temperature: temperature of the air [Celsius]
        - relative_humidity: relative humidity [%]
        - atmospheric_pressure: atmospheric pressure [Pa]
        - available_energy: net radiation or available energy at the canopy [W/m2]
        - aerodynamic_resistance: aerodynamic resistance of the canopy [s/m]

    Output:
        - canopy_conductance: canopy conductance [mm/s]

    Parameters
    ----------
    canopy_temperature
    air_temperature
    relative_humidity
    atmospheric_pressure
    available_energy
    aerodynamic_resistance

    Returns
    -------
    Union[array, float]: canopy conducntace [mm/s]

    """
    gam = psychrometric_gamma(air_temperature, relative_humidity, atmospheric_pressure)
    pCpp = rho_cp_moist(air_temperature, relative_humidity, atmospheric_pressure)
    tc_ta = canopy_temperature - air_temperature

    conductance = 1000.0 / (
        aerodynamic_resistance * (es_calc(canopy_temperature) - (es_calc(air_temperature) - ea_calc(air_temperature,
                                                                                                    relative_humidity))) / (
            gam * (aerodynamic_resistance * available_energy / pCpp - tc_ta)) - aerodynamic_resistance)
    return conductance


def canopy_conductance_2(canopy_temperature, air_temperature, relative_humidity, atmospheric_pressure, available_energy,
                         aerodynamic_resistance):
    """
    Function to calculate the canopy conductance from micrometeorological parameters

    Berni et al. (2009)

    Input:
        - canopy_temperature: temperature of the canopy [Celsius]
        - air_temperature: temperature of the air [Celsius]
        - relative_humidity: relative humidity [%]
        - atmospheric_pressure: atmospheric pressure [Pa]
        - available_energy: net radiation or available energy at the canopy [W/m2]
        - aerodynamic_resistance: aerodynamic resistance of the canopy [s/m]

    Output:
        - canopy_conductance: canopy conductance [mm/s]

    Parameters
    ----------
    canopy_temperature
    air_temperature
    relative_humidity
    atmospheric_pressure
    available_energy
    aerodynamic_resistance

    Returns
    -------
    Union[array, float]: canopy conducntace [mm/s]

    """
    rho = atmospheric_pressure / (287.07 * (air_temperature + 273.15))
    rho_cp = (rho * (1002.5 + 275e-6 * np.power(air_temperature + 273.15 - 200.0, 2)))
    le = (available_energy - sensible_heat_flux_2(air_temperature, relative_humidity,
                                                  atmospheric_pressure,
                                                  canopy_temperature, aerodynamic_resistance))
    gam = psychrometric_gamma(air_temperature, relative_humidity, atmospheric_pressure)
    gam_le = (gam * le)
    #gam_le[gam_le < 0.01] = 0.01

    rc = rho_cp * (
        es_calc(canopy_temperature) - ea_calc(air_temperature, relative_humidity)) / gam_le - aerodynamic_resistance

    return 1000.0 / rc


def sensible_heat_flux_2(air_temperature, relative_humidity, atmospheric_pressure, canopy_temperature,
                         aerodynamic_resistance):
    rho = atmospheric_pressure / (287.07 * (air_temperature + 273.15))
    rho_cp = (rho * (1002.5 + 275e-6 * np.power(air_temperature + 273.15 - 200.0, 2)))

    aerodynamic_resistance[aerodynamic_resistance < 0.01] = 0.01
    tc_ta = canopy_temperature - air_temperature
    return rho_cp * tc_ta / aerodynamic_resistance


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
        Canopy conductance [mm/s]
    aerodynamic_resistance: Union[array, float]
        Aerodynamic resistance [s/m]
    longwave: Union[array, float]
        Downwelling longwave radiation [W/m2]. If not provider will be calculated using air temperature and humidity
        assuming clear skies.
    emissivity: Union[array, float]
        Canopy emissivity. Default 0.98

    Returns
    -------
    Union[array, float]: simulated canopy temperature [Celsius]

    """
    # Start with arbitrary value
    canopy_temperature = air_temperature + 3.0

    gam = psychrometric_gamma(air_temperature, relative_humidity, atmospheric_pressure)
    delt = delta_calc(air_temperature)
    pCpp = rho_cp_moist(air_temperature, relative_humidity, atmospheric_pressure)
    vpd = vpd_calc(air_temperature, relative_humidity)
    canopy_resistance = 1000.0 / gc

    # If longwave radiation is not provided we calculate it from using the model and defaults for clear skies
    if longwave is None:
        longwave = longwave_downwelling(air_temperature, relative_humidity)

    # Iterate for implicit canopy temperature
    for i in range(5):
        rn = available_energy(shortwave, albedo, longwave, canopy_temperature, emissivity)
        tc_ta = (aerodynamic_resistance * rn * gam * (1 + canopy_resistance / aerodynamic_resistance) / pCpp - vpd) / (
            delt + gam * (1 + canopy_resistance / aerodynamic_resistance))
        canopy_temperature = tc_ta + air_temperature

    return canopy_temperature
