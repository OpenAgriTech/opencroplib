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
from .atmophere import delta_calc, psychrometric_gamma, vpd_calc
from .aerodynamic import adjust_wind_height

__author__ = "Jose A. Jimenez-Berni  <berni@ias.csic.es>"
__version__ = "0.0.1"
__date__ = "July 2019"


def et0_FAO56(air_temperature, relative_humidity, air_pressure, wind_speed, solar_radiation, net_radiation, z_wind=2.0,
              soil_heat_flux=None):
    """
    Estimate ETo using FAO56.
    Calculates ETo using FAO56 for the time steps provided in the input dataframes.
    Daily values can be calculated using time aggregation using
    ``eto_FAO56_daily()``.
    Based on FAO equation 39 in Allen et al (1998).
    :param air_temperature: Air temperature [degrees Celsius]
    :param relative_humidity: Air relative humidity [%]
    :param air_pressure: Barometric pressure [Pa]
    :param solar_radiation: Solar radiation [W m-2]
    :param net_radiation: Net radiation [W m-2]
    :param wind_speed: Wind speed [m s-1]
    :param z_wind: Height of the wind speed measurements [m]. Default 2m
    :param soil_heat_flux: If G is measured, it can be provided here. If None, it will be calculated as
    G=0.1Rn during the day and G=0.5Rn at night time
    :return: ETo [mm h-1]
    :rtype: float
    """

    # Correct wind speed if measurement is not 2m
    if z_wind != 2.0:
        wind_speed = adjust_wind_height(wind_speed, 2.0, z_wind)

    delta = delta_calc(air_temperature)
    gamma = psychrometric_gamma(air_temperature, relative_humidity, air_pressure)
    vpd = vpd_calc(air_temperature, relative_humidity)
    net_radiation_flux_density = net_radiation*0.0036
    if soil_heat_flux is None:
        if np.isscalar(solar_radiation):
            soil_heat_flux = 0.0
        else:
            soil_heat_flux = np.where(solar_radiation > 0, net_radiation_flux_density*0.1,
                                      net_radiation_flux_density*0.5)

    soil_heat_flux_density = soil_heat_flux*0.0036

    et0 = (0.408*delta*(net_radiation_flux_density-soil_heat_flux_density)+gamma*9.4/(air_temperature+273.0) *
           wind_speed*vpd/1000.0)/(delta+gamma*(1+0.34*wind_speed))

    return et0


def et0_FAO56_daily(air_temperature, relative_humidity, air_pressure, wind_speed, solar_radiation, net_radiation,
                    z_wind=2.0, soil_heat_flux=None):
    """
    Estimate ETo using FAO56.
    Calculates ETo using FAO56 for the time steps provided in the input dataframes.
    Daily values can be calculated using time aggregation using
    ``eto_FAO56_daily()``.
    Based on FAO equation 39 in Allen et al (1998).
    :param air_temperature: Air temperature [degrees Celsius]
    :param relative_humidity: Air relative humidity [%]
    :param air_pressure: Barometric pressure [Pa]
    :param solar_radiation: Solar radiation [MJ day-2]
    :param net_radiation: Net radiation [MJ day-2]
    :param wind_speed: Wind speed [m s-1]
    :param z_wind: Height of the wind speed measurements [m]. Default 2m
    :param soil_heat_flux: If G is measured, it can be provided here. If None, it will be calculated as
    G=0.1Rn during the day and G=0.5Rn at night time
    :return: ETo [mm h-1]
    :rtype: float
    """

    # Correct wind speed if measurement is not 2m
    if z_wind != 2.0:
        wind_speed = adjust_wind_height(wind_speed, 2.0, z_wind)

    delta = delta_calc(air_temperature)
    gamma = psychrometric_gamma(air_temperature, relative_humidity, air_pressure)
    vpd = vpd_calc(air_temperature, relative_humidity)
    if soil_heat_flux is None:
        # Assume null soil heat flux if not provided
        if np.isscalar(solar_radiation):
            soil_heat_flux = 0.0
        else:
            soil_heat_flux = 0.0

    et0 = (0.408*delta*(net_radiation-soil_heat_flux)+gamma*900.0/(air_temperature+273.0) *
           wind_speed*vpd/1000.0)/(delta+gamma*(1+0.34*wind_speed))

    return et0


def et0_ASCE(air_temperature, relative_humidity, air_pressure, wind_speed, solar_radiation, net_radiation, z_wind=2.0,
             soil_heat_flux=None):
    """
    Estimate ETo using FAO56.
    Calculates ETo using FAO56 for the time steps provided in the input dataframes.
    Daily values can be calculated using time aggregation using
    ``eto_FAO56_daily()``.
    Based on FAO equation 39 in Allen et al (1998).
    :param air_temperature: Air temperature [degrees Celsius]
    :param relative_humidity: Air relative humidity [%]
    :param air_pressure: Barometric pressure [Pa]
    :param solar_radiation: Solar radiation [W m-2]
    :param net_radiation: Net radiation [W m-2]
    :param wind_speed: Wind speed [m s-1]
    :param z_wind: Height of the wind speed measurements [m]. Default 2m
    :param soil_heat_flux: If G is measured, it can be provided here. If None, it will be calculated as
    G=0.1Rn during the day and G=0.5Rn at night time
    :return: ETo [mm h-1]
    :rtype: float
    """

    # Correct wind speed if measurement is not 2m
    if z_wind != 2.0:
        wind_speed = adjust_wind_height(wind_speed, 2.0, z_wind)

    delta = delta_calc(air_temperature)
    gamma = psychrometric_gamma(air_temperature, relative_humidity, air_pressure)
    vpd = vpd_calc(air_temperature, relative_humidity)
    net_radiation_flux_density = net_radiation*0.0036
    if soil_heat_flux is None:
        if np.isscalar(solar_radiation):
            soil_heat_flux = 0.0
        else:
            soil_heat_flux = np.where(solar_radiation > 0, net_radiation_flux_density*0.1,
                                      net_radiation_flux_density*0.5)

    soil_heat_flux_density = soil_heat_flux*0.0036

    if np.isscalar(net_radiation_flux_density):
        cd = 0.24 if net_radiation_flux_density >= 0 else 0.96
    else:
        cd = np.where(net_radiation_flux_density >= 0, 0.24, 0.96)

    et0 = (0.408*delta*(net_radiation_flux_density-soil_heat_flux_density)+gamma*9.4/(air_temperature+273.0) *
           wind_speed*vpd/1000.0)/(delta+gamma*(1+cd*wind_speed))

    return et0
