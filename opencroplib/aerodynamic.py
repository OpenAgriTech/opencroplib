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

import numpy

__author__ = "Jose A. Jimenez-Berni  <berni@ias.csic.es>"
__version__ = "0.1"
__date__ = "July 2017"


def aerodynamic_resistance_viney(wind_speed, air_temperature, canopy_temperature, canopy_height, measurement_height,
                                 d_ratio=2/3.0, z0m_ratio=0.13, z0h_ratio=0.13):
    """
    Calculates aerodynamic resistance based on Viney 1991

    Reference: Viney, Neil R. 1991. “An Empirical Expression for Aerodynamic Resistance in the Unstable Boundary Layer.”
        Boundary-Layer Meteorology 56 (4). Kluwer Academic Publishers: 381–93.


    Parameters
    ----------
    d_ratio : float
        Zero-plane displacement (d) = d_ratio * canopy_height
    z0m_ratio : float
        Roughness length (z0) for momentum = z0m_ratio * canopy_height
    z0h_ratio : float
        Roughness length (z0) for heat transfer = z0h_ratio * canopy_height
    wind_speed : Union[ndarray, float]
        Wind speed measured at a given height (measurement_height) [m/s]
    air_temperature : Union[ndarray, float]
        Air temperature [Celsius]
    canopy_temperature : Union[ndarray, float]
        Canopy temperature, used to estimate buoyancy effects [Celsius]
    canopy_height : Union[ndarray, float]
        Canopy height [m]
    measurement_height : Union[ndarray, float]
        Height at which wind speed measurements are done [m]


    Returns
    -------
    Union[ndarray, float] : aerodynamic resistance [s/m]
    """

    # von Karman constant
    k = 0.4
    d = d_ratio * canopy_height
    zom = z0m_ratio * canopy_height
    zoh = z0h_ratio * canopy_height
    logaritmic_profile_ym = numpy.log((measurement_height - d) / zom)

    # Set a minimum speed of 0.1 to avoid math issues
    if numpy.isscalar(wind_speed):
        if wind_speed == 0:
            wind_speed = 0.1
    else:
        wind_speed[wind_speed == 0] = 0.1

    wind_speed_canopy = 1.82 * wind_speed * (numpy.log((measurement_height - d) / zom)) / numpy.log((100.0 - d) / zom)
    a = 1.0591 - 0.0552 * numpy.log(1.72 + numpy.power(4.03 - logaritmic_profile_ym, 2))
    b = 1.9117 - 0.2237 * numpy.log(1.86 + numpy.power(2.12 - logaritmic_profile_ym, 2))
    c = 0.8437 - 0.1243 * numpy.log(3.49 + numpy.power(2.79 - logaritmic_profile_ym, 2))
    ra_prime = numpy.log((measurement_height - d) / zoh) * numpy.log((measurement_height - d) / zom) / (
        k * k * wind_speed_canopy)
    richardson_number = (measurement_height - d) * 9.81 * (air_temperature - canopy_temperature) / (
        wind_speed_canopy * wind_speed_canopy * (air_temperature + 273.15))
    if numpy.isscalar(richardson_number):
        if richardson_number > 0:
            return ra_prime
        else:
            return ra_prime / (a + b * numpy.power(-richardson_number, c))
    else:
        ra_prime[richardson_number <= 0] = ra_prime / (
            a + b * numpy.power(-richardson_number[richardson_number <= 0], c))
        return ra_prime


def aerodynamic_resistance_fao(wind_speed, canopy_height, measurement_height,
                               d_ratio=2/3.0, z0m_ratio=0.123, z0h_ratio=0.0123):
    """
    Function to aerodynamic resistance based on FAO's Penman-Monteith
    Reference: http://www.fao.org/docrep/X0490E/x0490e06.htm#aerodynamic resistance (ra)

    Parameters
    ----------
    d_ratio : float
        Zero-plane displacement (d) = d_ratio * canopy_height
    z0m_ratio : float
        Roughness length (z0) for momentum = z0m_ratio * canopy_height
    z0h_ratio : float
        Roughness length (z0) for heat transfer = z0h_ratio * canopy_height
    wind_speed : Union[ndarray, float]
        Wind speed measured at a given height (measurement_height) [m/s]
    canopy_height : Union[ndarray, float]
        Canopy height [m]
    measurement_height : Union[ndarray, float]
        Height at which wind speed measurements are done [m]

    Returns
    -------
    Union[ndarray, float] : aerodynamic resistance [s/m]
    """

    # von Karman constant
    k = 0.41

    d = d_ratio * canopy_height
    zom = z0m_ratio * canopy_height
    zoh = z0h_ratio * canopy_height


    # Set a minimum speed of 0.1 to avoid math issues
    if numpy.isscalar(wind_speed):
        if wind_speed == 0:
            wind_speed = 0.1
    else:
        wind_speed[wind_speed == 0] = 0.1

    wind_speed_canopy = 1.82 * wind_speed * (numpy.log((measurement_height - d) / zom)) / numpy.log((100.0 - d) / zom)
    ra_s = numpy.log((measurement_height - d) / zoh) * numpy.log((measurement_height - d) / zom) / (
        numpy.power(k, 2) * wind_speed_canopy)
    return ra_s


def raupach_z0_d(canopy_height, area_index, frontal_leaf_area_factor=0.36):
    """
    Calculate z0 and d as a function of height and frontal leaf area, following Raupach (1994)

    Reference: Raupach, M. R. 1994. “Simplified Expressions for Vegetation Roughness Length and Zero-Plane Displacement
    as Functions of Canopy Height and Area Index.” Boundary-Layer Meteorology 71 (1-2): 211–16.


    Parameters
    ----------
    canopy_height : Union[ndarray, float]
        Canopy height [m]
    area_index : Union[ndarray, float]
        Area index
    frontal_leaf_area_factor : float
        Scale factor for converting leaf area index to frontal area index

    Returns
    -------
    tuple: (z0, d) roughness length, zero-plane displacement

    """
    # von Karman constant
    k = 0.41

    c1 = 0.37
    # Free parameter in eq (8)
    cd1 = 7.5
    cw = 2
    cr = 0.35
    cs = 0.003
    gamma_prime_max = 0.3

    # Conversion to frontal leaf area
    frontal_area_index = frontal_leaf_area_factor * area_index

    # Roughness-sublayer influence function. Default value is 0.913.
    psi_h = numpy.log(cw) - 1 + 1 / cw

    gamma_prime = numpy.sqrt(cs + frontal_area_index * cr)

    if numpy.isscalar(gamma_prime):
        if gamma_prime < gamma_prime_max:
            gamma_prime = gamma_prime
        else:
            gamma_prime = gamma_prime_max

    else:
        gamma_prime[gamma_prime > gamma_prime_max] = gamma_prime_max

    d_h = 1 - (1 - numpy.exp(-numpy.sqrt(cd1 * 2 * frontal_area_index))) / numpy.sqrt(cd1 * 2 * frontal_area_index)
    z0mh = (1 - d_h) * numpy.exp(-k * 1/gamma_prime - psi_h)

    raupach_z0 = canopy_height * z0mh
    raupach_d = canopy_height * d_h

    return raupach_z0, raupach_d


def adjust_wind_height(wind_speed, z, measurement_height=2.0, d_ratio=2/3.0, z0m_ratio=0.13, z0h_ratio=0.13):
    """
        Function to aerodynamic resistance based on FAO's Penman-Monteith
        Reference: http://www.fao.org/docrep/X0490E/x0490e06.htm#aerodynamic resistance (ra)

        Parameters
        ----------
        d_ratio : float
            Zero-plane displacement (d) = d_ratio * canopy_height
        z0m_ratio : float
            Roughness length (z0) for momentum = z0m_ratio * canopy_height
        z0h_ratio : float
            Roughness length (z0) for heat transfer = z0h_ratio * canopy_height
        wind_speed : Union[ndarray, float]
            Wind speed measured at a given height (measurement_height) [m/s]
        canopy_height : Union[ndarray, float]
            Canopy height [m]
        measurement_height : Union[ndarray, float]
            Height at which wind speed measurements are done [m]

        Returns
        -------
        Union[ndarray, float] : aerodynamic resistance [s/m]
        """

    # von Karman constant
    k = 0.4
    d = d_ratio * z
    zom = z0m_ratio * z
    zoh = z0h_ratio * z
    logaritmic_profile_ym = numpy.log((measurement_height - d) / zom)

    # Set a minimum speed of 0.1 to avoid math issues
    if numpy.isscalar(wind_speed):
        if wind_speed == 0:
            wind_speed = 0.1
    else:
        wind_speed[wind_speed == 0] = 0.1

    wind_speed_canopy = 1.82 * wind_speed * (numpy.log((measurement_height - d) / zom)) / numpy.log((100.0 - d) / zom)

    return wind_speed_canopy
