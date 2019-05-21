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

__author__ = "Jose A. Jimenez-Berni  <berni@ias.csic.es>"
__version__ = "0.1"
__date__ = "July 2017"

"""
Based on the atmospheric functions from Meteolib:
A library with Python functions for calculations of micro-meteorological parameters.
Waterloo, M. J. 2014. “Meteorology and Evaporation Modules Documentation.”

Release. http://python.hydrology-amsterdam.nl/moduledoc/Meteorologyandevaporationmodules.pdf

The original code has been adapted to be compatible with numpy and pandas.

"""


def es_calc(air_temperature):
    """
    Function to calculate saturated vapour pressure from temperature.
    For T<0 C:  Saturation vapour pressure equation for ice: Goff, J.A.,and S.
    Gratch, Low-pressure properties of water from \-160 to 212 F.
    Transactions of the American society of
    heating and ventilating engineers, pp 95-122, presented
    at the 52nd annual meeting of the American society of
    heating and ventilating engineers, New York, 1946.

    For T>=0 C: Goff, J. A. Saturation pressure of water on the new Kelvin
    temperature scale, Transactions of the American
    society of heating and ventilating engineers, pp 347-354,
    presented at the semi-annual meeting of the American
    society of heating and ventilating engineers, Murray Bay,
    Quebec. Canada, 1957.


    Examples:
        >>> es_calc(30.0)
        >>> 4242.7259946566316
        >>> air_temperature = [20, 25]
        >>> es_calc(air_temperature)
        >>> array([ 2337.080198,  3166.824419])

    Parameters
    ----------
    air_temperature : Union[ndarray, float]
        Measured air temperature [Celsius]

    Returns
    -------
    Union[ndarray, float] : saturated vapour pressure [Pa]
    """

    # Determine length of array
    # Check if we have a scalar value or an array
    if np.isscalar(air_temperature):
        # Calculate saturated vapour pressures, distinguish between water/ice
        if air_temperature < 0:
            # Calculate saturation vapour pressure for ice
            log_pi = - 9.09718 * (273.16 / (air_temperature + 273.15) - 1.0) \
                     - 3.56654 * np.log10(273.16 / (air_temperature + 273.15)) \
                     + 0.876793 * (1.0 - (air_temperature + 273.15) / 273.16) \
                     + np.log10(6.1071)
            es = np.power(10, log_pi)
        else:
            # Calculate saturation vapour pressure for water
            log_pw = 10.79574 * (1.0 - 273.16 / (air_temperature + 273.15)) \
                     - 5.02800 * np.log10((air_temperature + 273.15) / 273.16) \
                     + 1.50475E-4 * (1 - np.power(10, (-8.2969 * ((air_temperature + 273.15) / 273.16 - 1.0)))) + \
                     0.42873E-3 * (np.power(10, (+4.76955 * (1.0 - 273.16 / (air_temperature + 273.15)))) - 1) + \
                     0.78614
            es = np.power(10, log_pw)
    else:  # Dealing with an array

        # Calculate saturated vapour pressures, distinguish between water/ice
        es = np.zeros(air_temperature.shape)

        # Saturation vapour pressure equation for ice
        log_pi = - 9.09718 * (273.16 / (air_temperature + 273.15) - 1.0) \
                 - 3.56654 * np.log10(273.16 / (air_temperature + 273.15)) \
                 + 0.876793 * (1.0 - (air_temperature + 273.15) / 273.16) \
                 + np.log10(6.1071)

        # Calculate saturation vapour pressure for water
        log_pw = 10.79574 * (1.0 - 273.16 / (air_temperature + 273.15)) \
                 - 5.02800 * np.log10((air_temperature + 273.15) / 273.16) \
                 + 1.50475E-4 * (1 - np.power(10, (-8.2969 * ((air_temperature + 273.15) / 273.16 - 1.0)))) + \
                 0.42873E-3 * (np.power(10, (+4.76955 * (1.0 - 273.16 / (air_temperature + 273.15)))) - 1) + \
                 0.78614

        es[np.array(air_temperature < 0, dtype=bool)] = np.power(10, log_pi[np.array(air_temperature < 0, dtype=bool)])
        es[np.array(air_temperature >= 0, dtype=bool)] = np.power(10, log_pw[np.array(air_temperature >= 0, dtype=bool)])
    # Convert from hPa to Pa
    es *= 100.0
    return es


def es_buck(air_temperature):
    """

    Parameters
    ----------
    air_temperature : Union[array, float]
        Moist air temperature [Celsius]. Note: assumes air temperature > 0C

    Returns
    -------
    Union[array, float] : saturated vapour pressure [Pa]

    """
    return 611.21*np.exp((18.678-air_temperature/234.5)*(air_temperature/(257.14+air_temperature)))

def delta_calc(air_temperature):
    """
    Function to calculate the slope of the temperature - vapour pressure curve
    (Delta) from air temperatures. Source: Technical regulations 49, World
    Meteorological Organisation, 1984. Appendix A. 1-Ap-A-3.


    Examples:
        >>> delta_calc(30.0)
        >>> 243.34309166827097
        >>> x = [20, 25]
        >>> delta_calc(x)
        >>> array([ 144.665841,  188.625046])

    Parameters
    ----------
    air_temperature : Union[ndarray, float]
        Air temperature [Celsius]

    Returns
    -------
    ndarray
    """

    # calculate vapour pressure in Pa
    es = es_calc(air_temperature)
    # Convert es (Pa) to kPa
    es /= 1000.0
    # Calculate Delta
    Delta = es * 4098.0 / np.power((air_temperature + 237.3), 2) * 1000.0

    return Delta


def ea_calc(air_temperature, relative_humidity):
    """
    Function to calculate actual saturation vapour pressure.


    Examples:
        >>> ea_calc(25, 60)
        >>> 1900.0946514729308

    Parameters
    ----------
    air_temperature : Union[ndarray, float]
        Air temperature [Celsius]
    relative_humidity : Union[ndarray, float]
        Relative humidity [%]

    Returns
    -------
    Union[ndarray, float] : actual vapour pressure [Pa]
    """

    # Calculate saturation vapour pressures
    es = es_calc(air_temperature)
    # Calculate actual vapour pressure
    eact = relative_humidity / 100.0 * es

    return eact  # in Pa


def vpd_calc(air_temperature, relative_humidity):
    """
    Function to calculate vapour pressure deficit.


    Examples:
        >>> vpd_calc(30,60)
        >>> 1697.0903978626527
        >>> T=[20,25]
        >>> RH=[50,100]
        >>> vpd_calc(T,RH)
        >>> array([ 1168.540099,   0. ])

    Parameters
    ----------
    air_temperature : Union[ndarray, float]
        Air temperature [Celsius]
    relative_humidity : Union[ndarray, float]
        Relative humidity [%]

    Returns
    -------
    Union[ndarray, float] : Vapour pressure deficit [Pa]
    """

    # Calculate saturation vapour pressures
    es = es_calc(air_temperature)
    eact = ea_calc(air_temperature, relative_humidity)
    # Calculate vapour pressure deficit
    return es - eact  # in Pa


def latent_heat_lambda(air_temperature):
    """
    Function to calculate the latent heat of vapourisation,
    lambda, from air temperature. Source: J. Bringfelt. Test of a forest
    evapotranspiration model. Meteorology and Climatology Reports 52,
    SMHI, Norrkopping, Sweden, 1986.

    Input:
        - air_temperature: (array of) air temperature [Celsius]

    Output:
        - L: (array of) lambda [J kg-1 K-1]

    Examples:
        >>> latent_heat_lambda(25)
        >>> 2440883.8804624998
        >>> t=[10, 20, 30]
        >>> latent_heat_lambda(t)
        >>> array([ 2476387.3842125,  2452718.3817125,  2429049.3792125])

    Parameters
    ----------
    air_temperature

    Returns
    -------
    ndarray
    """

    # Calculate lambda
    return 4185.5 * (751.78 - 0.5655 * (air_temperature + 273.15))  # in J/kg


def specific_heat_cp(air_temperature, relative_humidity, air_pressure):
    """
    Function to calculate the specific heat of air, c_p, from air temperatures, relative humidity and air pressure.

    Examples:
        >>> specific_heat_cp(25,60,101300)
        1014.0749457208065
        >>> t=[10, 20, 30]
        >>> rh=[10, 20, 30]
        >>> air_pressure=[100000, 101000, 102000]
        >>> specific_heat_cp(t,rh,air_pressure)
        array([ 1005.13411289,  1006.84399787,  1010.83623841])

    Parameters
    ----------
    air_temperature : Union[ndarray, float]
        Air temperature [Celsius]
    relative_humidity : Union[ndarray, float]
        Relative humidity [%]
    air_pressure : Union[ndarray, float]
        Air pressure [Pa]

    Returns
    -------
    Union[ndarray, float] : specific heat of air c_p [J/kg/K]
    """

    # calculate vapour pressures
    eact = ea_calc(air_temperature, relative_humidity)
    # Calculate cp
    cp = 0.24 * 4185.5 * (1 + 0.8 * (0.622 * eact / (air_pressure - eact)))

    return cp  # in J/kg/K


def psychrometric_gamma(air_temperature, relative_humidity, air_pressure):
    """
    Function to calculate the psychrometric constant gamma.
    Source: J. Bringfelt. Test of a forest evapotranspiration model.
    Meteorology and Climatology Reports 52, SMHI, Norrköpping, Sweden,
    1986.

    Examples:
        >>> psychrometric_gamma(10,50,101300)
        66.263433186572274
        >>> t=[10, 20, 30]
        >>> rh=[10, 20, 30]
        >>> air_pressure=[100000, 101000, 102000]
        >>> psychrometric_gamma(t,rh,air_pressure)
        array([ 65.255188,  66.656958,  68.242393])

    Parameters
    ----------
    air_temperature : Union[ndarray, float]
        Air temperature [Celsius]
    relative_humidity : Union[ndarray, float]
        Relative humidity [%]
    air_pressure : Union[ndarray, float]
        Air pressure [Pa]

    Returns
    -------
    Union[ndarray, float] : psychrometric constant gamma [Pa/K]
    """

    cp = specific_heat_cp(air_temperature, relative_humidity, air_pressure)
    lamb = latent_heat_lambda(air_temperature)
    # Calculate gamma
    gamma = cp * air_pressure / (0.622 * lamb)

    return gamma  # in Pa/K


def air_density_rho(air_temperature, relative_humidity, air_pressure):
    """
    Function to calculate the density of air, rho, from air
    temperatures, relative humidity and air pressure.


    Examples:
        >>> t=[10, 20, 30]
        >>> rh=[10, 20, 30]
        >>> air_pressure=[100000, 101000, 102000]
        >>> air_density_rho(t,rh,air_pressure)
        array([ 1.22948419,  1.19787662,  1.16635358])
        >>> air_density_rho(10,50,101300)
        1.2431927125520903

    Parameters
    ----------
    air_temperature : Union[ndarray, float]
        Air temperature [Celsius]
    relative_humidity : Union[ndarray, float]
        Relative humidity [%]
    air_pressure : Union[ndarray, float]
        Air pressure [Pa]

    Returns
    -------
    Union[ndarray, float] : density of air, rho [kg/m3]
    """

    eact = ea_calc(air_temperature, relative_humidity)
    rho = 1.201 * (290.0 * (air_pressure - 0.378 * eact)) / (1000.0 * (air_temperature + 273.15)) / 100.0

    return rho  # in kg/m3


def potential_temperature_theta(air_temperature, relative_humidity, air_pressure):
    """
    Function to calculate the potential temperature air, theta, from air
    temperatures, relative humidity and air pressure. Reference pressure
    1000 hPa.

    Examples:
        >>> t=[5, 10, 20]
        >>> rh=[45, 65, 89]
        >>> air_pressure=[101300, 102000, 99800]
        >>> potential_temperature_theta(t,rh,air_pressure)
        array([  3.97741582,   8.40874555,  20.16596828])
        >>> potential_temperature_theta(5,45,101300)
        3.977415823848844

    Parameters
    ----------
    air_temperature : Union[ndarray, float]
        Air temperature [Celsius]
    relative_humidity : Union[ndarray, float]
        Relative humidity [%]
    air_pressure : Union[ndarray, float]
        Air pressure [Pa]

    Returns
    -------
    Union[ndarray, float] : potential temperature of air, theta [Celsius]
    """

    # Determine cp
    cp = specific_heat_cp(air_temperature, relative_humidity, air_pressure)
    theta = (air_temperature + 273.15) * np.power((100000.0 / air_pressure), (287.0 / cp)) - 273.15

    return theta  # in degrees celsius


def dew_point(air_temperature, relative_humidity):
    """
    Calculates dew point from air temperature and relative humidity
    Reference: TODO

    Parameters
    ----------
    relative_humidity : ndarray
        Relative humidity [%]
    air_temperature: ndarray
        Dry bulb air temperature [Celsius]

    Returns
    -------
    ndarray: dew point [Celsius]
    """
    a = 17.27
    b = 237.7
    a1 = 17.27
    b1 = 237.7
    # TODO: Check if we can use gamma function instead
    gamma = a1 * air_temperature / (b1 + air_temperature) + np.log(relative_humidity / 100.0)

    dew = b * gamma / (a - gamma)
    return dew


def rho_cp_dry(air_temperature, air_pressure):
    """

    Parameters
    ----------
    air_temperature
    air_pressure

    Returns
    -------

    """
    rho = air_temperature / (287.07 * (air_temperature + 273.15))
    return rho * (1002.5 + 275e-6 * np.power(air_temperature + 273.15 - 200.0, 2))


def rho_cp_moist(air_temperature, relative_humidity, air_pressure):
    """

    Parameters
    ----------
    air_temperature
    relative_humidity
    air_pressure

    Returns
    -------

    """
    rho = air_density_rho(air_temperature, relative_humidity, air_pressure)
    return rho * specific_heat_cp(air_temperature, relative_humidity, air_pressure)
