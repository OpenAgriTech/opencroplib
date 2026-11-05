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

import numpy
from .atmophere import dew_point, ea_calc

__author__ = "Jose A. Jimenez-Berni  <berni@ias.csic.es>"
__version__ = "0.1.3"
__date__ = "July 2019"

BOLTZMANN = 5.670373e-8  # [W m-2 K-4]
BOLTZMANN_HOURLY = 2.043e-10  # [MJ K-4 m-2 hour-1] (FAO-56: 4.903e-9 / 24)
STEFAN_BOLTZMANN_CONSTANT = 0.000000004903  # [MJ K-4 m-2 day-1]
SUPPORTED_ATMOSPHERIC_EMISSIVITY_MODELS = ['B75', 'Swinbank', 'BM']
#: Solar constant [ MJ m-2 min-1]
SOLAR_CONSTANT = 0.0820


def planck(temperature_kelvin, wavelength=10.25):
    """
    Calculate Planck's black body thermal radiance for a given temperature and central wavelength
    Reference: https://ncc.nesdis.noaa.gov/data/planck.html

    Parameters
    ----------
    wavelength: float
        Central wavelength used for the calculation [micrometers]. Default is 10.25um
    temperature_kelvin: float
        Object temperature [Kelvin]

    Returns
    -------
    float
        Spectral radiance of the black body with a given temperature at a given wavelength [W m-2 sr-1 um-1].

    """

    c1 = 1.1910427e-22
    c2 = 1.4387752e-2
    cwl = wavelength * 1e-6

    return c1 / numpy.power(cwl, 5) / (numpy.exp(c2 / temperature_kelvin / cwl) - 1.0)


def iplanck(radiance, wavelength=10.25):
    """
    Calculate the inverse of Planck's function for a given thermal radiance and central wavelength
    Reference: https://ncc.nesdis.noaa.gov/data/planck.html

    Parameters
    ----------
    wavelength : float
        Central wavelength used for the calculation [micrometers]. Default is 10.25um
    radiance : float
        Spectral radiance of the black body with a given temperature at a given wavelength [W m-2 sr-1 um-1].

    Returns
    -------
    float: Object temperature [Kelvin]

    """

    h = 6.6260755e-34
    c = 2.9979246e8
    k = 1.380658e-23

    cwl = wavelength * 1e-6

    return h * c / (k * cwl) / numpy.log((2.0 * h * c * c) / (radiance * 1e6 * numpy.power(cwl, 5)) + 1.0)


def longwave_downwelling(air_temperature, relative_humidity, cloudiness_factor=0, emissivity_model='B75'):
    """
    Function to calculate the longwave downwelling radiation from the air temperature and relative humidity
    includes a cloudiness factor for non-clear-sky conditions

    Reference: Crawford, Todd M., and Claude E. Duchon. 1999. “An Improved Parameterization for Estimating Effective
    Atmospheric Emissivity for Use in Calculating Daytime Downwelling Longwave Radiation.”
    Journal of Applied Meteorology 38 (4). American Meteorological Society: 474–80.

    See also for the formulations of the emissivity models:
    Alados, I., I. Foyo-Moreno, and L. Alados-Arboledas. 2012. “Estimation of Downwelling Longwave Irradiance under
    All-Sky Conditions.” International Journal of Climatology 32 (5). John Wiley & Sons, Ltd.: 781–93.

    Parameters
    ----------
    cloudiness_factor : Union[ndarray, float]
        Cloudiness factor [0-1] for non-clear sky conditions. 0=clear sky
    emissivity_model : str
        Choose between ['B75', 'Swinbank', 'BM']
    relative_humidity : Union[ndarray, float]
        Relative humidity [%]
    air_temperature : Union[ndarray, float]
        Air temperature [Celsius]

    Returns
    -------
    Union[ndarray, float]: longwave downwelling radiation [W/m2]

    """

    if emissivity_model not in SUPPORTED_ATMOSPHERIC_EMISSIVITY_MODELS:
        raise ValueError('{foo} model not supported. Choose between {supported}'.format(foo=repr(emissivity_model),
                                                                                        supported=SUPPORTED_ATMOSPHERIC_EMISSIVITY_MODELS))

    ea_t = ea_calc(air_temperature, relative_humidity) * 0.01

    atmospheric_emissivity = 1

    if emissivity_model == 'B75':
        atmospheric_emissivity = (cloudiness_factor + (1 - cloudiness_factor) * (
                1.24 * numpy.power(ea_t / (air_temperature + 273.15), 1.0 / 7.0)))
    elif emissivity_model == 'Swinbank':
        atmospheric_emissivity = (
                cloudiness_factor + (1 - cloudiness_factor) * (9.36 * 1e-6 * numpy.power(air_temperature + 273.15, 2)))
    elif emissivity_model == 'BM':
        dew_t = dew_point(air_temperature, relative_humidity) / 100
        a = 0.711
        b = 0.56
        c = 0.73
        atmospheric_emissivity = (
                cloudiness_factor + (1 - cloudiness_factor) * (a + b * dew_t + c * numpy.power(dew_t, 2)))

    lwd_sim = atmospheric_emissivity * BOLTZMANN * numpy.power(air_temperature + 273.15, 4)

    return lwd_sim


def available_energy(shortwave, albedo, longwave, canopy_temperature, emissivity=0.98):
    """
    Function to calculate the available energy based on the balance between shortwave and longwave

    Reference...


    Parameters
    ----------
    emissivity : float
        Object emissivity (e.g. 0.98 for vegetation). It is assumed that emissivity and absorptivity are equal.
    canopy_temperature : Union[ndarray, float]
        Canopy temperature used to calculate upwelling longwave radiation [Celsius]
    longwave : Union[ndarray, float]
        Downwelling longwave radiation (see longwave_downwelling functions) [W/m2]
    albedo : float
        Surface albedo
    shortwave : Union[ndarray, float]
        Downwelling shortwave radiation (from pyranometer)

    Returns
    -------
    Union[ndarray, float]: total available energy or net radiation [W/m2]

    """
    return (1.0 - albedo) * shortwave + emissivity * (
                longwave - BOLTZMANN * numpy.power(canopy_temperature + 273.15, 4))


def solar_declination(doy):
    """
    Calculates solar declination from the day of the year using the Fourier series of Spencer (1971), as given in
    Iqbal (1983). Maximum error ~0.035 deg.

    The previous formula, arcsin(0.409 sin(2 pi (doy - 82) / 365)), used the obliquity in radians (0.409) as the
    amplitude of sin(declination), giving a maximum declination of 24.1 deg instead of 23.44 deg.

    References:
        Spencer, J.W. 1971. Fourier series representation of the position of the sun. Search 2 (5): 172.
        Iqbal, M. 1983. An Introduction to Solar Radiation. Academic Press, Eq. 1.3.1.

    Parameters
    ----------
    doy : Union[ndarray, float]
        day of the year [1-366]

    Returns
    -------
    Union[ndarray, float] : solar declination [radians]

    """
    day_angle = 2.0 * numpy.pi * (doy - 1.0) / 365.0
    return (0.006918 - 0.399912 * numpy.cos(day_angle) + 0.070257 * numpy.sin(day_angle)
            - 0.006758 * numpy.cos(2.0 * day_angle) + 0.000907 * numpy.sin(2.0 * day_angle)
            - 0.002697 * numpy.cos(3.0 * day_angle) + 0.00148 * numpy.sin(3.0 * day_angle))



def solar_declination_fao(doy):
    """
    Calculates solar declination from the day of the year according to FAO equation 24 in Allen et al (1998).

    Parameters
    ----------
    doy : Union[ndarray, float]
        day of the year

    Returns
    -------
    Union[ndarray, float] : solar declination [radians]

    """

    return 0.409 * numpy.sin(((2.0 * numpy.pi / 365.0) * doy - 1.39))


def solar_elevation(latitude, doy, solar_h):
    """
    Calculates the (geometric, not refracted) solar elevation for a given latitude, day of the year and solar time:

        sin(h) = cos(lat) cos(dec) cos(w) + sin(lat) sin(dec)

    with the declination from solar_declination and the hour angle w = 15 deg (solar_h - 12). Compared with the
    NREL SPA algorithm the error is ~0.2 deg RMS (< 0.6 deg) when solar_h comes from time_utils.get_solar_time.

    Reference: Iqbal, M. 1983. An Introduction to Solar Radiation. Academic Press, Eq. 1.5.1.

    Parameters
    ----------
    latitude : Union[ndarray, float]
        Latitude in degrees (+ North / - South)
    doy : Union[ndarray, float]
        Day of the year [0-365]
    solar_h : Union[ndarray, float]
        Decimal hour of the day [0.0-24.0]

    Returns
    -------
    Union[ndarray, float] : solar elevation [radians]

    """
    lat_r = numpy.radians(latitude)
    dec = solar_declination(doy)
    solar_h_r = numpy.radians((solar_h - 12.0) * 15.0)
    sin_h = numpy.cos(lat_r) * numpy.cos(dec) * numpy.cos(solar_h_r) + numpy.sin(lat_r) * numpy.sin(dec)
    return numpy.arcsin(sin_h)


def solar_azimuth(latitude, doy, solar_h):
    """
    Calculates the solar azimuth for a given latitude, day of the year and solar time.

    Convention (Duffie and Beckman): angle from due south, negative towards the east (morning) and positive towards
    the west (afternoon), in [-pi, pi]. The equivalent azimuth clockwise from north is solar_azimuth + pi.
    The previous implementation returned the unsigned angle (arccos), so morning and afternoon positions could not
    be distinguished.

        tan(azimuth) = sin(w) / (sin(lat) cos(w) - cos(lat) tan(dec))

    Reference: Duffie, J.A., Beckman, W.A. 2013. Solar Engineering of Thermal Processes, 4th ed. Wiley, Eq. 1.6.6.

    Parameters
    ----------
    latitude : Union[ndarray, float]
        Latitude in degrees (+ North / - South)
    doy : Union[ndarray, float]
        Day of the year [0-365]
    solar_h : Union[ndarray, float]
        Decimal hour of the day [0.0-24.0]

    Returns
    -------
    Union[ndarray, float] : solar azimuth from south, west positive [radians]

    """
    lat_r = numpy.radians(latitude)
    dec = solar_declination(doy)
    solar_h_r = numpy.radians((solar_h - 12.0) * 15.0)
    west = numpy.cos(dec) * numpy.sin(solar_h_r)
    south = numpy.sin(lat_r) * numpy.cos(dec) * numpy.cos(solar_h_r) - numpy.cos(lat_r) * numpy.sin(dec)
    return numpy.arctan2(west, south)



def toa_irradiance(doy):
    """
    Calculates top of atmosphere irradiance for the day of the year

    Reference...


    Parameters
    ----------
    doy : Union[ndarray, float]
        Day of the year [0-365]

    Returns
    -------
    Union[ndarray, float] : top of atmosphere terrestrial irradiance [W/m2]

    """
    solar_constant = 1367.0
    jp = 2.0 * numpy.pi * doy / 365.25
    ep = 1.0 + 0.03344 * numpy.cos(jp - 0.048869)
    return solar_constant * ep


def _float_if_scalar(values):
    return float(values) if numpy.ndim(values) == 0 else values


def clear_sky_beam(latitude, doy, solar_h, altitude=0, t_lk=1):
    """
    Clear-sky beam (direct) irradiance on a horizontal surface from the ESRA clear-sky model (Rigollier et al. 2000)
    as implemented in GRASS r.sun (Hofierka and Suri 2002):

        B = I0 * sin(h0) * exp(-0.8662 * TL * m * dR(m))

    with the refraction-corrected solar altitude, the Kasten and Young (1989) relative optical air mass m corrected
    for the site altitude, and the Rayleigh optical thickness dR(m).

    References:
        Rigollier, C., Bauer, O., Wald, L. 2000. On the clear sky model of the ESRA - European Solar Radiation Atlas -
        with respect to the Heliosat method. Solar Energy 68 (1): 33-48.
        Hofierka, Jaroslav, and Marcel Suri. 2002. “The Solar Radiation Model for Open Source GIS: Implementation and
        Applications.” In Proceedings of the Open Source GIS - GRASS Users Conference.

    Parameters
    ----------
    latitude : Union[ndarray, float]
        Latitude in degrees (+ North / - South)
    doy : Union[ndarray, float]
        Day of the year [0-365]
    solar_h : Union[ndarray, float]
        Decimal hour of the day [0.0-24.0]
    altitude : Union[ndarray, float]
        Ground altitude [m]. Default is sea level (0m)
    t_lk : Union[ndarray, float]
        Linke turbidity factor for an air mass of 2 [-]. Default 1 (clean, dry atmosphere)

    Returns
    -------
    Union[ndarray, float] : clear-sky beam irradiance on a horizontal surface [W/m2]. 0 when the sun is below the
    horizon.
    """
    s_elev = solar_elevation(latitude, doy, solar_h)
    ext_rad = toa_irradiance(doy)

    # Atmospheric refraction correction; solar altitude and correction in radians (ESRA / r.sun)
    dh = 0.061359 * (0.1594 + 1.123 * s_elev + 0.065656 * s_elev * s_elev) / (
            1.0 + 28.9344 * s_elev + 277.3971 * s_elev * s_elev)
    href = s_elev + dh
    sun_up = numpy.asarray(s_elev > 0)
    # Only evaluated where the sun is up; the clip avoids invalid powers of negative numbers below the horizon
    href_d = numpy.where(sun_up, numpy.degrees(href), 90.0)

    # Relative optical air mass (Kasten and Young 1989), corrected for the site altitude
    p_p0 = numpy.exp(-altitude / 8434.5)
    m = p_p0 / (numpy.sin(numpy.radians(href_d)) + 0.50572 * numpy.power(href_d + 6.07995, -1.6364))

    # Rayleigh optical thickness at air mass m
    dr = numpy.where(m <= 20.0,
                     1.0 / (6.6296 + 1.7513 * m - 0.1202 * m ** 2 + 0.0065 * m ** 3 - 0.00013 * m ** 4),
                     1.0 / (10.4 + 0.718 * m))

    beam = numpy.where(sun_up, ext_rad * numpy.exp(-0.8662 * t_lk * m * dr) * numpy.sin(s_elev), 0.0)
    return _float_if_scalar(beam)


def clear_sky_diffuse(latitude, doy, solar_h, t_lk=1):
    """
    Clear-sky diffuse irradiance on a horizontal surface from the ESRA clear-sky model (Rigollier et al. 2000) as
    implemented in GRASS r.sun (Hofierka and Suri 2002):

        D = I0 * Trd(TL) * Fd(h0, TL)

    with the diffuse transmission function at zenith Trd(TL) and the diffuse angular function
    Fd = A1 + A2 sin(h0) + A3 sin(h0)^2.

    References: see clear_sky_beam.

    Parameters
    ----------
    latitude : Union[ndarray, float]
        Latitude in degrees (+ North / - South)
    doy : Union[ndarray, float]
        Day of the year [0-365]
    solar_h : Union[ndarray, float]
        Decimal hour of the day [0.0-24.0]
    t_lk : Union[ndarray, float]
        Linke turbidity factor for an air mass of 2 [-]. Default 1 (clean, dry atmosphere)

    Returns
    -------
    Union[ndarray, float] : clear-sky diffuse irradiance on a horizontal surface [W/m2]. 0 when the sun is below the
    horizon.
    """
    s_elev = solar_elevation(latitude, doy, solar_h)
    sin_h = numpy.sin(s_elev)
    ext_rad = toa_irradiance(doy)

    trd = -1.5843e-2 + 3.0543e-2 * t_lk + 3.797e-4 * t_lk ** 2
    a1_prime = 2.6463e-1 - 6.1581e-2 * t_lk + 3.1408e-3 * t_lk ** 2
    a1 = numpy.where(a1_prime * trd < 2.2e-3, 2.2e-3 / trd, a1_prime)
    a2 = 2.04020 + 1.8945e-2 * t_lk - 1.1161e-2 * t_lk ** 2
    a3 = -1.3025 + 3.9231e-2 * t_lk + 8.5079e-3 * t_lk ** 2
    fd = a1 + a2 * sin_h + a3 * sin_h ** 2

    diffuse = numpy.where(numpy.asarray(s_elev > 0), ext_rad * trd * fd, 0.0)
    return _float_if_scalar(diffuse)


def potential_irradiance(latitude, doy, solar_h, altitude=0, t_lk=1):
    """
    Calculates the clear-sky global (beam + diffuse) irradiance on a horizontal surface for a given latitude, day of
    the year, solar time, altitude and Linke turbidity, using the ESRA clear-sky model (Rigollier et al. 2000) as
    implemented in GRASS r.sun (Hofierka and Suri 2002). See clear_sky_beam and clear_sky_diffuse for the
    components.

    References:
        Rigollier, C., Bauer, O., Wald, L. 2000. On the clear sky model of the ESRA - European Solar Radiation Atlas -
        with respect to the Heliosat method. Solar Energy 68 (1): 33-48.
        Hofierka, Jaroslav, and Marcel Suri. 2002. “The Solar Radiation Model for Open Source GIS: Implementation and
        Applications.” In Proceedings of the Open Source GIS - GRASS Users Conference.
        http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.19.9831.


    Parameters
    ----------
    altitude : Union[ndarray, float]
        Ground altitude [m]. Default is sea level (0m). Affects the beam component only.
    t_lk : Union[ndarray, float]
        Linke turbidity factor for an air mass of 2 [-]. Default 1 (clean, dry atmosphere); typically 2-3 for clear
        rural skies and 4-6 for turbid or polluted skies.
    latitude : Union[ndarray, float]
        Latitude in degrees (+ North / - South)
    doy : Union[ndarray, float]
        Day of the year [0-365]
    solar_h : Union[ndarray, float]
        Decimal hour of the day [0.0-24.0]

    Returns
    -------
    Union[ndarray, float] : clear-sky global irradiance on a horizontal surface [W/m2]. 0 when the sun is below the
    horizon.

    """
    return _float_if_scalar(clear_sky_beam(latitude, doy, solar_h, altitude, t_lk) +
                            clear_sky_diffuse(latitude, doy, solar_h, t_lk))


def daily_potential_toa(latitude, day_of_year):
    """
    Estimate daily extraterrestrial radiation (*Ra*, 'top of the atmosphere
    radiation').
    Based on equation 21 in Allen et al (1998). If monthly mean radiation is
    required make sure *sol_dec*. *sha* and *irl* have been calculated using
    the day of the year that corresponds to the middle of the month.
    **Note**: From Allen et al (1998): "For the winter months in latitudes
    greater than 55 degrees (N or S), the equations have limited validity.
    Reference should be made to the Smithsonian Tables to assess possible
    deviations."
    :param latitude: Latitude [degrees] (+ North / - South)
    :param day_of_year: Day of the year

    :return: Daily extraterrestrial radiation [MJ m-2 day-1]
    :rtype: float
    """

    lat_r = numpy.radians(latitude)

    sol_dec = solar_declination_fao(day_of_year)
    sha = sunset_hour_angle(latitude, sol_dec)
    ird = inv_rel_dist_earth_sun(day_of_year)

    tmp1 = (24.0 * 60.0) / numpy.pi
    tmp2 = sha * numpy.sin(lat_r) * numpy.sin(sol_dec)
    tmp3 = numpy.cos(lat_r) * numpy.cos(sol_dec) * numpy.sin(sha)
    return tmp1 * SOLAR_CONSTANT * ird * (tmp2 + tmp3)


def daily_clear_sky_irradiance(altitude, et_rad):
    """
    Estimate clear sky radiation from altitude and extraterrestrial radiation.
    Based on equation 37 in Allen et al (1998) which is recommended when
    calibrated Angstrom values are not available.
    :param altitude: Elevation above sea level [m]
    :param et_rad: Extraterrestrial radiation [MJ m-2 day-1]. Can be
        estimated using ``et_rad()``.
    :return: Clear sky radiation [MJ m-2 day-1]
    :rtype: float
    """
    return (0.00002 * altitude + 0.75) * et_rad


def inv_rel_dist_earth_sun(day_of_year):
    """
    Calculate the inverse relative distance between earth and sun from
    day of the year.
    Based on FAO equation 23 in Allen et al (1998).
    :param day_of_year: Day of the year [1 to 366]
    :return: Inverse relative distance between earth and the sun
    :rtype: float
    """

    return 1 + (0.033 * numpy.cos((2.0 * numpy.pi / 365.0) * day_of_year))


def sunset_hour_angle(latitude, sol_dec):
    """
    Calculate sunset hour angle (*Ws*) from latitude and solar
    declination.
    Based on FAO equation 25 in Allen et al (1998).
    :param latitude: Latitude [degrees]. Note: *latitude* should be negative
        if it in the southern hemisphere, positive if in the northern
        hemisphere.
    :param sol_dec: Solar declination [radians]. Can be calculated using
        ``sol_dec()``.
    :return: Sunset hour angle [radians].
    :rtype: float
    """

    lat_r = numpy.radians(latitude)

    cos_sha = -numpy.tan(lat_r) * numpy.tan(sol_dec)
    # If tmp is >= 1 there is no sunset, i.e. 24 hours of daylight
    # If tmp is <= 1 there is no sunrise, i.e. 24 hours of darkness
    # See http://www.itacanet.org/the-sun-as-a-source-of-energy/
    # part-3-calculating-solar-angles/
    # Domain of acos is -1 <= x <= 1 radians (this is not mentioned in FAO-56!)
    return numpy.arccos(numpy.clip(cos_sha, -1.0, 1.0))


def net_out_lw_daily(tmin, tmax, sol_rad, cs_rad, avp):
    """
    Estimate net outgoing longwave radiation.
    This is the net longwave energy (daily net energy flux) leaving the
    earth's surface. It is proportional to the absolute temperature of
    the surface raised to the fourth power according to the Stefan-Boltzmann
    law. However, water vapour, clouds, carbon dioxide and dust are absorbers
    and emitters of longwave radiation. This function corrects the Stefan-
    Boltzmann law for humidity (using actual vapor pressure) and cloudiness
    (using solar radiation and clear sky radiation). The concentrations of all
    other absorbers are assumed to be constant.
    The output can be converted to equivalent evaporation [mm day-1] using
    ``energy2evap()``.
    Based on FAO equation 39 in Allen et al (1998).
    :param tmin: Absolute daily minimum temperature [degrees Celsius]
    :param tmax: Absolute daily maximum temperature [degrees Celsius]
    :param sol_rad: Solar radiation [MJ m-2 day-1]. If necessary this can be
        estimated using ``sol_rad()``.
    :param cs_rad: Clear sky radiation [MJ m-2 day-1]. Can be estimated using
        ``cs_rad()``.
    :param avp: Actual vapour pressure [kPa]. Can be estimated using functions
        with names beginning with 'avp_from'.
    :return: Net outgoing longwave radiation [MJ m-2 day-1]
    :rtype: float
    """

    tmp1 = (STEFAN_BOLTZMANN_CONSTANT *
            ((numpy.power(tmax+273.15, 4) + numpy.power(tmin+273.15, 4)) / 2))
    tmp2 = (0.34 - (0.14 * numpy.sqrt(avp)))
    tmp3 = 1.35 * (sol_rad / cs_rad) - 0.35
    return tmp1 * tmp2 * tmp3


def net_out_lw_hourly(t_mean, relative_humidity, daily_clearness):
    """
    Estimate net outgoing longwave radiation.
    This is the net longwave energy (hourly net energy flux) leaving the
    earth's surface. It is proportional to the absolute temperature of
    the surface raised to the fourth power according to the Stefan-Boltzmann
    law. However, water vapour, clouds, carbon dioxide and dust are absorbers
    and emitters of longwave radiation. This function corrects the Stefan-
    Boltzmann law for humidity (using actual vapor pressure) and cloudiness
    (using solar radiation and clear sky radiation). The concentrations of all
    other absorbers are assumed to be constant.
    The output can be converted to equivalent evaporation [mm day-1] using
    ``energy2evap()``.
    Based on FAO equation 39 in Allen et al (1998), with the hourly Stefan-Boltzmann constant
    (2.043e-10 MJ K-4 m-2 hour-1) as in FAO-56 Example 19.
    :param t_mean: Hourly mean temperature [degrees Celsius]
    :param relative_humidity: Hourly mean relative humidity [%]
    :param daily_clearness: clearness index calculated as Rs/Rs_clear_sky.
    :return: Net outgoing longwave radiation [MJ m-2 hour-1]
    :rtype: float
    """

    # Actual vapour pressure in kPa (ea_calc returns Pa)
    avp = ea_calc(t_mean, relative_humidity) / 1000.0
    tmp1 = (BOLTZMANN_HOURLY * numpy.power((t_mean+273.15), 4))
    tmp2 = (0.34 - (0.14 * numpy.sqrt(avp)))
    tmp3 = 1.35 * daily_clearness - 0.35
    return tmp1 * tmp2 * tmp3


def net_in_sol_rad(sol_rad, albedo=0.23):
    """
    Calculate net incoming solar (or shortwave) radiation from gross
    incoming solar radiation, assuming a grass reference crop.
    Net incoming solar radiation is the net shortwave radiation resulting
    from the balance between incoming and reflected solar radiation. The
    output can be converted to equivalent evaporation [mm day-1] using
    ``energy2evap()``.
    Based on FAO equation 38 in Allen et al (1998).
    :param sol_rad: Gross incoming solar radiation [MJ m-2 day-1]. If
        necessary this can be estimated using functions whose name
        begins with 'sol_rad_from'.
    :param albedo: Albedo of the crop as the proportion of gross incoming solar
        radiation that is reflected by the surface. Default value is 0.23,
        which is the value used by the FAO for a short grass reference crop.
        Albedo can be as high as 0.95 for freshly fallen snow and as low as
        0.05 for wet bare soil. A green vegetation over has an albedo of
        about 0.20-0.25 (Allen et al, 1998).
    :return: Net incoming solar (or shortwave) radiation [MJ m-2 day-1].
    :rtype: float
    """
    return (1 - albedo) * sol_rad


def net_rad(ni_sw_rad, no_lw_rad):
    """
    Calculate daily net radiation at the crop surface, assuming a grass
    reference crop.
    Net radiation is the difference between the incoming net shortwave (or
    solar) radiation and the outgoing net longwave radiation. Output can be
    converted to equivalent evaporation [mm day-1] using ``energy2evap()``.
    Based on equation 40 in Allen et al (1998).
    :param ni_sw_rad: Net incoming shortwave radiation [MJ m-2 day-1]. Can be
        estimated using ``net_in_sol_rad()``.
    :param no_lw_rad: Net outgoing longwave radiation [MJ m-2 day-1]. Can be
        estimated using ``net_out_lw_rad()``.
    :return: Daily net radiation [MJ m-2 day-1].
    :rtype: float
    """
    return ni_sw_rad - no_lw_rad


# Below, alternative formulations for the downwelling longwave radiation

def longwave_downwelling_1(air_temperature, relative_humidity):
    """
    Function to calculate the longwave downwelling radiation from the air temperature and relative humidity, with
    the Brunt-type clear-sky emissivity implied by FAO-56:

        emissivity = 0.64 + 0.044 * sqrt(ea),  ea in hPa

    FAO-56 Eq. 39 writes the net emissivity between the surface and the atmosphere as (0.34 - 0.14 sqrt(ea)),
    ea in kPa. With a surface emissivity of 0.98 (the default of available_energy), the atmospheric emissivity is
    0.98 - 0.34 + 0.14 sqrt(ea[kPa]) = 0.64 + 0.0443 sqrt(ea[hPa]). (Identified from the coefficients; the original
    implementation gave no reference.)

    References:
        Allen, R.G., Pereira, L.S., Raes, D., Smith, M. 1998. Crop evapotranspiration. FAO Irrigation and Drainage
        Paper 56, Eq. 39.
        Brunt, D. 1932. Notes on radiation in the atmosphere. Quarterly Journal of the Royal Meteorological Society
        58: 389-420 (functional form; original coefficients 0.52 and 0.065).

    Parameters
    ----------
    relative_humidity : float
     Relative humidity [%]
    air_temperature : float
     Air temperature [Celsius]

    Returns
    -------
    float: longwave downwelling radiation [W/m2]

    """

    ea_t = ea_calc(air_temperature, relative_humidity) * 0.01

    lwd_sim = (0.64 + 0.044 * numpy.sqrt(ea_t)) * \
              BOLTZMANN * numpy.power(air_temperature + 273.15, 4)
    return lwd_sim


def longwave_downwelling_2(air_temperature, relative_humidity):
    """
    Function to calculate the longwave downwelling radiation from the air temperature and relative humidity, with
    the clear-sky atmospheric emissivity of Idso (1981): 0.70 + 5.95e-5 ea exp(1500 / T), ea in hPa, T in K.

    Reference: Idso, S.B. 1981. A set of equations for full spectrum and 8- to 14-um and 10.5- to 12.5-um thermal
    radiation from cloudless skies. Water Resources Research 17 (2): 295-304.

    Parameters
    ----------
    relative_humidity : float
        Relative humidity [%]
    air_temperature : float
        Air temperature [Celsius]

    Returns
    -------
    float: longwave downwelling radiation [W/m2]

    """

    ea_t = ea_calc(air_temperature, relative_humidity) * 0.01

    lwd_sim = (0.7 + 0.0000595 * ea_t * numpy.exp(1500.0 / (air_temperature + 273.15))) * \
              BOLTZMANN * numpy.power(air_temperature + 273.15, 4)
    return lwd_sim


def longwave_downwelling_3(air_temperature, relative_humidity):
    """
    Function to calculate the longwave downwelling radiation from the air temperature and relative humidity, for a
    sky with a fixed cloudiness of 0.3:

        emissivity = 0.3 + 0.7 * 1.24 * (ea / T) ** (1/7),  ea in hPa, T in K

    This is the Brutsaert (1975) clear-sky emissivity with the Crawford and Duchon (1999) cloud correction,
    emissivity = clf + (1 - clf) * emissivity_clear, for clf = 0.3. It is identical to
    longwave_downwelling(air_temperature, relative_humidity, cloudiness_factor=0.3, emissivity_model='B75').
    (Identified from the formula; the original implementation gave no reference.)

    References: see longwave_downwelling.

    Parameters
    ----------
    relative_humidity : float
        Relative humidity [%]
    air_temperature : float
        Air temperature [Celsius]

    Returns
    -------
    float: longwave downwelling radiation [W/m2]

    """

    ea_t = ea_calc(air_temperature, relative_humidity) * 0.01

    lwd_sim = (0.3 + 0.7 * 1.24 * numpy.power(ea_t / (air_temperature + 273.15), 1.0 / 7.0)) * \
              BOLTZMANN * numpy.power(air_temperature + 273.15, 4)
    return lwd_sim