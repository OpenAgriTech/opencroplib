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

__author__ = "Jose A. Jimenez-Berni  <berni@ias.csic.es>"
__version__ = "0.1"
__date__ = "July 2017"

# von Karman constant, used by every function in this module (FAO-56 value)
VON_KARMAN = 0.41

# Default FAO-56 reference grass and blending height used to transfer a wind speed measured over grass to the crop
FAO_GRASS_HEIGHT = 0.12
FAO_GRASS_MEASUREMENT_HEIGHT = 2.0
BLENDING_HEIGHT = 100.0
FAO_GRASS_D_RATIO = 2 / 3.0
FAO_GRASS_Z0M_RATIO = 0.123

# Wind speeds below this value are clamped to avoid division by (near) zero [m/s]
MIN_WIND_SPEED = 0.1


def _check_above_roughness(height, d, z0m, name):
    """Raise ValueError when the log-profile at `height` is undefined (height - d <= z0m)."""
    if numpy.any(numpy.asarray(height - d) <= numpy.asarray(z0m)):
        raise ValueError("{name} - d must be greater than z0m for the logarithmic wind profile "
                         "(is the height inside or too close to the canopy?)".format(name=name))


def grass_to_blending_wind_factor(grass_height=FAO_GRASS_HEIGHT,
                                  grass_measurement_height=FAO_GRASS_MEASUREMENT_HEIGHT,
                                  blending_height=BLENDING_HEIGHT,
                                  d_ratio=FAO_GRASS_D_RATIO, z0m_ratio=FAO_GRASS_Z0M_RATIO):
    """
    Ratio between the wind speed at the blending height and the wind speed measured over reference grass,
    assuming a neutral logarithmic profile over the grass:

        u(blending_height) / u(grass_measurement_height) = ln((zb - d_g) / z0m_g) / ln((zm_g - d_g) / z0m_g)

    With the FAO-56 defaults (h = 0.12 m, measured at 2 m, blending height 100 m) this is ~1.81.

    Parameters
    ----------
    grass_height : float
        Height of the grass over which wind is measured [m]
    grass_measurement_height : float
        Height of the anemometer over the grass [m]
    blending_height : float
        Height at which the wind is assumed to be unaffected by the surface [m]
    d_ratio : float
        Zero-plane displacement of the grass = d_ratio * grass_height
    z0m_ratio : float
        Roughness length for momentum of the grass = z0m_ratio * grass_height

    Returns
    -------
    float : wind speed ratio [-]
    """
    d = d_ratio * grass_height
    z0m = z0m_ratio * grass_height
    _check_above_roughness(grass_measurement_height, d, z0m, "grass_measurement_height")
    _check_above_roughness(blending_height, d, z0m, "blending_height")
    return numpy.log((blending_height - d) / z0m) / numpy.log((grass_measurement_height - d) / z0m)


def wind_speed_2m(wind_speed, measurement_height, grass_height=FAO_GRASS_HEIGHT, d_ratio=FAO_GRASS_D_RATIO,
                  z0m_ratio=FAO_GRASS_Z0M_RATIO):
    """
    Converts a wind speed measured at `measurement_height` over reference grass to the standard 2 m height, assuming
    a neutral logarithmic profile over the grass:

        u2 = uz * ln((2 - d) / z0m) / ln((z - d) / z0m)

    With the FAO-56 defaults (h = 0.12 m, d = 2/3 h, z0m = 0.123 h) this is FAO-56 Eq. 47,
    u2 = uz * 4.87 / ln(67.8 z - 5.42), to within 0.1 %.

    Reference: Allen, R.G., Pereira, L.S., Raes, D., Smith, M., 1998. Crop evapotranspiration. FAO Irrigation and
    Drainage Paper 56, Eq. 47.

    Parameters
    ----------
    wind_speed : Union[ndarray, float]
        Wind speed measured at measurement_height over grass [m/s]
    measurement_height : Union[ndarray, float]
        Anemometer height over the grass [m]
    grass_height : float
        Height of the reference grass [m]. Default 0.12 (FAO-56)
    d_ratio : float
        Zero-plane displacement of the grass = d_ratio * grass_height
    z0m_ratio : float
        Roughness length for momentum of the grass = z0m_ratio * grass_height

    Returns
    -------
    Union[ndarray, float] : wind speed at 2 m over grass [m/s]

    Raises
    ------
    ValueError
        If measurement_height - d <= z0m (anemometer inside or too close to the grass).
    """
    d = d_ratio * grass_height
    z0m = z0m_ratio * grass_height
    _check_above_roughness(measurement_height, d, z0m, "measurement_height")
    return wind_speed * numpy.log((2.0 - d) / z0m) / numpy.log((measurement_height - d) / z0m)


def _wind_at_reference_height(wind_speed, canopy_height, measurement_height, d_ratio, z0m_ratio, wind_over_grass,
                              grass_height, grass_measurement_height, blending_height, min_wind_speed):
    """
    Wind speed at `measurement_height` over the crop, plus d and z0m of the crop.

    The input is never modified: speeds below `min_wind_speed` are clamped on a copy.
    """
    d = d_ratio * canopy_height
    zom = z0m_ratio * canopy_height
    _check_above_roughness(measurement_height, d, zom, "measurement_height")

    # numpy.maximum returns a new object (and keeps pandas Series as Series); NaN is propagated
    wind_speed = numpy.maximum(wind_speed, min_wind_speed)

    if wind_over_grass:
        _check_above_roughness(blending_height, d, zom, "blending_height")
        factor = grass_to_blending_wind_factor(grass_height, grass_measurement_height, blending_height)
        wind_speed = factor * wind_speed * numpy.log((measurement_height - d) / zom) / numpy.log(
            (blending_height - d) / zom)
    return wind_speed, d, zom


def aerodynamic_resistance_viney(wind_speed, air_temperature, canopy_temperature, canopy_height, measurement_height,
                                 d_ratio=2 / 3.0, z0m_ratio=0.13, z0h_ratio=0.13, wind_over_grass=True,
                                 grass_height=FAO_GRASS_HEIGHT, grass_measurement_height=FAO_GRASS_MEASUREMENT_HEIGHT,
                                 blending_height=BLENDING_HEIGHT, min_wind_speed=MIN_WIND_SPEED):
    """
    Calculates aerodynamic resistance based on Viney 1991

    Reference: Viney, Neil R. 1991. “An Empirical Expression for Aerodynamic Resistance in the Unstable Boundary Layer.”
        Boundary-Layer Meteorology 56 (4). Kluwer Academic Publishers: 381–93.

    Viney (1991) is an expression for the unstable boundary layer only. The stability correction is applied when the
    bulk Richardson number Ri < 0 (canopy warmer than the air). In neutral (Ri = 0) and stable (Ri > 0) conditions no
    correction is applied and the neutral aerodynamic resistance is returned. Where Ri is NaN (e.g. missing air or
    canopy temperature) the result is NaN.

    Parameters
    ----------
    d_ratio : float
        Zero-plane displacement (d) = d_ratio * canopy_height
    z0m_ratio : float
        Roughness length (z0) for momentum = z0m_ratio * canopy_height
    z0h_ratio : float
        Roughness length (z0) for heat transfer = z0h_ratio * canopy_height
    wind_speed : Union[ndarray, float]
        Wind speed [m/s]. See `wind_over_grass` for where it is assumed to be measured.
    air_temperature : Union[ndarray, float]
        Air temperature [Celsius]
    canopy_temperature : Union[ndarray, float]
        Canopy temperature, used to estimate buoyancy effects [Celsius]
    canopy_height : Union[ndarray, float]
        Canopy height [m]
    measurement_height : Union[ndarray, float]
        Reference height above the ground at which ra is computed [m]
    wind_over_grass : bool
        If True (default), `wind_speed` was measured at `grass_measurement_height` over reference grass of height
        `grass_height`. It is extrapolated to `blending_height` over the grass and then brought down to
        `measurement_height` over the crop with a neutral log profile. In this case `measurement_height` is the
        reference height above the crop, not the anemometer height.
        If False, `wind_speed` was measured at `measurement_height` over the crop and is used directly.
    grass_height : float
        Height of the reference grass [m]. Only used if wind_over_grass is True. Default 0.12 (FAO-56)
    grass_measurement_height : float
        Anemometer height over the grass [m]. Only used if wind_over_grass is True. Default 2.0
    blending_height : float
        Blending height [m]. Only used if wind_over_grass is True. Default 100.0
    min_wind_speed : float
        Wind speeds below this value are set to this value (on a copy) to avoid divisions by zero [m/s]

    Raises
    ------
    ValueError
        If measurement_height - d <= z0m, i.e. the reference height is inside or too close to the canopy.

    Returns
    -------
    Union[ndarray, float] : aerodynamic resistance [s/m]
    """

    zoh = z0h_ratio * canopy_height
    wind_speed_canopy, d, zom = _wind_at_reference_height(wind_speed, canopy_height, measurement_height, d_ratio,
                                                          z0m_ratio, wind_over_grass, grass_height,
                                                          grass_measurement_height, blending_height, min_wind_speed)
    logarithmic_profile_ym = numpy.log((measurement_height - d) / zom)

    a = 1.0591 - 0.0552 * numpy.log(1.72 + numpy.power(4.03 - logarithmic_profile_ym, 2))
    b = 1.9117 - 0.2237 * numpy.log(1.86 + numpy.power(2.12 - logarithmic_profile_ym, 2))
    c = 0.8437 - 0.1243 * numpy.log(3.49 + numpy.power(2.79 - logarithmic_profile_ym, 2))
    ra_prime = numpy.log((measurement_height - d) / zoh) * logarithmic_profile_ym / (
            VON_KARMAN * VON_KARMAN * wind_speed_canopy)
    richardson_number = (measurement_height - d) * 9.81 * (air_temperature - canopy_temperature) / (
            wind_speed_canopy * wind_speed_canopy * (air_temperature + 273.15))

    # Unstable (Ri < 0): divide by the Viney correction. Neutral or stable (Ri >= 0): neutral ra (correction = 1).
    # NaN Ri (missing temperatures): NaN.
    # -Ri is clipped at 0 so that the power is never evaluated on a negative base in the stable branch.
    richardson_number = numpy.asarray(richardson_number, dtype=float)
    correction = numpy.where(richardson_number < 0,
                             a + b * numpy.power(numpy.maximum(-richardson_number, 0.0), c),
                             numpy.where(numpy.isnan(richardson_number), numpy.nan, 1.0))
    return ra_prime / correction


def aerodynamic_resistance_fao(wind_speed, canopy_height, measurement_height,
                               d_ratio=2 / 3.0, z0m_ratio=0.123, z0h_ratio=0.0123, wind_over_grass=True,
                               grass_height=FAO_GRASS_HEIGHT, grass_measurement_height=FAO_GRASS_MEASUREMENT_HEIGHT,
                               blending_height=BLENDING_HEIGHT, min_wind_speed=MIN_WIND_SPEED):
    """
    Function to aerodynamic resistance based on FAO's Penman-Monteith (neutral conditions)
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
        Wind speed [m/s]. See `wind_over_grass` for where it is assumed to be measured.
    canopy_height : Union[ndarray, float]
        Canopy height [m]
    measurement_height : Union[ndarray, float]
        Reference height above the ground at which ra is computed [m]
    wind_over_grass : bool
        If True (default), `wind_speed` was measured at `grass_measurement_height` over reference grass of height
        `grass_height`. It is extrapolated to `blending_height` over the grass and then brought down to
        `measurement_height` over the crop with a neutral log profile. In this case `measurement_height` is the
        reference height above the crop, not the anemometer height.
        If False, `wind_speed` was measured at `measurement_height` over the crop and is used directly.
    grass_height : float
        Height of the reference grass [m]. Only used if wind_over_grass is True. Default 0.12 (FAO-56)
    grass_measurement_height : float
        Anemometer height over the grass [m]. Only used if wind_over_grass is True. Default 2.0
    blending_height : float
        Blending height [m]. Only used if wind_over_grass is True. Default 100.0
    min_wind_speed : float
        Wind speeds below this value are set to this value (on a copy) to avoid divisions by zero [m/s]

    Returns
    -------
    Union[ndarray, float] : aerodynamic resistance [s/m]

    Raises
    ------
    ValueError
        If measurement_height - d <= z0m, i.e. the reference height is inside or too close to the canopy.
    """

    zoh = z0h_ratio * canopy_height
    wind_speed_canopy, d, zom = _wind_at_reference_height(wind_speed, canopy_height, measurement_height, d_ratio,
                                                          z0m_ratio, wind_over_grass, grass_height,
                                                          grass_measurement_height, blending_height, min_wind_speed)
    ra_s = numpy.log((measurement_height - d) / zoh) * numpy.log((measurement_height - d) / zom) / (
            VON_KARMAN * VON_KARMAN * wind_speed_canopy)
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
    z0mh = (1 - d_h) * numpy.exp(-VON_KARMAN / gamma_prime - psi_h)

    raupach_z0 = canopy_height * z0mh
    raupach_d = canopy_height * d_h

    return raupach_z0, raupach_d


def adjust_wind_height(wind_speed, z, measurement_height=2.0, d_ratio=2 / 3.0, z0m_ratio=0.13, z0h_ratio=0.13,
                       wind_over_grass=True, grass_height=FAO_GRASS_HEIGHT,
                       grass_measurement_height=FAO_GRASS_MEASUREMENT_HEIGHT, blending_height=BLENDING_HEIGHT,
                       min_wind_speed=MIN_WIND_SPEED):
    """
    Estimates the wind speed at `measurement_height` above a surface (crop) of height `z`, from a wind speed measured
    over reference grass, assuming neutral logarithmic wind profiles.

    The wind measured at `grass_measurement_height` over grass is extrapolated to `blending_height` (where it is
    assumed to be independent of the surface) and then brought down to `measurement_height` over the crop:

        u(zm) = u_grass * F_grass * ln((zm - d) / z0m) / ln((zb - d) / z0m)

    where F_grass = grass_to_blending_wind_factor(...) (~1.81 for the FAO-56 defaults) and d and z0m are those of
    the crop.

    Parameters
    ----------
    wind_speed : Union[ndarray, float]
        Wind speed [m/s]. See `wind_over_grass` for where it is assumed to be measured.
    z : Union[ndarray, float]
        Height of the crop (surface) over which the wind speed is wanted [m]
    measurement_height : Union[ndarray, float]
        Height above the ground at which the wind speed is wanted [m]. Default 2.0
    d_ratio : float
        Zero-plane displacement of the crop (d) = d_ratio * z
    z0m_ratio : float
        Roughness length for momentum of the crop (z0m) = z0m_ratio * z
    z0h_ratio : float
        Not used. Kept for backwards compatibility.
    wind_over_grass : bool
        If True (default), `wind_speed` was measured at `grass_measurement_height` over reference grass of height
        `grass_height`. It is extrapolated to `blending_height` over the grass and then brought down to
        `measurement_height` over the crop with a neutral log profile. In this case `measurement_height` is the
        reference height above the crop, not the anemometer height.
        If False, `wind_speed` was measured at `measurement_height` over the crop and is used directly.
    grass_height : float
        Height of the reference grass [m]. Only used if wind_over_grass is True. Default 0.12 (FAO-56)
    grass_measurement_height : float
        Anemometer height over the grass [m]. Only used if wind_over_grass is True. Default 2.0
    blending_height : float
        Blending height [m]. Only used if wind_over_grass is True. Default 100.0
    min_wind_speed : float
        Wind speeds below this value are set to this value (on a copy) to avoid divisions by zero [m/s]

    Returns
    -------
    Union[ndarray, float] : wind speed at measurement_height over the crop [m/s]

    Raises
    ------
    ValueError
        If measurement_height - d <= z0m, i.e. the reference height is inside or too close to the canopy.
    """
    wind_speed_canopy, _, _ = _wind_at_reference_height(wind_speed, z, measurement_height, d_ratio, z0m_ratio,
                                                        wind_over_grass, grass_height, grass_measurement_height,
                                                        blending_height, min_wind_speed)
    return wind_speed_canopy
