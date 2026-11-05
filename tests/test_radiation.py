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
import datetime as dt
import math
from unittest import TestCase

import numpy
import pytest
import pandas
import pytz

from opencroplib.radiation import SUPPORTED_ATMOSPHERIC_EMISSIVITY_MODELS
from opencroplib.radiation import planck, iplanck, potential_irradiance, longwave_downwelling, solar_elevation, \
    toa_irradiance, clear_sky_beam, clear_sky_diffuse, solar_declination, solar_azimuth, daily_potential_toa, \
    sunset_hour_angle, solar_declination_fao, net_out_lw_daily, net_out_lw_hourly
from opencroplib.time_utils import get_day_of_year, get_solar_time


class TestRadiationFunctions(TestCase):
    def setUp(self):
        self.dt_utc = dt.datetime(2017, 7, 15, 12, 0, tzinfo=pytz.utc)
        self.dt_madrid = pytz.timezone("Europe/Madrid").localize(dt.datetime(2017, 7, 15, 14, 0))
        self.dt_naive = dt.datetime(2017, 7, 15, 12, 0)
        self.long_utc = 0.0
        self.long_spain = -4.0
        self.latitude_spain = 37.0
        self.canopy_temperature = 29.0
        self.air_temperature = 30.0
        self.relative_humidity = 20.0

    def test_planck_functions(self):
        # Make sure direct and inverse functions return the same
        self.assertAlmostEqual(self.canopy_temperature + 273.15, iplanck(planck(self.canopy_temperature + 273.15)), 2)
        self.assertAlmostEqual(10, planck(iplanck(10)), 2)

    def test_potential_irradiance(self):
        # General test: ESRA global (beam + diffuse) clear-sky irradiance, TL = 1. Reference value from GRASS r.sun
        # (rsun_beam + rsun_diffuse below) at the same solar altitude. The previous implementation returned 1176 W/m2
        # (beam only, miscoded ESRA Rayleigh thickness and refraction, and a declination 0.7 deg too large).
        doy = get_day_of_year(self.dt_utc)
        solar_h = get_solar_time(self.long_spain, self.dt_utc)
        alt = solar_elevation(self.latitude_spain, doy, solar_h)
        self.assertAlmostEqual(potential_irradiance(self.latitude_spain, doy, solar_h),
                               rsun_beam(doy, alt, 0.0, 1.0) + rsun_diffuse(doy, alt, 1.0), delta=1e-9)
        self.assertAlmostEqual(potential_irradiance(self.latitude_spain, doy, solar_h), 1161.4, delta=0.1)
        # There should be more irradiance as we climb
        self.assertGreater(potential_irradiance(self.latitude_spain,
                                                get_day_of_year(self.dt_utc),
                                                get_solar_time(self.long_spain, self.dt_utc),
                                                1000),
                           potential_irradiance(self.latitude_spain,
                                                get_day_of_year(self.dt_utc),
                                                get_solar_time(self.long_spain, self.dt_utc),
                                                0))
        # Test that it works for a range of hours, against the r.sun reference
        doy = get_day_of_year(self.dt_utc)
        hours = numpy.arange(0, 24, 1)
        self.assertAlmostEqual(potential_irradiance(self.latitude_spain, doy, hours).sum(),
                               sum(rsun_beam(doy, a, 0.0, 1.0) + rsun_diffuse(doy, a, 1.0)
                                   for a in solar_elevation(self.latitude_spain, doy, hours)),
                               delta=1e-6)

        # Earlier in the summer there should be more radiation
        self.assertLess(potential_irradiance(self.latitude_spain,
                                                get_day_of_year(self.dt_utc),
                                                numpy.arange(0, 24, 1)).sum(),
                           potential_irradiance(self.latitude_spain,
                                                get_day_of_year(self.dt_utc)-2,
                                                numpy.arange(0, 24, 1)).sum())

        # Aerosols should reduce radiation
        self.assertGreater(potential_irradiance(self.latitude_spain,
                                                get_day_of_year(self.dt_utc),
                                                numpy.arange(0, 24, 1)).sum(),
                           potential_irradiance(self.latitude_spain,
                                                get_day_of_year(self.dt_utc),
                                                numpy.arange(0, 24, 1),
                                                t_lk=2).sum())

    def test_longwave_downwelling(self):
        # Test supported models
        with self.assertRaises(ValueError) as context:
            longwave_downwelling(self.air_temperature, self.relative_humidity,
                                                           emissivity_model='foo')
        self.assertTrue('model not supported' in str(context.exception))

        # Test the available models
        for model in SUPPORTED_ATMOSPHERIC_EMISSIVITY_MODELS:
            lw = longwave_downwelling(self.air_temperature, self.relative_humidity, emissivity_model=model)
            self.assertTrue(300 <= lw <= 500, msg="{model}: {lw}".format(model=model, lw=lw))


def rsun_beam(doy, alt, z, linke):
    """
    Independent transcription of GRASS GIS r.sun (raster/r.sun/rsunlib.c: com_sol_const and brad) for a horizontal
    surface, kept in r.sun's own form: ESRA clear-sky beam irradiance [W/m2] for solar altitude `alt` [rad], site
    elevation `z` [m] and Linke turbidity `linke`. r.sun only computes it when the solar altitude is > 0.
    """
    if alt <= 0.0:
        return 0.0
    i0 = 1367.0 * (1 + 0.03344 * math.cos(2 * math.pi * doy / 365.25 - 0.048869))
    elevation_corr = math.exp(-z / 8434.5)
    temp1 = 0.1594 + alt * (1.123 + 0.065656 * alt)
    temp2 = 1. + alt * (28.9344 + 277.3971 * alt)
    h0refract = alt + 0.061359 * temp1 / temp2  # in radians
    m = elevation_corr / (math.sin(h0refract) + 0.50572 * (math.degrees(h0refract) + 6.07995) ** -1.6364)
    if m <= 20.:
        rayl = 1. / (6.6296 + m * (1.7513 + m * (-0.1202 + m * (0.0065 - m * 0.00013))))
    else:
        rayl = 1. / (10.4 + 0.718 * m)
    return i0 * math.sin(alt) * math.exp(-rayl * m * 0.8662 * linke)


def rsun_diffuse(doy, alt, linke):
    """
    Independent transcription of GRASS GIS r.sun (raster/r.sun/rsunlib.c: com_sol_const and drad) for a horizontal
    surface: ESRA clear-sky diffuse irradiance [W/m2]. r.sun only computes it when the solar altitude is > 0.
    """
    if alt <= 0.0:
        return 0.0
    i0 = 1367.0 * (1 + 0.03344 * math.cos(2 * math.pi * doy / 365.25 - 0.048869))
    tn = -0.015843 + linke * (0.030543 + 0.0003797 * linke)
    a1b = 0.26463 + linke * (-0.061581 + 0.0031408 * linke)
    a1 = 0.0022 / tn if a1b * tn < 0.0022 else a1b
    a2 = 2.04020 + linke * (0.018945 - 0.011161 * linke)
    a3 = -1.3025 + linke * (0.039231 + 0.0085079 * linke)
    sin_alt = math.sin(alt)
    return i0 * (a1 + a2 * sin_alt + a3 * sin_alt * sin_alt) * tn


@pytest.mark.parametrize("latitude", [-35.0, 0.0, 37.0, 66.0])
@pytest.mark.parametrize("doy", [1, 172, 266])
@pytest.mark.parametrize("altitude, t_lk", [(0.0, 1.0), (0.0, 3.0), (1000.0, 3.0), (3000.0, 6.0), (0.0, 8.0)])
def test_clear_sky_matches_rsun(latitude, doy, altitude, t_lk):
    hours = numpy.arange(0.0, 24.0, 0.1)
    alts = solar_elevation(latitude, doy, hours)
    beam = [rsun_beam(doy, a, altitude, t_lk) for a in alts]
    diffuse = [rsun_diffuse(doy, a, t_lk) for a in alts]
    numpy.testing.assert_allclose(clear_sky_beam(latitude, doy, hours, altitude, t_lk), beam, rtol=1e-12, atol=1e-9)
    numpy.testing.assert_allclose(clear_sky_diffuse(latitude, doy, hours, t_lk), diffuse, rtol=1e-12, atol=1e-9)
    numpy.testing.assert_allclose(potential_irradiance(latitude, doy, hours, altitude, t_lk),
                                  numpy.add(beam, diffuse), rtol=1e-12, atol=1e-9)


def test_potential_irradiance_low_sun_and_scalars():
    # Air mass > 20 (sun ~3 deg): the m > 20 branch of the Rayleigh thickness. The old code had the branches swapped.
    hours = numpy.linspace(4.0, 7.0, 3001)
    alts = solar_elevation(37.0, 196, hours)
    for a_deg in (1.0, 2.0, 3.0, 4.0):
        h = hours[numpy.argmin(numpy.abs(numpy.degrees(alts) - a_deg))]
        assert clear_sky_beam(37.0, 196, h, 0, 3) == pytest.approx(
            rsun_beam(196, solar_elevation(37.0, 196, h), 0, 3), rel=1e-12)
    # Night: exactly 0, scalar in -> float out, no warnings from negative bases
    with numpy.errstate(all="raise"):
        assert potential_irradiance(37.0, 196, 0.0) == 0.0
        assert isinstance(potential_irradiance(37.0, 196, 12.0), float)
        assert numpy.all(potential_irradiance(37.0, 196, numpy.arange(0.0, 24.0, 0.5)) >= 0.0)


def test_toa_irradiance():
    # r.sun com_sol_const: 1367 * (1 + 0.03344 cos(2 pi doy / 365.25 - 0.048869)); perihelion early January
    for doy in (1, 100, 196, 300):
        assert toa_irradiance(doy) == pytest.approx(
            1367.0 * (1 + 0.03344 * math.cos(2 * math.pi * doy / 365.25 - 0.048869)), rel=1e-12)
    assert toa_irradiance(3) > toa_irradiance(185)


@pytest.fixture(scope="module")
def spa_reference():
    """Solar position from the NREL SPA (pvlib) and opencroplib, 2017-2020 hourly, five latitudes."""
    pvlib = pytest.importorskip("pvlib")
    times = pandas.date_range("2017-01-01", "2020-12-31 23:00", freq="h", tz="UTC")
    doy = numpy.asarray(times.day_of_year, float)
    rows = []
    for lat, lon in ((-35.0, 150.0), (0.0, -60.0), (37.4, -5.9), (52.0, 5.0), (66.0, 25.0)):
        spa = pvlib.solarposition.spa_python(times, lat, lon)
        solar_h = numpy.asarray(get_solar_time(lon, times), float)
        rows.append(dict(elevation=numpy.degrees(solar_elevation(lat, doy, solar_h)),
                         azimuth_north=(numpy.degrees(solar_azimuth(lat, doy, solar_h)) + 180.0) % 360.0,
                         spa_elevation=spa["elevation"].values, spa_azimuth=spa["azimuth"].values))
    return {k: numpy.concatenate([r[k] for r in rows]) for k in rows[0]}


def test_solar_declination_spencer():
    pvlib = pytest.importorskip("pvlib")
    doy = numpy.arange(1, 367, dtype=float)
    numpy.testing.assert_allclose(solar_declination(doy), pvlib.solarposition.declination_spencer71(doy), atol=1e-12)
    # The old formula reached 24.1 deg; the real maximum is the obliquity, 23.44 deg
    assert numpy.degrees(solar_declination(doy)).max() == pytest.approx(23.44, abs=0.05)


def test_solar_elevation_against_spa(spa_reference):
    ref = spa_reference
    up = ref["spa_elevation"] > 5
    err = (ref["elevation"] - ref["spa_elevation"])[up]
    assert numpy.sqrt(numpy.mean(err ** 2)) < 0.25
    assert numpy.max(numpy.abs(err)) < 0.65


def test_solar_azimuth_against_spa(spa_reference):
    ref = spa_reference
    # Azimuth is ill-conditioned near the zenith (a small elevation error becomes a large azimuth error)
    sel = (ref["spa_elevation"] > 5) & (ref["spa_elevation"] < 70)
    err = ((ref["azimuth_north"] - ref["spa_azimuth"] + 180.0) % 360.0 - 180.0)[sel]
    assert numpy.sqrt(numpy.mean(err ** 2)) < 0.35
    assert numpy.max(numpy.abs(err)) < 2.0


def test_solar_azimuth_sign_convention():
    # Duffie and Beckman: from south, east (morning) negative, west (afternoon) positive
    assert solar_azimuth(37.0, 196, 9.0) < 0.0 < solar_azimuth(37.0, 196, 15.0)
    assert solar_azimuth(37.0, 196, 12.0) == pytest.approx(0.0, abs=1e-12)
    assert solar_azimuth(37.0, 196, 9.0) == pytest.approx(-solar_azimuth(37.0, 196, 15.0))
    # In the tropics in July the noon sun is to the north
    assert abs(solar_azimuth(5.0, 196, 12.0)) == pytest.approx(math.pi)


def test_clear_sky_against_ineichen():
    # Plausibility against a different clear-sky model with the same Linke turbidity scale (Ineichen and Perez 2002).
    # The two models are known to differ by a few percent; this only guards against gross errors.
    pvlib = pytest.importorskip("pvlib")
    location = pvlib.location.Location(37.4, -5.9, altitude=0)
    for month, t_lk in ((1, 2.0), (4, 3.0), (7, 3.0), (7, 4.5), (10, 3.0)):
        times = pandas.date_range("2018-{:02d}-15".format(month), periods=24 * 12, freq="5min", tz="UTC")
        ineichen = location.get_clearsky(times, model="ineichen", linke_turbidity=t_lk)
        doy = numpy.asarray(times.day_of_year, float)
        solar_h = numpy.asarray(get_solar_time(-5.9, times), float)
        ours = potential_irradiance(37.4, doy, solar_h, 0.0, t_lk)
        assert ours.sum() / ineichen["ghi"].sum() == pytest.approx(1.0, abs=0.10)


def test_fao56_example_8_extraterrestrial():
    # FAO-56 Example 8: 20 S, 3 September (J = 246): ws = 1.527 rad, Ra = 32.2 MJ m-2 day-1
    assert sunset_hour_angle(-20.0, solar_declination_fao(246)) == pytest.approx(1.527, abs=5e-4)
    assert daily_potential_toa(-20.0, 246) == pytest.approx(32.2, abs=0.05)


def test_fao56_net_longwave():
    # FAO-56 Example 10 (daily): Tmax 25.1, Tmin 19.1 C, ea 2.1 kPa, Rs/Rso = 14.5/18.8 -> Rnl = 3.5 MJ m-2 day-1
    assert net_out_lw_daily(19.1, 25.1, 14.5, 18.8, 2.1) == pytest.approx(3.5, abs=0.05)
    # FAO-56 Example 19 (hourly, 14-15 h): T 38 C, RH 52 % (ea 3.445 kPa), Rs 2.450, Rso 2.658 MJ m-2 h-1.
    # Rnl = Rns - Rn = 0.77 * 2.450 - 1.749 = 0.138 MJ m-2 h-1. It used ea in Pa and a sigma 4x too small.
    assert net_out_lw_hourly(38.0, 52.0, 2.450 / 2.658) == pytest.approx(0.77 * 2.450 - 1.749, abs=0.003)
