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
from unittest import TestCase

import numpy
import pandas

from opencroplib.aerodynamic import aerodynamic_resistance_viney, aerodynamic_resistance_fao, adjust_wind_height, \
    grass_to_blending_wind_factor, wind_speed_2m, VON_KARMAN, MIN_WIND_SPEED
from opencroplib.plant import sensible_heat_flux_2


class TestAerodynamicFunctions(TestCase):
    def setUp(self):
        self.air_temperature = 30.0
        self.canopy_height = 0.5
        self.measurement_height = 2.0

    def test_inputs_not_modified(self):
        wind = numpy.array([0.0, 0.05, 2.0])
        wind_series = pandas.Series(wind.copy())
        tc = numpy.array([32.0, 32.0, 32.0])
        for w in (wind, wind_series):
            original = numpy.array(w, copy=True)
            aerodynamic_resistance_viney(w, self.air_temperature, tc, self.canopy_height, self.measurement_height)
            aerodynamic_resistance_fao(w, self.canopy_height, self.measurement_height)
            adjust_wind_height(w, self.canopy_height, self.measurement_height)
            numpy.testing.assert_array_equal(numpy.asarray(w), original)

        ra = numpy.array([0.0, 20.0])
        sensible_heat_flux_2(30.0, 40.0, 101325.0, 32.0, ra)
        numpy.testing.assert_array_equal(ra, [0.0, 20.0])

    def test_scalars_and_small_wind(self):
        # Scalars work, and small speeds (not only exact zeros) are clamped to the minimum
        ra_zero = aerodynamic_resistance_fao(0.0, self.canopy_height, self.measurement_height)
        ra_small = aerodynamic_resistance_fao(0.05, self.canopy_height, self.measurement_height)
        ra_min = aerodynamic_resistance_fao(MIN_WIND_SPEED, self.canopy_height, self.measurement_height)
        self.assertTrue(numpy.isfinite(ra_zero))
        self.assertEqual(ra_zero, ra_min)
        self.assertEqual(ra_small, ra_min)
        self.assertTrue(numpy.isfinite(aerodynamic_resistance_viney(0.0, 30.0, 32.0, 0.5, 2.0)))
        self.assertTrue(numpy.isfinite(adjust_wind_height(0.0, 0.5, 2.0)))
        # sensible_heat_flux_2 used item assignment, failing for a scalar resistance
        self.assertAlmostEqual(sensible_heat_flux_2(30.0, 40.0, 101325.0, 32.0, 0.0),
                               sensible_heat_flux_2(30.0, 40.0, 101325.0, 32.0, 0.01))

    def test_grass_wind_factor(self):
        self.assertAlmostEqual(grass_to_blending_wind_factor() / 1.82, 1.0, delta=0.01)

    def test_reference_height_inside_canopy_raises(self):
        # Anemometer at 3 m in a 4 m orchard: d = 2.67 m, z0m = 0.49 m -> z - d <= z0m
        with self.assertRaises(ValueError):
            aerodynamic_resistance_fao(2.0, 4.0, 3.0)
        with self.assertRaises(ValueError):
            aerodynamic_resistance_viney(2.0, 30.0, 32.0, 4.0, 3.0)
        with self.assertRaises(ValueError):
            adjust_wind_height(2.0, 4.0, 3.0)
        with self.assertRaises(ValueError):
            aerodynamic_resistance_fao(2.0, 4.0, 3.0, wind_over_grass=False)
        with self.assertRaises(ValueError):
            aerodynamic_resistance_fao(numpy.array([2.0, 2.0]), numpy.array([1.0, 4.0]), 3.0)

    def test_wind_over_grass_false(self):
        u = numpy.array([1.0, 3.0])
        numpy.testing.assert_allclose(adjust_wind_height(u, 4.0, 6.0, wind_over_grass=False), u)
        d, z0m, z0h = 2 / 3.0 * 4.0, 0.123 * 4.0, 0.0123 * 4.0
        expected = numpy.log((6.0 - d) / z0h) * numpy.log((6.0 - d) / z0m) / (VON_KARMAN ** 2 * u)
        numpy.testing.assert_allclose(aerodynamic_resistance_fao(u, 4.0, 6.0, wind_over_grass=False), expected)

    def test_viney(self):
        wind = numpy.array([1.0, 2.0, 4.0])
        # Scalar heights with an array of wind speeds
        neutral = aerodynamic_resistance_fao(wind, self.canopy_height, self.measurement_height,
                                             z0m_ratio=0.13, z0h_ratio=0.13)
        unstable = aerodynamic_resistance_viney(wind, self.air_temperature, self.air_temperature + 5.0,
                                                self.canopy_height, self.measurement_height)
        stable = aerodynamic_resistance_viney(wind, self.air_temperature, self.air_temperature - 5.0,
                                              self.canopy_height, self.measurement_height)
        self.assertEqual(unstable.shape, wind.shape)
        self.assertTrue(numpy.all(unstable < neutral))
        # Same von Karman constant in both functions -> identical neutral ra
        numpy.testing.assert_allclose(stable, neutral)
        # Mixed stable / unstable array
        mixed = aerodynamic_resistance_viney(wind, self.air_temperature, numpy.array([35.0, 25.0, 35.0]),
                                             self.canopy_height, self.measurement_height)
        numpy.testing.assert_allclose(mixed, [unstable[0], stable[1], unstable[2]])
        # Scalar in, scalar out
        self.assertTrue(numpy.ndim(aerodynamic_resistance_viney(2.0, 30.0, 35.0, 0.5, 2.0)) == 0)

    def test_viney_neutral_and_missing(self):
        wind = numpy.array([1.0, 2.0, 4.0])
        neutral = aerodynamic_resistance_fao(wind, self.canopy_height, self.measurement_height,
                                             z0m_ratio=0.13, z0h_ratio=0.13)
        # Ri = 0 (Tc == Ta) is neutral: no Viney correction (it used to return ra / a)
        numpy.testing.assert_allclose(
            aerodynamic_resistance_viney(wind, 30.0, 30.0, self.canopy_height, self.measurement_height), neutral)
        # Missing canopy temperature -> NaN, not the neutral ra
        ra = aerodynamic_resistance_viney(wind, 30.0, numpy.array([numpy.nan, 35.0, 25.0]),
                                          self.canopy_height, self.measurement_height)
        self.assertTrue(numpy.isnan(ra[0]))
        self.assertTrue(numpy.all(numpy.isfinite(ra[1:])))
        self.assertTrue(numpy.isnan(aerodynamic_resistance_viney(2.0, 30.0, numpy.nan, 0.5, 2.0)))

    def test_wind_speed_2m(self):
        # FAO-56 Eq. 47: u2 = uz * 4.87 / ln(67.8 z - 5.42)
        for z in (3.0, 5.0, 10.0):
            self.assertAlmostEqual(wind_speed_2m(3.0, z) / (3.0 * 4.87 / numpy.log(67.8 * z - 5.42)), 1.0, delta=1e-3)
        self.assertAlmostEqual(wind_speed_2m(3.0, 2.0), 3.0)
        # Wind measured higher up is faster, so it must be reduced
        self.assertLess(wind_speed_2m(3.0, 10.0), 3.0)
        with self.assertRaises(ValueError):
            wind_speed_2m(3.0, 0.05)
