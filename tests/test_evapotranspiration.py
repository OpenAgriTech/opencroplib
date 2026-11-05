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

from opencroplib.aerodynamic import wind_speed_2m
from opencroplib.atmophere import delta_calc, psychrometric_gamma, vpd_calc
from opencroplib.evapotranspiration import et0_FAO56, et0_FAO56_daily, et0_ASCE


class TestEvapotranspirationFunctions(TestCase):
    def test_wind_height_conversion(self):
        # ETo with wind measured at 10 m must equal ETo with the equivalent 2 m wind (it used to go the other way)
        u10 = numpy.array([1.0, 3.0, 6.0])
        args = (numpy.full(3, 30.0), 30.0, 101325.0)
        for et0 in (et0_FAO56, et0_ASCE):
            numpy.testing.assert_allclose(et0(*args, u10, 800.0, 500.0, z_wind=10.0),
                                          et0(*args, wind_speed_2m(u10, 10.0), 800.0, 500.0))
        numpy.testing.assert_allclose(et0_FAO56_daily(*args, u10, 25.0, 15.0, z_wind=10.0),
                                      et0_FAO56_daily(*args, wind_speed_2m(u10, 10.0), 25.0, 15.0))

    def test_fao56_example_19_hourly(self):
        # FAO-56 Example 19 (N'Diaye, Senegal, 8 m a.s.l.), default G = 0.1 Rn (day) / 0.5 Rn (night).
        # Wind coefficient was 9.4 instead of 37 and the default G was converted to MJ twice (G ~ 0).
        pressure = 101.3 * ((293.0 - 0.0065 * 8.0) / 293.0) ** 5.26 * 1000.0
        day = et0_FAO56(38.0, 52.0, pressure, 3.3, 500.0, 1.749 / 0.0036)
        night = et0_FAO56(28.0, 90.0, pressure, 1.9, 0.0, -0.100 / 0.0036)
        self.assertAlmostEqual(day, 0.63, delta=0.01)
        self.assertAlmostEqual(night, 0.0, delta=0.01)
        # Same result for arrays (scalar inputs used to get G = 0)
        numpy.testing.assert_allclose(et0_FAO56(numpy.array([38.0, 28.0]), numpy.array([52.0, 90.0]), pressure,
                                                numpy.array([3.3, 1.9]), numpy.array([500.0, 0.0]),
                                                numpy.array([1.749, -0.100]) / 0.0036), [day, night])

    def test_hourly_equation(self):
        # ASCE-EWRI (2005) Eq. 1 with Cn = 37 and Cd = 0.34 (FAO-56) or 0.24 / 0.96 (ASCE, day / night),
        # in kPa and MJ m-2 h-1, using opencroplib's own delta, gamma and VPD
        t = numpy.array([10.0, 25.0, 35.0, 30.0])
        rh = numpy.array([70.0, 40.0, 20.0, 80.0])
        u2 = numpy.array([1.0, 3.0, 6.0, 2.0])
        rn_w = numpy.array([300.0, 600.0, 700.0, -60.0])
        pressure = 101325.0
        rn = rn_w * 0.0036
        g = numpy.where(rn_w > 0, 0.1, 0.5) * rn
        delta = delta_calc(t) / 1000.0
        gamma = psychrometric_gamma(t, rh, pressure) / 1000.0
        vpd = vpd_calc(t, rh) / 1000.0

        def eq1(cd):
            return (0.408 * delta * (rn - g) + gamma * 37.0 / (t + 273.0) * u2 * vpd) / (delta + gamma * (1 + cd * u2))

        solar = numpy.where(rn_w > 0, 1.0, 0.0)
        numpy.testing.assert_allclose(et0_FAO56(t, rh, pressure, u2, solar, rn_w), eq1(0.34), rtol=1e-10)
        numpy.testing.assert_allclose(et0_ASCE(t, rh, pressure, u2, solar, rn_w),
                                      eq1(numpy.where(rn_w >= 0, 0.24, 0.96)), rtol=1e-10)
        # Hot and dry: more wind -> more ETo (with the 9.4 coefficient it decreased)
        self.assertTrue(numpy.all(numpy.diff(et0_FAO56(30.0, 30.0, pressure, numpy.array([0.5, 1, 2, 4, 8]), 1.0,
                                                       500.0)) > 0))
