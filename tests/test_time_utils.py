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
from unittest import TestCase

import numpy
import pandas
import pytest
import pytz

from opencroplib.time_utils import get_day_of_year, get_solar_time


class TestTimeFunctions(TestCase):
    def setUp(self):
        self.dt_utc = dt.datetime(2017, 7, 15, 12, 0, tzinfo=pytz.utc)
        self.dt_madrid = pytz.timezone("Europe/Madrid").localize(dt.datetime(2017, 7, 15, 14, 0))
        self.dt_naive = dt.datetime(2017, 7, 15, 12, 0)
        self.long_utc = 0.0
        self.long_spain = -4.0

    def test_get_day_of_year(self):
        # 15 July 2017 is day 196 (1 January = 1). opencroplib < 0.2.0 returned 195 for a single datetime.
        self.assertEqual(get_day_of_year(self.dt_utc), 196)
        self.assertEqual(get_day_of_year(self.dt_madrid), 196)
        self.assertEqual(get_day_of_year(self.dt_naive), 196)
        self.assertEqual(get_day_of_year(dt.datetime(2024, 1, 1)), 1)
        self.assertEqual(get_day_of_year(dt.datetime(2024, 12, 31)), 366)

    def test_solar_time(self):
        self.assertAlmostEqual(get_solar_time(self.long_utc, self.dt_utc), 11.908, 2)
        self.assertAlmostEqual(get_solar_time(self.long_utc, self.dt_madrid), 11.908, 2)
        self.assertAlmostEqual(get_solar_time(self.long_utc, self.dt_naive), 11.908, 2)
        self.assertAlmostEqual(get_solar_time(self.long_spain, self.dt_utc), 11.641, 2)
        self.assertAlmostEqual(get_solar_time(self.long_spain, self.dt_madrid), 11.641, 2)
        self.assertAlmostEqual(get_solar_time(self.long_spain, self.dt_naive), 11.641, 2)


@pytest.mark.parametrize("tz", [None, "UTC", "Europe/Madrid", "America/Los_Angeles", "Australia/Sydney"])
def test_scalar_and_index_agree(tz):
    # A single datetime and a DatetimeIndex must give the same day of year and solar time
    index = pandas.date_range("2019-12-30 00:30", "2021-01-02 23:30", freq="7h", tz=tz)
    scalars = [t.to_pydatetime() for t in index]
    numpy.testing.assert_array_equal([get_day_of_year(t) for t in scalars], numpy.asarray(get_day_of_year(index)))
    for longitude in (-120.0, -4.8, 0.0, 150.0):
        solar_h = numpy.array([get_solar_time(longitude, t) for t in scalars])
        numpy.testing.assert_allclose(solar_h, numpy.asarray(get_solar_time(longitude, index), float), atol=1e-12)
        # Always in [0, 24): the scalar path used to return negative hours west of Greenwich
        assert numpy.all((solar_h >= 0.0) & (solar_h < 24.0))


def test_equation_of_time_against_spencer():
    pvlib = pytest.importorskip("pvlib")
    from opencroplib.time_utils import equation_of_time
    doy = numpy.arange(1, 366, dtype=float)
    assert numpy.max(numpy.abs(equation_of_time(doy) - pvlib.solarposition.equation_of_time_spencer71(doy))) < 1.2
