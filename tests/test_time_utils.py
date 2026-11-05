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
        self.assertEqual(get_day_of_year(self.dt_utc), 195)
        self.assertEqual(get_day_of_year(self.dt_madrid), 195)
        self.assertEqual(get_day_of_year(self.dt_naive), 195)

    def test_solar_time(self):
        self.assertAlmostEqual(get_solar_time(self.long_utc, self.dt_utc), 11.908, 2)
        self.assertAlmostEqual(get_solar_time(self.long_utc, self.dt_madrid), 11.908, 2)
        self.assertAlmostEqual(get_solar_time(self.long_utc, self.dt_naive), 11.908, 2)
        self.assertAlmostEqual(get_solar_time(self.long_spain, self.dt_utc), 11.641, 2)
        self.assertAlmostEqual(get_solar_time(self.long_spain, self.dt_madrid), 11.641, 2)
        self.assertAlmostEqual(get_solar_time(self.long_spain, self.dt_naive), 11.641, 2)
