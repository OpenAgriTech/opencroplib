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

from opencroplib.aerodynamic import aerodynamic_resistance_viney
from opencroplib.plant import canopy_conductance, canopy_conductance_2, canopy_temperature_sim, sensible_heat_flux
from opencroplib.radiation import available_energy, longwave_downwelling


class TestPlantFunctions(TestCase):
    def setUp(self):
        self.shortwave = 800.0
        self.albedo = 0.23
        self.emissivity = 0.98
        self.air_pressure = 101325.0

    def _round_trip(self, gc, air_temperature, relative_humidity, aerodynamic_resistance):
        tc = canopy_temperature_sim(self.shortwave, air_temperature, relative_humidity, self.air_pressure, gc,
                                    aerodynamic_resistance, albedo=self.albedo, emissivity=self.emissivity)
        longwave = longwave_downwelling(air_temperature, relative_humidity)
        net_radiation = available_energy(self.shortwave, self.albedo, longwave, tc, self.emissivity)
        ra = aerodynamic_resistance(tc) if callable(aerodynamic_resistance) else aerodynamic_resistance
        gc_back = canopy_conductance(tc, air_temperature, relative_humidity, self.air_pressure, net_radiation, ra)
        return tc, gc_back

    def test_round_trip(self):
        # Simulate Tc for known gc, invert and recover gc within 0.1 %
        gc = numpy.array([0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0])
        for air_temperature, relative_humidity, ra in ((15.0, 70.0, 40.0), (30.0, 30.0, 20.0), (38.0, 15.0, 10.0)):
            tc, gc_back = self._round_trip(gc, air_temperature, relative_humidity, ra)
            numpy.testing.assert_allclose(gc_back, gc, rtol=1e-3)
            # More conductance -> cooler canopy
            self.assertTrue(numpy.all(numpy.diff(tc) < 0))
        # canopy_conductance_2 is the same inversion
        tc, _ = self._round_trip(gc, 30.0, 30.0, 20.0)
        rn = available_energy(self.shortwave, self.albedo, longwave_downwelling(30.0, 30.0), tc, self.emissivity)
        numpy.testing.assert_allclose(canopy_conductance_2(tc, 30.0, 30.0, self.air_pressure, rn, 20.0), gc,
                                      rtol=1e-3)

    def test_round_trip_callable_ra(self):
        # Stability-dependent ra(Tc) for a 4 m orchard, reference height 6 m
        def ra(tc):
            return aerodynamic_resistance_viney(2.0, 32.0, tc, 4.0, 6.0)

        gc = numpy.array([1.0, 5.0, 15.0])
        _, gc_back = self._round_trip(gc, 32.0, 25.0, ra)
        numpy.testing.assert_allclose(gc_back, gc, rtol=1e-3)

    def test_non_transpiring_limit(self):
        ra = 20.0
        tc = canopy_temperature_sim(self.shortwave, 30.0, 30.0, self.air_pressure, numpy.array([0.0, 0.5]), ra)
        self.assertGreater(tc[0], tc[1])
        # gc = 0 -> LE = 0 -> Rn = H
        rn = available_energy(self.shortwave, self.albedo, longwave_downwelling(30.0, 30.0), tc[0], self.emissivity)
        self.assertAlmostEqual(rn, sensible_heat_flux(30.0, 30.0, self.air_pressure, tc[0], ra), places=4)

    def test_canopy_conductance_invalid(self):
        # LE = Rn - H <= 0 -> NaN (Rn = 50 W m-2, H ~ 600 W m-2)
        self.assertTrue(numpy.isnan(canopy_conductance(40.0, 30.0, 30.0, self.air_pressure, 50.0, 20.0)))
        # e*(Tc) < ea (condensation) -> rc < 0 -> NaN
        self.assertTrue(numpy.isnan(canopy_conductance(15.0, 30.0, 80.0, self.air_pressure, 500.0, 20.0)))
        # Vectorised, never negative or infinite
        gc = canopy_conductance(numpy.array([40.0, 15.0, 31.0]), 30.0, numpy.array([30.0, 80.0, 30.0]),
                                self.air_pressure, numpy.array([50.0, 500.0, 500.0]), 20.0)
        self.assertTrue(numpy.isnan(gc[0]) and numpy.isnan(gc[1]))
        self.assertTrue(numpy.isfinite(gc[2]) and gc[2] > 0)

    def test_pandas_series(self):
        index = pandas.date_range("2009-07-15 12:00", periods=3, freq="h")
        tc = pandas.Series([31.0, 32.0, 33.0], index=index)
        gc = canopy_conductance(tc, 30.0, 30.0, self.air_pressure, 500.0, 20.0)
        self.assertIsInstance(gc, pandas.Series)
        self.assertTrue(gc.index.equals(index))
        tc_sim = canopy_temperature_sim(self.shortwave, pandas.Series([30.0, 31.0, 32.0], index=index), 30.0,
                                        self.air_pressure, 5.0, 20.0)
        self.assertIsInstance(tc_sim, pandas.Series)
