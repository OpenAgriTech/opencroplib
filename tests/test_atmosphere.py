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

from opencroplib.atmophere import rho_cp_dry, es_calc, delta_calc, vpd_calc, latent_heat_lambda, specific_heat_cp, \
    psychrometric_gamma, air_density_rho, potential_temperature_theta


class TestAtmosphereFunctions(TestCase):
    def test_rho_cp_dry(self):
        # rho = 1.204 kg m-3 and cp ~ 1005 J kg-1 K-1 at 20 C and 101325 Pa (density used T instead of P)
        self.assertAlmostEqual(rho_cp_dry(20.0, 101325.0) / (1.204 * 1005.0), 1.0, delta=0.01)

    def test_es_calc_list(self):
        # The array branch used .shape, failing for plain lists
        numpy.testing.assert_allclose(es_calc([20, 25]), [es_calc(20.0), es_calc(25.0)])
        numpy.testing.assert_allclose(es_calc([20, 25]), [2337.08, 3166.82], rtol=1e-5)

    def test_list_inputs(self):
        # Docstring examples with plain lists (delta_calc etc. used list + float)
        t, rh, p = [10, 20, 30], [10, 20, 30], [100000, 101000, 102000]
        ta, rha, pa = numpy.array(t, float), numpy.array(rh, float), numpy.array(p, float)
        numpy.testing.assert_allclose(delta_calc([20, 25]), [144.665841, 188.625046], rtol=1e-6)
        numpy.testing.assert_allclose(vpd_calc([20, 25], [50, 100]), [1168.540099, 0.0], atol=1e-5)
        numpy.testing.assert_allclose(latent_heat_lambda(t), latent_heat_lambda(ta))
        numpy.testing.assert_allclose(specific_heat_cp(t, rh, p), specific_heat_cp(ta, rha, pa))
        numpy.testing.assert_allclose(psychrometric_gamma(t, rh, p), [65.255188, 66.656958, 68.242393], rtol=1e-6)
        numpy.testing.assert_allclose(air_density_rho(t, rh, p), [1.22948419, 1.19787662, 1.16635358], rtol=1e-6)
        numpy.testing.assert_allclose(potential_temperature_theta([5, 10, 20], [45, 65, 89], [101300, 102000, 99800]),
                                      [3.97741582, 8.40874555, 20.16596828], rtol=1e-6)
