# -*- coding: utf-8 -*-
"""
Cross-validation of the ETo functions against established reference implementations:

- refet (https://github.com/WSWUP/RefET): ASCE-EWRI (2005) standardized reference ET, maintained by the
  Western States Water Use Program / Desert Research Institute.
- pyet (https://github.com/pyet-org/pyet): Vremec, Collenteur & Birk (2024), Geosci. Model Dev. 17, 7083-7103.

Two levels are checked:
1. The equation itself: refet's ASCE Eq. 1 fed with opencroplib's own delta, gamma and VPD must give the same ETo
   to machine precision.
2. The full calculation: with refet / pyet computing their own thermodynamics (Tetens es and delta,
   gamma = 0.000665 P) the results must agree within a tolerance. The only remaining difference is that opencroplib
   uses Bringfelt's psychrometric constant (humidity-dependent cp, temperature-dependent lambda), which differs from
   0.000665 P by about -2 % at 5 C and +3 % at 42 C.
"""
import numpy
import pandas
import pytest

from opencroplib.aerodynamic import wind_speed_2m
from opencroplib.atmophere import delta_calc, psychrometric_gamma, vpd_calc
from opencroplib.evapotranspiration import et0_FAO56, et0_FAO56_daily, et0_ASCE

refet_calcs = pytest.importorskip("refet.calcs")
pyet = pytest.importorskip("pyet")

W_TO_MJ_H = 0.0036
W_TO_MJ_D = 0.0864
PRESSURE = 101325.0


@pytest.fixture(scope="module")
def grid():
    t, rh, u2, rn = numpy.meshgrid([5.0, 15.0, 25.0, 35.0, 42.0], [10.0, 30.0, 50.0, 70.0, 95.0],
                                   [0.5, 1.0, 2.0, 4.0, 8.0], [-80.0, 50.0, 200.0, 400.0, 700.0], indexing="ij")
    t, rh, u2, rn = (x.ravel() for x in (t, rh, u2, rn))
    daytime = rn > 0
    rn_mj = rn * W_TO_MJ_H
    es = refet_calcs.sat_vapor_pressure(t)
    return dict(
        t=t, rh=rh, u2=u2, rn=rn, daytime=daytime, solar=numpy.where(daytime, 1.0, 0.0),
        rn_mj=rn_mj, g_mj=numpy.where(daytime, 0.1, 0.5) * rn_mj,
        # opencroplib thermodynamics [kPa]
        delta_o=delta_calc(t) / 1000.0, gamma_o=psychrometric_gamma(t, rh, PRESSURE) / 1000.0,
        vpd_o=vpd_calc(t, rh) / 1000.0,
        # refet thermodynamics [kPa]
        delta_r=refet_calcs.es_slope(t, method="refet"), gamma_r=0.000665 * PRESSURE / 1000.0,
        vpd_r=es * (1.0 - rh / 100.0),
    )


def _rel_diff(ours, ref, min_value):
    mask = numpy.abs(ref) > min_value
    return numpy.abs(ours[mask] - ref[mask]) / numpy.abs(ref[mask])


@pytest.mark.parametrize("method", ["fao56", "asce"])
def test_hourly_equation_matches_refet(grid, method):
    g = grid
    if method == "fao56":
        ours, cd = et0_FAO56(g["t"], g["rh"], PRESSURE, g["u2"], g["solar"], g["rn"]), 0.34
    else:
        ours, cd = et0_ASCE(g["t"], g["rh"], PRESSURE, g["u2"], g["solar"], g["rn"]), numpy.where(g["rn"] >= 0, 0.24,
                                                                                                    0.96)
    ref = refet_calcs.etsz(g["rn_mj"], g["g_mj"], g["t"], g["u2"], g["vpd_o"], g["delta_o"], g["gamma_o"], 37.0, cd)
    numpy.testing.assert_allclose(ours, ref, rtol=1e-10, atol=1e-12)


def test_daily_equation_matches_refet(grid):
    g = grid
    rn = g["rn"].clip(0) * W_TO_MJ_D
    ours = et0_FAO56_daily(g["t"], g["rh"], PRESSURE, g["u2"], 1.0, rn)
    ref = refet_calcs.etsz(rn, 0.0, g["t"], g["u2"], g["vpd_o"], g["delta_o"], g["gamma_o"], 900.0, 0.34)
    numpy.testing.assert_allclose(ours, ref, rtol=1e-10, atol=1e-12)


@pytest.mark.parametrize("method", ["fao56", "asce"])
def test_hourly_agrees_with_refet_and_pyet(grid, method):
    g = grid
    if method == "fao56":
        ours, cd = et0_FAO56(g["t"], g["rh"], PRESSURE, g["u2"], g["solar"], g["rn"]), 0.34
    else:
        ours, cd = et0_ASCE(g["t"], g["rh"], PRESSURE, g["u2"], g["solar"], g["rn"]), numpy.where(g["rn"] >= 0, 0.24,
                                                                                                    0.96)
    ref = refet_calcs.etsz(g["rn_mj"], g["g_mj"], g["t"], g["u2"], g["vpd_r"], g["delta_r"], g["gamma_r"], 37.0, cd)
    index = pandas.date_range("2020-07-01", periods=g["t"].size, freq="h")
    series = lambda x: pandas.Series(x, index=index)
    ref_pyet = pyet.pm_asce(series(g["t"]), series(g["u2"]), rn=series(g["rn_mj"]), g=series(g["g_mj"]),
                            rh=series(g["rh"]), pressure=PRESSURE / 1000.0, cn=37.0,
                            cd=series(numpy.broadcast_to(cd, g["t"].shape)), clip_zero=False).values
    numpy.testing.assert_allclose(ref_pyet, ref, rtol=1e-10)  # the two reference libraries agree

    rel = _rel_diff(ours, ref, 0.05)
    assert rel.mean() < 0.01
    assert rel.max() < 0.06
    day = g["daytime"]
    assert _rel_diff(ours[day], ref[day], 0.05).max() < 0.025


def test_daily_agrees_with_pyet(grid):
    g = grid
    rn = g["rn"].clip(0) * W_TO_MJ_D
    ours = et0_FAO56_daily(g["t"], g["rh"], PRESSURE, g["u2"], 1.0, rn)
    index = pandas.date_range("2020-01-01", periods=g["t"].size, freq="D")
    series = lambda x: pandas.Series(x, index=index)
    ref = pyet.pm_fao56(series(g["t"]), series(g["u2"]), rn=series(rn), g=0, rh=series(g["rh"]),
                        pressure=PRESSURE / 1000.0, clip_zero=False).values
    rel = _rel_diff(ours, ref, 0.1)
    assert rel.mean() < 0.01
    assert rel.max() < 0.05


def test_fao56_example_19_matches_refet_with_book_inputs():
    # FAO-56 Example 19 with the book's intermediate values (delta, gamma, es, ea, G)
    pressure = refet_calcs.air_pressure(8.0)[0] * 1000.0
    book_day = refet_calcs.etsz(1.749, 0.175, 38.0, 3.3, 6.625 - 3.445, 0.358, 0.0673, 37.0, 0.34)
    book_night = refet_calcs.etsz(-0.100, -0.050, 28.0, 1.9, 3.780 - 3.402, 0.220, 0.0673, 37.0, 0.34)
    assert et0_FAO56(38.0, 52.0, pressure, 3.3, 500.0, 1.749 / W_TO_MJ_H) == pytest.approx(book_day, abs=0.005)
    assert et0_FAO56(28.0, 90.0, pressure, 1.9, 0.0, -0.100 / W_TO_MJ_H) == pytest.approx(book_night, abs=0.002)


def test_wind_height_conversion_matches_refet():
    z = numpy.array([2.5, 3.0, 5.0, 10.0, 20.0])
    numpy.testing.assert_allclose(wind_speed_2m(1.0, z), refet_calcs.wind_height_adjust(1.0, z), rtol=5e-4)
