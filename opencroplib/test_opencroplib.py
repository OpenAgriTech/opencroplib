from unittest import TestCase
import datetime as dt

import pytz

from opencroplib.time_utils import get_day_of_year, get_solar_time
from opencroplib.radiation import planck, iplanck, potential_irradiance, longwave_downwelling
from opencroplib.radiation import SUPPORTED_ATMOSPHERIC_EMISSIVITY_MODELS

import numpy


# Test units for time.py functions
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
        # General test
        self.assertAlmostEqual(potential_irradiance(self.latitude_spain,
                                                    get_day_of_year(self.dt_utc),
                                                    get_solar_time(self.long_spain, self.dt_utc)),
                               1173.0, 1)
        # There should be more irradiance as we climb
        self.assertGreater(potential_irradiance(self.latitude_spain,
                                                get_day_of_year(self.dt_utc),
                                                get_solar_time(self.long_spain, self.dt_utc),
                                                1000),
                           potential_irradiance(self.latitude_spain,
                                                get_day_of_year(self.dt_utc),
                                                get_solar_time(self.long_spain, self.dt_utc),
                                                0))
        # Test that it works for a range of hours
        self.assertGreater(potential_irradiance(self.latitude_spain,
                                                get_day_of_year(self.dt_utc),
                                                numpy.arange(0, 24, 1)).sum(),
                           10000)

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
        self.assertTrue('model not supported' in context.exception.message)

        # Test the available models
        for model in SUPPORTED_ATMOSPHERIC_EMISSIVITY_MODELS:
            lw = longwave_downwelling(self.air_temperature, self.relative_humidity, emissivity_model=model)
            self.assertTrue(300 <= lw <= 500, msg="{model}: {lw}".format(model=model, lw=lw))

if __name__ == '__main__':
    unittest.main()
