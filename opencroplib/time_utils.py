# -*- coding: utf-8 -*-
#
#    Copyright 2017 Jose A. Jimenez-Berni
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.
#
import datetime as dt
import numpy
import pytz
import pandas as pd

__author__ = "Jose A. Jimenez-Berni  <berni@ias.csic.es>"
__version__ = "0.1"
__date__ = "July 2017"


def get_day_of_year(datetime):
    """
    Calculate the day of the year (DOY) for a given datetime: 1 for 1 January, as expected by the radiation
    functions (FAO-56 J, Spencer 1971). The calendar date of the datetime is used: the local date for timezone-aware
    datetimes, the date itself for naive ones.

    opencroplib < 0.2.0 returned the number of days *after* 1 January (0 for 1 January) for a single datetime, one
    day less than for a pd.DatetimeIndex.

    Parameters
    ----------
    datetime : Union[pd.DatetimeIndex, dt.datetime, pd.Timestamp]
        datetime representing the date

    Returns
    -------
    Union[pd.Index, int] : day of the year [1-366]

    """
    if isinstance(datetime, pd.DatetimeIndex):
        return datetime.day_of_year
    return datetime.timetuple().tm_yday


def get_solar_time(longitude_deg, datetime):
    """
    Calculates the local apparent solar time (decimal hour, 0-24) for a given location and datetime:

        solar_h = UTC time + 4 min * longitude + equation of time

    The equation of time is evaluated for the calendar date of `datetime` (see get_day_of_year). Pass timezone-aware
    local times so that this date, and the day of the year used with the result, is the local day.

    Parameters
    ----------
    longitude_deg : Union[ndarray, float]
        Geographic longitude [degrees]
    datetime : Union[pd.DatetimeIndex, dt.datetime]
        Datetime with the date and time. If the datetime has not timezone it is assumed that is UTC

    Returns
    -------
    Union[ndarray, float]: decimal hour of the day

    """
    if isinstance(datetime, pd.DatetimeIndex):
        if datetime.tzinfo is None:
            # If the datetime has not timezone it is assumed that is UTC
            utc_datetime = datetime.tz_localize(pytz.utc)
        else:
            utc_datetime = datetime.tz_convert(pytz.utc)
        doy = get_day_of_year(datetime)
        solar_h = ((utc_datetime.hour * 60.0) + utc_datetime.minute + utc_datetime.second / 60.0 + (
            4.0 * longitude_deg) + equation_of_time(doy)) / 60.0
        return (solar_h + 24.0) % 24.0

    if datetime.tzinfo is None:
        # If the datetime has not timezone it is assumed that is UTC
        utc_datetime = datetime.replace(tzinfo=pytz.utc)
    else:
        utc_datetime = datetime.astimezone(pytz.utc)
    doy = get_day_of_year(datetime)
    solar_h = ((utc_datetime.hour * 60.0) + utc_datetime.minute + utc_datetime.second / 60.0 + (
        4.0 * longitude_deg) + equation_of_time(doy)) / 60.0
    # Same 0-24 h range as for a DatetimeIndex (the scalar path could return negative hours or hours > 24)
    return (solar_h + 24.0) % 24.0


def equation_of_time(day_of_the_year):
    """
    Equation of time (apparent minus mean solar time):

        EoT = 9.87 sin(2B) - 7.53 cos(B) - 1.5 sin(B),  B = 2 pi (doy - 81) / 364

    Within ~1 min of the Fourier series of Spencer (1971).

    Parameters
    ----------
    day_of_the_year : Union[ndarray, float]
        Day of the year [1-366]

    Returns
    -------
    Union[ndarray, float] : equation of time [minutes]

    """
    b = (2.0 * numpy.pi / 364.0) * (day_of_the_year - 81.0)
    return (9.87 * numpy.sin(2.0 * b)) - (7.53 * numpy.cos(b)) - (1.5 * numpy.sin(b))
