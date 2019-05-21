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
    Calculate day of the year (DOY) for a given datetime object

    Parameters
    ----------
    datetime : Union[pd.DatetimeIndex, dt.datetime]
        datetime representing the date

    Returns
    -------
    object : day of the year as days after first of January

    """
    if isinstance(datetime, pd.DatetimeIndex):
        return datetime.dayofyear
    else:
        year_start = dt.datetime(datetime.year, 1, 1, tzinfo=datetime.tzinfo)
        delta = (datetime - year_start)
        return delta.days


def get_solar_time(longitude_deg, datetime):
    """
    Calculates the decimal hour for a given location and datetime

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
    return ((utc_datetime.hour * 60.0) + utc_datetime.minute + utc_datetime.second / 60.0 + (
        4.0 * longitude_deg) + equation_of_time(doy)) / 60.0


def equation_of_time(day_of_the_year):
    """

    Parameters
    ----------
    day_of_the_year

    Returns
    -------

    """
    b = (2.0 * numpy.pi / 364.0) * (day_of_the_year - 81.0)
    return (9.87 * numpy.sin(2.0 * b)) - (7.53 * numpy.cos(b)) - (1.5 * numpy.sin(b))
