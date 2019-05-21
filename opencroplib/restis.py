def longwave_downwelling_1(air_temperature, relative_humidity):
    """
    Function to calculate the longwave downwelling radiation from the air temperature and relative humidity
    Reference:


    Parameters
    ----------
    relative_humidity : float
     Relative humidity [%]
    air_temperature : float
     Air temperature [Celsius]

    Returns
    -------
    float: longwave downwelling radiation [W/m2]

    """

    ea_t = ea_calc(air_temperature, relative_humidity) * 0.01

    lwd_sim = (0.64 + 0.044 * numpy.sqrt(ea_t)) * \
              BOLTZMANN * numpy.power(air_temperature + 273.15, 4)
    return lwd_sim


def longwave_downwelling_2(air_temperature, relative_humidity):
    """
    Function to calculate the longwave downwelling radiation from the air temperature and relative humidity

    Reference...

    Parameters
    ----------
    relative_humidity : float
        Relative humidity [%]
    air_temperature : float
        Air temperature [Celsius]

    Returns
    -------
    float: longwave downwelling radiation [W/m2]

    """

    ea_t = ea_calc(air_temperature, relative_humidity) * 0.01

    lwd_sim = (0.7 + 0.0000595 * ea_t * numpy.exp(1500.0 / (air_temperature + 273.15))) * \
              BOLTZMANN * numpy.power(air_temperature + 273.15, 4)
    return lwd_sim


def longwave_downwelling_3(air_temperature, relative_humidity):
    """
    Function to calculate the longwave downwelling radiation from the air temperature and relative humidity

    Reference...

    Parameters
    ----------
    relative_humidity : float
        Relative humidity [%]
    air_temperature : float
        Air temperature [Celsius]

    Returns
    -------
    float: longwave downwelling radiation [W/m2]

    """

    ea_t = ea_calc(air_temperature, relative_humidity) * 0.01

    lwd_sim = (0.3 + 0.7 * 1.24 * numpy.power(ea_t / (air_temperature + 273.15), 1.0 / 7.0)) * \
              BOLTZMANN * numpy.power(air_temperature + 273.15, 4)
    return lwd_sim
