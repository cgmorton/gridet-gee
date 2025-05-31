import math

import ee

from . import utils


def CalculateHourlyASCEReferenceET(
        elevation,
        temperature,
        relative_humidity,
        wind_speed,
        solar_radiation,
        extraterrestrial_radiation,
        anemometer_height,
        reference_type,
        pressure,
):
    # <summary>
    # Calculates hourly reference evapotranspiration from the ASCE standardized equation.
    # </summary>
    # <param name="elevation">Height Above Mean Sea Level (Feet)</param>
    # <param name="temperature">Average Hourly Temperature (Fahrenheit)</param>
    # <param name="relative_humidity">Average Hourly Relative Humidity (Percent)</param>
    # <param name="solar_radiation">Total Hourly Incoming Solar Radiation (Langleys)</param>
    # <param name="extraterrestrial_radiation">Total Hourly Extraterrestrial Solar Radiation (Langleys)</param>
    # <param name="wind_speed">Average Hourly Wind Speed (Miles/Hour)</param>
    # <param name="anemometer_height">Height of Anemometer from the Ground (Feet)</param>
    # <param name="reference_type">Short(Grass) or Long(Alfalfa) Reference Height</param>
    # <param name="pressure">Average Hourly Air Pressure (kiloPascals)</param>
    # <returns>Estimated Reference Evapotranspiration (Inches/Hour)</returns>
    # <remarks>Source: ASCE Standardized Reference Evapotranspiration Equation (2005)</remarks>

    # Mean Air Temperature (°C)
    t = utils.ToCelsius(temperature)

    # Elevation Above Mean Sea Level (m)
    z = utils.ToMeters(elevation)

    # Psychrometric Constant (kPa/°C)
    psi = pressure.multiply(0.000665)

    # Slope of the Saturation Vapor Pressure-Temperature Curve (kPa/°C)
    # delta = t.expression('2503 * exp(17.27 * t / (t + 237.3)) / (t + 237.3) ** 2', {'t': t})
    delta = t.add(237.3).pow(-1).multiply(t).multiply(17.27).exp().multiply(2503).divide(t.add(237.3).pow(2))

    # Saturation Vapor Pressure (kPa)
    # es = T.expression('0.6108 * exp(17.27 * T / (T + 237.3))', {'T': T})
    es = t.add(237.3).pow(-1).multiply(t).multiply(17.27).exp().multiply(0.6108)

    # Actual Vapor Pressure (kPa)
    # Conversion from Percent
    ea = relative_humidity.multiply(es).multiply(0.01)

    # Extraterrestrial Radiation (MJ/(m^2*h))
    ra = utils.ToMJperM2(extraterrestrial_radiation)

    # Calculated Clear-sky Radiation (MJ/(m^2*h))
    rso = z.multiply(0.00002).add(0.75).multiply(ra)

    # Measured Solar Radiation (MJ/(m^2*h))
    rs = utils.ToMJperM2(solar_radiation)

    # Albedo (unitless)
    albedo = 0.23

    # Net Outgoing Long-Wave Radiation (MJ/(m^2*h))
    rns = rs.multiply(1 - albedo)

    # Cloudiness Function (unitless)
    # TODO: Check if this is working correctly
    fcd = rs.divide(rso).clamp(0.3, 1).multiply(1.35).subtract(0.35)
    fcd = rs.multiply(0).add(0.7).where(rso.gt(0), fcd)

    # Stefan-Boltzmann Constant (MJ/(K^4*m^2*h))
    sigma = 0.0000000002042

    # Net Shortwave Radiation (MJ/(m^2*h))
    # rnl = t.expression(
    #   'sigma * fcd * (0.34 - 0.14 * sqrt(ea)) * pow(t, 4)',
    #   {'sigma': sigma, 'fcd': fcd, 'ea': ea, 't': utils.ToKelvin(t)}
    # )
    rnl = (
        ea.sqrt().multiply(-0.14).add(0.34).multiply(utils.ToKelvin(t).pow(4))
        .multiply(fcd).multiply(sigma)
    )

    # Net Radiation (MJ/(m^2*h))
    rn = rns.subtract(rnl)

    # TODO: Check this GEE syntax and double check the value mappings
    # Assignment of Reference Values
    # Soil Heat Flux Density (MJ/(m^2*h))
    # Numerator Constant that Changes with Reference Type and Calculation Time Step (K*mm*s^3/(Mg*h))
    # Denominator Constant that Changes with Reference Type and Calculation Time Step (s/m)
    # CGM - Assuming tall/alfalfa reference unless type is set to short
    cd_day = 0.25
    cd_night = 1.7
    g_rn_day = 0.04
    g_rn_night = 0.2
    cn = 66.0
    if reference_type.lower() in ['short', 'grass']:
        cd_day = 0.24
        cd_night = 0.96
        g_rn_day = 0.1
        g_rn_night = 0.5
        cn = 37.0

    cd = ee.Image.constant(cd_day).where(rn.lt(0), cd_night)
    g = rn.multiply(ee.Image.constant(g_rn_day).where(rn.lt(0), g_rn_night))

    # Height of Wind Measurement Above Ground Surface (m)
    zw = utils.ToMeters(anemometer_height)

    # Measured Wind Speed at zw (m/s)
    uz = utils.ToMetersPerSecond(wind_speed)

    # Adjusted Wind Speed at 2 m Above Ground Surface (m/s)
    # u2 = uz.expression('4.87 * uz / log(67.8 * zw - 5.42)', {'uz': uz, 'zw': zw})
    u2 = uz.multiply(4.87).divide(zw.multiply(67.8).subtract(5.42).log())

    # Inverse Latent Heat of Vaporization (kg/MJ)
    lhv = 0.408

    # ASCE Standardized Reference Evapotranspiration (mm/hr)
    # et = t.expression(
    #     '(lhv * delta * (rn - g) + psi * cn * u2 * (es - ea) / (t + 273)) / '
    #     '(delta + psi * (1 + cd * u2))',
    #     {
    #         'lhv': lhv, 'delta': delta, 'rn': rn, 'g': g, 'psi': psi,
    #         'cn': cn, 'cd': cd, 'u2': u2, 'es': es, 'ea': ea, 't': t,
    #     }
    # )
    et = (
        rn.subtract(g).multiply(delta).multiply(lhv)
        .add(es.subtract(ea).multiply(u2).multiply(cn).multiply(psi).divide(t.add(273)))
        .divide(u2.multiply(cd).add(1).multiply(psi).add(delta))
    )

    return utils.ToInches(et).rename(['et_reference'])


def CalculateRelativeHumidity(specific_humidity, pressure, temperature):
    # <summary>
    # Calculates relative humidity from specific humidity, pressure, and temperature.
    # </summary>
    # <param name="specific_humidity">unitless</param>
    # <param name="pressure">Pa</param>
    # <param name="temperature">C</param>
    # <returns></returns>
    # <remarks></remarks>

    return temperature.expression(
        '81.80628272 * SpecificHumidity * P / exp(35.34 * T / (2 * T + 487)) / (189 * SpecificHumidity + 311)',
        {'SpecificHumidity': specific_humidity, 'P': pressure, 'T': temperature}
    )
    # return (
    #     temperature.multiply(2).add(487).pow(-1).multiply(temperature).multiply(35.34).exp()
    #     .pow(-1).multiply(pressure).multiply(specific_humidity).multiply(81.80628272)
    #     .divide(specific_humidity.multiply(189).add(311))
    # )


def CalculateWindSpeedAdjustmentFactor(zw, h=0.12/0.3048, zu=2/0.3048):
    # <summary>
    # Adjusts wind speed from one height to another using a logarithmic profile.
    # </summary>
    # <param name="zw">Wind reference height (ft)</param>
    # <param name="h">Vegetation height (ft)</param>
    # <param name="zu">Desired wind reference height (ft)</param>
    # <returns>Factor to adjust wind speed</returns>
    # <remarks></remarks>

    # zero plane displacement height
    d = 0.67 * h

    # aerodynamic roughness length
    zom = 0.123 * h

    # return zw.expression(
    #     'log((zu - d) / zom) / log((zw - d) / zom)',
    #     {'d': d, 'zom': zom, 'zu': zu, 'zw': zw}
    # )
    return (
        zw.subtract(d).divide(zom).log().pow(-1)
        .multiply(ee.Number(zu).subtract(d).divide(zom).log())
    )


def AdjustAirPressure(pressure, temperature, delta_z):
    # <summary>
    # Adjusts pressure for change in elevation.
    # </summary>
    # <param name="pressure">Average Air Pressure at Referenced Site (Any Unit--Multiplier Relative)</param>
    # <param name="temperature">Average Air Temperature at Interpolated Site (Fahrenheit)</param>
    # <param name="delta_z">Difference in Elevation of Referenced Site from Interpolated Site (Feet)</param>

    # Acceleration of Gravity on Earth (ft/s^2)
    g = 32.174

    # Gas Constant for Air (ft-lb/(slug-R))
    r = 1716.0

    return pressure.expression(
        'P / exp(G * delta_z / (R * T))',
        {'P': pressure, 'T': utils.ToRankine(temperature), 'delta_z': delta_z, 'G': g, 'R': r}
    )
    # return pressure.divide(utils.ToRankine(temperature).multiply(r).pow(-1).multiply(delta_z).multiply(g).exp())


def CalculateDailyHargreavesReferenceET(
    minimum_temperature,
    average_temperature,
    maximum_temperature,
    extraterrestrial_radiation,
):
    # <summary>
    # Calculates daily reference evapotranspiration from the Hargreaves-Samani equation.
    # </summary>
    # <param name="minimum_temperature">Maximum Daily Temperature (Fahrenheit)</param>
    # <param name="average_temperature">Average Daily Temperature (Fahrenheit)</param>
    # <param name="maximum_temperature">Minimum Daily Temperature (Fahrenheit)</param>
    # <param name="extraterrestrial_radiation">Total Daily Extraterrestrial Solar Radiation (Langleys)</param>
    # <returns>Estimated Reference Evapotranspiration (Inches/Day)</returns>
    # <remarks>
    # Source: Hargreaves, G. H., Samani, Z. A. (1982). Estimating potential evapotranspiration.
    # Journal of the Irrigation and Drainage Division, 108(3), 225-230.
    # </remarks>

    # return average_temperature * (maximum_temperature - minimum_temperature) ^ 0.5 * extraterrestrial_radiation / 800000  # 1340000
    return (
        maximum_temperature.subtract(minimum_temperature).pow(2)
        .multiply(average_temperature)
        .multiply(extraterrestrial_radiation).divide(800000)
        .rename('hargreaves_reference_et')
    )
