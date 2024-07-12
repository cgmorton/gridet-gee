import math

import ee

from . import utils


def CalculateInstantaneousRa(
        longitude, latitude, elevation, slope, azimuth, temperature, pressure, solar_position
):
    # <summary>
    # Calculates instantaneous incoming extraterrestrial radiation for a given location
    # and time period on the earth for a sloped surface.
    # </summary>
    # <param name="longitude">(Decimal Degrees)</param>
    # <param name="latitude">(Decimal Degrees)</param>
    # <param name="elevation">Height Above Mean Sea Level (Feet)</param>
    # <param name="slope">Positive Inclination from Horizontal (Decimal Degrees)</param>
    # <param name="azimuth">Positive Clockwise Slope Orientation from North(0) [East(90), South(180), West(270)] (Decimal Degrees)</param>
    # <param name="temperature">Air Temperature (Fahrenheit)</param>
    # <param name="pressure">Air Pressure (kPa)</param>
    # <param name="solar_position">Earth-Sun Orientation and Distance</param>
    # <returns>Estimated Extraterrestrial Radiation for a Sloped Surface (Langleys/Hour)</returns>
    # <remarks>
    # Adapted from: Reda, I., Andreas, A. (2004). Solar position algorithm for solar radiation applications.
    # Solar energy, 76(5), 577-589.
    # </remarks>

    # Unpack the solar_position bands/variables
    # alpha - Geocentric Sun Right Ascension (radians)
    # delta - Geocentric Sun Declination (radians)
    # R - Earth-Sun Distance (astronomical units)
    # v - Apparent Sidereal Time at Greenwich (radians)
    # CGM - The order of the Solar Position list variables [R, delta, alpha, v]
    sp_r = ee.Number(solar_position[0])
    sp_delta = ee.Number(solar_position[1])
    sp_alpha = ee.Number(solar_position[2])
    sp_v = ee.Number(solar_position[3])
    # sp_r = solar_position['R']
    # sp_alpha = solar_position['alpha']
    # sp_delta = solar_position['delta']
    # sp_v = solar_position['v']
    # sp_r = solar_position.select(['R'])
    # sp_alpha = solar_position.select(['alpha'])
    # sp_delta = solar_position.select(['delta'])
    # sp_v = solar_position.select(['v'])
    # print('sp_r', sp_r.getInfo())
    # print('sp_delta', sp_delta.getInfo())
    # print('sp_alpha', sp_alpha.getInfo())
    # print('sp_v', sp_v.getInfo())

    # Variable Conversions
    sigma = utils.ToRadians(longitude)
    phi = utils.ToRadians(latitude)
    z = utils.ToMeters(elevation)
    omega = utils.ToRadians(slope)
    gamma = utils.ToRadians(azimuth.subtract(180))
    # CGM - Keeping a copy for now to keep track of the original variable renaming
    # σ As Double = ToRadians(longitude)
    # ϕ As Double = ToRadians(Latitude)
    # Z As Double = ToMeters(Elevation)
    # ω As Double = ToRadians(Slope)
    # γ As Double = ToRadians(Azimuth - 180)

    # Observer Local Hour Angle (radians)
    h = utils.LimitAngle(sigma.add(sp_v).subtract(sp_alpha))

    # Equatorial Horizontal Parallax of the Sun (radians)
    # "8.794 / (3600 * R)"
    xi = utils.ToRadians(sp_r.multiply(3600).pow(-1).multiply(8.794))

    # Parallax in the Sun Right Ascension
    # u = atan(0.99664719 * tan(phi))
    # x = cos(u) + z * cos(phi) / 6378140
    # y = 0.99664719 * sin(u) + z * sin(phi) / 6378140
    # delta_a = atan2(-x * sin(xi) * sin(h), cos(sp_delta) - x * sin(xi) * cos(h))
    u = phi.tan().multiply(0.99664719).atan()
    x = z.multiply(phi.cos()).divide(6378140).add(u.cos())
    y = u.sin().multiply(0.99664719).add(z.multiply(phi.sin()).divide(6378140))
    # CGM - Swapped x and y in GEE atan2 function
    delta_a = (
        x.multiply(xi.sin()).multiply(h.cos()).multiply(-1).add(sp_delta.cos())
        .atan2(x.multiply(-1).multiply(xi.sin()).multiply(h.sin()))
    )
    # delta_a = (
    #     x.multiply(-1).multiply(xi.sin()).multiply(h.sin())
    #     .atan2(sp_delta.cos().subtract(x.multiply(xi.sin()).multiply(h.cos())))
    # )

    # CGM - This was commented out in the source
    # Topocentric Sun Right Ascension (radians)
    # αp = SolarPosition.α.add(delta_a)

    # Topocentric Sun Declination (radians)
    # "atan2((sin(sp_delta) - y * sin(xi)) * cos(delta_a), cos(sp_delta) - x * sin(xi) * cos(h))"
    # CGM - Swapped x and y in GEE atan2 function
    delta_p = (
        x.multiply(xi.sin()).multiply(h.cos()).multiply(-1).add(sp_delta.cos())
        .atan2(y.multiply(xi.sin()).multiply(-1).add(sp_delta.sin()).multiply(delta_a.cos()))
    )

    # Topocentric Local Hour Angle (radians)
    hp = h.subtract(delta_a)

    # Topocentric Elevation Angle Without Atmospheric Refraction (radians)
    # "asin(sin(phi) * sin(delta_p) + cos(phi) * cos(delta_p) * cos(hp))"
    e0 = (
        phi.sin().multiply(delta_p.sin())
        .add(phi.cos().multiply(delta_p.cos()).multiply(hp.cos()))
        .asin()
    )

    # Atmospheric Refraction Correction (radians)
    # "(P / 101) * (283 / (273 + T)) * 1.02 / (60 * tan((e0_deg + 10.3 / (e0_deg + 5.11)) * (math.pi / 180)))"
    e0_deg = utils.ToDegrees(e0)
    delta_e = (
        e0_deg.add(5.11).pow(-1).multiply(10.3).add(e0_deg)
        .multiply(math.pi / 180).tan().multiply(60).pow(-1)
        .multiply(pressure).divide(101)
        .divide(utils.ToCelsius(temperature).add(273)).multiply(283)
        .multiply(1.02)
    )
    delta_e = utils.ToRadians(delta_e)

    # Topocentric Elevation Angle (radians)
    e = e0.add(delta_e)

    # Zenith Angle (radians)
    # "math.pi / 2 - e"
    theta = e.multiply(-1).add(math.pi / 2)

    # Topocentric Astronomers Azimuth Angle (radians)
    # "atan2(sin(hp), cos(hp) * sin(phi) - tan(delta_p) * cos(phi))"
    # CGM - Swapped x and y in GEE atan2 function
    gamma_c = hp.cos().multiply(phi.sin()).subtract(delta_p.tan().multiply(phi.cos())).atan2(hp.sin())

    # Topocentric Azimuth Angle (radians)
    # phi_c = utils.LimitAngle(gamma_c.add(math.pi))

    # Incidence Angle for a Surface Oriented in Any Direction
    # "acos(cos(theta) * cos(omega) + sin(omega) * sin(theta) * cos(gamma_c - gamma))"
    i = (
        theta.cos().multiply(omega.cos())
        .add(omega.sin().multiply(theta.sin()).multiply(gamma_c.subtract(gamma).cos()))
        .acos()
    )

    # Solar Constant (W/m^2)
    gsc = 1367.0

    # Extraterrestrial Radiation (W/m^2)
    # Updated calculation
    # Dim Ra As Double = If(θ >= 0 AndAlso θ <= π / 2 AndAlso I <= π / 2, Gsc * Math.Cos(I) / SolarPosition.R ^ 2, 0)
    ra = i.cos().multiply(gsc).divide(sp_r.pow(2))
    ra = ra.where(theta.lt(0).Or(theta.gt(math.pi / 2)).Or(i.gt(math.pi / 2)), 0)
    # # Original calculation
    # # Dim Ra As Double = If e > 0 And π / 2 - I > 0 Then Ra = Gsc * Math.Cos(I) / SolarPosition.R ^ 2
    # ra = ra.where(e.lte(0).Or(i.multiply(-1).add(math.pi / 2).lte(0)), 0)

    return utils.ToLangleysPerHour(ra).rename(['ra'])


def CalculateSolarPosition(record_date):
    # <summary>
    # Calculates the geocentric sun angles and distance at a given instant in time.
    # </summary>
    # <param name="record_date">Event Time</param>
    # <returns>Earth-Sun Orientation and Distance</returns>
    # <remarks>
    # Adapted from: Reda, I., Andreas, A. (2004).
    # Solar position algorithm for solar radiation applications. Solar energy, 76(5), 577-589.
    # </remarks>

    # Julian Day
    month = ee.Date(record_date).get('month')
    year = ee.Date(record_date).get('year')
    # CGM - Correct year before month since month is being used as a "mask"
    year = year.add(month.lt(3).multiply(-1))
    month = month.add(month.lt(3).multiply(12))
    # jd = Int(365.25 * (year + 4716)) + Int(30.6001 * (month + 1)) +
    #         RecordDate.Day + RecordDate.TimeOfDay.TotalHours / 24 - 1524.5
    jd = (
        year.add(4716).multiply(365.25).int()
        .add(month.add(1).multiply(30.6001).int())
        .add(ee.Date(record_date).get('day'))
        .add(ee.Date(record_date).getFraction('day'))
        .subtract(1524.5)
    )
    a = year.divide(100).int()
    jd = jd.add(jd.gt(2299160).multiply(a.divide(4).int().subtract(a).add(2)))
    # print('jd', jd.getInfo())

    # Approximate Difference between the Earth Rotation Time and Terrestrial Time (seconds)
    delta_t = LookupDeltaT(ee.Date(record_date))

    # Julian Day, Ephemeris
    jde = delta_t.divide(86400).add(jd)  # seconds to day fraction

    # Julian Century, Ephemeris
    jc = jd.subtract(2451545).divide(36525)
    jce = jde.subtract(2451545).divide(36525)

    # Julian Millennium, Ephemeris
    jme = jce.divide(10)

    # Earth Heliocentric Longitude (radians)
    l = utils.LimitAngle(CalculateEarthHeliocentricLongitude(jme))

    # Earth Heliocentric Latitude (radians)
    b = CalculateEarthHeliocentricLatitude(jme)

    # Earth-Sun Distance (astronomical units, au)
    r = CalculateEarthHeliocentricRadius(jme)

    # Geocentric Longitude (radians)
    theta = utils.LimitAngle(l.add(math.pi))

    # Geocentric Latitude (radians)
    beta = b.multiply(-1)

    # Mean Elongation of the Moon from the Sun (degrees)
    x = [0, 0, 0, 0, 0]
    x[0] = CalculateSolarPolynomial(jce, 1.0 / 189474, -0.0019142, 445267.11148, 297.85036)

    # Mean Anomaly of the Sun-Earth (degrees)
    x[1] = CalculateSolarPolynomial(jce, 1.0 / -300000, -0.0001603, 35999.05034, 357.52772)

    # Mean Anomaly of the Moon (degrees)
    x[2] = CalculateSolarPolynomial(jce, 1.0 / 56250, 0.0086972, 477198.867398, 134.96298)

    # Moon# s Argument of Latitude  (degrees)
    x[3] = CalculateSolarPolynomial(jce, 1.0 / 327270, -0.0036825, 483202.017538, 93.27191)

    # Longitude of the Ascending Node of the Moon’s Mean Orbit on the Ecliptic,
    # Measured from the Mean Equinox of the Date (degrees)
    x[4] = CalculateSolarPolynomial(jce, 1.0 / 450000, 0.0020708, -1934.136261, 125.04452)

    # Nutation in Longitude (radians)
    delta_psi = CalculateNutation(jce, x, 'Longitude')
    # print('delta_psi (deg)', delta_psi.multiply(180 / math.pi).getInfo(), '-0.00399840')

    # Nutation In Obliquity (radians)
    delta_epsilon = CalculateNutation(jce, x, 'Obliquity')
    # print('delta_epsilon (deg)', delta_epsilon.multiply(180 / math.pi).getInfo(), '0.00166657')

    # Mean Obliquity of the Ecliptic (arc seconds)
    u = jme.divide(10)
    epsilon0 = (
        u.multiply(2.45).add(5.79).multiply(u).add(27.87).multiply(u).add(7.12)
        .multiply(u).add(-39.05).multiply(u).add(-249.67).multiply(u).add(-51.38)
        .multiply(u).add(1999.25).multiply(u).add(-1.55).multiply(u).add(-4680.96)
        .multiply(u).add(84381.448)
    )
    # esilon0 = (
    #     84381.448 + u * (-4680.96 + u * (-1.55 + u * (1999.25 + u * (-51.38 + u * (-249.67 +
    #     u * (-39.05 + u * (7.12 + u * (27.87 + u * (5.79 + u * 2.45)))))))))

    # True Obliquity of the Ecliptic (radians)
    epsilon = utils.ToRadians(epsilon0.divide(3600)).add(delta_epsilon)

    # Aberration Correction (radians)
    delta_tau = utils.ToRadians(r.multiply(3600).pow(-1).multiply(-20.4898))

    # Apparent Sun Longitude (radians)
    # CGM - Lambda is a reserved name in python
    lam = theta.add(delta_psi).add(delta_tau)

    # Mean Sidereal Time at Greenwich (radians)
    # ν0 = 280.46061837 + 360.98564736629 * (jd - 2451545) + jc * jc * (0.000387933 - jc / 38710000)
    v0 = (
        jc.divide(38710000).multiply(-1).add(0.000387933).multiply(jc.pow(2))
        .add(jd.subtract(2451545).multiply(360.98564736629))
        .add(280.46061837)
    )
    v0 = utils.LimitAngle(utils.ToRadians(v0))

    # Apparent Sidereal Time at Greenwich (radians)
    v = delta_psi.multiply(epsilon.cos()).add(v0)

    # Geocentric Sun Right Ascension (radians)
    # CGM - Swapped x and y in atan2 function
    alpha_x = lam.sin().multiply(epsilon.cos()).subtract(beta.tan().multiply(epsilon.sin()))
    alpha = utils.LimitAngle(lam.cos().atan2(alpha_x))
    # alpha = utils.LimitAngle(atan2(sin(lam) * cos(epsilon) - tan(beta) * sin(epsilon), cos(lam)))

    # Geocentric Sun Declination (radians)
    # delta = asin(sin(beta) * cos(epsilon) + cos(beta) * sin(epsilon) * sin(lambda))
    delta = (
        beta.sin().multiply(epsilon.cos())
        .add(beta.cos().multiply(epsilon.sin()).multiply(lam.sin()))
        .asin()
    )

    # TODO: Decide on what the format/structure of the return object(s)
    return r, delta, alpha, v
    # return ee.List([r, delta, alpha, v])
    # return {'R': r, 'delta': delta, 'alpha': alpha, 'v': v}
    # return ee.Image([r, delta, alpha, v]).rename(['R', 'delta', 'alpha', 'v'])


# TODO: Update values with data for more recent years (commented out below)
def LookupDeltaT(record_date):
    # <summary>
    # Returns a yearly average difference between universal time and actual terrestrial time (1620-2014).
    # </summary>
    # <param name="record_date">Event Time</param>
    # <returns>Time Difference (seconds)</returns>
    # <remarks>
    # Source: http://asa.usno.navy.mil/SecK/DeltaT.html (link is dead)
    # Source: ftp://maia.usno.navy.mil/ser7/deltat.data
    # Source: https://maia.usno.navy.mil/products/deltaT
    # </remarks>

    delta_t = ee.List([
        # 1620
        124, 119, 115, 110, 106, 102, 98, 95, 91, 88,
        85, 82, 79, 77, 74, 72, 70, 67, 65, 63,
        62, 60, 58, 57, 55, 54, 53, 51, 50, 49,
        48, 47, 46, 45, 44, 43, 42, 41, 40, 38,
        37, 36, 35, 34, 33, 32, 31, 30, 28, 27,
        26, 25, 24, 23, 22, 21, 20, 19, 18, 17,
        16, 15, 14, 14, 13, 12, 12, 11, 11, 10,
        10, 10, 9, 9, 9, 9, 9, 9, 9, 9,
        # 1700
        9, 9, 9, 9, 9, 9, 9, 9, 10, 10,
        10, 10, 10, 10, 10, 10, 10, 11, 11, 11,
        11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
        11, 11, 11, 11, 12, 12, 12, 12, 12, 12,
        12, 12, 12, 12, 13, 13, 13, 13, 13, 13,
        13, 14, 14, 14, 14, 14, 14, 14, 15, 15,
        15, 15, 15, 15, 15, 16, 16, 16, 16, 16,
        16, 16, 16, 16, 16, 17, 17, 17, 17, 17,
        17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
        17, 17, 16, 16, 16, 16, 15, 15, 14, 14,
        # 1800
        13.7, 13.4, 13.1, 12.9, 12.7, 12.6, 12.5, 12.5, 12.5, 12.5,
        12.5, 12.5, 12.5, 12.5, 12.5, 12.5, 12.5, 12.4, 12.3, 12.2,
        12.0, 11.7, 11.4, 11.1, 10.6, 10.2, 9.6, 9.1, 8.6, 8.0,
        7.5, 7.0, 6.6, 6.3, 6.0, 5.8, 5.7, 5.6, 5.6, 5.6,
        5.7, 5.8, 5.9, 6.1, 6.2, 6.3, 6.5, 6.6, 6.8, 6.9,
        7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.7, 7.8, 7.8,
        7.88, 7.82, 7.54, 6.97, 6.40, 6.02, 5.41, 4.10, 2.92, 1.82,
        1.61, 0.10, -1.02, -1.28, -2.69, -3.24, -3.64, -4.54, -4.71, -5.11,
        -5.40, -5.42, -5.20, -5.46, -5.46, -5.79, -5.63, -5.64, -5.80, -5.66,
        -5.87, -6.01, -6.19, -6.64, -6.44, -6.47, -6.09, -5.76, -4.66, -3.74,
        # 1900
        -2.72, -1.54, -0.02, 1.24, 2.64, 3.86, 5.37, 6.14, 7.75, 9.13,
        10.46, 11.53, 13.36, 14.65, 16.01, 17.20, 18.24, 19.06, 20.25, 20.95,
        21.16, 22.25, 22.41, 23.03, 23.49, 23.62, 23.86, 24.49, 24.34, 24.08,
        24.02, 24.00, 23.87, 23.95, 23.86, 23.93, 23.73, 23.92, 23.96, 24.02,
        24.33, 24.83, 25.30, 25.70, 26.24, 26.77, 27.28, 27.78, 28.25, 28.71,
        29.15, 29.57, 29.97, 30.36, 30.72, 31.07, 31.35, 31.68, 32.18, 32.68,
        33.15, 33.59, 34.00, 34.47, 35.03, 35.73, 36.54, 37.43, 38.29, 39.20,
        40.18, 41.17, 42.23, 43.37, 44.49, 45.48, 46.46, 47.52, 48.53, 49.59,
        50.54, 51.38, 52.17, 52.96, 53.79, 54.34, 54.87, 55.32, 55.82, 56.30,
        56.86, 57.57, 58.31, 59.12, 59.99, 60.78, 61.63, 62.30, 62.97, 63.47,
        # 2000
        63.83, 64.09, 64.30, 64.47, 64.57, 64.69, 64.85, 65.15, 65.46, 65.78,
        66.20, 66.45, 66.74, 67.09, 67.45, 67.84, 68.35, 68.78, 69.09, 66.20,
        66.45, 66.74, 67.09, 67.45, 67.84, 68.35, 68.78, 69.09,
        # # CGM - Updated values from new DeltaT link starting in 1974
        # #   and predicted values for 2025-2033
        # # 1900
        # -2.72, -1.54, -0.02, 1.24, 2.64, 3.86, 5.37, 6.14, 7.75, 9.13,
        # 10.46, 11.53, 13.36, 14.65, 16.01, 17.20, 18.24, 19.06, 20.25, 20.95,
        # 21.16, 22.25, 22.41, 23.03, 23.49, 23.62, 23.86, 24.49, 24.34, 24.08,
        # 24.02, 24.00, 23.87, 23.95, 23.86, 23.93, 23.73, 23.92, 23.96, 24.02,
        # 24.33, 24.83, 25.30, 25.70, 26.24, 26.77, 27.28, 27.78, 28.25, 28.71,
        # 29.15, 29.57, 29.97, 30.36, 30.72, 31.07, 31.349, 31.677, 32.166, 32.671,
        # 33.150, 33.584, 33.992, 34.466, 35.030, 35.738, 36.546, 37.429, 38.291, 39.204,
        # 40.18, 41.17, 42.23, 43.37, 44.4841, 45.4761, 46.4567, 47.5214, 48.5344, 49.5861,
        # 50.5387, 51.3808, 52.1668, 52.9565, 53.7882, 54.3427, 54.8713, 55.3222, 55.8197, 56.3000,
        # 56.8553, 57.5653, 58.3092, 59.1218, 59.9845, 60.7853, 61.6287, 62.2950, 62.9659, 63.4673,
        # # 2000
        # 63.8285, 64.0908, 64.2998, 64.4734, 64.5736, 64.6876, 64.8452, 65.1464, 65.4573, 65.7768,
        # 66.0699, 66.3246, 66.6030, 66.9069, 67.2810, 67.6439, 68.1024, 68.5927, 68.9676, 69.2202,
        # 69.3612, 69.3594, 69.2945, 69.2039, 69.1752, 69.04, 69.05, 69.14, 69.34, 69.63,
        # 69.97, 70.32, 70.62, 70.98,
    ])

    y = record_date.get('year').int().subtract(1620).max(0).min(delta_t.size().subtract(1))

    return ee.Number(delta_t.get(y))


def CalculateEarthHeliocentricLongitude(jme):
    # <summary>
    # Calculates earth heliocentric longitude, latitude, or earth-sun distance.
    # </summary>
    # <param name="jme">Julian Millennium</param>
    # <remarks>
    # Adapted from: Reda, I., Andreas, A. (2004). Solar position algorithm for solar radiation applications.
    # Solar energy, 76(5), 577-589.
    # </remarks>

    # Swapped the arrays to make the indexing easier
    abc = [
        [
            [
                175347046, 3341656, 34894, 3497, 3418, 3136, 2676, 2343, 1324, 1273,
                1199, 990, 902, 857, 780, 753, 505, 492, 357, 317,
                284, 271, 243, 206, 205, 202, 156, 132, 126, 115,
                103, 102, 102, 99, 98, 86, 85, 85, 80, 79,
                75, 74, 74, 70, 62, 61, 57, 56, 56, 52, 52, 51, 49, 41, 41, 39, 37, 37, 36, 36, 33, 30, 30, 25
            ],
            [
                0, 4.6692568, 4.6261, 2.7441, 2.8289, 3.6277, 4.4181, 6.1352, 0.7425, 2.0371,
                1.1096, 5.233, 2.045, 3.508, 1.179, 2.533, 4.583, 4.205, 2.92, 5.849,
                1.899, 0.315, 0.345, 4.806, 1.869, 2.458, 0.833, 3.411, 1.083, 0.645,
                0.636, 0.976, 4.267, 6.21, 0.68, 5.98, 1.3, 3.67, 1.81, 3.04,
                1.76, 3.5, 4.68, 0.83, 3.98, 1.82, 2.78, 4.39, 3.47, 0.19,
                1.33, 0.28, 0.49, 5.37, 2.4, 6.17, 6.04, 2.57, 1.71, 1.78,
                0.59, 0.44, 2.74, 3.16
            ],
            [
                0, 6283.07585, 12566.1517, 5753.3849, 3.5231, 77713.7715, 7860.4194, 3930.2097, 11506.7698, 529.691,
                1577.3435, 5884.927, 26.298, 398.149, 5223.694, 5507.553, 18849.228, 775.523, 0.067, 11790.629,
                796.298, 10977.079, 5486.778, 2544.314, 5573.143, 6069.777, 213.299, 2942.463, 20.775, 0.98,
                4694.003, 15720.839, 7.114, 2146.17, 155.42, 161000.69, 6275.96, 71430.7, 17260.15, 12036.46,
                5088.63, 3154.69, 801.82, 9437.76, 8827.39, 7084.9, 6286.6, 14143.5, 6279.55, 12139.55,
                1748.02, 5856.48, 1194.45, 8429.24, 19651.05, 10447.39, 10213.29, 1059.38, 2352.87, 6812.77,
                17789.85, 83996.85, 1349.87, 4690.48
            ],
        ],
        [
            [
                628331966747, 206059, 4303, 425, 119, 109, 93, 72, 68, 67,
                59, 56, 45, 36, 29, 21, 19, 19, 17, 16,
                16, 15, 12, 12, 12, 12, 11, 10, 10, 9,
                9, 8, 6, 6
            ],
            [
                0, 2.678235, 2.6351, 1.59, 5.796, 2.966, 2.59, 1.14, 1.87, 4.41,
                2.89, 2.17, 0.4, 0.47, 2.65, 5.34, 1.85, 4.97, 2.99, 0.03,
                1.43, 1.21, 2.83, 3.26, 5.27, 2.08, 0.77, 1.3, 4.24, 2.7,
                5.64, 5.3, 2.65, 4.67
            ],
            [
                0, 6283.07585, 12566.1517, 3.523, 26.298, 1577.344, 18849.23, 529.69, 398.15, 5507.55,
                5223.69, 155.42, 796.3, 775.52, 7.11, 0.98, 5486.78, 213.3, 6275.96, 2544.31,
                2146.17, 10977.08, 1748.02, 5088.63, 1194.45, 4694, 553.57, 6286.6, 1349.87, 242.73,
                951.72, 2352.87, 9437.76, 4690.48
            ],
        ],
        [
            [
                52919, 8720, 309, 27, 16, 16, 10, 9, 7, 5,
                4, 4, 3, 3, 3, 3, 3, 3, 2, 2
            ],
            [
                0, 1.0721, 0.867, 0.05, 5.19, 3.68, 0.76, 2.06, 0.83, 4.66,
                1.03, 3.44, 5.14, 6.05, 1.19, 6.12, 0.31, 2.28, 4.38, 3.75
            ],
            [
                0, 6283.0758, 12566.152, 3.52, 26.3, 155.42, 18849.23, 77713.77, 775.52, 1577.34,
                7.11, 5573.14, 796.3, 5507.55, 242.73, 529.69, 398.15, 553.57, 5223.69, 0.98
            ],
        ],
        [
            [289, 35, 17, 3, 1, 1, 1],
            [5.844, 0, 5.49, 5.2, 4.72, 5.3, 5.97],
            [6283.076, 0, 12566.15, 155.42, 3.52, 18849.23, 242.73],
        ],
        [[114, 8, 1], [3.142, 4.13, 3.84], [0, 6283.08, 12566.15]],
        [[1], [3.14], [0]],
    ]
    abc_n = [64, 34, 20, 7, 3, 1]

    l0 = (
        ee.Array(ee.List.repeat(jme, abc_n[0]))
        .multiply(ee.Array(abc[0][2])).add(ee.Array(abc[0][1])).cos()
        .multiply(ee.Array(abc[0][0]))
        .reduce(ee.Reducer.sum(), [0]).get([0])
    )
    l1 = (
        ee.Array(ee.List.repeat(jme, abc_n[1]))
        .multiply(ee.Array(abc[1][2])).add(ee.Array(abc[1][1])).cos()
        .multiply(ee.Array(abc[1][0]))
        .reduce(ee.Reducer.sum(), [0]).get([0])
    )
    l2 = (
        ee.Array(ee.List.repeat(jme, abc_n[2]))
        .multiply(ee.Array(abc[2][2])).add(ee.Array(abc[2][1])).cos()
        .multiply(ee.Array(abc[2][0]))
        .reduce(ee.Reducer.sum(), [0]).get([0])
    )
    l3 = (
        ee.Array(ee.List.repeat(jme, abc_n[3]))
        .multiply(ee.Array(abc[3][2])).add(ee.Array(abc[3][1])).cos()
        .multiply(ee.Array(abc[3][0]))
        .reduce(ee.Reducer.sum(), [0]).get([0])
    )
    l4 = (
        ee.Array(ee.List.repeat(jme, abc_n[4]))
        .multiply(ee.Array(abc[4][2])).add(ee.Array(abc[4][1])).cos()
        .multiply(ee.Array(abc[4][0]))
        .reduce(ee.Reducer.sum(), [0]).get([0])
    )
    l5 = (
        ee.Array(ee.List.repeat(jme, abc_n[5]))
        .multiply(ee.Array(abc[5][2])).add(ee.Array(abc[5][1])).cos()
        .multiply(ee.Array(abc[5][0]))
        .reduce(ee.Reducer.sum(), [0]).get([0])
    )

    EarthHeliocentricValue = (
        l0.add(jme.multiply(l1))
        .add(jme.pow(2).multiply(l2))
        .add(jme.pow(3).multiply(l3))
        .add(jme.pow(4).multiply(l4))
        .add(jme.pow(5).multiply(l5))
        .divide(100000000)
    )

    return utils.LimitAngle(EarthHeliocentricValue)


def CalculateEarthHeliocentricLatitude(jme):
    # <summary>
    # Calculates earth heliocentric latitude.
    # </summary>
    # <param name="jme">Julian Millennium</param>
    # <remarks>
    # Adapted from: Reda, I., Andreas, A. (2004). Solar position algorithm for solar radiation applications.
    # Solar energy, 76(5), 577-589.
    # </remarks>

    # Swapped the arrays to make the indexing easier
    abc = [
        [
            [280.0, 102.0, 80.0, 44.0, 32.0],
            [3.199, 5.422, 3.88, 3.7, 4.0],
            [84334.662, 5507.553, 5223.69, 2352.87, 1577.34],
        ],
        [
            [9.0, 6.0],
            [3.9, 1.73],
            [5507.55, 5223.69],
        ],
    ]
    abc_n = [5, 2]

    # CGM - There has to be a cleaner way to do this...
    #   Can this be built in a for loop?
    b0 = (
        ee.Array(ee.List.repeat(jme, abc_n[0]))
        .multiply(ee.Array(abc[0][2])).add(ee.Array(abc[0][1])).cos()
        .multiply(ee.Array(abc[0][0]))
        .reduce(ee.Reducer.sum(), [0]).get([0])
    )
    b1 = (
        ee.Array(ee.List.repeat(jme, abc_n[1]))
        .multiply(ee.Array(abc[1][2])).add(ee.Array(abc[1][1])).cos()
        .multiply(ee.Array(abc[1][0]))
        .reduce(ee.Reducer.sum(), [0]).get([0])
    )
    # print('b0', b0.getInfo())
    # print('b1', b1.getInfo())

    return b1.multiply(jme).add(b0).divide(100000000)


def CalculateEarthHeliocentricRadius(jme):
    # <summary>
    # Calculates earth heliocentric earth-sun distance.
    # </summary>
    # <param name="jme">Julian Millennium</param>
    # <param name="Type">Earth Heliocentric Property Type</param>
    # <remarks>
    # Adapted from: Reda, I., Andreas, A. (2004). Solar position algorithm for solar radiation applications.
    # Solar energy, 76(5), 577-589.
    # </remarks>

    # Transposed the arrays from the original source to make the indexing easier
    abc = [
        [
            [
                100013989, 1670700, 13956, 3084, 1628, 1576, 925, 542, 472, 346,
                329, 307, 243, 212, 186, 175, 110, 98, 86, 86,
                65, 63, 57, 56, 49, 47, 45, 43, 39, 38,
                37, 37, 36, 35, 33, 32, 32, 28, 28, 26
            ],
            [
                0, 3.0984635, 3.05525, 5.1985, 1.1739, 2.8469, 5.453, 4.564, 3.661, 0.964,
                5.9, 0.299, 4.273, 5.847, 5.022, 3.012, 5.055, 0.89, 5.69, 1.27,
                0.27, 0.92, 2.01, 5.24, 3.25, 2.58, 5.54, 6.01, 5.36, 2.39,
                0.83, 4.9, 1.67, 1.84, 0.24, 0.18, 1.78, 1.21, 1.9, 4.59
            ],
            [
                0, 6283.07585, 12566.1517, 77713.7715, 5753.3849, 7860.4194, 11506.77, 3930.21, 5884.927, 5507.553,
                5223.694, 5573.143, 11790.629, 1577.344, 10977.079, 18849.228, 5486.778, 6069.78, 15720.84, 161000.69,
                17260.15, 529.69, 83996.85, 71430.7, 2544.31, 775.52, 9437.76, 6275.96, 4694, 8827.39,
                19651.05, 12139.55, 12036.46, 2942.46, 7084.9, 5088.63, 398.15, 6286.6, 6279.55, 10447.39
            ],
        ],
        [
            [103019, 1721, 702, 32, 31, 25, 18, 10, 9, 9],
            [1.10749, 1.0644, 3.142, 1.02, 2.84, 1.32, 1.42, 5.91, 1.42, 0.27],
            [6283.07585, 12566.1517, 0, 18849.23, 5507.55, 5223.69, 1577.34, 10977.08, 6275.96, 5486.78],
        ],
        [
            [4359, 124, 12, 9, 6, 3],
            [5.7846, 5.579, 3.14, 3.63, 1.87, 5.47],
            [6283.0758, 12566.152, 0, 77713.77, 5573.14, 18849.23],
        ],
        [
            [145, 7], [4.273, 3.92], [6283.076, 12566.15],
        ],
        [
            [4], [2.56], [6283.08],
        ],
    ]
    # Hardcoded lengths
    abc_n = [40, 10, 6, 2, 1]

    r0 = (
        ee.Array(ee.List.repeat(jme, abc_n[0]))
        .multiply(ee.Array(abc[0][2])).add(ee.Array(abc[0][1])).cos()
        .multiply(ee.Array(abc[0][0]))
        .reduce(ee.Reducer.sum(), [0]).get([0])
    )
    r1 = (
        ee.Array(ee.List.repeat(jme, abc_n[1]))
        .multiply(ee.Array(abc[1][2])).add(ee.Array(abc[1][1])).cos()
        .multiply(ee.Array(abc[1][0]))
        .reduce(ee.Reducer.sum(), [0]).get([0])
    )
    r2 = (
        ee.Array(ee.List.repeat(jme, abc_n[2]))
        .multiply(ee.Array(abc[2][2])).add(ee.Array(abc[2][1])).cos()
        .multiply(ee.Array(abc[2][0]))
        .reduce(ee.Reducer.sum(), [0]).get([0])
    )
    r3 = (
        ee.Array(ee.List.repeat(jme, abc_n[3]))
        .multiply(ee.Array(abc[3][2])).add(ee.Array(abc[3][1])).cos()
        .multiply(ee.Array(abc[3][0]))
        .reduce(ee.Reducer.sum(), [0]).get([0])
    )
    r4 = (
        ee.Array(ee.List.repeat(jme, abc_n[4]))
        .multiply(ee.Array(abc[4][2])).add(ee.Array(abc[4][1])).cos()
        .multiply(ee.Array(abc[4][0]))
        .reduce(ee.Reducer.sum(), [0]).get([0])
    )

    return (
        r0.add(r1.multiply(jme))
        .add(r2.multiply(jme.pow(2)))
        .add(r3.multiply(jme.pow(3)))
        .add(r4.multiply(jme.pow(4)))
        .divide(100000000)
    )


def CalculateSolarPolynomial(x, a, b, c, d):
    # <summary>
    # Calculates a third order polynomial for earth solar calculations (A*V^3+B*V^2+C*V+D).
    # </summary>
    # <param name="x">Input Variable Raised to a Power</param>
    # <param name="a">3rd Order Constant</param>
    # <param name="b">2nd Order Constant</param>
    # <param name="c">1st Order Constant</param>
    # <param name="d">Constant</param>
    # <returns>Calculated Value</returns>

    return utils.LimitAngleDeg(x.multiply(a).add(b).multiply(x).add(c).multiply(x).add(d))


def CalculateNutation(jce, x, nutation_type):
    # <summary>
    # Calculates nutation in longitude or obliquity for an input Julian Century.
    # </summary>
    # <param name="jce">Julian Century</param>
    # <param name="x">Earth-Moon-Sun Properties Array</param>
    # <param name="NutationType">Nutation Type</param>
    # <returns>Nutation</returns>
    # <remarks>
    # Adapted from: Reda, I., Andreas, A. (2004). Solar position algorithm for solar radiation applications.
    # Solar energy, 76(5), 577-589.
    # </remarks>

    y = ee.Array([
        [
            0, -2, 0, 0, 0, 0, -2, 0, 0, -2, -2, -2, 0, 2, 0, 2, 0, 0, -2, 0,
            2, 0, 0, -2, 0, -2, 0, 0, 2, -2, 0, -2, 0, 0, 2, 2, 0, -2, 0, 2,
            2, -2, -2, 2, 2, 0, -2, -2, 0, -2, -2, 0, -1, -2, 1, 0, 0, -1, 0, 0,
            2, 0, 2
        ],
        [
            0, 0, 0, 0, 1, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 1, 0, -1, 0, 0, 0, 1, 1, -1, 0,
            0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, -1, 1, -1,
            -1, 0, -1
        ],
        [
            0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, -1, 0, 1, -1, -1, 1, 2, -2,
            0, 2, 2, 1, 0, 0, -1, 0, -1, 0, 0, 1, 0, 2, -1, 1, 0, 1, 0, 0,
            1, 2, 1, -2, 0, 1, 0, 0, 2, 2, 0, 1, 1, 0, 0, 1, -2, 1, 1, 1,
            -1, 3, 0
        ],
        [
            0, 2, 2, 0, 0, 0, 2, 2, 2, 2, 0, 2, 2, 0, 0, 2, 0, 2, 0, 2,
            2, 2, 0, 2, 2, 2, 2, 0, 0, 2, 0, 0, 0, -2, 2, 2, 2, 0, 2, 2,
            0, 2, 2, 0, 0, 0, 2, 0, 2, 0, 2, -2, 0, 0, 0, 2, 2, 0, 0, 2,
            2, 2, 2
        ],
        [
            1, 2, 2, 2, 0, 0, 2, 1, 2, 2, 0, 1, 2, 0, 1, 2, 1, 1, 0, 1,
            2, 2, 0, 2, 0, 0, 1, 0, 1, 2, 1, 1, 1, 0, 1, 2, 2, 0, 2, 1,
            0, 2, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 2,
            2, 2, 2
        ],
    ])
    # y = ee.Array([
    #   [0, 0, 0, 0, 1], [-2, 0, 0, 2, 2], [0, 0, 0, 2, 2], [0, 0, 0, 0, 2],
    #   [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [-2, 1, 0, 2, 2], [0, 0, 0, 2, 1],
    #   [0, 0, 1, 2, 2], [-2, -1, 0, 2, 2], [-2, 0, 1, 0, 0], [-2, 0, 0, 2, 1],
    #   [0, 0, -1, 2, 2], [2, 0, 0, 0, 0], [0, 0, 1, 0, 1], [2, 0, -1, 2, 2],
    #   [0, 0, -1, 0, 1], [0, 0, 1, 2, 1], [-2, 0, 2, 0, 0], [0, 0, -2, 2, 1],
    #   [2, 0, 0, 2, 2], [0, 0, 2, 2, 2], [0, 0, 2, 0, 0], [-2, 0, 1, 2, 2],
    #   [0, 0, 0, 2, 0], [-2, 0, 0, 2, 0], [0, 0, -1, 2, 1], [0, 2, 0, 0, 0],
    #   [2, 0, -1, 0, 1], [-2, 2, 0, 2, 2], [0, 1, 0, 0, 1], [-2, 0, 1, 0, 1],
    #   [0, -1, 0, 0, 1], [0, 0, 2, -2, 0], [2, 0, -1, 2, 1], [2, 0, 1, 2, 2],
    #   [0, 1, 0, 2, 2], [-2, 1, 1, 0, 0], [0, -1, 0, 2, 2], [2, 0, 0, 2, 1],
    #   [2, 0, 1, 0, 0], [-2, 0, 2, 2, 2], [-2, 0, 1, 2, 1], [2, 0, -2, 0, 1],
    #   [2, 0, 0, 0, 1], [0, -1, 1, 0, 0], [-2, -1, 0, 2, 1], [-2, 0, 0, 0, 1],
    #   [0, 0, 2, 2, 1], [-2, 0, 2, 0, 1], [-2, 1, 0, 2, 1], [0, 0, 1, -2, 0],
    #   [-1, 0, 1, 0, 0], [-2, 1, 0, 0, 0], [1, 0, 0, 0, 0], [0, 0, 1, 2, 0],
    #   [0, 0, -2, 2, 2], [-1, -1, 1, 0, 0], [0, 1, 1, 0, 0], [0, -1, 1, 2, 2],
    #   [2, -1, -1, 2, 2], [0, 0, 3, 2, 2], [2, -1, 0, 2, 2]
    # ])
    a = ee.Array([
        -171996, -13187, -2274, 2062, 1426, 712, -517, -386, -301, 217,
        -158, 129, 123, 63, 63, -59, -58, -51, 48, 46,
        -38, -31, 29, 29, 26, -22, 21, 17, 16, -16,
        -15, -13, -12, 11, -10, -8, 7, -7, -7, -7,
        6, 6, 6, -6, -6, 5, -5, -5, -5, 4,
        4, 4, -4, -4, -4, 3, -3, -3, -3, -3,
        -3, -3, -3
    ])
    b = ee.Array([
        -174.2, -1.6, -0.2, 0.2, -3.4, 0.1, 1.2, -0.4, 0, -0.5,
        0, 0.1, 0, 0, 0.1, 0, -0.1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, -0.1, 0, 0.1,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0
    ])
    c = ee.Array([
        92025, 5736, 977, -895, 54, -7, 224, 200, 129, -95,
        0, -70, -53, 0, -33, 26, 32, 27, 0, -24,
        16, 13, 0, -12, 0, 0, -10, 0, -8, 7,
        9, 7, 6, 0, 5, 3, -3, 0, 3, 3,
        0, -3, -3, 3, 3, 0, 3, 3, 3, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0
    ])
    d = ee.Array([
        8.9, -3.1, -0.5, 0.5, -0.1, 0, -0.6, 0, -0.1, 0.3,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0
    ])

    x_array = ee.Array(x).repeat(1, 63)
    # x_array = ee.Array(x).repeat(1, 63).transpose()
    sigma_xy = ee.Array(y.multiply(x_array).reduce(ee.Reducer.sum(), [0]).toList().get(0))

    if nutation_type.lower() == 'longitude':
        # CGM - Switching to list to drop the empty dimension
        #   How do you do this with gee arrays?
        delta_psi = (
            sigma_xy.multiply(math.pi / 180).sin()
            .multiply(b.multiply(jce).add(a))
            .reduce(ee.Reducer.sum(), [0]).get([0])
        )
        return utils.ToRadians(delta_psi.divide(36000000))
    elif nutation_type.lower() == 'obliquity':
        delta_e = (
            sigma_xy.multiply(math.pi / 180).cos()
            .multiply(d.multiply(jce).add(c))
            .reduce(ee.Reducer.sum(), [0]).get([0])
        )
        return utils.ToRadians(delta_e.divide(36000000))
    else:
        raise ValueError(f'Unsupported nutation_type: {nutation_type}')


def TranslateRs(rs_source, ra_source, ra_destination, slope):
    # <summary>
    # Translates incoming solar radiation from one site to another taking into account slope differences
    # (Source interpolation site assumed nearly horizontal).
    # </summary>
    # <param name="rs_source">Incoming Solar Radiation at Source Interpolation Site (Same Units as Other Input Radiations)</param>
    # <param name="ra_source">Extraterrestrial Solar Radiation at Source Interpolation Site (Same Units as Other Input Radiations)</param>
    # <param name="ra_destination">Extraterrestrial Solar Radiation at Source Interpolation Site (Same Units as Other Input Radiations)</param>
    # <param name="slope">Positive Inclination from Horizontal (Decimal Degrees)</param>
    # <returns>Incoming Solar Radiation at Destination Translation Site (Same Units as Input Radiations)</returns>
    # <remarks>
    # Adapted from: Allen, R. G., Trezza, R., Tasumi, M. (2006).
    # Analytical integrated functions for daily solar radiation on slopes.
    # Agricultural and Forest Meteorology, 139(1), 55-73.
    # </remarks>

    # TODO: Need to check for zeros or check for this somehow
    # if rs_source <= 0 or ra_source <= 0:
    #     return 0
    # else:

    # Wherever Ra is <=0, set Rs to 0, then set Ra to 0.01 to avoid divide by zero
    rs_source = rs_source.where(ra_source.lte(0), 0).max(0)
    ra_source = ra_source.max(0.01)

    # Convert to Radians
    slope = utils.ToRadians(slope)

    # Actual Atmospheric Transmissivity (Unitless)
    tau_sw = rs_source.divide(ra_source)

    # Source Actual Direct Radiation Index (Unitless)
    # CGM - The 0.28 coefficient should have a negative sign
    #   but leaving as positive to match GridET code for now
    kb = tau_sw.multiply(0.765).add(0.828).multiply(tau_sw).add(0.28).multiply(tau_sw).add(0.022)
    # kb = tau_sw.multiply(0.765).add(0.828).multiply(tau_sw).subtract(0.28).multiply(tau_sw).add(0.022)
    # kb = tau_sw.expression("((0.765 * b() + 0.828) * b() - 0.28) * b() + 0.022")
    kb = kb.where(tau_sw.gte(0.42), tau_sw.multiply(1.56).subtract(0.55))
    kb = kb.where(tau_sw.lte(0.175), 0.175)

    # Source Actual Diffuse Radiation Index (Unitless)
    kd = tau_sw.subtract(kb)

    # Direct Radiation Fraction (Unitless)
    fb = ra_destination.divide(ra_source)

    # Reflectance Factor (Unitless)
    # 0.75 + 0.25 * cos(slope) - 0.5 * slope / math.pi
    fi = slope.cos().multiply(0.25).subtract(slope.multiply(0.5).divide(math.pi)).add(0.75)

    # Sky-view Factor (Unitless)
    fia = kd.expression(
        '(1 - kb) * (1 + ((kb / (kb + kd)) ** 0.5) * ((sin(slope / 2)) ** 3)) * fi + (fb * kb)',
        {'fi': fi, 'fb': fb, 'kb': kb, 'kd': kd, 'slope': slope}
    )

    # Incoming Solar Radiation at Interpolated Site (Source Site Units)
    rs_destination = rs_source.expression(
        'rs_source * ((fb * kb + fia * kd) / tau_sw + 0.23 * (1 - fi))',
        {'fb': fb, 'fi': fi, 'fia': fia, 'kb': kb, 'kd': kd, 'rs_source': rs_source, 'tau_sw': tau_sw}
    )
    # rs_destination = (
    #     kb.multiply(fb).add(fia.multiply(kd)).divide(tau_sw)
    #     .add(fi.multiply(-1).add(1).multiply(0.23)).multiply(rs_source)
    # )

    return rs_destination
