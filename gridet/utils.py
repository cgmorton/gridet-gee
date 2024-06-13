import math

import ee


def LimitAngle(value):
    # <summary>
    # Redefines an angle between 0 to 2π.
    # </summary>
    # <param name="value">Input angle (Radians [or User-Defined])</param>
    # <returns>Basic Positive Angle</returns>
    upper_bounds = 2 * math.pi
    output = value.mod(upper_bounds)
    return output.add(output.lt(0).multiply(upper_bounds))

    # CGM - This approach will handle an arbitrary number of negative "wraps" but is probably overkill
    # lower_bounds = 0
    # upper_bounds = 2 * math.pi
    # limit_range = upper_bounds - lower_bounds
    # return balue.subtract(lower_bounds).mod(limit_range).add(limit_range).mod(limit_range).add(lower_bounds)


# TODO: This could be combined with function above by setting UpperBounds as a function default
def LimitAngleDeg(value):
    # <summary>
    # Redefines an angle between 0 to 360.
    # </summary>
    # <param name="value">Input angle (Radians [or User-Defined])</param>
    # <returns>Basic Positive Angle</returns>
    upper_bounds = 360
    output = value.mod(upper_bounds)
    return output.add(output.lt(0).multiply(upper_bounds))


# TODO: Implement a more dynamic approach
def SumProduct(Array1, Array2):
    return (
        ee.Number(Array1[0]).multiply(Array2[0])
        .add(Array1[1].multiply(Array2[1]))
        .add(Array1[2].multiply(Array2[2]))
        .add(Array1[3].multiply(Array2[3]))
        .add(Array1[4].multiply(Array2[4]))
    )


# Unit Conversions
def ToFahrenheit(Celsius):
    return Celsius.multiply(1.8).add(32)


def ToCelsius(Fahrenheit):
    return Fahrenheit.subtract(32).divide(1.8)


def ToRadians(Degrees):
    return Degrees.multiply(math.pi / 180)


def ToDegrees(Radians):
    return Radians.multiply(180 / math.pi)


# CGM - Feet to Meters conversion should be 0.3048 exactly
def ToMeters(Feet):
    return Feet.multiply(0.304800609601219)


def ToFeet(Meters):
    return Meters.divide(0.304800609601219)


def ToInches(Millimeters):
    return Millimeters.divide(25.4)


def ToMJperM2(Langleys):
    return Langleys.multiply(0.04184)


def ToKelvin(Celsius):
    return Celsius.add(273.16)


def ToRankine(Fahrenheit):
    return Fahrenheit.add(459.67)


def ToLangleysPerHour(WattsPerSquaredMeter):
    return WattsPerSquaredMeter.divide(11.63)


def ToMilesPerHour(MetersPerSecond):
    return MetersPerSecond.divide(0.44704)


def ToMetersPerSecond(MilesPerHour):
    return MilesPerHour.multiply(0.44704)
