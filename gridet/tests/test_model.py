import math
import pprint

import ee
import pytest

from gridet import model as model
from gridet import utils as utils


@pytest.mark.parametrize(
    'specific_humidity, pressure, temperature, expected',
    [
        # Units for this function are unitless, Pa, Celsius
        [0.00369, 85979, 30.15, 11.884309351100159],
    ]
)
def test_CalculateRelativeHumidity(specific_humidity, pressure, temperature, expected):
    output = model.CalculateRelativeHumidity(
        ee.Image(specific_humidity), ee.Image(pressure), ee.Image(temperature)
    )
    constant_geom = ee.Geometry.Rectangle([0, 0, 10, 10], 'EPSG:32613', False)
    output = output.reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)
    assert abs(float(output.getInfo()['constant']) - expected) < 0.000001


@pytest.mark.parametrize(
    'zw, expected',
    [
        # Should be 0.3048
        [10 / 0.3048, 0.747725544465109],
        [2 / 0.3048, 1.00000042805903],
    ]
)
def test_CalculateWindSpeedAdjustmentFactor(zw, expected):
    output = model.CalculateWindSpeedAdjustmentFactor(ee.Number(zw))
    assert abs(float(output.getInfo()) - expected) < 0.000001


@pytest.mark.parametrize(
    'pressure, temperature, delta_z, expected',
    [
        [85979, 86.27, 0, 85979],
        [85.979, 86.27, 0, 85.979],
        [85.979, 86.27, -10, 86.008533],
    ]
)
def test_AdjustAirPressure(pressure, temperature, delta_z, expected):
    output = model.AdjustAirPressure(ee.Image(pressure), ee.Number(temperature), delta_z)
    constant_geom = ee.Geometry.Rectangle([0, 0, 10, 10], 'EPSG:32613', False)
    output = output.reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)
    assert abs(float(output.getInfo()['constant']) - expected) < 0.0001


@pytest.mark.parametrize(
    'elev, temp, rh, wind, rs, ra, zw, type, pres, expected',
    [
        [4584, 86.27, 11.88, 0.786, 75, 105, 2 / 0.3048, 'grass', 85.979, 0.024805971],
    ]
)
def test_CalculateHourlyASCEReferenceET(elev, temp, rh, wind, rs, ra, zw, type, pres, expected):
    output = model.CalculateHourlyASCEReferenceET(
        ee.Image(elev), ee.Image(temp), ee.Image(rh), ee.Image(wind),
        ee.Image(rs), ee.Image(ra), ee.Number(zw), type, ee.Image(pres)
    )
    constant_geom = ee.Geometry.Rectangle([0, 0, 10, 10], 'EPSG:32613', False)
    output = output.reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)
    assert abs(float(output.getInfo()['et_reference']) - expected) < 0.0001
