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
        [0.0036495, 83615.86, 27.439999999999998, 13.377897386726037],
    ]
)
def test_CalculateRelativeHumidity(specific_humidity, pressure, temperature, expected):
    output = model.CalculateRelativeHumidity(
        ee.Image(specific_humidity), ee.Image(pressure), ee.Image(temperature)
    )
    constant_geom = ee.Geometry.Rectangle([0, 0, 10, 10], 'EPSG:32613', False)
    output = output.reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)
    assert abs(float(output.getInfo()['constant']) - expected) < 0.000000000001


@pytest.mark.parametrize(
    'zw, expected',
    [
        # Should be 0.3048
        # [10 / 0.3048, 0.747725544465109],
        # [2 / 0.3048, 1.00000042805903],
        [2 / 0.304800609601219, 1.0000004280590313],
        [10 / 0.304800609601219, 0.747725544465109],
    ]
)
def test_CalculateWindSpeedAdjustmentFactor(zw, expected):
    output = model.CalculateWindSpeedAdjustmentFactor(ee.Number(zw))
    assert abs(float(output.getInfo()) - expected) < 0.000000000001


@pytest.mark.parametrize(
    'pressure, temperature, delta_z, expected',
    [
        # [85979, 86.27, 0, 85979],
        # [85.979, 86.27, 0, 85.979],
        # [85.979, 86.27, -10, 86.008533],
        #
        [83.61585999999998, 79.83377333896347, -50.028510677078884, 83.7613647326532],
    ]
)
def test_AdjustAirPressure(pressure, temperature, delta_z, expected):
    output = model.AdjustAirPressure(ee.Image(pressure), ee.Number(temperature), delta_z)
    constant_geom = ee.Geometry.Rectangle([0, 0, 10, 10], 'EPSG:32613', False)
    output = output.reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)
    assert abs(float(output.getInfo()['constant']) - expected) < 0.000000000001


@pytest.mark.parametrize(
    'elev, temp, rh, wind, rs, ra, zw, type, pres, expected',
    [
        #
        [
            5421.52587890625, 79.00013596481276, 19.82057137642429, 4.218869373485886,
            42.367205269590094, 67.70220269275103, 6.561666666666671, 'grass', 83.6761147644775,
            0.016221255353246
        ],
        [
            5421.52587890625, 77.89298241625437, 23.927610818177225, 2.6562464950566755,
            64.7354304608692, 97.88128445826067, 6.561666666666671, 'grass', 83.55629280487929,
            0.0221100611749851
        ],
        [
            5421.52587890625, 84.19681253900778, 18.318759434188863, 3.0152093309544292,
            70.78540815587768, 107.9367407529658, 6.561666666666671, 'grass', 83.52770220243846,
            0.0256871607772338
        ],
    ]
)
def test_CalculateHourlyASCEReferenceET(elev, temp, rh, wind, rs, ra, zw, type, pres, expected):
    output = model.CalculateHourlyASCEReferenceET(
        ee.Image(elev), ee.Image(temp), ee.Image(rh), ee.Image(wind),
        ee.Image(rs), ee.Image(ra), ee.Number(zw), type, ee.Image(pres)
    )
    constant_geom = ee.Geometry.Rectangle([0, 0, 10, 10], 'EPSG:32613', False)
    output = output.reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)
    assert abs(float(output.getInfo()['et_reference']) - expected) < 0.000000000001
