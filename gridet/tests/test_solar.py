import math
import pprint

import ee
import pytest

from gridet import solar as solar
from gridet import utils as utils


@pytest.mark.parametrize(
    'date_str, expected',
    [
        ['2003-10-17', 64.47],
        ['1600-01-01', 124],
        ['1620-01-01', 124],
        ['2000-01-01', 63.83],
        ['2014-01-01', 67.43],
        ['2015-01-01', 67.43],
        ['2022-01-01', 67.43],
    ]
)
def test_LookupDeltaT(date_str, expected):
    assert solar.LookupDeltaT(ee.Date(date_str)).getInfo() == expected


@pytest.mark.parametrize(
    'value, expected',
    [
        [0.0037927819922933584, 24.0182616917],
    ]
)
def test_CalculateEarthHeliocentricLongitude(value, expected):
    output = solar.CalculateEarthHeliocentricLongitude(ee.Number(value))
    assert abs(float(output.multiply(180 / math.pi).getInfo()) - expected) < 0.0000000001


@pytest.mark.parametrize(
    'value, expected',
    [
        [0.0037927819922933584, -0.0001011219],
    ]
)
def test_CalculateEarthHeliocentricLatitude(value, expected):
    output = solar.CalculateEarthHeliocentricLatitude(ee.Number(value))
    assert abs(float(output.multiply(180 / math.pi).getInfo()) - expected) < 0.0000000001


@pytest.mark.parametrize(
    'value, expected',
    [
        [0.0037927819922933584, 0.9965422974],
    ]
)
def test_CalculateEarthHeliocentricRadius(value, expected):
    output = solar.CalculateEarthHeliocentricRadius(ee.Number(value))
    assert abs(float(output.getInfo()) - expected) < 0.0000000001


@pytest.mark.parametrize(
    'value, nutation_type, expected',
    [
        [0.037927819922933585, 'Longitude', -0.00399840],
        [0.037927819922933585, 'Obliquity', 0.00166657],
    ]
)
def test_CalculateNutation(value, nutation_type, expected):
    x = [
        ee.Number(265.86117906490836),
        ee.Number(282.89321846136477),
        ee.Number(234.07570261126602),
        ee.Number(60.07101228227839),
        ee.Number(51.686951165383405)
    ]
    output = solar.CalculateNutation(ee.Number(value), x, nutation_type)
    assert abs(float(output.multiply(180 / math.pi).getInfo()) - expected) < 0.00000001


# TODO: Test that ValueError is raised for invalid nutation type values


@pytest.mark.parametrize(
    'x, a, b, c, d, expected',
    [
        [0.2249555313, 1.0 / 189474, -0.0019142, 445267.11148, 297.85036, 23.14989],
    ]
)
def test_CalculateSolarPolynomial(x, a, b, c, d, expected):
    output = solar.CalculateSolarPolynomial(ee.Number(x), a, b, c, d)
    assert abs(float(output.getInfo()) - expected) < 0.0001


# TODO: Maybe move jd calculation to a separate def maybe to support testing
# # Test values for jd from here:
# #   https:#ssd.jpl.nasa.gov/tools/jdc/#/cd
# print(CalculateSolarPosition(ee.Date('2022-07-01')), '2459761.5000000')
# print(CalculateSolarPosition(ee.Date('2022-07-01').advance(1, 'hour')), '2459761.5416667')
# print(CalculateSolarPosition(ee.Date('2022-07-01').advance(19.5, 'hour')), '2459762.3125000')


@pytest.mark.parametrize(
    'date, expected',
    [
        ['2003-10-17T19:30:30', [0.9965423054291, -0.16256571547167, 3.52953362654697, 5.55908153353535]],
        ['2017-07-01T00:00:00', [1.0166579614132, 0.40324232094146, 1.74866608416633, 4.87361879669341]],
        ['2017-07-01T18:00:00', [1.0166660750518, 0.40233145063824, 1.76219620718522, 3.31572436194462]],
        ['2017-07-01T20:00:00', [1.0166668153730, 0.40222616971529, 1.76369882812728, 3.84075668424932]],
    ]
)
def test_CalculateSolarPosition(date, expected):
    R, delta, alpha, v = solar.CalculateSolarPosition(date)
    assert abs(float(R.getInfo()) - expected[0]) < 0.0000000001
    assert abs(float(delta.getInfo()) - expected[1]) < 0.0000000001
    assert abs(float(alpha.getInfo()) - expected[2]) < 0.0000000001
    assert abs(float(v.getInfo()) - expected[3]) < 0.0000000001


@pytest.mark.parametrize(
    'date, lon, lat, elev, slope, azimuth, temp, pressure, expected',
    [
        # ['2003-10-17T19:30:30', -105.1786, 39.742476, 1830.14, 30, 170, 11.0, 820, 107.1048697],
        ['2017-07-01T18:00:00', -112.688, 39.313, 4584, 0.04, 275.6, 86.27, 85.979, 102.328243428094],
        ['2017-07-01T20:00:00', -112.688, 39.313, 4584, 0.04, 275.6, 86.27, 85.979, 108.68201013007],
        ['2017-07-01T00:00:00', -113.43322516890046, 41.81669468299659, 5654.087675195883,
         2.1556215286254883, 192.59628295898438, 68.22186585091022, 83.05958497082149, 62.5452024483686],
    ]
)
def test_CalculateInstantaneousRa(date, lon, lat, elev, slope, azimuth, temp, pressure, expected):
    solar_position = solar.CalculateSolarPosition(date)
    output = solar.CalculateInstantaneousRa(
        ee.Image(lon),          # deg
        ee.Image(lat),          # deg
        ee.Image(elev),         # ft
        ee.Image(slope),        # deg
        ee.Image(azimuth),      # deg
        ee.Image(temp),         # F
        ee.Image(pressure),     # kpa
        solar_position
    )
    constant_geom = ee.Geometry.Rectangle([0, 0, 10, 10], 'EPSG:32613', False)
    output = output.reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)

    assert abs(float(output.getInfo()['ra']) - expected) < 0.00000001


@pytest.mark.parametrize(
    'rs_source, ra_source, ra_destination, slope, expected',
    [
        [75, 105, 105, 0, 75],
        [50, 105, 105, 0, 50],
        [100, 105, 105, 0, 100],
        [75, 105, 105, 1, 75.0292537253124],
        [50, 105, 105, 1, 49.9647767536174],
        [100, 105, 105, 1, 100.064447932062],
        [75, 105, 110, 0, 78.2446428571429],
        [50, 105, 110, 0, 51.2375],
        [100, 105, 110, 0, 104.756547619048],
        [75, 105, 100, 0, 71.7553571428571],
        [50, 105, 100, 0, 48.7625],
        [100, 105, 100, 0, 95.2434523809524],
        [75, 105, 100, 1, 71.7846108681696],
        #
        [0, 105, 105, 0, 0],
        [75, 0, 105, 0, 0],
        [0, 0, 105, 0, 0],
        [-1, 105, 105, 0, 0],
        [75, -1, 105, 0, 0],
        [-1, -1, 105, 0, 0],
        #
        [0.42061625, 0.0567711, 0.38273549857, 2.26058959961, -8.91727921654038],
        [0.42061625, 0.0567711, 0, 2.26058959961, 2.03225504837265],
    ]
)
def test_TranslateRs(rs_source, ra_source, ra_destination, slope, expected):
    output = solar.TranslateRs(
        ee.Image(rs_source), ee.Image(ra_source), ee.Image(ra_destination), ee.Image(slope)
    )
    constant_geom = ee.Geometry.Rectangle([0, 0, 10, 10], 'EPSG:32613', False)
    output = output.reduceRegion(ee.Reducer.first(), geometry=constant_geom, scale=1)

    assert abs(float(output.getInfo()['constant']) - expected) < 0.01
