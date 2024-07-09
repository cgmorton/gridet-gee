import math
import pprint

import ee
import pytest

from gridet import utils as utils


def test_ToFahrenheit():
    assert float(utils.ToFahrenheit(ee.Number(0)).getInfo()) == 32
    assert float(utils.ToFahrenheit(ee.Number(20)).getInfo()) == 68
    assert float(utils.ToFahrenheit(ee.Number(100)).getInfo()) == 212


def test_ToRadians():
    assert float(utils.ToRadians(ee.Number(0)).getInfo()) == 0
    assert abs(float(utils.ToRadians(ee.Number(180)).getInfo()) - math.pi) < 0.0001
    assert abs(float(utils.ToRadians(ee.Number(360)).getInfo()) - (2 * math.pi)) < 0.0001


def test_ToFeet():
    # Note, the feet converion in GridET is for survey feet
    assert float(utils.ToFeet(ee.Number(2)).getInfo()) == 6.561666666666671
    # assert float(utils.ToFeet(ee.Number(2)).getInfo()) == 2 / 0.3048


@pytest.mark.parametrize(
    'angle, expected',
    [
        # All seem reasonable except that the range is [0, 2*PI)
        #   so 2*PI is converted to 0 (is that okay?)
        [0, 0],
        [math.pi, math.pi],
        [2 * math.pi - 0.1, 2 * math.pi - 0.1],
        [2 * math.pi, 0],
        [2 * math.pi + 0.1, 0.1],
        [-0.1, 2 * math.pi - 0.1],
    ]
)
def test_LimitAngle(angle, expected, tol=0.01):
    assert abs(float(utils.LimitAngle(ee.Number(angle)).getInfo()) - expected) < tol


@pytest.mark.parametrize(
    'angle, expected',
    [
        # All seem reasonable except that the range is [0, 2*PI)
        #   so 2*PI is converted to 0 (is that okay?)
        [0, 0],
        [180, 180],
        [359.9, 359.9],
        [360, 0],
        [360.1, 0.1],
        [-1, 359],
    ]
)
def test_LimitAngleDeg(angle, expected, tol=0.01):
    assert abs(float(utils.LimitAngleDeg(ee.Number(angle)).getInfo()) - expected) < tol


def test_SumProduct(a=[1, 2, 3, 4, 5], b=[6, 7, 8, 9, 10], expected=130, tol=0.01):
    # The items in the first array must be ee.Number() or ee.Image() objects
    assert abs(float(utils.SumProduct([ee.Number(x) for x in a], b).getInfo()) - expected) < tol
