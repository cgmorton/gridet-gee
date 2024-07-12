import math
import pprint

import ee
import pytest

from gridet import nldas as nldas
#from gridet import model as model
#from gridet import utils as utils


@pytest.mark.parametrize(
    'date_str, surface, xy, expected',
    [
        ['2017-07-01T00:00:00', 'grass', [-113.94379, 40.06814], 0.0162212552250172],
        ['2017-07-01T18:00:00', 'grass', [-113.94379, 40.06814], 0.0221100612997057],
        ['2017-07-01T20:00:00', 'grass', [-113.94379, 40.06814], 0.0256871606629055],
    ]
)
def test_etsz(date_str, surface, xy, expected):
    nldas_img = (
        ee.ImageCollection('NASA/NLDAS/FORA0125_H002')
        .filterDate(date_str, ee.Date(date_str).advance(1, 'hour'))
        .first()
    )
    mask_img = ee.Image('projects/openet/assets/reference_et/utah/gridet/ancillary/mask')
    rr_args = {
        'reducer': ee.Reducer.first(), 'geometry': ee.Geometry.Point(xy),
        'crs':  mask_img.projection().wkt(),
        'crsTransform': ee.List(ee.Dictionary(ee.Algorithms.Describe(mask_img.projection())).get('transform'))
    }

    # Build a new NLDAS image with the test point values burned in
    nldas_values = nldas_img.reduceRegion(**rr_args).getInfo()
    nldas_img = ee.Image([
        nldas_img.select(['temperature']).multiply(0).add(nldas_values['temperature']),
        nldas_img.select(['pressure']).multiply(0).add(nldas_values['pressure']),
        nldas_img.select(['specific_humidity']).multiply(0).add(nldas_values['specific_humidity']),
        nldas_img.select(['shortwave_radiation']).multiply(0).add(nldas_values['shortwave_radiation']),
        nldas_img.select(['wind_u']).multiply(0).add(nldas_values['wind_u']),
        nldas_img.select(['wind_v']).multiply(0).add(nldas_values['wind_v']),
    ]).set({'system:time_start': nldas_img.get('system:time_start')})

    output = nldas.etsz(nldas_img, surface).reduceRegion(**rr_args)

    assert abs(float(output.getInfo()['et_reference']) - expected) < 0.0001
