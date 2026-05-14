#import argparse
from datetime import datetime, timedelta
import logging
#import math
import pprint

import ee

import gridet.nldas as nldas

ee.Initialize(project='ee-cmorton')

nldas_coll = ee.ImageCollection('NASA/NLDAS/FORA0125_H002')

# Load the GridET ancillary assets (elevation, lat, lon, slope, aspect)
mask_img = ee.Image('projects/openet/assets/reference_et/utah/gridet/ancillary/mask')
mask_info = mask_img.getInfo()
mask_shape = mask_info['bands'][0]['dimensions']
mask_crs = mask_info['bands'][0]['crs']
mask_transform = mask_info['bands'][0]['crs_transform']

# test_point = ee.Geometry.Point(-113.31386, 39.56243)
# test_point = ee.Geometry.Point(-113.4343, 39.5636)
# test_point = ee.Geometry.Point(-113.3114, 39.1846)
# test_point = ee.Geometry.Point(-113.4315, 39.1894)
# test_point = ee.Geometry.Point(-113.5606, 39.1888)

test_point = ee.Geometry.Point(-110.061, 40.0659)


project_rr = {
    'reducer': ee.Reducer.first(), 'geometry': test_point,
    'crs': mask_crs, 'crsTransform': mask_transform,
}
nldas_rr = {
    'reducer': ee.Reducer.first(), 'geometry': test_point,
    'crs': 'EPSG:4326', 'crsTransform': [0.125, 0, -125, 0, -0.125, 53],
}

export_dt = datetime(2017, 7, 1)
print(f'{export_dt.strftime("%Y-%m-%d")}')

# nldas_ns = 13914.936349159196
# nldas_ns = 13914.936349159196 / 2
# interp_offsets = [
#     [0, 0],
#     # [-nldas_ns, 0],
#     [-nldas_ns, nldas_ns],
#     # [0, nldas_ns],
#     [nldas_ns, nldas_ns],
#     # [nldas_ns, 0],
#     [nldas_ns, -nldas_ns],
#     # [0, -nldas_ns],
#     [-nldas_ns, -nldas_ns],
# ]


# Pull the original GridET value at the test point to compare against
gridet_img = ee.ImageCollection('projects/openet/assets/reference_et/utah/gridet/daily/v0')\
    .filterDate(export_dt.strftime("%Y-%m-%d")).first().select(['eto'])
pprint.pprint('GridE Old: {}'.format(gridet_img.reduceRegion(**project_rr).getInfo()['eto']))


# for interp_offset in interp_offsets:
#     print(interp_offset)

# Compute each hour image separately
hourly_eto_list = []
for h in range(0, 24):
    print(h)
    hour_dt = ee.Date(export_dt).advance(h, 'hour')
    output = (
        nldas.etsz(nldas_coll.filterDate(hour_dt, hour_dt.advance(1, 'hour')).first(), 'grass')
        .reduceRegion(**project_rr).getInfo()['et_reference']
    )
    hourly_eto_list.append(output)
print(f' {sum(hourly_eto_list) * 25.4}')


