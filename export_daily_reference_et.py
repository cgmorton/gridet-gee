import argparse
from datetime import datetime, timedelta
import logging
import pprint
import time

import ee

logging.getLogger('earthengine-api').setLevel(logging.INFO)
logging.getLogger('googleapiclient').setLevel(logging.INFO)
logging.getLogger('requests').setLevel(logging.INFO)
logging.getLogger('urllib3').setLevel(logging.INFO)


def main(
        start_dt,
        end_dt,
        delay_time=0,
        project_id=None,
        overwrite_flag=False,
        reverse_flag=False,
        noclamp_flag=True,
        round_flag=True,
):
    """"""

    logging.info('\nBuild GridET daily reference ET assets')
    logging.info(f'Start: {start_dt.strftime("%Y-%m-%d")}')
    logging.info(f'End:   {end_dt.strftime("%Y-%m-%d")}\n')

    if project_id is None:
        ee.Initialize(ee.ServiceAccountCredentials(
            '_', key_file='/Users/Charles.Morton@dri.edu/Projects/keys/openet-dri-gee.json')
        )
    else:
        ee.Initialize(project=project_id)

    if noclamp_flag:
        import gridet.nldas_noclamp as nldas
        output_coll_id = 'projects/openet/assets/reference_et/utah/gridet/daily/v0_windspeed'
        export_id_fmt = 'gridet_daily_{date}_noclamp'

    else:
        import gridet.nldas as nldas
        output_coll_id = 'projects/openet/assets/reference_et/utah/gridet/daily/v0_test'
        export_id_fmt = 'gridet_daily_{date}',

    if not ee.data.getInfo(output_coll_id):
        logging.info('\nBuilding output image collection')
        ee.data.createAsset({'type': 'IMAGE_COLLECTION'}, output_coll_id)

    # Load the GridET ancillary assets (elevation, lat, lon, slope, aspect)
    mask_img = ee.Image('projects/openet/assets/reference_et/utah/gridet/ancillary/mask')
    mask_info = mask_img.getInfo()
    mask_shape = mask_info['bands'][0]['dimensions']
    mask_crs = mask_info['bands'][0]['crs']
    mask_transform = mask_info['bands'][0]['crs_transform']

    # Process each date
    for export_dt in sorted(date_range(start_dt, end_dt), reverse=reverse_flag):
        logging.info(f'{export_dt.strftime("%Y-%m-%d")}')

        nldas_coll = (
            ee.ImageCollection('NASA/NLDAS/FORA0125_H002')
            .filterDate(ee.Date(export_dt), ee.Date(export_dt).advance(1, 'day'))
        )

        # Compute the hourly image collection
        # Define the day as the UTC day (0 UTC to 0 UTC)
        hourly_eto_coll = nldas_coll.map(nldas.eto)
        hourly_etr_coll = nldas_coll.map(nldas.etr)

        # Compute daily and convert to mm (from inches)
        output_img = (
            hourly_eto_coll.sum().multiply(25.4).rename('eto')
            .addBands(hourly_etr_coll.sum().multiply(25.4).rename('etr'))

        )
        # pprint.pprint(hourly_eto_coll.first().reduceRegion(**project_rr).getInfo())
        # pprint.pprint('GridET New: {}'.format(output_img.reduceRegion(**project_rr).getInfo()['eto']))

        if round_flag:
            output_img = output_img.multiply(10000).round().divide(10000).float()

        output_img = output_img.set({
            'system:time_start': ee.Date(export_dt).millis(),
            'units_eto': 'mm',
            'units_etr': 'mm',
        })

        export_id = f'{output_coll_id}/{export_dt.strftime("%Y%m%d")}'
        if ee.data.getInfo(export_id):
            if overwrite_flag:
                logging.info('  Removing existing asset')
                ee.data.deleteAsset(export_id)
            else:
                logging.debug('  Asset already exists, skipping')
                continue

        try:
            task = ee.batch.Export.image.toAsset(
                image=output_img,
                description=export_id_fmt.format(date=export_dt.strftime("%Y%m%d")),
                assetId=export_id,
                dimensions=mask_shape,
                crs=mask_crs,
                crsTransform=mask_transform,
                maxPixels=10000000000,
            )
            logging.info('  Starting export task')
            task.start()
        except:
            logging.info('  Unhandled exception starting export task, skipping')

        time.sleep(delay_time)

    # # Code for debugging
    # test_point = ee.Geometry.Point(-113.5613, 37.19)
    # test_point = ee.Geometry.Point(-113.5625, 37.1875)
    # project_rr = {
    #     'reducer': ee.Reducer.first(), 'geometry': test_point,
    #     'crs': mask_crs, 'crsTransform': mask_transform,
    # }
    # nldas_rr = {
    #     'reducer': ee.Reducer.first(), 'geometry': test_point,
    #     'crs': 'EPSG:4326', 'crsTransform': [0.125, 0, -125, 0, -0.125, 53],
    # }

    # def pprint_point_values(image):
    #     xy = [-113.94379, 40.06814]
    #     dx, dy = 0.006, 0.005
    #     points = ee.FeatureCollection([
    #         ee.Feature(ee.Geometry.Point(xy), None) for xy in [
    #             [xy[0] - dx, xy[1] + dy], [xy[0] + 0, xy[1] + dy], [xy[0] + dx, xy[1] + dy],
    #             [xy[0] - dx, xy[1] + 0],  [xy[0] + 0, xy[1] + 0],  [xy[0] + dx, xy[1] + 0],
    #             [xy[0] - dx, xy[1] - dy], [xy[0] + 0, xy[1] - dy], [xy[0] + dx, xy[1] - dy]]
    #     ])
    #     point_values = (
    #         image
    #         .reduceRegions(**{
    #             'collection': points, 'reducer': ee.Reducer.first(),
    #             'crs': mask_crs, 'crsTransform': mask_transform})
    #         .getInfo()
    #     )
    #     point_values = [ftr['properties']['first'] for ftr in point_values['features']]
    #     print(f'{point_values[0]:0.06f} {point_values[1]:0.06f} {point_values[2]:0.06f}')
    #     print(f'{point_values[3]:0.06f} {point_values[4]:0.06f} {point_values[5]:0.06f}')
    #     print(f'{point_values[6]:0.06f} {point_values[7]:0.06f} {point_values[8]:0.06f}')
    #
    #
    # elevation = ee.Image('projects/openet/assets/reference_et/utah/gridet/ancillary/elevation')
    # aspect = ee.Image('projects/openet/assets/reference_et/utah/gridet/ancillary/aspect')
    # slope = ee.Image('projects/openet/assets/reference_et/utah/gridet/ancillary/slope')
    # latitude = ee.Image('projects/openet/assets/reference_et/utah/gridet/ancillary/latitude')
    # longitude = ee.Image('projects/openet/assets/reference_et/utah/gridet/ancillary/longitude')
    #
    # print('\nElevation')
    # pprint_point_values(elevation)
    #
    # print('\nAspect')
    # pprint_point_values(aspect)
    #
    # print('\nSlope')
    # pprint_point_values(slope)
    #
    # print('\nLat/Lon')
    # pprint_point_values(latitude)
    # pprint_point_values(longitude)

    # for export_dt in date_range(start_dt, end_dt):
    #     logging.info(f'{export_dt.strftime("%Y-%m-%d")}')
    #
    #     # CGM - Compute each hour image separately
    #     # hourly_eto_list = []
    #     sum_img = mask_img.multiply(0)
    #     for h in range(0, 24):
    #         print(h)
    #         hour_dt = ee.Date(export_dt).advance(h, 'hour')
    #         hour_img = nldas.eto(nldas_coll.filterDate(hour_dt, hour_dt.advance(1, 'hour')).first())
    #         hour_img = hour_img.multiply(25.4)
    #         sum_img = sum_img.add(hour_img)
    #         pprint_point_values(hour_img)
    #         # pprint_point_values(sum_img)
    #         # output = hour_img.reduceRegion(**project_rr).getInfo()['eto']
    #         # hourly_eto_list.append(output)
    #
    #     # print('Total:', sum(hourly_eto_list) * 25.4)
    #     print('\nSum')
    #     pprint_point_values(sum_img)
    #
    #
    #     # Pull the original GridET value at the test point to compare against
    #     gridet_img = ee.ImageCollection('projects/openet/assets/reference_et/utah/gridet/daily/v0')\
    #         .filterDate(export_dt.strftime("%Y-%m-%d")).first().select(['eto'])
    #     # pprint.pprint('GridET Old: {}'.format(gridet_img.reduceRegion(**project_rr).getInfo()['eto']))
    #     print('\nGridET')
    #     pprint_point_values(gridet_img)


def date_range(start_dt, end_dt, days=1, skip_leap_days=False):
    """Generate dates within a range (inclusive)

    Parameters
    ----------
    start_dt : datetime
        start date.
    end_dt : datetime
        end date.
    days : int, optional
        Step size (the default is 1).
    skip_leap_days : bool, optional
        If True, skip leap days while incrementing (the default is False).

    Yields
    ------
    datetime

    """
    import copy
    curr_dt = copy.copy(start_dt)
    while curr_dt <= end_dt:
        if not skip_leap_days or curr_dt.month != 2 or curr_dt.day != 29:
            yield curr_dt
        curr_dt += timedelta(days=days)


def arg_valid_date(input_date):
    """Check that a date string is ISO format (YYYY-MM-DD)"""
    try:
        return datetime.strptime(input_date, '%Y-%m-%d')
    except ValueError:
        msg = f'Not a valid date: "{input_date}".'
        raise argparse.ArgumentTypeError(msg)


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='GridET asset export',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--start', type=arg_valid_date, metavar='DATE', default=None,
        help='Start date (format YYYY-MM-DD)')
    parser.add_argument(
        '--end', type=arg_valid_date, metavar='DATE', default=None,
        help='End date (format YYYY-MM-DD)')
    parser.add_argument(
        '--project', type=str,
        # required=True
        help='Earth Engine Project ID')
    parser.add_argument(
        '--delay', default=0, type=float,
        help='Delay (in seconds) between each export tasks')
    parser.add_argument(
        '--overwrite', default=False, action='store_true',
        help='Force overwrite of existing files')
    parser.add_argument(
        '--reverse', default=False, action='store_true',
        help='Process dates in reverse order')
    parser.add_argument(
        '--debug', default=logging.INFO, const=logging.DEBUG,
        help='Debug level logging', action='store_const', dest='loglevel'
    )
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = arg_parse()
    logging.basicConfig(level=args.loglevel, format='%(message)s')

    main(
        start_dt=args.start,
        end_dt=args.end,
        delay_time=args.delay,
        project_id=args.project,
        overwrite_flag=args.overwrite,
        reverse_flag=args.reverse,
    )
