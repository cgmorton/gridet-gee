import argparse
from datetime import datetime, timedelta
import logging
import pprint
import time

import ee

import gridet.nldas as nldas
import gridet.utils

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
):
    """"""

    logging.info('\nBuild GridET daily meteorology assets')
    logging.info(f'Start: {start_dt.strftime("%Y-%m-%d")}')
    logging.info(f'End:   {end_dt.strftime("%Y-%m-%d")}\n')

    if project_id is None:
        ee.Initialize(ee.ServiceAccountCredentials(
            '_', key_file='/Users/Charles.Morton@dri.edu/Projects/keys/openet-dri-gee.json')
        )
    else:
        ee.Initialize(project=project_id)

    output_coll_id = 'projects/openet/assets/meteorology/gridet/daily'
    export_id_fmt = 'gridet_meteorology_daily_{date}'

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
        hourly_coll = nldas_coll.map(nldas.adjust)

        # Compute the daily aggregations
        # Convert the units to match the default NLDAS units
        output_img = (
            ee.Image([
                gridet.utils.ToCelsius(hourly_coll.select(['temperature']).reduce(ee.Reducer.min())),
                gridet.utils.ToCelsius(hourly_coll.select(['temperature']).reduce(ee.Reducer.max())),
                hourly_coll.select(['relative_humidity']).reduce(ee.Reducer.mean()),
                gridet.utils.ToMetersPerSecond(hourly_coll.select(['wind_speed']).reduce(ee.Reducer.mean())),
                hourly_coll.select(['shortwave_radiation']).reduce(ee.Reducer.sum()).multiply(11.6300),
                hourly_coll.select(['extraterrestrial_radiation']).reduce(ee.Reducer.sum()).multiply(11.6300),
                hourly_coll.select(['pressure']).reduce(ee.Reducer.mean()).multiply(1000),
            ])
            .rename([
                'temperature_minimum',
                'temperature_maximum',
                'relative_humidity',
                'wind_speed',
                'shortwave_radiation',
                'extraterrestrial_radiation',
                'pressure',
            ])
            .set({
                'system:time_start': ee.Date(export_dt).millis(),
                'units_temperature_minimum': 'C',
                'units_temperature_maximum': 'C',
                'units_relative_humidity': 'percent',
                'units_wind_speed': 'm s-1',
                'units_shortwave_radiation': 'W m-2',
                'units_extraterrestrial_radiation': 'W m-2',
                'units_pressure': 'Pa',
                # 'units_temperature_minimum': 'Fahrenheit',
                # 'units_temperature_maximum': 'Fahrenheit',
                # 'units_relative_humidity': 'Percent',
                # 'units_wind_speed': 'Miles/Hour',
                # 'units_shortwave_radiation': 'Langleys',
                # 'units_extraterrestrial_radiation': 'Langleys',
                # 'units_pressure': 'KiloPascals',
            })
        )

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
        except Exception as e:
            logging.info('  Unhandled exception starting export task, skipping')
            logging.debug(f'  {e}')

        time.sleep(delay_time)


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
        '--project', type=str, help='Earth Engine Project ID')
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
