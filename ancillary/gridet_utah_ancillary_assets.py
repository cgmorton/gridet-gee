import argparse
import logging

import ee


def main(project_id, overwrite_flag=False):
    """"""

    logging.info('Build/Ingest GridET Ancillary Assets')

    ee.Initialize(project=project_id)

    mask_id = 'projects/openet/assets/reference_et/utah/gridet/ancillary/mask'

    # The 10m 3DEP image will not reduce in a single call
    elev_src_img = ee.Image('NASA/NASADEM_HGT/001').select('elevation')
    # elev_src_img = ee.Image('USGS/3DEP/10m')

    # TODO: Get the folder ID from the mask ID instead?
    folder_id = 'projects/openet/assets/reference_et/utah/gridet/ancillary'
    elevation_id = f'{folder_id}/elevation'
    latitude_id = f'{folder_id}/latitude'
    longitude_id = f'{folder_id}/longitude'
    slope_id = f'{folder_id}/slope'
    aspect_id = f'{folder_id}/aspect'
    # azimuth_id = f'{folder_id}/azimuth'

    mask_img = ee.Image(mask_id).unmask(1)

    shape = mask_img.getInfo()['bands'][0]['dimensions']
    crs = mask_img.projection().wkt()
    transform = ee.List(ee.Dictionary(ee.Algorithms.Describe(mask_img.projection())).get('transform'))
    # extent = ee.Geometry.Rectangle([-116, 35.5, -107, 43.5], 'EPSG:4326', False)

    # The multiply(1) is to drop the properties on the source elevation image
    elev_img = (
        elev_src_img
        .reduceResolution(ee.Reducer.mean().unweighted(), False, 65535)
        .reproject(crs, transform)
        .multiply(1).rename('elevation')
    )

    if overwrite_flag or not ee.data.getInfo(elevation_id):
        logging.info('elevation')
        if ee.data.getInfo(elevation_id):
            ee.data.deleteAsset(elevation_id)

        task = ee.batch.Export.image.toAsset(
            image=elev_img,
            description='gridet_elevation_asset',
            assetId=elevation_id,
            dimensions=shape,
            crs=crs,
            crsTransform=transform,
            maxPixels=10000000000,
            # pyramidingPolicy=,
            # shardSize=,
        )
        task.start()

    if overwrite_flag or not ee.data.getInfo(latitude_id):
        logging.info('latitude')
        if ee.data.getInfo(latitude_id):
            ee.data.deleteAsset(latitude_id)

        task = ee.batch.Export.image.toAsset(
            image=mask_img.multiply(0).add(ee.Image.pixelLonLat().select('latitude')).rename(['latitude']),
            description='gridet_latitude_asset',
            assetId=latitude_id,
            dimensions=shape,
            crs=crs,
            crsTransform=transform,
            maxPixels=10000000000,
            # pyramidingPolicy=,
        )
        task.start()

    if overwrite_flag or not ee.data.getInfo(longitude_id):
        logging.info('longitude')
        if ee.data.getInfo(longitude_id):
            ee.data.deleteAsset(longitude_id)

        task = ee.batch.Export.image.toAsset(
            image=mask_img.multiply(0).add(ee.Image.pixelLonLat().select('longitude')).rename(['longitude']),
            description='gridet_longitude_asset',
            assetId=longitude_id,
            dimensions=shape,
            crs=crs,
            crsTransform=transform,
            maxPixels=10000000000,
            # pyramidingPolicy=,
        )
        task.start()

    if overwrite_flag or not ee.data.getInfo(slope_id):
        logging.info('slope')
        if ee.data.getInfo(slope_id):
            ee.data.deleteAsset(slope_id)

        # Slope units are degrees, range is [0,90).
        slope_img = ee.Terrain.slope(elev_img)

        task = ee.batch.Export.image.toAsset(
            image=slope_img,
            description='gridet_slope_asset',
            assetId=slope_id,
            dimensions=shape,
            crs=crs,
            crsTransform=transform,
            maxPixels=10000000000,
            # pyramidingPolicy=,
        )
        task.start()

    if not ee.data.getInfo(aspect_id) or overwrite_flag:
        logging.info('aspect')
        if ee.data.getInfo(aspect_id):
            ee.data.deleteAsset(aspect_id)

        # Aspect units are degrees where 0=N, 90=E, 180=S, 270=W.
        aspect_img = ee.Terrain.aspect(elev_img)

        task = ee.batch.Export.image.toAsset(
            image=aspect_img,
            description='gridet_aspect_asset',
            assetId=aspect_id,
            dimensions=shape,
            crs=crs,
            crsTransform=transform,
            maxPixels=10000000000,
            # pyramidingPolicy=,
        )
        task.start()


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='Ingest Utah GridET ancillary assets into Earth Engine',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # parser.add_argument(
    #     '--zero', default=False, action='store_true',
    #     help='Set elevation nodata values to 0')
    parser.add_argument(
        '--overwrite', default=False, action='store_true',
        help='Force overwrite of existing files')
    parser.add_argument('--project', help='Earth Engine Project ID')
    parser.add_argument(
        '--debug', default=logging.INFO, const=logging.DEBUG,
        help='Debug level logging', action='store_const', dest='loglevel')

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = arg_parse()
    logging.basicConfig(level=args.loglevel, format='%(message)s')

    main(project_id=args.project, overwrite_flag=args.overwrite)
