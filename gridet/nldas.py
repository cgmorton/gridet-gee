import math

import ee

import gridet.model as model
import gridet.solar as solar
import gridet.utils as utils


def eto(nldas_img):
    """Compute hourly grass reference ET from an hourly NLDAS image

    Parameters
    ----------
    nldas_img : ee.Image

    Returns
    -------
    ee.Image

    """
    return (
        etsz(nldas_img, surface='grass')
        .rename('eto')
        .set('system:time_start', nldas_img.get('system:time_start'))
    )


def etr(nldas_img):
    """Compute hourly alfalfa reference ET from an hourly NLDAS image

    Parameters
    ----------
    nldas_img : ee.Image

    Returns
    -------
    ee.Image

    """
    return (
        etsz(nldas_img, surface='alfalfa')
        .rename('etr')
        .set('system:time_start', nldas_img.get('system:time_start'))
    )


def etsz(nldas_img, surface):
    """Compute hourly reference ET from an hourly NLDAS image

    Parameters
    ----------
    nldas_img : ee.Image
    surface : {'grass', 'alfalfa'}

    Returns
    -------
    ee.Image

    """
    if surface.lower() not in ['grass', 'alfalfa']:
        raise ValueError('surface must be "grass" or "alfalfa"')

    # Apply the GridET adjustments/corrections to the NLDAS hourly image
    # This function will also compute the variables needed the reference ET calculation
    hourly_img = adjust(nldas_img)

    # TODO: Pass in as input parameters
    # Load the GridET ancillary assets (elevation, lat, lon, slope, aspect)
    elevation = ee.Image('projects/openet/assets/reference_et/utah/gridet/ancillary/elevation')

    # CGM - Commenting out for now since GridET elevation assets are already in feet
    # Convert elevations to feet
    # elevation = utils.ToFeet(elevation)

    # NLDAS wind speed is at 10m but is corrected down to 2m when read in below
    # Note, the Reference ET function in this module is expecting the height in feet
    anemometer_height = utils.ToFeet(ee.Number(2))

    # Compute hourly reference ET
    return model.CalculateHourlyASCEReferenceET(
        elevation=elevation,
        temperature=hourly_img.select(['temperature']),
        relative_humidity=hourly_img.select(['relative_humidity']),
        wind_speed=hourly_img.select(['wind_speed']),
        solar_radiation=hourly_img.select(['shortwave_radiation']),
        extraterrestrial_radiation=hourly_img.select(['extraterrestrial_radiation']),
        anemometer_height=anemometer_height,
        reference_type=surface,
        pressure=hourly_img.select(['pressure']),
    ).rename('et_reference').set({'system:time_start': nldas_img.get('system:time_start')})


def adjust(nldas_img):
    """Apply GridET adjustments/corrections to the hourly NLDAS variables

    Parameters
    ----------
    nldas_img : ee.Image

    Returns
    -------
    ee.Image

    """

    record_date = ee.Date(nldas_img.get('system:time_start'))

    nldas_source = (
        ee.Image(nldas_img)
        .select(['temperature', 'specific_humidity', 'pressure', 'wind_u', 'wind_v', 'shortwave_radiation'])
    )

    # Load the GridET NLDAS ancillary assets (note elevation is in feet)
    nldas_elevation = ee.Image('projects/openet/assets/reference_et/utah/gridet/ancillary/nldas2_elevation')
    nldas_aspect = ee.Image('projects/openet/assets/reference_et/utah/gridet/ancillary/nldas2_aspect')
    nldas_slope = ee.Image('projects/openet/assets/reference_et/utah/gridet/ancillary/nldas2_slope')
    nldas_latitude = ee.Image('projects/openet/assets/reference_et/utah/gridet/ancillary/nldas2_latitude')
    nldas_longitude = ee.Image('projects/openet/assets/reference_et/utah/gridet/ancillary/nldas2_longitude')
    # nldas_elevation = ee.Image('projects/openet/assets/meteorology/nldas/ancillary/elevation')
    # nldas_aspect = ee.Image('projects/openet/assets/meteorology/nldas/ancillary/aspect')
    # nldas_slope = ee.Image('projects/openet/assets/meteorology/nldas/ancillary/slope')
    # nldas_latitude = ee.Image('projects/openet/assets/meteorology/nldas/ancillary/latitude')
    # nldas_longitude = ee.Image('projects/openet/assets/meteorology/nldas/ancillary/longitude')

    # TODO: Pass in as input parameters
    # Load the GridET ancillary assets (elevation, lat, lon, slope, aspect)
    elevation = ee.Image('projects/openet/assets/reference_et/utah/gridet/ancillary/elevation')
    aspect = ee.Image('projects/openet/assets/reference_et/utah/gridet/ancillary/aspect')
    slope = ee.Image('projects/openet/assets/reference_et/utah/gridet/ancillary/slope')
    latitude = ee.Image('projects/openet/assets/reference_et/utah/gridet/ancillary/latitude')
    longitude = ee.Image('projects/openet/assets/reference_et/utah/gridet/ancillary/longitude')

    mask_img = ee.Image('projects/openet/assets/reference_et/utah/gridet/ancillary/mask')
    mask_crs = mask_img.projection().wkt()
    mask_transform = ee.List(ee.Dictionary(ee.Algorithms.Describe(mask_img.projection())).get('transform'))

    # CGM - Commenting out for now since GridET elevation assets are already in feet
    # Convert elevations to feet
    # elevation = utils.ToFeet(elevation)
    # nldas_elevation = utils.ToFeet(nldas_elevation)

    # Elevation difference between GridET elevation and interpolated NLDAS elevation
    # CGM - This is slightly different then the approach in GridET
    #   where DeltaZ is computed for the NLDAS grid cells then interpolated
    # CGM - The bilinear resampling was not working correctly in initial testing
    delta_z = (
        elevation
        .subtract(nldas_elevation.resample('bilinear'))
        #.subtract(nldas_elevation.resample('bilinear').reproject(mask_crs, mask_transform))
        .rename('delta_z')
    )

    # GLOBALS
    # Derived from Strong, C., Khatri, K. B., Kochanski, A. K., Lewis, C. S., & Allen, L. N. (2017).
    # Reference evapotranspiration from coarse-scale and dynamically downscaled data in complex terrain:
    # Sensitivity to interpolation and resolution. Journal of Hydrology, 548, 406-418.
    # Fahrenheit/Foot
    # CGM - Modified to add the December value at the beginning and the January value at the end
    #   so the month value can be used directly as the index
    monthly_lapse_rate = ee.List([
        -0.002132,
        -0.001742, -0.002699, -0.003064, -0.003432, -0.003262, -0.00319,
        -0.003046, -0.002941, -0.002659, -0.002622, -0.002247, -0.002132,
        -0.001742
    ])
    AirTemperatureCorrectionCoefficients = [1.58, 0.59, -1.53, -3.73, 1.4, 0.0551]
    HumidityCorrectionCoefficients = [-21.9, 0.78, 3.55, 11.6, -5.05, 0.274]

    # Lapse Rate [Fahrenheit/Foot]
    # Get the bracketing dates (for the 15th of the month)
    # Move the start date forward a month if on or past the 15th
    # Move the end date back a month if the before the 15th
    m1 = record_date.update(day=15).advance(-1, 'month').advance(record_date.get('day').gte(15), 'month')
    m2 = record_date.update(day=15).advance(1, 'month').advance(record_date.get('day').lt(15).multiply(-1), 'month')
    fraction = record_date.difference(m1, 'days').divide(m2.difference(m1, 'days'))
    # Using the month as the index since the table was padded with the Dec and Jan values
    m1_lapse = ee.Number(monthly_lapse_rate.get(m1.get('month')))
    m2_lapse = ee.Number(monthly_lapse_rate.get(m2.get('month')))
    lapse_rate = m2_lapse.subtract(m1_lapse).multiply(fraction).add(m1_lapse)

    # The getRelative() function returns values 0-364, so subtract 1 is not needed
    function_doy = record_date.getRelative('day', 'year').divide(365).multiply(2 * math.pi)
    function_hour = record_date.get('hour').divide(23).multiply(2 * math.pi)
    time_variables = [1, function_doy.cos(), function_doy.sin(), function_hour.cos(), function_hour.sin()]

    # Bilinearly interpolate to the output resolution
    nldas_interp = nldas_source.resample('bilinear')
    #nldas_interp = nldas_source.resample('bilinear').reproject(mask_crs, mask_transform)

    # Correct for Topographical Differences
    # CGM - This is slightly different than GridET approach since
    #   DeltaZ is computed directly at the GridET scale
    #   instead of being interpolated from the NLDAS grid
    # The X * (1 - coeff) is equivalent to (X - X * coeff)
    temperature = (
        utils.ToFahrenheit(nldas_interp.select(['temperature']))
        .add(delta_z.multiply(lapse_rate))
        .multiply(1 - AirTemperatureCorrectionCoefficients[5])
        .subtract(utils.SumProduct(time_variables, AirTemperatureCorrectionCoefficients[:5]))
        # TODO: Check if this is needed and what min and max values are used
        # .clamp(MinimumAirTemperature, MaximumAirTemperature)
    )

    # GridET uses the NLDAS air pressure instead of computing from elevation
    #   so that it can be adjusted for elevation in the same way as temperature
    # CGM - This is slightly different than GridET approach since
    #   DeltaZ is computed directly at the GridET scale
    pressure = model.AdjustAirPressure(
        pressure=nldas_interp.select(['pressure']).divide(1000),
        temperature=temperature,
        delta_z=delta_z,
    )

    # Extraterrestrial Solar Radiation (Langley/hr)
    # Compute Ra at each 15 minute interval bracketing the image and then average
    # First compute Ra at the GridET resolution using interpolated inputs
    project_ra = solar.CalculateInstantaneousRa(
        longitude=longitude,
        latitude=latitude,
        elevation=elevation,
        slope=slope,
        azimuth=aspect,
        temperature=temperature,
        pressure=pressure,
        solar_position=solar.CalculateSolarPosition(record_date),
    )
    for minute_offset in [-30, -15, 15, 30]:
        temp_ra = solar.CalculateInstantaneousRa(
            longitude=longitude,
            latitude=latitude,
            elevation=elevation,
            slope=slope,
            azimuth=aspect,
            temperature=temperature,
            pressure=pressure,
            solar_position=solar.CalculateSolarPosition(record_date.advance(minute_offset, 'minute')),
        )
        project_ra = project_ra.add(temp_ra)
    project_ra = project_ra.divide(5).rename('extraterrestrial_radiation')

    nldas_source_ra = solar.CalculateInstantaneousRa(
        longitude=nldas_longitude,
        latitude=nldas_latitude,
        elevation=nldas_elevation,
        slope=nldas_slope,
        azimuth=nldas_aspect,
        temperature=utils.ToFahrenheit(nldas_source.select(['temperature'])),
        pressure=nldas_source.select(['pressure']).divide(1000),
        solar_position=solar.CalculateSolarPosition(record_date),
    )
    for minute_offset in [-30, -15, 15, 30]:
        temp_ra = solar.CalculateInstantaneousRa(
            longitude=nldas_longitude,
            latitude=nldas_latitude,
            elevation=nldas_elevation,
            slope=nldas_slope,
            azimuth=nldas_aspect,
            temperature=utils.ToFahrenheit(nldas_source.select(['temperature'])),
            pressure=nldas_source.select(['pressure']).divide(1000),
            solar_position=solar.CalculateSolarPosition(record_date.advance(minute_offset, 'minute')),
        )
        nldas_source_ra = nldas_source_ra.add(temp_ra)
    nldas_source_ra = nldas_source_ra.divide(5)

    nldas_source_rs = utils.ToLangleysPerHour(nldas_source.select(['shortwave_radiation']))

    # Bilinearly interpolate NLDAS Ra and Rs to the GridET grid before TranslateRs call
    # CGM - This is slightly different than GridET approach since
    #   TranslateRs is being applied after interpolating the NLDAS Ra and Rs
    nldas_interp_ra = nldas_source_ra.resample('bilinear')
    nldas_interp_rs = nldas_source_rs.resample('bilinear')

    # Solar Radiation (Langley/hr)
    rs = solar.TranslateRs(
        rs_source=nldas_interp_rs,
        ra_source=nldas_interp_ra,
        ra_destination=project_ra,
        slope=slope,
    )

    # Wind speed (mph)
    # In GridET, the wind vectors are bilinearly interpolated to the new grid/point,
    #   then adjusted down to 2m and converted to mph
    # The wind speed is then computed from the vectors and clamped to 5.5 mph (2.45872 m/s)
    u2 = model.CalculateWindSpeedAdjustmentFactor(zw=utils.ToFeet(ee.Number(10)))
    wind_speed = (
        utils.ToMilesPerHour(nldas_interp.select(['wind_u']).multiply(u2)).pow(2)
        .add(utils.ToMilesPerHour(nldas_interp.select(['wind_v']).multiply(u2)).pow(2))
        .sqrt()
        .clamp(0, 5.5)
        .rename('wind_speed')
    )

    # Relative Humidity (%)
    # Compute relative humidity at the coarse scale (before interpolating)
    # CalculateRelativeHumidity is expecting pressure as Pa and temperature as C
    #   so don't convert NLDAS input values here
    # Clamp Relative Humidity to 7-100%
    nldas_source_rh = model.CalculateRelativeHumidity(
        specific_humidity=nldas_source.select(['specific_humidity']),
        pressure=nldas_source.select(['pressure']),
        temperature=nldas_source.select(['temperature']),
    )
    relative_humidity = (
        nldas_source_rh
        .resample('bilinear')
        #.reproject(mask_crs, mask_transform)
        .multiply(1 - HumidityCorrectionCoefficients[5])
        .subtract(utils.SumProduct(time_variables, HumidityCorrectionCoefficients[:5]))
        .clamp(7, 100)
        .rename('relative_humidity')
    )

    return (
        ee.Image([temperature, relative_humidity, wind_speed, rs, project_ra, pressure])
        .set({'system:time_start': record_date.millis()})
    )
