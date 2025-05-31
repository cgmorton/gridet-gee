# gridet-gee

Google Earth Engine implementation of the GridET model (https://github.com/claytonscottlewis/GridET) reference ET calculations.

## Supported Calculations and Datasets

Currently, there is only support for computing hourly grass and alfalfa reference ET from NLDAS hourly images using the "eto" and "etr" functions respectively.  These functions take a single NLDAS image as input and can be mapped over the NLDAS image collection (see the example below) or called directly.  Support for other meteorology image collections or reference ET calculations may be added in the future.

## Example

Example of computing NLDAS based grass reference ET (ETo) for a single day (0-23 UTC) by mapping the "gridet.nldas.eto" function over the NLDAS image collection.

```
import gridet.nldas

eto_img = (
    ee.ImageCollection('NASA/NLDAS/FORA0125_H002')
    .filterDate('2017-07-01', '2027-07-02')
    .map(gridet.nldas.eto)
    .sum()
    .rename(['eto'])
    .set({'system:time_start': ee.Date('2017-07-01').millis()})
)
``` 

## Ancillary Assets

The GEE GridET code is currently hardcoded to be run for the same spatial resolution and study area as the original GridET model using the following ancillary assets.  These data have a spatial resolution of 1760 ft (~536 m), are in a NAD83 Lambert Conformal Conic projection, and were generated from the USGS NED 10m dataset.

| Type      | Asset ID                                                            |
|-----------|---------------------------------------------------------------------|
| elevation | projects/openet/assets/reference_et/utah/gridet/ancillary/elevation |
| aspect    | projects/openet/assets/reference_et/utah/gridet/ancillary/aspect    |
| slope     | projects/openet/assets/reference_et/utah/gridet/ancillary/slope     |
| latitude  | projects/openet/assets/reference_et/utah/gridet/ancillary/latitude  |
| longitude | projects/openet/assets/reference_et/utah/gridet/ancillary/longitude |

NLDAS2 Ancillary Assets

| Type      | Asset ID                                                                   |
|-----------|----------------------------------------------------------------------------|
| elevation | projects/openet/assets/reference_et/utah/gridet/ancillary/nldas2_elevation |
| aspect    | projects/openet/assets/reference_et/utah/gridet/ancillary/nldas2_aspect    |
| slope     | projects/openet/assets/reference_et/utah/gridet/ancillary/nldas2_slope     |
| latitude  | projects/openet/assets/reference_et/utah/gridet/ancillary/nldas2_latitude  |
| longitude | projects/openet/assets/reference_et/utah/gridet/ancillary/nldas2_longitude |

## GEE Implementation Differences

The code is currently separated into four modules (model, nldas, solar, and utils), but this was mostly done to simplify testing and may be consolidated or further split in the future.

Most of the internal variables were renamed to follow the Python PEP8 convention of using all lower case for variable.  Variables that used greek letter symbols were spelled out.  The function names are still the original camel case, but will likely be modified in the future.  The function doc strings have not been modified, but will likely be reformated in the future.  As much as possible, the original comments have been left in place.   

The GEE implementation has a few differences with the original VB.net version that are all some variation of not being able to mimic the exact flow of calculation at the NLDAS scale, interpolation to the GridET scale, and then additional calculation.  In initial testing these code differences do not seem to significantly change the output values, but additional testing is planned.

* TranslateRs is being applied after spatially interpolating the NLDAS grid cell Ra and Rs to the GridET grid.
* Air temperature is computed directly at the GridET scale using the spatially interpolated air temperature and delta Z, instead of computed for each NLDAS grid cell and then interpolating.
* Pressure is computed directly at the GridET scale from the spatially interpolated delta Z, instead of being computed for each NLDAS grid cell and then interpolating.

##  References

Lewis, C. & Allen, L. (2017). Potential crop evapotranspiration and surface evaporation estimates via a gridded weather forcing dataset. Journal of Hydrology, 546, 450-463. https://doi.org/10.1016/j.jhydrol.2016.11.055

Lewis, C., Geli, H., Neale, C. (2014). Comparison of the NLDAS Weather Forcing Model to Agrometeorological Measurements in the western United States. Journal of Hydrology, 510 385-392. https://doi.org/10.1016/j.jhydrol.2013.12.040
