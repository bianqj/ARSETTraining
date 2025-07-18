import xarray as xr
import requests
import netCDF4
import boto3
from botocore import UNSIGNED
from botocore.config import Config
from requests.exceptions import ChunkedEncodingError, ConnectionError
import time
import os
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import warnings
import ssl
import urllib3

warnings.filterwarnings('ignore')
# Disable SSL warnings
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


def get_s3_keys(bucket, s3_client, prefix=''):
    """
    Generate the keys in an S3 bucket.

    :param bucket: Name of the S3 bucket.
    :param prefix: Only fetch keys that start with this prefix (optional).
    """
    kwargs = {'Bucket': bucket}

    if isinstance(prefix, str):
        kwargs['Prefix'] = prefix

    while True:
        resp = s3_client.list_objects_v2(**kwargs)
        for obj in resp['Contents']:
            key = obj['Key']
            if key.startswith(prefix):
                yield key

        try:
            kwargs['ContinuationToken'] = resp['NextContinuationToken']
        except KeyError:
            break


def plot_goes_aod(aod_filename, latlon_filename, output_path, map_extent=None,
                  clim_range=(0, 1), title='GOES Aerosol Optical Depth'):
    """
    Plot GOES AOD data with geographic features.

    Parameters:
    -----------
    aod_filename : str
        Path to the AOD NetCDF file
    latlon_filename : str
        Path to the lat/lon NetCDF file
    output_path : str
        Path where the plot will be saved
    map_extent : list, optional
        Map extent as [lon_min, lon_max, lat_min, lat_max]
        Default: [-125, -66, 24, 50] (CONUS)
    clim_range : tuple, optional
        Color scale limits (min, max). Default: (0, 1)
    title : str, optional
        Plot title. Default: 'GOES Aerosol Optical Depth'

    Returns:
    --------
    bool : True if successful, False otherwise
    """
    try:
        # Clear any existing plots
        plt.close('all')

        # Set default map extent if not provided
        if map_extent is None:
            map_extent = [-125, -66, 24, 50]  # CONUS extent

        # Read AOD data
        with netCDF4.Dataset(aod_filename, 'r') as ds:
            AOD_raw = ds.variables['AOD'][:]
            print(f"AOD data shape: {AOD_raw.shape}")
            print(f"AOD data type: {type(AOD_raw)}")

        # Read lat/lon data
        with netCDF4.Dataset(latlon_filename, 'r') as ds_latlon:
            lat_AOD_raw = ds_latlon.variables['latitude'][:]
            lon_AOD_raw = ds_latlon.variables['longitude'][:]
            print(f"Lat/Lon data shapes: {lat_AOD_raw.shape}, {lon_AOD_raw.shape}")

        # Create the plot
        fig = plt.figure(figsize=(14, 10))
        ax = plt.axes(projection=ccrs.Mercator())
        ax.set_extent(map_extent, crs=ccrs.PlateCarree())

        # Plot the data
        im = ax.pcolormesh(lon_AOD_raw, lat_AOD_raw, AOD_raw,
                           transform=ccrs.PlateCarree(),
                           shading='auto',
                           cmap='viridis')

        # Add geographic features
        ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
        ax.add_feature(cfeature.BORDERS, linewidth=0.8)
        ax.add_feature(cfeature.STATES, linewidth=0.5, edgecolor='black', facecolor='none')
        ax.add_feature(cfeature.LAND, alpha=0.1)
        ax.add_feature(cfeature.OCEAN, alpha=0.1)

        # Add gridlines
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.5, color='gray', alpha=0.7, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 16}
        gl.ylabel_style = {'size': 16}

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, orientation='vertical', pad=0.02, shrink=0.8)
        cbar.set_label('Aerosol Optical Depth (AOD)', fontsize=16)
        cbar.ax.tick_params(labelsize=14)

        # Set color limits and title
        im.set_clim(clim_range[0], clim_range[1])
        plt.title(title, fontsize=18, pad=20)
        plt.tight_layout()

        # Save the plot
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"Plot saved successfully to: {output_path}")
        return True

    except Exception as e:
        print(f"Error creating plot: {e}")
        return False


def download_and_process_goes_data(bucket_name='noaa-goes16', product_name='ABI-L2-AODF',
                                   start_year=2020, end_year=2022, start_day=1, end_day=365,
                                   start_hour=0, end_hour=23, output_base_dir='training_sample\\',
                                   apply_subset=False, subset_coords=None, create_plots=False,
                                   latlon_file=None, plot_output_dir=None):
    """
    Download GOES data and optionally create plots.

    Parameters:
    -----------
    bucket_name : str
        S3 bucket name. Default: 'noaa-goes16'
    product_name : str
        GOES product name. Default: 'ABI-L2-AODF'
    start_year : int
        Starting year for download. Default: 2020
    end_year : int
        Ending year for download (exclusive). Default: 2022
    start_day : int
        Starting day of year (1-366). Default: 1
    end_day : int
        Ending day of year (1-366, inclusive). Default: 365
    start_hour : int
        Starting hour (0-23). Default: 0
    end_hour : int
        Ending hour (0-23, inclusive). Default: 23
    output_base_dir : str
        Base directory for downloaded files
    apply_subset : bool
        Whether to apply spatial subsetting. Default: False (saves full file)
    subset_coords : dict, optional
        Subset coordinates with keys 'x' and 'y', each containing slice objects
        Only used if apply_subset=True
        Examples:
        - GOES East: {'x': slice(986, 1161), 'y': slice(1007, 1108)}
        - GOES West: {'x': slice(3475, 3676), 'y': slice(963, 1052)}
    create_plots : bool
        Whether to create plots for downloaded files. Default: False
    latlon_file : str, optional
        Path to lat/lon file (required if create_plots=True)
    plot_output_dir : str, optional
        Directory to save plots (required if create_plots=True)
    """

    # Set default subset coordinates if subsetting is enabled but no coords provided
    if apply_subset and subset_coords is None:
        subset_coords = {'x': slice(986, 1161), 'y': slice(1007, 1108)}  # Default GOES East
        print("Using default GOES East subset coordinates")

    # Initialize s3 client with SSL verification disabled
    s3_client = boto3.client('s3',
                             config=Config(signature_version=UNSIGNED),
                             use_ssl=False)  # Disable SSL to avoid certificate issues

    # Create output directories
    for year in range(start_year, end_year):
        year_dir = os.path.join(output_base_dir, str(year))
        os.makedirs(year_dir, exist_ok=True)

    if create_plots and plot_output_dir:
        os.makedirs(plot_output_dir, exist_ok=True)

    subset_msg = "with subsetting" if apply_subset else "full files"
    print(
        f"Downloading {subset_msg} for years {start_year}-{end_year - 1}, days {start_day}-{end_day}, hours {start_hour}-{end_hour}")

    for year in range(start_year, end_year):
        # Handle leap years for day range validation
        max_day = 366 if year % 4 == 0 and (year % 100 != 0 or year % 400 == 0) else 365
        actual_end_day = min(end_day, max_day)

        for day_of_year in range(start_day, actual_end_day + 1):
            for hour in range(start_hour, end_hour + 1):
                keys1 = get_s3_keys(bucket_name,
                                    s3_client,
                                    prefix=f'{product_name}/{year}/{day_of_year:03.0f}/{hour:02.0f}/OR_{product_name}-M6')
                try:
                    for key1 in keys1:
                        retries = 10
                        for i in range(retries):
                            try:
                                # Use HTTP instead of HTTPS and disable SSL verification
                                resp1 = requests.get(f'http://{bucket_name}.s3.amazonaws.com/{key1}',
                                                     timeout=30, verify=False)
                                file_name1 = key1.split('/')[-1].split('.')[0]

                                try:
                                    nc4_ds1 = netCDF4.Dataset(file_name1, memory=resp1.content)
                                    store1 = xr.backends.NetCDF4DataStore(nc4_ds1)
                                    DS1 = xr.open_dataset(store1)

                                    # Apply subset only if requested
                                    if apply_subset and subset_coords:
                                        DS_processed = DS1.isel(x=subset_coords['x'], y=subset_coords['y'])
                                        filename_prefix = 'subset_'
                                        print(f"Applied subset: x{subset_coords['x']}, y{subset_coords['y']}")
                                    else:
                                        DS_processed = DS1
                                        filename_prefix = 'full_'
                                        print("Saving full file (no subset applied)")

                                    # Save the file
                                    output_file = os.path.join(output_base_dir, str(year),
                                                               f'{filename_prefix}{file_name1[-72:]}.nc')
                                    DS_processed.to_netcdf(output_file)
                                    print(f"Saved: {output_file}")

                                    # Create plot if requested
                                    if create_plots and latlon_file and plot_output_dir:
                                        plot_filename = os.path.join(plot_output_dir, f'plot_{file_name1[-72:]}.png')
                                        plot_success = plot_goes_aod(
                                            aod_filename=output_file,
                                            latlon_filename=latlon_file,
                                            output_path=plot_filename,
                                            title=f'GOES AOD - {file_name1[-72:]}'
                                        )
                                        if plot_success:
                                            print(f"Plot saved: {plot_filename}")

                                    break

                                except Exception as e:
                                    print(f'Error processing file {file_name1}: {e}')
                                    continue

                            except (ChunkedEncodingError, ConnectionError) as e:
                                print(f"Attempt {i + 1} failed with error: {e}")
                                time.sleep(2)
                            except requests.exceptions.ReadTimeout:
                                print(f"Read timeout occurred. Attempt {i + 1} of {retries}. Retrying in 2 seconds...")
                                time.sleep(2)

                except KeyError:
                    print(
                        f'{product_name}/{year}/{day_of_year:03.0f}/{hour:02.0f}/OR_{product_name}-M6 not found in the resp dictionary.')
                    continue


if __name__ == "__main__":
    # Example usage 1: Download full GOES-16 files (default behavior)
    # print("Starting GOES-16 full file download...")
    # download_and_process_goes_data(
    #     bucket_name='noaa-goes16',
    #     start_year=2020,
    #     end_year=2021,
    #     start_day=1,
    #     end_day=2,  # Just 2 days for testing
    #     start_hour=12,
    #     end_hour=13,
    #     output_base_dir='training_sample\\goes16_full\\',
    #     apply_subset=False,  # Save full files (default)
    #     create_plots=False
    # )

    # Example usage 2: Download with regional subset (smaller files)
    """
    print("Starting GOES-16 subset download...")
    download_and_process_goes_data(
        bucket_name='noaa-goes16',
        start_year=2020,
        end_year=2021,
        start_day=1,
        end_day=5,
        start_hour=15,
        end_hour=17,
        output_base_dir='training_sample\\goes16_subset\\',
        apply_subset=True,  # Apply subsetting
        subset_coords={'x': slice(986, 1161), 'y': slice(1007, 1108)},  # GOES East region
        create_plots=False
    )
    """

    # Example usage 3: Download GOES-17 West Coast subset
    """
    print("Starting GOES-17 West Coast subset download...")
    download_and_process_goes_data(
        bucket_name='noaa-goes17',
        start_year=2020,
        end_year=2021,
        start_day=1,
        end_day=3,
        start_hour=18,
        end_hour=20,
        output_base_dir='training_sample\\goes17_west\\',
        apply_subset=True,
        subset_coords={'x': slice(3475, 3676), 'y': slice(963, 1052)},  # GOES West region
        create_plots=False
    )
    """

    # Example usage 4: Full files with plotting

    print("Starting full file download with plotting...")
    download_and_process_goes_data(
        bucket_name='noaa-goes18',
        start_year=2025,
        end_year=2026,
        start_day=11,
        end_day=11,     # Just one day
        start_hour=15,
        end_hour=18,    # Just one hour
        output_base_dir='training_sample\\goes18_full\\',
        apply_subset=False,  # Full files
        create_plots=True,
        latlon_file='goes18_abi_full_disk_lat_lon.nc', # available at https://www.star.nesdis.noaa.gov/atmospheric-composition-training/satellite_data_goes_imager_projection.php#lat_lon_files
        plot_output_dir='plots\\goes18\\'
    )


    # Example usage 5: Comparison of full vs subset files
    """
    # Download same data both ways for size comparison
    datasets = [
        {'name': 'full', 'subset': False, 'coords': None},
        {'name': 'regional', 'subset': True, 'coords': {'x': slice(986, 1161), 'y': slice(1007, 1108)}}
    ]

    for dataset in datasets:
        print(f"Downloading {dataset['name']} dataset...")
        download_and_process_goes_data(
            bucket_name='noaa-goes16',
            start_year=2020,
            end_year=2021,
            start_day=1,
            end_day=1,  # Just one day
            start_hour=15,
            end_hour=15,  # Just one hour
            output_base_dir=f'comparison\\{dataset["name"]}\\',
            apply_subset=dataset['subset'],
            subset_coords=dataset['coords'],
            create_plots=False
        )
    """
