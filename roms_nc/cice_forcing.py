"""
Interpolate forcing data to the roms grid. Cice requires forcing data to match grid.
to run interactively on Fram, for example:
1. Install a conda environment into any folder with enough space and files
   number quota (home directory does not satisfy these requirements)
2. install xarray into a conda environment
3. activate the conda environment
4. salloc --nodes=1 --time=00:30:00 --qos=short --account=nn8103k
5. python cice_forcing.py
"""
import os
import glob

import xarray as xr


def interpolate_cice(
    filepath: str, var_name: str, time_name: str,
    roms_grid_ds: xr.Dataset
    ):
    """
    Interpolates var_name variables for the certain time_name timesteps
    to the coordinates from roms_grid_ds

    Returns:
        A dataset for a grid from roms_grid_ds
    """
    ds = xr.open_dataset(filepath)

    coords = {
        time_name: ds[var_name][time_name],
        'lon': roms_grid_ds.lon_rho,
        'lat': roms_grid_ds.lat_rho,
    }
    return ds[var_name].interp(coords)


def main():
    roms_grid_ds = xr.open_dataset('/cluster/projects/nn8103k/A20/Grid/A20niva_grd_v1.nc')
    data_folder = '/cluster/projects/nn8103k/A20/FORCING/HINDCAST/v3/ERA5/6hourly/'
    cice_folder = '/cluster/projects/nn8103k/A20/FORCING/HINDCAST/v3/ERA5/cice/'
    globs_names = (
        ('era5_q_year_*_roms_halo.nc', 'Qair', 'qair_time'),
        ('era5_tp_year_*_roms_halo.nc', 'rain', 'rain_time'),
        ('era5_msl_year_*_roms_halo.nc', 'Pair', 'pair_time'),
        ('era5_t2m_year_*_roms_halo.nc', 'Tair', 'Tair_time'),
        ('era5_u10_year_*_roms_halo.nc', 'Uwind', 'wind_time'),
        ('era5_v10_year_*_roms_halo.nc', 'Vwind', 'wind_time'),
        ('era5_tcc_year_*_roms_halo.nc', 'cloud', 'cloud_time'),
    )

    # choose a parameter to interpolate
    glob_str, var_name, time_name = globs_names[0]
    filepaths = sorted(glob.glob(data_folder+glob_str))
    for filepath in filepaths:
        ds_cice = interpolate_cice(filepath, var_name, time_name, roms_grid_ds)
        filename = os.path.basename(filepath)
        ds_cice.to_netcdf(cice_folder+filename)
        print(f"{filename} is written.")


if __name__ == "__main__":
    main()

