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

import numpy as np
import xarray as xr

GRID_FOLDER = '/cluster/projects/nn8103k/A20/Grid/A20niva_grd_v1.nc'
DATA_FOLDER = '/cluster/projects/nn8103k/A20/FORCING/HINDCAST/v3/ERA5/6hourly/'
CICE_FOLDER = '/cluster/projects/nn8103k/A20/FORCING/HINDCAST/v3/ERA5/cice/'
CICE_FOLDER_DUMP = '/cluster/projects/nn8103k/A20/FORCING/HINDCAST/v3/ERA5/cice/dump/'
GLOBS_NAMES = (
    ('era5_q_year_*_roms_halo.nc', 'Qair', 'qair_time'),
    ('era5_tp_year_*_roms_halo.nc', 'rain', 'rain_time'),
    ('era5_msl_year_*_roms_halo.nc', 'Pair', 'pair_time'),
    ('era5_t2m_year_*_roms_halo.nc', 'Tair', 'Tair_time'),
    ('era5_u10_year_*_roms_halo.nc', 'Uwind', 'wind_time'),
    ('era5_v10_year_*_roms_halo.nc', 'Vwind', 'wind_time'),
    ('era5_tcc_year_*_roms_halo.nc', 'cloud', 'cloud_time'),
)


class CiceHandler:
    """
    Prepares some cice forcing files.
    """
    def __init__(self, path_roms_grid):
        self.grid_ds = xr.open_dataset(path_roms_grid)

    def interpolate(
            self,
            filepath: str, var_name: str, time_name: str,
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
            'lon': self.grid_ds.lon_rho,
            'lat': self.grid_ds.lat_rho,
        }
        return ds[var_name].interp(coords)

    def rotate(
            self,
            u_filepath, u_name, v_filepath, v_name
            ):
        u_ds = xr.open_dataset(u_filepath)
        v_ds = xr.open_dataset(v_filepath)
        u_east, v_north = u_ds[u_name], v_ds[v_name]

        # angle of rotation from east to x
        alpha = self.grid_ds['angle']
        cos_alpha = np.cos(alpha)
        sin_alpha = np.sin(alpha)

        u_x = u_east * cos_alpha + v_north * sin_alpha
        v_y = v_north * cos_alpha - u_east * sin_alpha

        u_ds[u_name] = u_x
        v_ds[v_name] = v_y
        return u_ds, v_ds


def interpolate_roms_to_cice(parameter: int):
    handler = CiceHandler(GRID_FOLDER)
    glob_str, var_name, time_name = GLOBS_NAMES[parameter]
    filepaths = sorted(glob.glob(DATA_FOLDER+glob_str))
    for filepath in filepaths:
        ds_cice = handler.interpolate(filepath, var_name, time_name)
        filename = os.path.basename(filepath)
        ds_cice.to_netcdf(CICE_FOLDER+filename)
        print(f"{filename} is written.")


def rotate_winds():
    handler = CiceHandler(GRID_FOLDER)
    glob_str_u, var_name_u, _ = GLOBS_NAMES[4]
    glob_str_v, var_name_v, _ = GLOBS_NAMES[5]
    filepaths_u = sorted(glob.glob(CICE_FOLDER_DUMP+glob_str_u))
    filepaths_v = sorted(glob.glob(CICE_FOLDER_DUMP+glob_str_v))
    for filepath_u, filepath_v in zip(filepaths_u, filepaths_v):
        u_ds, v_ds = handler.rotate(filepath_u, var_name_u, filepath_v, var_name_v)
        filename_u = os.path.basename(filepath_u)
        filename_v = os.path.basename(filepath_v)
        u_ds.to_netcdf(CICE_FOLDER+filename_u)
        v_ds.to_netcdf(CICE_FOLDER+filename_v)
        print(f"{filename_u} is written.")
        print(f"{filename_v} is written.")


if __name__ == "__main__":
    rotate_winds()
