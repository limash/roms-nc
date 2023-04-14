"""
Rewrites cice forcing files of choice with a constant value
"""
import os
import glob

import xarray as xr


def rewrite(filepath: str, var_name: str):
    """
    Puts a constant value to the cice forcing dataset
    """
    ds = xr.open_dataset(filepath)
    ds[var_name].values[:, :, :] = 0
    return ds


def main():
    data_folder = '/cluster/projects/nn8103k/A20/FORCING/HINDCAST/v3/ERA5/cice/dump/'
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

    glob_str, var_name, _ = globs_names[1]
    filepaths = sorted(glob.glob(data_folder+glob_str))
    for filepath in filepaths:
        ds_cice = rewrite(filepath, var_name)
        filename = os.path.basename(filepath)
        ds_cice.to_netcdf(cice_folder+filename)
        print(f"{filename} is written.")


if __name__ == "__main__":
    main()
