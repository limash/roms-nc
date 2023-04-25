"""
Extracts a variable from several netcdf files.
"""
import xarray as xr


def write_variable(var_name: str, dir_name: str):
    print(f"Start: {var_name}")
    ds = xr.open_mfdataset(dir_name + 'roho160_avg_000*.nc')
    ds[var_name].to_netcdf(dir_name + f'roho160_{var_name}.nc', 'w')
    ds.close()
    print(f"Finish: {var_name}")


def main():
    wdir = '/cluster/projects/nn9297k/ROHO160+/OutputData/s_layers_25/2_dec2018-sep2019/'
    variables = 'temp', 'salt'
    [write_variable(variable, wdir) for variable in variables]


if __name__ == "__main__":
    main()

