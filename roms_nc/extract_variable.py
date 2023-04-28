"""
Extracts a variable from several netcdf files.
"""
import xarray as xr


def write_variable(var_name: str, dir_name: str, number: int):
    print(f"Start: {var_name}")
    ds = xr.open_mfdataset(dir_name + f'roho160_his_00{number}*.nc')
    ds[var_name].to_netcdf(dir_name + f'roho160_{var_name}_{number}.nc', 'w')
    ds.close()
    print(f"Finish: {var_name}")


def write_variables(var_names: tuple, dir_name: str, number: int):
    ds = xr.open_mfdataset(dir_name + f'roho160_his_00{number}*.nc')
    for var_name in var_names:
        print(f"Start: {var_name}")
        ds[var_name].to_netcdf(dir_name + f'roho160_{var_name}_{number}.nc', 'w')
        print(f"Finish: {var_name}")
    ds.close()


def main():
    variables = ('temp', 'salt')
    wdir = '/cluster/projects/nn9297k/ROHO160+/OutputData/s_layers_25/1_dec2017-dec2018/'
    number = 3
    # [write_variable(variable, wdir, number) for variable in variables]  # pylint: disable=expression-not-assigned
    write_variables(variables, wdir, number)


if __name__ == "__main__":
    main()

