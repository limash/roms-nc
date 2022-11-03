"""
Function for ROMS output files.
"""
import numpy as np


def tranform_to_z(ds):
    """
    Transforms s coordinate to depth (z).
    """
    if ds.Vtransform == 1:
        zo_rho = ds.hc * (ds.s_rho - ds.Cs_r) + ds.Cs_r * ds.h
        z_rho = zo_rho + ds.zeta * (1 + zo_rho / ds.h)
    elif ds.Vtransform == 2:
        zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
        z_rho = ds.zeta + (ds.zeta + ds.h) * zo_rho
    else:
        raise ValueError
    return z_rho.transpose()


def roho160_tranform_to_z(ds):
    """
    Transforms s coordingate to z with Vtransform = 2
    """
    zo_rho = (ds.hc * ds.s_rho + ds.Cs_rho * ds.h) / (ds.hc + ds.h)
    z_rho = ds.zeta + (ds.zeta + ds.h) * zo_rho
    return z_rho.transpose()


def check_rivers_runoff(ds, verbose=False, runoff=0):
    """
    Returns coordinates and direction for
    rivers with runoff >= threshold
    """
    river_xi_eta_dir = []
    for i in range(ds.dims['river']):
        river = ds.isel(river=i, river_time=0)
        assert i+1 == int(river.river)
        cell_runoff = (abs(river.river_transport) * river.river_Vshape).max().values
        if cell_runoff >= runoff:
            if verbose:
                print(f"River {i+1} max runoff: {cell_runoff}")
                print(f"Coordinates: {river.river_Xposition}; {river.river_Eposition}")
                print(f"Direction: {river.river_direction} \n")
            river_xi_eta_dir.append((
                int(river.river),
                int(river.river_Xposition.values),
                int(river.river_Eposition.values),
                int(river.river_direction.values)
            ))
    return river_xi_eta_dir


def gen_rivers_map(ds, grid_ds):
    """
    Generates a numpy array with NaN values at the places of rivers.
    """
    map_shape = grid_ds.dims['xi_rho'], grid_ds.dims['eta_rho']
    riv_map = np.ones(map_shape)
    for river_id in ds.river.values:
        river_id = int(river_id-1)
        x_pos = int(ds.river_Xposition[river_id].values)
        y_pos = int(ds.river_Eposition[river_id].values)
        riv_map[x_pos, y_pos] = np.NaN

    return riv_map


if __name__ == "__main__":
    print("This is a module file.")
