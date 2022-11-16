"""
Function for ROMS output files.
"""
import numpy as np
import scipy as sp


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


class GridHandler:
    """
    Provides some method to fix ROMS grid.
    """
    def __init__(self, grid_ds):
        self.grid_ds = grid_ds

    @staticmethod
    def update_u_v_psi_masks_from_rho_mask(grid_ds):
        """
        Rewrite other masks
        Andre matlab code:
        % mask at u, v and psi points
        mask_u = mask_rho(:,1:end-1).*mask_rho(:,2:end);
        mask_v = mask_rho(1:end-1,:).*mask_rho(2:end,:);
        mask_psi = ...
            mask_rho(1:end-1,1:end-1).*mask_rho(1:end-1,2:end).*...
            mask_rho(2:end,1:end-1).*mask_rho(2:end,2:end);
        """
        left_mask = grid_ds.mask_rho.isel(xi_rho=slice(1, None)).values  # left border
        right_mask = grid_ds.mask_rho.isel(xi_rho=slice(None, -1)).values  # Right border
        u_mask = left_mask * right_mask
        grid_ds.mask_u.values = u_mask

        bottom_mask = grid_ds.mask_rho.isel(eta_rho=slice(1, None)).values  # bottom border
        upper_mask = grid_ds.mask_rho.isel(eta_rho=slice(None, -1)).values  # upper border
        v_mask = bottom_mask * upper_mask
        grid_ds.mask_v.values = v_mask

        rho_mask = grid_ds.mask_rho.values
        psi_mask = (
            rho_mask[:-1, :-1] * rho_mask[:-1, 1:] *
            rho_mask[1:, :-1] * rho_mask[1:, 1:]
        )
        grid_ds.mask_psi.values = psi_mask

        return grid_ds

    @staticmethod
    def to_netcdf(grid_ds, name):
        grid_ds.to_netcdf(path=f'fram_data/norfjords_160m_grid_{name}.nc', format='NETCDF4')

    @staticmethod
    def filter_gaussian(grid_ds, sigma=1):
        """
        Applies gaussian filter to ds.h.values
        """
        new_values = sp.ndimage.gaussian_filter(grid_ds.h.values, sigma, mode='nearest')
        grid_ds.h.values = new_values

        return grid_ds

    @staticmethod
    def filter_min_depth(grid_ds, min_depth=10):
        grid_ds['h'] = grid_ds.h.where(grid_ds.h > min_depth, min_depth)

        return grid_ds

    def find_matches(self, border_position: str):
        """
        find indices of matches
        """
        up_eta_rho = self.grid_ds.dims['eta_rho']
        up_xi_rho = self.grid_ds.dims['xi_rho']

        border = {
            'west': self.grid_ds.mask_rho.isel(xi_rho=slice(0, 2)),
            'north': self.grid_ds.mask_rho.isel(eta_rho=slice(up_eta_rho-2, up_eta_rho)),  # grid_ds.mask_rho[-2:, :]
            'east': self.grid_ds.mask_rho.isel(xi_rho=slice(up_xi_rho-2, up_xi_rho)),
            'south': self.grid_ds.mask_rho.isel(eta_rho=slice(0, 2)),
        }

        pattern = {
            'west': np.array([[1., 0.]]),
            'north': np.array([[0.], [1.]]),
            'east': np.array([[0., 1.]]),
            'south': np.array([[1.], [0.]]),
        }

        if border_position in ('north', 'south'):
            idx = np.argwhere(np.asarray(border[border_position].values == pattern[border_position],
                              dtype=np.int32).transpose() @ np.array([1, 1]) > 1)
        elif border_position in ('west', 'east'):
            idx = np.argwhere(np.asarray(border[border_position].values == pattern[border_position],
                              dtype=np.int32) @ np.array([1, 1]) > 1)
        else:
            raise ValueError

        print(f"{border_position} indices: {idx.transpose()[0]}")

        return idx

    def get_values(self, eta_rho, xi_rho):
        try:
            return self.grid_ds.mask_rho.isel(eta_rho=eta_rho, xi_rho=xi_rho).values
        except IndexError:

            return np.NAN

    def print_values(self, west_idx, north_idx, east_idx, south_idx):
        up_eta_rho = self.grid_ds.dims['eta_rho']
        up_xi_rho = self.grid_ds.dims['xi_rho']

        print(f"West values: {self.get_values(eta_rho=west_idx.squeeze(), xi_rho=0)}")
        print(f"North values: {self.get_values(eta_rho=up_eta_rho-1, xi_rho=north_idx.squeeze())}")
        print(f"East values: {self.get_values(eta_rho=east_idx.squeeze(), xi_rho=up_xi_rho-1)}")
        print(f"South values: {self.get_values(eta_rho=0, xi_rho=south_idx.squeeze())}")


if __name__ == "__main__":
    print("This is a module file.")
