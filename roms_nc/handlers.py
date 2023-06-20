"""
ROMS data processing
"""
from typing import Optional
import cftime

import numpy as np
import scipy as sp
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("darkgrid")


def make_fake_std(ds):
    # variables = ['salt', 'temp', 'u', 'v', 'ubar', 'vbar']
    # ds[variables] = ds[variables].where(ds[variables] == 1, other = 1)
    ds['salt'] = ds['salt'].where(ds['salt'] == 10, other = 10)
    ds['temp'] = ds['temp'].where(ds['temp'] == 10, other = 10)
    ds['u'] = ds['u'].where(ds['u'] == 1, other = 1)
    ds['v'] = ds['v'].where(ds['v'] == 1, other = 1)
    ds['ubar'] = ds['ubar'].where(ds['ubar'] == 1, other = 1)
    ds['vbar'] = ds['vbar'].where(ds['vbar'] == 1, other = 1)
    ds['zeta'] = ds['zeta'].where(ds['zeta'] == 1, other = 1)
    return ds


def extrapolate_fill(ds, parameter_tuple):
    """
    ds: boundary dataset
    parameter_tuple: temp_tuple | salt_tuple
    for example: par_tuple = ('temp_north', 'salt_north', 'temp_south', 'salt_south')
    """
    for name, dim in parameter_tuple:
        ds[name] = ds[name].interpolate_na(dim=dim, method="linear", fill_value="extrapolate")
    return ds


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


def assign_new_time(ds, path):
    time_cftime = cftime.datetime(2017, 1, 15)
    ds = ds.assign_coords(ocean_time=[time_cftime])
    ds.to_netcdf(path=path, format='NETCDF4')


def roho160_tranform_to_z(ds):
    """
    Transforms s coordingate to z with Vtransform = 2
    """
    zo_rho = (ds.hc * ds.s_rho + ds.Cs_rho * ds.h) / (ds.hc + ds.h)
    z_rho = ds.zeta + (ds.zeta + ds.h) * zo_rho
    return z_rho.transpose()


class GridHandler:
    """
    Provides some method to fix ROMS grid.
    """
    def __init__(self, grid_ds):
        self.up_xi = grid_ds.dims['xi_rho'] - 1
        self.up_eta = grid_ds.dims['eta_rho'] - 1
        self.grid_ds = grid_ds

    def get_xi_slice(self, xi, delta):
        return slice(max(min(xi - delta, self.up_xi), 0), max(min(xi + delta + 1, self.up_xi), 0))

    def get_eta_slice(self, eta, delta):
        return slice(max(min(eta - delta, self.up_eta), 0), max(min(eta + delta + 1, self.up_eta), 0))

    def plot_h(self, xi, eta, delta):

        plt.figure(figsize=(5, 5))

        xi_slice = self.get_xi_slice(xi, delta)
        eta_slice = self.get_eta_slice(eta, delta)

        h = self.grid_ds.h.isel(xi_rho=xi_slice, eta_rho=eta_slice) * \
            self.grid_ds.mask_rho.isel(xi_rho=xi_slice, eta_rho=eta_slice)

        plt.subplot(111)
        h.plot()

    def plot_masks(self, xi, eta, delta=10):

        plt.figure(figsize=(15, 4))

        xi_slice = self.get_xi_slice(xi, delta)
        eta_slice = self.get_eta_slice(eta, delta)

        plt.subplot(121)
        self.grid_ds.mask_u.isel(xi_u=xi_slice, eta_u=eta_slice).plot()

        plt.subplot(122)
        self.grid_ds.mask_v.isel(xi_v=xi_slice, eta_v=eta_slice).plot()

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
        new_values = sp.ndimage.gaussian_filter(grid_ds.h.values, sigma, mode='nearest') # type: ignore
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


class RiversHandler:
    """
    Provides rivers forcing related methods.
    """
    def __init__(self, grid_ds):
        self.up_xi = grid_ds.dims['xi_rho'] - 1
        self.up_eta = grid_ds.dims['eta_rho'] - 1
        self.grid_ds = grid_ds

    def get_xi_slice(self, xi, delta):
        return slice(max(min(xi - delta, self.up_xi), 0), max(min(xi + delta + 1, self.up_xi), 0))

    def get_eta_slice(self, eta, delta):
        return slice(max(min(eta - delta, self.up_eta), 0), max(min(eta + delta + 1, self.up_eta), 0))

    def plot_river_estuary_u(self, river, river_xi, river_eta, river_dir, roff):

        assert river_dir == 0
        mask_u = self.grid_ds.mask_u.copy()
        plt.figure(figsize=(15, 4))

        xi_slice = self.get_xi_slice(river_xi, 4)
        eta_slice = self.get_eta_slice(river_eta, 4)

        plt.subplot(121)
        mask_u.isel(
            xi_u=xi_slice,
            eta_u=eta_slice,
            ).plot()

        mask_u[dict(xi_u=river_xi, eta_u=river_eta)] = np.nan

        plt.subplot(122)
        mask_u.isel(
            xi_u=xi_slice,
            eta_u=eta_slice,
            ).plot()
        # to show the fortran locations
        plt.suptitle(f"River: {river}; Xi: {river_xi+1}; Eta: {river_eta}; Direction: {river_dir}, Runoff: {roff}.")

    def plot_river_estuary_v(self, river, river_xi, river_eta, river_dir, roff):

        assert river_dir == 1
        mask_v = self.grid_ds.mask_v.copy()
        plt.figure(figsize=(15, 4))

        xi_slice = self.get_xi_slice(river_xi, 4)
        eta_slice = self.get_eta_slice(river_eta, 4)

        plt.subplot(121)
        mask_v.isel(
            xi_v=xi_slice,
            eta_v=eta_slice,
            ).plot()

        mask_v[dict(xi_v=river_xi, eta_v=river_eta)] = np.nan

        plt.subplot(122)
        mask_v.isel(
            xi_v=xi_slice,
            eta_v=eta_slice,
            ).plot()
        # to show the fortran locations
        plt.suptitle(f"River: {river}; Xi: {river_xi}; Eta: {river_eta+1}; Direction: {river_dir}, Runoff: {roff}.")

    def plot_river_estuary_rho(self, river, river_xi, river_eta, river_dir):

        mask_rho = self.grid_ds.mask_rho.copy()
        plt.figure(figsize=(15, 4))

        plt.subplot(121)
        mask_rho.isel(
            xi_rho=slice(river_xi - 4, river_xi + 5),
            eta_rho=slice(river_eta - 4, river_eta + 5)
            ).plot()

        mask_rho[dict(xi_rho=river_xi, eta_rho=river_eta)] = np.nan

        plt.subplot(122)
        mask_rho.isel(
            xi_rho=slice(river_xi - 4, river_xi + 5),
            eta_rho=slice(river_eta - 4, river_eta + 5)
            ).plot()
        # to show the fortran locations
        plt.suptitle(f"River: {river}; Xi: {river_xi+1}; Eta: {river_eta+1}; Direction: {river_dir}.")

    def plot_rivers_luv(self, river_xi_eta_dir_roff):
        """
        Probably, not sure: Due to roms notation where rho points start from zero,
        buy u and v start from 1, to show them on map, where all point start from 0 (in python),
        we need for rivers in 0 (u) direction subtract 1 from xi coordinate
                              1 (v)                           eta
        https://www.myroms.org/wiki/Grid_Generation
        https://www.myroms.org/wiki/River_Runoff

        Args:
            river_xi_eta_dir: a list of river, xi, eta, dir
        """
        for river, xi, eta, direction, roff in river_xi_eta_dir_roff:
            if direction == 0:
                xi, eta = xi - 1, eta  # pylint: disable=self-assigning-variable
                self.plot_river_estuary_u(river, xi, eta, direction, roff)
            else:
                xi, eta = xi, eta - 1  # pylint: disable=self-assigning-variable
                self.plot_river_estuary_v(river, xi, eta, direction, roff)

    def plot_rivers_lw(self, river_xi_eta_dir):
        for river, xi, eta, direction in river_xi_eta_dir:
            self.plot_river_estuary_rho(river, xi-1, eta-1, direction)

    @staticmethod
    def to_netcdf(ds, name):
        ds.to_netcdf(path=f'fram_data/norfjords_160m_river_{name}.nc', format='NETCDF4')

    @staticmethod
    def add_s_layers(river_ds):
        # Creating new array to fill the gaps
        v_data = np.zeros((5, 175))
        v_add = xr.DataArray(v_data, dims=river_ds['river_Vshape'].dims, coords=river_ds['river_Vshape'].coords)
        temp_data = np.ones((4019, 5, 175)) * 5
        t_add = xr.DataArray(temp_data, dims=river_ds['river_temp'].dims, coords=river_ds['river_temp'].coords)
        salt_data = np.ones((4019, 5, 175))
        s_add = xr.DataArray(salt_data, dims=river_ds['river_salt'].dims, coords=river_ds['river_salt'].coords)

        # New ready data arrays
        river_vshape = xr.concat([v_add, river_ds['river_Vshape']], dim="s_rho")
        river_temp = xr.concat([t_add, river_ds['river_temp']], dim="s_rho")
        river_salt = xr.concat([s_add, river_ds['river_salt']], dim="s_rho")

        # Remove old and add new variables with new s_rho=40 dimension
        river_ds = river_ds.drop_dims('s_rho')
        river_ds = river_ds.assign(river_Vshape=river_vshape, river_temp=river_temp, river_salt=river_salt)

        return river_ds

    @staticmethod
    def stretch_s_layers(river_ds, stretch_to = 10):
        n_rivers = river_ds.dims['river']
        for n_river in range(1, n_rivers+1):
            na_river = river_ds.river_Vshape.sel(river=n_river).values
            n_points = np.count_nonzero(na_river)
            if n_points < stretch_to:
                x = np.arange(0, n_points)
                y = na_river[-n_points:]
                f = sp.interpolate.interp1d(x, y) # type: ignore
                xnew = np.linspace(0, n_points-1, num=stretch_to)
                ynew = f(xnew)
                ynew_unity = ynew / np.sum(ynew)
                river_ds.river_Vshape.loc[-stretch_to:, n_river] = ynew_unity

        return river_ds

    @staticmethod
    def check_coords(ds, xi, eta):
        """
        Locate rivers around the provided coordinates.
        """
        xis = [xi + i for i in range(-10, 11)]
        etas =  [eta + i for i in range(-10, 11)]
        for river_id_orig in ds.river.values:
            river_id = int(river_id_orig-1)
            x_pos = int(ds.river_Xposition[river_id].values)
            y_pos = int(ds.river_Eposition[river_id].values)
            dirct = int(ds.river_direction[river_id].values)
            roff = ds.river_transport.isel(river=river_id, river_time=0).values
            if x_pos in xis and y_pos in etas:
                print(f"River: {int(river_id_orig)}; Xi: {x_pos}; Eta: {y_pos}; Dir: {dirct}; Runoff: {roff}")

    @classmethod
    def check_rivers_runoff(cls, river_ds, verbose=False, runoff=0):
        """
        Returns coordinates and direction for rivers with runoff >= threshold
        By default it will return all rivers (runoff > 0)
        """
        river_xi_eta_dir_roff = []
        for i in range(river_ds.dims['river']):
            river = river_ds.isel(river=i, river_time=0)
            assert i+1 == int(river.river)
            cell_runoff = (abs(river.river_transport) * river.river_Vshape).max().values
            if cell_runoff >= runoff:
                if verbose:
                    print(f"River {i+1} max runoff: {cell_runoff}")
                    print(f"Direction: {river.river_direction} \n")
                    print(f"Coordinates: {river.river_Xposition}; {river.river_Eposition}")
                river_xi_eta_dir_roff.append((
                    int(river.river),
                    int(river.river_Xposition.values),
                    int(river.river_Eposition.values),
                    int(river.river_direction.values),
                    cell_runoff
                ))

        return river_xi_eta_dir_roff

    @classmethod
    def check_rivers_entire_runoff(
        cls,
        river_ds,
        time_slice = slice('2017-01-15','2019-12-31'),
        verbose=False,
        runoff=0
        ):
        """
        Returns coordinates and direction for rivers with runoff >= threshold
        By default it will return all rivers (runoff > 0)
        """
        river_xi_eta_dir_roff = []
        for river_number in river_ds.coords['river']:
            river_max_runoff = abs(river_ds.river_transport.sel(
                river_time=time_slice,
                river=river_number
                )).max().values
            river = river_ds.sel(river=river_number)
            if river_max_runoff >= runoff:
                if verbose:
                    print(f"River {river_number} max runoff: {river_max_runoff}")
                    print(f"Coordinates: {river.river_Xposition}; {river.river_Eposition}")
                    print(f"Direction: {river.river_direction} \n")
                river_xi_eta_dir_roff.append((
                    int(river_number),
                    int(river.river_Xposition.values),
                    int(river.river_Eposition.values),
                    int(river.river_direction.values),
                    int(river_max_runoff)
                ))
        return river_xi_eta_dir_roff

    @classmethod
    def get_bad_rivers(cls, river_ds, bad_rivers_numbers):
        river_xi_eta_dir = RiversHandler.check_rivers_runoff(river_ds, verbose=False)
        bad_rivers = [river_xi_eta_dir[coord] for coord in bad_rivers_numbers]
        print(bad_rivers)
        return bad_rivers

    @staticmethod
    def limit_transport(river_ds, min_val=-10, max_val=10):
        river_ds['river_transport'] = river_ds.river_transport.where(river_ds.river_transport < max_val, max_val)
        river_ds['river_transport'] = river_ds.river_transport.where(river_ds.river_transport > min_val, min_val)

        return river_ds

    @staticmethod
    def gen_rivers_map(river_ds, grid_ds):
        """
        Generates a numpy array with NaN values at the places of rivers.
        It puts them to the rho coordinates, so the coordinates are approximate.
        """
        map_shape = grid_ds.dims['xi_rho'], grid_ds.dims['eta_rho']
        riv_map = np.ones(map_shape)
        for river_id in river_ds.river.values:
            river_id = int(river_id-1)
            x_pos = int(river_ds.river_Xposition[river_id].values)
            y_pos = int(river_ds.river_Eposition[river_id].values)
            riv_map[x_pos, y_pos] = np.NaN

        return riv_map

    @staticmethod
    def mask_all_rivers(river_ds):
        river_ds['river_transport'].values = np.zeros((4019, 175))

        return river_ds

    @staticmethod
    def mask_some_rivers(river_ds, rivers):
        """
        Assignment uses i locations (implicit)
        """

        try:
            river_ds.river_transport[dict(river=rivers)] = np.zeros((4019, 1))
        except ValueError:
            river_ds.river_transport[dict(river=rivers)] = np.zeros((4019))

        return river_ds

    @classmethod
    def split_river(
        cls,
        ds,
        number: float,
        dir_swap: bool = False,
        flux_swap: bool = False,
        f_xi: Optional[float] = None,
        f_eta: Optional[float] = None,
        s_xi: Optional[float] = None,
        s_eta: Optional[float] = None,
        ):
        """
        Decrease flux twice at the current position.
        Adds a new river at the provided position.

        Args:
            ds: Entire dataset to change.
            number: river coordinate number
            dir_swap: river ROMS direction
            flux_swap: multiplies by -1 if True
            f_xi: current river new xi coordinate
            f_eta: current river eta new coordinate
            s_xi: new river xi coordinate
            s_eta: new river eta coordinate

        Returns:
            ds: a new dataset with an additional 'river'.
        """
        directions = {
            0: 1,
            1: 0,
        }

        river_ds = ds.sel(river=number)

        if dir_swap:
            direction = river_ds.river_direction.values
            river_ds.river_direction.values = directions[int(direction)]

        if_swap = -1 if flux_swap else 1
        river_ds.river_transport.values = if_swap * 1/2 * river_ds.river_transport.values

        if f_xi is not None:
            river_ds.river_Xposition.values = f_xi
        if f_eta is not None:
            river_ds.river_Eposition.values = f_eta

        ds.loc[dict(river=number)] = river_ds

        if s_xi is not None:
            river_ds.river_Xposition.values = s_xi
        if s_eta is not None:
            river_ds.river_Eposition.values = s_eta

        river_ds.river.values = ds.river[-1].values + 1
        ds = xr.concat([ds, river_ds], dim="river")

        return ds

    def split_report(
        self,
        river_ds,
        number: float,
        dir_swap: bool = False,
        flux_swap: bool = False,
        f_xi: Optional[float] = None,
        f_eta: Optional[float] = None,
        s_xi: Optional[float] = None,
        s_eta: Optional[float] = None,
        ):

        river_ds = RiversHandler.split_river(river_ds, number,
                                             dir_swap=dir_swap,
                                             flux_swap=flux_swap,
                                             f_xi=f_xi, f_eta=f_eta,
                                             s_xi=s_xi, s_eta=s_eta)
        oldriv = river_ds.sel(river=number)
        oldriv_c = (int(oldriv.river),
                    int(oldriv.river_Xposition.values),
                    int(oldriv.river_Eposition.values),
                    int(oldriv.river_direction.values),
                    oldriv.river_transport.isel(river_time=0).values),
        newriv = river_ds.isel(river=-1)
        newriv_c = (int(newriv.river),
                    int(newriv.river_Xposition.values),
                    int(newriv.river_Eposition.values),
                    int(newriv.river_direction.values),
                    newriv.river_transport.isel(river_time=0).values),
        self.plot_rivers_luv(oldriv_c)
        self.plot_rivers_luv(newriv_c)

        return river_ds


class OutputHandler:
    """
    Provides some method to fix ROMS grid.
    """
    def __init__(self, grid_ds):
        self.up_xi = grid_ds.dims['xi_rho'] - 1
        self.up_eta = grid_ds.dims['eta_rho'] - 1
        self.grid_ds = grid_ds

    def get_xi_slice(self, xi, delta):
        return slice(max(min(xi - delta, self.up_xi), 0), max(min(xi + delta + 1, self.up_xi), 0))

    def get_eta_slice(self, eta, delta):
        return slice(max(min(eta - delta, self.up_eta), 0), max(min(eta + delta + 1, self.up_eta), 0))

    def plot_masks(self, xi, eta, delta=10):

        plt.figure(figsize=(15, 4))

        xi_slice = self.get_xi_slice(xi, delta)
        eta_slice = self.get_eta_slice(eta, delta)

        plt.subplot(121)
        self.grid_ds.mask_u.isel(xi_u=xi_slice, eta_u=eta_slice).plot()

        plt.subplot(122)
        self.grid_ds.mask_v.isel(xi_v=xi_slice, eta_v=eta_slice).plot()

    def plot_u_v(self, ds, xi, eta, s=-1, ocean_time=-1, delta=5):

        plt.figure(figsize=(15, 4))

        xi_slice = self.get_xi_slice(xi, delta)
        eta_slice = self.get_eta_slice(eta, delta)

        u = ds.u.isel(xi_u=xi_slice, eta_u=eta_slice, s_rho=s, ocean_time=ocean_time) #  * \
            # self.grid_ds.mask_u.isel(xi_u=xi_slice, eta_u=eta_slice)
        plt.subplot(121)
        u.plot.pcolormesh(cmap='brg')

        v = ds.v.isel(xi_v=xi_slice, eta_v=eta_slice, s_rho=s, ocean_time=ocean_time) # * \
            # self.grid_ds.mask_v.isel(xi_v=xi_slice, eta_v=eta_slice)
        plt.subplot(122)
        v.plot.pcolormesh(cmap='brg')

    def plot_w(self, ds, xi, eta, s=-1, ocean_time=-1, delta=5):

        plt.figure(figsize=(15, 4))

        xi_slice = self.get_xi_slice(xi, delta)
        eta_slice = self.get_eta_slice(eta, delta)

        w = ds.w.isel(xi_rho=xi_slice, eta_rho=eta_slice, s_w=s, ocean_time=ocean_time) #  * \
            # self.grid_ds.mask_u.isel(xi_u=xi_slice, eta_u=eta_slice)
        plt.subplot(111)
        w.plot.pcolormesh(cmap='brg')

    def plot_ubar_vbar(self, ds, xi, eta, ocean_time=-1, delta=5):

        plt.figure(figsize=(15, 4))

        xi_slice = self.get_xi_slice(xi, delta)
        eta_slice = self.get_eta_slice(eta, delta)

        u = ds.ubar.isel(xi_u=xi_slice, eta_u=eta_slice, ocean_time=ocean_time) #  * \
            # self.grid_ds.mask_u.isel(xi_u=xi_slice, eta_u=eta_slice)
        plt.subplot(121)
        u.plot()

        v = ds.vbar.isel(xi_v=xi_slice, eta_v=eta_slice, ocean_time=ocean_time) # * \
            # self.grid_ds.mask_v.isel(xi_v=xi_slice, eta_v=eta_slice)
        plt.subplot(122)
        v.plot()

    def plot_zeta(self, ds, xi, eta, ocean_time=-1, delta=5):

        plt.figure(figsize=(7, 4))

        xi_slice = self.get_xi_slice(xi, delta)
        eta_slice = self.get_eta_slice(eta, delta)

        (self.grid_ds.mask_rho.isel(xi_rho=xi_slice, eta_rho=eta_slice) *
         ds.zeta.isel(ocean_time=ocean_time,
                      xi_rho=xi_slice, eta_rho=eta_slice)).plot()

    def plot_temp_salt(self, ds, xi, eta, s=-1, ocean_time=-1, delta=5):

        plt.figure(figsize=(15, 4))

        xi_slice = self.get_xi_slice(xi, delta)
        eta_slice = self.get_eta_slice(eta, delta)

        salt = ds.salt.isel(ocean_time=ocean_time,
                     xi_rho=xi_slice,
                     eta_rho=eta_slice,
                     s_rho=s)
        plt.subplot(121)
        salt.plot.pcolormesh(cmap='brg')

        temp = ds.temp.isel(ocean_time=ocean_time,
                     xi_rho=xi_slice,
                     eta_rho=eta_slice,
                     s_rho=s)
        plt.subplot(122)
        temp.plot.pcolormesh(cmap='brg')

    @staticmethod
    def check_zeta(ds):
        """
        da_stacked containes values of zeta < -1 and their coordinates
        """
        da = ds.zeta.isel(ocean_time=-1)
        da_stacked = da.where(da<-1).stack(x=['eta_rho','xi_rho'])
        da_stacked = da_stacked[da_stacked.notnull()]
        return da_stacked

    @staticmethod
    def check_temp(ds):
        """
        Finds lacations (eta, xi, s) of the points with high temperature
        """
        temp_da = ds.temp.isel(ocean_time=-1)
        temp_da_stacked = temp_da.where(temp_da>12).stack(x=['eta_rho', 'xi_rho', 's_rho'])
        temp_da_stacked = temp_da_stacked[temp_da_stacked.notnull()]
        return temp_da_stacked

    @staticmethod
    def check_v(ds):
        """
        Finds lacations (eta, xi, s) of the points with high v
        """
        da = ds.v.isel(ocean_time=-1)
        da_stacked = da.where(da>4).stack(x=['eta_v', 'xi_v', 's_rho'])
        da_stacked = da_stacked[da_stacked.notnull()]
        return da_stacked

    @staticmethod
    def plot_bad_temp(ds, temp_da_point):
        bad_temp_point = temp_da_point
        eta, xi, s = bad_temp_point.eta_rho.values, bad_temp_point.xi_rho.values, bad_temp_point.s_rho.values
        print(f"Xi: {xi}; Eta: {eta}")
        ds.temp.isel(ocean_time=-1,
                     eta_rho=slice(eta-5,eta+5),
                     xi_rho=slice(xi-5,xi+5)).sel(s_rho=s).plot(figsize=(10, 5))

    @staticmethod
    def plot_entire_u(ds):
        ds.u.isel(ocean_time=-1, s_rho=-1).plot()

    @staticmethod
    def plot_entire_v(ds, s=-1, ocean_time=-1):
        ds.v.isel(ocean_time=ocean_time, s_rho=s).plot()

    def plot_entire_masked_zeta(self, ds):
        (self.grid_ds.mask_rho * ds.zeta.isel(ocean_time=-1)).plot(figsize=(14, 7))


if __name__ == "__main__":
    print("This is a module file.")
