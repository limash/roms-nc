import glob

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from PIL import Image


def plot_temp_salt_contour(time_str, temp, salt, save_dir=None, counter=None):

    temp_levs = np.linspace(-4, 15, 60)
    salt_levs = np.linspace(0, 38, 100)

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(14, 7))
    cf_temp = axs[0].contourf(temp[:,:], cmap='coolwarm', levels=temp_levs)
    cf_salt = axs[1].contourf(salt[:,:], cmap='PuBuGn', levels=salt_levs)

    names = [f'Temperature, {time_str}', 'Salinity']
    cf_top = [cf_temp, cf_salt]
    for i, ax in enumerate(axs.ravel()):
        ax.set_title(names[i])
    for cf, ax in zip(cf_top, axs):
        fig.colorbar(cf, ax=ax, location='bottom', pad=0.05)
    fig.tight_layout()

    if save_dir is not None:
        plt.savefig(save_dir + f'temp_salt_contours/{counter:03}', bbox_inches='tight', dpi=150)
        plt.close()
        print(f"Counter: {counter}")


def make_gif(frame_folder, gifname):
    frames = [Image.open(image) for image in sorted(glob.glob(f"{frame_folder}/*.png"))]
    frame_one = frames[0]
    frame_one.save(frame_folder + gifname, format="GIF", append_images=frames, save_all=True, duration=200, loop=0)


def save_plots():
    folder = '/cluster/projects/nn9297k/ROHO160+/OutputData/s_layers_25/1_dec2017-dec2018/'
    ds_temp = xr.open_dataset(folder + 'roho160_temp_0.nc')
    ds_salt = xr.open_dataset(folder + 'roho160_salt_0.nc')

    for i in range(ds_temp['ocean_time'].shape[0]):
        plot_temp_salt_contour(
            time_str=str(ds_temp['ocean_time'].isel(ocean_time=i).values),
            temp=ds_temp['temp'].isel(ocean_time=i, s_rho=-1).values,
            salt=ds_salt['salt'].isel(ocean_time=i, s_rho=-1).values,
            save_dir=folder,
            counter=i
            )

def save_gif():
    dir_name = '/cluster/projects/nn9297k/ROHO160+/OutputData/s_layers_25/1_dec2017-dec2018/temp_salt_contours'
    gif_name = "/temp_salt.gif"
    make_gif(dir_name, gif_name)


if __name__ == "__main__":
    save_plots()
    save_gif()

