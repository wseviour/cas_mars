"""
Script to plot PV and tracer side by side
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import cartopy.crs as ccrs
import cartopy.util

kt = ['0.0', '1.0', '2.0']
#kt = ['1.0']
tr = 'hat'

for ikt in kt:
    PATH = '../model_output/netcdf/ann57.-70.-nu4-urlx-kt%s-%s.c-0020.T85.nc' % (ikt,tr)

    ds = xr.open_dataset(PATH)

    lats = ds.coords['latitude'].data
    lons = ds.coords['longitude'].data

    days = [90,110,140,190]

    fig = plt.figure(figsize=(6,12))

    for iday in days:
        ax1 = fig.add_subplot(len(days),2,days.index(iday)*2+1, projection=ccrs.NorthPolarStereo())
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax1.set_boundary(circle, transform=ax1.transAxes)
        ax1.set_title('PV: day %03d' % iday)
        ax1.set_extent([-180, 180, 20, 90], ccrs.PlateCarree())
        cyclic_data, cyclic_lons = cartopy.util.add_cyclic_point(ds.q[iday*4,:].data, coord=lons)
        ax1.contourf(cyclic_lons, lats, cyclic_data, np.linspace(0,1.5,30),cmap='Blues',extend='both', transform=ccrs.PlateCarree())

        ax2 = fig.add_subplot(len(days),2,days.index(iday)*2+1+1, projection=ccrs.NorthPolarStereo())
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax2.set_boundary(circle, transform=ax2.transAxes)
        ax2.set_title('Tracer: day %03d' % iday)
        ax2.set_extent([-180, 180, 20, 90], ccrs.PlateCarree())
        cyclic_data, cyclic_lons = cartopy.util.add_cyclic_point(ds.s[iday*4,:].data, coord=lons)
        ax2.contourf(cyclic_lons, lats, cyclic_data,np.linspace(0,0.9,30),cmap='Reds',extend='both', transform=ccrs.PlateCarree())

    plt.tight_layout()
    plt.savefig('../plots/tracer_%s_pv_kt%s.pdf' % (tr,ikt))
    #plt.show()



