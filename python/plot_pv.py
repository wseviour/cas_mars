'''
Script to make PV animation and plot from BOB binary output
'''

import numpy as np
import matplotlib
matplotlib.use('qt5agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.path as mpath
from read_BOB import ds_from_BOB
import cartopy.util


PATH = '../model_output/'
ext = 'res_test57.-70.-nu4-urlx-kt4.0.c-0020sat200.0.T85'

ds = ds_from_BOB(PATH+ext, ['q','u','v','h'], 85, time_step=0.5)

days=[160,180,200,220,240,260]
anim=False
plot=True
# No need to change below here
#-------------------------#


lats = ds.coords['latitude'].data
lons = ds.coords['longitude'].data
omega=2*np.pi/86165.0
f = 2*omega*np.sin(np.deg2rad(lats))


if anim:

    for iday in range(200,220):

        prin(iday)
        fig = plt.figure()
        ax1 = fig.add_subplot(1,1,1, projection=ccrs.NorthPolarStereo())
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax1.set_boundary(circle, transform=ax1.transAxes)
        ax1.set_title('day %03d' % iday)
        ax1.set_extent([-180, 180, 40, 90], ccrs.PlateCarree())
        plt.contourf(lons, lats, ds['q'][iday,:].data,np.linspace(0,1.5,30),cmap='viridis',extend='both', transform=ccrs.PlateCarree())
        #plt.colorbar()
        plt.savefig('../plots/tmp/day_%03d.png' % iday)
        plt.close()


if plot:

    fig = plt.figure(figsize=(7,4.7))
    for day in days:
        print(day)
        ax = fig.add_subplot(2,3,days.index(day)+1, projection=ccrs.NorthPolarStereo())
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax.set_boundary(circle, transform=ax.transAxes)
        ax.set_extent([-180, 180,40, 90], ccrs.PlateCarree())
        con = ax.contourf(lons, lats, ds['q'][day,:].data,np.linspace(0,1.6,30),cmap='viridis',extend='both', transform=ccrs.PlateCarree())
        for c in con.collections:
            c.set_rasterized(True)
        ax.text(-0.05,0.9, 'sol %s' % day, transform=ax.transAxes, fontsize=9)
        ax.gridlines(color='white',ylocs=[0,20,40,60,80,90])

    plt.tight_layout()
    #plt.savefig('../plots/PV_'+ext+'.pdf')
    plt.show()
