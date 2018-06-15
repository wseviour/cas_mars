'''
Script to make PV animation and plot from BOB binary output
'''

import numpy as np
import iris
import glob
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.path as mpath
import iris.plot as iplt
from read_BOB import cube_from_BOB
import cartopy.util

#PATH = '/home/local/WIN/wseviou1/data/swvac/'
PATH = '../model_output/'
ext = 'res_test57.-70.-nu4-urlx-kt4.0.c-0020sat200.0.T85'
var = 'q'
var_name = 'potential_vorticity_of_atmosphere_layer'
res = 85
anim=True
plot=False

days=[160,180,200,220,240,260]
#days=[0,10,20,30,40,49]
#days=[0,13,15,17,19,21,23,25,100]
#days=[20,22,24,26,28,30,35,200]

# No need to change below here
#-------------------------#

files = sorted(glob.glob(PATH+ext+'/'+var+'.?????'))

cube = cube_from_BOB(files, res, var_name, new_format=False)

lats = cube.coord('latitude').points
lons = cube.coord('longitude').points
omega=2*np.pi/86165.0
f = 2*omega*np.sin(np.deg2rad(lats))


if anim:

    for iday in range(200,220):

        print iday
        fig = plt.figure()
        ax1 = fig.add_subplot(1,1,1, projection=ccrs.NorthPolarStereo())
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax1.set_boundary(circle, transform=ax1.transAxes)
        ax1.set_title('day %03d' % iday)
        ax1.set_extent([-180, 180, 40, 90], ccrs.PlateCarree())
        cyclic_data, cyclic_lons = cartopy.util.add_cyclic_point(cube[iday,:,:].data, coord=lons)
        plt.contourf(cyclic_lons, lats, cyclic_data,np.linspace(0,1.5,30),cmap='viridis',extend='both', transform=ccrs.PlateCarree())
        #plt.colorbar()
        plt.savefig('../plots/tmp/day_%03d.png' % iday)
        plt.close()


if plot:

    fig = plt.figure(figsize=(7,4.7))
    for day in days:
        print day
        ax = fig.add_subplot(2,3,days.index(day)+1, projection=ccrs.NorthPolarStereo())
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax.set_boundary(circle, transform=ax.transAxes)
        ax.set_extent([-180, 180,40, 90], ccrs.PlateCarree())
        cyclic_data, cyclic_lons = cartopy.util.add_cyclic_point(cube[day,:,:].data, coord=lons)
        con = ax.contourf(cyclic_lons, lats, cyclic_data,np.linspace(0,1.6,30),cmap='viridis',extend='both', transform=ccrs.PlateCarree())
        for c in con.collections:
            c.set_rasterized(True)
        ax.text(-0.05,0.9, 'sol %s' % day, transform=ax.transAxes, fontsize=9)
        ax.gridlines(color='white',ylocs=[0,20,40,60,80,90])

    plt.tight_layout()
    #plt.savefig('../plots/PV_'+ext+'.pdf')
    plt.show()
