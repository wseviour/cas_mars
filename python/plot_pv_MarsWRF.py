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


PATH = '../MarsWRF_data/'
file = 'isentropic_surface_data_ls=215_to_ls=240.nc'
cube = iris.load_cube(PATH+file, 'EPV')
L_S = iris.load_cube(PATH+file, 'L_S')
lat = iris.load_cube(PATH+file, 'LAT')
lon = iris.load_cube(PATH+file, 'LON')
theta = iris.load_cube(PATH+file, 'THETA')

lon = lon.data[0,:]
for ilon in range(lon.shape[0]):
    if lon[ilon] < 0:
        lon[ilon] += 360
lon = np.append(lon[90:],lon[:90])

cube = cube.data[:50,:]
cube = np.append(cube[:,:,:,90:],cube[:,:,:,:90], axis=3)

L_S = iris.coords.DimCoord(L_S.data[:50], standard_name='angle_of_incidence', units='degree', var_name='L_S')
lon = iris.coords.DimCoord(lon,standard_name='longitude',units='degrees_E',var_name='lon',circular=True)
lat = iris.coords.DimCoord(lat.data[:,0],standard_name='latitude',units='degrees_N', var_name='lat')
theta = iris.coords.DimCoord(theta.data,standard_name='air_potential_temperature', units='K', var_name='T')

cube = iris.cube.Cube(cube,long_name='Ertel Potential Voriticity',var_name='EPV',
                    dim_coords_and_dims=[(L_S,0),(theta,1),(lat,2),(lon,3)])

level = iris.Constraint(air_potential_temperature = 300)
cube = cube.extract(level)


# #PATH = '/home/local/WIN/wseviou1/data/swvac/'
# PATH = '../output/'
# ext = 'res_test57.-70.-nu4-urlx-kt4.0.c-0020sat200.0.T85'
# var = 'q'
# var_name = 'potential_vorticity_of_atmosphere_layer'
# res = 85
anim=False
plot=True
#
days=[5,10]
# #days=[0,10,20,30,40,49]
# #days=[0,13,15,17,19,21,23,25,100]
# #days=[20,22,24,26,28,30,35,200]
#
# # No need to change below here
# #-------------------------#
#
# files = sorted(glob.glob(PATH+ext+'/'+var+'.?????'))
#
# cube = cube_from_BOB(files, res, var_name, new_format=False)
#
lats = cube.coord('latitude').points
lons = cube.coord('longitude').points
# omega=2*np.pi/86165.0
# f = 2*omega*np.sin(np.deg2rad(lats))


if anim:

    for iday in range(200,220):

        print iday
        fig = plt.figure()
        ax1 = fig.add_subplot(1,1,1, projection=ccrs.NorthPolarStereo())
        ang = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(ang), np.cos(ang)]).T
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
        ang = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(ang), np.cos(ang)]).T
        circle = mpath.Path(verts * radius + center)
        ax.set_boundary(circle, transform=ax.transAxes)
        ax.set_extent([-180, 180,45, 90], ccrs.PlateCarree())
        cyclic_data, cyclic_lons = cartopy.util.add_cyclic_point(cube[day,:,:].data, coord=lons)
        con = ax.contourf(cyclic_lons, lats, cyclic_data,np.linspace(0,0.003,30),cmap='viridis',extend='both', transform=ccrs.PlateCarree())
        for c in con.collections:
            c.set_rasterized(True)
        ax.text(-0.05,0.9, 'sol %s' % day, transform=ax.transAxes, fontsize=9)
        ax.gridlines(color='white',ylocs=[0,20,40,60,80,90])

    plt.tight_layout()
    #plt.savefig('../plots/PV_'+ext+'.pdf')
    plt.show()
