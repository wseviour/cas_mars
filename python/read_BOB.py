'''
Script to make a netCDF file from BOB binary output
'''

import numpy as np
import iris
import glob
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.path as mpath
import iris.plot as iplt


def cube_from_BOB(files, res, var_name, new_format=False, time_step=1.,units=None):
    """ Returns a single cube containing output of `files', a
    list of BOB binary files
    """
    # res, nlons
    nlons = {2730 : 8192,
             1365 : 4096,
             682 : 2058,
             341 : 1024,
             170 : 512,
             85 : 256,
             42 : 128,
             21 : 64
             }
    try:
        nlon = nlons[res]
        nlat = nlon/2
    except KeyError:
        print "resolution not found"

    if new_format:
        a = [np.reshape(np.fromfile(f, dtype='float32')[1:-1], (nlat,nlon)) for f in files]
    else:
        a = [np.reshape(np.fromfile(f, dtype='float32'), (nlat,nlon)) for f in files]
    arr = np.stack(a, axis=0)

    lats = np.genfromtxt('../swbob/grids/GRID.T%s' % res, max_rows=nlat/2)
    lats = np.rad2deg(np.append(lats[::-1,1],-lats[:,1]))

    #lats = np.linspace(90,-90,nlat)
    lon_spacing = 360./nlons[res]
    lons = np.arange(0,360,lon_spacing)
    time = np.arange(len(files))/time_step

    lats_coord = iris.coords.DimCoord(lats, standard_name='latitude', units='degrees_N')
    lons_coord = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees_E')
    time_coord = iris.coords.DimCoord(time, standard_name='time', units='days')

    cube = iris.cube.Cube(arr, standard_name=var_name)
    cube.add_dim_coord(time_coord,0)
    cube.add_dim_coord(lons_coord,2)
    cube.add_dim_coord(lats_coord,1)
    if units:
        cube.units = units

    return cube


if __name__ == '__main__':

    PATH = '/home/local/WIN/wseviou1/data/swvac/ann50-nu4-urlx.c00sat50.0.T42'
    var = 'q'
    var_name = 'potential_vorticity_of_atmosphere_layer'
    res = 42
    anim=False
    plot=True

    #days=[6,10,12,15,18,21,23,27,30]
    days=[11,13,15,17,19,21,23,25,100]

    # No need to change below here
    #-------------------------#

    files = sorted(glob.glob(PATH+'/'+var+'.?????'))

    cube = cube_from_BOB(files)

    if anim:

        for iday in range(0,200):

            fig = plt.figure()
            ax1 = fig.add_subplot(1,1,1, projection=ccrs.NorthPolarStereo())
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)
            ax1.set_boundary(circle, transform=ax1.transAxes)
            ax1.set_extent([-180, 180, 50, 90], ccrs.PlateCarree())
            iplt.contourf(cube[iday,:,:],np.linspace(-0.01,0.02,30),cmap='Greys',extend='both')
            plt.colorbar()
            plt.savefig('/home/local/WIN/wseviou1/plots/day%03d.png' % iday)
            plt.close()


    if plot:

        fig = plt.figure(figsize=(8,8))
        for day in days:
            print day
            ax = fig.add_subplot(3,3,days.index(day)+1, projection=ccrs.NorthPolarStereo())
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)
            ax.set_boundary(circle, transform=ax.transAxes)
            ax.set_extent([-180, 180, 50, 90], ccrs.PlateCarree())
            iplt.contourf(cube[day,:,:],np.linspace(-0.01,0.02,30),cmap='Greys',extend='both')

        plt.tight_layout()
