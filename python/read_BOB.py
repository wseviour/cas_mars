'''
Script to make a netCDF file from BOB binary output
'''

import numpy as np
import glob
# import matplotlib
# matplotlib.use('qt5agg')
# import matplotlib.pyplot as plt
import xarray as xr
# import matplotlib.path as mpath
import cartopy.crs as ccrs

def ds_from_BOB(run_dir, vars, res, time_step=0.25,units=None):
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
        nlat = int(nlon/2)
    except KeyError:
        print("resolution not found")

    tmp = np.rad2deg(np.genfromtxt(run_dir+'/GRID.T%s' % res, max_rows=nlat/2)[:,1])
    lats = np.append(90-tmp,-90+np.flip(tmp,axis=0))
    #lats = np.rad2deg(np.append(lats[::-1,1],-lats[:,1]))
    lon_spacing = 360./float(nlons[res])
    lons = np.arange(0,360,lon_spacing)
    ds = xr.Dataset()

    for var in vars:
        print('* Reading %s' % var)
        files = sorted(glob.glob(run_dir+'/'+var+'.?????'))[1:]
        a = [np.reshape(np.fromfile(f, dtype='float32'), (nlat,nlon)) for f in files]
        arr = np.stack(a, axis=0)
        time = np.arange(len(files))*time_step
        arr = xr.DataArray(arr,dims=('time','latitude','longitude'),coords=[time,lats,lons])
        if units:
            arr.attrs['units'] = units[var]
        ds[var] = arr

    return ds


if __name__ == '__main__':

    PATH = '../model_output/res_test57.-70.-nu4-urlx-kt4.0.c-0020sat200.0.T85'
    res = 85

    ds = ds_from_BOB(PATH, ['q','u','v','h'], 85, time_step=0.5)

    #ds.to_netcdf(PATH+'/'+PATH.split('/')[-1]+'.nc')

    # ax = plt.axes(projection=ccrs.NorthPolarStereo())
    # ax.set_extent([-180,180,45,90], ccrs.PlateCarree())
    # ax.contourf(ds.coords['longitude'].data,ds.coords['latitude'].data, ds['q'].isel(time=300).data, 21, transform=ccrs.PlateCarree())
    #
    # plt.show()
