'''
Script to interpolate winds and ouput text file for use in CAS script
'''

import numpy as np
import iris
import glob
from read_BOB import cube_from_BOB
import csv
import os

ndays = 90
start_time=10
theta_lev = 300

PATH = '../MarsWRF_data/'
file = 'isentropic_surface_data_ls=240_to_ls=266.nc'
wind_file = '../MarsWRF_data/MarsWRF_winds_%sK/%s' % (theta_lev,file.split('.nc')[0])
if not os.path.isdir(wind_file):
    os.system("mkdir -p %s" % wind_file)

prefix_u = wind_file+'/u'
prefix_v = wind_file+'/v'

u = iris.load_cube(PATH+file, 'U')
v = iris.load_cube(PATH+file, 'V')
L_S = iris.load_cube(PATH+file, 'L_S')
lat = iris.load_cube(PATH+file, 'LAT')
lon = iris.load_cube(PATH+file, 'LON')
theta = iris.load_cube(PATH+file, 'THETA')

lon = lon.data[0,:]
for ilon in range(lon.shape[0]):
    if lon[ilon] < 0:
        lon[ilon] += 360
lon = np.append(lon[90:],lon[:90])

# Regrid from -180-180 to 0-360
u = u.data[:100,:]
u = np.append(u[:,:,:,90:],u[:,:,:,:90], axis=3)
v = v.data[:100,:]
v = np.append(v[:,:,:,90:],v[:,:,:,:90], axis=3)

L_S = iris.coords.DimCoord(L_S.data[:100], standard_name='angle_of_incidence', units='degree', var_name='L_S')
lon = iris.coords.DimCoord(lon,standard_name='longitude',units='degrees_E',var_name='lon',circular=True)
lat = iris.coords.DimCoord(lat.data[:,0],standard_name='latitude',units='degrees_N', var_name='lat')
theta = iris.coords.DimCoord(theta.data,standard_name='air_potential_temperature', units='K', var_name='T')

u_cube = iris.cube.Cube(u,long_name='zonal wind',var_name='U',
                    dim_coords_and_dims=[(L_S,0),(theta,1),(lat,2),(lon,3)])
v_cube = iris.cube.Cube(v,long_name='meridional wind',var_name='V',
                    dim_coords_and_dims=[(L_S,0),(theta,1),(lat,2),(lon,3)])

level = iris.Constraint(air_potential_temperature = theta_lev)
u_cube = u_cube.extract(level)
v_cube = v_cube.extract(level)

nlon = 144
nlat = 92
remainder = 4

new_grid = [('longitude', np.arange(0,360,360./nlon)),
            ('latitude',  np.arange(0, 90, 90./nlat)[::-1])]

v_interp = v_cube.interpolate(new_grid, iris.analysis.Linear())
u_interp = u_cube.interpolate(new_grid, iris.analysis.Linear())

for iday in range(start_time,start_time+ndays):
    with open(prefix_u+'%05d' % iday, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ')
        for ilat in range(nlat):
            for iline in range(nlon/5):
                writer.writerow(u_interp[iday,ilat,iline*5:5+iline*5].data)
            writer.writerow(u_interp[iday,ilat,-remainder:].data)

for iday in range(start_time,start_time+ndays):
    with open(prefix_v+'%05d' % iday, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ')
        for ilat in range(nlat):
            for iline in range(nlon/5):
                writer.writerow(v_interp[iday,ilat,iline*5:5+iline*5].data)
            writer.writerow(v_interp[iday,ilat,-remainder:].data)
