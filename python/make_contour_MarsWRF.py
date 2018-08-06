'''
Script for producing contour file for cas
'''

import numpy as np
import iris
import glob
from read_BOB import cube_from_BOB
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.util
import matplotlib.path as mpath
import csv
import os
from math import radians, cos, sin, asin, sqrt

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    return c * r

def calc_con_len(con):
    con_len = 0
    for ipt in range(1,con.shape[0]-1):
        con_len += haversine(con[ipt-1,0],con[ipt-1,1],con[ipt,0],con[ipt,1])
    return con_len

PATH = '../MarsWRF_data/'
file = 'isentropic_surface_data_ls=240_to_ls=266.nc'
#file = 'isentropic_surface_data_ls=215_to_ls=240.nc'
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

lev = 300
level = iris.Constraint(air_potential_temperature = lev)
cube = cube.extract(level)

pv_cons = [0.0015,0.002,0.0022,0.0025,0.0027,0.003,0.0029,0.0027,0.0025]
#pv_cons = [0.0010,0.0012,0.0015,0.002,0.0022,0.0025,0.0027,0.003]
day=10
ndays=60
factor=1
lats = cube.coord('latitude').points
lons = cube.coord('longitude').points

count=0
for pv_con in pv_cons:

    print pv_con
    inner = False
    if count > 0:
        if pv_cons[count-1] > pv_con:
            inner = True

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection=ccrs.NorthPolarStereo())
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)
    ax.set_extent([-180, 180,20, 90], ccrs.PlateCarree())
    ax.gridlines()
    #cyclic_data, cyclic_lons = cartopy.util.add_cyclic_point(cube[day,:,:].data, coord=lons)
    con1 = ax.contourf(lons, lats, cube[day,:,:].data*factor,np.linspace(0,0.003,30),cmap='viridis', extend='both',transform=ccrs.PlateCarree())
    con = ax.contour(lons, lats, cube[day,:,:].data*factor,[pv_con],colors='k', transform=ccrs.PlateCarree())
    if len(con.allsegs[0]) == 1:
        a = con.allsegs[0][0]
    else:
        lens = np.zeros(len(con.allsegs[0]))
        for icon in range(len(con.allsegs[0])):
            lens[icon] = calc_con_len(con.allsegs[0][icon])
        if inner:
            # 2nd longest contour
            a = con.allsegs[0][np.where(lens == np.sort(lens)[-2])[0][0]]
        else:
            a = con.allsegs[0][np.argmax(lens)]
    plt.plot(a[:,0],a[:,1], transform=ccrs.Geodetic(), color='red')
    plt.tight_layout()
    plt.show()

    if a[:,0][0] > a[:,1][1]:
        a = a[::-1,:]

    try:
        os.remove('pv%s.in' % pv_con)
    except OSError:
        pass

    if inner:
        filename = '../cas/input_contours_MarsWRF/%sK/pv%s_tstep_%s_%s_%sK_inner.in' % (lev,pv_con,day,file.split('.nc')[0],lev)
    else:
        filename = '../cas/input_contours_MarsWRF/%sK/pv%s_tstep_%s_%s_%sK.in' % (lev,pv_con,day,file.split('.nc')[0],lev)

    with open(filename, "a") as csvfile:
        csvfile.write("Contour Advection with Surgery\n")
        csvfile.write("PV %s contour\n" % pv_con)
        csvfile.write("\n")
        csvfile.write("%s  24  0.25000000  0.25000000  0.1000000  0.0000000\n" % (ndays))
        csvfile.write("1 %s 0.00000\n" % a.shape[0])
        csvfile.write("%s %d %d 1.00000\n" % (a.shape[0], a[0,0], a[0,1]))

    with open(filename, 'ab') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ')
        for irow in range(a.shape[0]):
            writer.writerow(a[irow,:])

    count +=1
