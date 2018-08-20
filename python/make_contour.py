'''
Script for producing contour file for cas
'''

import numpy as np
import glob
from read_BOB import ds_from_BOB
import matplotlib
matplotlib.use('qt5agg')
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
    r = 3389 # Radius of mars in kilometers.
    return c * r

def calc_con_len(con):
    con_len = 0
    for ipt in range(1,con.shape[0]-1):
        con_len += haversine(con[ipt-1,0],con[ipt-1,1],con[ipt,0],con[ipt,1])
    return con_len


# day of simulation for contour
ndays = 50
start_time = 200
PATH = '../model_output/'
ext = 'res_test57.-70.-nu4-urlx-kt4.0.c-0020sat200.0.T85/'
res = 85

ds = ds_from_BOB(PATH+ext, ['q'], 85, time_step=0.5)

day = 200
factor=1
pv_cons = [0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.35,1.3,1.2,1.1]
PATH = '../model_output/'
ext = 'res_test57.-70.-nu4-urlx-kt4.0.c-0020sat200.0.T85'
var = 'q'
var_name = 'potential_vorticity_of_atmosphere_layer'
res = 85


files = sorted(glob.glob(PATH+ext+'/'+var+'.?????'))
cube = cube_from_BOB(files, res, var_name)
lats = cube.coord('latitude').points
lons = cube.coord('longitude').points
count=0

for pv_con in pv_cons:

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
    con1 = ax.contourf(lons, lats, cube[day,:,:].data*factor,np.linspace(0,1.6,11),cmap='viridis', transform=ccrs.PlateCarree())
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
        filename = '../cas/input_contours/pv%s_tstep_%s_%s_inner.in' % (pv_con,day,ext)
    else:
        filename = '../cas/input_contours/pv%s_tstep_%s_%s.in' % (pv_con,day,ext)

    with open(filename, "a") as csvfile:
        csvfile.write("Contour Advection with Surgery\n")
        csvfile.write("PV %s contour\n" % pv_con)
        csvfile.write("\n")
        csvfile.write("%s  24  0.5000000  0.5000000  0.1000000  0.0000000\n" % (ndays))
        csvfile.write("1 %s 0.00000\n" % a.shape[0])
        csvfile.write("%s %d %d 1.00000\n" % (a.shape[0], a[0,0], a[0,1]))

    with open(filename, 'ab') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ')
        for irow in range(a.shape[0]):
            writer.writerow(a[irow,:])

    count +=1
