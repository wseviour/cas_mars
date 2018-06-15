'''
Script to interpolate winds and ouput text file for use in CAS script
'''

import numpy as np
import iris
import glob
from read_BOB import cube_from_BOB
import csv
import os

ndays = 50
start_time = 200
prefix_u = 'u'
prefix_v = 'v'
PATH = '../model_output/'
ext = 'res_test57.-70.-nu4-urlx-kt4.0.c-0020sat200.0.T85/'
res = 85

os.system('mkdir ../cas/winds_'+ext)

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

var = 'u'
u_files = sorted(glob.glob(PATH+ext+'/'+var+'.?????'))
u_cube = cube_from_BOB(u_files, res, 'eastward_wind')
# restrict to NH
u_cube = u_cube[:,:nlat/2,:]

var = 'v'
v_files = sorted(glob.glob(PATH+ext+'/'+var+'.?????'))
v_cube = cube_from_BOB(v_files, res, 'northward_wind')
# restrict to NH
v_cube = v_cube[:,:nlat/2,:]

#nlon=72
#nlat = 46
#remainder =2

nlon = 144
nlat = 92
remainder = 4

new_grid = [('longitude', np.arange(0,360,360./nlon)),
            ('latitude',  np.arange(0, 90, 90./nlat)[::-1])]

v_interp = v_cube.interpolate(new_grid, iris.analysis.Linear())
u_interp = u_cube.interpolate(new_grid, iris.analysis.Linear())

for iday in range(start_time,start_time+ndays):
    with open('../cas/winds_'+ext+prefix_u+'%05d' % iday, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ')
        for ilat in range(nlat):
            for iline in range(nlon/5):
                writer.writerow(u_interp[iday,ilat,iline*5:5+iline*5].data)
            writer.writerow(u_interp[iday,ilat,-remainder:].data)

for iday in range(start_time,start_time+ndays):
    with open('../cas/winds_'+prefix_v+'%05d' % iday, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ')
        for ilat in range(nlat):
            for iline in range(nlon/5):
                writer.writerow(v_interp[iday,ilat,iline*5:5+iline*5].data)
            writer.writerow(v_interp[iday,ilat,-remainder:].data)
