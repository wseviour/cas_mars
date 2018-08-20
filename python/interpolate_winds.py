'''
Script to interpolate winds and ouput text file for use in CAS script
'''

import numpy as np
import glob
from read_BOB import ds_from_BOB
import csv
import os

ndays = 50
start_time = 200
prefix_u = 'u'
prefix_v = 'v'
PATH = '../model_output/'
ext = 'res_test57.-70.-nu4-urlx-kt4.0.c-0020sat200.0.T85/'
res = 85

if not os.path.isdir(PATH+ext+'winds'):
    os.system('mkdir '+PATH+ext+'winds')

ds = ds_from_BOB(PATH+ext, ['u','v'], 85, time_step=0.5)
# Restrict to NH
nlat = int(len(ds.coords['latitude'])/2)
ds = ds.isel(latitude=slice(None,nlat))

# Resolution for cas winds
nlon = 144
nlat = 92
remainder = 4  # remainder after writing u 5 lon points at a time

# Interpolate to new grid
new_lon = np.arange(0,360,360./nlon)
new_lat = np.arange(0, 90, 90./nlat)[::-1]
ds_cas = ds.interp(latitude=new_lat,longitude=new_lon, kwargs={'fill_value':None})

for iday in range(start_time,start_time+ndays):
    with open(PATH+ext+'winds/u'+'%05d' % iday, 'w') as csvfile:
        print('Writing '+PATH+ext+'winds/u'+'%05d' % iday)
        writer = csv.writer(csvfile, delimiter=' ')
        for ilat in range(nlat):
            for iline in range(nlon//5):
                writer.writerow(ds_cas['u'][iday,ilat,iline*5:5+iline*5].data)
            writer.writerow(ds_cas['u'][iday,ilat,-remainder:].data)

for iday in range(start_time,start_time+ndays):
    with open(PATH+ext+'winds/v'+'%05d' % iday, 'w') as csvfile:
        print('Writing '+PATH+ext+'winds/v'+'%05d' % iday)
        writer = csv.writer(csvfile, delimiter=' ')
        for ilat in range(nlat):
            for iline in range(nlon//5):
                writer.writerow(ds_cas['v'][iday,ilat,iline*5:5+iline*5].data)
            writer.writerow(ds_cas['v'][iday,ilat,-remainder:].data)
