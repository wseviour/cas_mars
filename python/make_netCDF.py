'''
Script to make a netCDF file from BOB output
'''

from read_BOB import ds_from_BOB
import glob
import numpy as np
import matplotlib.pyplot as plt
import sys

if sys.argv:
    print(sys.argv[1])
    ext = sys.argv[1]
else:
    ext = 'ann57.-70.-nu4-urlx-kt2.0-sinlat.c-0020.T85'

PATH = '../model_output/raw/'

res = int(ext.split('T')[-1])

ds = ds_from_BOB(PATH+ext, ['q','s','u','v','h'], res, time_step=0.25) 

ds.to_netcdf(PATH.split('raw')[0]+'netcdf/'+ext+'.nc')

