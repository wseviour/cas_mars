'''
Script to make a netCDF file from BOB output
'''

from read_BOB import cube_from_BOB
import glob
import numpy as np
import matplotlib.pyplot as plt
import iris

PATH = '../model_output/'
ext = 'res_test57.-70.-nu4-urlx-kt4.0.c-0020sat200.0.T85'
var = ['q','u','v']
var_name = ['potential_vorticity_of_atmosphere_layer','eastward_wind','northward_wind']
units = [None,'m s-1','m s-1']
res = 85
cubes = iris.cube.CubeList()
for i in range(len(var)):
    files = sorted(glob.glob(PATH+ext+'/'+var[i]+'.?????'))
    cube = cube_from_BOB(files, res, var_name[i],time_step=2.,units=units[i])
    cubes.append(cube)

iris.save(cubes, PATH+ext+'.nc')
