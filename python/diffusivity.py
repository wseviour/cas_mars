"""
Script to calculate diffusivity from swbob tracer output
"""

import numpy as np
import matplotlib
matplotlib.use('qt5agg')
import matplotlib.pyplot as plt
import xarray as xr
from windspharm.xarray import VectorWind
from scipy import signal

def calc_Koc(ds, detrend=False, plot=False):

    if detrend:
        ds['s'].data = signal.detrend(ds['s'].data, axis=0) + ds['s'].mean(dim='time').data

    s_mean = ds['s'].mean(dim='time')
    s_anom = ds['s'] - s_mean

    w = VectorWind(ds['u'],ds['v'])

    grad_s_mean = w.gradient(s_mean)
    grad_s_mean = VectorWind(grad_s_mean[0],grad_s_mean[1])
    abs_grad_s_mean = grad_s_mean.magnitude()
    sq_abs_grad_s_mean = abs_grad_s_mean * abs_grad_s_mean

    grad_s_anom = w.gradient(s_anom)
    grad_s_anom = VectorWind(grad_s_anom[0],grad_s_anom[1])
    abs_grad_s_anom = grad_s_anom.magnitude()
    sq_abs_grad_s_anom = abs_grad_s_anom * abs_grad_s_anom
    sq_abs_grad_s_anom = sq_abs_grad_s_anom.mean(dim='time')

    sq_abs_grad_s_mean = sq_abs_grad_s_mean.mean(dim='longitude')
    sq_abs_grad_s_anom = sq_abs_grad_s_anom.mean(dim='longitude')


    Koc = sq_abs_grad_s_anom / sq_abs_grad_s_mean

    if plot:

        plt.ion()
        fig1, ax1 = plt.subplots()
        ax1.plot(lats,sq_abs_grad_s_mean.data, label='mean')
        ax1.plot(lats,sq_abs_grad_s_anom.data, label='anom')
        ax11 = ax1.twinx()
        ax11.plot(lats,ds['q'].mean(dim=['time','longitude']).data)


        fig2, ax2 = plt.subplots()
        ax2.plot(lats,Koc, label='Koc')
        ax21 = ax2.twinx()
        ax21.plot(lats,ds['q'].mean(dim=['time','longitude']).data)

        plt.show()

    return Koc



RUNS = '../output/ann57.-70.-nu4-urlx-kt4.0-sinlat.c-0020.T85/ann57.-70.-nu4-urlx-kt4.0-sinlat.c-0020.T85.nc'

#ds = xr.open_dataset(RUNS).isel(time=slice(1000,3000))
ds = xr.open_dataset(RUNS).isel(time=slice(200,2000))
lats = ds.coords['latitude'].data

Koc = calc_Koc(ds,detrend=True,plot=True)
