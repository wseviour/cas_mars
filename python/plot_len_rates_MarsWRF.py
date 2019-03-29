from setup_cas import Cas, read_MarsWRF
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os

ndays = 10 # number days for CAS calculation
#start_time = 250 # starting time step of CAS
time_step= 0.25 # time step in model (should be 0.25 for MarsWRF)

start_times = [80,190,200,210,220,230,300,310,320]#range(250,301,10)

levels = range(4,25)

thetas = []
lin_lats = np.arange(45,86,2)
lin_len_rates = np.zeros((len(start_times),len(levels),len(lin_lats))) 

ds = read_MarsWRF()

for start_time in start_times:
    for level in levels:


        ds_lev = ds.isel(theta=level)

        working_dir = '../MarsWRF_data/cas/Z%.2d/' % level
        if os.path.isdir(working_dir):
            pass
        else:
            os.system('mkdir '+working_dir)



        CA = Cas(ds_lev, working_dir, start_time, ndays, time_step)
        CA.interpolate_winds()
        CA.make_contours2(lats=np.arange(50,86,2),plot=False)
        CA.run_cas()
        len_rates = CA.len_rate_eqlat()
        #CA.plot_cas(0)

        lin_len_rates[start_times.index(start_time),levels.index(level),:] = np.interp(lin_lats,len_rates[:,0],len_rates[:,1])
        if start_times.index(start_time) == 0:
            thetas.append(float(ds_lev.theta.data))


    q = ds.q.mean(dim='longitude').isel(time=start_time).isel(theta=slice(levels[0],levels[-1]+1)).sel(latitude=slice(90,40))

    scaling = (np.array(thetas)/200.)**(-1*(1+4.0))
    scaling = np.tile(scaling[:,np.newaxis],(1,len(q.latitude)))

    qs = q * scaling

    plt.contourf(lin_lats, thetas, lin_len_rates[start_times.index(start_time),:],21)
    plt.contour(qs.latitude.data,qs.theta.data,qs.data,colors='k')
    plt.savefig('../plots/len_rate_vs_theta_start_time_%s.pdf' % start_time)
    plt.close()
    #plt.show()
