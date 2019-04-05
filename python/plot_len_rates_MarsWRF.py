from setup_cas import Cas, read_MarsWRF
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
import save_restore

ndays = 10 # number days for CAS calculation
#start_time = 250 # starting time step of CAS
time_step= 0.25 # time step in model (should be 0.25 for MarsWRF)

start_times = [240]#range(70,400,8)#[192,200,208]#range(250,301,10)

levels = range(4,25)


lin_lats = np.arange(45,86,2)
lin_len_rates = np.zeros((len(start_times),len(levels),len(lin_lats))) 

#ds = read_MarsWRF()
ds = xr.open_dataset('../MarsWRF_data/isentropic_surface_data_wjms.nc')
                     
# for start_time in start_times:
#     for level in levels:


#         ds_lev = ds.isel(theta=level)

#         working_dir = '../MarsWRF_data/cas/Z%.2d/' % level
#         if os.path.isdir(working_dir):
#             pass
#         else:
#             os.system('mkdir '+working_dir)



#         CA = Cas(ds_lev, working_dir, start_time, ndays, time_step)
#         if start_times.index(start_time) == 0:
#             CA.interpolate_winds(length=start_times[-1]-start_times[0]+40)
#         CA.make_contours2(lats=np.arange(50,86,2),plot=False)
#         CA.run_cas()
#         len_rates = CA.len_rate_eqlat()
#         #CA.plot_cas(0)

#         lin_len_rates[start_times.index(start_time),levels.index(level),:] = np.interp(lin_lats,len_rates[:,0],len_rates[:,1])
#         # if start_times.index(start_time) == 0:
#         #     thetas.append(float(ds_lev.theta.data))



a= save_restore.restore('len_rate_vs_theta_start_time_200.pypic')
lin_len_rates = a['lin_len_rates']

thetas = ds.theta[4:25].data

q = ds.q.mean(dim='longitude').isel(time=start_times[0]).isel(theta=slice(levels[0],levels[-1]+1)).sel(latitude=slice(90,40))

scaling = (np.array(thetas)/200.)**(-1*(1+4.0))
scaling = np.tile(scaling[:,np.newaxis],(1,len(q.latitude)))

qs = q * scaling

l = plt.pcolormesh(lin_lats, thetas, np.mean(lin_len_rates[0:2,:],axis=0))
c = plt.contour(qs.latitude.data,qs.theta.data,qs.data*1e5,15,colors='k')
plt.clabel(c,fmt='%.0f')
plt.ylabel('Isentropic level [K]')
plt.xlabel('Equivalent latitude')
plt.xlim(45,90)
plt.colorbar(l,label= 'stretching rate [sol$^{-1}$]')
plt.savefig('../plots/for_paper/len_rate_vs_theta_start_time_240.pdf')
plt.show()

    
    
#save_restore.save('len_rate_vs_time_theta_%d.pypic' % level,  lin_len_rates = lin_len_rates, lin_lats=lin_lats, start_times=start_times, theta=ds.theta[level].data)

# a = save_restore.restore('len_rate_vs_time_theta_6.pypic')
# lin_len_rates = a['lin_len_rates']
# lin_lats = a['lin_lats']
# start_times = a['start_times']

# q = ds.q.mean(dim='longitude').isel(theta=level).isel(time=slice(70,400))
# scaling = (300/200.)**(-1*(1+4.0))
# q = q*scaling*1e5

# t_to_ls = 0.1346

# fig = plt.figure()
# ax1 = fig.add_subplot(111)
# lens = ax1.pcolormesh(190+(t_to_ls*np.arange(70,400,8)),lin_lats, lin_len_rates[:,0,:].T)
# pv = ax1.contour(190+(t_to_ls*q.time.data), q.latitude.data, q.data.T, np.arange(0,50,5),colors='k')
# ax1.clabel(pv, fmt='%.0f')
# ax1.set_ylim(45,86)
# plt.colorbar(lens, label='stretching rate [sol$^{-1}$]')
# ax1.set_ylabel('equivalent latitude')
# ax1.set_xlabel('time [approx. $L_S$]')
# plt.savefig('../plots/tmp.pdf')
# plt.show()

