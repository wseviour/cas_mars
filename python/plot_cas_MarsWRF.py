from setup_cas import Cas, read_MarsWRF
import cartopy.crs as ccrs
import cartopy.util
import matplotlib.path as mpath
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os


#break_times = [247,227]

# 190 - annulus with barrier
# 73 - patch 

ndays = 10 # number days for CAS calculation
start_time = 190 # starting time step of CAS
time_step= 0.25 # time step in model (should be 0.25 for MarsWRF)
level = 18

ds = read_MarsWRF()

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
cons_dic, lens_dic = CA.read_cons_lens()
len_rates = CA.len_rate_eqlat()


keys = sorted(cons_dic.keys())


lats = ds.coords['latitude'].data
lons = ds.coords['longitude'].data

plot_days = [0,3,6] # days to plot

fig = plt.figure(figsize=(16,4))
for iday in plot_days:
    ax = fig.add_subplot(1,len(plot_days),1+plot_days.index(iday), projection=ccrs.NorthPolarStereo())
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)
    ax.set_extent([-180, 180,30, 90], ccrs.PlateCarree())
    ax.gridlines()

    scaling = (ds.theta[level]/200.)**(-1*(1+4.0))
    cyclic_data, cyclic_lons = cartopy.util.add_cyclic_point(ds_lev.q[start_time+int(iday/time_step),:].data*scaling.data, coord = lons) ##
    pv = ax.contourf(cyclic_lons, lats,cyclic_data,cmap='Blues', transform=ccrs.PlateCarree(),extend='both')
    for key in keys:
        try:
            con = cons_dic[key][int(iday/time_step)]
            con = np.append(con[-1,:][np.newaxis,:],con,axis=0)
            if con.shape[0] < 5000:
                ax.plot(con[:,0], con[:,1], transform=ccrs.Geodetic(),color='k', label=key,lw=0.5) #color=colors[keys.index(key)]
        except IndexError:
            pass 
    ax.set_title('$\Theta = %d K$, $L_s = %.1f$' % (ds.theta[level], ds.time[start_time+int(iday/time_step)]))

    # ax = fig.add_subplot(1,len(plot_days)+1,len(plot_days)+1)
    # len_rates = len_rates[len_rates[:,0]<89]
    # ax.plot(len_rates[:,0], len_rates[:,1],marker='.',color='blue')
    # ax.set_ylabel('Stretching rate (sol$^{-1}$)')
    # ax.set_xlabel('Equivalent latitude')
    # ax.tick_params('y', colors='blue')
    # ax.set_xlim(40,90)
    # ax.set_ylim(0,0.5)
    # ax2 = ax.twinx()
    # ax2.plot(lats,ds.q[start_time,level,:].mean(dim='longitude').data,color='red')
    # ax2.set_ylabel('PV')
    # ax2.tick_params('y', colors='red')


plt.tight_layout()
plt.savefig('../plots/MarsWRF_CAS_theta_%d_Ls_%.1f.pdf' % (ds.theta[level], ds.time[start_time]))
plt.show()
