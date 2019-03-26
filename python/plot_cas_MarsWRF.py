from setup_cas import Cas
import cartopy.crs as ccrs
import cartopy.util
import matplotlib.path as mpath
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

def read_MarsWRF(run_dir, level):

    # Resolution for cas winds
    nlon = 144
    nlat = 92
    
    ds_raw = xr.open_mfdataset(run_dir+'*.nc', concat_dim='Time')
    ds_raw = ds_raw.where(ds_raw.L_S < 300, drop=True) # Drop missing values in MarsWRF files

    lat = ds_raw.LAT[0,:,0]
    lon = ds_raw.LON[0,0,:] % 360 # make sure in 0-360 format, not -180-180
    lon = np.append(lon[90:],lon[:90])
    theta = ds_raw.THETA[0,:]
    l_s = ds_raw.L_S

    U = ds_raw.U.data
    V = ds_raw.V.data
    EPV = ds_raw.EPV.data
    U = np.append(U[:,:,:,90:],U[:,:,:,:90], axis=3)
    V = np.append(V[:,:,:,90:],V[:,:,:,:90], axis=3)
    EPV = np.append(EPV[:,:,:,90:],EPV[:,:,:,:90], axis=3)            

    ds = xr.Dataset({'u':(['l_s','theta','latitude','longitude'], U),
                     'v':(['l_s','theta','latitude','longitude'], V),
                     'q':(['l_s','theta','latitude','longitude'], EPV)},
                     coords = {'l_s':('l_s',l_s),
                               'theta':('theta',theta),
                               'latitude':('latitude',lat),
                               'longitude':('longitude',lon)})

    # ds = ds.isel(theta=level) # Select given level
    # ds = ds.sel(latitude=slice(0,90))
    # new_lon = np.arange(0,360,360./nlon)
    # new_lat = np.arange(0, 90, 90./nlat)[::-1]
    # ds_cas = ds.interp(latitude=new_lat,longitude=new_lon, kwargs={'fill_value':None}) 
    

    return ds

start_times = [85]#85[130,260]
levels = [4,6]

for start_time in start_times:
    for level in levels:
#start_time = 20 #260 # starting time step


        time_step = 0.25 # time step size in data
        ndays = 10 # number of days for CAS run
        #level = 6 #6 # vertical level

        try:
            del CA
        except NameError:
            pass

        CA = Cas('../MarsWRF_data/',start_time,ndays,time_step,model_type='MarsWRF',level=level)
        CA.interpolate_winds()
        CA.make_contours(con_var='q',lats=np.arange(50,86,2),plot=False)
        CA.run_cas()
        cons_dic, lens_dic = CA.read_cons_lens()
        len_rates = CA.len_rate_eqlat()

        keys = sorted(cons_dic.keys())
        indices = [x for x, s in enumerate(keys) if 'inner' in s]
        if len(indices) > 0:
            keys[indices[0]:] = sorted(keys[indices[0]:],reverse=True)
        colors = plt.cm.rainbow(np.linspace(0,1,len(keys)))

        ds = read_MarsWRF('../MarsWRF_data/', level)
        lats = ds.coords['latitude'].data
        lons = ds.coords['longitude'].data

        plot_days = [0,2,4] # days to plot

        fig = plt.figure(figsize=(16,4))
        for iday in plot_days:
            ax = fig.add_subplot(1,len(plot_days)+1,1+plot_days.index(iday), projection=ccrs.NorthPolarStereo())
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)
            ax.set_boundary(circle, transform=ax.transAxes)
            ax.set_extent([-180, 180,30, 90], ccrs.PlateCarree())
            ax.gridlines()
            if level == 4:
                clevs = np.linspace(0,0.001,21)
            else:
                clevs = np.linspace(0,0.003,21)
            cyclic_data, cyclic_lons = cartopy.util.add_cyclic_point(ds.q[start_time+int(iday/time_step),level,:].data, coord = lons) ##
            pv = ax.contourf(cyclic_lons, lats,cyclic_data, clevs,cmap='Blues', transform=ccrs.PlateCarree(),extend='both')
            for key in keys:
                try:
                    con = cons_dic[key][int(iday/time_step)]
                    con = np.append(con[-1,:][np.newaxis,:],con,axis=0)
                    if con.shape[0] < 5000:
                        ax.plot(con[:,0], con[:,1], transform=ccrs.Geodetic(),color='k', label=key,lw=0.5) #color=colors[keys.index(key)]
                except IndexError:
                    pass 
            ax.set_title('$\Theta = %d K$, $L_s = %.1f$' % (ds.theta[level], ds.l_s[start_time+int(iday/time_step)]))

        ax = fig.add_subplot(1,len(plot_days)+1,len(plot_days)+1)
        len_rates = len_rates[len_rates[:,0]<89]
        ax.plot(len_rates[:,0], len_rates[:,1],marker='.',color='blue')
        ax.set_ylabel('Stretching rate (sol$^{-1}$)')
        ax.set_xlabel('Equivalent latitude')
        ax.tick_params('y', colors='blue')
        ax.set_xlim(40,90)
        ax.set_ylim(0,0.5)
        ax2 = ax.twinx()
        ax2.plot(lats,ds.q[start_time,level,:].mean(dim='longitude').data,color='red')
        ax2.set_ylabel('PV')
        ax2.tick_params('y', colors='red')


        plt.tight_layout()
        plt.savefig('../plots/MarsWRF_CAS_theta_%d_Ls_%.1f.pdf' % (ds.theta[level], ds.l_s[start_time]))
        plt.show()
