from setup_cas import Cas, read_MarsWRF
import cartopy.crs as ccrs
import cartopy.util
import matplotlib.path as mpath
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os

fig1plot = False
fig2plot = False
fig3plot = False
fig4plot = False
fig5plot = True


#break_times = [247,227]

# 190 - annulus with barrier
# 73 - patch 
if fig1plot:
    ndays = 10 # number days for CAS calculation
    start_time = 240 # starting time step of CAS
    time_step= 0.25 # time step in model (should be 0.25 for MarsWRF)
    level = 6

    #ds = read_MarsWRF()
    ds = xr.open_dataset('../MarsWRF_data/isentropic_surface_data_wjms.nc')
    
    ds_lev = ds.isel(theta=level)


    working_dir = '../MarsWRF_data/cas/Z%.2d/' % level
    if os.path.isdir(working_dir):
        pass
    else:
        os.system('mkdir '+working_dir)


    CA = Cas(ds_lev, working_dir, start_time, ndays, time_step)
    CA.interpolate_winds()
    CA.make_contours2(lats=np.arange(50,86,4),plot=False)
    CA.run_cas()
    cons_dic, lens_dic = CA.read_cons_lens()
    len_rates = CA.len_rate_eqlat()


    keys = sorted(cons_dic.keys())


    lats = ds.coords['latitude'].data
    lons = ds.coords['longitude'].data

    plot_days = [0,4] # days to plot

    alphabet = ['a','b','c']

    eqlats = {'q_0.0008713_tstep_240.in.cas':52.7,
              'q_0.0019731_tstep_240.in.cas':63.1,
              'q_0.0026796_tstep_240_inner.in.cas':81.7}
    
    fig = plt.figure(figsize=(12,4))
    for iday in plot_days:
        ax = fig.add_subplot(1,len(plot_days)+1,plot_days.index(iday)+1, projection=ccrs.NorthPolarStereo())
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax.set_boundary(circle, transform=ax.transAxes)
        ax.set_extent([-180, 180,40, 90], ccrs.PlateCarree())
        ax.gridlines(linestyle='--',alpha=0.5)

        scaling = (ds.theta[level]/200.)**(-1*(1+4.0))
        cyclic_data, cyclic_lons = cartopy.util.add_cyclic_point(ds_lev.q[start_time+int(iday/time_step),:].data*scaling.data, coord = lons) ##
        clevs = np.arange(8,41,2)
        pv = ax.contourf(cyclic_lons, lats,cyclic_data*1e5,clevs,cmap='Greys', transform=ccrs.PlateCarree(),extend='both')
        for key in keys:
            try:
                if key in eqlats.keys():
                    con = cons_dic[key][int(iday/time_step)]
                    con = np.append(con[-1,:][np.newaxis,:],con,axis=0)
                    ax.plot(con[:,0], con[:,1], transform=ccrs.Geodetic(),color='r', label=key,lw=2) #color=colors[keys.index(key)]
                else:
                    con = cons_dic[key][int(iday/time_step)]
                    con = np.append(con[-1,:][np.newaxis,:],con,axis=0)
                    ax.plot(con[:,0], con[:,1], transform=ccrs.Geodetic(),color='r', label=key,lw=0.5) #color=colors[keys.index(key)]
            except IndexError:
                pass 
        ax.set_title('(%s) $t = %d$ sols' % (alphabet[plot_days.index(iday)], iday))

    ax = fig.add_subplot(1,len(plot_days)+1,len(plot_days)+1)
    for key in keys:
        if key in eqlats.keys():
            ax.plot(np.arange(0,10,0.25),np.array(lens_dic[key])/lens_dic[key][0],color='k')
            ax.text(8,np.array(lens_dic[key])[-2]/lens_dic[key][0]+1,'$\phi_e = %.1f^{\circ}$' % eqlats[key] )
    ax.set_ylabel('length($t$) / length(0)')
    ax.set_xticks([0,2,4,6,8,10])
    ax.set_yticks([1,3,5,7,9,11,13])
    ax.set_xlabel('$t$ [sols]')
    ax.hlines(1,0,10,linestyles='--')
    ax.set_xlim(0,11)
    ax.set_ylim(0,15)
    ax.set_title('(c)')


    cax = fig.add_axes([0.33333-0.1,0.12,0.2,0.05])
    cbar = fig.colorbar(pv, cax=cax, orientation='horizontal', label='PV [$10^{-5}$ Km$^2$ kg$^{-1}$ s$^{-1}]$')

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
    plt.savefig('../plots/for_paper/MarsWRF_CAS_theta_%d_time_%d.pdf' % (ds.theta[level], ds.time[start_time]))
    plt.show()


if fig2plot:

    ndays = 10 # number days for CAS calculation
    start_time = 240 # starting time step of CAS
    time_step= 0.25 # time step in model (should be 0.25 for MarsWRF)
    level = 6

    #ds = read_MarsWRF()
    ds = xr.open_dataset('../MarsWRF_data/isentropic_surface_data_wjms.nc')
    
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

    lats = ds.coords['latitude'].data
    lons = ds.coords['longitude'].data

    
    fig2 = plt.figure(figsize=(4,4))

    ax1 = fig2.add_subplot(111)
    ax1.plot(len_rates[:,0], len_rates[:,1],marker='.',color='k')
    ax1.set_ylabel('Stretching rate (sol$^{-1}$)')
    ax1.set_xlabel('Equivalent latitude')
    #ax1.tick_params('y', colors='blue')
    ax1.set_xlim(50,90)
    ax1.set_ylim(0,0.5)
    ax2 = ax1.twinx()
    ax2.plot(lats,ds_lev.q[start_time,:].mean(dim='longitude').data*1e5,color='k',linestyle='--')
    ax2.set_ylabel('PV [$10^{-5}$ Km$^2$ kg$^{-1}$ s$^{-1}$]')
    ax2.set_ylim(50,300)
    #ax1.tick_params('y', colors='red')


    # ax3 = fig2.add_subplot(212)
    # ax3.plot(lats, ds_lev.u[start_time,:].mean(dim='longitude').data, color='k')
    # ax3.set_xlim(40,90)


    plt.tight_layout()
    plt.savefig('../plots/for_paper/MarsWRF_lenrates_theta_%d_Time_%d.pdf' % (ds.theta[level], ds.time[start_time]))
    plt.show()



if fig3plot:
    ndays = 10 # number days for CAS calculation
    start_times = [200,320]#300 # starting time step of CAS
    time_step= 0.25 # time step in model (should be 0.25 for MarsWRF)
    level = 6

    #ds = read_MarsWRF()
    ds = xr.open_dataset('../MarsWRF_data/isentropic_surface_data_wjms.nc')

    ds_lev = ds.isel(theta=level)

    fig = plt.figure(figsize=(12,8))
    
    for start_time in start_times:
    
    
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

        plot_days = [0,4] # days to plot

        alphabet = ['a','b','c','d','e','f']

        
        for iday in plot_days:
            ax = fig.add_subplot(2,len(plot_days)+1,plot_days.index(iday)+1+start_times.index(start_time)*(len(plot_days)+1), projection=ccrs.NorthPolarStereo())
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)
            ax.set_boundary(circle, transform=ax.transAxes)
            ax.set_extent([-180, 180,40, 90], ccrs.PlateCarree())
            ax.gridlines(linestyle='--',alpha=0.5)

            scaling = (ds.theta[level]/200.)**(-1*(1+4.0))
            cyclic_data, cyclic_lons = cartopy.util.add_cyclic_point(ds_lev.q[start_time+int(iday/time_step),:].data*scaling.data, coord = lons) ##
            clevs = np.arange(8,41,2)
            pv = ax.contourf(cyclic_lons, lats,cyclic_data*1e5,clevs,cmap='Greys', transform=ccrs.PlateCarree(),extend='both')
            for key in keys[::2]:
                try:
                    con = cons_dic[key][int(iday/time_step)]
                    con = np.append(con[-1,:][np.newaxis,:],con,axis=0)
                    ax.plot(con[:,0], con[:,1], transform=ccrs.Geodetic(),color='r', label=key,lw=0.5) #color=colors[keys.index(key)]
                except IndexError:
                    pass 
            ax.set_title('(%s) $t = %d$ sols' % (alphabet[plot_days.index(iday)+start_times.index(start_time)*(len(plot_days)+1)],iday))

        len_rates = len_rates[len_rates[:,0] < 90,:]

        ax1 = fig.add_subplot(2,len(plot_days)+1,len(plot_days)+1+start_times.index(start_time)*(len(plot_days)+1))
        ax1.plot(len_rates[:,0], len_rates[:,1],marker='.',color='k')
        ax1.set_ylabel('Stretching rate (sol$^{-1}$)')
        ax1.set_xlabel('Equivalent latitude')
        ax1.set_title('(%s)' % alphabet[len(plot_days)+start_times.index(start_time)*(len(plot_days)+1)])
        #ax1.tick_params('y', colors='blue')
        ax1.set_xlim(50,90)
        ax1.set_ylim(0,0.5)
        ax2 = ax1.twinx()
        ax2.plot(lats,ds_lev.q[start_time,:].mean(dim='longitude').data*1e5,color='k',linestyle='--')
        ax2.set_ylabel('PV [$10^{-5}$ Km$^2$ kg$^{-1}$ s$^{-1}$]')
        #ax2.set_ylim(50,300)

    cax = fig.add_axes([0.33333-0.1,0.04,0.2,0.02])
    cbar = fig.colorbar(pv, cax=cax, orientation='horizontal', label='PV [$10^{-5}$ Km$^2$ kg$^{-1}$ s$^{-1}]$')


    plt.tight_layout()
    plt.savefig('../plots/for_paper/MarsWRF_varyingtimes.pdf',bbox_inches='tight')
    plt.show()



if fig4plot:
    ndays = 10 # number days for CAS calculation
    start_time = 200 # starting time step of CAS
    time_step= 0.25 # time step in model (should be 0.25 for MarsWRF)
    level = 16

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

    plot_days = [0,4] # days to plot

    alphabet = ['a','b','c']

    fig = plt.figure(figsize=(12,4))
    for iday in plot_days:
        ax = fig.add_subplot(1,len(plot_days)+1,plot_days.index(iday)+1, projection=ccrs.NorthPolarStereo())
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax.set_boundary(circle, transform=ax.transAxes)
        ax.set_extent([-180, 180,40, 90], ccrs.PlateCarree())
        ax.gridlines(linestyles='--')

        scaling = (ds.theta[level]/200.)**(-1*(1+4.0))
        cyclic_data, cyclic_lons = cartopy.util.add_cyclic_point(ds_lev.q[start_time+int(iday/time_step),:].data*scaling.data, coord = lons) ##
        clevs = np.arange(8,61,2)
        pv = ax.contourf(cyclic_lons, lats,cyclic_data*1e5,clevs,cmap='Blues', transform=ccrs.PlateCarree(),extend='both')
        for key in keys:
            try:
                con = cons_dic[key][int(iday/time_step)]
                con = np.append(con[-1,:][np.newaxis,:],con,axis=0)
                ax.plot(con[:,0], con[:,1], transform=ccrs.Geodetic(),color='k', label=key,lw=0.5) #color=colors[keys.index(key)]
            except IndexError:
                pass 
        ax.set_title('(%s) $t = %d$ sols, $L_s = %.1f^{\circ}$' % (alphabet[plot_days.index(iday)],iday, ds.time[start_time+int(iday/time_step)]))

        ax1 = fig.add_subplot(1,3,3)
        ax1.plot(len_rates[:,0], len_rates[:,1],marker='.',color='k')
        ax1.set_ylabel('Stretching rate (sol$^{-1}$)')
        ax1.set_xlabel('Equivalent latitude')
        ax1.set_title('(%s)' % alphabet[2])
        #ax1.tick_params('y', colors='blue')
        ax1.set_xlim(50,90)
        ax1.set_ylim(0,1.0)
        ax2 = ax1.twinx()
        ax2.plot(lats,ds_lev.q[start_time,:].mean(dim='longitude').data*scaling.data*1e5,color='k',linestyle='--')
        ax2.set_ylabel('PV [$10^{-5}$ Km$^2$ kg$^{-1}$ s$^{-1}$]')
        #ax2.set_ylim(50,300)


    cax = fig.add_axes([0.33333-0.12,0.12,0.2,0.05])
    cbar = fig.colorbar(pv, cax=cax, orientation='horizontal', label='PV [$10^{-5}$ Km$^2$ kg$^{-1}$ s$^{-1}]$')
    plt.tight_layout()
    plt.savefig('../plots/for_paper/MarsWRF_lenrates_theta_%d_Ls_%.1f.pdf' % (ds.theta[level], ds.time[start_time]))
    plt.show()



if fig5plot:
    ndays = 10 # number days for CAS calculation
    start_time = 240 # starting time step of CAS
    time_step= 0.25 # time step in model (should be 0.25 for MarsWRF)
    level = 6

    #ds = read_MarsWRF()
    ds = xr.open_dataset('../MarsWRF_data/isentropic_surface_data_wjms.nc')
    
    ds = ds.isel(theta=level).isel(time=slice(start_time,start_time+40)).mean(dim='time').mean(dim='longitude')

    shear = ds.u.differentiate('latitude')
    dqdl =  ds.q.differentiate('latitude')
    
    
    kt = 2.0
    start_time = 500
    FILE = '/home/bridge/kz18101/Research/cas_mars/model_output/netcdf/ann57.-70.-nu4-urlx-kt%s-sinlat-trfac1.c-0020.T85.nc' % kt
    run_name = FILE.split('/')[-1]
    ds_swm = xr.open_dataset(FILE)
    ds_swm = ds_swm.isel(time=slice(start_time,start_time+40)).mean(dim='time').mean(dim='longitude')

    shear_swm = ds.u.differentiate('latitude')
    dqdl_swm =  ds.q.differentiate('latitude')
    

    fig = plt.figure(figsize=(8,7))
    ax1 = fig.add_subplot(221)
    ax1.plot(ds.latitude.data, ds.u.data,color='k')
    ax12 = ax1.twinx()
    ax12.plot(ds_swm.latitude, ds_swm.u.data,color='r')
    ax1.set_xlim(40,90)
    ax1.vlines([60, 66], 0, 1, transform=ax1.get_xaxis_transform(), colors='r',linestyles='--',alpha=0.5)
    ax1.vlines([59, 70], 0, 1, transform=ax1.get_xaxis_transform(), colors='k',linestyles='--',alpha=0.5)
    ax12.tick_params('y', colors='red')
    ax1.set_title('(a) $\overline{u}$')
    ax1.set_ylabel('m s$^{-1}$')
    ax12.set_ylabel('m s$^{-1}$')

    
    ax2 = fig.add_subplot(222)
    ax2.plot(ds.latitude.data, ds.q.data*1e5,color='k')
    ax22 = ax2.twinx()
    ax22.plot(ds_swm.latitude, ds_swm.q.data,color='r')
    ax2.set_xlim(40,90)
    ax2.vlines([60, 66], 0, 1, transform=ax2.get_xaxis_transform(), colors='r',linestyles='--',alpha=0.5)
    ax2.vlines([59, 70], 0, 1, transform=ax2.get_xaxis_transform(), colors='k',linestyles='--',alpha=0.5)
    ax22.set_ylim(0,1.5)
    ax22.tick_params('y', colors='red')
    ax2.set_title('(b) $\overline{q}$')
    ax2.set_ylabel('$10^{-5}$ Km$^2$ kg$^{-1}$ s$^{-1}$')
    ax22.set_ylabel('$2\Omega/H$')

    
    ax3 = fig.add_subplot(223)
    ax3.plot(ds.latitude.data, ds.u.differentiate('latitude').data,color='k')
    ax32 = ax3.twinx()
    ax32.plot(ds_swm.latitude, ds_swm.u.differentiate('latitude').data,color='r')
    ax3.set_xlim(40,90)
    ax3.vlines([60, 66], 0, 1, transform=ax3.get_xaxis_transform(), colors='r',linestyles='--',alpha=0.5)
    ax3.vlines([59, 70], 0, 1, transform=ax3.get_xaxis_transform(), colors='k',linestyles='--',alpha=0.5)
    ax3.set_ylim(-5,4)
    ax32.set_ylim(-5,4)
    ax32.tick_params('y', colors='red')
    ax3.set_title('(c) $\mathrm{d}\overline{u}/\mathrm{d}\phi$')
    ax3.set_ylabel('m s$^{-1}$ deg$^{-1}$')
    ax32.set_ylabel('m s$^{-1}$ deg$^{-1}$')
    ax3.set_xlabel('latitude')
    
    ax4 = fig.add_subplot(224)
    ax4.plot(ds.latitude.data, ds.q.differentiate('latitude').data*1e5,color='k')
    ax42 = ax4.twinx()
    ax42.plot(ds_swm.latitude, ds_swm.q.differentiate('latitude').data,color='r')
    ax4.set_xlim(40,90)
    ax4.vlines([60, 66], 0, 1, transform=ax4.get_xaxis_transform(), colors='r',linestyles='--',alpha=0.5)
    ax4.vlines([59, 70], 0, 1, transform=ax4.get_xaxis_transform(), colors='k',linestyles='--',alpha=0.5)
    ax42.tick_params('y', colors='red')
    ax4.set_title('(d) $\mathrm{d}\overline{q}/\mathrm{d}\phi$')
    ax4.set_ylabel('$10^{-5}$ Km$^2$ kg$^{-1}$ s$^{-1}$ deg$^{-1}$')
    ax42.set_ylabel('$2\Omega/H$ deg$^{-1}$')
    ax4.set_xlabel('latitude')
    ax4.set_ylim(-6,15)

    plt.tight_layout()
    plt.savefig('../plots/for_paper/zonal_diags_comparison.pdf')
    plt.show()


    


    
    
