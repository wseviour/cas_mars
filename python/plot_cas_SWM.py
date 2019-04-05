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
fig4plot = True

if fig1plot:

    ndays = 10 # number days for CAS calculation
    #start_time = 250 # starting time step of CAS
    time_step= 0.25 # time step in model (should be 0.25 for MarsWRF)

    lin_lats = np.arange(45,86,2)

    kt = 2.0

    start_time = 500

    fig = plt.figure(figsize=(12,4))


    FILE = '/home/bridge/kz18101/Research/cas_mars/model_output/netcdf/ann57.-70.-nu4-urlx-kt%s-sinlat-trfac1.c-0020.T85.nc' % kt
    run_name = FILE.split('/')[-1]

    ds = xr.open_dataset(FILE)

    working_dir = '../model_output/cas/'+run_name+'/'
    if os.path.isdir(working_dir):
        pass
    else:
        os.system('mkdir '+working_dir)

    CA = Cas(ds, working_dir, start_time, ndays, time_step)
    CA.interpolate_winds()
    CA.make_contours2(lats=np.arange(50,86,2),plot=False)
    CA.run_cas()
    cons_dic, lens_dic = CA.read_cons_lens()
    len_rates = CA.len_rate_eqlat()

    keys = sorted(cons_dic.keys())
    
    plot_days = [0,4] # days to plot

    alphabet = ['a','b','c','d','e','f']


    lats = ds.coords['latitude'].data
    lons = ds.coords['longitude'].data


    for iday in plot_days:
        ax = fig.add_subplot(1,len(plot_days)+1,plot_days.index(iday)+1, projection=ccrs.NorthPolarStereo())
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax.set_boundary(circle, transform=ax.transAxes)
        ax.set_extent([-180, 180,40, 90], ccrs.PlateCarree())
        ax.gridlines(linestyle='--',alpha=0.5)

        cyclic_data, cyclic_lons = cartopy.util.add_cyclic_point(ds.q[start_time+int(iday/time_step),:].data, coord = lons) ##
        clevs = np.arange(0,1.5,0.1)
        pv = ax.contourf(cyclic_lons, lats,cyclic_data,clevs,cmap='Greys', transform=ccrs.PlateCarree(),extend='both')
        for key in keys[::2]:
            try:
                con = cons_dic[key][int(iday/time_step)]
                con = np.append(con[-1,:][np.newaxis,:],con,axis=0)
                ax.plot(con[:,0], con[:,1], transform=ccrs.Geodetic(),color='r', label=key,lw=0.5) #color=colors[keys.index(key)]
            except IndexError:
                pass 
        ax.set_title('(%s) $t = %d$ sols' % (alphabet[plot_days.index(iday)],iday))

    len_rates = len_rates[len_rates[:,0] < 90,:]

    ax1 = fig.add_subplot(1,len(plot_days)+1,len(plot_days)+1)
    ax1.plot(len_rates[:,0], len_rates[:,1],marker='.',color='k')
    ax1.set_ylabel('Stretching rate (sol$^{-1}$)')
    ax1.set_xlabel('Equivalent latitude')
    ax1.set_title('(%s)' % alphabet[len(plot_days)])
    #ax1.tick_params('y', colors='blue')
    ax1.set_xlim(50,90)
    ax1.set_ylim(0,0.5)
    ax2 = ax1.twinx()
    ax2.plot(lats,ds.q[start_time,:].mean(dim='longitude').data,color='k',linestyle='--')
    ax2.set_ylabel('PV [$2\Omega/H$]')
    ax2.set_ylim(0.2,1.6)

    cax = fig.add_axes([0.33333-0.12,0.12,0.2,0.05])
    cbar = fig.colorbar(pv, cax=cax, orientation='horizontal', label='PV [$2\Omega/H$]')


    plt.tight_layout()
    plt.savefig('../plots/for_paper/%s_start_%s.pdf' % (run_name, start_time))
    plt.show()

if fig2plot:

    ndays = 10 # number days for CAS calculation
    #start_time = 250 # starting time step of CAS
    time_step= 0.25 # time step in model (should be 0.25 for MarsWRF)

    lin_lats = np.arange(45,86,2)

    start_time = 500

    fig = plt.figure(figsize=(12,4))


    PATH = '/home/bridge/kz18101/Research/cas_mars/model_output/netcdf/'

    FILES = ['ann57.-70.-nu4-urlx-kt4.0-ring.c-0020.T85.nc',
             'ann57.-70.-nu4-urlx-kt1.5-sinlat.c-0020.T85.nc',
             'ann57.-70.-nu4-urlx-kt0.0-sinlat-trfac1.c-0020.T85.nc']

    titles = ['$t_r = 0.25$ sols','$t_r = 0.7$ sols', r"$t_r \rightarrow \infty$"]

    fig = plt.figure(figsize=(8,12))
    
    for ifile in FILES:
        run_name = ifile

        ds = xr.open_dataset(PATH+ifile)

        working_dir = '../model_output/cas/'+run_name+'/'
        if os.path.isdir(working_dir):
            pass
        else:
            os.system('mkdir '+working_dir)

        CA = Cas(ds, working_dir, start_time, ndays, time_step)
        CA.interpolate_winds()
        CA.make_contours2(lats=np.arange(50,86,2),plot=False)
        CA.run_cas()
        cons_dic, lens_dic = CA.read_cons_lens()
        len_rates = CA.len_rate_eqlat()

        keys = sorted(cons_dic.keys())

        plot_days = [4] # days to plot

        alphabet = ['a','b','c','d','e','f']


        lats = ds.coords['latitude'].data
        lons = ds.coords['longitude'].data


        for iday in plot_days:
            ax = fig.add_subplot(len(FILES),len(plot_days)+1,2*FILES.index(ifile)+1, projection=ccrs.NorthPolarStereo())
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)
            ax.set_boundary(circle, transform=ax.transAxes)
            ax.set_extent([-180, 180,40, 90], ccrs.PlateCarree())
            ax.gridlines(linestyle='--',alpha=0.5)

            cyclic_data, cyclic_lons = cartopy.util.add_cyclic_point(ds.q[start_time+int(iday/time_step),:].data, coord = lons) ##
            clevs = np.arange(0,1.5,0.1)
            pv = ax.contourf(cyclic_lons, lats,cyclic_data,clevs,cmap='Greys', transform=ccrs.PlateCarree(),extend='both')
            for key in keys[::2]:
                try:
                    con = cons_dic[key][int(iday/time_step)]
                    con = np.append(con[-1,:][np.newaxis,:],con,axis=0)
                    ax.plot(con[:,0], con[:,1], transform=ccrs.Geodetic(),color='r', label=key,lw=0.5) #color=colors[keys.index(key)]
                except IndexError:
                    pass 
            ax.set_title('(%s) %s' % (alphabet[2*FILES.index(ifile)],titles[FILES.index(ifile)]))

        len_rates = len_rates[len_rates[:,0] < 90,:]

        ax1 = fig.add_subplot(len(FILES),len(plot_days)+1,2*FILES.index(ifile)+2)
        ax1.plot(len_rates[:,0], len_rates[:,1],marker='.',color='k')
        ax1.set_ylabel('Stretching rate (sol$^{-1}$)')
        ax1.set_xlabel('Equivalent latitude')
        ax1.set_title('(%s)' % alphabet[1+2*FILES.index(ifile)])
        #ax1.tick_params('y', colors='blue')
        ax1.set_xlim(40,90)
        ax1.set_ylim(0,0.5)
        ax2 = ax1.twinx()
        ax2.plot(lats,ds.q[start_time,:].mean(dim='longitude').data,color='k',linestyle='--')
        ax2.set_ylabel('PV [$2\Omega/H$]')
        ax2.set_ylim(0.2,1.8)

        #cax = fig.add_axes([0.33333-0.12,0.12,0.2,0.05])
        #cbar = fig.colorbar(pv, cax=cax, orientation='horizontal', label='PV [$2\Omega/H$]')


    plt.tight_layout()
    plt.savefig('../plots/for_paper/SWM_changing_relax.pdf')
    plt.show()


if fig3plot:

    ndays = 10 # number days for CAS calculation
    #start_time = 250 # starting time step of CAS
    time_step= 0.25 # time step in model (should be 0.25 for MarsWRF)

    lin_lats = np.arange(45,86,2)

    start_times = [500,550,550]

    fig = plt.figure(figsize=(12,4))


    PATH = '/home/bridge/kz18101/Research/cas_mars/model_output/netcdf/'

    FILES = ['ann57.-70.-nu4-urlx-kt2.0-sinlat-trfac1-dq20.2.c-0020sat200.0.T85.nc',
             #'ann62.-70.-nu4-urlx-kt2.0-sinlat-trfac1-dq20.5.c-0020sat200.0.T85.nc'
             'ann57.-70.-nu4-urlx-kt2.0-sinlat-trfac1-dq21.0.c-0020sat200.0.T85.nc',
             'ann52.-65.-nu4-urlx-kt2.0-sinlat-trfac1.c-0020sat200.0.T85.nc',
            ]
             

    titles = ['$q_p = 0.2 \cdot 2\Omega/H$','$q_p = 2\Omega/H$', r"$\theta_1 = 55^\circ, \theta_2 = 65^\circ$"]

    fig = plt.figure(figsize=(8,12))
    
    for ifile in FILES:

        start_time = start_times[FILES.index(ifile)]
        
        run_name = ifile

        ds = xr.open_dataset(PATH+ifile)

        working_dir = '../model_output/cas/'+run_name+'/'
        if os.path.isdir(working_dir):
            pass
        else:
            os.system('mkdir '+working_dir)

        CA = Cas(ds, working_dir, start_time, ndays, time_step)
        CA.interpolate_winds()
        CA.make_contours2(lats=np.arange(50,86,2),plot=False)
        CA.run_cas()
        cons_dic, lens_dic = CA.read_cons_lens()
        len_rates = CA.len_rate_eqlat()

        keys = sorted(cons_dic.keys())

        plot_days = [4] # days to plot

        alphabet = ['a','b','c','d','e','f']


        lats = ds.coords['latitude'].data
        lons = ds.coords['longitude'].data


        for iday in plot_days:
            ax = fig.add_subplot(len(FILES),len(plot_days)+1,2*FILES.index(ifile)+1, projection=ccrs.NorthPolarStereo())
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)
            ax.set_boundary(circle, transform=ax.transAxes)
            ax.set_extent([-180, 180,40, 90], ccrs.PlateCarree())
            ax.gridlines(linestyle='--',alpha=0.5)

            cyclic_data, cyclic_lons = cartopy.util.add_cyclic_point(ds.q[start_time+int(iday/time_step),:].data, coord = lons) ##
            clevs = np.arange(0.2,2.0,0.1)
            pv = ax.contourf(cyclic_lons, lats,cyclic_data,clevs,cmap='Greys', transform=ccrs.PlateCarree(),extend='both')
            for key in keys[::2]:
                try:
                    con = cons_dic[key][int(iday/time_step)]
                    con = np.append(con[-1,:][np.newaxis,:],con,axis=0)
                    ax.plot(con[:,0], con[:,1], transform=ccrs.Geodetic(),color='r', label=key,lw=0.5) #color=colors[keys.index(key)]
                except IndexError:
                    pass 
            ax.set_title('(%s) %s' % (alphabet[2*FILES.index(ifile)],titles[FILES.index(ifile)]))

        len_rates = len_rates[len_rates[:,0] < 88,:]

        ax1 = fig.add_subplot(len(FILES),len(plot_days)+1,2*FILES.index(ifile)+2)
        ax1.plot(len_rates[:,0], len_rates[:,1],marker='.',color='k')
        ax1.set_ylabel('Stretching rate (sol$^{-1}$)')
        ax1.set_xlabel('Equivalent latitude')
        ax1.set_title('(%s)' % alphabet[1+2*FILES.index(ifile)])
        #ax1.tick_params('y', colors='blue')
        ax1.set_xlim(40,90)
        ax1.set_ylim(0,0.8)
        ax2 = ax1.twinx()
        ax2.plot(lats,ds.q[start_time,:].mean(dim='longitude').data,color='k',linestyle='--')
        ax2.set_ylabel('PV [$2\Omega/H$]')
        ax2.set_ylim(0.2,1.8)

        #cax = fig.add_axes([0.33333-0.12,0.12,0.2,0.05])
        #cbar = fig.colorbar(pv, cax=cax, orientation='horizontal', label='PV [$2\Omega/H$]')


    plt.tight_layout()
    plt.savefig('../plots/for_paper/SWM_changing_annulus.pdf')
    plt.show()


if fig4plot:

             
    PATH = '/home/bridge/kz18101/Research/cas_mars/model_output/netcdf/'

    FILES = ['ann57.-70.-nu4-urlx-kt2.0-sinlat-trfac1.c-0020.T85.nc',
             'ann57.-70.-nu4-urlx-kt2.0-sinlat-trfac1-dq20.2.c-0020sat200.0.T85.nc',
             'ann57.-70.-nu4-urlx-kt2.0-sinlat-trfac1-dq21.0.c-0020sat200.0.T85.nc',
             'ann52.-65.-nu4-urlx-kt2.0-sinlat-trfac1.c-0020sat200.0.T85.nc',
            ]
             

    titles = ['$q_p = 0.7 \cdot 2\Omega/H$','$q_p = 0.2 \cdot 2\Omega/H$','$q_p = 2\Omega/H$', r"$\theta_1 = 55^\circ, \theta_2 = 65^\circ$"]
    
    ls = ['-','--',':','-.']
    
    fig = plt.figure(figsize=(5,5))
    ax1 = fig.add_subplot(111)
    
    for ifile in FILES:
        run_name = ifile

        ds = xr.open_dataset(PATH+ifile)
        lats = ds.latitude.data

        
        ax1.plot(lats, ds.q[0,:,0].data, label = titles[FILES.index(ifile)],color='k',ls=ls[FILES.index(ifile)])

    ax1.set_xlim(40,90)
    ax1.set_ylim(0,1.8)
    ax1.legend()
    ax1.set_xlabel('latitude')
    ax1.set_ylabel('PV [$2\Omega/H$]')
    plt.tight_layout()
    plt.savefig('../plots/for_paper/pv_profiles.pdf')
    plt.show()
