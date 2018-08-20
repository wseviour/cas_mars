import numpy as np
import matplotlib
matplotlib.use('qt5agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.util
import matplotlib.path as mpath
from setup_cas import Cas
from read_BOB import ds_from_BOB
import save_restore

if __name__ == '__main__':

    calc_from_scratch = False
    plot_len_rates = False
    plot_contours = True

    rlax_vals = ['5.0','4.0','1.0','0.0']
    #rlax_vals=['2.0']

    ndays = 10 # number days for CAS calculation
    start_time = 800 # starting time step of CAS
    time_step= 0.25 # time step in SWM
    res = 85 # resolution of SWM

    if calc_from_scratch:

        len_rates = {}
        cons = {}
        for rlax in rlax_vals:

            PATH = '../model_output/'
            ext = 'ann57.-70.-nu4-urlx-kt%s-sinlat.c-0020.T85/' % rlax
            Cas_rlax = Cas(PATH+ext, start_time, ndays, time_step)
            #Cas_rlax.interpolate_winds()
            #Cas_rlax.make_contours(lats=np.arange(55,89,1.5))
            #Cas_rlax.run_cas()
            #len_rates[rlax] = Cas_rlax.len_rate_eqlat()
            icons, ilens = Cas_rlax.read_cons_lens()
            cons[rlax] = icons
        save_restore.save('len_rates.pypic', len_rates=len_rates)
        save_restore.save('cons.pypic', cons=cons)

    elif plot_len_rates:
        len_rates = save_restore.restore('len_rates.pypic')['len_rates']

        fig = plt.figure(figsize=(7,9))
        for rlax in rlax_vals:
            PATH = '../model_output/'
            ext = 'ann57.-70.-nu4-urlx-kt%s-sinlat.c-0020.T85/' % rlax
            ds = ds_from_BOB(PATH+ext, 'q', 85, time_step=time_step)
            q_zm = ds['q'].isel(time=slice(start_time,start_time+int(ndays/time_step)))
            q_zm = q_zm.mean(dim=['longitude','time'])
            q_zm_lats = q_zm.coords['latitude'].data
            try:
                tstd = 1/float(rlax)
            except ZeroDivisionError:
                tstd = 0
            ax1 = fig.add_subplot(2,2,rlax_vals.index(rlax)+1)
            ax1.set_title("$t_r = %.2f$ sols" % tstd if tstd != 0 else r"$t_r \rightarrow \infty$")
            ax1.plot(len_rates[rlax][:,0],len_rates[rlax][:,1],color='b',marker='.')
            ax1.set_ylabel('lengthening rate (1/sol)')
            ax1.tick_params('y', colors='b')
            ax1.set_xlabel('Equivalent latitude')
            ax1.set_xlim(55,90)
            ax1.set_ylim(-0.01,0.45)
            ax2 = ax1.twinx()
            ax2.plot(q_zm_lats,q_zm.data,color='r')
            ax2.set_xlim(55,90)
            ax2.set_ylabel('PV')
            ax2.tick_params('y', colors='r')
            ax2.set_ylim(0,2)
            ax2.set_xticks([55,60,65,70,75,80,85,90])
        plt.tight_layout()
        plt.savefig("../plots/lenrate_rlax.pdf",bbox_inches='tight')
        plt.show()

    elif plot_contours:

        fig = plt.figure(figsize=(7,9))
        for rlax in rlax_vals:
            cons = save_restore.restore('cons.pypic')['cons'][rlax]
            try:
                tstd = 1/float(rlax)
            except ZeroDivisionError:
                tstd = 0
            colors = plt.cm.rainbow(np.linspace(0,1,len(cons.keys())))

            ax = fig.add_subplot(2,2,rlax_vals.index(rlax)+1, projection=ccrs.NorthPolarStereo())
            theta = np.linspace(0, 2*np.pi, 200)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)
            ax.set_boundary(circle, transform=ax.transAxes)
            ax.set_extent([-180, 180,20, 90], ccrs.PlateCarree())
            ax.gridlines()
            colorcount = 0
            for key in cons.keys():
                con = cons[key][20]
                ax.plot(con[:,0], con[:,1], transform=ccrs.Geodetic(),color=colors[colorcount], label=key, lw=0.5)
                colorcount += 1
            ax.set_title("5 sol CA, $t_r = %.2f$ sols" % tstd if tstd != 0 else r"5 sol CA, $t_r \rightarrow \infty$",fontsize=10)
        plt.tight_layout()
        plt.savefig("../plots/contours_rlax.pdf",bbox_inches='tight')
        plt.show()
