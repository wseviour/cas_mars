from pyproj import Proj
from shapely.geometry import shape
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
import pylab

def pstere_proj(ds,timestep=600):

    d = ds.s.isel(time=timestep).sel(latitude = slice(90,10))
    pa = Proj("+proj=stere +lat_0=90",preserve_units=True)
    lonv, latv = np.meshgrid(d.longitude.data, d.latitude.data)
    x, y = pa(lonv,latv)
    reg_x = np.linspace(np.min(x),np.max(x),100)
    reg_y = np.linspace(np.min(y),np.max(y),100)    
    xi, yi = np.meshgrid(reg_x, reg_y)
    d2 = mlab.griddata(x.flatten(),y.flatten(),d.data.flatten(),xi,yi,interp='linear')

    return d2, xi, yi

def calc_keff(c, dx, dy, rac, timestep, ds, lats=np.linspace(30,88,2)):
   
    R = 6371.e3

    A = 2*np.pi*(R**2)*(1-np.sin(np.deg2rad(lats)))
    Ncon = len(lats)
    cons = ds.s.isel(time = timestep).mean(dim='longitude').interp(latitude=lats).data

    grad_c = np.gradient(c)
    grad_c_y = grad_c[0] / dy
    grad_c_x = grad_c[1] / dx
    # we need to integrate |grad c|^2, this makes it easy
    grad_c2xRAC = (grad_c_x**2 + grad_c_y**2) * rac.filled(0.)


    A_C_lin = np.zeros(Ncon) # the area inside each contour
    grad_c2_C_lin = np.zeros(Ncon) # integral of c^2 inside contour
    for n in range(Ncon):
        mask = np.ma.masked_where(c <= cons[n], c).mask
        A_C_lin[n] = np.ma.masked_array(rac, mask=mask).sum()
        #plt.imshow(np.ma.masked_array(rac, mask=mask))
        #plt.show()
        grad_c2_C_lin[n] = np.ma.masked_array(grad_c2xRAC, mask=mask).sum()

    #C = np.interp(A, A_C_lin, cons)
    #X = np.interp(A, A_C_lin, grad_c2_C_lin)

    #dCdA = np.gradient(C) / np.gradient(A)
    #dXdA = np.gradient(X) / np.gradient(A)
    
    dCdA = np.gradient(cons, A_C_lin)
    dXdA = np.gradient(grad_c2_C_lin, A_C_lin)

    Le2 = dXdA / (dCdA**2)

    Eqlat = np.rad2deg(np.arcsin(1 - A_C_lin/(2*np.pi*(R**2))))

    keff = Le2/((2*np.pi*R*np.cos(np.deg2rad(Eqlat)))**2)

    return Eqlat, keff


if __name__ == '__main__':

    kt = ['2.0']#['0.0','0.4','1.0','2.0']#['0.0','1.0','2.0']
    dq = ['','-dq20.2']
    lat = [65]
    res = 170
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    #ax2 = fig.add_subplot(212)

    FILES = ['ann57.-70.-nu4-urlx-kt2.0-sinlat-trfac1-dq20.2.c-0020sat200.0.T85.nc',
             'ann57.-70.-nu4-urlx-kt2.0-sinlat-trfac1.c-0020.T85.nc',
             #'ann57.-70.-nu4-urlx-kt2.0-sinlat-trfac1-dq21.0.c-0020sat200.0.T85.nc'
            ]

    titles = ['$q_p = 0.2 \cdot 2\Omega/H$','$q_p = 0.7 \cdot 2\Omega/H$','$q_p = 2\Omega/H$']

    ls = ['--','-',':']
    
    for ifile in FILES:
    
        FILE = '/home/bridge/kz18101/Research/cas_mars/model_output/netcdf/'+ifile

        timesteps = range(500,550)
        #timesteps=[450]
        lats = np.arange(30,88,2)
        #lats = [55]
        keffi = np.zeros((len(lats), len(timesteps)))

        R = 6371e3

        ds = xr.open_dataset(FILE)

        for timestep in timesteps:

            print(timestep)
            c, xi, yi = pstere_proj(ds,timestep=timestep)

            dx = np.ones_like(c) * np.diff(xi)[0,0]
            dy = np.ones_like(c) * np.diff(xi)[0,0]
            rac = np.ones_like(c) * np.diff(xi)[0,0] * np.diff(xi)[0,0]

            Eqlat, keff = calc_keff(c, dx, dy, rac, timestep, ds, lats=lats)

            keffi[:,timesteps.index(timestep)] = np.interp(lats,Eqlat,keff)


        keffi = np.ma.masked_array(keffi, mask = np.isnan(keffi))
        ax1.plot(lats, np.ma.mean(keffi,axis=1),color='k',ls=ls[FILES.index(ifile)],label=titles[FILES.index(ifile)])
        ax2.plot(ds.coords['latitude'].data, ds.q.mean(dim='longitude')[timesteps[0]:timesteps[-1],:].mean(dim='time').data,color='r',ls=ls[FILES.index(ifile)])
        
    ax1.set_ylabel('$\kappa_{\mathrm{eff}}/\kappa$',fontsize=12)
    ax1.set_xlabel('$\phi_e$',fontsize=12)
    ax1.set_xlim(45,90)
    #ax2.set_xlim(40,90)
    ax1.set_ylim(1,30)
    ax2.set_ylim(0.0,1.9)
    ax2.tick_params('y',colors='red')
    ax2.set_ylabel('PV [$2\Omega/H$]')
    ax1.legend()
    plt.savefig('../plots/for_paper/keff.pdf')
    plt.show()

