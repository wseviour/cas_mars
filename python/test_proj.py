import numpy as np
from pyproj import Proj
import xarray as xr
#from scipy.interpolate import griddata  
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.util
import matplotlib.path as mpath

d = xr.open_dataset('/home/bridge/kz18101/Research/cas_mars/model_output/netcdf/ann57.-70.-nu4-urlx-kt5.0-hat.c-0020.T85.nc').q
d = d.sel(latitude = slice(90,20))

pa = Proj("+proj=stere +lat_0=90 +lon_0=0")
lonv, latv = np.meshgrid(d.longitude.data, d.latitude.data)
x, y = pa(lonv,latv)

xi, yi = np.meshgrid(np.linspace(np.min(x),np.max(x),1000), np.linspace(np.min(y),np.max(y),1000))



d2 = mlab.griddata(x.flatten(),y.flatten(),d.data[100,:].flatten(),xi,yi,interp='linear')

lats = d.coords['latitude'].data
lons = d.coords['longitude'].data

fig = plt.figure(figsize=(10,5))
ax1 = fig.add_subplot(1,2,1)
ax1.contourf(np.linspace(np.min(x),np.max(x),1000),np.linspace(np.min(y),np.max(y),1000),d2)
xycon = ax1.contour(np.linspace(np.min(x),np.max(x),1000),np.linspace(np.min(y),np.max(y),1000),d2,[0.7],colors='k')

ax2 = fig.add_subplot(1,2,2,projection=ccrs.NorthPolarStereo())
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
ax2.set_boundary(circle, transform=ax2.transAxes)
ax2.set_extent([-180, 180,20, 90], ccrs.PlateCarree())
ax2.gridlines()
cyclic_data, cyclic_lons = cartopy.util.add_cyclic_point(d[100,:].data, coord = lons) ##
con1 = ax2.contourf(cyclic_lons, lats, cyclic_data, transform=ccrs.PlateCarree())
con = ax2.contour(cyclic_lons, lats, cyclic_data,[0.7],colors='k', transform=ccrs.PlateCarree())


plt.show()

xycon = xycon.allsegs[0][0]
con_map = con.allsegs[0][0]

lon, lat = pa(xycon[:,1],xycon[:,0],inverse=True)

print(con_map)

