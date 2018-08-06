'''
Script to make plot for Darryn's proposal
'''

import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from itertools import islice
import matplotlib.path as mpath
import glob
from read_BOB import cube_from_BOB
import cartopy.util
from math import radians, cos, sin, asin, sqrt
from scipy.optimize import curve_fit
import iris

def func(x, a, b, c):
    return a * np.exp(b * x) + c

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    return c * r

def calc_con_len(con):
    con_len = 0
    for ipt in range(1,con.shape[0]-1):
        con_len += haversine(con[ipt-1,0],con[ipt-1,1],con[ipt,0],con[ipt,1])
    return con_len

# file_circ1 = '../sandbox/2deg/circle_sat0_kt0_pv0.9.cas'
# file_circ2 = '../sandbox/2deg/circle_sat0_kt0_pv1.4.cas'
#
# file_ann1 = '../sandbox/2deg/circle_sat0_kt1_pv0.9.cas'
# file_ann2 = '../sandbox/2deg/circle_sat0_kt1_pv1.4.cas'

PATH = '../MarsWRF_data/'
#file = 'isentropic_surface_data_ls=215_to_ls=240.nc'
file = 'isentropic_surface_data_ls=240_to_ls=266.nc'
cube = iris.load_cube(PATH+file, 'EPV')
L_S = iris.load_cube(PATH+file, 'L_S')
lat = iris.load_cube(PATH+file, 'LAT')
lon = iris.load_cube(PATH+file, 'LON')
theta = iris.load_cube(PATH+file, 'THETA')
lon = lon.data[0,:]
for ilon in range(lon.shape[0]):
    if lon[ilon] < 0:
        lon[ilon] += 360
lon = np.append(lon[90:],lon[:90])
cube = cube.data[:50,:]
cube = np.append(cube[:,:,:,90:],cube[:,:,:,:90], axis=3)
L_S = iris.coords.DimCoord(L_S.data[:50], standard_name='angle_of_incidence', units='degree', var_name='L_S')
lon = iris.coords.DimCoord(lon,standard_name='longitude',units='degrees_E',var_name='lon',circular=True)
lat = iris.coords.DimCoord(lat.data[:,0],standard_name='latitude',units='degrees_N', var_name='lat')
theta = iris.coords.DimCoord(theta.data,standard_name='air_potential_temperature', units='K', var_name='T')
cube = iris.cube.Cube(cube,long_name='Ertel Potential Voriticity',var_name='EPV',
                    dim_coords_and_dims=[(L_S,0),(theta,1),(lat,2),(lon,3)])
theta_lev = 300
level = iris.Constraint(air_potential_temperature = theta_lev)
cube = cube.extract(level)
lats = cube.coord('latitude').points
lons = cube.coord('longitude').points

fileskt1 = sorted(glob.glob("../cas/cas_output_MarsWRF/%sK/*%s*.cas" % (theta_lev,file.split('.nc')[0])))

nsteps=60
factor = 1
initday = 10

fig = plt.figure(figsize=(12,4))

rows = [4,5,6,7]


cons_dic = {}
lens_dic = {}


for ifile in fileskt1:

    with open(ifile) as lines:
        params = np.genfromtxt(islice(lines, 3, 4))
    ndays = params[0] # number of days
    ntstep = params[1] # timesteps per day
    tgrid = params[2] # time between gridded data (days)
    tsave = params[3] # time between output data (days)
    amu = params[4]
    dm = params[5]

    with open(ifile) as lines:
        tot_params = np.genfromtxt(islice(lines, 4, 5))
    ncont = tot_params[0] # number of contours
    npts_tot = tot_params[1]

    with open(ifile) as lines:
        con_params = np.genfromtxt(islice(lines, 5, 5+ncont))

    cons = []
    with open(ifile) as lines:
        con = np.genfromtxt(islice(lines, 5+ncont, 5+ncont+con_params[0]))
    cons.append(con)


    for istep in range(1,nsteps):

        with open(ifile) as lines:
            tot_params = np.genfromtxt(islice(lines, 5+(ncont*istep)+npts_tot+istep-1, 5+(ncont*istep)+npts_tot+(istep-1)+1))
        npts_tot2 = tot_params[1]

        with open(ifile) as lines:
            con_params = np.genfromtxt(islice(lines, 5+1+(ncont*istep)+npts_tot+istep-1, 5+1+ncont+(ncont*istep)+npts_tot+istep-1))

        with open(ifile) as lines:
            con = np.genfromtxt(islice(lines, 5+1+ncont+(ncont*istep)+npts_tot+istep-1, 5+1+ncont+(ncont*istep)+npts_tot+istep-1+con_params[0]))
        cons.append(con)

        npts_tot = npts_tot+npts_tot2

    lens = []
    for icon in cons:
        lens.append(calc_con_len(icon))


    conname = ifile.split('_tstep')[0]
    print conname
    if 'inner' in ifile:
        conname = 'inner'+conname

    cons_dic[conname] = cons
    lens_dic[conname] = lens


keys = sorted(cons_dic.keys())
indices = [x for x, s in enumerate(keys) if 'inner' in s]
if len(indices) > 0:
    keys[indices[0]:] = sorted(keys[indices[0]:],reverse=True)

colors = plt.cm.rainbow(np.linspace(0,1,len(keys)))

ax = fig.add_subplot(1,4,1, projection=ccrs.NorthPolarStereo())
ang = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(ang), np.cos(ang)]).T
circle = mpath.Path(verts * radius + center)
ax.set_boundary(circle, transform=ax.transAxes)
ax.set_extent([-180, 180,45, 90], ccrs.PlateCarree())
ax.gridlines()
for key in keys:
    con = cons_dic[key][0]
    con = np.append(con[-1,:][np.newaxis,:],con,axis=0)
    ax.plot(con[:,0], con[:,1], transform=ccrs.Geodetic(),color=colors[keys.index(key)], label=key)

cyclic_data, cyclic_lons = cartopy.util.add_cyclic_point(cube[initday,:,:].data, coord=lons)
con = ax.contourf(cyclic_lons, lats, cyclic_data*factor,np.linspace(0,0.003,30),cmap='gray_r',extend='both', transform=ccrs.PlateCarree())
ax.set_title('initial')

ax = fig.add_subplot(1,4,2, projection=ccrs.NorthPolarStereo())
ang = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(ang), np.cos(ang)]).T
circle = mpath.Path(verts * radius + center)
ax.set_boundary(circle, transform=ax.transAxes)
ax.set_extent([-180, 180,45, 90], ccrs.PlateCarree())
ax.gridlines()
for key in keys:
    con = cons_dic[key][12]
    if (con.shape[0] > 0):
        if (con[:,0].shape[0] > 20):
            con = np.append(con[-1,:][np.newaxis,:],con,axis=0)
            ax.plot(con[:,0], con[:,1], transform=ccrs.Geodetic(),color=colors[keys.index(key)])

cyclic_data, cyclic_lons = cartopy.util.add_cyclic_point(cube[initday+12,:,:].data, coord=lons)
con = ax.contourf(cyclic_lons, lats, cyclic_data*factor,np.linspace(0,0.003,30),cmap='gray_r',extend='both', transform=ccrs.PlateCarree())
ax.set_title('Day 3')
#ax.text(-0.09,0.55, 'Day 5', transform=ax.transAxes, fontsize=14,rotation=90)


ax = fig.add_subplot(1,4,3, projection=ccrs.NorthPolarStereo())
ang = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(ang), np.cos(ang)]).T
circle = mpath.Path(verts * radius + center)
ax.set_boundary(circle, transform=ax.transAxes)
ax.set_extent([-180, 180,45, 90], ccrs.PlateCarree())
ax.gridlines()

for key in keys:
    con = cons_dic[key][28]
    con = np.append(con[-1,:][np.newaxis,:],con,axis=0)
    ax.plot(con[:,0], con[:,1], transform=ccrs.Geodetic(),color=colors[keys.index(key)])
cyclic_data, cyclic_lons = cartopy.util.add_cyclic_point(cube[initday+28,:,:].data, coord=lons)
con = ax.contourf(cyclic_lons, lats, cyclic_data*factor,np.linspace(0,0.003,30),cmap='gray_r',extend='both', transform=ccrs.PlateCarree())
ax.set_title('Day 7')

#ax.text(-0.09,0.55, 'Day 10', transform=ax.transAxes, fontsize=14,rotation=90)


ax = fig.add_subplot(1,4,4)
for key in keys:
    ax.plot(np.arange(0,len(lens_dic[key])/4.,0.25),np.array(lens_dic[key])/lens_dic[key][0],color=colors[keys.index(key)])
ax.set_ylabel('length/original length')
ax.set_xlabel('days')
ax.set_ylim(0,30)
ax.set_xlim(0,15)

plt.tight_layout()
#plt.subplots_adjust(wspace=0.1, hspace=0.1)
plt.savefig('../plots/'+file+'.pdf', bbox_inches='tight')
plt.show()
