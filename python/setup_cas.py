import numpy as np
import glob
from read_BOB import ds_from_BOB
import csv
import os
import matplotlib
matplotlib.use('qt5agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import cartopy.crs as ccrs
import cartopy.util
import matplotlib.path as mpath
from math import radians, cos, sin, asin, sqrt, pi
from itertools import islice
from pyproj import Proj
from shapely.geometry import shape
from scipy.optimize import curve_fit
import xarray as xr


class Cas:
    """
    Perform CAS calculations
    """
    def __init__(self, ds, working_dir, start_time, ndays, time_step, model_type='SWM'):
        self.ds = ds
        self.working_dir = working_dir
        self.start_time = start_time
        self.ndays = ndays
        self.time_step = time_step
        self.model_type = model_type

    def exp_func(self, x, a, b):
        """
        Exponential function for fitting growth rate
        """
        return a * np.exp(b * x)

    def haversine(self, lon1, lat1, lon2, lat2):
        """
        Calculate the great circle distance between two points
        on a sphere (specified in decimal degrees)
        """
        # convert decimal degrees to radians
        lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

        # haversine formula
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
        c = 2 * asin(sqrt(a))
        r = 3389 # Radius of mars in kilometers.
        return c * r

    def calc_con_len(self, con):
        """
        Calculate the length of a contour on a sphere
        """
        con_len = 0
        for ipt in range(1,con.shape[0]-1):
            con_len += self.haversine(con[ipt-1,0],con[ipt-1,1],con[ipt,0],con[ipt,1])
        return con_len

    def eqlat_from_contour(self, conlats, conlons):
        """
        Calculate equivalent latitude of a contour
        """
        pa = Proj("+proj=stere +lat_0=90")
        x, y = pa(conlons, conlats)
        cop = {"type": "Polygon", "coordinates": [zip(x, y)]}
        area = shape(cop).area
        radius = 6371000 #Earth radius for use with Proj
        eqlat = np.rad2deg(asin(1-(area/(2*pi*radius**2))))
        return eqlat

    def calc_lengthening_rate(self, lens, time_step):
        """
        Given length time series, return lengthening rate
        """
        try:
            popt,pcov = curve_fit(self.exp_func, np.arange(len(lens[:30])), lens[:30], p0=(0.1,0.01))
            rate = popt[1]/time_step
            if rate < 0:
                rate = 0.0
            return rate
        except RuntimeError:
            print('no exponential fit')
            return 0.0

    def interpolate_winds(self):
        """
        Interpolate winds and write to text file for input into CAS routine
        """
        nsteps = int(self.ndays/self.time_step)

        # Resolution for cas winds
        nlon = 144
        nlat = 92
        remainder = 4  # remainder after writing u 5 lon points at a time
        
        # First check if wind files already exist
        file_list = [self.working_dir+'winds/u%.5d' % istep for istep in range(self.start_time,self.start_time+nsteps)]
        file_test = [os.path.isfile(f) for f in file_list]
        if all(file_test):
            print("winds already exist, no need to interpolate them")
            return
        else:
            print("not all wind file available, interpolating...")

            if os.path.isdir(self.working_dir+'winds'):
                os.system('rm -f '+self.working_dir+'winds/*')
            else:
                os.system('mkdir '+self.working_dir+'winds')

            # Interpolate to new grid for CAS
            new_lon = np.arange(0,360,360./nlon)
            new_lat = np.arange(0, 90, 90./nlat)[::-1]
            ds_cas = self.ds.interp(latitude=new_lat,longitude=new_lon, kwargs={'fill_value':None})

            for istep in range(self.start_time,self.start_time+nsteps):
                with open(self.working_dir+'winds/u'+'%05d' % istep, 'w') as csvfile:
                    print('Writing '+self.working_dir+'winds/u'+'%05d' % istep)
                    writer = csv.writer(csvfile, delimiter=' ')
                    for ilat in range(nlat):
                        for iline in range(nlon//5):
                            writer.writerow(ds_cas['u'][istep,ilat,iline*5:5+iline*5].data)
                        writer.writerow(ds_cas['u'][istep,ilat,-remainder:].data)

            for istep in range(self.start_time,self.start_time+nsteps):
                with open(self.working_dir+'winds/v'+'%05d' % istep, 'w') as csvfile:
                    print('Writing '+self.working_dir+'winds/v'+'%05d' % istep)
                    writer = csv.writer(csvfile, delimiter=' ')
                    for ilat in range(nlat):
                        for iline in range(nlon//5):
                            writer.writerow(ds_cas['v'][istep,ilat,iline*5:5+iline*5].data)
                        writer.writerow(ds_cas['v'][istep,ilat,-remainder:].data)


    def make_contours(self, con_var='q', lats=np.arange(50,86,2), plot=False):

        if os.path.isdir(self.working_dir+'contours'):
            try:
                os.system('rm -f '+self.working_dir+'contours/*.in')
            except OSError:
                pass
        else:
            os.system('mkdir '+self.working_dir+'contours')
        
        # Only use 90 - 20 latitude
        d = self.ds[con_var].sel(latitude = slice(90,20))[self.start_time,:]
        
        #cons = d.mean(dim = 'longitude').interp(latitude=lats).data
        cons = d.sel({'latitude':lats,'longitude':0}, method='nearest').data
        
        lats = d.coords['latitude'].data
        lons = d.coords['longitude'].data
        count=0

        for icon in cons:
            print('contour: '+str(icon))
            inner = False
            if (count > 0) and (cons[count-1] > icon):
                inner = True
                print('inner')
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1, projection=ccrs.NorthPolarStereo())
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)
            ax.set_boundary(circle, transform=ax.transAxes)
            ax.set_extent([-180, 180,20, 90], ccrs.PlateCarree())
            ax.gridlines()
            cyclic_data, cyclic_lons = cartopy.util.add_cyclic_point(d.data, coord = lons) ##
            con1 = ax.contourf(cyclic_lons, lats, cyclic_data,cmap='viridis', transform=ccrs.PlateCarree())
            con = ax.contour(cyclic_lons, lats, cyclic_data,[icon],colors='k', transform=ccrs.PlateCarree())

            fig2 = plt.figure()
            p = d.roll(longitude=0).plot.contour(levels=[icon])
            plt.close()
            
            if len(p.allsegs[0]) == 1:
                a = p.allsegs[0][0]
            else:
                lens = [p.allsegs[0][i].shape[0] for i in range(len(p.allsegs[0]))] 
                if inner:
                    # 2nd longest contour
                    a = p.allsegs[0][np.where(lens == np.sort(lens)[-2])[0][0]]
                else:
                    a = p.allsegs[0][np.argmax(lens)]
                #a = a[1:,:]
            a = a[a[:,0]<360]
            ax.plot(a[:,0],a[:,1], transform=ccrs.Geodetic(), color='red')
            plt.tight_layout()
            if plot:
                plt.show()
            plt.close()

            if a[:,0][0] > a[:,0][1]:
                a = a[::-1,:]

           
            
            if inner:
                filename = self.working_dir+'contours/%s_%.4f_tstep_%s_inner.in' % (con_var,icon,self.start_time)
            else:
                filename = self.working_dir+'contours/%s_%.4f_tstep_%s.in' % (con_var,icon,self.start_time)

            with open(filename, "w") as csvfile:
                csvfile.write("Contour Advection with Surgery\n")
                csvfile.write("%s %.4f contour\n" % (con_var,icon))
                csvfile.write("\n")
                csvfile.write("%s  24  %.7f  %.7f  0.1000000  0.0000000\n" % (self.ndays,self.time_step,self.time_step))
                csvfile.write("1 %s 0.00000\n" % a.shape[0])
                csvfile.write("%s %d %d 1.00000\n" % (a.shape[0], a[0,0], a[0,1]))

            with open(filename, "a") as csvfile:
                writer = csv.writer(csvfile, delimiter=' ')
                for irow in range(a.shape[0]):
                    writer.writerow(a[irow,:])

            count +=1



    def run_cas(self):
        """
        Run CAS over all contours
        """
        cons = sorted([os.path.basename(x) for x in glob.glob(self.working_dir+"contours/*.in")])

        if os.path.isdir(self.working_dir+'cas_output'):
            try:
                os.system('rm -f '+self.working_dir+'cas_output/*.cas')
            except OSError:
                pass
        else:
            os.system('mkdir '+self.working_dir+'cas_output')

        os.chdir('../cas')

        with open('run_params.in', "w") as inputfile:
            inputfile.write("con.in\n")
            inputfile.write("output.cas\n")
            inputfile.write("output.log\n")
            inputfile.write(str(self.start_time))
            inputfile.write("\n")

        for icon in cons:
            print("**** %s ****" % icon)
            os.system("rm -r winds")
            os.system("cp -R %s/winds winds" % self.working_dir)
            os.system("cp %s/contours/%s con.in" % (self.working_dir,icon))
            os.system("./cas < run_params.in")
            os.system("cp output.cas %s/cas_output/%s.cas" % (self.working_dir, icon))

        os.chdir('../python')

    def read_cons_lens(self):
        """
        Read CAS output and return dictionary of contour coords and lengths
        """
        cas_files = sorted(glob.glob(self.working_dir+"cas_output/*.cas"))

        cons_dic = {}
        lens_dic = {}

        for ifile in cas_files:

            f = open(ifile,'r')
            lines = f.readlines()

            params = np.fromstring(lines[3], sep=' ')
            ndays = params[0] # number of days
            ntstep = params[1] # timesteps per day
            tgrid = params[2] # time between gridded data (days)
            tsave = params[3] # time between output data (days)
            amu = params[4]
            dm = params[5]

            tot_params = np.fromstring(lines[4], sep=' ')
            ncont = tot_params[0] # number of contours
            npts_tot = tot_params[1]

            con_params = np.empty(0)
            for iline in range(5,5+int(ncont)):
                con_params = np.append(con_params, np.fromstring(lines[iline], sep=' '))

            cons = []
            con = np.zeros((int(con_params[0]),2))
            for iline in range(int(5+ncont), int(5+ncont+con_params[0])):
                con[iline-int(5+ncont),:] = np.fromstring(lines[iline], sep=' ')
            cons.append(con)

            f.close()

            for istep in range(1,int(ndays/tsave)):
                try:
                    f = open(ifile,'r')
                    lines = f.readlines()

                    tot_params = np.fromstring(lines[int(5+(ncont*istep)+npts_tot+istep-1)], sep=' ')
                    npts_tot2 = tot_params[1]

                    con_params = np.empty(0)
                    for iline in range(int(5+1+(ncont*istep)+npts_tot+istep-1), int(5+1+ncont+(ncont*istep)+npts_tot+istep-1)):
                        con_params = np.append(con_params, np.fromstring(lines[iline], sep=' '))

                    con = np.zeros((int(con_params[0]),2))
                    for iline in range(int(5+1+ncont+(ncont*istep)+npts_tot+istep-1), int(5+1+ncont+(ncont*istep)+npts_tot+istep-1+con_params[0])):
                        con[iline-int(5+1+ncont+(ncont*istep)+npts_tot+istep-1),:] = np.fromstring(lines[iline], sep=' ')
                    cons.append(con)

                    npts_tot = npts_tot+npts_tot2
                    f.close()

                except IndexError:
                    f.close()

            lens = []
            for icon in cons:
                lens.append(self.calc_con_len(icon))

            conname = ifile.split('/')[-1]

            cons_dic[conname] = cons
            lens_dic[conname] = lens

        return cons_dic, lens_dic


    def len_rate_eqlat(self):
        """
        Retuns array of equivalent latitudes and their lengthening rates
        """
        cons, lens = self.read_cons_lens()
        len_rates = np.zeros((len(cons),2))
        count=0
        for icon in cons.keys():
            len_rates[count,0] = self.eqlat_from_contour(cons[icon][0][:,1], cons[icon][0][:,0])
            len_rates[count,1] = self.calc_lengthening_rate(lens[icon], self.time_step)
            count+=1
        len_rates = len_rates[len_rates[:,0].argsort()] # sort by eqlat

        return len_rates

    
    def plot_cas(self, plot_day):

        cons, lens = self.read_cons_lens()
        lats = self.ds.coords['latitude'].data
        lons = self.ds.coords['longitude'].data
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection=ccrs.NorthPolarStereo())
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax.set_boundary(circle, transform=ax.transAxes)
        ax.set_extent([-180, 180,30, 90], ccrs.PlateCarree())
        ax.gridlines()
        cyclic_data, cyclic_lons = cartopy.util.add_cyclic_point(self.ds.q[start_time+int(plot_day/self.time_step),:].data, coord = lons)
        pv = ax.contourf(cyclic_lons, lats,cyclic_data, cmap='Blues', transform=ccrs.PlateCarree(),extend='both')
        for key in cons.keys():
            con = cons[key][int(plot_day/time_step)]
            ax.plot(con[:,0], con[:,1], transform=ccrs.Geodetic(),color='k')
        ax.set_title('day '+str(plot_day))
        plt.show()
        
        
def read_MarsWRF(level=0):
    """
    Read MarsWRF data at a given theta level, return a xarray Dataset
    """
    ds_raw = xr.open_mfdataset('../MarsWRF_data/*.nc', concat_dim='Time')
    ds_raw = ds_raw.where(ds_raw.L_S < 300, drop=True) # Drop missing values in MarsWRF files
    lat = ds_raw.LAT[0,:,0]
    lon = ds_raw.LON[0,0,:] % 360 # make sure in 0-360 format, not -180-180
    lon = np.append(lon[90:],lon[:90])
    theta = ds_raw.THETA[0,:]
    l_s = ds_raw.L_S

    print('Reading MarsWRF data at theta = '+str(theta[level].data))
    
    U = ds_raw.U.data
    V = ds_raw.V.data
    EPV = ds_raw.EPV.data
    U = np.append(U[:,:,:,90:],U[:,:,:,:90], axis=3)
    V = np.append(V[:,:,:,90:],V[:,:,:,:90], axis=3)
    EPV = np.append(EPV[:,:,:,90:],EPV[:,:,:,:90], axis=3)            

    ds = xr.Dataset({'u':(['time','theta','latitude','longitude'], U),
                     'v':(['time','theta','latitude','longitude'], V),
                     'q':(['time','theta','latitude','longitude'], EPV)},
                     coords = {'time':('time',l_s),
                               'theta':('theta',theta),
                               'latitude':('latitude',lat),
                               'longitude':('longitude',lon)})

    ds = ds.isel(theta=level) # Select given level
    ds = ds.sel(latitude=slice(0,90))

    return ds


if __name__ == '__main__':


    ndays = 10 # number days for CAS calculation
    start_time = 200 # starting time step of CAS
    time_step= 0.25 # time step in SWM


    ############# SWM test ############
    run_name = 'ann57.-70.-nu4-urlx-kt2.0-ring.c-0020.T85.nc'

    ds = xr.open_dataset('../model_output/netcdf/'+run_name)

    working_dir = '../model_output/cas/'+run_name+'/'
    if os.path.isdir(working_dir):
        pass
    else:
        os.system('mkdir '+working_dir)

    ############# MarsWRF test ###############
    # level = 10
    # ds = read_MarsWRF(level=level)

    # working_dir = '../MarsWRF_data/cas/Z%.2d/' % level
    # if os.path.isdir(working_dir):
    #     pass
    # else:
    #     os.system('mkdir '+working_dir)  
        

    CA = Cas(ds, working_dir, start_time, ndays, time_step)
    #CA.interpolate_winds()
    CA.make_contours(lats=np.arange(50,86,2),plot=False)
    CA.run_cas()

    CA.plot_cas(0)
    # Cas = Cas(PATH+ext, start_time, ndays, time_step)
    # Cas.interpolate_winds()
    # Cas.make_contours(con_var='q',plot=True)
    # Cas.run_cas()
    # len_rates = Cas.len_rate_eqlat()

    #Cas = Cas('../MarsWRF_data/',100,10,0.25,model_type='MarsWRF',level=10)
    #Cas.interpolate_winds()
    #Cas.make_contours(con_var='q',lats=[65],plot=True)
    #Cas.run_cas()

    
    # interpolate_winds(PATH+ext, start_time, ndays, time_step)
    # make_contours(PATH+ext, start_time, ndays, con_var='q')
    # run_cas(PATH+ext,start_time)
    #
    # len_rates = len_rate_eqlat(PATH+ext, time_step)


        # def make_contours(self, con_var='q', lats=np.arange(50,86,2), plot=False):

        # if os.path.isdir(self.working_dir+'contours'):
        #     try:
        #         os.system('rm -f '+self.working_dir+'contours/*.in')
        #     except OSError:
        #         pass
        # else:
        #     os.system('mkdir '+self.working_dir+'contours')
        
        # # Only use 90 - 20 latitude
        # d = self.ds[con_var].sel(latitude = slice(90,20))[self.start_time,:]
        
        # # Find contour levels by interpolation at lon=0
        # cons = d.mean(dim = 'longitude').interp(latitude=lats).data

        # pa = Proj("+proj=stere +lat_0=90",preserve_units=True)
        # lonv, latv = np.meshgrid(d.longitude.data, d.latitude.data)
        # x, y = pa(lonv,latv)
        # reg_x = np.linspace(np.min(x),np.max(x),500)
        # reg_y = np.linspace(np.min(y),np.max(y),500)    
        # xi, yi = np.meshgrid(reg_x, reg_y)
        # d2 = mlab.griddata(x.flatten(),y.flatten(),d.data.flatten(),xi,yi,interp='linear')

        # count = 0
        # for icon in cons:

        #     inner = False
        #     if count > 0:
        #         if cons[count-1] > icon:
        #             inner = True
            
        #     fig = plt.figure(figsize=(10,5))
        #     ax1 = fig.add_subplot(1,2,1)
        #     ax1.contourf(reg_x,reg_y,d2)
        #     xycon = ax1.contour(reg_x,reg_y,d2,[icon],colors='k')
            
        #     latloncon = []
        #     for icon in range(len(xycon.allsegs[0])):
        #         ilons, ilats = pa(xycon.allsegs[0][icon][:,0],xycon.allsegs[0][icon][:,1],inverse=True)
        #         ilons = ilons % 360
        #         latloncon.append(np.vstack((ilons,ilats)).T)
                
            
        #     if len(latloncon) == 1:
        #         a = latloncon[0]
        #     else:
        #         print('more than one contour')
        #         lens = np.zeros(len(latloncon))
        #         for iicon in range(len(latloncon)):
        #             lens[iicon] = self.calc_con_len(latloncon[iicon])
        #         if inner:
        #             # 2nd longest contour
        #             a = latloncon[np.where(lens == np.sort(lens)[-2])[0][0]]
        #         else:
        #             a = latloncon[np.argmax(lens)]

        #     a = a[:-1,:]
        #     #lon, lat = pa(a[:,0],a[:,1],inverse=True)
        #     #lon = lon % 360

        #     #a[:,0] = lon
        #     #a[:,1] = lat
        #     # if a[0,0] > a[1,0]:
        #     #     a = a[::-1,:]
            
        #     if plot:
        #         lats = d.coords['latitude'].data
        #         lons = d.coords['longitude'].data
        #         #ax1.plot(a[:,0],a[:,1],color='red')
        #         ax2 = fig.add_subplot(1,2,2,projection=ccrs.NorthPolarStereo())
        #         theta = np.linspace(0, 2*np.pi, 100)
        #         center, radius = [0.5, 0.5], 0.5
        #         verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        #         circle = mpath.Path(verts * radius + center)
        #         ax2.set_boundary(circle, transform=ax2.transAxes)
        #         ax2.set_extent([-180, 180,20, 90], ccrs.PlateCarree())
        #         ax2.gridlines()
        #         cyclic_data, cyclic_lons = cartopy.util.add_cyclic_point(d.data, coord = lons) ##
        #         con1 = ax2.contourf(cyclic_lons, lats, cyclic_data, transform=ccrs.PlateCarree())
        #         con = ax2.contour(cyclic_lons, lats, cyclic_data,[icon],colors='k', transform=ccrs.PlateCarree())
        #         ax2.plot(a[:,0],a[:,1], transform=ccrs.Geodetic(), color='red')
        #         plt.show()
        #     else:
        #         plt.close()

        #     if inner:
        #         filename = self.working_dir+'contours/%s_%.7f_tstep_%s_inner.in' % (con_var,icon,self.start_time)
        #     else:
        #         filename = self.working_dir+'contours/%s_%.7f_tstep_%s.in' % (con_var,icon,self.start_time)

        #     with open(filename, "w") as csvfile:
        #         csvfile.write("Contour Advection with Surgery\n")
        #         csvfile.write("%s %.4f contour\n" % (con_var,icon))
        #         csvfile.write("\n")
        #         csvfile.write("%s  24  %.7f  %.7f  0.10000000  0.0000000\n" % (self.ndays,self.time_step,self.time_step))
        #         csvfile.write("1 %s 0.00000\n" % a.shape[0])
        #         csvfile.write("%s %d %d 1.00000\n" % (a.shape[0], a[0,0], a[0,1]))

        #     with open(filename, "a") as csvfile:
        #         writer = csv.writer(csvfile, delimiter=' ')
        #         for irow in range(a.shape[0]):
        #             writer.writerow(a[irow,:])
            
        #     count += 1
