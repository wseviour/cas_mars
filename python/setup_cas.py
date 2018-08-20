import numpy as np
import glob
from read_BOB import ds_from_BOB
import csv
import os
import matplotlib
matplotlib.use('qt5agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.util
import matplotlib.path as mpath
from math import radians, cos, sin, asin, sqrt, pi
from itertools import islice
from pyproj import Proj
from shapely.geometry import shape
from scipy.optimize import curve_fit


class Cas:
    """
    Perform CAS calculations
    """
    def __init__(self, run_dir, start_time, ndays, time_step):
        self.run_dir = run_dir
        self.start_time = start_time
        self.ndays = ndays
        self.time_step = time_step

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
            return rate
        except RuntimeError:
            print('no exponential fit')
            return 0.0

    def interpolate_winds(self):
        """
        Interpolate winds and write to text file for input into CAS routine
        """
        nsteps = int(self.ndays/self.time_step)

        # First check if wind files already exist
        file_list = [self.run_dir+'winds/u%.5d' % istep for istep in range(self.start_time,self.start_time+nsteps)]
        file_test = [os.path.isfile(f) for f in file_list]
        if all(file_test):
            print("winds already exist, no need to interpolate them")
            return
        else:
            print("not all wind file available, interpolating...")

            ds = ds_from_BOB(self.run_dir, ['u','v'], 85, time_step=self.time_step)

            if os.path.isdir(self.run_dir+'winds'):
                os.system('rm -f '+self.run_dir+'winds/*')
            else:
                os.system('mkdir '+self.run_dir+'winds')

            # Restrict to NH
            nlat = int(len(ds.coords['latitude'])/2)
            ds = ds.isel(latitude=slice(None,nlat))

            # Resolution for cas winds
            nlon = 144
            nlat = 92
            remainder = 4  # remainder after writing u 5 lon points at a time

            # Interpolate to new grid
            new_lon = np.arange(0,360,360./nlon)
            new_lat = np.arange(0, 90, 90./nlat)[::-1]
            ds_cas = ds.interp(latitude=new_lat,longitude=new_lon, kwargs={'fill_value':None})

            for istep in range(self.start_time,self.start_time+nsteps):
                with open(self.run_dir+'winds/u'+'%05d' % istep, 'w') as csvfile:
                    print('Writing '+self.run_dir+'winds/u'+'%05d' % istep)
                    writer = csv.writer(csvfile, delimiter=' ')
                    for ilat in range(nlat):
                        for iline in range(nlon//5):
                            writer.writerow(ds_cas['u'][istep,ilat,iline*5:5+iline*5].data)
                        writer.writerow(ds_cas['u'][istep,ilat,-remainder:].data)

            for istep in range(self.start_time,self.start_time+nsteps):
                with open(self.run_dir+'winds/v'+'%05d' % istep, 'w') as csvfile:
                    print('Writing '+self.run_dir+'winds/v'+'%05d' % istep)
                    writer = csv.writer(csvfile, delimiter=' ')
                    for ilat in range(nlat):
                        for iline in range(nlon//5):
                            writer.writerow(ds_cas['v'][istep,ilat,iline*5:5+iline*5].data)
                        writer.writerow(ds_cas['v'][istep,ilat,-remainder:].data)

    def make_contours(self, con_var='q',lats=np.arange(50,86,2)):
        """
        Produce contours for input into CAS
        """
        if os.path.isdir(self.run_dir+'contours'):
            try:
                os.system('rm -f '+self.run_dir+'contours/*.in')
            except OSError:
                pass
        else:
            os.system('mkdir '+self.run_dir+'contours')

        ds = ds_from_BOB(self.run_dir, [con_var], 85, time_step=self.time_step)

        ds_zm = ds[con_var].mean(dim='longitude')[self.start_time,:]

        lats = np.arange(50,86,2) # approx lats for PV contours

        cons = ds_zm.interp(latitude=lats).data
        #cons = [0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.35,1.3,1.2,1.1]

        lats = ds.coords['latitude'].data
        lons = ds.coords['longitude'].data
        count=0

        for icon in cons:

            inner = False
            if count > 0:
                if cons[count-1] > icon:
                    inner = True

            fig = plt.figure()
            ax = fig.add_subplot(1,1,1, projection=ccrs.NorthPolarStereo())
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)
            ax.set_boundary(circle, transform=ax.transAxes)
            ax.set_extent([-180, 180,20, 90], ccrs.PlateCarree())
            ax.gridlines()
            con1 = ax.contourf(lons, lats, ds[con_var][self.start_time,:,:].data,np.linspace(0,1.6,11),cmap='viridis', transform=ccrs.PlateCarree())
            con = ax.contour(lons, lats, ds[con_var][self.start_time,:,:].data,[icon],colors='k', transform=ccrs.PlateCarree())
            if len(con.allsegs[0]) == 1:
                a = con.allsegs[0][0]
            else:
                lens = np.zeros(len(con.allsegs[0]))
                for iicon in range(len(con.allsegs[0])):
                    lens[iicon] = self.calc_con_len(con.allsegs[0][iicon])
                if inner:
                    # 2nd longest contour
                    a = con.allsegs[0][np.where(lens == np.sort(lens)[-2])[0][0]]
                else:
                    a = con.allsegs[0][np.argmax(lens)]
            plt.plot(a[:,0],a[:,1], transform=ccrs.Geodetic(), color='red')
            plt.tight_layout()
            plt.close()

            if a[:,0][0] > a[:,1][1]:
                a = a[::-1,:]

            if inner:
                filename = self.run_dir+'contours/%s_%.4f_tstep_%s_inner.in' % (con_var,icon,self.start_time)
            else:
                filename = self.run_dir+'contours/%s_%.4f_tstep_%s.in' % (con_var,icon,self.start_time)

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
        cons = sorted([os.path.basename(x) for x in glob.glob(self.run_dir+"contours/*.in")])

        if os.path.isdir(self.run_dir+'cas_output'):
            try:
                os.system('rm -f '+self.run_dir+'cas_output/*.cas')
            except OSError:
                pass
        else:
            os.system('mkdir '+self.run_dir+'cas_output')

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
            os.system("cp -R %s/winds winds" % self.run_dir)
            os.system("cp %s/contours/%s con.in" % (self.run_dir,icon))
            os.system("./cas < run_params.in")
            os.system("cp output.cas %s/cas_output/%s.cas" % (self.run_dir, icon))

        os.chdir('../python')

    def read_cons_lens(self):
        """
        Read CAS output and return dictionary of contour coords and lengths
        """
        cas_files = sorted(glob.glob(self.run_dir+"cas_output/*.cas"))

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

            conname = ifile.split('/')[-1].split('_')[1][:6]

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



if __name__ == '__main__':

    PATH = '../model_output/'
    #ext = 'res_test57.-70.-nu4-urlx-kt4.0.c-0020sat200.0.T85/'
    ext = 'ann57.-70.-nu4-urlx-kt2.0-sinlat.c-0020.T85/'
    ndays = 20 # number days for CAS calculation
    start_time = 200 # starting time step of CAS
    time_step= 0.5 # time step in SWM
    res = 85 # resolution of SWM

    Cas = Cas(PATH+ext, start_time, ndays, time_step)
    Cas.interpolate_winds()
    Cas.make_contours(con_var='q')
    Cas.run_cas()
    len_rates = Cas.len_rate_eqlat()

    # interpolate_winds(PATH+ext, start_time, ndays, time_step)
    # make_contours(PATH+ext, start_time, ndays, con_var='q')
    # run_cas(PATH+ext,start_time)
    #
    # len_rates = len_rate_eqlat(PATH+ext, time_step)
