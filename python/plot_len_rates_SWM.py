from setup_cas import Cas
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
import save_restore

ndays = 10 # number days for CAS calculation
#start_time = 250 # starting time step of CAS
time_step= 0.25 # time step in model (should be 0.25 for MarsWRF)

lin_lats = np.arange(45,86,2)

kt = ['0.0','0.4','1.0','2.0']


#start_times = range(500,520,8)
start_times = [500]

lin_len_rates = np.zeros((len(kt),len(start_times),len(lin_lats))) 



for ikt in kt:
    for start_time in start_times:
        FILE = '/home/bridge/kz18101/Research/cas_mars/model_output/netcdf/ann57.-70.-nu4-urlx-kt%s-sinlat-trfac1.c-0020sat200.0.T170.nc' % ikt
        run_name = FILE.split('/')[-1]
        
        ds = xr.open_dataset(FILE)
    
        working_dir = '../model_output/cas/'+run_name+'/'
        if os.path.isdir(working_dir):
            pass
        else:
            os.system('mkdir '+working_dir)

        CA = Cas(ds, working_dir, start_time, ndays, time_step)
        if start_times.index(start_time) == 0:
            CA.interpolate_winds(length=start_times[-1]-start_times[0]+40)
        CA.make_contours2(lats=np.arange(50,86,2),plot=False)
        CA.run_cas()
        len_rates = CA.len_rate_eqlat()
        #CA.plot_cas(0)

        lin_len_rates[kt.index(ikt),start_times.index(start_time),:] = np.interp(lin_lats,len_rates[:,0],len_rates[:,1])


save_restore.save('SWM_CAS.pypic', lin_len_rates=lin_len_rates, lin_lats=lin_lats, kt =kt)
