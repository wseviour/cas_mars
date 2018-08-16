#ifndef MODELTIME_H_
#define MODELTIME_H_

#include <type.h>

      INT_TYPE  ntimestep        ! number of timesteps to integrate
      INT_TYPE  nspday           ! number of timesteps per day
      INT_TYPE  itime           ! keeps track of the numer of steps
      INT_TYPE  itime0           ! starting time index (0 or 1)

      REAL_TYPE time_current     ! the current  model time
      REAL_TYPE time_start       ! the starting model time for this run
      REAL_TYPE time_stop        ! the stopping model time for this run
      REAL_TYPE time_advance     ! the advance of model time for this run
      REAL_TYPE time_out         ! interval between output of 2d fields
      REAL_TYPE time_zon         ! interval between output of zonal fields

      common /modeltime/time_current, time_start, time_stop, time_advance,
     $                  time_out, time_zon, 
     $                  ntimestep, nspday, itime, itime0

#endif
