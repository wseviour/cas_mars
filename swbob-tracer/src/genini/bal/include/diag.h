#ifndef DIAG_H_
#define DIAG_H_

#include <dims.h>

      INT_TYPE nzavg
      parameter(nzavg=6)

      ! zonally averaged fields

      REAL_TYPE zon(plat,plev,nzavg) ! instantaneous zonal averages

      common /zonalavg/ zon

#endif
