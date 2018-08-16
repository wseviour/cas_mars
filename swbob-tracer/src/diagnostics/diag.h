#ifndef DIAG_H_
#define DIAG_H_

#include <dims.h>

      INT_TYPE nzavg
      parameter(nzavg=7)

      ! zonally averaged fields

      REAL_TYPE zon(plat,plev,nzavg) ! instantaneous zonal averages
      REAL_TYPE zon_node_acc(nlat2*2*plev*nzavg,jelemd)

      common /zonalavg/ zon,zon_node_acc

#endif
