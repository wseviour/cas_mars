#ifndef FORDAMP_H_
#define FORDAMP_H_

#include <dims.h>

      REAL_TYPE kvcoef(nlat2,0:1,plev,jelemd)
      REAL_TYPE ktcoef(nlat2,0:1,plev,jelemd)
      REAL_TYPE phitopo(nlon,nlat2,0:1,plev,ielemd,jelemd)

      common /forcing/ kvcoef,ktcoef,phitopo

#endif
