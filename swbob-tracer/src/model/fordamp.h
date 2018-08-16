#ifndef FORDAMP_H_
#define FORDAMP_H_

#include <dims.h>

      REAL_TYPE kvcoef(nlat2,0:1,plev,jelemd)
      REAL_TYPE ktcoef(nlat2,0:1,plev,jelemd)
      REAL_TYPE ftopo(nlon,nlat2,0:1,plev,ielemd,jelemd)
      REAL_TYPE ftopokt(nlon,nlat2,0:1,plev,ielemd,jelemd)
      REAL_TYPE fthr(nlon,nlat2,0:1,plev,ielemd,jelemd)
      REAL_TYPE fmechu(nlon,nlat2,0:1,plev,ielemd,jelemd)
      REAL_TYPE fmechv(nlon,nlat2,0:1,plev,ielemd,jelemd)
      REAL_TYPE fum1n(nlat2,0:1,plev,jelemd)
      REAL_TYPE fum1t(nlat2,0:1,plev,jelemd)

      common /forcing/ kvcoef,ktcoef,ftopo,ftopokt,fthr,fmechu,fmechv,fum1n,fum1t

#endif
