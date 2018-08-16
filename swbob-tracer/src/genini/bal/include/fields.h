#ifndef FIELDS_H_
#define FIELDS_H_

#include <dims.h>

      REAL_TYPE nlt(nlon,nlat2,2*nnlt*plev,ielemd,jelemd)  ! physical space non linear terms
      REAL_TYPE phy(nlon,nlat2,2*nphy*plev,ielemd,jelemd)  ! physical space state variables
      REAL_TYPE vomega(nlon,nlat2,2*(plev+1),ielemd,jelemd) ! vertical velocity

      REAL_TYPE afc(abuflen,0:M_NODE-1)        ! fourier coefficients of non linear terms
      REAL_TYPE afcx(abuflen,0:P_NODE-1)       ! transposed fourier coefficients of non linear terms
      REAL_TYPE sfc(sbuflen,0:M_NODE-1)        ! fourier coefficients of state variables
      REAL_TYPE sfcx(sbuflen,0:P_NODE-1)       ! transposed fourier coefficients of state variables

      common /fields/ nlt,phy,vomega,afc,afcx,sfc,sfcx

#endif
