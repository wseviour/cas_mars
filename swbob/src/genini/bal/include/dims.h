#ifndef DIMS_H_
#define DIMS_H_

#include <type.h>
#include <params.h>

      INT_TYPE plon,plat,plev

      parameter(plon=PLON)
      parameter(plev=PLEV)
      parameter(plat=plon/2)

      INT_TYPE plat2d
      parameter(plat2d=((plat/2)-1)/P_NODE + 1)

      INT_TYPE nlon,ielemd
      parameter(nlon=8,ielemd=(plon/nlon))      ! Note: assumes nlon evenly divides plon

      INT_TYPE nlat2,jelemg,jelemd
      parameter(nlat2=4,jelemg=(plat/2)/nlat2)  ! Note: plat/2 assumed divisible by nlat2
      parameter(jelemd=(jelemg-1)/P_NODE + 1)

      INT_TYPE  mm,nn
      INT_TYPE mm2d

      parameter(mm=(plon-1)/3,nn=mm)
      parameter(mm2d=mm/(2*M_NODE)+1)

      INT_TYPE nscg  ! global spectral coefficient counts
      parameter(nscg=(nn+1)*(mm+1)-((mm+1)*mm)/2)

      INT_TYPE nvsc,nsc  ! local paired coefficient counts
      parameter(nvsc=mm2d*(nn+4),nsc=mm2d*(nn+2))

      INT_TYPE nphy,nnlt,nsp
      parameter(nphy=6,nnlt=5,nsp=3)

      INT_TYPE nsyn,nanal
      parameter(nsyn=6,nanal=5)

      INT_TYPE sbuflen
      parameter(sbuflen=2*nlat2*nsyn*2*plev*mm2d*2*jelemd)

      INT_TYPE abuflen
      parameter(abuflen=2*nlat2*nanal*2*plev*mm2d*2*jelemd)


#endif


