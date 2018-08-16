#ifndef PHYSICAL_H_
#define PHYSICAL_H_

#include <dims.h>

      INTEGER   ielem(jelemg)       ! number of longitude elements for a given latitude chunk
      REAL_TYPE lon(nlon,ielemd,jelemd)  ! nodal longitudes in rads
      REAL_TYPE lat  (nlat2,jelemg) ! latitudes of gaussian grids
      REAL_TYPE cslat(nlat2,jelemg) ! cosine of latitude
      REAL_TYPE snlat(nlat2,jelemg) ! sine of latitude
      REAL_TYPE da   (nlat2,jelemg) ! area weight

      common /physical/ ielem,lat,cslat,snlat,da,lon
      
#endif
