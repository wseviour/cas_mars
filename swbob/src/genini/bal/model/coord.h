#ifndef COORD_H_
#define COORD_H_

#include <dims.h>

      REAL_TYPE Alpha_coord           ! angle between coordinate and rotational poles
      parameter(Alpha_coord=0.0)    

      INT_TYPE  plonmax(jelemg)       ! max longitude for each tile

      REAL_TYPE colatitude(plat/2)    ! colatitudes 
      REAL_TYPE gausswt(plat/2)       ! gaussian weights 

      REAL_TYPE coslat(plat)          ! cosine of latitudes

      common /coord/ plonmax,colatitude,gausswt,coslat

#endif
