       subroutine sc_xfer_layer(sc3d,sc2d,m,nlev,k)
       implicit none

#include <dims.h>


       INT_TYPE  m
       INT_TYPE  nlev
       REAL_TYPE sc3d(0:1,m:nn,nlev)
       REAL_TYPE sc2d(0:1,m:nn)
       INT_TYPE  k

       INT_TYPE n

       do n=m,nn
          sc2d(0,n)=sc3d(0,n,k)
          sc2d(1,n)=sc3d(1,n,k)
       end do

       end
