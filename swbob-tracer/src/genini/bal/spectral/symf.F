      subroutine symf(fc,gw)
      implicit none

#include <dims.h>

      REAL_TYPE fc(0:1,nlat2,nnlt,0:1)
      REAL_TYPE gw(nlat2,nnlt)

      REAL_TYPE fc_re,fc_im
      INT_TYPE j,k

      do k=1,nnlt
         do j=1,nlat2
            fc_re=fc(0,j,k,0)
            fc(0,j,k,0)=gw(j,k)*(fc_re+fc(0,j,k,1))
            fc(0,j,k,1)=gw(j,k)*(fc_re-fc(0,j,k,1))

            fc_im=fc(1,j,k,0)
            fc(1,j,k,0)=gw(j,k)*(fc_im+fc(1,j,k,1))
            fc(1,j,k,1)=gw(j,k)*(fc_im-fc(1,j,k,1))
         end do
      end do

      end
