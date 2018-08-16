      subroutine gen_symf(fc,gw)
      implicit none

#include <dims.h>

      REAL_TYPE fc(0:1,nlat2,0:1)
      REAL_TYPE gw(nlat2)

      REAL_TYPE fc_re,fc_im
      INT_TYPE j

      do j=1,nlat2
         fc_re=fc(0,j,0)
         fc(0,j,0)=gw(j)*(fc_re+fc(0,j,1))
         fc(0,j,1)=gw(j)*(fc_re-fc(0,j,1))

         fc_im=fc(1,j,0)
         fc(1,j,0)=gw(j)*(fc_im+fc(1,j,1))
         fc(1,j,1)=gw(j)*(fc_im-fc(1,j,1))
      end do

      end
