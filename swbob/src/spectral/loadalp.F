       subroutine loadalp(Pseed,P,m)
       implicit none

#include <dims.h>

       INT_TYPE  m
       REAL_TYPE Pseed(nlat2,0:nn+1)
       REAL_TYPE P(nlat2,0:nn+1)

       INT_TYPE n,j

       do n=m,nn+1
          do j=1,nlat2
             P(j,n)=Pseed(j,n)
          end do
       end do

       end
  


