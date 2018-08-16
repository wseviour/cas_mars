! genalp assumes you have (or have loaded) the previous row of 
! ALP's in slot P(1-ix), it uses these to compute P(ix).

      subroutine genalp(P,
     $                  genp,
     $                  ix,
     $                  m)
      implicit none

#include <dims.h>

      INT_TYPE  m
      REAL_TYPE P(nlat2,0:nn+1,0:1)
      REAL_TYPE genp(3,m:nn+1)
      INT_TYPE  ix

      INT_TYPE n,j
      REAL_TYPE c1,c2,c3

      do n=m,nn+1
         c1=genp(1,n)
         c2=genp(2,n)
         c3=genp(3,n)

         do j=1,nlat2
            p(j,n,ix)=p(j,n-2,1-ix)*c1 +
     $                p(j,n  ,1-ix)*c2 +
     $                p(j,n-2,ix  )*c3
         end do
      end do

      end





