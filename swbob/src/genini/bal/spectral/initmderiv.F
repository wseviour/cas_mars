      subroutine initmderiv(eps,m)
      implicit none

#include <dims.h>

      INT_TYPE  m
      REAL_TYPE eps(m:nn+1)

      INT_TYPE n

      ! ==================================================
      !
      ! See equations 2.52 and 2.53 in Rudy Jacob's thesis
      ! NCAR/CT-144, 1993 for Hmn coefs below.
      ! ==================================================
      
      eps(m)=0.0d0
      do n=m+1,nn+1
         eps(n)=DSQRT(((n+m)*(n-m))/((2.0d0*n+1.0d0)*(2.0d0*n-1.0d0)))
      end do

      end



