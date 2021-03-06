      subroutine initspecuv(Alpha_coord,Radius,Omega)
      implicit none

c------------------------------------------------------------------
c      initialize coriolis factors in spectral space
c      see equation B.18 in NCAR Technical Note TN-343+STR
c      entitled "Description of a Global Shallow Water Model Based
c      on the Spectral Transform Method" by James J. Hack
c      and Ruediger Jakob (Feb 1992).
c-----------------------------------------------------------------

#include <dims.h>
#include <specuv.h>
c      REAL_TYPE fcor(0:1)         ! coriolis factor spectral coefficients m=1,n=0,1
c      REAL_TYPE annp1(0:nn+1)     ! Radius/(n*(n+1))
c      common /specuv/fcor,annp1

      REAL_TYPE Alpha_coord,Radius,Omega
      INT_TYPE n

c-----------------------------------------------------------------
c     Executable statements
c-----------------------------------------------------------------

      fcor(0) =  DSQRT(8.0d0/3.0d0)*Omega*DCOS(Alpha_coord)
      fcor(1) = -DSQRT(4.0d0/3.0d0)*Omega*DSIN(Alpha_coord)

      annp1(0)=0
      do n=1,nn+1
         annp1(n)=Radius/(n*(n+1))
      end do

      end 
