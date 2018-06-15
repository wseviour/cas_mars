      subroutine initspecuv(Alpha_coord,Radius,Omega)
      implicit none

#include <dims.h>
#include <specuv.h>

      REAL_TYPE Alpha_coord,Radius,Omega
      INT_TYPE n

      ! ----------------------------------------------------
      !
      ! initialize coreolis factors in spectral space
      ! see equation B.18 in NCAR Technical Note TN-343+STR
      ! entitled "Description of a Global Shallow Water Model Based
      ! on the Spectral Transform Method" by James J. Hack
      ! and Ruediger Jakob (Feb 1992).
      !
      ! ----------------------------------------------------

      fcor(0) =  DSQRT(8.0d0/3.0d0)*Omega*DCOS(Alpha_coord)
      fcor(1) = -DSQRT(4.0d0/3.0d0)*Omega*DSIN(Alpha_coord)

      annp1(0)=0
      do n=1,nn+1
         annp1(n)=Radius/(n*(n+1))
      end do

      end 
