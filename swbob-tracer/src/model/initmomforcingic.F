#include <type.h>
      subroutine initmomforcingic(ufor,vfor,
     $                         ueq,
     $                         pkvcoef,
     $                         coslat,
     $                         k)
      implicit none

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      Include Files
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include <dims.h>
#include <planet.h>
#include <fordamp0.h>

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      REAL_TYPE ufor(nlon,nlat2,0:1),vfor(nlon,nlat2,0:1)
      REAL_TYPE ueq(nlat2,0:1)
      REAL_TYPE pkvcoef(nlat2,0:1)
      REAL_TYPE coslat(nlat2)

      INT_TYPE i,j,ins,k

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     ueq is equilibrium velocity, read in from ic.dat in initforsc
c     ufor and vfor are added to the terms B and A in bldnlt.F

      do ins=0,1
         do j=1,nlat2
            do i=1,nlon
               ufor(i,j,ins) = ueq(j,ins)*pkvcoef(j,ins)*coslat(j)
               vfor(i,j,ins) = 0.d0
            enddo
         enddo
      enddo

      end


        



