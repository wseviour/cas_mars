#include <type.h>
      subroutine inittheforcingic(thfor,
     $                         theq,
     $                         pktcoef,
     $                         k)
      implicit none

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      Include Files
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include <dims.h>
#include <planet.h>
#include <fordamp0.h>

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      REAL_TYPE thfor(nlon,nlat2,0:1)
      REAL_TYPE theq(nlat2,0:1)
      REAL_TYPE pktcoef(nlat2,0:1)

      INT_TYPE i,j,ins,k

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     theq is read in from ic.dat in initforsc

      do ins=0,1
         do j=1,nlat2
            do i=1,nlon
               thfor(i,j,ins) = theq(j,ins)*pktcoef(j,ins)
            enddo
         enddo
      enddo

      end


        



