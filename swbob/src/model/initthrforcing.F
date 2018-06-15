#include <type.h>
      subroutine initthrforcing(fthr0,
     $                          fthe0,
     $                          fthp0,
     $                          k)
      implicit none

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      Include Files
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include <dims.h>
#include <planet.h>
#include <fordamp0.h>

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      REAL_TYPE fthr0(nlon,nlat2,0:1)
      REAL_TYPE fthe0(nlon,nlat2,0:1)
      REAL_TYPE fthp0(nlon,nlat2,0:1)

      INT_TYPE i,j,ins,k

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     computes fthr=fthe+fthp (needed to build r term in bldnlt)

      do ins=0,1
         do j=1,nlat2
            do i=1,nlon
               fthr0(i,j,ins)=0.d0
#ifdef FRLAX
               fthr0(i,j,ins)=fthr0(i,j,ins)+fthe0(i,j,ins)
#endif
#ifdef FHEAT
               fthr0(i,j,ins)=fthr0(i,j,ins)+fthp0(i,j,ins)
#endif
c               fthr0(i,j,ins)=fthe0(i,j,ins)+fthp0(i,j,ins)
            enddo
         enddo
      enddo

      end
