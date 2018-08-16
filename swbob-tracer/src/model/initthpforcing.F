#include <type.h>
      subroutine initthpforcing(thpfor,
     $                         coslat,
     $                         sinlat,
     $                         rlat,
     $                         lon,
     $                         k)
      implicit none

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      Include Files
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include <dims.h>
#include <planet.h>
#include <fordamp0.h>

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      REAL_TYPE thpfor(nlon,nlat2,0:1)

      REAL_TYPE coslat(nlat2)
      REAL_TYPE sinlat(nlat2)
      REAL_TYPE rlat(nlat2)

      INT_TYPE i,j,ins,k
      REAL_TYPE snlat,cslat
      REAL_TYPE lon(nlon)
      REAL_TYPE londep(nlon), latdep(nlat2,0:1)
      REAL_TYPE degrad,rad0,drad

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     input parameters: lat0, wavenumber, amp

      degrad=Pi/180.d0
      rad0=lath*degrad
      drad=dlath*degrad

      do j=1,nlat2
         latdep(j,0)=exp(-(( rlat(j)-rad0)/drad)**2)
         latdep(j,1)=exp(-((-rlat(j)-rad0)/drad)**2)

c     SP2012 exoplanet:
c         latdep(j,0)=dcos(rlat(j))
c         latdep(j,1)=dcos(rlat(j))

      enddo
      do i=1,nlon
         londep(i)=dcos(waveh*lon(i))
      enddo

      do ins=0,1
         do j=1,nlat2
            do i=1,nlon
               thpfor(i,j,ins)=amph*londep(i)*latdep(j,ins)
            end do
         end do
      end do

      end


        




