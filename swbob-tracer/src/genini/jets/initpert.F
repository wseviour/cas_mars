#include <type.h>
      subroutine initpert(phi,
     $                   rlat,
     $                   coslat,
     $                   sinlat,
     $                   lon,
     $                   k)
 
      implicit none

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Include Files
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#include <dims.h>
#include <constants.h>
#include <planet.h>

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      REAL_TYPE phi(nlon,nlat2,0:1)
      REAL_TYPE coslat(nlat2),sinlat(nlat2),rlat(nlat2)
      REAL_TYPE lon(nlon)
      REAL_TYPE cslat,snlat

      REAL_TYPE lonc,latc,hw0,hw1,pertamp
      REAL_TYPE londep1(nlon), londep2(nlon)
      REAL_TYPE londep(nlon), latdep(nlat2,0:1)
      REAL_TYPE th2t

      INT_TYPE i,j,k,ins

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     bump is centred in longitude on 0E
      hw0  = 1.d0/3.d0

      latc = 45.d0
      hw1  = 1.d0/15.d0
      pertamp = 00.d0

      do i=1,nlon
         if(lon(i).le.pi)
     :        londep(i)=1.d0/(dcosh( (lon(i)-0.d0)/hw0 ))**2
         if(lon(i).gt.pi)
     :        londep(i)=1.d0/(dcosh( (2.d0*pi-lon(i))/hw0 ))**2
      enddo
      do j=1,nlat2
         latdep(j,0)=1.d0/(dcosh( ( rlat(j)-latc*pi/180.d0)/hw1 ))**2
         latdep(j,1)=1.d0/(dcosh( (-rlat(j)-latc*pi/180.d0)/hw1 ))**2
      enddo
        
      do ins=0,1
         do j=1,nlat2
            do i=1,nlon
               phi(i,j,ins)=phi(i,j,ins)+pertamp*londep(i)*latdep(j,ins)
            end do
         end do
      end do

      end

