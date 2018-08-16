#include <type.h>
      subroutine initpert(vor,
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

      REAL_TYPE vor(nlon,nlat2,0:1)
      REAL_TYPE coslat(nlat2),sinlat(nlat2),rlat(nlat2)
      REAL_TYPE lon(nlon)
      REAL_TYPE cslat,snlat

      REAL_TYPE lonc,latc,hw0,hw1,pertamp
      REAL_TYPE londep1(nlon), londep2(nlon)
      REAL_TYPE londep(nlon), latdep(nlat2,0:1)
      REAL_TYPE th2t

      REAL_TYPE r0

      INT_TYPE i,j,k,ins

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     bump is centred in longitude on 0E
      hw0  = 1.d0/3.d0

      latc = 45.d0
      hw1  = 1.d0/15.d0
c      pertamp = 800.d0
      pertamp = 1.d0/1000000.d0
      
c     do i=1,nlon
c         if(lon(i).le.pi)
c     :        londep(i)=1.d0/(dcosh( (lon(i)-0.d0)/hw0 ))**2
c         if(lon(i).gt.pi)
c     :        londep(i)=1.d0/(dcosh( (2.d0*pi-lon(i))/hw0 ))**2
c      enddo
c      do j=1,nlat2
c         latdep(j,0)=1.d0/(dcosh( ( rlat(j)-latc*pi/180.d0)/hw1 ))**2
c         latdep(j,1)=1.d0/(dcosh( (-rlat(j)-latc*pi/180.d0)/hw1 ))**2
c      enddo
c        
c      do ins=0,1
c         do j=1,nlat2
c            do i=1,nlon
c               vor(i,j,ins)=vor(i,j,ins)+pertamp*londep(i)*latdep(j,ins)
c            end do
c         end do
c      end do

      do ins=0,1
         do j=1,nlat2
            do i=1,nlon
               vor(i,j,ins)=vor(i,j,ins)+pertamp*(rand()-0.5d0)
            end do
         end do
      end do


      end

