#include <type.h>
      subroutine inituth(uini,vorini,thini)
      implicit none

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      Include Files
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include <coord.h>
#include <dims.h>
#include <constants.h>
#include <planet.h>

      REAL_TYPE uini(plat),u1(plat),u2(plat)
      REAL_TYPE vorini(plat)
      REAL_TYPE thini(plat)
      REAL_TYPE thinimean

      REAL_TYPE sinlat(plat), rlat(plat), dmu(plat)
      INT_TYPE i,j,ins
      REAL_TYPE snlat,cslat,deglat

      REAL_TYPE umax
      REAL_TYPE dthdphi(plat)
      REAL_TYPE hw,two,dg0,ph0

      do j=1,plat/2
         rlat(j)=Pi/2-colatitude(j)
         rlat(plat+1-j)=-rlat(j)
         sinlat(j)=DSIN(Pi/2-colatitude(j))
         sinlat(plat+1-j)=-sinlat(j)
      end do

c     sech^2 jet at 45N:
      umax=0.1d0
      hw=1.d0/15.d0
      dg0=20.d0
      ph0=dg0*pi/180.d0
      two=2.d0

      do j=1,plat
         u1(j) = umax / ( dcosh((rlat(j)-ph0)/hw) )**2
         u2(j) = umax / ( dcosh((rlat(j)+ph0)/hw) )**2
         uini(j) = u1(j) + u2(j)
      enddo

c     vorticity = -1/acosphi*d(cosphi*uini)/dphi
      do j=1,plat
         snlat=sinlat(j)
         cslat=coslat(j)
         vorini(j) = (1.d0/radius) * ( uini(j)*snlat/cslat
     :        + (two/hw)*dtanh((rlat(j)-ph0)/hw)*u1(j)
     :        + (two/hw)*dtanh((rlat(j)+ph0)/hw)*u2(j) )
     :        + two*Omega*snlat
      enddo

c     (af+u*tanhphi)*u=dthini/dphi
      do j=1,plat
         snlat=sinlat(j)
         cslat=coslat(j)
         dthdphi(j) = - ( 2*omega*snlat*radius + uini(j)*snlat/cslat )*uini(j)
      enddo

c     integrate w.r.t phi
      thini(1)=0.d0
      do j=2,plat
         thini(j)=thini(j-1)+
     :           0.5d0*(dthdphi(j)+dthdphi(j-1))*(rlat(j)-rlat(j-1))
      enddo

c     Subtract mean thini from thini at each level
c     (mu goes from 1 to -1)

      do j=2,plat-2
         dmu(j)=0.5d0*(sinlat(j-1)-sinlat(j+1))
      enddo
      dmu(1)=0.5d0*(1.d0-sinlat(2))
      dmu(plat)=0.5d0*(sinlat(plat-1)-(-1.d0))
      thinimean=0.d0
      do j=1,plat
         thinimean=thinimean+thini(j)*dmu(j)*0.5d0
      enddo
      do j=1,plat
         thini(j)=thini(j)-thinimean
      enddo

      write(*,*) ' -------------- vorini --------------- '
      write(*,*) vorini

      write(*,*) ' --------------  uini  --------------- '
      write(*,*) uini

      write(*,*) ' --------------  thini --------------- '
      write(*,*) thini

      return
      end

