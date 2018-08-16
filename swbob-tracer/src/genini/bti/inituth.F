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

      REAL_TYPE uini(plat)
      REAL_TYPE vorini(plat)
      REAL_TYPE thini(plat)
      REAL_TYPE thinimean

      REAL_TYPE sinlat(plat), rlat(plat), dmu(plat)
      INT_TYPE i,j,ins
      REAL_TYPE snlat,cslat,deglat

      REAL_TYPE umax
      REAL_TYPE dthdphi(plat)
      REAL_TYPE hw,ffd,two,hwc

      do j=1,plat/2
         rlat(j)=Pi/2-colatitude(j)
         rlat(plat+1-j)=-rlat(j)
         sinlat(j)=DSIN(Pi/2-colatitude(j))
         sinlat(plat+1-j)=-sinlat(j)
      end do

c     uini is a barotropically unstable state (sech^2 jet at 45N)

      umax=80.d0
      hw=1.d0/15.d0
c      hw=1.d0/18.d0
      hwc=pi/18.d0
      ffd=pi/4.d0     !forty five degrees
      two=2.d0

#if 0
      do j=1,plat/2
c         uini(j)=umax * dcos(two*(rlat(j)-ffd)) / ( dcosh((rlat(j)-ffd)/hw) )**2
         uini(j)=umax / ( dcosh((rlat(j)-ffd)/hw) )**2
      enddo
      do j=plat/2+1,plat
         uini(j)=0.d0
      enddo

c     vorini is the corresponding vorticity = -1/acosphi*d(cosphi*uini)/dphi

      do j=1,plat
         snlat=sinlat(j)
         cslat=coslat(j)
c         vorini(j) = uini(j)/radius * ( snlat/cslat
c     :        + two*dsin(two*(rlat(j)-ffd))/dcos(two*(rlat(j)-ffd))
c     :        + two/hw*dtanh((rlat(j)-ffd)/hw) )
c     :        + 2.0d0*Omega*snlat
         vorini(j) = uini(j)/radius * ( snlat/cslat
     :        + two/hw*dtanh((rlat(j)-ffd)/hw) )
     :        + 2.0d0*Omega*snlat
      enddo
#endif
#if 1
      do j=1,plat
         if(dabs(rlat(j)-ffd).le.hwc) then
            uini(j)=umax * ( dcos(pi/2.d0*(rlat(j)-ffd)/hwc) )**2
         else
            uini(j)=0.d0
         endif
      enddo
      do j=plat/2+1,plat
         uini(j)=0.d0
      enddo

c     vorini is the corresponding vorticity = -1/acosphi*d(cosphi*uini)/dphi

      do j=1,plat
         snlat=sinlat(j)
         cslat=coslat(j)
         if(dabs(rlat(j)-ffd).le.hwc) then
            vorini(j) = 1.d0/radius * ( 
     :           snlat/cslat*uini(j)
     :           + umax*pi/hwc*dcos(pi/2.d0*(rlat(j)-ffd)/hwc)
     :                        *dsin(pi/2.d0*(rlat(j)-ffd)/hwc) )
     :           + 2.0d0*Omega*snlat
         else
            vorini(j) = 2.0d0*Omega*snlat
         endif
      enddo
#endif




c     thini is the height in gradient wind balance with uini
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

