#include <type.h>
      subroutine inituthq(uini,vorini,thini)
      implicit none

c     find u,zeta,theta in balance with a specified zonal PV field
c     by iteration

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      Include Files
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include <coord.h>
#include <dims.h>
#include <constants.h>
#include <planet.h>

      REAL_TYPE uini(plat),vorini(plat),thini(plat)
      REAL_TYPE q0(plat),z0(plat),u0(plat)
      REAL_TYPE h0(plat),h1(plat),z1(plat),q1(plat)
      REAL_TYPE dhdphi(plat),dhdmu(plat)

      REAL_TYPE sinlat(plat), rlat(plat), dmu(plat), da(plat)
      REAL_TYPE snlat,cslat,deglat,mu

      REAL_TYPE two,dg0,ph0,m0,dm,hhp,hhm,datot,uc
      REAL_TYPE grav,hbar,htmp,htmp1,fac,herr,tol,qtmp,qtmp1,qtmp2,ztmp

      INT_TYPE i,j,ins,jc

      two=2.d0

      open(20,file='q0',form='formatted')
      open(21,file='z0',form='formatted')
      open(22,file='u0',form='formatted')
      open(23,file='h0',form='formatted')


c     hbar and gravity:
      hbar=1.d0
      grav=phibar/hbar

      write(*,*) 'radius, gravity = ',radius, grav

      do j=1,plat/2
         rlat(j)=Pi/2-colatitude(j)
         rlat(plat+1-j)=-rlat(j)
         sinlat(j)=DSIN(Pi/2-colatitude(j))
         sinlat(plat+1-j)=-sinlat(j)
         da(j)=gausswt(j)
         da(plat+1-j)=gausswt(j)
      end do

      dg0=45.d0
      ph0=dg0*pi/180.d0
      m0=dsin(ph0)
      dm=0.1

      m0=0.3
      hhm=0.12
      hhp=0.1

      do j=1,plat
         jc=plat+1-j
         mu=sinlat(j)
c         q0(j)=two*omega*(mu+(mu-m0)*exp(-((mu-m0)/dm)**2))
         q0(j)=two*omega*(mu+(mu-m0)*exp(-((mu-m0)/dm)**2)
     .                      +(mu+m0)*exp(-((mu+m0)/dm)**2))
c         q0(jc)=-q0(j)
      enddo
c      do j=1,plat/2
c         jc=plat+1-j
c         mu=sinlat(j)
c         if (mu.lt.(m0-hhm)) then 
c            q0(j)=two*omega*mu
c         else if (mu.lt.(m0)) then
c            q0(j)=two*omega*(m0-hhm)
c         else if (mu.lt.(m0+hhp)) then
c            q0(j)=two*omega*(m0+hhp)
c         else
c            q0(j)=two*omega*mu
c         endif
c         q0(jc)=-q0(j)
c      enddo

c     check global PV=0
      qtmp=0.d0
      do j=1,plat
         qtmp=qtmp+q0(j)*da(j)
      enddo
      write(*,*) 'global PV = ',qtmp*0.5d0

      tol=1e-15
      do j=1,plat
         h0(j)=0.d0
      enddo

 1    continue

      do j=1,plat
         z0(j)=q0(j)*(hbar+h0(j))-two*omega*sinlat(j)
      enddo

c     set global zeta to 0
      ztmp=0.d0
      do j=2,plat
         ztmp=ztmp+z0(j)*da(j)
      enddo
      ztmp=ztmp*0.5d0
      write(*,*) 'global zeta = ',ztmp
      do j=1,plat
         z0(j)=z0(j)-ztmp
      enddo

c     find u0
      u0(1)=0
      do j=2,plat
         u0(j)=u0(j-1)+z0(j)*da(j)
      enddo
      do j=1,plat
         u0(j)=u0(j)*radius/coslat(j)
      enddo

c     (af+u*tanhphi)*u=dthini/dphi=g*dhdphi
      do j=1,plat
         uc=u0(j)/coslat(j)
c         dhdphi(j) = - ( 2*omega*snlat*radius + u0(j)*snlat/cslat )*u0(j)/grav
         dhdmu(j) = ( 2*omega*radius + uc )*uc*sinlat(j)/grav
      enddo
c     integrate w.r.t phi
      h1(1)=0.d0
      do j=2,plat
         h1(j)=h1(j-1)+dhdmu(j)*da(j)
      enddo

c     make h even about mu=0:
      do j=1,plat/2
         jc=plat+1-j
         htmp=0.5d0*(h1(j)+h1(jc))
         h1(j)=htmp
         h1(jc)=htmp
      enddo

c     Subtract mean h from h at each level
      htmp=0.d0
      do j=1,plat
         htmp=htmp+h1(j)*da(j)
      enddo
      htmp=htmp*0.5d0
      do j=1,plat
         h1(j)=h1(j)-htmp
      enddo

c     mean square diff between h0 and h1:
      herr=0.d0
      do j=1,plat
         herr=herr+(h0(j)-h1(j))**2*da(j)*0.5d0
      enddo
      write(*,*) '  ****  herr = ', herr, '  ****  '

c      do j=1,plat
c         q1(j)=(z0(j)+two*omega*sinlat(j))/(hbar+h1(j))
c      enddo
c      write(20,*) q1

      write(21,*) z0
      write(22,*) u0
      write(23,*) h0

      do j=1,plat
         h0(j)=h1(j)
      enddo

      if (herr.ge.tol) goto 1

c      do j=1,plat
c         z1(j)=q0(j)*(hbar+h1(j))-two*omega*sinlat(j)
c      enddo

      write(*,*) ' ---------------- q0 ----------------- '
c      write(*,*) q0
      write(*,*) ' ---------------- z0 ----------------- '
c      write(*,*) z0
      write(*,*) ' ---------------- u0 ----------------- '
c      write(*,*) u0
      write(*,*) ' ---------------- h1 ----------------- '
c      write(*,*) h1
      write(*,*) ' ---------------- z1 ----------------- '
c      write(*,*) z1

      do j=1,plat
         vorini(j) = z0(j)+two*omega*sinlat(j)
         uini(j)   = u0(j)
         thini(j)  = grav*h0(j)
      enddo

c      write(*,*) ' ---------------- q0 ----------------- '
c      write(*,*) q0

c      write(*,*) ' -------------- vorini --------------- '
c      write(*,*) vorini

c      write(*,*) ' --------------  uini  --------------- '
c      write(*,*) uini

c      write(*,*) ' --------------  thini --------------- '
c      write(*,*) thini

      return
      end

