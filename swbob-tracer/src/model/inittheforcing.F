#include <type.h>
      subroutine inittheforcing(thfor,
     $                         theq,
     $                         pktcoef,
     $                         coslat,
     $                         sinlat,
     $                         rlat,   !da,thsum,
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

      REAL_TYPE coslat(nlat2)
      REAL_TYPE sinlat(nlat2)
      REAL_TYPE rlat(nlat2)
!      REAL_TYPE da(nlat2)

      INT_TYPE i,j,ins,k
      REAL_TYPE snlat,cslat

      REAL_TYPE wst,degrad,radj,drad
      REAL_TYPE snlt0,ueq0

      REAL_TYPE p0,p1,p2,p3,p4,mu  !,thsum

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     computes thfor=kt*theta(eq) not t(eq)

      degrad=Pi/180.d0
      radj=latj*degrad
      drad=dlat*degrad

c     polvor = strength of polar vortex, or Hadley cell, or ...
c     latj = latitude of jet, or Hadley cell extent, or ...
c     (both in fordamp0.h)


#if 1
c     polar vortex:
      do ins=0,1
         do j=1,nlat2
            cslat=coslat(j)
            if (ins.eq.0) then
               wst=dtanh((rlat(j)-radj)/drad)
            else
               wst=dtanh((-rlat(j)-radj)/drad)
            end if

            theq(j,ins) = -0.5d0*(1.d0+wst)*polvor

            do i=1,nlon
               thfor(i,j,ins) = theq(j,ins)*pktcoef(j,ins)
            enddo
         enddo
      enddo


      if (dlat.eq.99.0) then 
c     polar vortex, similar to Rong+Waugh profile, using Legendre polynomials
c         write(*,*) "Using Rong+Waugh like profile"
      do ins=0,1
         do j=1,nlat2
            mu=(-1.0)**ins*sinlat(j)
            p1=mu
            p2=0.5*(3.0*mu**2-1.0)
            p3=0.5*(5.0*mu**3-3.0*mu)
            p4=0.125*(35.0*mu**4-30.0*mu**2+3.0)
c            theq(j,ins) = -0.125*(3.0*p1 + 2.0*p2 + 2.0*p3 + p4)*polvor
c            if (ins.eq.0) theq(j,ins) = -(0.5*p2 + (1./3.)*p4 + 0.125)*polvor
c            if (ins.eq.1) theq(j,ins) = (0.25*p2 + 0.125)*polvor
            if (ins.eq.0) then 
c               theq(j,ins) = -polvor*(sin((pi/2.)*(rlat(j)-radj)/(pi/2.-radj)))**2
               theq(j,ins) = -polvor*(sin(rlat(j)))**3
               if (rlat(j).lt.radj) theq(j,ins) = 0.d0
            endif
            if (ins.eq.1) theq(j,ins) = 0.d0

!            thsum=thsum+theq(j,ins)*da(j)

c     add contribution corresponding to solid body rotation, usb*cosphi:
c            theq(j,ins)=theq(j,ins)+usb*(2.d0*omega*radius+usb)*(coslat(j))**2

            do i=1,nlon
               thfor(i,j,ins) = theq(j,ins)*pktcoef(j,ins)
            enddo
         enddo
      enddo
      endif
#endif

#if 0
c     Hadley cell, following Shell & Held, 2004 (symmetric about equator):
      ueq0=polvor
      snlt0=sin(radj)

      do ins=0,1
         do j=1,nlat2
            snlat=sinlat(j)
            if (abs(rlat(j)).le.radj) then
               theq(j,ins) =  -ueq0*radius*omega*snlat**2
            else
               theq(j,ins) =  -ueq0*radius*omega*snlt0**2
            endif
c     modified SH04 - sine profile extending all the way to the pole:
c     (set phibar=0.5*ueq0max*radius*omega, and allow 0<ueq0<ueq0max)
c            theq(j,ins)= ueq0*radius*omega*(0.5-snlat**2)

c     changed (nov.2012): now allow dphi=ueq0*radius*omega<=3/2*phibar
            theq(j,ins)= ueq0*radius*omega*(1.d0/3.d0-snlat**2)

c            
            do i=1,nlon
               thfor(i,j,ins) = theq(j,ins)*pktcoef(j,ins)
            enddo
         enddo
      enddo
#endif

      end


        



