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

      INT_TYPE  navu
      external  navu

      REAL_TYPE uini(plat), vorini(plat), thini(plat)

      REAL_TYPE sinlat(plat), rlat(plat), da(plat)

      REAL_TYPE qq(plat), qqn(plat), zz(plat), zz0(plat), ff(plat)
      REAL_TYPE uu(plat), hh(plat), dhdmu(plat)
      REAL_TYPE qq0(plat)


      REAL_TYPE lat0,lat1,dq0,dq1,q0,q1,dtor,rlat0,rlat1

      REAL_TYPE ctop,cbot,coff,hbar,havg
      REAL_TYPE phibari,twom,twoma
      REAL_TYPE alpha,tol,error

      INT_TYPE i,j,iter,iunit
      INT_TYPE hw,iw

      REAL_TYPE wd,wde,fac,sl0,sl1,sli,slj,dy

      INT_TYPE ny

      parameter(ny=100000)
      parameter(tol=1.0d-10,dtor=pi/180.d0)

      REAL_TYPE mu(0:ny),fh(0:ny),qh(0:ny)

      do j=1,plat/2
         rlat(j)=Pi/2-colatitude(j)
         rlat(plat+1-j)=-rlat(j)
         sinlat(j)=DSIN(Pi/2-colatitude(j))
         sinlat(plat+1-j)=-sinlat(j)
         da(j)=gausswt(j)
         da(plat+1-j)=da(j)
      end do

c     test da:
c      ctop=0.d0
c      do j=1,plat
c         ctop=ctop+da(j)
c      enddo
c      write(*,*) 'integral over sphere = ',ctop
c      stop

      phibari=1.d0/phibar

c      phibari=0.d0

      twom=2.d0*omega
      twoma=twom*radius
      hbar=1.d0

c     read inputs.  q should be scaled by hbar/fpole
      iunit=navu()
      open(iunit,file='genx.in',form='formatted',status='old')
      read(iunit,*) lat0,lat1,dq0,dq1
      read(iunit,*) alpha,wd,wde
      close(12)

      rlat0=lat0*dtor
      rlat1=lat1*dtor
      q0=(dsin(rlat0)+dq0)*twom
      q1=q0+dq1*twom

c      write(*,*) 'q0,q1',q0,q1
c      write(*,*) 'rlat0,rlat1',rlat0,rlat1

      do j=1,plat
         ff(j)=twom*sinlat(j)
         hh(j)=hbar
         if (rlat(j).lt.rlat0) then
            qq(j) = ff(j)
         else if (rlat(j).lt.rlat1) then
            qq(j)=q0
         else
            qq(j)=q1
         endif
      enddo

#if 0
c     apply a weak smoothing to eliminate Gibbs fringes in subsequent evolution
      hw=0
      if(hw.ge.1) then
      iw=2*hw+1
      do j=1,plat
         qq0(j)=qq(j)
      enddo
      do j=1+hw,plat-hw
         do i=1,hw
            qq(j)=qq(j)+qq0(j-i)+qq0(j+i)
         enddo
         qq(j)=qq(j)/dble(iw)
      enddo
      endif
#endif


#if 1
c     alternative smoothing (convolve with gaussian kernel):
c      wd=0.01   ! value used for rbsplit
c      wde=10.0d0*wd
c      wd=0.04   ! value used for rbvac
c      wde=3.0d0*wd

c     first make high resolution qq profile:
      sl0=dsin(rlat0)
      sl1=dsin(rlat1)
      do i=0,ny
         mu(i)=1.d0-dble(2*i)/ny
         fh(i)=twom*mu(i)
         if (mu(i).lt.sl0) then
            qh(i) = fh(i)
         else if (mu(i).lt.sl1) then
            qh(i)=q0
         else
            qh(i)=q1
         endif
      enddo
      dy=2.d0/dble(ny)
      fac=1.d0/(sqrt(2.d0*pi)*wd)
      do j=1,plat
         qq0(j)=qq(j)
         slj=dsin(rlat(j))
         if ( (abs(slj-sl0).le.wde) .or. (abs(slj-sl1).le.wde) ) then 
            qq0(j)=0.d0
            do i=0,ny
               sli=mu(i)
               qq0(j)=qq0(j)+qh(i)*exp(-(slj-sli)**2/(2.d0*wd**2))*dy
            enddo
            write(*,'(i3,3(1x,1p,e16.8))') j,slj,qq(j),fac*qq0(j)
            qq(j)=fac*qq0(j)
         endif
      enddo
#endif

c     iterate:
      iter=0
 1    continue

c     ensure global integral of zeta = 0
      ctop=0.d0
      cbot=0.d0
      do j=1,plat
         ctop=ctop+qq(j)*hh(j)*da(j)
         cbot=cbot+hh(j)*da(j)
      enddo

      write(*,*) " ctop, cbot = ", ctop, cbot

      coff=-ctop/cbot

      do j=1,plat
         if (iter.ge.1) zz0(j)=zz(j)
         zz(j)=(coff+qq(j))*hh(j) - ff(j)
         if (iter.ge.1) zz(j)=alpha*zz(j)+(1.d0-alpha)*zz0(j)
      enddo

c      uu(1)=0.5d0*(zz(2)+zz(1))*(rlat(1)-rlat(0))
      uu(1)=zz(1)*da(1)
      do j=2,plat
c         uu(j)=uu(j-1)+0.5d0*(zz(j)+zz(j-1))*(rlat(j)-rlat(j-1))
         uu(j)=uu(j-1)+zz(j)*da(j)
      enddo
c     +ve sign for uu, dhdmu since j goes from NP to SP and da>0
      do j=1,plat
         uu(j)=uu(j)*radius/coslat(j)**2
         dhdmu(j)=phibari*(uu(j)+twoma)*uu(j)*sinlat(j)
      enddo
      do j=1,plat
         hh(j)=0.d0
      enddo
      hh(1)=dhdmu(1)*da(1)
      do j=2,plat
         hh(j)=hh(j-1)+dhdmu(j)*da(j)
      enddo
      havg=0.d0
      do j=1,plat
         havg=havg+hh(j)*da(j)
      enddo
      havg=havg/2.d0

      do j=1,plat
         hh(j)=hh(j)-havg+hbar
      enddo

      iter=iter+1

      error=0.d0
      do j=1,plat
         qqn(j)=(ff(j)+zz(j))/hh(j) - coff
         error=error+sqrt((qqn(j)-qq(j))**2)
      enddo
      write(*,*) iter, error

      if (error.lt.tol) goto 2
      goto 1

 2    continue

c      do j=1,plat
c         write(*,*) qq(j),qqn(j)
c         write(*,*) uu(j)*coslat(j),hh(j)
c      enddo
c      stop

c     now have the total height hh and relative vorticity zz,
c     get height perturbation and absolute vorticity:

      do j=1,plat
         thini(j)=(hh(j)-hbar)*phibar
         vorini(j)=zz(j)+ff(j)
         uini(j)=uu(j)*coslat(j)
      enddo


c      write(*,*) ' -------------- vorini --------------- '
c      write(*,*) ' --------------  uini  --------------- '
c      write(*,*) ' --------------  thini --------------- '

      open(12,file='ic.dat',status='new',form='formatted')
      write(12,*) omega,radius,phibar,coff
      write(12,*) rlat
      write(12,*) qq
      write(12,*) vorini
      write(12,*) thini
      write(12,*) uini
      close(12)

c     tmp: resting atmosphere  (to test spin up from ic.dat thini/uini profile)
c      do j=1,plat
c         thini(j)=0.d0
c         vorini(j)=ff(j)
c         uini(j)=0.d0
c      enddo

      return
      end

