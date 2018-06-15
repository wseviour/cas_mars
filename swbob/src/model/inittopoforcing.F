#include <type.h>
      subroutine inittopoforcing(phitopo,phitopokt,
     $                         pktcoef,
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
      REAL_TYPE phitopo(nlon,nlat2,0:1)
      REAL_TYPE phitopokt(nlon,nlat2,0:1)
      REAL_TYPE pktcoef(nlat2,0:1)

      REAL_TYPE coslat(nlat2)
      REAL_TYPE sinlat(nlat2)
      REAL_TYPE rlat(nlat2)

      INT_TYPE i,j,ins,k
      REAL_TYPE snlat,cslat
      REAL_TYPE lon(nlon)
      REAL_TYPE londep(nlon), latdep(nlat2,0:1)

      REAL_TYPE degrad,rad0,drad

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     input parameters: lat0, wavenumber, amp, (a1,a2,a3)

      degrad=Pi/180.d0
      rad0=latt*degrad
      drad=dlatt*degrad

c      a1=1.d0/sqrt(3.d0)
c      a2=1.d0/sqrt(3.d0)
c      a3=1.d0/sqrt(3.d0)
c      a1=1.d0/sqrt(7.d0)

#if 1
      do j=1,nlat2
         if(rlat(j).ge.rad0) then
            latdep(j,0)=( dsin((rlat(j)-rad0)/(pi/2.d0-rad0)*pi) )**2
         else
            latdep(j,0)=0.d0
         endif
         latdep(j,1)=0.d0
      enddo
#endif

      if (latt.lt.0.d0) then
         write(*,*) "global topography"
      do j=1,nlat2
         latdep(j,0)=dcos(rlat(j))
         latdep(j,1)=dcos(rlat(j))
      enddo
      endif

      if (wavet.lt.0.0) then
         write(*,*) "localized topography"
      do j=1,nlat2
         latdep(j,0)=exp(-(( rlat(j)-rad0)*wavet*2.d0/pi)**2)
         latdep(j,1)=exp(-((-rlat(j)-rad0)*wavet*2.d0/pi)**2)
      enddo
      endif

#if 0
      do j=1,nlat2
         latdep(j,0)=exp(-(( rlat(j)-rad0)/drad)**2)
c         latdep(j,1)=exp(-((-rlat(j)-rad0)/drad)**2)
         latdep(j,1)=latdep(j,0) ! symmetric N/S

c     SP2012 exoplanet:
c         latdep(j,0)=dcos(rlat(j))
c         latdep(j,1)=dcos(rlat(j))

      enddo
#endif

      do i=1,nlon
         if (wavet.gt.0.) londep(i)=dcos(wavet*lon(i))
         if (wavet.lt.0.) londep(i)=exp(-((lon(i)-pi)*wavet*2.d0/pi)**2)

c         londep(i)=a1*dcos(lon(i))+a2*dcos(2.d0*lon(i))+a3*dcos(3.d0*lon(i))

c         londep(i)=a1*(dcos(lon(i))+dcos(2.d0*lon(i))+dcos(3.d0*lon(i))
c     .        +dcos(4.d0*lon(i))+dcos(5.d0*lon(i))+dcos(6.d0*lon(i))+dcos(7.d0*lon(i)))
c         londep(i)=exp(-((lon(i)-pi)*wavet*2.d0/pi)**2)
      enddo

      do ins=0,1
         do j=1,nlat2
            do i=1,nlon
               phitopo(i,j,ins)=ampt*londep(i)*latdep(j,ins)
               phitopokt(i,j,ins)=phitopo(i,j,ins)*pktcoef(j,ins)
             end do
         end do
      end do

      end
