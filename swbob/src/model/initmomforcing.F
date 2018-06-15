#include <type.h>
c      subroutine initvorforcing(vorfor,
      subroutine initmomforcing(foru,forv,
     $                         coslat,
     $                         sinlat,
     $                         rlat,
     $                         rlon,
     $                         k)
      implicit none

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      Include Files
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include <dims.h>
#include <planet.h>
#include <fordamp0.h>

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      REAL_TYPE vorfor(nlon,nlat2,0:1)
      REAL_TYPE foru(nlon,nlat2,0:1)
      REAL_TYPE forv(nlon,nlat2,0:1)

      REAL_TYPE coslat(nlat2)
      REAL_TYPE sinlat(nlat2)
      REAL_TYPE rlat(nlat2)
      REAL_TYPE rlon(nlon)

      INT_TYPE i,j,ins,k
      REAL_TYPE snlat,cslat

      REAL_TYPE mu0,fhat,fhatp,fhatz,arg,degrad
      REAL_TYPE fu,fv

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     lat0,wavenumber,amp in model/fordamp0.h
c     a1 is the width of the forcing profile (in mu)

c     equatorial momentum forcing a al Shell & Held, 2004:
c     (parameters: amp,wavenumber)
#if 0
      do ins=0,1
         do j=1,nlat2
            cslat=coslat(j)
            fu=ampm*cslat**wavem
            fv=ampm*cslat**wavem
            do i=1,nlon
               foru(i,j,ins) = fu
               forv(i,j,ins) = fv    ????
            enddo
         enddo            
      enddo
#endif

c     solid body rotation:
#if 0
      do ins=0,1
         do j=1,nlat2
            cslat=coslat(j)
            fu=ampm*cslat
            do i=1,nlon
               foru(i,j,ins) = fu
               forv(i,j,ins) = 0.0
            enddo
         enddo            
      enddo
#endif

c forcing a midlatitude jet (compact divF):
c (need to update to new implementation of FMECH in bldnlt (just need F, not curl(F))
#if 0
      degrad=Pi/180.d0
      mu0=dsin(lat0*degrad)

      do ins=0,1
         do j=1,nlat2
            cslat=coslat(j)
            if (ins.eq.0) then
               snlat=sinlat(j)
            else
               snlat=-sinlat(j)
            endif

            if (abs(snlat-mu0).le.a1) then
               arg=pi*(snlat-mu0)/(2.d0*a1)
               fhat=amp*(dcos(arg))**4
               fhatp=-amp*4.d0*(pi/(2.d0*a1))*dsin(arg)*(dcos(arg))**3
               fhatz=(-fhatp*cslat+fhat*snlat/cslat)/radius
            else
               fhatz=0.d0
            endif

            do i=1,nlon
               vorfor(i,j,ins) = fhatz
            enddo
         enddo
      enddo
#endif

      end


        



