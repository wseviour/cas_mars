#include <type.h>
      subroutine initvorth(vorini,thini)
      implicit none

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      Include Files
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include <coord.h>
#include <dims.h>
#include <constants.h>
#include <planet.h>

      REAL_TYPE vorini(plat)
      REAL_TYPE thini(plat)
      REAL_TYPE uini(plat)

      REAL_TYPE sinlat(plat), rlat(plat), dmu(plat)
      INT_TYPE i,j,ins
      REAL_TYPE snlat,cslat

      do j=1,plat/2
         rlat(j)=Pi/2-colatitude(j)
         rlat(plat+1-j)=-rlat(j)
         sinlat(j)=DSIN(Pi/2-colatitude(j))
         sinlat(plat+1-j)=-sinlat(j)
      end do

      do j=1,plat
         snlat=sinlat(j)
         cslat=coslat(j)
         vorini(j) = 2.0d0*Omega*snlat
      enddo

      do j=1,plat
         thini(j)=0.0d0
         uini(j)=0.0d0
      enddo

c     uini need in model/initmomforic.F:
      open(12,file='ic.dat',status='new',form='formatted')
      write(12,*) omega,radius,phibar
      write(12,*) rlat
      write(12,*) vorini
      write(12,*) vorini
      write(12,*) thini
      write(12,*) uini
      close(12)


      return
      end

