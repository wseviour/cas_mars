#include <type.h>
      subroutine inituth(vorini,thini,sigini)
      implicit none

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      Include Files
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include <coord.h>
#include <dims.h>
#include <constants.h>
#include <planet.h>

c     rks - added initial tracer distribution (sigini)
      
      REAL_TYPE vorini(plat)
      REAL_TYPE thini(plat)
      REAL_TYPE sigini(plat)

      REAL_TYPE sinlat(plat), rlat(plat), dmu(plat)
      INT_TYPE i,j,ins
      REAL_TYPE snlat,cslat,deglat

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

c     rks - tracer=1 between 30S and 30N, tracer=0 elsewhere:
      do j=1,plat
         snlat=sinlat(j)
         sigini(j) = 0.d0
         if (abs(snlat).le.0.5d0) sigini(j) = 1.d0
      enddo
      
      do j=1,plat
         thini(j)=0.0d0
      enddo

      write(*,*) ' -------------- vorini --------------- '
      write(*,*) vorini

      write(*,*) ' --------------  thini --------------- '
      write(*,*) thini

      write(*,*) ' --------------  sigini --------------- '
      write(*,*) sigini

      return
      end

