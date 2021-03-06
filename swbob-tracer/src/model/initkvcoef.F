#include <type.h>
      subroutine initkvcoef(kvcoef,
     $     coslat,
     $     sinlat,
     $     rlat,
     $     k)
      implicit none

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      Include Files
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include <dims.h>
#include <planet.h>
#include <fordamp0.h>
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      REAL_TYPE kvcoef(nlat2,0:1)

      REAL_TYPE coslat(nlat2)
      REAL_TYPE sinlat(nlat2)
      REAL_TYPE rlat(nlat2)

      INT_TYPE i,j,ins,k
      REAL_TYPE env,degrad,dradm,rlat0
      REAL_TYPE env1

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      degrad=Pi/180.d0
      rlat0=latm*degrad
      dradm=dlatm*degrad

      do ins=0,1
         do j=1,nlat2
c old:
c            cslat=coslat(j)
c            if (ins.eq.0) then
c               wst=0.5d0*(1.d0+dtanh((rlat(j)-lat0)/dlat))
c            else
c               wst=0.5d0*(1.d0+dtanh((-rlat(j)-lat0)/dlat))
c            end if
c            kvcoef(j,ins) = kv*(1-wst)

c            env=exp(-(rlat(j)/dradm)**2)
            env1=exp(-(rlat(j)/(15.*degrad))**2)

            if (ins.eq.0) then
               env=0.5d0*(1.d0-dtanh((rlat(j)-rlat0)/abs(dradm)))
            else
               env=0.5d0*(1.d0-dtanh((-rlat(j)-rlat0)/abs(dradm)))
            end if
            kvcoef(j,ins) = kv*env
c            kvcoef(j,ins) = kv

c            write(*,'(i2,i2,2x,f5.2,2x,f6.4,2x,f6.4)') ins,j,rlat(j)/degrad,env,env1

         end do
      end do

      end
