      subroutine ell_bld(phy,lon,lat,da,mm2,mm4,area,qc)

c     CALLED BY: fld.F

      implicit none

#include <dims.h>            ! nlat2
#include <stptrs.h>          ! u_p
#include <diagptrs.h>
#include <planet.h>          ! now includes phibar (7/30/05)

      REAL_TYPE phy(nlon,nlat2,nphy,0:1,plev)
      REAL_TYPE lon(nlon),lat(nlat2),da(nlat2)
      REAL_TYPE za,h,q,qc,x,y,hpi,colat
      REAL_TYPE area,mm2(3),mm4(3)

      INT_TYPE i,j,k,ins

      hpi=pi/2.d0
      do k=1,plev
         do ins=0,1
            do j=1,nlat2
               colat=hpi-lat(j)*(-1.d0)**ins
               do i=1,nlon

                  za   =phy(i,j,vor_p,ins,k)
                  h    =phy(i,j,theta_p,ins,k)/PhiBar
                  q=za/(1.d0+h)

                  if (q.ge.qc) then
c                  if (colat.le.(pi/6.)) then

c     Lambert projection:
                     x=2.d0*dsin(colat/2.d0)*dcos(lon(i))
                     y=2.d0*dsin(colat/2.d0)*dsin(lon(i))

                     area=area+da(j)
                     mm2(1)=mm2(1)+x**2*da(j)
                     mm2(2)=mm2(2)+x*y*da(j)
                     mm2(3)=mm2(3)+y**2*da(j)
                     mm4(1)=mm4(1)+x**4*da(j)
                     mm4(2)=mm4(2)+x**2*y**2*da(j)
                     mm4(3)=mm4(3)+y**4*da(j)

                  endif

               end do
            end do
         end do
      end do

      end 

















