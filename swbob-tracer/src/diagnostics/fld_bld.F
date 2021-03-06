      subroutine fld_bld(phy,
     $     phy_u,phy_v,phy_z,phy_d,phy_h,phy_q,phy_s,cslat,snlat)

c     CALLED BY: fld.F

      implicit none

#include <dims.h>            ! nlat2
#include <stptrs.h>          ! u_p
#include <diagptrs.h>
#include <planet.h>          ! now includes phibar (7/30/05)

      real*4 phy_z(nlon,nlat2,0:1,plev)
      real*4 phy_d(nlon,nlat2,0:1,plev)
      real*4 phy_u(nlon,nlat2,0:1,plev)
      real*4 phy_v(nlon,nlat2,0:1,plev)
      real*4 phy_h(nlon,nlat2,0:1,plev)
      real*4 phy_q(nlon,nlat2,0:1,plev)
c     rks - tracer field:
      real*4 phy_s(nlon,nlat2,0:1,plev)

      REAL_TYPE phy(nlon,nlat2,nphy,0:1,plev)
      REAL_TYPE cslat(nlat2),snlat(nlat2)
      real*4 u,v,h,za,zr,fcor

      INT_TYPE i,j,k,ins

      do k=1,plev
         do ins=0,1
            do j=1,nlat2
               fcor=(-1.0)**ins*real(2.0d0*Omega*snlat(j))
               do i=1,nlon
c
                  u    =real(phy(i,j,u_p    ,ins,k))/real(cslat(j))
                  v    =real(phy(i,j,v_p    ,ins,k))/real(cslat(j))
                  za   =real(phy(i,j,vor_p,ins,k))
                  zr   =za-fcor
                  h    =real(phy(i,j,theta_p,ins,k)/PhiBar)

                  phy_z(i,j,ins,k) = zr
                  phy_u(i,j,ins,k) = u
                  phy_v(i,j,ins,k) = v
                  phy_d(i,j,ins,k) = real(phy(i,j,div_p,ins,k))
                  phy_h(i,j,ins,k) = h
c                  phy_q(i,j,ins,k) = za/(1.0+h)
                  phy_q(i,j,ins,k) = za/((1.0+h)*(2.d0*Omega))

                  phy_s(i,j,ins,k) = real(phy(i,j,sig_p,ins,k))

               end do
            end do
         end do
      end do

      end
