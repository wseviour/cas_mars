      subroutine bldsf(phy,
     $     physf_node_z,physf_node_d,physf_node_h,
     $     cslat,snlat,da)

c     CALLED BY: sf_diag.F

      implicit none

#include <dims.h>            ! nlat2
#include <stptrs.h>          ! u_p
#include <diagptrs.h>
#include <planet.h>          ! now includes phibar (7/30/05)

      real*4 physf_node_z(nlon,nlat2,0:1,plev)
      real*4 physf_node_d(nlon,nlat2,0:1,plev)
      real*4 physf_node_h(nlon,nlat2,0:1,plev)

      REAL_TYPE phy(nlon,nlat2,nphy,0:1,plev)
      REAL_TYPE cslat(nlat2),snlat(nlat2),da(nlat2)
      real*4 u,v,h,za,zr,fcor

      INT_TYPE i,j,k,ins

      do k=1,plev
         do ins=0,1
            do j=1,nlat2
               fcor=(-1.0)**ins*real(2.0d0*Omega*snlat(j))
               do i=1,nlon
c                  
c                  u    =real(phy(i,j,u_p    ,ins,k))/real(cslat(j))
c                  v    =real(phy(i,j,v_p    ,ins,k))/real(cslat(j))
                  za   =real(phy(i,j,vor_p,ins,k))
                  zr   =za-fcor
                  h    =real(phy(i,j,theta_p,ins,k))/PhiBar
c
                  physf_node_z(i,j,ins,k) = zr
                  physf_node_d(i,j,ins,k) = real(phy(i,j,div_p,ins,k))
                  physf_node_h(i,j,ins,k) = h
c                  physf_node_q(i,j,ins,k) = za/(1.d0+h)

               end do
            end do
         end do
      end do

      end 

















