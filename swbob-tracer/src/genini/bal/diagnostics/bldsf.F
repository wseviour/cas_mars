      subroutine bldsf(phy,
     $     physf_node_z,physf_node_d,physf_node_h,
     $     cslat,snlat)

c     CALLED BY: sf_diag.F

      implicit none

#include <diag.h>            ! nzavg
#include <dims.h>            ! nlat2
#include <stptrs.h>          ! u_p
#include <diagptrs.h>
#include <planet.h>

      real*4 physf_node_z(nlon,nlat2,0:1,plev)
      real*4 physf_node_d(nlon,nlat2,0:1,plev)
      real*4 physf_node_h(nlon,nlat2,0:1,plev)

      REAL_TYPE phy(nlon,nlat2,nphy,0:1,plev)
      REAL_TYPE cslat(nlat2),snlat(nlat2)
      real*4 u,v,avor,div,phip,fcor,fac

      INT_TYPE i,j,k,ins

      fac=secperday*phibar

      do k=1,plev
         do ins=0,1
            do j=1,nlat2
               fcor=(-1.0)**ins*real(2.0d0*Omega*snlat(j))
               do i=1,nlon
c                  
                  avor=real(phy(i,j,vor_p,ins,k))
                  phip=real(phy(i,j,theta_p,ins,k))
                  div=real(phy(i,j,div_p,ins,k))
c
                  physf_node_z(i,j,ins,k) = avor-fcor
c                  physf_node_z(i,j,ins,k) = avor
                  physf_node_d(i,j,ins,k) = div
                  physf_node_h(i,j,ins,k) = phip
c                  physf_node_q(i,j,ins,k) = fac*avor/(phibar+phip)
c
               end do
            end do
         end do
      end do

      end 

















