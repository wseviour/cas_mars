      subroutine znl_bld(phy,
     $                   zon_node,
     $                   cslat,snlat,
     $                   nzavg)
      implicit none

#include <dims.h>
#include <stptrs.h>
#include <diagptrs.h>
#include <planet.h>

      INT_TYPE  nzavg
      REAL_TYPE phy(nlon,nlat2,nphy,0:1,plev)
      REAL_TYPE zon_node(nlat2,0:1,plev,nzavg)
      REAL_TYPE cslat(nlat2),snlat(nlat2)

      REAL_TYPE fcor,avor,u,v,h,q

      INT_TYPE i,j,k,ins

      do k=1,plev
         do ins=0,1
            do j=1,nlat2
               fcor=(-1.0)**ins*2.0d0*Omega*snlat(j)
               do i=1,nlon
                  avor =phy(i,j,vor_p  ,ins,k)
                  u    =phy(i,j,u_p    ,ins,k)/cslat(j)
                  v    =phy(i,j,v_p    ,ins,k)/cslat(j)
                  h    =phy(i,j,theta_p,ins,k)/Phibar
                  q    =avor/(1.d0+h)
                   
                  zon_node(j,ins,k,uzon_p)   = zon_node(j,ins,k,uzon_p)   + u
                  zon_node(j,ins,k,vzon_p)   = zon_node(j,ins,k,vzon_p)   + v
                  zon_node(j,ins,k,zzon_p)   = zon_node(j,ins,k,zzon_p)   + (avor-fcor)
                  zon_node(j,ins,k,hzon_p)   = zon_node(j,ins,k,hzon_p)   + h
                  zon_node(j,ins,k,qzon_p)   = zon_node(j,ins,k,qzon_p)   + q
                  zon_node(j,ins,k,kzon_p)   = zon_node(j,ins,k,kzon_p)   + 0.50d0*(u**2 + v**2)
                  zon_node(j,ins,k,vqzon_p)  = zon_node(j,ins,k,vqzon_p)  + v*q
c                  zon_node(j,ins,k,ptflux_p) = zon_node(j,ins,k,ptflux_p) + v*theta
c                  zon_node(j,ins,k,amflux_p) = zon_node(j,ins,k,amflux_p) + u*v*cslat(j)

               end do
            end do
         end do
      end do

      end 

















