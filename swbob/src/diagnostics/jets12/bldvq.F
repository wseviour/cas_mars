c      subroutine bldvq(phy,
c     $                   zon_node,
c     $                   cslat,snlat,
c     $                   nzavg,rlonmax)
c      implicit none

      subroutine bldvq(zon_node,nzavg)
      implicit none

#include <dims.h>
#include <stptrs.h>
#include <diagptrs.h>
#include <planet.h>

      INT_TYPE  nzavg
c      REAL_TYPE phy(nlon,nlat2,nphy,0:1,plev)
      REAL_TYPE zon_node(nlat2,0:1,plev,nzavg)
c      REAL_TYPE cslat(nlat2),snlat(nlat2),rlonmax

c      REAL_TYPE fcor,avor,h,v,vp,qp

      INT_TYPE i,j,k,ins

      do k=1,plev
         do ins=0,1
            do j=1,nlat2
               zon_node(j,ins,k,vqzon_p) = zon_node(j,ins,k,vqzon_p)
     .                          -zon_node(j,ins,k,vzon_p)*zon_node(j,ins,k,qzon_p) 
            enddo
         enddo
      enddo

c     following was wrong (was called for each longitudinal tile)
#if 0
      do k=1,plev
         do ins=0,1
            do j=1,nlat2
               fcor=(-1.0)**ins*2.0d0*Omega*snlat(j)
               do i=1,nlon
                  avor =phy(i,j,vor_p  ,ins,k)
                  h    =phy(i,j,theta_p,ins,k)/Phibar
                  v    =phy(i,j,v_p    ,ins,k)/cslat(j)

                  qp   =avor/(1.d0+h)-zon_node(j,ins,k,qzon_p)
                  vp   =v            -zon_node(j,ins,k,vzon_p)
                   
                  zon_node(j,ins,k,vqzon_p)   = zon_node(j,ins,k,vqzon_p)   + vp*qp

               end do
            end do
         end do
      end do

      do k=1,plev
         do ins=0,1
            do j=1,nlat2
               zon_node(j,ins,k,vqzon_p)   = zon_node(j,ins,k,vqzon_p) * rlonmax
            end do
         end do
      end do

#endif

      end 

















