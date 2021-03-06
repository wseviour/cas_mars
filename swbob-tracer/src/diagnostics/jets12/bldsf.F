      subroutine bldsf(phy,phvor,
     $     physf_node_z,physf_node_d,physf_node_h,momq,momc,Eph,Zph,
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

      REAL_TYPE momq(3,3),momc(3,3),Eph,Zph,nspm

      REAL_TYPE phy(nlon,nlat2,nphy,0:1,plev)
      REAL_TYPE phvor(nlon,nlat2,nphy,0:1,plev)
      REAL_TYPE cslat(nlat2),snlat(nlat2),da(nlat2)
      real*4 u,v,h,za,zr,fcor

      INT_TYPE i,j,k,ins,im,is

      do k=1,plev
         do ins=0,1
            nspm=(-1.0)**ins
            do j=1,nlat2
               fcor=(-1.0)**ins*real(2.0d0*Omega*snlat(j))
               do i=1,nlon
c                  
                  u    =real(phy(i,j,u_p    ,ins,k))/real(cslat(j))
                  v    =real(phy(i,j,v_p    ,ins,k))/real(cslat(j))
                  za   =real(phy(i,j,vor_p,ins,k))
                  zr   =za-fcor
                  h    =real(phy(i,j,theta_p,ins,k))/PhiBar
c
                  physf_node_z(i,j,ins,k) = zr
                  physf_node_d(i,j,ins,k) = real(phy(i,j,div_p,ins,k))
                  physf_node_h(i,j,ins,k) = h
c                  physf_node_q(i,j,ins,k) = za/(1.d0+h)

                  do im=1,3
                  do is=1,3
                  momq(im,is)=momq(im,is)+phvor(i,j,is,ins,k)**im*da(j)
                  momc(im,is)=momc(im,is)+(nspm*phvor(i,j,is,ins,k))**im*da(j)
                  enddo
                  enddo
                  Eph=Eph+(u**2+v**2)*da(j)
                  Zph=Zph+zr**2*da(j)

               end do
            end do
         end do
      end do

      end 

















