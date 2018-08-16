      subroutine equibld(phy,
     $                   unode,vnode,hnode,znode,
     $                   cslat,snlat)

      implicit none

#include <dims.h>
#include <stptrs.h>
#include <planet.h>

      REAL_TYPE phy(nlon,nlat2,nphy,0:1,plev)
      REAL_TYPE unode(nlat2,0:1,plev)
      REAL_TYPE vnode(nlat2,0:1,plev)
      REAL_TYPE hnode(nlat2,0:1,plev)
      REAL_TYPE znode(nlat2,0:1,plev)
      REAL_TYPE cslat(nlat2),snlat(nlat2)

      REAL_TYPE fcor,u,v,h,z

      INT_TYPE i,j,k,ins

      do k=1,plev
         do ins=0,1
            do j=1,nlat2
               fcor=(-1.0)**ins*2.0d0*Omega*snlat(j)
               do i=1,nlon
                  u    =phy(i,j,u_p    ,ins,k)/cslat(j)
                  v    =phy(i,j,v_p    ,ins,k)/cslat(j)
                  h    =phy(i,j,theta_p,ins,k)/Phibar
                  z    =phy(i,j,vor_p  ,ins,k)-fcor
                   
                  unode(j,ins,k)   = unode(j,ins,k)   + u
                  vnode(j,ins,k)   = vnode(j,ins,k)   + v
                  hnode(j,ins,k)   = hnode(j,ins,k)   + h
                  znode(j,ins,k)   = znode(j,ins,k)   + z

               enddo
            enddo
         enddo
      enddo

      end 

















