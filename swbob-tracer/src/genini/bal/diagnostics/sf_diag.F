      subroutine sf_diag()
      implicit none

#include <fields.h>          ! phy
#include <physical.h>        ! cslat(nlat2,jelemg),ielem(jelemg)
#include <decomp.h>          ! jebeg(0:P_NODE-1),jeend(0:P_NODE-1)
#include <mpcluster.h>       ! rank

      real*4 physf_node_z(nlon*nlat2*2*plev,ielemd,jelemd)    
      real*4 physf_node_d(nlon*nlat2*2*plev,ielemd,jelemd)    
      real*4 physf_node_h(nlon*nlat2*2*plev,ielemd,jelemd)    
      real*4 phyg_z(plon,plat,plev)
      real*4 phyg_d(plon,plat,plev)
      real*4 phyg_h(plon,plat,plev)
      INT_TYPE  je,jeg,ie,icountp,j

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Executable statements
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1
         do ie=1,ielem(jeg)
            call bldsf(phy(1,1,1,ie,je),
     $           physf_node_z(1,ie,je),
     $           physf_node_d(1,ie,je),
     $           physf_node_h(1,ie,je),
     $                   cslat(1,jeg),
     $                   snlat(1,jeg))

         end do
      end do

      ! gather them together on the root node:

      call sf_pgather(physf_node_z, phyg_z, plev)
      call sf_pgather(physf_node_d, phyg_d, plev)
      call sf_pgather(physf_node_h, phyg_h, plev)

      if(rank.eq.root)then
c         print*,'calling congen...'
c         call congen(phyg_q,icountp)
c         print*,'done congen'
         call sfout(phyg_z,phyg_d,phyg_h)
      end if

      end


