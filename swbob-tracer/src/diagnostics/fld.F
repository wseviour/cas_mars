      subroutine fld(icountp)
      implicit none

#include <fields.h>          ! phy
#include <physical.h>        ! cslat(nlat2,jelemg),ielem(jelemg)
#include <decomp.h>          ! jebeg(0:P_NODE-1),jeend(0:P_NODE-1)
#include <mpcluster.h>       ! rank
#include <spectral.h>        ! sc(2*(plev*nsp)*nsc)
#include <modeltime.h>
#include <timestep.h>

      real*4 phy_z(nlon*nlat2*2*plev,ielemd,jelemd)
      real*4 phy_d(nlon*nlat2*2*plev,ielemd,jelemd)
      real*4 phy_u(nlon*nlat2*2*plev,ielemd,jelemd)
      real*4 phy_v(nlon*nlat2*2*plev,ielemd,jelemd)
      real*4 phy_h(nlon*nlat2*2*plev,ielemd,jelemd)
      real*4 phy_q(nlon*nlat2*2*plev,ielemd,jelemd)
c     rks- tracer field
      real*4 phy_s(nlon*nlat2*2*plev,ielemd,jelemd)
      real*4 phyg_q(plon,plat,plev)
      INT_TYPE  je,jeg,ie,icountp,j
      character*1 lfld

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Executable statements
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1
         do ie=1,ielem(jeg)
            call fld_bld(phy(1,1,1,ie,je),
     $           phy_u(1,ie,je),
     $           phy_v(1,ie,je),
     $           phy_z(1,ie,je),
     $           phy_d(1,ie,je),
     $           phy_h(1,ie,je),
     $           phy_q(1,ie,je),
     $           phy_s(1,ie,je),
     $                   cslat(1,jeg),
     $                   snlat(1,jeg))
         enddo
      enddo

c      lfld='z'
c      call fld_out(phy_z,lfld,icountp)
c      lfld='d'
c      call fld_out(phy_d,lfld,icountp)
      lfld='h'
      call fld_out(phy_h,lfld,icountp)
      lfld='q'
      call fld_out(phy_q,lfld,icountp)
      lfld='s'
      call fld_out(phy_s,lfld,icountp)
      lfld='u'
      call fld_out(phy_u,lfld,icountp)
      lfld='v'
      call fld_out(phy_v,lfld,icountp)

c     gather PV fields on the root node and contour:
      call sf_pgather(phy_q, phyg_q, plev)
      if(rank.eq.root)then
         write(*,*),'calling congen...'
         call congen(phyg_q,icountp)
         write(*,*),'done congen'
      endif



      end
