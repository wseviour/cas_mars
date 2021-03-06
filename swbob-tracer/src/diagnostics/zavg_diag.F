      subroutine zavg_diag
      implicit none

#include <diag.h>
#include <decomp.h>
#include <fields.h>
#include <physical.h>
#include <mpcluster.h>
#include <fft.h>
#include <modeltime.h>
#include <timestep.h>

      REAL_TYPE zon_node  (nlat2*2*plev*nzavg,jelemd)    ! nodal zonal averages
      INT_TYPE  je,jeg,ie
      REAL_TYPE rlonmax
      REAL_TYPE accinv
      INT_TYPE itime,icount,nacc

c     zero out nodal zonal average arrays:
      call zeroreal(zon_node  ,jelemd*nlat2*2*plev*nzavg)

c     construct operators and zonally average:
      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1
         do ie=1,ielem(jeg)
            call bldzavg(phy(1,1,1,ie,je),
     $                   zon_node(1,je),
     $                   cslat(1,jeg),snlat(1,jeg),
     $                   nzavg)
         end do
      end do

c     scale zon_node by 1/lonmax:
      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1
         rlonmax=1.0d0/lonmax(je)
         call scal(nlat2*2*plev*nzavg, rlonmax, zon_node(1,je), 1)
      end do

c     calculate \bar{v'q'}=\bar{vq}-\bar{v}\bar{q}:  (needs fixing)
c      do jeg=jebeg(rank),jeend(rank)
c         je=jeg-jebeg(rank)+1
c         call bldvq(zon_node(1,je),nzavg)
c      end do

c     gather them together on the root node and output:
      call zavg_pgather(zon_node, zon, plev*nzavg)
      if (rank.eq.root) then
         call diagout
      endif

      end


