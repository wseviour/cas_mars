      subroutine zavg_diag(itime,icount,nacc)
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

c     calculate \bar{v'q'}=\bar{vq}-\bar{v}\bar{q}:
c      do jeg=jebeg(rank),jeend(rank)
c         je=jeg-jebeg(rank)+1
c         call bldvq(zon_node(1,je),nzavg)
c      end do

c     calculate \bar{v'q'}:   WRONG - can't divide by rlonmax inside ie loop
c     
c      do jeg=jebeg(rank),jeend(rank)
c         je=jeg-jebeg(rank)+1
c         rlonmax=1.0d0/lonmax(je)
c         do ie=1,ielem(jeg)
c            call bldvq(phy(1,1,1,ie,je),
c     $                   zon_node(1,je),
c     $                   cslat(1,jeg),snlat(1,jeg),
c     $                   nzavg,rlonmax)
c         end do
c      end do

c     accumulate:
      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1
         call axpy(nlat2*2*plev*nzavg,1.d0,zon_node(1,je),1,zon_node_acc(1,je),1)
      end do

#ifdef ADAPT
      if (mod(time_current,time_out).lt.(dt2/2.d0)) then
#else
      if (mod(itime,zavgfreq).eq.0) then
#endif
c     gather them together on the root node and output:
         accinv=1.d0/float(nacc)
         do jeg=jebeg(rank),jeend(rank)
            je=jeg-jebeg(rank)+1
            call scal(nlat2*2*plev*nzavg, accinv, zon_node_acc(1,je), 1)
         end do
         call zavg_pgather(zon_node_acc, zon, plev*nzavg)
         if (rank.eq.root) call diagout
         call zeroreal(zon_node_acc ,jelemd*nlat2*2*plev*nzavg)
      endif

      end


