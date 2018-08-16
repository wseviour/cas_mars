#define SF_TAG_BASE 1

      subroutine sf_pgather(physf_node, phyg, lplev)

      implicit none

#include <dims.h>         ! P_NODE from params.h
#include <mpcluster.h>    ! rank,root
#include <io.h>           ! MPI_UNIT

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Local Variables
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      INT_TYPE  lplev
      INT_TYPE  ipe
      INT_TYPE  rtag,stag,je,j
      real*4 physf_node (nlon*nlat2*2*lplev,ielemd,jelemd)    
      real*4 phyg (plon,plat,lplev)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Executable statemnets
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (rank.eq.root) then

         do ipe=0,P_NODE-1
            if (ipe.ne.root) then
               rtag=ipe+SF_TAG_BASE*P_NODE ! P_NODE is number of processes in params.h
               call recvspreal(physf_node,nlon*nlat2*2*lplev*ielemd*jelemd,ipe,rtag)
            end if
            call sf_unpack(physf_node,phyg,ipe,lplev)
         end do

      else

         stag=rank+SF_TAG_BASE*P_NODE
         print*,'stag = ', stag
         call sendspreal(physf_node,nlon*nlat2*2*lplev*ielemd*jelemd,root,stag)

      end if

      end
