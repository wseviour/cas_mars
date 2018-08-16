#define ZAVG_TAG_BASE 1

! assemble global zonal average from nodal zonal average, and scale it.

      subroutine znl_pgather(zon_node,
     $                        zon,
     $                        nlev)
      implicit none

#include <dims.h>
#include <mpcluster.h>

      INT_TYPE nlev
      REAL_TYPE zon_node(nlat2*2*nlev,jelemd)
      REAL_TYPE zon(plat,nlev)

      INT_TYPE  ipe

      INT_TYPE  rtag,stag

      ! send zon_node to root process

      if (rank.eq.root) then

         do ipe=0,P_NODE-1

            if (ipe.ne.root) then
               rtag=ipe+ZAVG_TAG_BASE*P_NODE
               call recvreal(zon_node,nlat2*2*jelemd*nlev,ipe,rtag)
            end if

            call znl_unpack(zon_node,zon,nlev,ipe)
            
         end do

      else

         stag=rank+ZAVG_TAG_BASE*P_NODE
         call sendreal(zon_node,nlat2*2*jelemd*nlev,root,stag)

      end if

      end


