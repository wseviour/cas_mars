      subroutine sf_unpack(physf_node,phyg,ipe,lplev)

#include <dims.h>
#include <coord.h>
#include <decomp.h>  !jebeg(0:P_NODE-1),jeend(0:P_NODE-1),jbeg(0:P_NODE-1),jend(0:P_NODE-1)
#include <physical.h> !ielem(jelemg)
#include <io.h>       ! MPI_UNIT

      INT_TYPE lplev

      real*4 phyg(plon,plat,lplev)
      real*4 physf_node(nlon,nlat2,0:1,lplev,ielemd,jelemd)
      INT_TYPE  ipe

      INT_TYPE je,jeg,j,i,ie,k
      INT_TYPE jptr,iptr

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Executable statemnets
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do k=1,lplev
         do jeg=jebeg(ipe),jeend(ipe) 
            je=jeg-jebeg(ipe)+1
            do ie=1,ielem(jeg)
               iptr=(ie-1)*nlon
               do j=1,nlat2
                  jptr=j+(jbeg(ipe)-1)+(je-1)*nlat2
                  do i=1,nlon
                     
                     phyg(iptr+i,jptr,k)=physf_node(i,j,0,k,ie,je)
                     phyg(iptr+i,plat + 1 -jptr,k)=physf_node(i,j,1,k,ie,je)
                     
                  end do
               end do
            end do
         end do
      end do
         
      end

      
