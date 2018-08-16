      subroutine zavg_unpack(zone_node,zon,nlev,ipe)

#include <dims.h>
#include <coord.h>
#include <decomp.h>

      INT_TYPE  nlev
      REAL_TYPE zon(plat,nlev)
      REAL_TYPE zone_node(nlat2,0:1,nlev,jelemd)
      INT_TYPE  ipe

      INT_TYPE je,jeg,j,k


      do k=1,nlev               
         jptr=jbeg(ipe)
         do jeg=jebeg(ipe),jeend(ipe) ! loop over elements
            je=jeg-jebeg(ipe)+1
            do j=1,nlat2
               zon(jptr,k)=zone_node(j,0,k,je)
               zon(plat + 1 -jptr,k)=zone_node(j,1,k,je)
               jptr=jptr+1
            end do
         end do
      end do
 
      end

      
