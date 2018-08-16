      subroutine znl_bldvq(zon_node,nzavg)
      implicit none

#include <dims.h>
#include <stptrs.h>
#include <diagptrs.h>
#include <planet.h>

      INT_TYPE  nzavg
      REAL_TYPE zon_node(nlat2,0:1,plev,nzavg)

      INT_TYPE i,j,k,ins

      do k=1,plev
         do ins=0,1
            do j=1,nlat2
               zon_node(j,ins,k,vqzon_p) = zon_node(j,ins,k,vqzon_p)
     .                          -zon_node(j,ins,k,vzon_p)*zon_node(j,ins,k,qzon_p) 
            enddo
         enddo
      enddo

      end 

















