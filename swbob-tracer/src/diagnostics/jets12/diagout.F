      subroutine diagout
      implicit none

#include <diag.h>
#include <diagptrs.h>
#include <coord.h>
#include <io.h>

      character*80 zavgfile
      character*60 stub
      character*4  tail,tail_lj

      INT_TYPE  ounit
      INT_TYPE  i,j,k
      real*4 zonsp (plat,plev,nzavg)
      
      do i=1,nzavg
         do k=1,plev
            do j=1,plat
               zonsp(j,k,i)=real(zon(j,k,i))
            enddo
         enddo
      enddo

c     Files units are opened in initio.F
      call writespreal(zonsp(1,1,uzon_p),plev*plat,26)
      call writespreal(zonsp(1,1,vzon_p),plev*plat,27)
      call writespreal(zonsp(1,1,zzon_p),plev*plat,28)
      call writespreal(zonsp(1,1,hzon_p),plev*plat,29)
      call writespreal(zonsp(1,1,qzon_p),plev*plat,30)
      call writespreal(zonsp(1,1,kzon_p),plev*plat,31)
      call writespreal(zonsp(1,1,vqzon_p),plev*plat,32)

      end
