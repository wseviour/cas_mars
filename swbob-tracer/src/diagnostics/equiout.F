      subroutine equiout
      implicit none

#include <coord.h>
#include <io.h>
#include <equilib.h>

      character*80 zavgfile
      character*60 stub
      character*4  tail,tail_lj

      INT_TYPE  ounit
      INT_TYPE  i,j,k
      real*4 ubarsp(plat,plev)
      real*4 vbarsp(plat,plev)
      real*4 hbarsp(plat,plev)
      real*4 zbarsp(plat,plev)
      
      do k=1,plev
         do j=1,plat
            ubarsp(j,k)=real(ubar(j,k))
            vbarsp(j,k)=real(vbar(j,k))
            hbarsp(j,k)=real(hbar(j,k))
            zbarsp(j,k)=real(zbar(j,k))
         enddo
      enddo

c     Files units are opened in initio.F
      call writespreal(ubarsp,plev*plat,36)
      call writespreal(vbarsp,plev*plat,37)
      call writespreal(hbarsp,plev*plat,38)
      call writespreal(zbarsp,plev*plat,39)

      do i=36,39
         call flush(i)
      enddo

      end
