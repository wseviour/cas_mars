#include <type.h>
      subroutine scdebug(sc,m)
      implicit none

!--------------------------------------------------------------------
!      Include Files
!--------------------------------------------------------------------
#include <dims.h>
#include <stptrs.h>

      INT_TYPE  m
      REAL_TYPE sc(0:1,nsp,m:nn,plev)

      INT_TYPE ic  ! re/im index
      INT_TYPE n   ! legendre order number
      INT_TYPE k   ! level index
      INT_TYPE l   ! state variable index

c----------------------------------------
c     Executable Statements
c----------------------------------------

      do k=1,plev
         do n=m,nn

            write(*,1) 'm,n,sc',m,n,sc(0,vor_p,n,k)

         end do ! Loop in n
      end do     ! end loop over levels

 1    format(a10,2(1x,i4),1x,e12.4)
      end




