#include <type.h>
      REAL_TYPE function eddyke()
      implicit none

#include <dims.h>
#include <diag.h>
#include <diagptrs.h>
#include <coord.h>

      INT_TYPE k,j,js
      REAL_TYPE dp

      dp=1.0D0/plev ! shallow water equations

      eddyke=0.0d0
      do k=1,plev
         do j=1,plat/2
            js=plat+1-j
            eddyke=eddyke + dp*gausswt(j)*(zon(j ,k,kezon_p)-0.50d0*(zon(j ,k,uzon_p)**2+zon(j ,k,vzon_p)**2)
     $                                   + zon(js,k,kezon_p)-0.50d0*(zon(js,k,uzon_p)**2+zon(js,k,vzon_p)**2))
         end do
      end do

      end
