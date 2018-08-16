      subroutine initdiag(navg)
      implicit none

#include <diag.h>

      INT_TYPE navg

      ! zero averaging counter

      navg=0

      ! zero accumulators

      call zeroreal(zon,plev*plat*nzavg)

      end
