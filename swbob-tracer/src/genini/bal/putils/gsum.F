!
!     gsum.F:: global sum one word in place
!
      subroutine gsum(x)
      implicit none

!--------------------------------------------------------------------
!     Include Files
!--------------------------------------------------------------------
#include <multinode.h>
!--------------------------------------------------------------------

      REAL_TYPE x

#if defined(MPI)

      REAL_TYPE y
      INT_TYPE info

      call MPI_allreduce(x,y,1,real_t,sum_oper,comm,info)
      x=y

#endif

      return
      end

