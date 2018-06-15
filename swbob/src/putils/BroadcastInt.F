      subroutine BroadcastInt(a,nwords)
      implicit none

!--------------------------------------------------------------------
!     Include Files
!--------------------------------------------------------------------
#include <multinode.h>
!--------------------------------------------------------------------


      INT_TYPE  nwords
      INT_TYPE  a(nwords)

#if defined(MPI)

      INT_TYPE info
      call MPI_Bcast(a,
     $               nwords,
     $               int_t,
     $               root,
     $               comm,
     $               info)
#endif

      end
