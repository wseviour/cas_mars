      subroutine pmaxreal(s,r,len)
      implicit none

!--------------------------------------------------------------------
!     Include Files
!--------------------------------------------------------------------

#include <type.h>
#include <multinode.h>

      INT_TYPE  len
      REAL_TYPE s(len)
      REAL_TYPE r(len)
 
!--------------------------------------------------------------------
!
!     Local Variables
!
!--------------------------------------------------------------------

#if defined(MPI)
      INT_TYPE ierr,i
      call MPI_reduce(s,
     $                r,
     $                len,
     $                real_t,
     $                max_oper,
     $                root,
     $                comm,
     $                ierr)
      do i=1,len
         s(i)=r(i)
      end do
#endif

         end 
