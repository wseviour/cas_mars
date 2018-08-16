        subroutine synch()
        implicit none

#include <multinode.h>

#if defined(MPI)
        INT_TYPE ierr
        call MPI_barrier(comm,ierr)
#endif

       end
