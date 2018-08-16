        subroutine recvspreal(buf,len,srcpe,tag)
        implicit none

#include <multinode.h>

        INT_TYPE  len
        real buf(len)
        INT_TYPE  srcpe       ! pe source of message
        INT_TYPE  tag

#if defined(MPI)

      INT_TYPE  ierr

#ifdef LAM_MPI
#include <mpif.h>       
      INT_TYPE status(MPI_STATUS_SIZE)
#else
      INT_TYPE status(stat_s)
#endif

      call MPI_recv(buf,
     $              len,
     $              realsp_t,
     $              srcpe,
     $              tag,
     $              comm,
     $              status,
     $              ierr)

#endif

       end
