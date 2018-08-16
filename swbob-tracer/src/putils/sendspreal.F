        subroutine sendspreal(buf,len,destpe,tag)
        implicit none

#include <multinode.h>

        INT_TYPE  len
        real buf(len)
        INT_TYPE  destpe         ! pe destination of send
        INT_TYPE  tag

#if defined(MPI)

        INT_TYPE  ierr

        call MPI_send(buf,
     $                len,
     $                realsp_t,
     $                destpe,
     $                tag,
     $                comm,
     $                ierr)

#endif

       end

