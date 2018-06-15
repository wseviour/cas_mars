        subroutine sendreal(buf,len,destpe,tag)
        implicit none

#include <multinode.h>

        INT_TYPE  len
        REAL_TYPE buf(len)
        INT_TYPE  destpe         ! pe destination of send
        INT_TYPE  tag

#if defined(MPI)

        INT_TYPE  ierr

        call MPI_send(buf,
     $                len,
     $                real_t,
     $                destpe,
     $                tag,
     $                comm,
     $                ierr)

#endif

       end
