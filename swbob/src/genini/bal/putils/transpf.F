c
c...transpf.F:: forward spectral transpose
c
c       fc(len,p) --> fcx(len,n)  where n >= p
c
      subroutine transpf(fc,fcx,len,p,n)
      implicit none

#include <multinode.h>
c
c...arguments
c
      INT_TYPE len      ! length of each slice
      INT_TYPE p        ! # of spectral processors
      INT_TYPE n        ! # of processors
      REAL_TYPE fc(len,0:p-1)    ! input array
      REAL_TYPE fcx(len,0:n-1)   ! output array
c
c...local variables      
c
      INT_TYPE i

#if defined(MPI)

      INT_TYPE stage
      INT_TYPE info

      INT_TYPE send,recv
      INT_TYPE pair,ceil2
      external pair,ceil2

#ifdef LAM_MPI
#include <mpif.h>       
      INT_TYPE status(MPI_STATUS_SIZE)
#else
      INT_TYPE status(stat_s)
#endif

#endif

      ! On processor part

      call synch()

      if(rank.lt.p) then
         do i=1,len
            fcx(i,rank) = fc(i,rank)
         enddo
      endif

#if defined(MPI)
#define SENDRECV_FUSED
      do stage=1,ceil2(n)-1
         send = pair(p,n,rank,stage)
         recv = pair(n,p,rank,stage)
#if defined(SENDRECV_FUSED)
         call mpi_sendrecv(fc(1,send),len,real_t,send,stage,
     .        fcx(1,recv),len,real_t,recv,stage,comm,status,info)
#endif
#if defined(SENDRECV_PAIR)
         if (rank.gt.send) then
            call MPI_send(fc(1,send),len,real_t,send,stage,comm,info)
            call MPI_recv(fcx(1,recv),len,real_t,recv,stage,comm,status,info)
         else
            call MPI_recv(fcx(1,recv),len,real_t,recv,stage,comm,status,info)
            call MPI_send(fc(1,send),len,real_t,send,stage,comm,info)
         end if
#endif
      enddo
#endif

      return
      end


      


