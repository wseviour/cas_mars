c
c...transpb.F:: backward spectral transpose
c
c       fcx(len,n) --> fc(len,p)  where n >= p
c
      subroutine transpb(fcx,fc,len,n,p)
      implicit none
#include <type.h>
#include <multinode.h>
c
c...arguments
c
      INT_TYPE len      ! length of each slice
      INT_TYPE p        ! # of spectral processors
      INT_TYPE n        ! # of processors
      REAL_TYPE fcx(len,0:n-1)   ! input array
      REAL_TYPE fc (len,0:p-1)   ! output array

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

      call synch()

      ! On processor part
           
      if(rank.lt.p) then
         do i=1,len
            fc(i,rank) = fcx(i,rank)
         enddo
      endif

#if defined(MPI)
#define SENDRECV_FUSED
      do stage=1,ceil2(n)-1
         send = pair(n,p,rank,stage)
         recv = pair(p,n,rank,stage)
#if defined(SENDRECV_FUSED)
         call mpi_sendrecv(fcx(1,send),len,real_t,send,stage,
     .        fc(1,recv),len,real_t,recv,stage,comm,status,info)
#endif

#if defined(SENDRECV_PAIR)
         if (rank.gt.send) then
            call MPI_send(fcx(1,send),len,real_t,send,stage,comm,info)
            call MPI_recv(fc(1,recv),len,real_t,recv,stage,comm,status,info)
         else
            call MPI_recv(fc(1,recv),len,real_t,recv,stage,comm,status,info)
            call MPI_send(fcx(1,send),len,real_t,send,stage,comm,info)
         end if
#endif

      enddo
#endif

      return
      end


      


