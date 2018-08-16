c
c ... haltmp.F:: graceful exit from message passing environment
c     Stands for halt mpi
c     MPI_finalize exit MPI gracefully

      subroutine haltmp(string)
      implicit none

#include <params.h>
#include <type.h>
#include <mpcluster.h>

#if defined(MPI)
      INT_TYPE info
#endif

      character*(*) string
      if (rank.eq.root) write(*,*) string

#if defined(MPI)
      call MPI_finalize(info)
#endif

      stop
      end

