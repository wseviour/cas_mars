c
c...initmp.F:: initialize message passing environment
c
      subroutine initmp(nproc_req, rank_ret, root_ret, nThreads)
      implicit none

#include <multinode.h>

      INT_TYPE nproc_req     ! number of processes requested
      INT_TYPE rank_ret      ! rank of this process
      INT_TYPE root_ret      ! root process 
      INT_TYPE nThreads      ! number of threads

#if defined(OMP)
      INT_TYPE omp_get_num_threads
      external omp_get_num_threads
#endif

      INT_TYPE nproc         ! number of processes started by MPI
      
#if defined(MPI)


      INT_TYPE info         
      logical running        ! state of MPI at beginning of initmp call

      !
      !     Basic MPI initialization
      !


#include <mpif.h>

      call MPI_initialized(running,info)

      if (.not.running) then
         call MPI_init(info)
      end if

      call MPI_comm_rank(MPI_COMM_WORLD,rank,info)
      call error_dup(rank)
      call MPI_comm_size(MPI_COMM_WORLD,nproc,info)
 
      !
      ! Check to see if MPI starts the same number of processes as
      ! The model has requested. Abort if they don't
      !

      root     = 0
      rank_ret = rank
      root_ret = root
#else
      nproc    = 1
      rank     = 0
      root     = 0 
      rank_ret = 0
      root_ret = 0
#endif

      if (nproc.ne.nproc_req) then
         call haltmp('initmp:: ABORT - requested and MPI process counts differ')
      endif

#if defined(MPI)

      !
      ! Set MPI Constants
      !

      comm = MPI_COMM_WORLD

#ifdef LAM_MPI
      real_t = MPI_DOUBLE_PRECISION
      realsp_t = MPI_REAL4
#else
      real_t = MPI_REAL8
      realsp_t = MPI_REAL4
#endif

      int_t  = MPI_INTEGER
      char_t = MPI_CHARACTER
      pack_t = MPI_PACKED
      bool_t = MPI_LOGICAL
      stat_s = MPI_STATUS_SIZE
      sum_oper = MPI_SUM
      max_oper = MPI_MAX
#endif

      !
      ! Find out how many threads there are
      !

#if defined(OMP)
!$OMP PARALLEL
      nThreads=omp_get_num_threads()
!$OMP END PARALLEL            
#else
      nThreads=1
#endif

      return
      end
