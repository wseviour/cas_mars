#ifndef MULTINODE_H_
#define MULTINODE_H_

#include <params.h>
#include <type.h>

c
c...multinode.h:: common for multinode MPI latitude decomposition
c
      INT_TYPE rank       ! node id
      INT_TYPE root       ! root node id
      INT_TYPE comm       ! mpi global communicator 
      INT_TYPE real_t     ! real mpi data type
      INT_TYPE realsp_t   ! real single precision formovies mpi data type
      INT_TYPE int_t      ! INT_TYPE mpi data type
      INT_TYPE char_t     ! character mpi data type
      INT_TYPE pack_t     ! packed mpi data type
      INT_TYPE bool_t     ! logical mpi data type
      INT_TYPE stat_s     ! status size
      INT_TYPE sum_oper     ! sum reduce operation
      INT_TYPE max_oper     ! max reduce operation
      common /multinode/ rank,root,comm,
     .     real_t,realsp_t,int_t,char_t,pack_t,bool_t,stat_s,sum_oper,max_oper
#endif      
