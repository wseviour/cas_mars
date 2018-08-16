/* 
 * pair.c:: send/receive pair generator
 *
 *    size:  number of nodes
 *    rmax:  rank upper bound for pairing (1..size)
 *    rank:  rank number (0..size-1)
 *    stage: pair index (1..size-1)
 */

#include <params.h>

#if defined(MPI)
#include <mpi.h>
#else
#define MPI_PROC_NULL 0
#endif

int pair(int *size,int *rmax,int *rank,int *stage){
  int partner;
  partner = (*rank)^(*stage);
  if((partner>*size-1)||(*rank>=*rmax))
    return(MPI_PROC_NULL);
  else
    return(partner);
}

int pair_(int *size,int *rmax,int *rank,int *stage){
    int partner;
    partner=pair(size,rmax,rank,stage);
    return(partner);
}
