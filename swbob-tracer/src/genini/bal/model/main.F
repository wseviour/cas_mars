
       program main
       implicit none

#include <params.h>
#include <mpcluster.h>
#include <spectral.h>
#include <fields.h>
#include <io.h>

      ! ===========================================================
      ! Start Multi-Threads Under Multi-Process code Parallel Region
      ! This *must* be the first thing we do.
      ! ===========================================================
      call initmp(P_NODE,rank,root,nThreads) 
      ! ===========================================================
      ! Initialize the model
      ! ===========================================================
      call initmodel()
      ! ===========================================================
      ! Perform one spherical harmonic transform cycle
      ! ===========================================================
      call sht()
      call output()
      ! ===========================================================
      ! End Multi-Threads Under Multi-Processes Parallel region.
      ! This *must* be the last thing we do.
      ! ===========================================================
      call haltmp('end of BOB')

      end










