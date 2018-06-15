      subroutine initio
      implicit none

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include <io.h>
#include <mpcluster.h>
#include <dims.h>
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      INT_TYPE  navu
      external  navu

      INT_TYPE  lenchr
      external  lenchr

      INT_TYPE  ierr        ! file io error flag

      character*80 stub
      character*5 tail,tail_lj
      character*80 proclogfile

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Executable Statements
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      stub="./proclog_"
      write(tail,1)rank
 1    format(i4)
      call ljust(tail,tail_lj)
      proclogfile=stub(1:lenchr(stub))//'.'//tail_lj(1:lenchr(tail_lj))
      open(unit=MPI_UNIT,
     $     file=proclogfile,
     $     form='formatted')

      write(MPI_UNIT,*) "Log File for Process",rank

      end
