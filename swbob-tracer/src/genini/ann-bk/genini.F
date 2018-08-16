      program genini

      implicit none

#include <coord.h>
#include <planet.h>
#include <mpcluster.h>
#include <io.h>

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      INT_TYPE   ierr
      REAL_TYPE  plon_file,plev_file,time_start
      INT_TYPE   ounit

      ! Returns index of next available fortran unit
      INT_TYPE navu
      external navu

      INT_TYPE  lenchr
      external  lenchr

      character*80 stub
      character*5 tail,tail_lj
      character*80 proclogfile
      character*80 filename

c-------------------------------------------------------------
c Begin executable statements
c-------------------------------------------------------------

      call initmp(P_NODE,rank,root,nThreads) ! initialize multiprocessor cluster info
      call gridfile()                        ! Read and broadcast the grid

c     Create name and open log file for each processor

      stub="./out_genini"
      write(tail,1)rank
 1    format(i4)
      call ljust(tail,tail_lj)               ! Remove leading blanks
      proclogfile=stub(1:lenchr(stub))//'.'//tail_lj(1:lenchr(tail_lj))

      open(unit=MPI_UNIT,
     $     file=proclogfile,
     $     form='formatted')

      write(MPI_UNIT,*) "---------------------------"
      write(MPI_UNIT,*) "Log File for Process",rank
      write(MPI_UNIT,*) "---------------------------"

c     Initialize the spectral transform apparatus
c     Decomp,FFTs, Coriolis factor, annp1(0:nn+1)

      call initsht(colatitude,
     $             gausswt,
     $             plonmax,
     $             Alpha_coord,
     $             Radius,
     $             Omega)

c     Rank 0 output the Initial condition in progsc.ini

      if (rank.eq.root) then
         filename='progsc.ini'
         ounit=navu()
         open(UNIT   = ounit,
     $        FILE   = filename,
     $        ACCESS = 'sequential',
     $        STATUS = 'unknown',
     $        FORM   = 'unformatted',
     $        IOSTAT = ierr)

         plon_file=PLON
         plev_file=PLEV
         time_start=0.0D0
         write(ounit)plon_file,plev_file,time_start
      end if

c     Compute the initial condition
c     Output spec coef of IC to progsc.ini unit ounit

      call initscg(ounit)

      if (rank.eq.root) then
         close(ounit)
      end if

c     Exit MPI gracefully (MPI_finalize)
      call haltmp('End of Genini')

      stop
      end
