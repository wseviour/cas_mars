      subroutine initio
      implicit none

#include <io.h>
#include <modeltime.h>
#include <mpcluster.h>
#include <dims.h>

      INT_TYPE  navu
      external  navu

      INT_TYPE  lenchr
      external  lenchr

      INT_TYPE  ierr        ! file io error flag
      character*80 msg      ! file io error message

      character*1 trunc
      character*80 stub
      character*5 tail,tail_lj
      character*80 proclogfile
      character*80 zavgfile

      ! scale I/O frequencies
      ! if rstfreq or zonfreq are negative, 
      ! interpret ABS value as number of days
      
      if (rstfreq.lt.0) then
         rstfreq = -nspday*rstfreq
      end if

      if (zonfreq.lt.0) then
         zonfreq = -nspday*zonfreq
         fldfreq = -nspday*fldfreq
      end if

      call broadcastint(rstfreq,1)
      call broadcastint(fldfreq,1)
      call broadcastint(zonfreq,1)
      call broadcastint(caseid,8)

c      stub="./proclog_"//caseid(1:lenchr(caseid))
      stub="./proclog_"
      write(tail,1)rank
 1    format(i4)
      call ljust(tail,tail_lj)
      proclogfile=stub(1:lenchr(stub))//'.'//tail_lj(1:lenchr(tail_lj))
      open(unit=MPI_UNIT,
     $     file=proclogfile,
     $     form='formatted')
      write(MPI_UNIT,*) "Log File for Process",rank

      write(MPI_UNIT,9000)
 9000 format("")
      write(MPI_UNIT,9010)
 9010 format("------------io control data---------------")
      write(MPI_UNIT,9030)rstfreq
 9030 format(" restart tape write frequency       = ",i8)
      write(MPI_UNIT,9040)zonfreq
 9040 format(" zoanl average write frequency       = ",i8)
      write(MPI_UNIT,9070)
 9070 format("-------------------------------------------")

      call flush(MPI_UNIT)

      if (rank.eq.root) then

c     open files for zonal averages:
c      stub=caseid(1:lenchr(caseid))//'.zavg'
      stub='zavg'
      zavgfile=stub(1:lenchr(stub))
      open(unit=26,
     $     file=zavgfile(1:lenchr(zavgfile))//'.u',
     $     form="UNFORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)

      open(unit=27,
     $     file=zavgfile(1:lenchr(zavgfile))//'.v',
     $     form="UNFORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)

      open(unit=28,
     $     file=zavgfile(1:lenchr(zavgfile))//'.z',
     $     form="UNFORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)

      open(unit=29,
     $     file=zavgfile(1:lenchr(zavgfile))//'.h',
     $     form="UNFORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)

      open(unit=30,
     $     file=zavgfile(1:lenchr(zavgfile))//'.q',
     $     form="UNFORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)

      open(unit=31,
     $     file=zavgfile(1:lenchr(zavgfile))//'.k',
     $     form="UNFORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)

      open(unit=32,
     $     file=zavgfile(1:lenchr(zavgfile))//'.vq',
     $     form="UNFORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)

c     open file for parameters:
      open(unit=81,
     $     file='paramt.dat',
     $     form="FORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)

c     open files for spectra:
      open(unit=85,
     $     file='Erot.dat',
     $     form="FORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)
      open(unit=86,
     $     file='Ediv.dat',
     $     form="FORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)
      open(unit=87,
     $     file='Etot.dat',
     $     form="FORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)
      open(unit=88,
     $     file='Ztot.dat',
     $     form="FORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)
      open(unit=89,
     $     file='urms.dat',
     $     form="FORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)
      open(unit=90,
     $     file='skwq.dat',
     $     form="FORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)
      open(unit=91,
     $     file='skwc.dat',
     $     form="FORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)
      open(unit=92,
     $     file='cycl.dat',
     $     form="FORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)
      open(unit=93,
     $     file='Eph.dat',
     $     form="FORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)
      open(unit=94,
     $     file='Zph.dat',
     $     form="FORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)
      open(unit=95,
     $     file='ell.dat',
     $     form="FORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)


#if 0
c     open file to dump the contours:
      open(33,file='pvocon.dat',status='unknown')
      write(33,'(a40)') caseid
      write(33,'(a40)') caseid
      write(33,*) nspday,ntimestep
      write(33,'(5(1x,f3.1))') 0.,0.,0.,0.,0.
      write(33,'(3(1x,f3.1))') 0.,0.,0.
#endif


#ifdef FMECH
      stub='equi'
      zavgfile=stub(1:lenchr(stub))
      open(unit=36,
     $     file=zavgfile(1:lenchr(zavgfile))//'.u',
     $     form="UNFORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)

      open(unit=37,
     $     file=zavgfile(1:lenchr(zavgfile))//'.v',
     $     form="UNFORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)

      open(unit=38,
     $     file=zavgfile(1:lenchr(zavgfile))//'.h',
     $     form="UNFORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)

      open(unit=39,
     $     file=zavgfile(1:lenchr(zavgfile))//'.z',
     $     form="UNFORMATTED",
     $     status="UNKNOWN",
     $     position="append",
     $     iostat=ierr)

#endif



      endif

      end
