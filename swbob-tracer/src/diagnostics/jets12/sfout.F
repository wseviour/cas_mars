      subroutine sfout(phyg_z,phyg_d,phyg_h,icountp)
      implicit none

#include <io.h>                 ! caseid
#include <dims.h>               ! plon,plat
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Local declarations
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      INT_TYPE  navu
      external  navu

      INT_TYPE  lenchr
      external  lenchr

      character*60 stub_z,sffile_z
      character*60 stub_d,sffile_d
      character*60 stub_h,sffile_h
      character*5  tail,tail_lj
      INT_TYPE  ounit,ierr,icountp
      real*4 phyg_z (plon,plat,plev)    
      real*4 phyg_d (plon,plat,plev)    
      real*4 phyg_h (plon,plat,plev)    

      INT_TYPE  i,j

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Executable statements
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ! Build file name caseid.x.iexp

c      stub_z=caseid(1:lenchr(caseid))//'.z.'
      stub_z='z.'
      stub_d='d.'
      stub_h='h.'
      write(tail,1)icountp
 1    format(i5)
      call ljust0(tail,tail_lj)
      sffile_z=stub_z(1:lenchr(stub_z))//tail_lj(1:lenchr(tail_lj))
      sffile_d=stub_d(1:lenchr(stub_d))//tail_lj(1:lenchr(tail_lj))
      sffile_h=stub_h(1:lenchr(stub_h))//tail_lj(1:lenchr(tail_lj))

      ounit=navu()
      open(unit=ounit,
     $     file=sffile_z,
     $     form="UNFORMATTED",
     $     status="UNKNOWN",
     $     iostat=ierr)
      call writespreal(phyg_z,plon*plat*plev,ounit)
      close(ounit)

      ounit=navu()
      open(unit=ounit,
     $     file=sffile_d,
     $     form="UNFORMATTED",
     $     status="UNKNOWN",
     $     iostat=ierr)
      call writespreal(phyg_d,plon*plat*plev,ounit)
      close(ounit)

      ounit=navu()
      open(unit=ounit,
     $     file=sffile_h,
     $     form="UNFORMATTED",
     $     status="UNKNOWN",
     $     iostat=ierr)
      call writespreal(phyg_h,plon*plat*plev,ounit)
      close(ounit)

      end
