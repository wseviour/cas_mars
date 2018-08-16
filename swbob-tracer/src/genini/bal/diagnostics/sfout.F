      subroutine sfout(phyg_z,phyg_d,phyg_h)
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

      character*60 sffile_z,sffile_d,sffile_h
      INT_TYPE  ounit,ierr
      real*4 phyg_z (plon,plat,plev)    
      real*4 phyg_d (plon,plat,plev)    
      real*4 phyg_h (plon,plat,plev)    

      INT_TYPE  i,j

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Executable statements
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ! Build file name caseid.x.iexp

      sffile_z='ic.z'
      sffile_d='ic.d'
      sffile_h='ic.h'

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
