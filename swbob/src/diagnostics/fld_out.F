      subroutine fld_out(phyo,lfld,icountp)
      implicit none

#include <io.h>
#include <dims.h>
#include <mpcluster.h>
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Local declarations
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      INT_TYPE  navu
      external  navu

      INT_TYPE  lenchr
      external  lenchr

      real*4 phyo(nlon,nlat2,0:1,plev,ielemd,jelemd)

      character*60 ofile
      character*5  tail
      character*1  lfld

      INT_TYPE  ounit,ierr,icountp
      INT_TYPE  i,j,ie,je,k,ins,ipe

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Executable statements
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      write(tail,1)icountp
 1    format(i5.5)
      ofile=lfld(1:lenchr(lfld))//"."//tail

      do ins=0,1

      do ipe=0,P_NODE-1
         call synch()
         if (rank.eq.ipe) then

            ounit=navu()
            open(unit=ounit,file=ofile,form='unformatted',
     .           status='unknown',access='stream',position='append')
            if ((rank.eq.0).and.(ins.eq.0)) rewind(ounit)

            if (ins.eq.0) then
               write(ounit) (((((phyo(i,j,0,k,ie,je)
     .              ,i=1,nlon),ie=1,ielemd)
     .              ,j=1,nlat2),je=1,jelemd),k=1,plev)
            endif
            if (ins.eq.1) then
               write(ounit) (((((phyo(i,j,1,k,ie,je)
     .              ,i=1,nlon),ie=1,ielemd)
     .              ,j=nlat2,1,-1),je=jelemd,1,-1),k=1,plev)
            endif

            close(ounit)

         endif

         call synch()
      enddo

      enddo

      end
