      subroutine testalp()
      implicit none

#include <alp.h>
#include <physical.h>
#include <gaussian.h>
#include <decomp.h>
#include <scptrs.h>
#include <mpcluster.h>
#include <io.h>

      REAL_TYPE psum(nvsc),ptmp
      REAL_TYPE err,val

      INT_TYPE mch,ioe,ix,jeg,j,n,m

      do mch=0,1
         do ioe=0,1

            ix=0          ! P memory pointer for compute on fly

            do jeg=1,jelemg
               call loadalp(Pseed(1,0,jeg,ioe,mch),P(1,0,ix,jeg),mbeg(mch,rank)+ioe)
            end do

            do m=mbeg(mch,rank)+ioe,mend(mch,rank),2


               do n=m,nn+1
                  psum(vscptr(m)+1+n-m)=0.0d0
               end do

               do jeg=1,jelemg

                  if (m.ge.mbeg(mch,rank)+2) then
                     call genalp(P(1,0,0,jeg),
     $                          genp(1,1+vscptr(m)),
     $                          ix,
     $                          m)
                  end if

                  do n=m,nn+1
                     ptmp=0.0d0
                     do j=1,nlat2
                        ptmp=ptmp+P(j,n,ix,jeg)*P(j,n,ix,jeg)*gen_gw(j,jeg)
                     end do
                     psum(vscptr(m)+1+n-m)=psum(vscptr(m)+1+n-m)+ptmp
                  end do

               end do

               ix=1-ix

            end do
         end do
      end do
#if 0
      write(MPI_UNIT,*)
      write(MPI_UNIT,*)"-----------------------------------------------------------"
      write(MPI_UNIT,*)"Maximum Deviations from Unity for Legendre Polynomial Norms"
      write(MPI_UNIT,*)"-----------------------------------------------------------"
      do mch=0,1
         do m=mbeg(mch,rank),mend(mch,rank)
            err=0.0d0
            do n=m,nn+1
               val=ABS(1.0d0-2.0*psum(vscptr(m)+1+n-m))
               if (err.lt.val) err=val
            end do
            write(MPI_UNIT,9010)m,err
 9010       format("wave no. ",i5," pnorm error= ",e22.15)
         end do
      end do
      write(MPI_UNIT,*)"-----------------------------------------------------------"
#endif
      end



