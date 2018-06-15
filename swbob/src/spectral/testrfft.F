      subroutine testrfft()
      implicit none

#include <dims.h>
#include <fft.h>
#include <decomp.h>
#include <mpcluster.h>
#include <io.h>

      REAL_TYPE ffterr,ffterr_max
      REAL_TYPE a(plon),b(plon)
      INT_TYPE  i,je,jeg

      ! Triangle wave test

#if 0
      write(MPI_UNIT,*)
      write(MPI_UNIT,9010)
 9010 format("--------------rfft test------------------")
#endif

      ffterr_max=-911.0       
      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1
         do i=1,lonmax(je)
            a(i)=i
            b(i)=i
         end do

         call rfftf(lonmax(je),a,fftwts(1,fftptr(je)))
         call rfftb(lonmax(je),a,fftwts(1,fftptr(je)))

         ffterr=0.0d0
         do i=1,lonmax(je)
            ffterr=ffterr+ABS((a(i)-b(i))**2)         
         end do
         ffterr=DSQRT(ffterr)

         ffterr_max=MAX(ffterr_max,ffterr)
      end do
#if 1
      write(MPI_UNIT,9020)ffterr_max
 9020 format(" MAX rfft fwd/back rms error = ",e12.6)
      write(MPI_UNIT,9030)
 9030 format("-----------------------------------------")
#endif

      end       






