      subroutine initsht(colatitude, 
     $                   gausswt, 
     $                   plonmax, 
     $                   Alpha_coord,
     $                   Radius, 
     $                   Omega)

c---------------------------------------------------------
c Description:
c     initialize decompositions, coordinate systems
c     mdecomp, jdecomp, initscptrs, initcoord
c     initialize numlons and fftptr
c     store max longitudes in node local lonmax,
c     initialize FFT weights,
c     initialize ALP related information
c     initialize and test the Associated Legendre Polynomial seeds
c     initialize coriolis factor and annp1(0:nn+1)
c---------------------------------------------------------

      implicit none

#include <fft.h>
#include <alp.h>
#include <mderiv.h>
#include <decomp.h>
#include <scptrs.h>

#include <mpcluster.h>

      REAL_TYPE colatitude(plat/2)    ! colatitudes 
      REAL_TYPE gausswt(plat/2)       ! gaussian weights 
      INT_TYPE  plonmax(jelemg)       ! max longitude for each latitude
      REAL_TYPE Alpha_coord           ! angle between coordinate and rotational pole

      REAL_TYPE Radius                ! Planetary Radius
      REAL_TYPE Omega                 ! Angular velocity of Planet rotation

      INT_TYPE m,mch,je,jeg

c----------------------------------------------------------
c Begin executable statements
c-------------------------------------------------------------

c-------------------------------------------------------------
c     Initialize decompositions, pointers, coordinate systems
c-------------------------------------------------------------
      call mdecomp()
      call jdecomp()
      call initscptrs()
      call initcoord(colatitude,gausswt,plonmax,Radius)

c-------------------------------------------------------------
c     Initialize numlons and fftptr (pointer to right fft in case
c     plonmax is different with latitudes)
c-------------------------------------------------------------
      numlons=1
      fftptr(1)=numlons
      do jeg=jebeg(rank)+1,jeend(rank)
         je=jeg-jebeg(rank)+1
         if (plonmax(jeg).ne.plonmax(jeg-1)) then
            numlons=numlons+1
         endif
         if (numlons.le.MAX_NUM_LONS) then
            fftptr(je)=numlons
         else
            call haltmp('ABORTING: too many DIFFERENT longitudes: check MAX_NUM_LONS')
         endif
      end do

c-------------------------------------------------------------
c     Store max longitudes in node local lonmax,
c     Initialize FFT weights (tablution of trigo fcts) (call rffti),
c-------------------------------------------------------------
      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1
         lonmax(je)=plonmax(jeg)
         call rffti(lonmax(je),fftwts(1,fftptr(je)))
      end do

#ifdef DEBUG_LEV1
      call testrfft()
#endif

c-------------------------------------------------------------
c     Initialize ALP related information
c-------------------------------------------------------------
      do mch=0,1
         do m=mbeg(mch,rank),mend(mch,rank)
            call initgenalp(genp(1,1+vscptr(m)),m)
            call initmderiv(eps(1+vscptr(m)),m)
         end do
      end do

c-------------------------------------------------------------
c     Initialize and test the Associated Legendre Polynomial seeds
c-------------------------------------------------------------
      call initalpseed(colatitude)

#ifdef DEBUG_LEV1
      call testalp()
#endif

c-------------------------------------------------------------
c     Initialize coriolis factor (fcor) and Radius/(n*(n+1)): annp1(0:nn+1)
c-------------------------------------------------------------
      call initspecuv(Alpha_coord,Radius,Omega)
      
      end     




