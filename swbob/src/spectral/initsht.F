
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
c     INT_TYPE  numlons  ! number of different longitude grids on a hemisphere
c     INT_TYPE  fftptr(jelemd) ! points to appropriate fft
c     REAL_TYPE fftwts(2*plon+15,MAX_NUM_LONS)  ! fftwts (tabulation of trigo fcts) needed for model
c     common /wts/ fftwts,fftptr,lonmax,numlons
#include <alp.h>
c     src_common/spectral/alp.h:      REAL_TYPE genp(3,nvsc) ! coef of Swarztrauber's reccurence relation
c     src_common/spectral/alp.h:      common /alp/ Pseed, genp, P
#include <mderiv.h>
c      REAL_TYPE eps(nvsc)
c      common /mderiv/ eps
#include <decomp.h>
#include <scptrs.h>
c     src_common/include/scptrs.h: INT_TYPE scgptr(0:mm),scptr(0:mm),vscptr(0:mm) ! pointer for the spectral space
c     src_common/include/scptrs.h: common /scptrs/scgptr,scptr,vscptr
#include <mpcluster.h>

c #define MAX_NUM_LONS 1 in params.h

      REAL_TYPE colatitude(plat/2)    ! colatitudes 
      REAL_TYPE gausswt(plat/2)       ! gaussian weights 
      INT_TYPE  plonmax(jelemg)       ! max longitude for each latitude, read from gridfile
      REAL_TYPE Alpha_coord           ! angle between coordinate and rotational pole

      REAL_TYPE Radius                ! Planetary Radius
      REAL_TYPE Omega                 ! Angular velocity of Planet rotation

      INT_TYPE m,mch,je,jeg

c-------------------------------------------------------------
c Begin executable statements
c-------------------------------------------------------------

c-------------------------------------------------------------
c     Initialize decompositions, pointers, coordinate systems
c-------------------------------------------------------------
      call mdecomp() ! decomposition of the spectral space
      call jdecomp() ! decomposition of the physical space
      call initscptrs() ! Initializes pointers 
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
            fftptr(je)=numlons ! pointer to point to appropriate fft
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
            call initgenalp(genp(1,1+vscptr(m)),m) ! builds genp
            call initmderiv(eps(1+vscptr(m)),m) ! builds eps
         end do
      end do

c-------------------------------------------------------------
c     Initialize and test the Associated Legendre Polynomial seeds
c-------------------------------------------------------------
      call initalpseed(colatitude) ! computes Pseed

#ifdef DEBUG_LEV1
      call testalp()
#endif

c-------------------------------------------------------------
c     Initialize coriolis factor (fcor) and Radius/(n*(n+1)): annp1(0:nn+1)
c-------------------------------------------------------------
      call initspecuv(Alpha_coord,Radius,Omega)
      
      end     




