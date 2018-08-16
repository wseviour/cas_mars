      subroutine spectral_to_physical(sc,phy)
      implicit none

#include <alp.h>
#include <mderiv.h>
#include <scptrs.h>
#include <decomp.h>
#include <mpcluster.h>
#include <physical.h>
#include <io.h>

#include <fft.h>
#if defined(MPI)
#ifdef LAM_MPI
#include <mpif.h>
#endif
#endif

      REAL_TYPE sc(2*(plev*nsp)*nsc)                  ! time level n state variable spectral coefficients
      REAL_TYPE sfcx(0:1,nlat2*nsyn*2*plev,mm2d,0:1,jelemd,0:P_NODE-1) ! transposed fourier coefficients of state variables
      REAL_TYPE sfc(0:1,nlat2*nsyn*2*plev,mm2d,0:1,jelemd,0:M_NODE-1) ! fourier coefficients of state variables
      REAL_TYPE phy(nlon,nlat2*nphy*2*plev,ielemd,jelemd)               ! physical space state variables

      INT_TYPE  m          ! the fourier wave number
      INT_TYPE  mch        ! there are two chunks of m for load balancing. mch is the chunk pointer
      INT_TYPE  mloc       ! local wave numbers index
      INT_TYPE  nlt_p      ! non linear spectral coefficients pointer
      INT_TYPE  sc_p       ! points to the working position in the state variable spectral coefficients,
                           ! based on the current working value of m

      INT_TYPE  ioe        ! odd even pointer
      INT_TYPE  ipe        ! indexes the process number
      INT_TYPE  ix         ! points to the working buffer of the ALP being computed on the fly
      
      INT_TYPE  je         ! local  latitude tile index
      INT_TYPE  jeg        ! global latitude tile index


      INT_TYPE  isyn
      INT_TYPE  i,ie
      INT_TYPE  iptr
      
      REAL_TYPE fft(plon)          ! fft temporary

      do ipe=0,P_NODE-1
         call zeroreal(sfcx(0,1,1,0,1,ipe),sbuflen)
      end do

      do mch=0,mchmax(rank)
         do ioe=0,1
            ix=0          ! P memory pointer for compute on fly
            do jeg=1,jelemg
               call loadalp(Pseed(1,0,jeg,ioe,mch),P(1,0,ix,jeg),mbeg(mch,rank)+ioe)
            end do
            mloc=ioe+1
            do m=mbeg(mch,rank)+ioe,mend(mch,rank),2
               nlt_p=1+2*plev*nnlt*vscptr(m)
               sc_p =1+2*plev*nsp*scptr(m)
               do ipe=0,P_NODE-1
                  do jeg=jebeg(ipe),jeend(ipe)
                     je=jeg-jebeg(ipe)+1
                     if (m.ge.mbeg(mch,rank)+2) then
                        call genalp(P(1,0,0,jeg),
     $                              genp(1,1+vscptr(m)),
     $                              ix,
     $                              m)
                     end if
                     call synthesis1(sc(sc_p),
     $                              sfcx(0,1,mloc,mch,je,ipe),
     $                              P(1,0,ix,jeg),
     $                              m)
                  end do
               end do
               ix=1-ix
               mloc=mloc+2
            end do
         end do
      end do

      call transpb(sfcx,sfc,sbuflen,P_NODE,M_NODE)   !...transpose sfcx->sfc


      ! fft synthesis phase

      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1
         ! apply OMP directive here
         do isyn=1,nlat2*nsyn*2*plev
            ! zero truncated part of fft work array
            do m=2*(mm+1),lonmax(je)
               fft(m)=0.0d0
            end do
            ! unpack sfc array into  fft work array
            do ipe=0,M_NODE-1
               do mch=0,mchmax(ipe)
                  m=mbeg(mch,ipe) 
                  mloc=1
                  if (m.eq.0) then
                     fft(1)=sfc(0,isyn,mloc,mch,je,ipe)
                  else
                     fft(2*m  )=sfc(0,isyn,mloc,mch,je,ipe)
                     fft(2*m+1)=sfc(1,isyn,mloc,mch,je,ipe)
                  end if
                  mloc=mloc+1
                  do m=mbeg(mch,ipe)+1,mend(mch,ipe)
                     fft(2*m  )=sfc(0,isyn,mloc,mch,je,ipe)
                     fft(2*m+1)=sfc(1,isyn,mloc,mch,je,ipe)
                     mloc=mloc+1
                  end do
               end do
            end do
            call rfftb(lonmax(je),fft,fftwts(1,fftptr(je)))
            iptr=1
            do ie=1,ielem(jeg)
               do i=1,nlon
                  phy(i,isyn,ie,je)=fft(iptr)
                  iptr=iptr+1
               end do
            end do
         end do   ! end synthesis multithreaded loop
      end do   ! end loop over jeg


      end
