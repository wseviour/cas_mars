      subroutine spectral_to_fourier(nltsc,sc,scm1,sfcx,sfc)
      implicit none

#include <alp.h>
#include <mderiv.h>
#include <scptrs.h>
#include <decomp.h>
#include <mpcluster.h>
#include <physical.h>
#include <io.h>

      REAL_TYPE nltsc(2*(plev*nnlt)*nvsc)             ! spectral space non linear terms
      REAL_TYPE sc(2*(plev*nsp)*nsc)                  ! time level n state variable spectral coefficients
      REAL_TYPE scm1(2*(plev*nsp)*nsc)                  ! time level n-1 state variable spectral coefficients
      REAL_TYPE sfcx(0:1,nlat2*nsyn*2*plev,mm2d,0:1,jelemd,0:P_NODE-1) ! transposed fourier coefficients of state variables
      REAL_TYPE sfc(0:1,nlat2*nsyn*2*plev,mm2d,0:1,jelemd,0:M_NODE-1) ! fourier coefficients of state variables

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

#if 0
      INT_TYPE fld_p,level
#endif

#if 0
      write(MPI_UNIT,*)"-----------------------------------------------------------------------------"
      write(MPI_UNIT,*)"INITIAL LOCAL SPECTRAL COEFFICIENTS"
      write(MPI_UNIT,*)"-----------------------------------------------------------------------------"
      write(MPI_UNIT,*)"  m   n        Abs Vor SC           Div SC            Phi SC"
      write(MPI_UNIT,*)"-----------------------------------------------------------------------------"

      do mch=0,mchmax(rank)
         do m=mbeg(mch,rank),mend(mch,rank)
c     serial output
            if(m.eq.0)then
c
               sc_p=1+2*nsp*plev*scptr(m)
               call scdebug(sc(sc_p),m)
            endif
         end do
      end do
      write(MPI_UNIT,*)"-----------------------------------------------------------------------------"
#endif

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

               ! compute u and v from vor and div

               call blduv3d(sc(sc_p),
     $                      nltsc(nlt_p),
     $                      eps(1+vscptr(m)),
     $                      m)

               do ipe=0,P_NODE-1
                  do jeg=jebeg(ipe),jeend(ipe)
                     je=jeg-jebeg(ipe)+1

                     if (m.ge.mbeg(mch,rank)+2) then
                        call genalp(P(1,0,0,jeg),
     $                              genp(1,1+vscptr(m)),
     $                              ix,
     $                              m)
                     end if

                     call synthesis(nltsc(nlt_p),
     $                              sc(sc_p),
     $                              scm1(sc_p), 
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

      end







