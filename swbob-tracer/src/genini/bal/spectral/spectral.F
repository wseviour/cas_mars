      subroutine spectral(afcx,sfcx)
      implicit none

#include <alp.h>
#include <mderiv.h>
#include <spectral.h>
#include <gaussian.h>
#include <scptrs.h>
#include <decomp.h>
#include <mpcluster.h>
c     Added by leo for output of sc
#include <io.h>
#include <dims.h>
#include <physical.h>

c----------------------------------------
c     Local Declarations
c----------------------------------------

      REAL_TYPE afcx(0:1,nlat2*nanal*2*plev,mm2d,0:1,jelemd,0:P_NODE-1) ! fourier coefficients of non linear terms
      REAL_TYPE sfcx(0:1,nlat2*nsyn*2*plev,mm2d,0:1,jelemd,0:P_NODE-1) ! fourier coefficients of non linear terms

      REAL_TYPE revorsc(0:1,nsc*plev) ! spec coef for restart
      REAL_TYPE redivsc(0:1,nsc*plev)
      REAL_TYPE rephisc(0:1,nsc*plev)

      INT_TYPE  m          ! the fourier wave number
      INT_TYPE  mch        ! there are two chunks of m for load balancing. mch is the chunk pointer
      INT_TYPE  mloc       ! wave numbers chunk based index
      INT_TYPE  nlt_p      ! points to the working position in the NLT spectral coefficients
                           ! based on the current working value of m
      INT_TYPE  sc_p       ! points to the working position in the state variable spectral coefficients,
                           ! based on the current working value of m
      INT_TYPE  for_p
      INT_TYPE  ioe        ! odd even pointer
      INT_TYPE  ipe        ! indexes the process number
      INT_TYPE  ix         ! points to the working buffer of the ALP being computed on the fly
      
      INT_TYPE  je         ! local latitude tile index
      INT_TYPE  jeg        ! global latitude tile index
#if 1
      INT_TYPE  k, recnum
      ! Used to output forurier coeffs to file 44
      INT_TYPE fld_p,level
#endif
      INT_TYPE  reounit    ! unit for ouput file restart
      character*80 refilename
      REAL_TYPE plon_file,plev_file,time_start
      ! Returns index of next available fortran unit
      INT_TYPE navu
      external navu


c----------------------------------------
c     Executable Statements
c----------------------------------------

      call zeroreal(nltsc,2*(plev*nnlt)*nvsc)
      call zeroreal(sfcx,sbuflen*P_NODE)

c      write(*,*)'---------------------------'
c      write(*,*)'Begin loop in m in spectral'
c      write(*,*)'---------------------------'

      do mch=0,mchmax(rank)
         do ioe=0,1

            mloc=ioe+1
            ix=0          ! P memory pointer for compute on fly

            do jeg=1,jelemg
               call loadalp(Pseed(1,0,jeg,ioe,mch),P(1,0,ix,jeg),mbeg(mch,rank)+ioe)
            end do

            do m=mbeg(mch,rank)+ioe,mend(mch,rank),2

               nlt_p=1+2*plev*nnlt*vscptr(m)
               sc_p =1+2*plev*nsp*scptr(m)
               for_p =1+2*plev*scptr(m)

               do ipe=0,P_NODE-1
                  do jeg=jebeg(ipe),jeend(ipe)
                     je=jeg-jebeg(ipe)+1

                     if (m.ge.mbeg(mch,rank)+2) then
                         call genalp(P(1,0,0,jeg),
     $                               genp(1,1+vscptr(m)),
     $                               ix,
     $                               m)
                     end if

                     call analysis(afcx(0,1,mloc,mch,je,ipe),
     $                             nltsc(nlt_p),
     $                             P(1,0,ix,jeg),
     $                             gw(1,1,jeg),
     $                             m)
                  end do
               end do

      ! compute rhs of eqns of motion and advance state variables
      ! one timestep, returns u and v coefficients in nltsc

               call advance(nltsc(nlt_p),
     $                      sc(sc_p),
     $                      eps(1+vscptr(m)),
     $                      m)

               do ipe=0,P_NODE-1
                  do jeg=jebeg(ipe),jeend(ipe)
                     je=jeg-jebeg(ipe)+1
                     
                     ! go to fourier space: sfcx
                     call synthesis(nltsc(nlt_p),
     $                              sc(sc_p),
     $                              sfcx(0,1,mloc,mch,je,ipe),
     $                              P(1,0,ix,jeg),
     $                              m)
                  end do
               end do

               mloc=mloc+2
               ix=1-ix

            end do ! End loop on m
         end do ! End loop on ioe
      end do ! End loop on mch
cccccccccccc     Save for initial condition  ccccccccccccccccccccccccc
#if 1
      if (rank.eq.root) then
         refilename='progsc.ini'
         reounit=navu()
         open(UNIT   = reounit,
     $        FILE   = refilename,
     $        ACCESS = 'sequential',
     $        STATUS = 'unknown',
     $        FORM   = 'unformatted')
            
         plon_file=PLON
         plev_file=PLEV
         time_start=0.0d0
         write(reounit)plon_file,plev_file,time_start
      end if

      do mch=0,mchmax(rank)
         do ioe=0,1
            do m=mbeg(mch,rank)+ioe,mend(mch,rank),2
               sc_p =1+2*plev*nsp*scptr(m)
               call transc(sc(sc_p),m,revorsc,redivsc,rephisc)
            end do              ! End loop on m
         end do                 ! End loop on ioe
      end do                    ! End loop on mch
      call out_resc(revorsc,redivsc,rephisc,reounit)
      call out_resc(revorsc,redivsc,rephisc,reounit)
      
      if (rank.eq.root) then
         close(reounit)
      end if

#endif

      end
