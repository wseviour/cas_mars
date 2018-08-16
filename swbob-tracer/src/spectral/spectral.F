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
#include <modeltime.h>
#include <dims.h>
#include <timestep.h>
#include <timeint.h>
#include <physical.h>
#include <stptrs.h>
#include <uspec.h>

c----------------------------------------
c     Local Declarations
c----------------------------------------

      REAL_TYPE afcx(0:1,nlat2*nanal*2*plev,mm2d,0:1,jelemd,0:P_NODE-1) ! fourier coefficients of non linear terms
      REAL_TYPE sfcx(0:1,nlat2*nsyn*2*plev,mm2d,0:1,jelemd,0:P_NODE-1) ! fourier coefficients of non linear terms

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

      REAL_TYPE revorsc(0:1,nsc*plev) ! spec coef for restart
      REAL_TYPE redivsc(0:1,nsc*plev) ! spec coef for restart
      REAL_TYPE rephisc(0:1,nsc*plev) ! spec coef for restart
      INT_TYPE  reounit    ! unit for ouput file restart
      character*80 refilename
      REAL_TYPE plon_file,plev_file,fld_file

      REAL_TYPE Erott(0:nn,4),Edivt(0:nn,4),Etott(0:nn,4),Ztott(0:nn,4)

      ! Returns index of next available fortran unit
      INT_TYPE navu
      external navu
      INT_TYPE  lenchr
      external  lenchr
      INT_TYPE n,i

      INT_TYPE iotoggle
      character*1 iostr
      common /rstrt0/ iotoggle

c----------------------------------------
c     Executable Statements
c----------------------------------------

      call zeroreal(nltsc,2*(plev*nnlt)*nvsc)
      call zeroreal(sfcx,sbuflen*P_NODE)
      call zeroreal(Erot,(nn+1)*4)
      call zeroreal(Ediv,(nn+1)*4)
      call zeroreal(Etot,(nn+1)*4)
      call zeroreal(Ztot,(nn+1)*4)

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
                  enddo
               enddo

      ! compute rhs of eqns of motion and advance state variables
      ! one timestep, returns u and v coefficients in nltsc

               call advance(nltsc(nlt_p),
     $                      sc(sc_p),
     $                      scm1(sc_p),
     $                      fthesc(for_p),
     $                      ftoposc(for_p),
     $                      ftopoktsc(for_p),
     $                      fvorsc(for_p),
     $                      fthpsc(for_p),
     $                      eps(1+vscptr(m)),
     $                      m)

               do ipe=0,P_NODE-1
                  do jeg=jebeg(ipe),jeend(ipe)
                     je=jeg-jebeg(ipe)+1
                     
                     ! go to fourier space: sfcx
                     call synthesis(nltsc(nlt_p),
     $                              sc(sc_p),
     $                              scm1(sc_p),
     $                              sfcx(0,1,mloc,mch,je,ipe),
     $                              P(1,0,ix,jeg),
     $                              m)
                  enddo
               enddo

               mloc=mloc+2
               ix=1-ix

            enddo ! End loop on m
         enddo ! End loop on ioe
      enddo ! End loop on mch

c     update random phases of spectral forcing:
      call forsc(fvorsc,sc)

c     sum KE spectrum (from blduv) across nodes:
      call psumreal(Erot,Erott,(nn+1)*4)
      call broadcastreal(Erott,(nn+1)*4)
      call copyreal(Erott,Erot,(nn+1)*4)
      call psumreal(Ediv,Edivt,(nn+1)*4)
      call broadcastreal(Edivt,(nn+1)*4)
      call copyreal(Edivt,Ediv,(nn+1)*4)
      call psumreal(Ztot,Ztott,(nn+1)*4)
      call broadcastreal(Ztott,(nn+1)*4)
      call copyreal(Ztott,Ztot,(nn+1)*4)
      do i=1,4
         do n=0,nn
            Etot(n,i)=Erot(n,i)+Ediv(n,i)
         enddo
      enddo

cccccccccccc     Save for restart  ccccccccccccccccccccccccc
c      if (itime.eq.ntimestep+1)then

c      if (time_current.gt.time_stop) then
      if ((time_current.gt.time_stop).or.
     .    (mod(time_current,time_zon).lt.(dt2/2.d0))) then 

         iotoggle=mod(iotoggle+1,2)
         write(iostr,'(i1)'),iotoggle

         if (rank.eq.root) then
c            refilename=caseid(1:lenchr(caseid))//'resta.ini'
c            refilename='resta.ini'
            refilename='resta.ini'//iostr
            reounit=navu()
            open(UNIT   = reounit,
     $           FILE   = refilename,
     $           ACCESS = 'sequential',
     $           STATUS = 'unknown',
     $           FORM   = 'unformatted')
            rewind reounit
            plon_file=PLON
            plev_file=PLEV
            time_start=time_current
            write(reounit)plon_file,plev_file,time_start
            write(reounit)timestep
         end if

c     SC coeffs:
         do mch=0,mchmax(rank)
            do ioe=0,1
               do m=mbeg(mch,rank)+ioe,mend(mch,rank),2
                  sc_p=1+2*plev*nsp*scptr(m)
                  call transc(sc(sc_p),m,revorsc,redivsc,rephisc)
               enddo
            enddo
         enddo
         fld_file=vor_p
         call out_resc(revorsc,reounit,fld_file)
         fld_file=div_p
         call out_resc(redivsc,reounit,fld_file)
         fld_file=phi_p
         call out_resc(rephisc,reounit,fld_file)

c     SCM1 coeffs:
         do mch=0,mchmax(rank)
            do ioe=0,1
               do m=mbeg(mch,rank)+ioe,mend(mch,rank),2
                  sc_p=1+2*plev*nsp*scptr(m)
                  call transc(scm1(sc_p),m,revorsc,redivsc,rephisc)
               enddo
            enddo
         enddo
         fld_file=vor_p
         call out_resc(revorsc,reounit,fld_file)
         fld_file=div_p
         call out_resc(redivsc,reounit,fld_file)
         fld_file=phi_p
         call out_resc(rephisc,reounit,fld_file)

c     forcing coeffs:
         fld_file=nsp+1
         do mch=0,mchmax(rank)
            do ioe=0,1
               do m=mbeg(mch,rank)+ioe,mend(mch,rank),2
                  sc_p=1+2*plev*scptr(m)
                  call transc1(fvorsc(sc_p),m,revorsc)
               enddo
            enddo
         enddo
         call out_resc(revorsc,reounit,fld_file)

         if (rank.eq.root) then
            close(reounit)
         end if

      endif

      end
