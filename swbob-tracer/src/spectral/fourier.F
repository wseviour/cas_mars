      subroutine fourier(sfc,phy,nlt,afc)
      implicit none

!--------------------------------------------------------------------
#include <physical.h>
#include <decomp.h>
#include <fft.h>
#include <mpcluster.h>
#include <io.h>
#include <fordamp.h>
#if defined(MPI)
#ifdef LAM_MPI
#include <mpif.h>
#endif
#endif
#include <planet.h>
#include <timestep.h>
#include <modeltime.h>

      REAL_TYPE sfc(0:1 ,nlat2*nsyn*2*plev,mm2d,0:1,jelemd,0:M_NODE-1)    ! fourier coefficients for synthesis
      REAL_TYPE phy(nlon,nlat2*nphy*2*plev,ielemd,jelemd)               ! physical space state variables
      REAL_TYPE nlt(nlon,nlat2*nnlt*2*plev,ielemd,jelemd)               ! physical space non linear terms
      REAL_TYPE afc(0:1 ,nlat2*nanal*2*plev,mm2d,0:1,jelemd,0:M_NODE-1)   ! fourier coefficients of non linear terms

      INT_TYPE  isyn,ianal
      INT_TYPE  i,ie,je,jeg

      INT_TYPE  mch,ipe
      INT_TYPE  m,mloc
      INT_TYPE  iptr
      INT_TYPE  nacc
      
      REAL_TYPE fft(plon)          ! fft temporary
      REAL_TYPE cfl,cflt,cflacc,umax,umaxt,rlonmax
      common /cfladapt/ cflacc,nacc

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

      !-----------------------------------------
      ! Build Non Linear Terms in Physical Space
      ! ----------------------------------------

      umax=0.0d0
      call zeroreal(fum1n  ,jelemd*nlat2*2*plev)

      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1
         do ie=1,ielem(jeg)
            call bldnlt(phy(1,1,ie,je),
     $                  nlt(1,1,ie,je),
     $                  kvcoef(1,0,1,je),
     $                  ktcoef(1,0,1,je),
     $                  fthr(1,1,0,1,ie,je),
     $                  fmechu(1,1,0,1,ie,je),
     $                  fmechv(1,1,0,1,ie,je),
c     $                  fum1(1,1,0,1,ie,je),
     $                  fum1n(1,0,1,je),
     $                  fum1t(1,0,1,je),
     $                  lat(1,jeg),
     $                  lon(1,ie,je),
     $                  umax)
         end do
      end do

      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1
         rlonmax=1.0d0/lonmax(je)
         call scal(nlat2*2*plev, rlonmax, fum1n(1,0,1,je), 1)
      enddo
      call copyreal(fum1n,fum1t,nlat2*2*plev*jelemd)


      cfl=timestep*umax*nn/Radius      

      call pmaxreal(cfl,cflt,1)
      call broadcastreal(cflt,1)
      cfl=cflt
#ifdef ADAPT
      reducedt=.false.
      cflacc=cflacc+cfl
      nacc=nacc+1

      if (time_current.eq.time_start) then
         cflacc=0.d0
         nacc=0
      endif
      if(cfl.ge.0.6d0) reducedt=.true.
      if ((mod(time_current,time_zon).lt.(dt2/2.d0)).and.(itime.gt.0)) then
#else
      reducedt=.false.
      if ( (mod(itime-1,zonfreq).eq.0).and.(itime.gt.0) )then
#endif
         if (rank.eq.root) then
            write(*,*)
            write(*,11)'t,cfl,dt:',time_current/SecPerDay,cfl,timestep/SecPerDay,'   '
         endif
         cflacc=0.d0
         nacc=0
      endif
c
      if(cfl.ge.1.5d0)then
         if (rank.eq.root) then
            write(*,*)
            write(*,*)'itime,cfl:',itime,cfl
         endif
         call haltmp('Stopped, cfl too big, greater than 1.5')
      end if
 10   format(a14,1x,f11.4,3x,f9.6,a3,$)
 11   format(a20,1x,f11.4,3x,f6.4,3x,f10.8,a3,$)

      ! fft analysis phase

      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1

         ! ----------------------------
         ! apply OMP directive here
         !-----------------------------

         do ianal=1,nlat2*nnlt*2*plev

            iptr=1
            do ie=1,ielem(jeg)
               do i=1,nlon
                  fft(iptr)=nlt(i,ianal,ie,je)
                  iptr=iptr+1
               end do
            end do

            call rfftf(lonmax(je),fft,fftwts(1,fftptr(je)))

            ! pack fft'd array into afc array
 
            do ipe=0,M_NODE-1
               do mch=0,mchmax(ipe)

                  m=mbeg(mch,ipe) 
                  mloc=1

                  if (m.eq.0) then
                     afc(0,ianal,mloc,mch,je,ipe)=fft(1)
                     afc(1,ianal,mloc,mch,je,ipe)=0.0d0
                  else
                     afc(0,ianal,mloc,mch,je,ipe)=fft(2*m  )
                     afc(1,ianal,mloc,mch,je,ipe)=fft(2*m+1)
                  end if
                  mloc=mloc+1

                  do m=mbeg(mch,ipe)+1,mend(mch,ipe)
                     afc(0,ianal,mloc,mch,je,ipe)=fft(2*m)
                     afc(1,ianal,mloc,mch,je,ipe)=fft(2*m+1)
                     mloc=mloc+1
                  end do

               end do
            end do

         end do   ! end multithreaded loop over ianal

      end do   ! end loop over jeg (latitude pair tile coord)

      end 


   
















