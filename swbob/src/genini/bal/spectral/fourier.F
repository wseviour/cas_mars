      subroutine fourier(sfc,phy,nlt,afc)
      implicit none

!--------------------------------------------------------------------
#include <physical.h>
#include <decomp.h>
#include <fft.h>
#include <mpcluster.h>
#include <io.h>
!--------------------------------------------------------------------
#if defined(MPI)
#ifdef LAM_MPI
#include <mpif.h>
#endif
#endif
#include <planet.h>
#include <fordamp.h>

      REAL_TYPE sfc(0:1 ,nlat2*nsyn*2*plev,mm2d,0:1,jelemd,0:M_NODE-1)    ! fourier coefficients for synthesis
      REAL_TYPE phy(nlon,nlat2*nphy*2*plev,ielemd,jelemd)               ! physical space state variables
      REAL_TYPE nlt(nlon,nlat2*nnlt*2*plev,ielemd,jelemd)               ! physical space non linear terms
      REAL_TYPE afc(0:1 ,nlat2*nanal*2*plev,mm2d,0:1,jelemd,0:M_NODE-1)   ! fourier coefficients of non linear terms

      INT_TYPE  isyn,ianal
      INT_TYPE  i,ie,je,jeg

      INT_TYPE  mch,ipe
      INT_TYPE  m,mloc
      INT_TYPE  iptr
      
      REAL_TYPE fft(plon)          ! fft temporary
      REAL_TYPE cfl,cflt,umax,umaxt

#if 0
      INT_TYPE fld_p,level
       
c     fld_p: 1=vor,2=div,3=theta/phi,4=U,5=V
      fld_p=5
      level=1
      write(*,*)"START OF FOURIER:"
      write(*,*)phy(1,1,1,1)
#endif

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

      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1

         ! ----------------------------
         ! apply OMP directive here
         !-----------------------------

         do ie=1,ielem(jeg)
            call bldnlt(phy(1,1,ie,je),
     $                  nlt(1,1,ie,je),
     $                  phitopo(1,1,0,1,ie,je),
     $                  lat(1,jeg),
     $                  lon(1,ie,je))

         end do
      end do

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
