      subroutine physical_to_spectral(phy,sc,nlev)
      implicit none

#include <dims.h>
#include <decomp.h>
#include <mpcluster.h>
#include <fft.h>
#include <alp.h>
#include <physical.h>
#include <scptrs.h>
#include <gaussian.h>
#include <io.h>

      INT_TYPE  nlev
      REAL_TYPE phy(nlon,nlat2*nlev*2,ielemd,jelemd)
      REAL_TYPE sc(0:1,nsc*nlev)  

      REAL_TYPE fc (0:1,nlat2*plev*2,mm2d,0:1,jelemd,0:M_NODE-1)  ! should be dynamically allocated
      REAL_TYPE fcx(0:1,nlat2*plev*2,mm2d,0:1,jelemd,0:P_NODE-1)  ! should be dynamically allocated
      REAL_TYPE fft(plon)          ! fft temporary

      INT_TYPE buflen
      INT_TYPE jeg,je
      INT_TYPE i,iptr,ianal,ie
      INT_TYPE ipe,ioe,mloc,mch,m,ix

      INT_TYPE sc_p

      buflen=2*nlat2*nlev*2*mm2d*2*jelemd

      ! fft analysis phase

      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1

         ! ----------------------------
         ! apply OMP directive here
         !-----------------------------

         do ianal=1,2*nlev*nlat2

            iptr=1
            do ie=1,ielem(jeg)
               do i=1,nlon
                  fft(iptr)=phy(i,ianal,ie,je)
                  iptr=iptr+1
               end do
            end do

            call rfftf(lonmax(je),fft,fftwts(1,fftptr(je)))

            ! pack fft'd array into fc array
 
            do ipe=0,M_NODE-1
               do mch=0,1            

                  mloc=1
                  m=mbeg(mch,ipe) 

                  if (m.eq.0) then
                     fc(0,ianal,mloc,mch,je,ipe)=fft(1)
                     fc(1,ianal,mloc,mch,je,ipe)=0.0d0
                  else
                     fc(0,ianal,mloc,mch,je,ipe)=fft(2*m  )
                     fc(1,ianal,mloc,mch,je,ipe)=fft(2*m+1)
                  end if
                  mloc=mloc+1

                  do m=mbeg(mch,ipe)+1,mend(mch,ipe)
                     fc(0,ianal,mloc,mch,je,ipe)=fft(2*m)
                     fc(1,ianal,mloc,mch,je,ipe)=fft(2*m+1)
                     mloc=mloc+1
                  end do

               end do
            end do   

         end do ! end multithreaded loop over ianal

      end do   ! end loop over jeg (latitude pair tile coord)

      
#if 0
      write(MPI_UNIT,*)"initial fourier coefficients"
      write(MPI_UNIT,*)"-----------------------------------------------------------------------------"
      write(MPI_UNIT,*)"  lat       m          FC"
      write(MPI_UNIT,*)"-----------------------------------------------------------------------------"
      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1
         call fcdebug(fc,lat(1,jeg),je,nlev)
      end do
      write(MPI_UNIT,*)"-----------------------------------------------------------------------------"
#endif

      call transpf(fc,fcx,buflen,M_NODE,P_NODE)   !...transpose fc->fcx

#if 0
      write(MPI_UNIT,*)"transposed fourier coefficients"
      write(MPI_UNIT,*)"-----------------------------------------------------------------------------"
      write(MPI_UNIT,*)"  lat       m          FC"
      write(MPI_UNIT,*)"-----------------------------------------------------------------------------"
      do ipe=0,P_NODE-1
         do jeg=jebeg(ipe),jeend(ipe)
            je=jeg-jebeg(ipe)+1
            call fcxdebug(fcx(0,1,1,0,je,ipe),lat(1,jeg),nlev)
         end do
      end do
      write(MPI_UNIT,*)"-----------------------------------------------------------------------------"
#endif

      call zeroreal(sc,2*nlev*nsc)

      do mch=0,1
         do ioe=0,1

            ix=0          ! P memory pointer for compute on fly

            do jeg=1,jelemg
               call loadalp(Pseed(1,0,jeg,ioe,mch),P(1,0,ix,jeg),mbeg(mch,rank)+ioe)
            end do

            mloc=ioe+1

            do m=mbeg(mch,rank)+ioe,mend(mch,rank),2

               sc_p=1+nlev*scptr(m)

               do ipe=0,P_NODE-1
                  do jeg=jebeg(ipe),jeend(ipe)
                     je=jeg-jebeg(ipe)+1

                     if (m.ge.mbeg(mch,rank)+2) then
                         call genalp(P(1,0,0,jeg),
     $                               genp(1,1+vscptr(m)),
     $                               ix,
     $                               m)
                     end if

                     call gen_analysis(fcx(0,1,mloc,mch,je,ipe),
     $                                 sc(0,sc_p),
     $                                 P(1,0,ix,jeg),
     $                                 gen_gw(1,jeg),
     $                                 lat(1,jeg),
     $                                 m,
     $                                 nn,
     $                                 nlev)

                  end do
               end do

               mloc=mloc+2
               ix=1-ix

            end do
         end do
      end do


      end




