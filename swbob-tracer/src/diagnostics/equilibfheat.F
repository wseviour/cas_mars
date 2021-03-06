      subroutine equilib
      implicit none

#include <diag.h>
#include <decomp.h>
#include <fields.h>
#include <physical.h>
#include <mpcluster.h>
#include <fft.h>
#include <modeltime.h>
#include <timestep.h>

#include <equilib.h>

      REAL_TYPE ubarnode(nlat2*2*plev,jelemd)
      REAL_TYPE vbarnode(nlat2*2*plev,jelemd)
      REAL_TYPE hbarnode(nlat2*2*plev,jelemd)
      REAL_TYPE zbarnode(nlat2*2*plev,jelemd)
      INT_TYPE  je,jeg,ie,k
      REAL_TYPE rlonmax
      REAL_TYPE accinv
      INT_TYPE icount,nacc

c     for sweeping out a series of forcing amplitudes:
#include <fordamp0.h>
#include <fordamp.h>
      REAL_TYPE udiff0,udiff00,udifft,updown
      common /fadjust/ updown
      REAL_TYPE dimfaca,dimfacu,forinc,ampmax,utol,u0
      parameter(dimfacu=86165./6371000.)
      parameter(dimfaca=86165.*86165./6371000.)
      parameter(forinc=2.0d-8*dimfaca,ampmax=1.2d-6*dimfaca,utol=0.001*dimfacu)

      amp=amph

c     zero out nodal zonal average arrays:
      call zeroreal(ubarnode  ,jelemd*nlat2*2*plev)
      call zeroreal(vbarnode  ,jelemd*nlat2*2*plev)
      call zeroreal(hbarnode  ,jelemd*nlat2*2*plev)
      call zeroreal(zbarnode  ,jelemd*nlat2*2*plev)

c     zonally average:
      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1
         do ie=1,ielem(jeg)
            call equibld(phy(1,1,1,ie,je),
     $                   ubarnode(1,je),vbarnode(1,je),
     $                   hbarnode(1,je),zbarnode(1,je),
     $                   cslat(1,jeg),snlat(1,jeg))
         enddo
      enddo
c     scale by 1/lonmax:
      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1
         rlonmax=1.0d0/lonmax(je)
         call scal(nlat2*2*plev, rlonmax, ubarnode(1,je), 1)
         call scal(nlat2*2*plev, rlonmax, hbarnode(1,je), 1)
         call scal(nlat2*2*plev, rlonmax, vbarnode(1,je), 1)
         call scal(nlat2*2*plev, rlonmax, zbarnode(1,je), 1)
      enddo

c     compute diff ubar-ubar0:
      udiff0=0.0
      udiff00=0.0
      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1
         do k=1,nlat2*2*plev
            udiff0=dmax1(udiff0,dabs(ubarnode(k,je)-ubarnode0(k,je)))
            udiff00=dmax1(udiff00,dabs(ubarnode(k,je)-ubarnode00(k,je)))
         enddo
      enddo
      call pmaxreal(udiff0,udifft,1)
      call broadcastreal(udifft,1)
      udiff0=udifft
      call pmaxreal(udiff00,udifft,1)
      call broadcastreal(udifft,1)
      udiff00=udifft

      call copyreal(ubarnode,ubarnode0,nlat2*2*plev*jelemd)

      if (mod(time_current,time_zon).lt.(dt2/2.d0)) then
         if (rank.eq.root) write(*,*) ' udiff0, udiff00 = ',udiff0/dimfacu,udiff00/dimfacu
      endif

c     if at equilibrium then update forcing and output:
      if ((udiff0.le.(utol*timestep)).and.(udiff00.ge.(100.0*utol))) then

c     extract ubar at equator:
         if (rank.eq.(P_NODE-1)) then 
            u0=ubarnode(nlat2,jeend(rank)-jebeg(rank)+1)
            call sendreal(u0,1,root,1)
         endif
         if (rank.eq.root) then
            call recvreal(u0,1,P_NODE-1,1)
         endif
         call broadcastreal(u0,1)

         if (rank.eq.root) then
            write(6,*)
            write(6,'(a,f11.4,a,f6.3,a,f8.3)') 
     .           ' equilibrium reached:  t = ',time_current,
     .           ';  famp = ',amp*1.0d6/dimfaca,
     .           ';  u0 = ',u0/dimfacu
            if (amp.eq.0.d0) open(55,file='equil.dat',form='formatted',position='append')
            write(55,'(f11.4,2x,f6.3,2x,f8.3)') time_current,amp*1.0d6/dimfaca,u0/dimfacu
            call flush(55)
         endif
         
c     increment forcing amplitude (increasing/decreasing):
         if (amp.eq.0.d0) updown=1.d0
         if (amp.ge.ampmax) updown=-1.d0
         amp=amp+forinc*updown

c     check for end of sweep:
         if ((amp.le.0.d0).and.(updown.eq.-1.d0)) then 
            if (rank.eq.root) close(55)
            time_current=time_stop
         endif

         do jeg=jebeg(rank),jeend(rank)
            je=jeg-jebeg(rank)+1
            do ie=1,ielem(jeg)
               do k=1,plev
                  call initmomforcing(fmechu(1,1,0,k,ie,je),fmechv(1,1,0,k,ie,je),
     $                 cslat(1,jeg),
     $                 snlat(1,jeg),
     $                 lat(1,jeg),
     $                 lon(1,ie,je),
     $                 k)
               enddo
            enddo
         enddo

         call copyreal(ubarnode,ubarnode00,nlat2*2*plev*jelemd)
         call zavg_pgather(ubarnode, ubar, plev)
         call zavg_pgather(vbarnode, vbar, plev)
         call zavg_pgather(hbarnode, hbar, plev)
         call zavg_pgather(zbarnode, zbar, plev)

         if (rank.eq.root) then
            call equiout
         endif

      endif




      end


