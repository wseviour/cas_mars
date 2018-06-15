      subroutine sf_diag(icountp,skwq,skwc,cycl,Eph,Zph)
      implicit none

#include <fields.h>          ! phy
#include <physical.h>        ! cslat(nlat2,jelemg),ielem(jelemg)
#include <decomp.h>          ! jebeg(0:P_NODE-1),jeend(0:P_NODE-1)
#include <mpcluster.h>       ! rank
#include <spectral.h>        ! sc(2*(plev*nsp)*nsc)
#include <modeltime.h>
#include <timestep.h>

      REAL_TYPE scvor(2*(plev*nsp)*nsc)
      REAL_TYPE phvor(nlon,nlat2,2*nphy*plev,ielemd,jelemd)
      real*4 physf_node_z(nlon*nlat2*2*plev,ielemd,jelemd)    
      real*4 physf_node_d(nlon*nlat2*2*plev,ielemd,jelemd)    
      real*4 physf_node_h(nlon*nlat2*2*plev,ielemd,jelemd)    
      real*4 phyg_z(plon,plat,plev)
      real*4 phyg_d(plon,plat,plev)
      real*4 phyg_h(plon,plat,plev)
      REAL_TYPE momq(3,3),momc(3,3),momt(3,3)
      REAL_TYPE skwq(3),skwc(3),cycl(3),Eph,Zph,Epht,Zpht
      INT_TYPE  je,jeg,ie,icountp,j,i,im

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Executable statements
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#if 1
c     vorticity in spectral sectors for skewness:
      call copyreal(scm1,scvor,2*(plev*nsp)*nsc)
      call zerovor(scvor)
      call spectral_to_physical(scvor,phvor)
      do im=1,3
      do i=1,3
         momq(im,i)=0.0
         momc(im,i)=0.0
      enddo
      enddo
      Eph=0.0
      Zph=0.0
#endif


      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1
         do ie=1,ielem(jeg)
            call bldsf(phy(1,1,1,ie,je),phvor(1,1,1,ie,je),
     $           physf_node_z(1,ie,je),
     $           physf_node_d(1,ie,je),
     $           physf_node_h(1,ie,je),momq,momc,Eph,Zph,
     $                   cslat(1,jeg),
     $                   snlat(1,jeg),
     $                   da(1,jeg))
         end do
      end do

      call psumreal(momq,momt,9)
      call broadcastreal(momt,9)
      call copyreal(momt,momq,9)
      call psumreal(momc,momt,9)
      call broadcastreal(momt,9)
      call copyreal(momt,momc,9)
      call psumreal(Eph,Epht,1)
      call broadcastreal(Epht,1)
      call copyreal(Epht,Eph,1)
      call psumreal(Zph,Zpht,1)
      call broadcastreal(Zpht,1)
      call copyreal(Zpht,Zph,1)

c      if (rank.eq.root) then 
c         write(*,*) momc(1,1),momc(2,1),momc(1,1)**2
c      endif

      Eph=0.5d0*Eph
      Zph=0.5d0*Zph
      do i=1,3
         skwq(i)=momq(3,i)/sqrt(momq(2,i)**3)
         skwc(i)=(momc(3,i)-3.0*momc(2,i)*momc(1,i)+2.0*momc(1,i)**3)
     .           /sqrt((momc(2,i)-momc(1,i)**2)**3)
         cycl(i)=momc(1,i)
      enddo
 
c     gather fields together on the root node:
      call sf_pgather(physf_node_z, phyg_z, plev)
      call sf_pgather(physf_node_d, phyg_d, plev)
      call sf_pgather(physf_node_h, phyg_h, plev)

      if(rank.eq.root)then
c         print*,'calling congen...'
c         call congen(phyg_q,icountp)
c         print*,'done congen'
         if (mod(time_current,time_out*100).lt.timestep)
     .        call sfout(phyg_z,phyg_d,phyg_h,icountp)
      end if

      end


