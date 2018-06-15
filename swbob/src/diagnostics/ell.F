      subroutine ell
      implicit none

#include <fields.h>          ! phy
#include <physical.h>        ! cslat(nlat2,jelemg),ielem(jelemg)
#include <decomp.h>          ! jebeg(0:P_NODE-1),jeend(0:P_NODE-1)
#include <mpcluster.h>       ! rank
#include <spectral.h>        ! sc(2*(plev*nsp)*nsc)
#include <modeltime.h>
#include <timestep.h>
#include <planet.h>

      INT_TYPE  navu
      external  navu

      REAL_TYPE area,areat,mm2(3),mm2t(3),mm4(3),mm4t(3)
      INT_TYPE  je,jeg,ie,i,j,iunit
      REAL_TYPE dummy,dq0,dq1,qc
      REAL_TYPE tdays


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Executable statements
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c     get midvalue of vortex edge from genx.in:
      iunit=navu()
      open(iunit,file='genx.in',form='formatted')
      read(iunit,*),dummy,dummy,dq0,dq1,qc
      close(iunit)
c      qc=0.5d0*(dq0+dq1)*2.d0*Omega
      qc=qc*2.d0*Omega

      area=0.d0
      do i=1,3
         mm2(i)=0.d0
         mm4(i)=0.d0
      enddo

      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1
         do ie=1,ielem(jeg)
            call ell_bld(phy(1,1,1,ie,je),
     $                   lon(1,ie,je),
     $                   lat(1,jeg),
     $                   da(1,jeg),
     $                   mm2,mm4,area,qc)
         enddo
      enddo

      call psumreal(area,areat,1)
      call broadcastreal(areat,1)
      call copyreal(areat,area,1)
      call psumreal(mm2,mm2t,3)
      call broadcastreal(mm2t,3)
      call copyreal(mm2t,mm2,3)
      call psumreal(mm4,mm4t,3)
      call broadcastreal(mm4t,3)
      call copyreal(mm4t,mm4,3)

      area=area*4.*pi
      do i=1,3
         mm2(i)=mm2(i)*4.*pi
         mm4(i)=mm4(i)*4.*pi
      enddo

      if (rank.eq.root) then
         tdays=time_current/SecPerDay
         write(95,8) tdays,area,(mm2(i),i=1,3),(mm4(i),i=1,3)
      endif
      call flush(95)

 8    format(f8.3,2x,7(1x,1p,e14.6))

      end


