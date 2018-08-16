      subroutine sp_diag
      implicit none

#include <fields.h>
#include <physical.h>
#include <decomp.h>
#include <mpcluster.h>
#include <spectral.h>
#include <uspec.h>
#include <dims.h>
#include <planet.h>
#include <fordamp0.h>
#include <modeltime.h>

      REAL_TYPE scvor(2*(plev*nsp)*nsc)
      REAL_TYPE phvor(nlon,nlat2,2*nphy*plev,ielemd,jelemd)
      REAL_TYPE momq(3,3),momc(3,3),momt(3,3)
      REAL_TYPE skwq(3),skwc(3),cycl(3),Eph,Zph,Epht,Zpht
      REAL_TYPE accinv,U1,E1,N1,Ro,Fr,Z1,tdays
      REAL_TYPE Eglb(4),Zglb(4),Ecnt(4),Enrm(4),Urms(4)
      INT_TYPE  je,jeg,ie,icountp,j,i,im,n,nfor



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

c     calculate |vorticity| moments and global energy/enstrophy
      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1
         do ie=1,ielem(jeg)
            call bldsp(phy(1,1,1,ie,je),phvor(1,1,1,ie,je),
     $                   momq,momc,Eph,Zph,
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
      Eph=0.5d0*Eph
      Zph=0.5d0*Zph

c     skewness etc:
      do i=1,3
         skwq(i)=momq(3,i)/sqrt(momq(2,i)**3)
         skwc(i)=(momc(3,i)-3.0*momc(2,i)*momc(1,i)+2.0*momc(1,i)**3)
     .           /sqrt((momc(2,i)-momc(1,i)**2)**3)
         cycl(i)=momc(1,i)
      enddo

c     global quantities from spectra:
      nfor=(nf1+nf0)/2
      do i=1,4
         Eglb(i)=0.d0
         Zglb(i)=0.d0
         Ecnt(i)=0.d0
         Enrm(i)=0.d0
         do n=0,nn
            Eglb(i)=Eglb(i)+Etot(n,i)
            Zglb(i)=Zglb(i)+Ztot(n,i)
#ifdef GAUSS
            Ecnt(i)=Ecnt(i)+Etot(n,i)*float(n)
            Enrm(i)=Enrm(i)+Etot(n,i)
#else
            if ((abs(n-nfor)).gt.(nf1-nf0)) Ecnt(i)=Ecnt(i)+Etot(n,i)*float(n)
            if ((abs(n-nfor)).gt.(nf1-nf0)) Enrm(i)=Enrm(i)+Etot(n,i)
#endif
         enddo
         Urms(i)=sqrt(2.0*Eglb(i))
         Ecnt(i)=Ecnt(i)/Enrm(i)
      enddo

c     write to file and stdout
      if (rank.eq.root) then
         tdays=time_current/SecPerDay
         write(85,*) tdays
         write(86,*) tdays
         write(87,*) tdays
         write(88,*) tdays
         do n=0,nn
            write(85,8) n,(0.5*Erot(n,i),i=1,4)
            write(86,8) n,(0.5*Ediv(n,i),i=1,4)
            write(87,8) n,(0.5*Etot(n,i),i=1,4)
            write(88,8) n,(0.5*Ztot(n,i),i=1,4)
         enddo
         write(89,9) tdays,(Urms(i),i=1,4)
         write(90,9) tdays,(skwq(i),i=1,3)
         write(91,9) tdays,(skwc(i),i=1,3)
         write(92,9) tdays,(cycl(i),i=1,3)
         write(93,9) tdays,Eph
         write(94,9) tdays,Zph
         do i=85,94
            call flush(i)
         enddo

         U1=Urms(1)
         E1=Eglb(1)
         Z1=Zglb(1)
         N1=Ecnt(1)
         Ro=U1/(2.d0*omega*radius)
         Fr=U1/(2.d0*omega*radius*Ldeform)
         write(*,10)  '   Ro = ',Ro
         write(81,'(f11.4,5(2x,f9.6))') tdays,Ro,Fr,U1,E1,N1
         open(21,file='param0.dat',form='formatted',status='unknown') 
         rewind(21)
         write(21,'(a64)'),
     .     'L_D        Ro         Fr         Urms       E          nmax'
         write(21,'(6(2x,f9.6))') Ldeform,Ro,Fr,U1,E1,N1
         close(21)
         call flush(81)
      endif
      
 8    format(i5,2x,4(1x,e14.6))
 9    format(f11.4,2x,4(1x,e14.6))
 10   format(a10,f9.6)

      end


