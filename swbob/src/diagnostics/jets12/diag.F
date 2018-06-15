      subroutine diag(itime,itime0)
      implicit none

#include <io.h>
#include <mpcluster.h>
#include <diag.h>
#include <diagptrs.h>
#include <modeltime.h>
#include <timestep.h>
#include <uspec.h>
#include <dims.h>
#include <planet.h>
#include <famp.h>

      REAL_TYPE skwq(3),skwc(3),cycl(3),Eph,Zph
      REAL_TYPE accinv,U1,E1,N1,Ro,Fr,Z1
      REAL_TYPE Eglb(4),Zglb(4),Ecnt(4),Enrm(4)
      INT_TYPE itime,itime0,icount,i,n,nacc,nfor
      common /diagctr/icount,nacc

ccccccccccccccccccccccccccccccccccccccccccccccc

      nfor=(nf1+nf0)/2

      do i=1,4
         Eglb(i)=0.d0
         Zglb(i)=0.d0
         Ecnt(i)=0.d0
         Enrm(i)=0.d0
         do n=0,nn
            Eglb(i)=Eglb(i)+Etot(n,i)
            Zglb(i)=Zglb(i)+Ztot(n,i)
            if ((abs(n-nfor)).gt.(nf1-nf0)) Ecnt(i)=Ecnt(i)+Etot(n,i)*float(n)
            if ((abs(n-nfor)).gt.(nf1-nf0)) Enrm(i)=Enrm(i)+Etot(n,i)
            Erotacc(n,i)=Erotacc(n,i)+Erot(n,i)
            Edivacc(n,i)=Edivacc(n,i)+Ediv(n,i)
            Etotacc(n,i)=Etotacc(n,i)+Etot(n,i)
            Ztotacc(n,i)=Ztotacc(n,i)+Ztot(n,i)
         enddo
         Eglbacc(i)=Eglbacc(i)+Eglb(i)
         Zglbacc(i)=Zglbacc(i)+Zglb(i)
         Urmsacc(i)=Urmsacc(i)+sqrt(2.0*Eglb(i))
         Ecntacc(i)=Ecntacc(i)+Ecnt(i)/Enrm(i)
      enddo
      nacc=nacc+1

c--------------------------------------
      if (itime.eq.itime0) then
         do i=1,4
            Urmsacc(i)=0.d0
            Eglbacc(i)=0.d0
            Zglbacc(i)=0.d0
            Ecntacc(i)=0.d0
            do n=0,nn
               Erotacc(n,i)=0.d0
               Edivacc(n,i)=0.d0
               Etotacc(n,i)=0.d0
               Ztotacc(n,i)=0.d0
            enddo
         enddo
         nacc=0
         if(itime0.eq.0)then
c            icount=0
c            call sf_diag(icount)   ! output IC
c            call zavg_diag(icount) ! calls bldzavg, scal, zavg_pgather
         else
            do i=1,4
               Eglb(i)=0.d0
               Zglb(i)=0.d0
               Ecnt(i)=0.d0
               Enrm(i)=0.d0
               do n=0,nn
                  Eglb(i)=Eglb(i)+Etot(n,i)
                  Zglb(i)=Zglb(i)+Ztot(n,i)
                  if ((abs(n-nfor)).gt.(nf1-nf0)) Ecnt(i)=Ecnt(i)+Etot(n,i)*float(n)
                  if ((abs(n-nfor)).gt.(nf1-nf0)) Enrm(i)=Enrm(i)+Etot(n,i)
                  Erotacc(n,i)=Erotacc(n,i)+Erot(n,i)
                  Edivacc(n,i)=Edivacc(n,i)+Ediv(n,i)
                  Etotacc(n,i)=Etotacc(n,i)+Etot(n,i)
                  Ztotacc(n,i)=Ztotacc(n,i)+Ztot(n,i)
               enddo
               Eglbacc(i)=Eglbacc(i)+Eglb(i)
               Zglbacc(i)=Zglbacc(i)+Zglb(i)
               Urmsacc(i)=Urmsacc(i)+sqrt(2.0*Eglb(i))
               Ecntacc(i)=Ecntacc(i)+Ecnt(i)/Enrm(i)
            enddo
            nacc=nacc+1
         endif
#ifdef ADAPT
         icount=nint(time_start+time_out)
#else
         icount=nint(time_start)+1
#endif
      else

         call zavg_diag(itime,icount,nacc)

#ifdef ADAPT
c         if (mod(time_current,time_out).lt.timestep) then
c         write(*,*) time_current,mod(time_current,time_out),dt2/2.d0
         if (mod(time_current,time_out).lt.(dt2/2.d0)) then
#else
         if (mod(itime,zavgfreq).eq.0) then
#endif

            call sf_diag(icount,skwq,skwc,cycl,Eph,Zph)
c            call zavg_diag(icount)

            accinv=1.d0/float(nacc)
            if (rank.eq.root) then
               write(85,*) time_current
               write(86,*) time_current
               write(87,*) time_current
               write(88,*) time_current
               do n=0,nn
                  write(85,8) n,(0.5*accinv*Erotacc(n,i),i=1,4)
                  write(86,8) n,(0.5*accinv*Edivacc(n,i),i=1,4)
                  write(87,8) n,(0.5*accinv*Etotacc(n,i),i=1,4)
                  write(88,8) n,(0.5*accinv*Ztotacc(n,i),i=1,4)
               enddo
               write(89,9) time_current,(accinv*Urmsacc(i),i=1,4)
               write(90,9) time_current,(skwq(i),i=1,3)
               write(91,9) time_current,(skwc(i),i=1,3)
               write(92,9) time_current,(cycl(i),i=1,3)
               write(93,9) time_current,Eph
               write(94,9) time_current,Zph
               U1=accinv*Urmsacc(1)
               E1=accinv*Eglbacc(1)
               Z1=accinv*Zglbacc(1)
               N1=accinv*Ecntacc(1)
               Ro=U1/(2.d0*omega*radius)
               Fr=U1/(2.d0*omega*radius*Ldeform)
               write(*,10)  '   Ro = ',Ro
               write(81,'(f11.4,5(2x,f9.6))') time_current,Ro,Fr,U1,E1,N1
            open(21,file='param0.dat',form='formatted',status='unknown') 
               rewind(21)
               write(21,'(a64)'),
     .    'L_D        Ro         Fr         Urms       E          nmax'
               write(21,'6(2x,f9.6)') Ldeform,Ro,Fr,U1,E1,N1
               close(21)
            endif

            do i=1,4
               do n=0,nn
                  Erotacc(n,i)=0.d0
                  Edivacc(n,i)=0.d0
                  Etotacc(n,i)=0.d0
                  Ztotacc(n,i)=0.d0
               enddo
               Eglbacc(i)=0.d0
               Zglbacc(i)=0.d0
               Ecntacc(i)=0.d0
               Urmsacc(i)=0.d0
            enddo
            nacc=0
            
c            icount=icount+1
            icount=icount+nint(time_out)

         endif

      endif

 8    format(i5,2x,4(1x,e14.6))
 9    format(f11.4,2x,4(1x,e14.6))
 10   format(a10,f9.6)

      end 
