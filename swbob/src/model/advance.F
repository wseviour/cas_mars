#include <type.h>
      subroutine advance(nltsc,
     $                   sc,
     $                   scm1,
     $                   fthesc,
     $                   ftoposc,
     $                   ftopoktsc,
     $                   fvorsc,
     $                   fthpsc,
     $                   eps,
     $                   m)
      implicit none

!--------------------------------------------------------------------
!      Include Files
!--------------------------------------------------------------------
#include <dims.h>
#include <planet.h>          ! now includes phibar (7/30/05)
#include <nltptrs.h>
#include <stptrs.h>
#include <timestep.h>
#include <timeint.h>
#include <modeltime.h>
#include <hyper.h>
#include <fordamp.h>     
#include <fordamp0.h>

      INT_TYPE  m
      REAL_TYPE nltsc(0:1,nnlt,m:nn+1,plev)
      REAL_TYPE sc(0:1,nsp,m:nn,plev)
      REAL_TYPE scm1(0:1,nsp,m:nn,plev)
      REAL_TYPE fthesc(0:1,m:nn,plev)
      REAL_TYPE ftoposc(0:1,m:nn,plev)
      REAL_TYPE ftopoktsc(0:1,m:nn,plev)
      REAL_TYPE fvorsc(0:1,m:nn,plev)
      REAL_TYPE fthpsc(0:1,m:nn,plev)
      REAL_TYPE eps(m:nn+1)

      INT_TYPE ic, n, k, l

      REAL_TYPE thterm(0:1)
      REAL_TYPE vorterm(0:1)
      REAL_TYPE divterm(0:1)
      REAL_TYPE scp1(0:1,nsp,0:nn)

      REAL_TYPE w, dt, dtinv
      REAL_TYPE flap, hdutmp, hdttmp, hducor, hdtcor
      REAL_TYPE Sinv,Sm
#ifdef FTOPO
      REAL_TYPE rampdays,ramp,delay,rmpc
      common /toporamp/ ramp,rmpc
#endif
      REAL_TYPE phtc,phts,phhc,phhs,targ

c----------------------------------------
c     Executable Statements
c----------------------------------------

      dt=dt2/2.d0
      dtinv=1.d0/dt
      w = m

c      if(m.eq.0) print*,'in advance; ',itime
      if (m.eq.0) then 
         write(6,'(a1,$)')'.'
         call flush(6)
      endif


#ifdef M2
c     all odd wavenumbers = 0 with wave-2 symmetry
      if (mod(m,2).eq.0) then
#endif

#ifdef FTOPO
      delay=0.d0*SecPerDay
c      rampdays=100.d0*SecPerDay
c changing to instantaneous topography - wjms
      rampdays=0.d0*SecPerDay
c      if (time_current.lt.rampdays) then
c         ramp=(dsin( (time_current/rampdays)*pi/2.d0 ) )**2
      if (time_current.lt.delay) then
         ramp=0.d0
      else if (time_current.lt.(delay+rampdays)) then
         ramp=(dsin( ((time_current-delay)/rampdays)*pi/2.d0 ) )**2
      else
         ramp=1.d0
      endif
      rmpc=1.d0-ramp

c linear in time:
c      ramp=time_current/time_stop

c tmp: level off at constant amplitude:  
c      ramp=min(ramp,2.d0*160.d0/1000.d0)

c tmp - pulsed topography:
c      ramp=(dsin( (time_current/rampdays)*pi/2.d0 ) )**2
      targ=freqt*time_current
      phtc=dcos(targ)
      phts=dsin(targ)
#endif
#ifdef FVORT
      targ=freqt*time_current
      phtc=dcos(targ)
      phts=dsin(targ)
#endif
#ifdef FHEAT
      targ=freqh*time_current
      phhc=dcos(targ)
      phhs=dsin(targ)
#endif
      do k=1,plev

         do n=m,nn

            flap        =  (n*(n+1))/(Radius*Radius)
            hdutmp      =  hdufac(n)
            hdttmp      =  hdtfac(n)
            hdtcor      =  1.D0/( 1.D0+dt2*hdttmp )
            hducor      =  1.D0/( 1.D0+dt2*hdutmp )
            Sinv        =  1.0d0 / ( 1.0d0 + hducor*hdtcor*flap*PhiBar*dt*dt )
            Sm          =  ( 1.0d0 - hducor*hdtcor*flap*PhiBar*dt*dt )

            if (n.eq.m) then
               thterm(0)  =  w*nltsc(1,c_p,n,k)-n*eps(n+1)*nltsc(0,d_p,n+1,k) - hdttmp*scm1(0,theta_p,n,k)
               thterm(1)  = -w*nltsc(0,c_p,n,k)-n*eps(n+1)*nltsc(1,d_p,n+1,k) - hdttmp*scm1(1,theta_p,n,k)

               divterm(0) = -w*nltsc(1,b_p,n,k)-n*eps(n+1)*nltsc(0,a_p,n+1,k) + flap*nltsc(0,e_p,n,k) - hdutmp*scm1(0,div_p,n,k)
c#ifdef FTOPO
c     $                                        + flap*ramp*ftoposc(0,n,k)
c#endif
               divterm(1) =  w*nltsc(0,b_p,n,k)-n*eps(n+1)*nltsc(1,a_p,n+1,k) + flap*nltsc(1,e_p,n,k) - hdutmp*scm1(1,div_p,n,k)
c#ifdef FTOPO
c     $                                        + flap*ramp*ftoposc(1,n,k)
c#endif

               vorterm(0) =  w*nltsc(1,a_p,n,k)-n*eps(n+1)*nltsc(0,b_p,n+1,k) - hdutmp*scm1(0,vor_p,n,k)
               vorterm(1) = -w*nltsc(0,a_p,n,k)-n*eps(n+1)*nltsc(1,b_p,n+1,k) - hdutmp*scm1(1,vor_p,n,k)

            else
               thterm(0)  =  w*nltsc(1,c_p,n,k)-n*eps(n+1)*nltsc(0,d_p,n+1,k) +
     $                                        (n+1)*eps(n)*nltsc(0,d_p,n-1,k) - hdttmp*scm1(0,theta_p,n,k)
               thterm(1)  = -w*nltsc(0,c_p,n,k)-n*eps(n+1)*nltsc(1,d_p,n+1,k) +
     $                                        (n+1)*eps(n)*nltsc(1,d_p,n-1,k) - hdttmp*scm1(1,theta_p,n,k)

               divterm(0) = -w*nltsc(1,b_p,n,k)-n*eps(n+1)*nltsc(0,a_p,n+1,k) +
     $                                        (n+1)*eps(n)*nltsc(0,a_p,n-1,k) + flap*nltsc(0,e_p,n,k) - hdutmp*scm1(0,div_p,n,k)
c#ifdef FTOPO
c     $                                        + flap*ramp*ftoposc(0,n,k)
c#endif
               divterm(1) =  w*nltsc(0,b_p,n,k)-n*eps(n+1)*nltsc(1,a_p,n+1,k) +
     $                                        (n+1)*eps(n)*nltsc(1,a_p,n-1,k) + flap*nltsc(1,e_p,n,k) - hdutmp*scm1(1,div_p,n,k)
c#ifdef FTOPO
c     $                                        + flap*ramp*ftoposc(1,n,k)
c#endif

               vorterm(0) =  w*nltsc(1,a_p,n,k)-n*eps(n+1)*nltsc(0,b_p,n+1,k) +
     $                                        (n+1)*eps(n)*nltsc(0,b_p,n-1,k) - hdutmp*scm1(0,vor_p,n,k)
               vorterm(1) = -w*nltsc(0,a_p,n,k)-n*eps(n+1)*nltsc(1,b_p,n+1,k) +
     $                                        (n+1)*eps(n)*nltsc(1,b_p,n-1,k) - hdutmp*scm1(1,vor_p,n,k)

            end if

c     forcing options:
c     ( NOTE - momentum forcing added to A,B terms in bldnlt )
#ifdef FDRAG
c     vorticity damping (linear drag):
#if 1
c     damping on all wavenumbers
            vorterm(0)= vorterm(0) - nltsc(0,f_p,n,k)
            vorterm(1)= vorterm(1) - nltsc(1,f_p,n,k)
#endif
#if 0
c     damping on m=0 only
            if (m.eq.0) then
               vorterm(0)= vorterm(0) - nltsc(0,f_p,n,k)
            endif
#endif
#endif
#ifdef FRLAX
c     thermal relaxation (Newtonian cooling):
c            if (m.eq.0) write(*,*)fthesc(0,n,k)
#if 1
c     relaxation on all wavenumbers
            thterm(0)= thterm(0) + fthesc(0,n,k) - nltsc(0,g_p,n,k)
            thterm(1)= thterm(1) + fthesc(1,n,k) - nltsc(1,g_p,n,k)
#endif
#if 0
c     relaxation on all m=0 only
            if (m.eq.0) then
               thterm(0)= thterm(0) + fthesc(0,n,k) - nltsc(0,g_p,n,k)
            endif
#endif
#if 0
c     relaxation on all m>0 only
            if (m.ne.0) then
               thterm(0)= thterm(0) + fthesc(0,n,k) - nltsc(0,g_p,n,k)
               thterm(1)= thterm(1) + fthesc(1,n,k) - nltsc(1,g_p,n,k)
            endif
#endif
#ifdef FRLAXZ
c     relax h+htopo to he
            thterm(0)= thterm(0) - ramp*(phtc*ftopoktsc(0,n,k) + phts*ftopoktsc(1,n,k))
            thterm(1)= thterm(1) - ramp*(phtc*ftopoktsc(1,n,k) - phts*ftopoktsc(0,n,k))
#endif
#endif
#ifdef FTOPO
c     topographic forcing:
            divterm(0)= divterm(0) + flap*ramp*(phtc*ftoposc(0,n,k) + phts*ftoposc(1,n,k))
            divterm(1)= divterm(1) + flap*ramp*(phtc*ftoposc(1,n,k) - phts*ftoposc(0,n,k))


c    !!! TEMPORARY !!!
c     divergence source (no laplacian):
c            divterm(0)= divterm(0) + ramp*(phtc*ftoposc(0,n,k) + phts*ftoposc(1,n,k))
c            divterm(1)= divterm(1) + ramp*(phtc*ftoposc(1,n,k) - phts*ftoposc(0,n,k))


#endif
#ifdef FHEAT
c     thermal forcing (asymmetric):
            thterm(0)= thterm(0) + phhc*fthpsc(0,n,k) + phhs*fthpsc(1,n,k)
            thterm(1)= thterm(1) + phhc*fthpsc(1,n,k) - phhs*fthpsc(0,n,k)
#endif
#ifdef FVORT
c     vorticity forcing (asymmetric):
            vorterm(0)= vorterm(0) + phtc*ftoposc(0,n,k) + phts*ftoposc(1,n,k)
            vorterm(1)= vorterm(1) + phtc*ftoposc(1,n,k) - phts*ftoposc(0,n,k)
#endif
#ifdef FSPEC
c     vorticity forcing (stochastic):
            vorterm(0)= vorterm(0) + fvorsc(0,n,k)*ramp
            vorterm(1)= vorterm(1) + fvorsc(1,n,k)*ramp
#endif

            scp1(0,vor_p,n)=scm1(0,vor_p,n,k)+dt2*hducor*vorterm(0)
            scp1(1,vor_p,n)=scm1(1,vor_p,n,k)+dt2*hducor*vorterm(1)

            scp1(0,div_p,n)=Sinv*(Sm*scm1(0,div_p,n,k)
     $           +dt2*hducor*(divterm(0) + flap*(hdtcor*dt*thterm(0) + scm1(0,theta_p,n,k) )))
            scp1(1,div_p,n)=Sinv*(Sm*scm1(1,div_p,n,k)
     $           +dt2*hducor*(divterm(1) + flap*(hdtcor*dt*thterm(1) + scm1(1,theta_p,n,k) )))

            scp1(0,theta_p,n)=Sinv*(Sm*scm1(0,theta_p,n,k)
     $           +dt2*hdtcor*(thterm(0) - PhiBar*(hducor*dt*divterm(0) + scm1(0,div_p,n,k) )))
            scp1(1,theta_p,n)=Sinv*(Sm*scm1(1,theta_p,n,k)
     $           +dt2*hdtcor*(thterm(1) - PhiBar*(hducor*dt*divterm(1) + scm1(1,div_p,n,k) )))

         end do ! Loop in n


c     time filter:
         do n=m,nn
            if (dofilter) then  
               do l=1,nsp
                  do ic=0,1
                     scm1(ic,l,n,k)=(1.0d0-2.d0*tfc)*sc(ic,l,n,k) + 
     &                    tfc*(scm1(ic,l,n,k) + scp1(ic,l,n))
                  end do
               end do
            end if
            sc(0,vor_p,n,k)=scp1(0,vor_p,n)
            sc(1,vor_p,n,k)=scp1(1,vor_p,n)
            sc(0,div_p,n,k)=scp1(0,div_p,n)
            sc(1,div_p,n,k)=scp1(1,div_p,n)
            sc(0,theta_p,n,k)=scp1(0,theta_p,n)
            sc(1,theta_p,n,k)=scp1(1,theta_p,n)
         end do

c     compute u and v from absvor and div:
         call blduv(sc(0,1,m,k),
     $              nltsc(0,1,m,k),
     $              eps,
     $              m)

      end do     ! end loop over levels


#ifdef M2
      else
c     odd zonal wavenumbers=0:
      do k=1,plev
         do n=m,nn
            sc(0,vor_p,n,k)=0.d0
            sc(1,vor_p,n,k)=0.d0
            sc(0,div_p,n,k)=0.d0
            sc(1,div_p,n,k)=0.d0
            sc(0,theta_p,n,k)=0.d0
            sc(1,theta_p,n,k)=0.d0
         enddo
         call blduv(sc(0,1,m,k),
     $              nltsc(0,1,m,k),
     $              eps,
     $              m)
      enddo
      endif
#endif


      end




