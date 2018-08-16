      subroutine forsc(fvorsc,sc)

c     called from spectral.F

      implicit none

#include <dims.h>
#include <stptrs.h>
#include <scptrs.h>
#include <decomp.h>
#include <mpcluster.h>
#include <planet.h>
#include <fordamp0.h>
#include <timestep.h>
#include <specuv.h>

      REAL_TYPE sc(0:1,nsp,nsc*plev)
      REAL_TYPE fvorsc(0:1,nsc*plev)
      REAL_TYPE ftmpsc(0:1,nsc*plev)
      REAL_TYPE r0,r1,rnorm
      REAL_TYPE v0,v1,vnorm,psi0,psi1,psif,psift,psistart
      REAL_TYPE fspec(0:nn)
      INT_TYPE m,mch,k,sc_p,n

      REAL_TYPE zero,one,small
      parameter(zero=0.d0,one=1.d0,small=1.0d-12)

      REAL_TYPE ffac

c     time-correlation (gamma=correlation radius):
      fcor0=1.d0-timestep/max(gamma,timestep)
      fcor1=sqrt(1.d0-fcor0**2)

c     spectrum of forcing
      do n=0,nn
#ifdef GAUSS         
         fspec(n)=exp(-(float(n-nf0)/float(nf1))**2)
#else
         fspec(n)=1.d0
#endif
      enddo

c     update forcing:

      if (gamma.gt.timestep) then
c     write(*,*) 'markovian'

      psif=zero
      psift=zero
      do mch=0,1
         do m=mbeg(mch,rank),mend(mch,rank)
            do k=1,plev
               sc_p=1+plev*scptr(m)+(nn+1-m)*(k-1)
#ifdef GAUSS
               do n=m,nn
#else
               do n=max(m,nf0),nf1
#endif
                  call myrand(r0)
                  call myrand(r1)
                  rnorm=one/sqrt(r0**2+r1**2)
                  r0=r0*rnorm
                  r1=r1*rnorm
                  
                  v0=fvorsc(0,sc_p+n-m)
                  v1=fvorsc(1,sc_p+n-m)
                  vnorm=one/sqrt(v0**2+v1**2+small)
                  v0=v0*vnorm
                  v1=v1*vnorm
c                  vfor(0,n,m)=fcor0*v0+fcor1*r0
c                  vfor(1,n,m)=fcor0*v1+fcor1*r1
                  ftmpsc(0,sc_p+n-m)=fcor0*v0+fcor1*r0
                  ftmpsc(1,sc_p+n-m)=fcor0*v1+fcor1*r1
                  if ((n.eq.1).and.(m.le.1)) then
                     psi0=annp1(n)*(sc(0,vor_p,sc_p+n-m)-fcor(m))
                  else
                     psi0=annp1(n)*sc(0,vor_p,sc_p+n-m)
                  endif
                  psi1=annp1(n)*sc(1,vor_p,sc_p+n-m)
c     account for -n<=m<=n modes:
                  if (m.eq.0) then 
                     psif=psif+ftmpsc(0,sc_p+n-m)*psi0+ftmpsc(1,sc_p+n-m)*psi1
                  else
                     psif=psif+2.d0*(ftmpsc(0,sc_p+n-m)*psi0+ftmpsc(1,sc_p+n-m)*psi1)
                  endif
               enddo
            enddo
         enddo
      enddo
      call psumreal(psif,psift,1)
      call broadcastreal(psift,1)
      psif=psift
      psistart=0.0001d0
      if (abs(psif).gt.psistart) then
c         if (rank.eq.root) write(*,*) psif
         do mch=0,1
            do m=mbeg(mch,rank),mend(mch,rank)
               do k=1,plev
                  sc_p=1+plev*scptr(m)+(nn+1-m)*(k-1)
#ifdef GAUSS
                  do n=m,nn
#else
                  do n=max(m,nf0),nf1
#endif
                     fvorsc(0,sc_p+n-m)=eps0*fspec(n)*ftmpsc(0,sc_p+n-m)/psif
                     fvorsc(1,sc_p+n-m)=eps0*fspec(n)*ftmpsc(1,sc_p+n-m)/psif
                  enddo
               enddo
            enddo
         enddo
      else
c         if (rank.eq.root) write(*,*) 'psif,psistart = ',psif,psistart
         do mch=0,1
            do m=mbeg(mch,rank),mend(mch,rank)
               do k=1,plev
                  sc_p=1+plev*scptr(m)+(nn+1-m)*(k-1)
#ifdef GAUSS
                  do n=m,nn
#else
                  do n=max(m,nf0),nf1
#endif
                     fvorsc(0,sc_p+n-m)=eps0*fspec(n)*ftmpsc(0,sc_p+n-m)/psistart
                     fvorsc(1,sc_p+n-m)=eps0*fspec(n)*ftmpsc(1,sc_p+n-m)/psistart
                  enddo
               enddo
            enddo
         enddo
      endif
c
      else
c      write(*,*) 'delta'

      do mch=0,1
         do m=mbeg(mch,rank),mend(mch,rank)
            do k=1,plev
               sc_p=1+plev*scptr(m)+(nn+1-m)*(k-1)
#ifdef GAUSS
               do n=m,nn
#else
               do n=max(m,nf0),nf1
#endif
                  call myrand(r0)
                  call myrand(r1)
                  rnorm=one/sqrt(r0**2+r1**2)
                  r0=r0*rnorm
                  r1=r1*rnorm
c                  ffac=sqrt((n+1)*eps0/(timestep*(nf1-nf0)))  ! no 2pi; diff from planar
#ifdef GAUSS
              stop "forsc.F: implement gaussian forcing distribution"
#else
                  ffac=sqrt(2.0*n*(n+1)*eps0/(timestep*(nf1-nf0+1)*(2*n+1)))
#endif
                  fvorsc(0,sc_p+n-m)=ffac*r0
                  fvorsc(1,sc_p+n-m)=ffac*r1
               enddo
            enddo
         enddo
      enddo

      endif



C       psif=zero
C       psift=zero
C       do mch=0,1
C          do m=mbeg(mch,rank),mend(mch,rank)
C             do k=1,plev
C                sc_p=1+plev*scptr(m)+(nn+1-m)*(k-1)
C                do n=max(m,nf0),nf1
C                   if ((n.eq.1).and.(m.le.1)) then
C                      psi0=annp1(n)*(sc(0,vor_p,sc_p+n-m)-fcor(m))
C                   else
C                      psi0=annp1(n)*sc(0,vor_p,sc_p+n-m)
C                   endif
C                   psi1=annp1(n)*sc(1,vor_p,sc_p+n-m)
C                   psif=psif+fvorsc(0,sc_p+n-m)*psi0+fvorsc(1,sc_p+n-m)*psi1
C                   psift=psift+annp1(n)*(fvorsc(0,sc_p+n-m)**2+fvorsc(1,sc_p+n-m)**2)
C                enddo
C             enddo
C          enddo
C       enddo
C       psift=psift*0.5d0*timestep
C       write(*,*) ' psif,psift' , psif/eps0,psift/eps0

c      do mch=0,1
c         do m=mbeg(mch,rank),mend(mch,rank)
c            do k=1,plev
c               sc_p=1+plev*scptr(m)+(nn+1-m)*(k-1)
c               write(*,*) m,sc_p,fvorsc(0,sc_p+nf-m),fvorsc(1,sc_p+nf-m)
c            enddo
c         enddo
c      enddo


      return
      end
