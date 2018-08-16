      subroutine bldrzd(sc,
     $                 nltsc,
     $                 eps,
     $                 m)
      implicit none

!--------------------------------------------------------------------
!      Include Files
!--------------------------------------------------------------------
#include <dims.h>
#include <nltptrs.h>
#include <stptrs.h>
#include <specuv.h>
#include <uspec.h>
!--------------------------------------------------------------------

      INT_TYPE  m
      REAL_TYPE sc(0:1,nsp,m:nn)
      REAL_TYPE nltsc(0:1,nnlt,m:nn+1)
      REAL_TYPE eps(m:nn+1)
  
      REAL_TYPE psi(0:1,0:nn)
      REAL_TYPE chi(0:1,0:nn)
      REAL_TYPE w
      REAL_TYPE zeta
      REAL_TYPE afac,zfac

      INT_TYPE  n,ic,i

      ! added by leo oct 30th 1999
      do n=m,nn+1
         do ic=0,1
            nltsc(ic,u_p,n)=0.0D0
            nltsc(ic,v_p,n)=0.0D0
         enddo
      enddo

      ! pretend like coreolis term doesn't exist (vor=zeta)

      do n=m,nn
         do ic=0,1
            zeta=sc(ic,vor_p,n)
            psi(ic,n) = -annp1(n)*zeta
            chi(ic,n) = -annp1(n)*sc(ic,div_p,n)
         end do
      end do
 
      ! actually vor and zeta differ only in coefficients m=1,n=0,1 
      ! where the coreolis factor contributes.
      ! Subtract fcor from vor for these wave numbers to obtain zeta, 
      ! then build psi with it

      n=1

      if (m.le.1) then
         zeta=sc(0,vor_p,n)-fcor(m) 
         psi(0,n) = -annp1(n)*zeta
      end if

      w=m

      ! Code according to eqs 3.11, 3.12 p33 in R. Jakob thesis

      do n=m,nn
         nltsc(0,u_p,n) = -w*chi(1,n)
         nltsc(1,u_p,n) =  w*chi(0,n)
         nltsc(0,v_p,n) = -w*psi(1,n)
         nltsc(1,v_p,n) =  w*psi(0,n)
      end do

      do n=m,nn-1
         nltsc(0,u_p,n) = nltsc(0,u_p,n)-(n+2)*eps(n+1)*psi(0,n+1)
         nltsc(1,u_p,n) = nltsc(1,u_p,n)-(n+2)*eps(n+1)*psi(1,n+1)
         nltsc(0,v_p,n) = nltsc(0,v_p,n)+(n+2)*eps(n+1)*chi(0,n+1)
         nltsc(1,v_p,n) = nltsc(1,v_p,n)+(n+2)*eps(n+1)*chi(1,n+1)
      end do

      do n=m+1,nn+1
         nltsc(0,u_p,n) = nltsc(0,u_p,n)+(n-1)*eps(n)*psi(0,n-1)
         nltsc(1,u_p,n) = nltsc(1,u_p,n)+(n-1)*eps(n)*psi(1,n-1)
         nltsc(0,v_p,n) = nltsc(0,v_p,n)-(n-1)*eps(n)*chi(0,n-1)
         nltsc(1,v_p,n) = nltsc(1,v_p,n)-(n-1)*eps(n)*chi(1,n-1)
      end do

      end
