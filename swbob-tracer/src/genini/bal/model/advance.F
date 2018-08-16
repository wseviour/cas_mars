#include <type.h>
      subroutine advance(nltsc,
     $                   sc,
     $                   eps,
     $                   m)
      implicit none

!--------------------------------------------------------------------
!      Include Files
!--------------------------------------------------------------------
#include <dims.h>
#include <planet.h>
#include <nltptrs.h>
#include <stptrs.h>
#include <specuv.h>

      INT_TYPE  m
      REAL_TYPE nltsc(0:1,nnlt,m:nn+1,plev)
      REAL_TYPE sc(0:1,nsp,m:nn,plev)
      REAL_TYPE eps(m:nn+1)

      INT_TYPE ic  ! re/im index
      INT_TYPE n   ! legendre order number
      INT_TYPE k   ! level index
      INT_TYPE l   ! state variable index

      REAL_TYPE divterm(0:1)

      REAL_TYPE w
      REAL_TYPE flap,flapinv

c----------------------------------------
c     Executable Statements
c----------------------------------------

      w = m

      do k=1,plev
         do n=m,nn

            flap        =  (n*(n+1))/(Radius*Radius)
            flapinv = (Radius*Radius)/(n*(n+1))

            if (n.eq.m) then
               divterm(0) = -w*nltsc(1,b_p,n,k)-n*eps(n+1)*nltsc(0,a_p,n+1,k) + flap*nltsc(0,e_p,n,k)
               divterm(1) =  w*nltsc(0,b_p,n,k)-n*eps(n+1)*nltsc(1,a_p,n+1,k) + flap*nltsc(1,e_p,n,k)
            else
               divterm(0) = -w*nltsc(1,b_p,n,k)-n*eps(n+1)*nltsc(0,a_p,n+1,k) +
     $              (n+1)*eps(n)*nltsc(0,a_p,n-1,k) + flap*nltsc(0,e_p,n,k)
               divterm(1) =  w*nltsc(0,b_p,n,k)-n*eps(n+1)*nltsc(1,a_p,n+1,k) +
     $              (n+1)*eps(n)*nltsc(1,a_p,n-1,k) + flap*nltsc(1,e_p,n,k)
            end if

            if(n.eq.0) then
               sc(0,theta_p,n,k)=0.d0
               sc(1,theta_p,n,k)=0.d0
            else
               sc(0,theta_p,n,k)=sc(0,theta_p,n,k)-flapinv*divterm(0)
               sc(1,theta_p,n,k)=sc(1,theta_p,n,k)-flapinv*divterm(1)
            endif

         end do ! Loop in n
      end do     ! end loop over levels

      end




