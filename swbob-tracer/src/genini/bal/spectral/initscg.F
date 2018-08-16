
      subroutine initscg(scg,ilev,ifld)

      implicit none

#include <dims.h>
#include <stptrs.h>
#include <scptrs.h>
#include <mpcluster.h>
#include <planet.h>
#include <io.h>
#include <vordist.h>
#include <specuv.h>

      REAL_TYPE scg(0:1,nscg)
      REAL_TYPE edensity(0:nn),vornsq(0:nn)
      REAL_TYPE rnda,rndb,an,uni,pm

      INT_TYPE  ifld,ilev,scg_p,mch,m,n
      INT_TYPE  iseed,i
      INT_TYPE  ounit               ! file output unit number

c-------------------------------------------------------------
c     Executable Statements
c-------------------------------------------------------------
c     
      if (rank.eq.root) then
         iseed=ilev
         call random_seed(iseed)
c
         if (ifld.eq.vor_p) then 
c
            do n=0,nn
               an=float(n)
               edensity(n)=amp*an**(gamma/2.d0)/(an+n0)**gamma
               vornsq(n)=an*(an+1)/(Radius**2*(2.d0*an+1.d0))*edensity(n)
            enddo
c
            do n=0,nn
               do m=0,n
                  call random_number(uni)
                  call random_number(pm)
                  rnda=(2.d0*uni-1.d0)
                  rndb=sign(1.d0,2.d0*pm-1.d0)*sqrt(1.d0-rnda**2)
                  scg_p=1+scgptr(m)
                  scg(0,scg_p+n-m)=rnda*vornsq(n)
                  scg(1,scg_p+n-m)=rndb*vornsq(n)
                  if (n.eq.1) scg(0,scg_p+n-m)=scg(0,scg_p+n-m)+fcor(m)
               enddo
            enddo

         endif

         if (ifld.eq.div_p) then 
            do n=1,nscg
               scg(0,n)=0.d0
               scg(1,n)=0.d0
            enddo
         endif
         if (ifld.eq.phi_p) then 
            do n=1,nscg
               scg(0,n)=0.d0
               scg(1,n)=0.d0
            enddo
         endif

      endif

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      end






