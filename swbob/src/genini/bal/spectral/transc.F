      subroutine transc(sc,m,
     $     revorsc,redivsc,rephisc)
      implicit none

#include <type.h>
#include <dims.h>
c      parameter(plev=PLEV)
c      parameter(mm=(plon-1)/3,nn=mm)
c      parameter(nphy=5,nnlt=6,nsp=3)
c      parameter(nvsc=mm2d*(nn+4),nsc=mm2d*(nn+2))
#include <scptrs.h>
c      INT_TYPE scgptr(0:mm),scptr(0:mm),vscptr(0:mm)
c      common /scptrs/scgptr,scptr,vscptr
#include <stptrs.h>
c      vor_p etc ...

c-----------------------------------------------
c     Local declaration
c-----------------------------------------------

      INT_TYPE k,m,n,ic
      INT_TYPE  sc_p            
      REAL_TYPE sc(0:1,nsp,m:nn,plev)
      REAL_TYPE revorsc(0:1,nsc*plev)
      REAL_TYPE redivsc(0:1,nsc*plev)
      REAL_TYPE rephisc(0:1,nsc*plev)

c-----------------------------------------------
c     executable statements
c-----------------------------------------------

      do k=1,plev
         do n=m,nn
            sc_p=1+plev*scptr(m)+(nn+1-m)*(k-1)+n-m
            do ic=0,1
               revorsc(ic,sc_p)=sc(ic,vor_p,n,k)
               redivsc(ic,sc_p)=sc(ic,div_p,n,k)
               rephisc(ic,sc_p)=sc(ic,phi_p,n,k)
c     revorsc(ic,sc_p)=0.d0
c     redivsc(ic,sc_p)=0.d0
c     rephisc(ic,sc_p)=0.d0
            end do
         end do
      end do

      end
