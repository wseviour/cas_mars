      subroutine zerofor(scflt)

c     zero vorticity spectral coefficients close to forcing range

      implicit none

#include <dims.h>
#include <stptrs.h>
#include <scptrs.h>
#include <decomp.h>
#include <mpcluster.h>
#include <fordamp0.h>

      REAL_TYPE scflt(0:1,nsp,nsc*plev)
      INT_TYPE m,mch,k,sc_p,n,nf00,nf11,fhwid

c      fhwid=(nf1-nf0)/2
      fhwid=(nf1-nf0)
      nf00=nf0-fhwid
      nf11=nf1+fhwid

#ifdef GAUSS
#else
      do mch=0,1
         do m=mbeg(mch,rank),mend(mch,rank)
            do k=1,plev
               sc_p=1+plev*scptr(m)+(nn+1-m)*(k-1)
               do n=max(m,nf00),nf11
                  scflt(0,vor_p,sc_p+n-m)=0.d0
                  scflt(1,vor_p,sc_p+n-m)=0.d0
               enddo
            enddo
         enddo
      enddo
#endif

      return
      end
