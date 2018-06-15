      subroutine initscptrs()
      implicit none

#include <dims.h>
#include <decomp.h>
#include <mpcluster.h>
#include <scptrs.h>

      INT_TYPE m,ipe

      scgptr(0)=0
      do m=1,mm
         scgptr(m)=scgptr(m-1) + (nn+1) - (m-1)
      end do

      do ipe=0,M_NODE-1
         vscptr(mbeg(0,ipe))=0
         scptr(mbeg(0,ipe))=0
         do m=mbeg(0,ipe)+1,mend(0,ipe)
            vscptr(m) = vscptr(m-1) + (nn+2) - (m-1)
            scptr(m) =  scptr(m-1) + (nn+1) - (m-1)
         end do

         vscptr(mbeg(1,ipe))=vscptr(mend(0,ipe)) + (nn+2) -(mend(0,ipe))
         scptr(mbeg(1,ipe))= scptr(mend(0,ipe)) + (nn+1) -(mend(0,ipe))

         do m=mbeg(1,ipe)+1,mend(1,ipe)
            vscptr(m) = vscptr(m-1) + (nn+2) - (m-1)
            scptr(m) =  scptr(m-1) + (nn+1) - (m-1)
         end do
      end do

      end 

