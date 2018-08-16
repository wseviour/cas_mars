#define SC_TAG_BASE 1
      subroutine sc_pgather(sc,scg,nlev,k)
      implicit none

#include <dims.h>
#include <scptrs.h>
#include <decomp.h>
#include <mpcluster.h>

      INT_TYPE  nlev
      INT_TYPE  k
      REAL_TYPE sc(0:1,nsc*plev)  
      REAL_TYPE scbuf(0:1,nsc)
      REAL_TYPE scg(0:1,nscg)

      INT_TYPE ipe,mch,m
      INT_TYPE sc_p
      INT_TYPE rtag,stag

      call synch()
      
      do mch=0,1
         do m=mbeg(mch,rank),mend(mch,rank)
            sc_p=1+nlev*scptr(m)
            call sc_xfer_layer(sc(0,sc_p),scbuf(0,scptr(m)+1),m,nlev,k)
         end do
      end do

      if (rank.eq.root) then

         do ipe=0,P_NODE-1

            if (ipe.ne.root) then
               rtag=ipe+SC_TAG_BASE*P_NODE
               call recvreal(scbuf,2*nsc,ipe,rtag)
            end if

            do mch=0,1
               do m=mbeg(mch,ipe),mend(mch,ipe)
                  call sc_xfer_layer(scbuf(0,scptr(m)+1),scg(0,scgptr(m)+1),m,1,1)
               end do
            end do
            
         end do

      else

         stag=rank+SC_TAG_BASE*P_NODE
         call sendreal(scbuf,2*nsc,root,stag)

      end if

      end

