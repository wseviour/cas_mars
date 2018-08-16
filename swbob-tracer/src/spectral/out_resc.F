      subroutine out_resc(lrevorsc,reounit,fld_file)

c     Calls: spectral/sc_pgather, utils/fldio/writereal
c------------------------------------------------
      implicit none
#include <type.h>
#include <dims.h>
#include <scptrs.h>
#include <stptrs.h>
#include <mpcluster.h>
c-----------------------------------------------
c     Local declaration
c-----------------------------------------------
      INT_TYPE k,m,n,ic,sc_p            
      REAL_TYPE scgrep(0:1,nscg)
      REAL_TYPE lrevorsc(0:1,nsc*plev)
      REAL_TYPE fld_file,level_file
      INT_TYPE reounit
c-----------------------------------------------
c     executable statements
c-----------------------------------------------
         do k=1,plev
            call sc_pgather(lrevorsc,scgrep,plev,k)
            level_file=k
            if (rank.eq.root) then
               write(reounit)level_file,fld_file
               call writereal(scgrep,2*nscg,reounit)
            endif
         enddo

      end
