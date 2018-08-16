      subroutine unpacksc(scbuf,sc,m,fld_p,level)
      implicit none

#include <dims.h>
#include <io.h>

      INT_TYPE  m
      REAL_TYPE scbuf(0:1,m:nn)
      REAL_TYPE sc(0:1,nsp,m:nn,plev)

      INT_TYPE fld_p,level

      INT_TYPE n,ic

      do n=m,nn
         do ic=0,1
            sc(ic,fld_p,n,level)=scbuf(ic,n)
         end do
      end do

      end
