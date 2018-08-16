
      subroutine unpacksc1(scbuf,sc,m,level)
      implicit none

#include <dims.h>
#include <io.h>

      INT_TYPE  m
      REAL_TYPE scbuf(0:1,m:nn)
      REAL_TYPE sc(0:1,m:nn,plev)

      INT_TYPE level

      INT_TYPE n,ic

      do n=m,nn
         do ic=0,1
            sc(ic,n,level)=scbuf(ic,n)
         enddo
      enddo

      return
      end
