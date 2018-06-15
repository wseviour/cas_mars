      subroutine blduv3d(sc,
     $                   nltsc,
     $                   eps,
     $                   m)
      implicit none

#include <dims.h>

      INT_TYPE  m
      REAL_TYPE nltsc(0:1,nnlt,m:nn+1,plev)
      REAL_TYPE sc(0:1,nsp,m:nn,plev)
      REAL_TYPE eps(m:nn+1)

      INT_TYPE k

      do k=1,plev
         call blduv(sc(0,1,m,k),
     $              nltsc(0,1,m,k),
     $              eps,
     $              m)
      end do

      end
