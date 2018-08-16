      subroutine copyreal(src,dest,len)
      implicit none

#include <type.h>

      INT_TYPE  len
      REAL_TYPE src(len)
      REAL_TYPE dest(len)

      INT_TYPE i

      do i=1,len
         dest(i)=src(i)
      end do

      end
