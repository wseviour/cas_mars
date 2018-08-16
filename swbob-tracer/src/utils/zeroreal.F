      subroutine zeroreal(dest,len)
      implicit none

#include <type.h>

      INT_TYPE  len
      REAL_TYPE dest(len)

      INT_TYPE i

      do i=1,len
         dest(i)=0.0d0
      end do

      end
