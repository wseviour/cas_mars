#include <type.h>
      REAL_TYPE function ave(fld,len)
      implicit none

#include <dims.h>

      INT_TYPE  len
      REAL_TYPE fld(len)

      INT_TYPE i

      ave=0
      do i=1,len
         ave=ave+fld(i)
      end do
 
      ave=ave/len

      end 
