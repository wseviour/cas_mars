#include <type.h>
      INT_TYPE function boundint(i,imin,imax,name)
      implicit none

      INT_TYPE i,imin,imax
      character*(*) name

      if (i.lt.imin.or.i.gt.imax) then
         print *,name," out of bounds"
         print *,name,imin," <= ",name," <= ",imax
         print *,"is not satisfied for ",name,"=",i
         boundint=1
      else
c         write(6,10)name,i
c 10      format("......",A8,"=",i3," bounds check PASSED")
         boundint=0
      end if

      end 




