#include <type.h>
      INT_TYPE function equalreal(r,rref,name)
      implicit none
      
      REAL_TYPE  r
      REAL_TYPE  rref
      character*(*) name

      REAL_TYPE tiny
      parameter(tiny=1.e-11)

      if (ABS(rref-r).lt.tiny) then
         print *,"......",name," comparison PASSED"
         equalreal=0
      else
         print *,"ERROR:",name," comparison FAILED"
         print *,"       file val=",r," model val=",rref
         equalreal=1
      end if  

      end

