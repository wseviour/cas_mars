#include <type.h>
      INT_TYPE function equalint(i,iref,name)
      implicit none
      
      INT_TYPE  i
      INT_TYPE  iref
      character*(*) name

      if (iref.eq.i) then
         print *,"......",name," comparison PASSED"
         equalint=0
      else
         print *,"ERROR:",name," comparison FAILED"
         print *,"       file val=",i," model val=",iref
         print *,"       This will cause an ABORT!"
         equalint=1
      end if  

      end 

