#include <type.h>
       INT_TYPE function log2(n)
       implicit none

       INT_TYPE n

       INT_TYPE m

       log2=0
       m=n

       do while(m.gt.1)
          m=m/2   
          log2=log2+1
       end do

       end 
