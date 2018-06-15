      subroutine decomp(beg,len,npe)
      implicit none

#include <type.h>

      INT_TYPE npe
      INT_TYPE beg(0:npe)
      INT_TYPE len

      INT_TYPE n

      beg(0)=1
      do n=1,npe-1
         if (n.le.mod(len,npe)) then
            beg(n)=beg(n-1)+(len-1)/npe+1
         else
            beg(n)=beg(n-1)+len/npe
         end if
      end do
      beg(npe)=len+1

      end
