#include <type.h>

       subroutine readreal(fld,len,iunit)
       implicit none

       INT_TYPE len
       REAL_TYPE fld(len)
       INT_TYPE iunit

       INT_TYPE i

       read(iunit)(fld(i),i=1,len)

       end

       subroutine writereal(fld,len,ounit)
       implicit none

       INT_TYPE len
       REAL_TYPE fld(len)
       INT_TYPE ounit

       INT_TYPE i

       write(ounit)(fld(i),i=1,len)

       end

       subroutine writespreal(fld,len,ounit)
       implicit none

       INT_TYPE len
       real*4 fld(len)
       INT_TYPE ounit

       INT_TYPE i

       write(ounit)(fld(i),i=1,len)

       end

       subroutine readint(fld,len,iunit)
       implicit none

       INT_TYPE len
       INT_TYPE fld(len)
       INT_TYPE iunit

       INT_TYPE i

       read(iunit)(fld(i),i=1,len)

       end

       subroutine writeint(fld,len,ounit)
       implicit none

       INT_TYPE len
       INT_TYPE fld(len)
       INT_TYPE ounit

       INT_TYPE i

       write(ounit)(fld(i),i=1,len)

       end
