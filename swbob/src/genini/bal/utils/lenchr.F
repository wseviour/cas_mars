#include <type.h>
      INT_TYPE function lenchr(chrstg)
      implicit none
C-----------------------------------------------------------------------
C
C Return position of right-most non-blank, non-null character
C in chrstg.
C
C---------------------------Code history--------------------------------
C
C Original version:  L. Bath, April 1992
C Standardized:      L. Bath, June 1992
C                    J. Rosinski April 1994
C                    T. Acker, March 1996
C
C-----------------------------------------------------------------------

C------------------------------Arguments--------------------------------
C
C Input arguments
C
C-----------------------------------------------------------------------

      character*(*) chrstg       !  Input character string
C
C--------------------------Local Variables------------------------------
C
      integer l     ! loop counter
C
C-----------------------------------------------------------------------
C
      lenchr = 0
      do l=len(chrstg),1,-1
C
C Find right-most non-blank character in string
C
         if (chrstg(l:l).ne.' ' .and. chrstg(l:l).ne.char(0)) then
            lenchr = l
            goto 10
         end if
      end do

   10 return

      end


