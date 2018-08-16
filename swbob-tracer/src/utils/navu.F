c ==================================================================
c
c     Utility function to return index of next available Fortran
c     unit number.
c
c ==================================================================

#include <type.h>

      INT_TYPE function navu()
      implicit none

c ==================================================================
c
c     Local Variables
c
c ==================================================================

      logical isopen          ! logical variable tested for by 
                              ! inquire function

      INT_TYPE iunit           ! Test unit

c ==================================================================
c
c     Code Body
c
c ==================================================================

      iunit=6
      isopen=.true.

      do while( isopen .and. iunit.lt.99)
         iunit=iunit+1
         inquire(unit=iunit,opened=isopen)
      end do

      if (isopen .and. iunit.eq.99) then
         navu=-1
      else
         navu=iunit
      end if

      end

