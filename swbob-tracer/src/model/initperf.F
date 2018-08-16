      subroutine initperf()
      implicit none

#include <type.h>
#include <timer.h>

      INT_TYPE k

      do k=1,4
         call timer_clear(k)
      end do

c following to time whole timestep loop:
      call timer_clear(5)

      end
