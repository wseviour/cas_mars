      subroutine timeadvance()
      implicit none

#include <modeltime.h>
#include <timestep.h>
#include <timeint.h>

      REAL_TYPE tiny
      parameter(tiny=1.e-11)

      ! advance the model time

      if (time_current.lt.tiny) then
         time_current = time_current + dt2
      else 
         time_current = time_current + 0.5*dt2
      end if

      if (time_current.lt.(timestep+tiny)) then
         ! initial toggle the filter (after first two timesteps)
         dt2     = 2*dt2
      end if

      if (time_current.lt.(timestep/2+tiny)) then
         ! initial: no filtering first couple of steps
         dofilter= .false.
      else
         ! restart always filters 
         dofilter= .true.
      end if

#ifdef ADAPT
      if (reducedt) then
         timestep=timestep*0.75d0
         dt2=2.d0*timestep
      endif
#endif

      end 
