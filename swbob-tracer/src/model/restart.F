      subroutine restart()
      implicit none

#include <type.h>
#include <io.h>
#include <modeltime.h>

      if (itime.gt.0) then

         if (mod(itime,rstfreq).eq.0) then
            print *," output of restart file at time=",time_current
         end if 

      end if

      end 
