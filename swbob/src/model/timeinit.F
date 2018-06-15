#include <type.h>
      ! returns timestep to begin on, sets up model time 
      ! and time integration scheme.

      INT_TYPE function timeinit()
      implicit none

#include <modeltime.h>
#include <timestep.h>
#include <timeint.h>
#include <io.h>
#include <planet.h>

      REAL_TYPE tiny
      parameter(tiny=1.e-11)

      INT_TYPE iotoggle
      common /rstrt0/ iotoggle

      iotoggle=0  ! for restart (see spectral.F)

      ! set times in modeltime.h

      time_current=time_start
      time_advance=timestep*ntimestep
      time_stop   =time_start+time_advance

      time_out    =float(fldfreq)*SecPerDay/float(nspday)
      time_zon    =float(zonfreq)*SecPerDay/float(nspday)

c new: change timestep to value read in initfile (checkheader) if restart:
      if (timesteprst.gt.0.d0) timestep=timesteprst


      ! set dt2 and filter flag before stepping begins
      
      if (time_start.lt.tiny) then  ! initial file
         dt2         = timestep/2
         dofilter    = .false.
         timeinit    = 0
      else                          ! restart file
         dt2         = 2*timestep
         dofilter    = .true.
         timeinit    = 1
      end if

      end
