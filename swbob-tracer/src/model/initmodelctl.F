      subroutine initmodelctl()
      implicit none

#include <planet.h>
#include <timestep.h>
#include <timeint.h>
#include <modeltime.h>
#include <io.h>

      ! scale/initialize timestep data

      timestep  = SecPerDay/nspday

      if (ntimestep.lt.0) then
         ntimestep = -nspday*ntimestep
      end if

      call broadcastint(ntimestep,1)
      call broadcastint(nspday,1)
      call broadcastreal(timestep,1)

      write(MPI_UNIT,9000)
 9000 format("")
      write(MPI_UNIT,9010)
 9010 format("---------model control data----------------")
      write(MPI_UNIT,9020)ntimestep
 9020 format(" total number of timesteps   = ",i8)
      write(MPI_UNIT,9030)nspday
 9030 format(" number of timesteps/day     = ",i8)
      write(MPI_UNIT,9040)timestep
 9040 format(" timestep                    = ",f12.6)
      write(MPI_UNIT,9050)tfc
 9050 format(" time filter coefficient     = ",f12.6)
      write(MPI_UNIT,9060)si
 9060 format(" semi implicit switch        = ",f12.6)
      write(MPI_UNIT,9080)
 9080 format("-------------------------------------------")

      end
