      subroutine initplanet()
      implicit none

#include <planet.h>
#include <io.h>

      write(MPI_UNIT,9000)
 9000 format("")
      write(MPI_UNIT,9010)
 9010 format("-----------planetary data-----------------")
      write(MPI_UNIT,9020)Radius
 9020 format(" Radius (m)          = ",e11.4)
      write(MPI_UNIT,9030)Gravity
 9030 format(" Gravity (m/s^2)     = ",e11.4)
      write(MPI_UNIT,9040)LengthOfDay
 9040 format(" Length of Day (hrs) = ",e11.4)
      write(MPI_UNIT,9050)
 9050 format("-------------------------------------------")

      end

 
