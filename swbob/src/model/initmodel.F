      subroutine initmodel()
      implicit none

#include <planet.h>
#include <coord.h>

      call getinputs()
      call gridfile()
      call initsht(colatitude,
     $             gausswt,
     $             plonmax,
     $             Alpha_coord,
     $             Radius,
     $             Omega)

      call initfile()

      end
