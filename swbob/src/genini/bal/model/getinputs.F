      subroutine getinputs()
      implicit none

#include <io.h>
#include <mpcluster.h>
#include <vordist.h>

      namelist /modelnl/ amp,
     $                   n0,
     $                   gamma

      if (rank.eq.root) then
         read (5,modelnl)
      end if

      call broadcastreal(amp,1)
      call broadcastreal(n0,1)
      call broadcastreal(gamma,1)

      call initio()
      call initplanet()      

      end 
