      subroutine getinputs()
      implicit none

#include <modelctl.h>
#include <io.h>
#include <mpcluster.h>
#include <fordamp0.h>
#include <hyper.h>

      namelist /modelnl/ ctitle,        ! model run title
     $                   caseid,        ! case id 
     $                   nspday,        ! number of steps per day
     $                   ntimestep,     ! number of timesteps to take
     $                   inifile,       ! initial file name
     $                   rstfreq,       ! restart file frequency
     $                   fldfreq,       ! 2d field output frequency
     $                   zonfreq,        ! zonal output frequency
     $                   delhi,         ! hyperdiffusion order
     $                   dello,         ! hypodiffusion order
     $                   nuhi,          ! hyperdiffusion coefficient
     $                   nulo           ! hypodiffusion coefficient

      if (rank.eq.root) then
         read (5,modelnl)
      end if

      ! scale/initialize planet data

      call broadcastreal(delhi,1)
      call broadcastreal(dello,1)
      call broadcastreal(nuhi,1)
      call broadcastreal(nulo,1)

      call initio()
      call initmodelctl()
      
      end 
