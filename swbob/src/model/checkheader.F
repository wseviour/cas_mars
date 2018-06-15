c ==================================================================
c
c     Utility function to return index of next available Fortran
c     unit number.
c
c ==================================================================

#include <type.h>
      INT_TYPE function checkheader(iunit)
      implicit none

#include <params.h>
#include <io.h>
#include <modeltime.h>
#include <timestep.h>

      INT_TYPE iunit           ! file unit

      INT_TYPE equalint        
      external equalint

      INT_TYPE equalreal        
      external equalreal
     
c ==================================================================
c
c     Local Variables
c
c ==================================================================

      REAL_TYPE tiny
      parameter(tiny=1.e-11)

      INT_TYPE ierr

c ==================================================================
c
c     Code Body
c
c ==================================================================

      REAL_TYPE   plon_file,plev_file

      read(iunit)plon_file,      
     $           plev_file,
     $           time_start

c new:
      timesteprst=0.d0
      if (time_start.gt.timestep) read(iunit) timesteprst    ! for restart file

      ierr=0
      ierr=ierr+equalint(NINT(plon_file),PLON,'plon')
      ierr=ierr+equalint(NINT(plev_file),PLEV,'plev')


      if (ierr.ne.0) ierr=HDR_CHECK_FAIL

      checkheader=ierr

      end

