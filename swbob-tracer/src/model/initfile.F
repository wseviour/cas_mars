      subroutine initfile()
      implicit none

#include <mpcluster.h>
#include <io.h>
#include <modeltime.h>
#include <timestep.h>

      INT_TYPE iunit,ierr

      INT_TYPE checkheader
      external checkheader

      INT_TYPE navu
      external navu

      ! open the file

      ierr=0

      if (rank.eq.root) then

         ! assign unit number
         iunit=navu()
         if (iunit.gt.0) then

            open(UNIT   = iunit,
     $           FILE   = inifile,
     $           ACCESS = 'sequential',
     $           STATUS = 'old',
     $           FORM   = 'unformatted',
     $           IOSTAT = ierr)

            if (ierr.eq.0) then        
               ! file opening was good
               print *,'opened startup file ',inifile,' on unit ',iunit,' successfully'
               ! check file parameters in header vs model parameters            
               ! if parameters are bad, set ierr to 1
               ierr=checkheader(iunit)
            else
               ierr=OPEN_FAIL
            endif

         else
            ierr = UNIT_ASGN_FAIL
         endif

      endif

      ! brodcast time_start for restart business

      call broadcastreal(time_start,1)
      call broadcastreal(timesteprst,1)

      ! broadcast ierr to the other processes

      call BroadcastInt(ierr,1)
      
      if (ierr.eq.OPEN_FAIL) then
         call haltmp('ABORT:error reading initial state file')
      else if (ierr.eq.UNIT_ASGN_FAIL) then
         call haltmp('ABORT:cannot assign unit for initial state file')
      else if (ierr.eq.HDR_CHECK_FAIL) then
         call haltmp('ABORT: header check fail for initial state file')
      else
         call startsht(iunit)
      endif

      if (rank.eq.root) then
         close(iunit)
      endif

      end




