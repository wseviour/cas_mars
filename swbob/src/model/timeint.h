#ifndef TIMEINT_H_
#define TIMEINT_H_

      REAL_TYPE si          ! semi-implicit switch      
      REAL_TYPE tfc         ! Robert time filter coefficient

      parameter(si     = 1.0D0,    ! semi-implicit switch
     $          tfc    = 0.05D0 )   ! time filtering coefficient

      LOGICAL  dofilter            ! switch to control time filtering
      common /timeint/ dofilter

#endif


