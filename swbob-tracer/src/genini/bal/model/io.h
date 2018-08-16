#ifndef IO_H_
#define IO_H_

! IO error messages

#define OPEN_FAIL      1
#define UNIT_ASGN_FAIL 2
#define HDR_CHECK_FAIL 3
#define MPI_UNIT       99

#include <type.h>

      character*16   caseid
      character*80  inifile
      character*80  rstfile
      INT_TYPE      rstfreq

      common /io/ caseid, inifile, rstfile, rstfreq

      INT_TYPE      zavgfreq      ! sampling frequency of zonal averages

      common /zavg/ zavgfreq

#endif      

