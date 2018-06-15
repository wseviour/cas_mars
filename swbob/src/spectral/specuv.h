#ifndef SPECUV_H_
#define SPECUV_H_

#include <dims.h>
      REAL_TYPE fcor(0:1)         ! coriolis factor spectral coefficients m=1,n=0,1
      REAL_TYPE annp1(0:nn+1)     ! Radius/(n*(n+1))

      common /specuv/fcor,annp1

#endif
