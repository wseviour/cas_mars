#ifndef DEL8_H_
#define DEL8_H_

#include <planet.h>

      REAL_TYPE d8damp, nudamp
      REAL_TYPE td8fac(0:nn),ud8fac(0:nn)
      REAL_TYPE td8cor(0:nn),ud8cor(0:nn)

      common /del8fac/ td8fac,ud8fac
      common /del8cor/ td8cor,ud8cor

      parameter( d8damp=1.0d0/(0.1d0*SecPerDay) )
      parameter( nudamp=1.517433222D-11 )

#endif


