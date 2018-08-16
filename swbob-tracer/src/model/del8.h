#ifndef DEL8_H_
#define DEL8_H_

#include <planet.h>

      REAL_TYPE d8damp, nudamp
      REAL_TYPE td8fac(0:nn),ud8fac(0:nn)
      REAL_TYPE td8cor(0:nn),ud8cor(0:nn)

      common /del8fac/ td8fac,ud8fac
      common /del8cor/ td8cor,ud8cor

      parameter( d8damp=24.0d0/(SecPerDay/10.d0) )
      parameter( nudamp=0.d0 )

#endif


