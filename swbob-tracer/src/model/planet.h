#ifndef PLANET_H_
#define PLANET_H_

#include <type.h>
#include <constants.h>

      REAL_TYPE Radius
      REAL_TYPE SecPerDay
      REAL_TYPE Omega 
      REAL_TYPE Ldeform
      REAL_TYPE PhiBar

      parameter(Radius      = 3396000.0,
     $          SecPerDay   = 88774.0,
     $          Ldeform     = 1774137.76411825273008285858,
     $          Omega       = .00007077731438461245)

      parameter(PhiBar=(2.d0*Omega*Ldeform)**2)

#endif
