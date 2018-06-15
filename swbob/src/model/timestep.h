#ifndef TIMESTEP_H_
#define TIMESTEP_H_

#include <type.h>


      REAL_TYPE timestep,timesteprst     ! input timestep
      REAL_TYPE dt2          ! 2*(current timestep)

      common /tstep/ timestep,timesteprst,dt2

      LOGICAL reducedt
      common /tsadapt/ reducedt

#endif
