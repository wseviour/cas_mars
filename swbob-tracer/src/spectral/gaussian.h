#ifndef GAUSSIAN_H_
#define GAUSSIAN_H_

#include <dims.h>

      REAL_TYPE gw(nlat2,nnlt,jelemg) ! distributed gaussian weights associated
                                      ! to the nonlinear terms
      REAL_TYPE gen_gw(nlat2,jelemg) ! distributed gaussian weights used in gen_analysis
      common /gaussian/gw,gen_gw

#endif
