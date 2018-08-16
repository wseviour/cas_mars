#ifndef ALP_H_
#define ALP_H_

#include <dims.h>

      REAL_TYPE Pseed(nlat2,0:nn+1,jelemg,0:1,0:1)
      REAL_TYPE genp(3,nvsc) ! coef of Swarztrauber's reccurence relation 
      REAL_TYPE P(nlat2,0:nn+1,0:1,jelemg)

      common /alp/ Pseed, genp, P

#endif
