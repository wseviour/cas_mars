#ifndef SPECTRAL_H_
#define SPECTRAL_H_

#include <dims.h>

      REAL_TYPE nltsc(2*(plev*nnlt)*nvsc)    ! spectral space non linear terms
      REAL_TYPE sc(2*(plev*nsp)*nsc)         ! time level n state variable spectral coefficients
      REAL_TYPE scm1(2*(plev*nsp)*nsc)       ! time level n-1 state variable spectral coefficients
      REAL_TYPE fthesc(2*plev*nsc)       
      REAL_TYPE fthpsc(2*plev*nsc)       
      REAL_TYPE ftoposc(2*plev*nsc)
      REAL_TYPE ftopoktsc(2*plev*nsc)
      REAL_TYPE fvorsc(2*plev*nsc)       

      common /spec/ nltsc,sc,scm1,fthesc,fthpsc,ftoposc,ftopoktsc,fvorsc

#endif