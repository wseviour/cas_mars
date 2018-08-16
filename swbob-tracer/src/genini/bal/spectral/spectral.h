#ifndef SPECTRAL_H_
#define SPECTRAL_H_

#include <dims.h>

      REAL_TYPE nltsc(2*(plev*nnlt)*nvsc)    ! spectral space non linear terms
      REAL_TYPE sc(2*(plev*nsp)*nsc)         ! time level n state variable spectral coefficients
      REAL_TYPE scm1(2*(plev*nsp)*nsc)       ! time level n-1 state variable spectral coefficients
      REAL_TYPE gpstar(2*nsc)                ! lower boundary geopotential
      REAL_TYPE forvorsc(2*plev*nsc)       
      REAL_TYPE fordivsc(2*plev*nsc)       
      REAL_TYPE forthesc(2*plev*nsc)       
      REAL_TYPE foreqsc(2*plev*nsc)       

      common /spec/ nltsc,sc,scm1,gpstar,forvorsc,fordivsc,forthesc,foreqsc

#endif
