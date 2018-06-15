#ifndef FFT_H_
#define FFT_H_

#include <dims.h>

      INT_TYPE  numlons                         ! number of different longitude grids on a hemisphere
      INT_TYPE  fftptr(jelemd)                  ! points to appropriate fft
      REAL_TYPE fftwts(2*plon+15,MAX_NUM_LONS)  ! fftwts needed for model (tabulation of trigo fcts)
      INT_TYPE  lonmax(jelemd)                  ! local max longitudes

      common /wts/ fftwts,fftptr,lonmax,numlons

#endif
