#ifndef DECOMP_H_
#define DECOMP_H_

#include <dims.h>

      INT_TYPE mbeg(0:1,0:M_NODE-1) 
      INT_TYPE mend(0:1,0:M_NODE-1)
      INT_TYPE mchmax(0:M_NODE-1)

      common /mdec/ mbeg,mend,mchmax

      INT_TYPE jbeg(0:P_NODE-1),jend(0:P_NODE-1)     
      INT_TYPE jebeg(0:P_NODE-1),jeend(0:P_NODE-1) ! bounds of global latitudinal index for a tile element on a given PE

      common /jdec/ jbeg,jend,jebeg,jeend

#endif

