      subroutine sht()
      implicit none

#include <fields.h>
#include <mpcluster.h>
#include <timer.h>

      call timer_start(1)
      call fourier(sfc,phy,nlt,afc)                  ! sfc -> afc
      call timer_stop(1)

      call timer_start(2)
      call transpf(afc,afcx,abuflen,M_NODE,P_NODE)   !...transpose afc->afcx
      call timer_stop(2)

      call timer_start(3)
      call spectral(afcx,sfcx)                       ! afcx -> sfcx
      call timer_stop(3)

      call timer_start(4)
      call transpb(sfcx,sfc,sbuflen,P_NODE,M_NODE)   !...transpose sfcx->sfc
      call timer_stop(4)

      end
