      subroutine perfstats()
      implicit none

#include <dims.h>
c      INT_TYPE plon,plat,plev
#include <timer.h>

#define REAL_BYTES 8

      real*8 legf_flops
      real*8 legb_flops

      real*8 xpf_bytes
      real*8 xpb_bytes

      real*8 fftf_flops
      real*8 fftb_flops

      real*8 time_fourier
      real*8 time_spectral
      real*8 time_xpf
      real*8 time_xpb

      real*8 time_timeloop

      INT_TYPE log2
      external log2

      time_fourier = timer_time(1)/timer_calls(1)
      time_xpf     = timer_time(2)/timer_calls(2)
      time_spectral= timer_time(3)/timer_calls(3)
      time_xpb     = timer_time(4)/timer_calls(4)

      time_timeloop= timer_time(5)

      fftf_flops=5*log2(plon/2)*(plon/2)*(plat*plev*nnlt)*(1.e-6)
      fftb_flops=5*log2(plon/2)*(plon/2)*(plat*plev*nsyn)*(1.e-6)

      legf_flops=plat*plev*nnlt*(mm+1)*(mm+4)*1.e-6
      legb_flops=plat*plev*(3*(mm+1)*(mm+2)+2*(mm+1)*(mm+4))*1.e-6

      xpf_bytes=(abuflen*REAL_BYTES)*P_NODE
      xpb_bytes=(sbuflen*REAL_BYTES)*P_NODE

      print *,"==================================================="
      print *,"FFT MFlop count    = ",fftf_flops+fftb_flops
      print *,"FFT elapse time    = ",time_fourier
      print *,"FFT MFlop/sec rate = ",(fftf_flops+fftb_flops)/time_fourier
      print *

      print *,"Transpose Fwd msg size (bytes) = ",xpf_bytes
      print *,"Transpose Fwd elapse time      = ",time_xpf
      print *,"Transpose Fwd MB/sec rate      = ",(xpf_bytes*1.e-6)/time_xpf
      print *

      print *,"Spectral Space MFlop count = ",legf_flops+legb_flops
      print *,"Spectral Space elapse time = ",time_spectral
      print *,"Spectral Space MFlop rate  = ",(legf_flops+legb_flops)/time_spectral
      print *

      print *,"Transpose Inv msg size (bytes) = ",xpb_bytes
      print *,"Transpose Inv elapse time      = ",time_xpb
      print *,"Transpose Inv MB/sec rate      = ",(xpb_bytes*1.e-6)/time_xpb
      print *
      print *,"==================================================="
      print *,"Total time/timestep=",time_fourier+time_xpf+time_spectral+time_xpb
      print *,"==================================================="

      print *
      print *,"Total time                     = ",time_timeloop
      print *

      end
