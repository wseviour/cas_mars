
How to run BOB
=============

The current set up only runs on Linux. There are some issues with the 
case insensitivity of the Mac system.

# `libtimer.so`

Before the model can be run, we need to compile `libtimer.so`. 
Navigate to `src/timer`, then delete any compiled `*.so` files. 
Then run:

     make timer_clear.o
     make timer_clock.o
     make timer_start.o
     make timer_stop.o
     make timer_time.o
     make libtimer.a


phibar = gH, defined as phibar = (Ld*2*Omega)**2 in `planet.h`. 