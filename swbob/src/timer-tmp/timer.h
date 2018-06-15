c
c            timer.h
c
c     timer routines - errors
c        0 = no error
c       -1 = timer_id is out of bounds
c
c**********************************************************************

      integer MAX_TIMER,MTP1
      parameter (MAX_TIMER = 64,MTP1=MAX_TIMER+1)

c**********************************************************************
c
c     btimers = the common area for timer information
c
      common /timer_com/ timer_assign,     timer_starts, 
     1                   timer_accum_time, timer_cur_time

c**********************************************************************
c
c     Keeps track of whether timer has been initialized already
c
      integer timer_assign(0:MAX_TIMER)

c**********************************************************************
c
c     timer_starts = number of times this clock has been started
c
      integer timer_starts(0:MAX_TIMER)

c**********************************************************************
c
c     timer_accum_time = time accumulated in each timer...
c
      real*8 timer_accum_time(0:MAX_TIMER)

c**********************************************************************
c
c     timer_cur_time = current time reading for each timer...
c
      real*8 timer_cur_time(0:MAX_TIMER)


c**********************************************************************


