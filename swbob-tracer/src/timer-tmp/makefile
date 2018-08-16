#
# Makefile for timer library
#

ARCH=$(shell uname -s)
LOCALPATH = $(HOME)/lib/$(ARCH)

ifeq ($(ARCH),Linux)
  AR       = ar ruv
#  FC       = pgf77
  FC       = g77
  CC       = gcc
#  FBASE    = -Mextend -I.
  FBASE    = -ffixed-line-length-none  -fno-second-underscore -I.
  FOPT     = -O2
  CPP      = /lib/cpp 
  CPPFLAG  = -DLinux -DLWR_CASE_USCORE -P
  CBASE    =  -DLinux -DLWR_CASE_USCORE -P
  COPT     = 
endif

ifeq ($(ARCH),AIX)
   AR        = ar ruv
   FC        = xlf90
   CC        = cc
   FBASE     = -qfixed=132
   FOPT      = -O3 -qarch=pwr4 -qtune=pwr4 -qcache=auto -Q+Anal
   FPREC     =
   CPP       = /lib/cpp
   CPPFLAG   =
   CBASE     = -DAIX
   COPT      = -O3 -qarch=pwr4 -qtune=pwr4 -qcache=auto -Q+Anal
endif

ifeq ($(ARCH),IRIX64)
  AR       = ar ruv
  FC       = f77
  CC       = cc
  FBASE    =
  FOPT     = -O3 -mips4 -64
  FPREC    = -r8 -i4 
  CPP      = /lib/cpp 
  CPPFLAG  = -D$(ARCH) -P
  CBASE    = -64 $(CPPFLAG)
  COPT     =  
endif

ifeq ($(ARCH),HP-UX)
   AR        = ar ruv
   FC        = fort77
   CC        = /opt/ansic/bin/cc
   CPP       = /opt/langtools/lbin/cpp
   CPPFLAG   = -DHP
   FBASE     = +Oinfo +es
   FOPT      = +O3 +DS2.0a +DA2.0N +Odataprefetch +Onolimit +Olibcalls \
	       +Ono_fail_safe +U77 +Onoinline 
   CBASE     = $(CPPFLAG) 
   COPT      = -O
endif

ifeq ($(ARCH),SunOS)
  AR       = ar ruv
  FC       = f77
  CC       = cc
  CPP      = cpp 
  CPPFLAG  = -D$(ARCH) 
  FBASE    = 
  FPREC    = 
  FOPT     = 
  CBASE    = $(CPPFLAG)
  COPT     = 
endif

FFLAG  = $(FBASE) $(FPREC) $(FOPT)
CFLAG  = $(CBASE)          $(COPT) 

LIBTIMER = libtimer.a
HDRS = timer.h

FSOURCES= timer_clear.F timer_start.F timer_stop.F timer_time.F\
          timer_calls.F

CSOURCES= timer_clock.c 

FOBJECTS = $(patsubst %.F, %.o, $(FSOURCES))
COBJECTS = $(patsubst %.c, %.o, $(CSOURCES))

libtimer.a: $(FOBJECTS) $(COBJECTS)
	$(AR) $(LIBTIMER) $(FOBJECTS) $(COBJECTS)
	mv libtimer.a $(LOCALPATH)/libtimer.a

%.o: %.F
	$(CPP) $(CPPFLAG) -P  $*.F >! $*.f
	$(FC) $(FFLAG) -c $*.f

%.o: %.c
	$(CC) $(CFLAG) -c $*.c

clean:
	rm -rf *.o *.f ! 


clobber:
	rm -rf *.o *.f ! $(LOCALPATH)/libtimer.a








