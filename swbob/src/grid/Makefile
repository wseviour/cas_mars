#
# Makefile to generate BOB's grid
#
ARCH=$(shell uname -s)
TOP=../src
BLDIR=.
EXE = gridx
VPATH=$(TOP)/grid:$(TOP)/include:$(TOP)/utils
GRIDPATH=-I$(TOP)/include

# If using a Mac, set ARCH to Linux - wjms
ifeq ($(ARCH),Darwin)
	ARCH=Linux
endif

ifeq ($(ARCH),Linux)
   ifeq ($(HOST),typhoon)
      FC        = gfortran
      CC        = cc
      FBASE     = -ffixed-line-length-none
   endif
   ifeq ($(HOST),lightning)
      FC        = pgf77
      CC        = pgcc
      FBASE     = -Mextend -r8
   endif
   ifeq ($(HOST),copson)
      FC        = pgf77
      CC        = pgcc
      FBASE     = -Mextend -r8
   endif
   CPPFLAGS  = $(GRIDPATH)
   FOPT      =
   FFLAGS    = $(FBASE) $(FOPT)
   CBASE     = -I/usr/local/lam/h
   COPT      =
   CFLAGS    = $(CBASE) $(COPT)
   LDFLAGS   = $(FOPT)
   MPILIB    =
   LOCALLIB  =
   LIBS      = $(LOCALLIB) $(MPILIB)
endif

ifeq ($(ARCH),IRIX64)
   FC        = f90
   CC        = cc
   CPPFLAGS  = $(GRIDPATH)
   FBASE     = -extend_source -cpp
   FOPT      = -O3 -64 -static
   FFLAGS    = $(FBASE) $(FOPT)
   CBASE     =
   COPT      = -O3 -64 -static
   CFLAGS    = $(CBASE) $(COPT)
   LDFLAGS   = $(FOPT)
   MPILIB    =
   LOCALLIB  =
   LIBS      = $(LOCALLIB) $(MPILIB)
endif

ifeq ($(ARCH),AIX)
   FC        = xlf90
   CC        = xlc
   CPPFLAGS  = $(GRIDPATH)
   FBASE     = -qfixed=132
   FOPT      = -O3 -qarch=pwr4 -qtune=pwr4 -qcache=auto -Q+Anal
   FFLAGS    = $(FBASE) $(FOPT)
   CBASE     =
   COPT      = -O3 -qarch=pwr4 -qtune=pwr4 -qcache=auto -Q+Anal
   CFLAGS    = $(CBASE) $(COPT)
   LDFLAGS   = $(FOPT)
   MPILIB    =
   LOCALLIB  =
   LIBS      = $(LOCALLIB) $(MPILIB)
endif

ifeq ($(ARCH),HP-UX)
   FC        = fort77
   CC        = /opt/ansic/bin/cc
   CPPFLAGS  = $(GRIDPATH) -I/opt/mpi/include
   FBASE     = +Oinfo +es
   FOPT      = +O3 +DS2.0a +DA2.0N +Odataprefetch +Onolimit +Olibcalls \
	       +Ono_fail_safe +U77 +Onoinline
   FFLAGS    = $(FBASE) $(FOPT)
   CBASE     =
   COPT      = -O
   CFLAGS    = $(CBASE) $(COPT)
   LDFLAGS   = $(FOPT)
   LOCALLIB  =
   MPILIB    =
   LIBS      = $(MPILIB) $(LOCALLIB)
endif

ifeq ($(ARCH),SunOS)
  FC       = f77
  CC       = cc
  CPPFLAGS  = $(GRIDPATH)
  FBASE    = -e -ftrap=%all,no%inexact -C -dalign -fsimple=0 -g
  FPREC    =
  FOPT     = -fast -xarch=v8
  FFLAGS   = $(FBASE) $(FOPT)
  CBASE    =
  COPT     = -O2 -xarch=v8
  CFLAGS   = $(CBASE) $(COPT)
  LDFLAGS  = $(FOPT)
  LOCALLIB =
  MPILIB   =
  LIBS     = $(MPILIB) $(LOCALLIB)
endif

GLOBAL_INC=dims.h type.h constants.h params.h mpcluster.h

MODEL_OBJ=GenGrid.o Pgaqd.o

UTILS_OBJ=fldio.o lenchr.o ljust.o navu.o copyreal.o zeroreal.o boundint.o \
equalint.o equalreal.o ave.o axpy.o scal.o

GRID_OBJ=$(MODEL_OBJ) $(UTILS_OBJ)

default: $(BLDIR)/$(EXE)

$(GRID_OBJ): $(GLOBAL_INC)
$(UTILS_OBJ): $(UTILS_INC)

.SUFFIXES:
.SUFFIXES: .F .o .c

$(BLDIR)/$(EXE): $(GRID_OBJ)
	$(FC) -o $@ $(GRID_OBJ) $(LDFLAGS) $(LIBS)

.F.o:
	$(FC) -c $(CPPFLAGS) $(FFLAGS) $<

.c.o:
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $<

clean:
	-rm -f *~ $(O_DIR)*.o $(BLDIR)/$(EXE)

clobber:
	-rm -f *~ $(O_DIR)*
