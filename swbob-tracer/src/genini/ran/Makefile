###
### GNU makefile for genini for SW 
###

ARCH=$(shell uname -s)
HOST=$(shell hostname -s)
LOCALPATH = $(HOME)/lib/$(ARCH)
TOP=../src
RUNDIR=.
GEN=genx

VPATH=$(TOP)/genini/zz:$(TOP)/include:$(TOP)/model:$(TOP)/putils:$(TOP)/rfft:$(TOP)/spectral:$(TOP)/utils
INCPATH=-I$(TOP)/genini/zz -I$(TOP)/include -I$(TOP)/putils -I$(TOP)/utils -I$(TOP)/model -I$(TOP)/spectral

#########################################
# Tested/supported architectures
# (22 February 2003, Rich Loft)
#
# 1) Linux/Portland Group compilers + MPI (optional, call make _MPI=1)
# 2) IBM AIX/xlf90 compilers        + MPI (optional, call make _MPI=1)
# 3) OSX (Darwin)/Absoft compilers
#########################################

ifeq ($(ARCH),Linux)
   ifeq ($(HOST),flare)
      FC        = pg77
      CC        = cc
      MPILIB    =
   endif
   ifeq ($(HOST),lightning)
      FC        = mpif77
      CC        = mpicc
      MPILIB    = -lmpich
   endif
   CPPFLAGS  = $(INCPATH) 
   FBASE     = $(ENVFBASE) -Mextend
   FOPT      = -O2 -Munroll=c:8,n:8 -Mdalign 
   FFLAGS    = $(FBASE) $(FOPT)
   CBASE     = #-DF77_UNDERSCORE
   COPT      = 
   CFLAGS    = $(CBASE) $(COPT)
   LDFLAGS   = $(FOPT)  
   LOCALLIB  = -L$(LOCALPATH) -ltimer
   LIBS      = $(LOCALLIB) $(MPILIB) 

   ifeq ($(_MPI),1)
      FFLAGS += -DMPI
   endif
endif

ifeq ($(ARCH),AIX)
   FC        = mpxlf90_r
   CC        = mpcc_r
   CPPFLAGS  = $(INCPATH)
   FBASE     = -qfixed=132
   FOPT      = -O3 -qarch=pwr4 -qtune=pwr4 -qcache=auto
   FFLAGS    = $(FBASE) $(FOPT)
   CBASE     = 
   COPT      = -O3 -qarch=pwr4 -qtune=pwr4 -qcache=auto
   CFLAGS    = $(CBASE) $(COPT) 
   LDFLAGS   = $(FOPT)
   MPILIB    = 
   LOCALLIB  = 
   LIBS      = $(LOCALLIB) $(MPILIB)

   ifeq ($(_MPI),1)
      FFLAGS += -DMPI
   endif

endif

ifeq ($(ARCH),Darwin)
   FC        = f90
   CC        = gcc
   CPPFLAGS  = $(BOBPATH)
   FBASE     = -D$(DELORDER) -D$(NU) -f fixed -W 132
   FOPT      = -O
   FFLAGS    = $(FBASE) $(FOPT)
   CBASE     =
   COPT      =
   CFLAGS    = $(CBASE) $(COPT)
   LDFLAGS   = $(FOPT)
   LOCALLIB  = -L$(LOCALPATH) -ltimer
   LIBS      = $(LOCALLIB)
endif

#########################################
# Untested architectures, which may have
# worked at one time.
# 1) IRIX
# 2) HP/UX
# 3) SunOS
#########################################

ifeq ($(ARCH),IRIX64)
   FC        = f90
   CC        = cc
   CPPFLAGS  = $(INCPATH)
   FBASE     = -extend_source -cpp
   FOPT      = -O3 -64 -mips4 -static
   FFLAGS    = $(FBASE) $(FOPT)
   CBASE     = -DF77_UNDERSCORE 
   COPT      = -O3 -64 -mips4 -static
   CFLAGS    = $(CBASE) $(COPT) 
   LDFLAGS   = $(FOPT)
   MPILIB    = -L$(LOCALPATH) -ltimer
   LOCALLIB  = -lmpi
   LIBS      = $(LOCALLIB) $(MPILIB)
endif


ifeq ($(ARCH),HP-UX)
   FC        = fort77
   CC        = /opt/ansic/bin/cc
   CPPFLAGS  = $(INCPATH) -I/opt/mpi/include 
   FBASE     = +Oinfo +es
   FOPT      = +O2 +Onofltacc +DS2.0a +DA2.0N +Odataprefetch +U77 
   FFLAGS    = $(FBASE) $(FOPT)
   CBASE     =  
   COPT      = -O
   LDFLAGS   = $(FOPT)
   LOCALLIB  = -L$(LOCALPATH) -ltimer
   MPILIB    = -L/opt/mpi/lib/pa1.1 -lmpi /lib/libail.sl -lcnx_syscall 
   LIBS      = $(MPILIB) $(LOCALLIB)
endif

ifeq ($(ARCH),SunOS)
  FC       = mpf77
  CC       = mpcc
  CPPFLAGS  = $(INCPATH) -I/opt/SUNWhpc/HPC3.1/include
  FBASE    = -$(ENVFBASE) -e -dalign
  FPREC    = 
  FOPT     = -fast -xarch=v8
  FFLAGS   = $(FBASE) $(FOPT)
  CBASE    = -I/opt/SUNWhpc/HPC3.1/include
  COPT     = -xO2 -xarch=v8
  CFLAGS   = $(CBASE) $(COPT)
  LDFLAGS  = $(FOPT)
  LOCALLIB = -L$(LOCALPATH) -ltimer
  MPILIB   = -lmpi
  LIBS     = $(MPILIB) $(LOCALLIB)
endif

GLOBAL_INC=dims.h params.h type.h mpcluster.h

GEN_OBJ=genini.o initscg.o inituth.o initpert.o

GEN_INC=planet.h coord.h mpcluster.h io.h dims.h stptrs.h scptrs.h decomp.h \
physical.h constants.h 

PUTILS_OBJ=BroadcastInt.o BroadcastReal.o gsum.o haltmp.o initmp.o recvreal.o \
sendreal.o transpb.o transpf.o ceil2.o error_dup.o pair.o synch.o

PUTILS_INC=multinode.h

UTILS_OBJ=fldio.o lenchr.o ljust.o navu.o copyreal.o zeroreal.o boundint.o \
equalint.o equalreal.o ave.o axpy.o scal.o

UTILS_INC=

RFFT_OBJ=radb2.o radb3.o radb4.o radb5.o radbg.o radf2.o radf3.o radf4.o radf5.o \
radfg.o rfftb.o rfftb1.o rfftf.o rfftf1.o rffti.o rffti1.o

SPECTRAL_OBJ=decomp.o dnlfk.o dnlft.o \
genalp.o initalpseed.o initcoord.o initgenalp.o initmderiv.o initscptrs.o initsht.o \
initspecuv.o jdecomp.o loadalp.o mdecomp.o \
testalp.o testrfft.o gen_analysis.o gen_symf.o physical_to_spectral.o \
sc_pgather.o sc_xfer_layer.o 

SPECTRAL_INC=physical.h alp.h decomp.h fields.h fft.h gaussian.h mderiv.h specuv.h \
spectral.h 

MODEL_OBJ=gridfile.o 

MODEL_INC=coord.h planet.h dims.h io.h

BOB_OBJ=$(GEN_OBJ) $(MODEL_OBJ) $(PUTILS_OBJ) $(RFFT_OBJ) $(SPECTRAL_OBJ) $(UTILS_OBJ)

default: $(RUNDIR)/$(GEN)

$(BOB_OBJ): $(GLOBAL_INC)
$(GEN_OBJ): $(GEN_INC)
$(UTILS_OBJ): $(UTILS_INC)
$(PUTILS_OBJ): $(PUTILS_INC)
$(SPECTRAL_OBJ): $(SPECTRAL_INC)

.SUFFIXES:
.SUFFIXES: .F .o .c

$(RUNDIR)/$(GEN): $(BOB_OBJ)
	$(FC) -o $@ $(BOB_OBJ) $(LDFLAGS) $(LIBS)

.F.o:
	$(FC) -c $(CPPFLAGS) $(FFLAGS) $<

.c.o:
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $<

clean:
	rm -f *.o $(RUNDIR)/$(GEN)
