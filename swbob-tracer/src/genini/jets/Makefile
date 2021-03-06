###
### GNU makefile for genini for SW 
###

ARCH=$(shell uname -s)
TOP=../src
RUNDIR=.
GEN=genx

VPATH=$(TOP)/genini/jets:$(TOP)/include:$(TOP)/model:$(TOP)/putils:$(TOP)/rfft:$(TOP)/spectral:$(TOP)/utils
INCPATH=-I$(TOP)/genini/jets -I$(TOP)/include -I$(TOP)/putils -I$(TOP)/utils -I$(TOP)/model -I$(TOP)/spectral
LOCALPATH = $(HOME)/lib/$(ARCH)

ifeq ($(ARCH),Linux)
   ifeq ($(HOST),localhost)
      COMP      = GNU
      FC        = g77
      CC        = cc
      FBASE     = -ffixed-line-length-none
      MPILIB    =
   endif
   ifeq ($(HOST),lightning)
      COMP      = MPG
      FC        = mpif90
      CC        = mpicc
      FBASE     = -Mextend -Mnosave
      MPILIB    = -lmpich
   endif
   ifeq ($(HOST),copson)
      COMP      = MPG
      FC        = mpif90
      CC        = mpicc
      FBASE     = -Mextend -Mnosave
      MPILIB    = -lmpich
   endif
   ifeq ($(HOST),mhdcluster)
      COMP      = MPG
      FC        = mpif90 -f90=pgf90
      CC        = mpicc -cc=pgcc
      FBASE     = -Mextend -Mnosave
      MPILIB    = -lmpich
   endif
   CPPFLAGS  = $(INCPATH) 
   FOPT      = -O3
   FFLAGS    = $(FBASE) $(FOPT) -D$(COMP)
   CBASE     = 
   COPT      = 
   CFLAGS    = $(CBASE) $(COPT)
   LDFLAGS   = $(FOPT)  
   MPILIB    = #-lmpich
   LOCALLIB  = -L$(LOCALPATH) -ltimer
   LIBS      = $(LOCALLIB) $(MPILIB) 

   ifeq ($(_MPI),1)
      FFLAGS += -DMPI -DLAM_MPI
   endif
   ifeq ($(FD),1)
      FFLAGS += -DFDRAG
   endif
   ifeq ($(FR),1)
      FFLAGS += -DFRLAX
   endif
   ifeq ($(FT),1)
      FFLAGS += -DFTOPO
   endif
   ifeq ($(FH),1)
      FFLAGS += -DFHEAT
   endif
   ifeq ($(FM),1)
      FFLAGS += -DFMECH
   endif
   ifeq ($(FS),1)
      FFLAGS += -DFSPEC
   endif
   ifeq ($(AD),1)
      FFLAGS += -DADAPT
   endif
   ifeq ($(HD),1)
      FFLAGS += -DHDISS
   endif
   ifeq ($(GS),1)
      FFLAGS += -DGAUSS
   endif
endif

ifeq ($(ARCH),IRIX64)
   LOCALPATH = $(HOME)/lib/$(ARCH)
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
endif

ifeq ($(ARCH),HP-UX)
   LOCALPATH = $(HOME)/lib/$(ARCH)
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
  LOCALPATH = $(HOME)/lib/$(ARCH)
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

GEN_OBJ=genini.o initscg.o inituth.o inituthq.o initpert.o

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
