###
### GNU makefile for genini for SW 
###

ARCH=$(shell uname -s)
TOP=../src
RUNDIR=.
GEN=genx

VPATH=$(TOP)/genini/zz:$(TOP)/include:$(TOP)/model:$(TOP)/putils:$(TOP)/rfft:$(TOP)/spectral:$(TOP)/utils
INCPATH=-I$(TOP)/genini/zz -I$(TOP)/include -I$(TOP)/putils -I$(TOP)/utils -I$(TOP)/model -I$(TOP)/spectral

ifeq ($(ARCH),Linux)
   ifeq ($(HOST),modon)    # and newer storms
      COMP      = GNU
      FC        = gfortran
      CC        = cc
      FBASE     = -ffixed-line-length-none -fno-second-underscore # -frecord-marker=4
      MPILIB    =
      ifeq ($(_MPI),1)
         FC        = /usr/lib64/openmpi/1.4-gcc/bin/mpif90
         MPILIB    = -lmpich
      endif
   endif
   ifeq ($(HOST),wha)     # and older storms
      COMP      = GNU
      FC        = gfortran
      CC        = cc
      FBASE     = -ffixed-line-length-none -fno-second-underscore
      MPILIB    =
   endif
   LOCALPATH = $(HOME)/lib/$(ARCH)
   CPPFLAGS  = $(INCPATH) 
   FOPT      = -O2
   FFLAGS    = $(FBASE) $(FOPT)
   CBASE     = # -DF77_UNDERSCORE
   COPT      = 
   CFLAGS    = $(CBASE) $(COPT)
   LDFLAGS   = $(FOPT)  
   MPILIB    = # -lmpich
   LOCALLIB  = -L$(LOCALPATH) -ltimer
   LIBS      = $(LOCALLIB) $(MPILIB) 
endif

GLOBAL_INC=dims.h params.h type.h mpcluster.h

GEN_OBJ=genini.o initscg.o initvorth.o

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
