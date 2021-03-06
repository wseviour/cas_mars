###
### GNU makefile for genini for SW 
###

ARCH=$(shell uname -s)
HOST=$(shell hostname -s)
LOCALPATH = $(HOME)/lib/$(ARCH)
TOP=../src/genini/bal
RUNDIR=.
GEN=genx

VPATH=$(TOP)/include:$(TOP)/model:$(TOP)/putils:$(TOP)/rfft:$(TOP)/spectral:$(TOP)/utils:$(TOP)/diagnostics
INCPATH=-I$(TOP)/include -I$(TOP)/model -I$(TOP)/putils -I$(TOP)/rfft -I$(TOP)/spectral

ifeq ($(ARCH),Linux)
#   ifeq ($(HOST),flare)
   ifeq ($(HOST),localhost)
      FC        = pg90
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
   CBASE     = 
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
   COPT      = -O4 -qarch=pwr4 -qtune=pwr4 -qcache=auto
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
   FC        = f77
   CC        = cc
   CPPFLAGS  = $(INCPATH)
   FBASE     = -extend_source -cpp -DSWE_TC_1
   FOPT      = -O3 -64 -mips4 -static
   FFLAGS    = $(FBASE) $(FOPT)
   CBASE     = 
   COPT      = -O3 -64 -mips4 -static
   CFLAGS    = $(CBASE) $(COPT) 
   LDFLAGS   = $(FOPT)
   MPILIB    = -lmpi
   LOCALLIB  = -L$(LOCALPATH) -ltimer
   LIBS      = $(LOCALLIB) $(MPILIB)
endif

ifeq ($(ARCH),HP-UX)
   FC        = fort77
   CC        = /opt/ansic/bin/cc
   CPPFLAGS  = $(INCPATH) -I/opt/mpi/include 
   FBASE     = +Oinfo +es -DSWE_TC_1
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
  FC        = mpf77
  CC        = mpcc
  CPPFLAGS  = $(INCPATH) -I/opt/SUNWhpc/HPC3.1/include
  FBASE     = -$(ENVFBASE) -e  -dalign 
  FPREC     = 
  FOPT      = -fast -xarch=v8
  FFLAGS    = $(FBASE) $(FOPT)
  CBASE     = -I/opt/SUNWhpc/HPC3.0/include
  COPT      =  -fast -xarch=v8
  CFLAGS    = $(CBASE) $(COPT)
  LDFLAGS   = $(FOPT)
  LOCALLIB  = -L$(LOCALPATH) -ltimer
  MPILIB    = -lmpi
  LIBS      = $(MPILIB) $(LOCALLIB)
endif


GLOBAL_INC=dims.h nltptrs.h type.h constants.h fields.h\
params.h stptrs.h mpcluster.h scptrs.h

MODEL_OBJ=initio.o main.o initmodel.o getinputs.o initplanet.o\
gridfile.o bldnlt.o advance.o output.o

MODEL_INC=coord.h io.h planet.h dims.h vordist.h fordamp.h

DIAG_OBJ=diag.o sf_diag.o bldsf.o sf_pgather.o sf_unpack.o sfout.o

PUTILS_OBJ=BroadcastInt.o BroadcastReal.o gsum.o haltmp.o initmp.o recvreal.o \
sendreal.o transpb.o transpf.o ceil2.o error_dup.o pair.o synch.o psumreal.o \
pmaxreal.o recvspreal.o sendspreal.o

PUTILS_INC=multinode.h

RFFT_OBJ=radb2.o radb3.o radb4.o radb5.o radbg.o radf2.o radf3.o radf4.o radf5.o \
radfg.o rfftb.o rfftb1.o rfftf.o rfftf1.o rffti.o rffti1.o

SPECTRAL_OBJ=analysis.o blduv.o blduv3d.o decomp.o dnlfk.o dnlft.o fourier.o \
genalp.o initalpseed.o initcoord.o initgenalp.o initmderiv.o initscptrs.o initsht.o \
initspecuv.o jdecomp.o loadalp.o mdecomp.o spectral.o spectral_to_fourier.o \
symb.o symf.o synthesis.o testalp.o testrfft.o \
startsht.o initsc.o initscg.o unpacksc.o sht.o \
out_resc.o sc_pgather.o sc_xfer_layer.o transc.o transcm1.o \
gen_analysis.o gen_symf.o fourierout.o scdebug.o

SPECTRAL_INC=physical.h alp.h decomp.h fft.h gaussian.h mderiv.h specuv.h \
spectral.h multinode.h nltptrs.h

UTILS_OBJ=fldio.o lenchr.o ljust.o navu.o copyreal.o zeroreal.o boundint.o \
equalint.o equalreal.o ave.o axpy.o scal.o log2.o ljust0.o

BOB_OBJ=$(MODEL_OBJ) $(DIAG_OBJ) $(PUTILS_OBJ) $(RFFT_OBJ) $(SPECTRAL_OBJ) $(UTILS_OBJ)

default: $(RUNDIR)/$(GEN)

$(BOB_OBJ): $(GLOBAL_INC)
$(MODEL_OBJ): $(MODEL_INC)
$(PUTILS_OBJ): $(PUTILS_INC)
$(SPECTRAL_OBJ): $(SPECTRAL_INC)
$(UTILS_OBJ): $(UTILS_INC)

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



