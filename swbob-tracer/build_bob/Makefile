###
### GNU makefile for BOB
###

ARCH=$(shell uname -s)
HOST=$(shell hostname -s)
LOCALPATH = $(HOME)/lib/$(ARCH)
TOP=../src
RUNDIR=.
MODEL=bobx

# If using a Mac, set ARCH to Linux - wjms
ifeq ($(ARCH),Darwin)
	ARCH=Linux
endif

VPATH=$(TOP)/include:$(TOP)/model:$(TOP)/putils:$(TOP)/rfft:$(TOP)/spectral:$(TOP)/utils:$(TOP)/diagnostics
BOBPATH=-I$(TOP)/include -I$(TOP)/model -I$(TOP)/putils -I$(TOP)/rfft -I$(TOP)/spectral -I$(TOP)/diagnostics

#########################################
# Tested/supported architectures
# (22 February 2003, Rich Loft)
#
# 1) Linux/Portland Group compilers + MPI (optional, call make _MPI=1)
# 2) IBM AIX/xlf90 compilers        + MPI (optional, call make _MPI=1)
# 3) OSX (Darwin)/Absoft compilers
#########################################

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
   ifeq ($(HOST),laptop)     # and older storms
      COMP      = GNU
      FC        = gfortran
      CC        = cc
      FBASE     = -ffixed-line-length-none -fno-second-underscore
      MPILIB    =
   endif
   CPPFLAGS  = $(BOBPATH)
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
   ifeq ($(RR),1)
      FFLAGS += -DRTERM
   endif
   ifeq ($(FD),1)
      FFLAGS += -DFDRAG
   endif
   ifeq ($(FR),1)
      FFLAGS += -DFRLAX
   endif
   ifeq ($(FRZ),1)
      FFLAGS += -DFRLAXZ
   endif
   ifeq ($(FT),1)
      FFLAGS += -DFTOPO
   endif
   ifeq ($(FH),1)
      FFLAGS += -DFHEAT
   endif
   ifeq ($(FV),1)
      FFLAGS += -DFVORT
   endif
   ifeq ($(FM),1)
      FFLAGS += -DFMECH
   endif
   ifeq ($(FS),1)
      FFLAGS += -DFSPEC
   endif
   ifeq ($(FU),1)
      FFLAGS += -DFURLX
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
   ifeq ($(M2),1)
      FFLAGS += -DM2
   endif
endif

ifeq ($(ARCH),AIX)
   FC        = mpxlf90_r
   CC        = mpcc_r
   CPPFLAGS  = $(BOBPATH)
   FBASE    = -qfixed=132 -WF,-D$(DELORDER) -WF,-D$(NU)
   FOPT      = -O3 -qarch=pwr4 -qtune=pwr4 -qcache=auto
   FFLAGS    = $(FBASE) $(FOPT)
   CBASE     =
   COPT      = -O4 -qarch=pwr4 -qtune=pwr4 -qcache=auto
   CFLAGS    = $(CBASE) $(COPT)
   LDFLAGS   = $(FOPT)
   MPILIB    =
   LOCALLIB  = -L$(LOCALPATH) -ltimer
   LIBS      = $(LOCALLIB) $(MPILIB)

   ifeq ($(_MPI),1)
      FFLAGS += -WF,-DMPI
   endif
   ifeq ($(RR),1)
      FFLAGS += -WF -DRTERM
   endif
   ifeq ($(FR),1)
      FFLAGS += -WF,-DFRLAX
   endif
   ifeq ($(FD),1)
      FFLAGS += -WF,-DFDRAG
   endif
   ifeq ($(FT),1)
      FFLAGS += -WF,-DFTOPO
   endif
   ifeq ($(FH),1)
      FFLAGS += -WF,-DFHEAT
   endif
   ifeq ($(FV),1)
      FFLAGS += -WF,-DFVORT
   endif
   ifeq ($(FM),1)
      FFLAGS += -WF,-DFMECH
   endif
   ifeq ($(FS),1)
      FFLAGS += -WF,-DFSPEC
   endif
   ifeq ($(AD),1)
      FFLAGS += -WF,-DADAPT
   endif
   ifeq ($(HD),1)
      FFLAGS += -WF,-DHDISS
   endif
   ifeq ($(GS),1)
      FFLAGS += -WF,-DGAUSS
   endif

endif

ifeq ($(ARCH),Darwin)
   FC        = f90
   CC        = gcc
   CPPFLAGS  = $(BOBPATH)
   FBASE     = -D$(DELORDER) -D$(NU) -D$(FT) -D$(FR) -f fixed -W 132
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
   CPPFLAGS  = $(BOBPATH)
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
   CPPFLAGS  = $(BOBPATH) -I/opt/mpi/include
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
  CPPFLAGS  = $(BOBPATH) -I/opt/SUNWhpc/HPC3.1/include
  FBASE     = -e -g
  FPREC     =
  FOPT      = -fast
  FFLAGS    = $(FBASE) $(FOPT)
  CBASE     = -I/opt/SUNWhpc/HPC3.1/include
  COPT      = -xO2
  CFLAGS    = $(CBASE) $(COPT)
  LDFLAGS   = $(FOPT)
  LOCALLIB  = -L$(LOCALPATH) -ltimer
  MPILIB    = -lmpi
  LIBS      = $(MPILIB) $(LOCALLIB)
endif

GLOBAL_INC=dims.h nltptrs.h type.h constants.h fields.h\
params.h stptrs.h mpcluster.h scptrs.h timer.h

MODEL_OBJ=main.o initmodel.o getinputs.o\
initmodelctl.o initio.o gridfile.o initfile.o checkheader.o \
bldnlt.o advance.o timeadvance.o output.o timeinit.o \
perfstats.o initperf.o hyper.o \
inittheforcing.o inittopoforcing.o initthpforcing.o initthrforcing.o \
initmomforcing.o initkvcoef.o initktcoef.o inittheforcingic.o initurlxforcing.o

MODEL_INC=coord.h io.h modelctl.h modeltime.h planet.h timeint.h\
timestep.h diag.h hyper.h fordamp.h fordamp0.h

PUTILS_OBJ=BroadcastInt.o BroadcastReal.o gsum.o haltmp.o initmp.o recvreal.o \
sendreal.o transpb.o transpf.o ceil2.o error_dup.o pair.o synch.o psumreal.o \
pmaxreal.o recvspreal.o sendspreal.o

PUTILS_INC=multinode.h

RFFT_OBJ=radb2.o radb3.o radb4.o radb5.o radbg.o radf2.o radf3.o radf4.o radf5.o \
radfg.o rfftb.o rfftb1.o rfftf.o rfftf1.o rffti.o rffti1.o

SPECTRAL_OBJ=analysis.o blduv.o blduv3d.o decomp.o dnlfk.o dnlft.o fourier.o \
genalp.o initalpseed.o initcoord.o initgenalp.o initmderiv.o initscptrs.o initsht.o \
initspecuv.o jdecomp.o loadalp.o mdecomp.o spectral.o spectral_to_fourier.o symb.o \
symf.o synthesis.o testalp.o testrfft.o initsc.o unpacksc.o unpacksc1.o startsht.o sht.o \
sc_pgather.o sc_xfer_layer.o \
initforsc.o forsc.o physical_to_spectral.o  gen_analysis.o gen_symf.o \
out_resc.o transc.o transc1.o spectral_to_physical.o synthesis1.o

SPECTRAL_INC=physical.h alp.h decomp.h fft.h gaussian.h mderiv.h specuv.h \
spectral.h fordamp.h fordamp0.h uspec.h modeltime.h

UTILS_OBJ=fldio.o lenchr.o ljust.o ljust0.o navu.o copyreal.o zeroreal.o boundint.o \
equalint.o equalreal.o ave.o axpy.o scal.o log2.o myrand.o

DIAG_OBJ=diag.o \
fld.o fld_bld.o fld_out.o \
spc.o spc_bld.o \
znl.o znl_pgather.o znl_unpack.o znl_out.o znl_bld.o znl_bldvq.o \
zacc.o ell.o ell_bld.o \
congen.o zerovor.o \
equilib.o equibld.o equiout.o \
sf_pgather.o sf_unpack.o

DIAG_INC=diag.h diagptrs.h congen.h planet.h io.h uspec.h modeltime.h timestep.h fordamp0.h equilib.h

BOB_OBJ=$(MODEL_OBJ) $(PUTILS_OBJ) $(RFFT_OBJ) $(SPECTRAL_OBJ) $(UTILS_OBJ) $(DIAG_OBJ)

default: $(RUNDIR)/$(MODEL)

$(BOB_OBJ): $(GLOBAL_INC)
$(MODEL_OBJ): $(MODEL_INC)
$(PUTILS_OBJ): $(PUTILS_INC)
$(SPECTRAL_OBJ): $(SPECTRAL_INC)
$(UTILS_OBJ): $(UTILS_INC)
$(DIAG_OBJ): $(DIAG_INC)

.SUFFIXES:
.SUFFIXES: .F .o .c

$(RUNDIR)/$(MODEL): $(BOB_OBJ)
	$(FC) -o $@ $(BOB_OBJ) $(LDFLAGS) $(LIBS)

.F.o:
	$(FC) -c $(CPPFLAGS) $(FFLAGS) $<

.c.o:
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $<

clean:
	rm -f *.o $(RUNDIR)/$(MODEL)
