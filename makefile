##############################################################################
#
#   Makefile for the program PB3D (Peeling Ballooning in 3D)
#   Toon Weyens
#
##############################################################################

##############################################################################
#   Directories
##############################################################################
# BLAS/LAPACK
BLASLAPACK_DIR=''# 1. XPS 9360
#BLASLAPACK_DIR=$(COMPILE_DIR)# 2. ITER
#BLASLAPACK_DIR=''# 3. GEORGE

# LIBSTELL (Note that by default, unlogically, everything is in bin!)
LIBSTELL_DIR=/opt/stellinstall/bin# 1. XPS 9360
#LIBSTELL_DIR=$(COMPILE_DIR)/bin# 2. ITER
#LIBSTELL_DIR=$(HOME)/bin# 3. GEORGE

# HDF5
# (from http://www.hdfgroup.org/ftp/HDF5/examples/howto/makefiles/Makefilef)
HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/openmpi# 1. XPS 9360
#HDF5_DIR=$(COMPILE_DIR)# 2. ITER
#HDF5_DIR_INC=/usr/include/hdf5/openmpi# 3. GEORGE
#HDF5_DIR_LIB=/usr/lib/x86_64-linux-gnu/hdf5/openmpi# 3. GEORGE

# NETCDF
NETCDFF_DIR=/opt/netcdf-fortran-4.4.4/4.4.4#  1. XPS 9360
#NETCDFF_DIR=$(COMPILE_DIR)#  2. ITER
#NETCDFF_DIR=$(HOME)#  3. GEORGE

# PETSC
#PETSC_ARCH = debug-complex
PETSC_ARCH = complex# 1. XPS 93600
#PETSC_ARCH = arch-linux2-c-opt# 2. ITER
PETSC_DIR = /opt/petsc-3.7.6# 1. XPS 9360
#PETSC_DIR=$(COMPILE_DIR)# 2. ITER
#PETSC_DIR = $(HOME)/Programs/petsc-3.6.4# 3. GEORGE

# SLEPC
SLEPC_DIR=/opt/slepc-3.7.4# 1. XPS 9360
#SLEPC_DIR=$(COMPILE_DIR)# 2. ITER
#SLEPC_DIR = $(HOME)/Programs/slepc-3.6.3# 3. GEORGE

# PB3D
PB3D_DIR = /opt/PB3D# 1. XPS 9360
#PB3D_DIR = $(HOME)/Programs_MPICH3.1.3/PB3D# 2. ITER
#PB3D_DIR = $(HOME)/Programs/PB3D# 3. GEORGE


##############################################################################
#   Include
##############################################################################
include  $(PETSC_DIR)/lib/petsc/conf/variables
include  $(SLEPC_DIR)/lib/slepc/conf/slepc_variables

INCLUDE = -I$(LIBSTELL_DIR)/libstell_dir \
  $(PETSC_FC_INCLUDES) \
  $(SLEPC_INCLUDE) \
  -I/usr/include/hdf5/openmpi \
  -I$(PB3D_DIR)/include#1. XPS 9360

#INCLUDE = -I$(LIBSTELL_DIR)/libstell_dir \
  #$(PETSC_FC_INCLUDES) \
  #$(SLEPC_INCLUDE) \
  #-I$(PB3D_DIR)/include#2. ITER

#INCLUDE = -I$(LIBSTELL_DIR)/libstell_dir \
  #-I$(HDF5_DIR_INC) \
  #$(PETSC_FC_INCLUDES) \
  #$(SLEPC_INCLUDE) \
  #-I$(PB3D_DIR)/include#3. GEORGE


##############################################################################
#   Link
##############################################################################
LINK = -llapack -lblas \
  $(LIBSTELL_DIR)/libstell.a \
  -L$(HDF5_DIR) -lhdf5_fortran -lhdf5 \
  -L$(NETCDFF_DIR)/lib -lnetcdff \
  -Wl,-R$(NETCDFF_DIR)/lib \
  $(PETSC_LIB) \
  $(SLEPC_LIB) \
  libdfftpack.a libfoul.a# 1. XPS 9360

#LINK = -L$(BLASLAPACK_DIR)/lib -lblas -llapack \
  #$(LIBSTELL_DIR)/libstell.a \
  #-L$(NETCDFF_DIR)/lib -lnetcdff \
  #-L$(HDF5_DIR)/lib -lhdf5_hl -lhdf5 -lhdf5_fortran -ldl -lm -lz \
  #$(PETSC_LIB) \
  #$(SLEPC_LIB) \
  #libdfftpack.a libfoul.a# 2. ITER

#LINK = -llapack -lblas \
  #$(LIBSTELL_DIR)/libstell.a \
  #-L$(HDF5_DIR_LIB) -lhdf5_fortran -lhdf5 \
  #-L$(NETCDFF_DIR)/lib -lnetcdff \
  #-Wl,-R$(NETCDFF_DIR)/lib \
  #$(PETSC_LIB) \
  #$(SLEPC_LIB) \
  #libdfftpack.a libfoul.a# 3. GEORGE


##############################################################################
#   Compiler
##############################################################################
COMPILER=mpifort
#COMPILER=/opt/scorep-3.0/bin/scorep mpifort


##############################################################################
#   Linker
##############################################################################
LINKER=mpifort
#LINKER=/opt/scorep-3.0/bin/scorep mpifort


##############################################################################
#   Compiler flags
#  	options (used with -D[name]):
# 		ldebug: debug
# 		lIB: infiniband
# 		lwith_gnu: use GNU compiler [default]
# 		lwith_intel: use INTEL compiler, (checked for version 12.0.2)
#   note: INTEL warning 6536 is suppressed, which informs about extra "USE".
#   note: INTEL warning 6843 is suppressed, which informs about empty intent(out) variables
##############################################################################
COMP_FLAGS = -g -Og -Wall -Wextra -pedantic -fimplicit-none -fbacktrace -fno-omit-frame-pointer -fcheck=all -cpp -Dldebug# debug, profiling with gprof2dot, GCC
#COMP_FLAGS = -O3 -fbacktrace -g -fimplicit-none -fno-omit-frame-pointer -cpp# optimized, GCC

#COMP_FLAGS = -O0 -DlIB -Dldebug -g -heap-arrays 100 -recursive -ftrapuv -check bounds -check uninit -traceback -implicitnone -fno-omit-frame-pointer -cpp -Dlwith_intel -diag-disable 6536 -diag-disable 6843# debug, profiling with gprof2dot, INTEL
#COMP_FLAGS = -O3 -DlIB -traceback -g -heap-arrays 100 -recursive -implicitnone -fno-omit-frame-pointer -cpp -Dlwith_intel -diag-disable 6536 -diag-disable 6843# optimized, INTEL

COMP_FLAGS_EX= -O2 -w

COMP_FLAGS_F= -O2 -funroll-loops -fexpensive-optimizations


##############################################################################
#   Link flags
##############################################################################
LINK_FLAGS = -fPIC


##############################################################################
#   Prepare
##############################################################################
# Add "Modules" to the search path for the prerequisites
VPATH = Modules:Libraries

# Contains list of source files (.o) and dependencies
DEPLIST = PB3D.dep
OBJLIST = ObjectList# defines "ObjectFiles"

# Includes source files and dependency list
include $(DEPLIST)# Dependencies of all the objects
include $(OBJLIST)# Names of all the objects


##############################################################################
#   Rules
##############################################################################
all:	PB3D POST

PB3D:	$(ObjectFiles) libdfftpack.a libfoul.a PB3D.o
	$(LINKER) -o $@ $(ObjectFiles) PB3D.o $(LINK) $(LINK_FLAGS)

POST:	$(ObjectFiles) libdfftpack.a libfoul.a POST.o
	$(LINKER) -o $@ $(ObjectFiles) POST.o $(LINK) $(LINK_FLAGS)

libdfftpack.a: 	dfft.o
	ar -rcs libdfftpack.a dfft.o

libfoul.a: 	foul.o
	ar -rcs libfoul.a foul.o

%.o: %.f90
	$(COMPILER) $(INCLUDE) $(COMP_FLAGS) -c $<

%.o: %.f
	$(COMPILER) $(COMP_FLAGS_F) -c $<

dfft.o: dfft.f
	$(COMPILER) $(COMP_FLAGS_EX) -c $<

foul.o: foul.f90
	$(COMPILER) $(COMP_FLAGS_EX) -c $<

clean:
	@rm -f *.o *.a *.mod *~ fort.* 

clean_all:
	@rm -f *.o *.mod *~ fort.* PB3D POST

code_stats:
	cloc .
	#@find . -name '*.f90' | xargs wc -l
