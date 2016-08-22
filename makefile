##############################################################################
#
#   Makefile for the program PB3D (Peeling Ballooning in 3D)
#   Toon Weyens
#
##############################################################################

##############################################################################
#   Paths
##############################################################################
# Home bin directory
HOME_BIN = $(HOME)/bin

# PB3D directory
PB3D_DIR = $(HOME)/Documents/PHD/PB3D# laptop
#PB3D_DIR = $(HOME)/Programs/PB3D# uranus

# Comiler
COMP_DIR = /usr/bin/mpifort# laptop openmpi 1.6.3
#COMP_DIR = /opt/openmpi/1.10.0/bin/mpifort# laptop custom openmpi 1.10.0
#COMP_DIR = /share/apps/openmpi-1.10.1/bin/mpifort# uranus

# Linker (needed for C++ preprocessing)
LINK_DIR = $(COMP_DIR)

# PETSC and SLEPC directories
# version 3.5.3
#PETSC_DIR = /opt/petsc-3.5.3
#PETSC_ARCH = debug-complex-seq
#SLEPC_DIR = /opt/slepc-3.5.3
#include  $(PETSC_DIR)/conf/variables
#include  $(SLEPC_DIR)/conf/slepc_variables
# version 3.6 on laptop
PETSC_DIR = /opt/petsc-3.6.4
SLEPC_DIR = /opt/slepc-3.6.3
# version 3.6 on uranus
#PETSC_DIR = $(HOME)/Programs/petsc-3.6.1
#SLEPC_DIR = $(HOME)/Programs/slepc-3.6.1
# version 3.6 for both
PETSC_ARCH = debug-complex
include  $(PETSC_DIR)/lib/petsc/conf/variables
include  $(SLEPC_DIR)/lib/slepc/conf/slepc_variables

# HDF5
# (from http://www.hdfgroup.org/ftp/HDF5/examples/howto/makefiles/Makefilef)
#HDF5_lib = /usr/lib/x86_64-linux-gnu/hdf5/openmpi# laptop repository
#HDF5_inc = /usr/include/hdf5/openmpi#laptop repository
HDF5_lib = /opt/hdf5-1.8.16/hdf5/lib# laptop 1.8.16
HDF5_inc = /opt/hdf5-1.8.16/hdf5/include/#laptop 1.8.16
#HDF5_lib = /opt/HDF5_1.8.15-patch1/lib# laptop 1.8.15 patched
#HDF5_inc = /opt/HDF5_1.8.15-patch1/include# laptop 1.8.15 patched
#HDF5_lib = $(HOME)/lib# uranus 1.8.15 patched
#HDF5_inc = $(HOME)/include# uranus 1.8.15 patched

# NETCDF
#NETCDF_lib = /usr/lib/# laptop repository
#NETCDF_inc = /usr/include/#laptop repository
NETCDF_lib = /opt/NetCDF-4.4.0/NetCDF/lib# laptop 4.4.0
NETCDF_inc = /opt/NetCDF-4.4.0/NetCDF/include/#laptop 4.4.0
#NETCDF_lib = $(HOME)/lib# uranus 4.4.4
#NETCDF_inc = $(HOME)/include# uranus 4.4.4

# Add "Modules" to the search path for the prerequisites
VPATH = Modules:Libraries

# Contains list of source files (.o) and dependencies
DEPLIST = PB3D.dep
OBJLIST = ObjectList# defines "ObjectFiles"

# Includes source files and dependency list
include $(DEPLIST)# Dependencies of all the objects
include $(OBJLIST)# Names of all the objects

##############################################################################
#   Compiler specifications
#  	options (used with -D[name]):
# 		ldebug: debug
##############################################################################
# compiler flags
#COMP_FLAGS = -g -O0 -Wall -Wextra -pedantic -fimplicit-none -fbacktrace -pg -fno-omit-frame-pointer -fcheck=bounds,array-temps,do,pointer,recursion -cpp -Dldebug# profiling with gprof2dot
COMP_FLAGS = -O3 -fimplicit-none -fno-omit-frame-pointer -cpp# optimized

# compiler include
COMP_INC = -I$(HDF5_inc) -I$(NETCDF_inc) -I$(HOME_BIN)/libstell_dir -I$(PB3D_DIR)/include #-I/opt/openmpi/1.10.0/include

# compiler command
COMPILE = $(COMP_DIR) $(COMP_INC) $(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE) $(COMP_FLAGS)

##############################################################################
#   Link specifications
##############################################################################
# link flags
LINK_FLAGS = -fPIC -pg

# libraries
LINK_LIB = $(HOME_BIN)/libstell.a libdfftpack.a -lgfortran -llapack -lblas \
	-L$(HDF5_lib) -lhdf5_fortran -lhdf5 -L$(NETCDF_lib) -lnetcdf -lnetcdff  \
	-Wl,-R$(NETCDF_lib) -lz -lpthread -ldl -lm# -Wl,-R[PATH] to set to default search path http://superuser.com/questions/192573/how-do-you-specify-the-location-of-libraries-to-a-binary-linux)

# link command
LINK    = $(LINK_DIR) $(LINK_FLAGS)

##############################################################################
#   Rules
##############################################################################
all:	PB3D POST

PB3D:	$(ObjectFiles) libdfftpack.a PB3D.o
	$(LINK) -o $@ $(ObjectFiles) PB3D.o $(LINK_LIB) $(PETSC_LIB) $(SLEPC_LIB)

POST:	$(ObjectFiles) libdfftpack.a POST.o
	$(LINK) -o $@ $(ObjectFiles) POST.o $(LINK_LIB) $(PETSC_LIB) $(SLEPC_LIB)

libdfftpack.a: 	$(ObjectFiles_dfftpack)
	ar -rcs libdfftpack.a $(ObjectFiles_dfftpack)

%.o : %.f90
	$(COMPILE) -c $<

%.o : %.f
	gfortran -O2 -funroll-loops -fexpensive-optimizations -c $<

clean:
	@rm -f *.o *.a *.mod *~ fort.* 

clean_all:
	@rm -f *.o *.mod *~ fort.* PB3D POST

code_stats:
	cloc .
	#@find . -name '*.f90' | xargs wc -l
