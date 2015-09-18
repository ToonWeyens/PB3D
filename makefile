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
#PB3D_DIR = $(HOME)/Programs/PB3D# quadrivium

# Comiler
COMP_DIR = /usr/bin/mpif90# laptop
#COMP_DIR = /share/apps/openmpi/1.6.5/gcc-4.6.4/bin/mpif90# quadrivium

# Linker (needed for C++ preprocessing)
LINK_DIR = /usr/bin/g++# laptop
#LINK_DIR = /share/apps/gcc/4.6.4/bin/g++# quadrivium

# PETSC and SLEPC directories
# version 3.5.3
#PETSC_DIR = /opt/petsc-3.5.3
#PETSC_ARCH = debug-complex-seq
#SLEPC_DIR = /opt/slepc-3.5.3
#include  $(PETSC_DIR)/conf/variables
#include  $(SLEPC_DIR)/conf/slepc_variables
# version 3.6.1 on laptop
PETSC_DIR = /opt/petsc-3.6.1
SLEPC_DIR = /opt/slepc-3.6.0
# version 3.6.1 on quadrivium
#PETSC_DIR = $(HOME)/Programs/petsc-3.6.1
#SLEPC_DIR = $(HOME)/Programs/slepc-3.6.1
# version 3.6.1 for both
PETSC_ARCH = debug-complex
include  $(PETSC_DIR)/lib/petsc/conf/variables
include  $(SLEPC_DIR)/lib/slepc/conf/slepc_variables

# HDF5
# (from http://www.hdfgroup.org/ftp/HDF5/examples/howto/makefiles/Makefilef)
HDF5_lib = /usr/lib/x86_64-linux-gnu# laptop repository
HDF5_inc = /usr/include/#laptop repository
#HDF5_lib = /opt/HDF5_1.8.15-patch1/lib# laptop 1.8.15 patched
#HDF5_inc = /opt/HDF5_1.8.15-patch1/include# laptop 1.8.15 patched
#HDF5_lib = $(HOME)/lib# quadrivium 1.8.15 patched
#HDF5_inc = $(HOME)/include# quadrivium 1.8.15 patched

# Add "Modules" to the search path for the prerequisites
VPATH = Modules

# Contains list of source files (.o) and dependencies
DEPLIST = PB3D.dep
OBJLIST = ObjectList# defines "ObjectFiles"

# Includes source files and dependency list
include $(DEPLIST)# Dependencies of all the objects
include $(OBJLIST)# Names of all the objects

##############################################################################
#   Compiler specifications
##############################################################################
# compiler flags
COMP_FLAGS = -g -O0 -Wall -Wextra -pedantic -fimplicit-none -fbacktrace -pg -fno-omit-frame-pointer -cpp -Dldebug# profiling
#COMP_FLAGS = -g -O0 -Wall -Wextra -pedantic -fimplicit-none -fbacktrace -cpp -Dldebug# no profiling

# compiler include
COMP_INC = -I$(HDF5_inc) -I$(HOME_BIN)/libstell_dir -I$(PB3D_DIR)/include

# compiler command
COMPILE = $(COMP_DIR) $(COMP_INC) $(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE) $(COMP_FLAGS)

##############################################################################
#   Link specifications
##############################################################################
# link flags
LINK_FLAGS = -fPIC -pg

# libraries
LINK_LIB = $(HOME_BIN)/libstell.a -lgfortran -llapack -lblas \
	$(HDF5_lib)/libhdf5_fortran.a $(HDF5_lib)/libhdf5.a \
	-lz -lpthread -ldl -lm


# link command
LINK    = $(LINK_DIR) $(LINK_FLAGS)

##############################################################################
#   Rules
##############################################################################
all:	PB3D PB3D_POST

PB3D:	$(ObjectFiles) PB3D.o
	$(LINK) -o $@ $(ObjectFiles) PB3D.o $(LINK_LIB) $(PETSC_LIB) $(SLEPC_LIB)

PB3D_POST:	$(ObjectFiles) PB3D_POST.o
	$(LINK) -o $@ $(ObjectFiles) PB3D_POST.o $(LINK_LIB) $(PETSC_LIB) $(SLEPC_LIB)

%.o : %.f90
	$(COMPILE) -c $<
clean:
	@rm -f *.o *.mod *~ fort.* 

code_stats:
	cloc .
	#@find . -name '*.f90' | xargs wc -l
