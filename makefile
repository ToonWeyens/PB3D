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
HOME_BIN = /home/toon/bin

# PB3D directory
PB3D_DIR = /home/toon/Documents/PHD/PB3D

# Comiler
#COMP_DIR = /usr/bin/gfortran # gfortran
COMP_DIR = /usr/bin/mpif90 # mpi fortran

# Linker
LINK_DIR = /usr/bin/g++ # g++ (needed for C++ preprocessing)

# PETSC directory
PETSC_DIR = /opt/petsc/petsc-3.4.5
PETSC_ARCH = debug-complex
#PETSC_ARCH = debug-complex-seq
SLEPC_DIR = /opt/slepc/slepc-3.4.4
include  $(SLEPC_DIR)/conf/slepc_variables
include  $(PETSC_DIR)/conf/variables

# Add "Modules" to the search path for the prerequisites
VPATH = Modules

# Contains list of source files (.o) and dependencies
DEPLIST = PB3D.dep
OBJLIST = ObjectList # defines "ObjectFiles"

# Includes source files and dependency list
include $(DEPLIST) # Dependencies of all the objects
include $(OBJLIST) # Names of all the objects

##############################################################################
#   Compiler specifications
##############################################################################
# compiler flags
COMP_FLAGS = -g -O0 -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace -pg -cpp -Dldebug
#COMP_FLAGS = -g -O0 -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace
#COMP_FLAGS = -O3

# compiler include
COMP_INC = -I/usr/include -I$(HOME_BIN)/libstell_dir -I$(PB3D_DIR)/include

# compiler command
COMPILE = $(COMP_DIR) $(COMP_INC) $(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE) $(COMP_FLAGS)

##############################################################################
#   Link specifications
##############################################################################
# link flags
LINK_FLAGS = -fPIC

# libraries
LINK_LIB = $(HOME_BIN)/libstell.a -lgfortran -lnetcdff -lnetcdf -llapack -lblas
#LINK_LIB = $(HOME_BIN)/libstell.a -L/usr/lib -lgfortran -lnetcdff -lnetcdf -llapack -lblas

# link command
LINK    = $(LINK_DIR) $(LINK_FLAGS)

##############################################################################
#   Rules
##############################################################################
PB3D:	$(ObjectFiles)
	$(LINK) -o $@ $(ObjectFiles) $(LINK_LIB) $(PETSC_LIB) $(SLEPC_LIB)

%.o : %.f90
	$(COMPILE) -c $<
clean:
	@rm -f *.o *.mod *~ gmon.out output.xdot fort.* results.txt PB3D_out* tempoutput.dat

output:
	@rm -f *.m *.nc Plots/*
