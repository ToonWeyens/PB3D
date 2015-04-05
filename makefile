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
PETSC_DIR = /opt/petsc-3.5.3
PETSC_ARCH = debug-complex
#PETSC_ARCH = debug-complex-seq
SLEPC_DIR = /opt/slepc-3.5.3
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
LINK_FLAGS = -fPIC -pg

# libraries (for HDF5, see http://hpc.ucla.edu/hoffman2/software/hdf.php#hdf5f90)
LINK_LIB = $(HOME_BIN)/libstell.a -lgfortran -llapack -lblas \
	   -L/usr/lib/x86_64-linux-gnu/libhdf5hl_fortran.a /usr/lib/x86_64-linux-gnu/libhdf5_hl.a \
           /usr/lib/x86_64-linux-gnu/libhdf5_fortran.a /usr/lib/x86_64-linux-gnu/libhdf5.a -Wl,-z,relro \
           -lpthread -lz -ldl -lm # -Wl,-Bsymbolic-functions -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu
# (To switch between own install of hdf5 and ubuntu packages, probably replace "/usr/lib/x86_64-linux-gnu" by "/opt/HDF5/lib/")
# (ALTERNATIVE from http://www.hdfgroup.org/ftp/HDF5/examples/howto/makefiles/Makefilef:
	   #-I/opt/HDF5/include /opt/HDF5/lib/libhdf5_fortran.a /opt/HDF5/lib/libhdf5.a)


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
	@rm -f *.o *.mod *~ gmon.out output.xdot fort.* results.txt PB3D_out* tempoutput.dat Run/*~

rm_output:
	@rm -f Run/PB3D_out* Run/Plots/* Run/Data/* Run/Scripts/*
code_stats:
	cloc .
	#@find . -name '*.f90' | xargs wc -l
