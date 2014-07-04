##############################################################################
#
#   Makefile for the program PB3D (Peeling Ballooning in 3D)
#   Toon Weyens
#
##############################################################################

##############################################################################
#   Paths
##############################################################################
HOME_BIN = /home/toon/bin

# Comiler
COMP_DIR = /usr/bin/gfortran # gfortran
#COMP_DIR = /usr/bin/f95 # f95
#COMP_DIR = /usr/local/solarisstudio12.3/bin/f95 # oracle f95

# Linker
LINK_DIR = /usr/bin/g++ # g++ (needed for C++ preprocessing)

# Add "Modules" to the search path for the prerequisites
VPATH = Modules

# Contains list of source files (.o) and dependencies
DEPLIST = PB3D.dep
OBJLIST = ObjectList

# Includes source files and dependency list
include $(DEPLIST) # Dependencies of all the objects
include $(OBJLIST) # Names of all the objects

##############################################################################
#   Compiler specifications
##############################################################################
# compiler flags
COMP_FLAGS = -g -O0 -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace -pg
#COMP_FLAGS = -O3

# compiler include
COMP_INC = -I/usr/include -I$(HOME_BIN)/libstell_dir

# compiler command
COMPILE = $(COMP_DIR) $(COMP_INC) $(COMP_FLAGS)

##############################################################################
#   Link specifications
##############################################################################
# link flags
LINK_FLAGS = -fPIC

# libraries
LINK_LIB = $(HOME_BIN)/libstell.a -L/usr/lib -lgfortran -lnetcdff -lnetcdf -llapack -lblas

# link command
LINK    = $(LINK_DIR) $(LINK_FLAGS) $(COMP_FLAGS)

##############################################################################
#   Rules
##############################################################################
PB3D:  $(ObjectFiles)
	$(LINK) -o $@ $(ObjectFiles) $(LINK_LIB)
%.o : %.f90
	$(COMPILE) -c $<
clean:
	- rm -f *.o *.mod *~ gmon.out output.xdot fort.* results.txt PB3D_out* tempoutput.dat

output:
	- rm -f *.m *.nc
