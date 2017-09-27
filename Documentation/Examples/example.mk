##############################################################################
#
#   Example makefile for the program PB3D (Peeling Ballooning in 3D)
#   \author Author: Toon Weyens
#
#   Don't forget to set the directories:
# 		- LIBSTELL_DIR
# 		- HDF5_DIR
# 		- NETCDFF_DIR (note: Fortran library)
# 		- PETSC_DIR
# 		- SLEPC_DIR
# 		- STRUMPACK_DIR
##############################################################################

##############################################################################
#   Include
##############################################################################
## [PETSc and SLEPc trick]
include  $(PETSC_DIR)/lib/petsc/conf/variables
include  $(SLEPC_DIR)/lib/slepc/conf/slepc_variables
## [PETSc and SLEPc trick]

## [PETSc and SLEPc trick inc]
INCLUDE = $(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
## [PETSc and SLEPc trick inc]
## [Libstell special]
INCLUDE += -I$(LIBSTELL_DIR)/libstell_dir
## [Libstell special]
INCLUDE += -I$(STRUMPACK_DIR)/include
## [PB3D include]
INCLUDE += -I$(PB3D_DIR)/include
## [PB3D include]
INCLUDE += -I/usr/include/hdf5/openmpi

##############################################################################
#   Link
##############################################################################
## [PB3D libraries]
LIB_INTERNAL = libdfftpack.a libfoul.a libbspline.a
## [PB3D libraries]

LINK := $(LIB_INTERNAL)

## [PETSc and SLEPc trick lib]
LINK += $(PETSC_LIB)
LINK += $(SLEPC_LIB)
## [PETSc and SLEPc trick lib]
LINK += $(LIBSTELL_DIR)/libstell.a
LINK += -L$(STRUMPACK_DIR)/lib -lstrumpack
LINK += -L$(HDF5_DIR) -lhdf5_fortran -lhdf5
LINK += -L$(NETCDFF_DIR)/lib -lnetcdff
LINK += -Wl,-R$(NETCDFF_DIR)/lib
LINK += -lscalapack -lblacs -lblas -lm
LINK += -lstdc++ -lmpi_cxx


##############################################################################
#   Compiler
##############################################################################
COMPILER=mpifort


##############################################################################
#   Linker
##############################################################################
LINKER=mpifort


##############################################################################
#   Compiler flags
#  	options (used with -D[name]):
# 		ldebug: debug
# 		lIB: infiniband
# 		lwith_gnu: use GNU compiler [default]
# 		lwith_intel: use INTEL compiler, (checked for version 12.0.2)
#   note: INTEL warning 6536 is suppressed, which informs about extra "USE".
#   note: INTEL warning 6843 is suppressed, which informs about empty
#    	intent(out) variables
##############################################################################
COMP_FLAGS = -finit-real=snan -g -Og -Wall -Wextra -pedantic \
	-fimplicit-none -fbacktrace -fno-omit-frame-pointer \
	-fcheck=all -cpp -Dldebug# debug, profiling with gprof2dot, GCC
#COMP_FLAGS = -O3 -fbacktrace -g -fimplicit-none -fno-omit-frame-pointer \
	#-cpp# optimized, GCC

#COMP_FLAGS = -O0 -DlIB -Dldebug -g -heap-arrays 100 -recursive \
	#-ftrapuv -check bounds -check uninit -traceback -implicitnone \
	#-fno-omit-frame-pointer -cpp -Dlwith_intel -diag-disable 6536 \
	#-diag-disable 6843# debug, profiling with gprof2dot, INTEL
#COMP_FLAGS = -O3 -DlIB -traceback -g -heap-arrays 100 -recursive \
	#-implicitnone -fno-omit-frame-pointer -cpp -Dlwith_intel \
	#-diag-disable 6536 -diag-disable 6843# optimized, INTEL

COMP_FLAGS_EX= -O2 -w

COMP_FLAGS_F= -O2 -funroll-loops -fexpensive-optimizations


##############################################################################
#   Link flags
##############################################################################
LINK_FLAGS = -fPIC -finit-real=snan# debug
#LINK_FLAGS = -fPIC# optimized


##############################################################################
#   Prepare
##############################################################################
# Add "Modules" and "Libraries" to the search path for the prerequisites
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

PB3D:	$(ObjectFiles) $(LIB_INTERNAL) PB3D.o
	$(LINKER) -o $@ $(ObjectFiles) PB3D.o $(LINK) $(LINK_FLAGS)

POST:	$(ObjectFiles) $(LIB_INTERNAL) POST.o
	$(LINKER) -o $@ $(ObjectFiles) POST.o $(LINK) $(LINK_FLAGS)

libdfftpack.a: 	dfft.o
	ar -rcs libdfftpack.a dfft.o

libfoul.a: 	foul.o
	ar -rcs libfoul.a foul.o

libbspline.a: 	bspline_sub_module.o
	ar -rcs libbspline.a bspline_sub_module.o

%.o: %.f90
	$(COMPILER) $(INCLUDE) $(COMP_FLAGS) -c $<

%.o: %.f
	$(COMPILER) $(COMP_FLAGS_F) -c $<

dfft.o: dfft.f
	$(COMPILER) $(COMP_FLAGS_EX) -c $<

foul.o: foul.f90
	$(COMPILER) $(COMP_FLAGS_EX) -c $<

bspline_sub_module.o: bspline_sub_module.f90
	$(COMPILER) $(COMP_FLAGS_EX) -c $<

clean:
	@rm -f *.o *.a *.mod *~ fort.* 

clean_all:
	@rm -f *.o *.mod *~ fort.* PB3D POST
