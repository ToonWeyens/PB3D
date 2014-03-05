# Provisional makefile for the program PB3D (Peeling Ballooning in 3D)

SHELL   = /bin/sh
#MYHOME  = $(HOME)
#PRECOMP = /lib/cpp -P -traditional -DLINUX -DNEED_BLAS -DSILO_AVAIL -DNETCDF
HOME_BIN= /home/toon/bin
VISIT_DIR = /usr/local/visit/2.7.1/linux-x86_64/libsim/V2
LIB_LINK= $(HOME_BIN)/libstell.a -L/usr/lib -lgfortran -lnetcdff -lnetcdf -llapack -lblas # $(VISIT_DIR)/lib/libsimV2.a $(VISIT_DIR)/lib/libsimV2f.a #-lgfortranbegin 
LIB     = libstell.a
LIB_DIR = /home/toon/Documents/School/PHD/Stellinstal/LIBSTELL
FLAGS = -g -O0 -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace # for debugging
#FLAGS = -O3 ! for optimization

COMPILE = /usr/bin/gfortran -I/usr/include -I$(HOME_BIN)/libstell_dir # -I$(VISIT_DIR)/include #-I/usr/lib/fortran/modules/plplot -ffixed-form
#COMPILE_FREE = gfortran -I/usr/include -I/usr/lib/fortran/modules/plplot -ffree-form
FFILE   = '$*''.f'
CFILE   = '$*''.c'
F90FILE = '$*''.f90'
#LINK    = g++ -fPIC $(FLAGS) $(SFLAGS) -o
LINK    = g++ -fPIC $(FLAGS) $(SFLAGS) -o
#MOD_PATH= -I
SPATH   = modules

#Contains list of source files (.o) and dependencies
DEPLIST = PB3D.dep
OBJLIST = ObjectList

#Includes source files and dependency list
include $(DEPLIST) # Dependencies of all the objects
include $(OBJLIST) # Names of all the objects
VPATH = $(SPATH)


PB3D:  $(ObjectFiles)
	$(LINK) $@ $(ObjectFiles) $(LIB_LINK)
%.o : %.f90
	$(COMPILE) $(FLAGS) -c $<
$(LIB) :
#	@cd $(LIB_DIR); make release # If you want the library to be updated. 

clean:
	- rm -f *.o *.mod




#.SUFFIXES :
#.SUFFIXES : .f .f90 .o
#PB3D:  $(LIB) $(ObjectFiles)
#	$(LINK) $@ $(ObjectFiles) $(LIB_LINK)
##Compile object files defined in OBJLIST.
#.f.o :
#	@if grep -q '^#if' $<; \
#      then \
#         cp $< $(CFILE); \
#         echo '$(PRECOMP) $<'; $(PRECOMP) $(CFILE) $(FFILE); \
#         rm -f $(CFILE); echo '$(COMPILE) $(FLAGS) $(MOD_PATH).. -c $<'; \
#         $(COMPILE) $(FLAGS) $(MOD_PATH).. -c $(FFILE); \
#      else \
#         echo '$(COMPILE) $(FLAGS) $(MOD_PATH). -c $<'; \
#         $(COMPILE) $(FLAGS) $(MOD_PATH). -c $<; \
#      fi
#
#.f90.o :
#	@if grep -q '^#if' $<; \
#      then \
#         cp $< $(CFILE); \
#         echo '$(PRECOMP) $<'; $(PRECOMP) $(CFILE) $(F90FILE); \
#         rm -f $(CFILE); echo '$(COMPILE_FREE) $(FLAGS) $(MOD_PATH).. -c $<'; \
#        $(COMPILE_FREE) $(FLAGS) $(MOD_PATH).. -c $(F90FILE); \
#      else \
#         echo '$(COMPILE_FREE) $(FLAGS) $(MOD_PATH). -c $<'; \
#         $(COMPILE_FREE) $(FLAGS) $(MOD_PATH).. -c $<; \
#      fi
