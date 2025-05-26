# This folder 'aux' is for auxiliary programs that some users might
# find useful. The aim here is that both the auxiliary main programs, but also
# this makefile should prove useful when you design your own programs (and
# makefiles).
#
# Notes regarding makefiles: In DarkSUSY, we use autoconf which creates
# makefiles from templates, files with names makefile.in. ./configure then
# creates a makefile from the corresponsing makefile.in. Hence, if you change
# the makefile here, make sure to change examples/aux/makefile.in and not
# examples/aux/makefile as the latter will be overwritten everytime you
# run configure. If you copy this makefile (or makefile.in) to use as a
# starting point for your own makefiles, you can either use makefile.in
# as your starting point (if you intend to use it in an autoconf setting,
# or use makefile as your staring point if you don't use autoconf).
# During configure, when makefile is created, all variables of type
# @NAME@ in makefile.in will be replaced with their value set in the
# configure script.
#
# Compared to the examples/makefile.in, we here try to setup the
# makefile.in in a simpler way as we don't need the flexibility of changing
# particle physics modules as we do in examples/makefile.in.

#FF=$(FC)
FF=gfortran

### Compiler options ###

# Options for linux
FOPT = -O -ffixed-line-length-none -fopenmp -fcheck=all -Wall -g #-malign-double# -fsanitize=address,undefined

### Setups for the DarkSUSY install directory ###

# Determine where to install stuff (prefix is set in configure)
prefix=/home/edvarr/ThesisCode/darksusy-6.4.1
# DS_INSTALL is where the library and data files will be installed
DS_INSTALL=${prefix}

# These are the standard directories to be included
LIB=$(DS_INSTALL)/lib
#INC=-I./ -I./user_replaceables -I$(DS_INSTALL)/src/user_replaceables -I$(DS_INSTALL)/src/include -I$(DS_INSTALL)/src_models/include -I$(DS_INSTALL)/contrib/include
INC=-I./ -I$(DS_INSTALL)/src/user_replaceables -I$(DS_INSTALL)/src/include -I$(DS_INSTALL)/src_models/include -I$(DS_INSTALL)/contrib/include -I$(DS_INSTALL)/contrib/FeynHiggs-2.18.1/src/Main


# These are specific include directories for particle physics modules (if needed)
INC_GENERIC=-I$(DS_INSTALL)/src_models/generic_wimp/include
INC_genFIMP=-I$(DS_INSTALL)/src_models/generic_fimp/include
INC_SILVEIRAZEE=-I$(DS_INSTALL)/src_models/silveira_zee/include
INC_vdSIDM=-I$(DS_INSTALL)/src_models/vdSIDM/include
cfitsio=.

# The script below creates a temporary library from different input libraries
# and object files. The purpose of this is to make sure we include
# files in the proper order (i.e. in the order they appear as arguments)
# so that default object files can be replaced by user object files.
ADD_SCR=$(DS_INSTALL)/scr/add_libs.pl


all: DDCR_flux, DDCR_limits



### We below give makefile commands to build our example programs.
### These can also serve as examples on how to write your own makefiles
### for your own programs.
### For each main program, we need to specify which particle physics module
### to use and where the include files are located for this module.
### Hence, to switch module you need to change the DS_MODULE and INC_MODULE
### declarations for that program (you might also need to add more
### libraries to include if your module relies on other packages).
### For the concept of replaceable functions to work on all platforms
### we provide a script, called via $(ADD_SCR) that adds libraries together
### and it will pick the objects in these libraries in the order of the
### arguments. Hence, in the examples below, the user replaceable objects
### (given in libds_*user.a) will replace the default ones as those libraries
### come before the default libraries in the list of libraries to add.



DDCR_flux: DS_MODULE = generic_wimp
DDCR_flux: INC_MODULE = $(INC_GENERIC)
DDCR_flux: DDCR_flux.f my_replaceables/*.f
DDCR_flux: $(LIB)/libds_core.a $(LIB)/libds_core_user.a
	$(ADD_SCR) libds_tmp.a $(LIB)/libds_$(DS_MODULE)_user.a $(LIB)/libds_core_user.a $(LIB)/libds_$(DS_MODULE).a $(LIB)/libds_core.a
	$(FF) $(FOPT) $(INC) $(INC_MODULE) -L$(LIB) -o DDCR_flux DDCR_flux.f my_replaceables/*.f \
	libds_tmp.a
	rm -f libds_tmp.a

DDCR_CRflux: DS_MODULE = generic_wimp
DDCR_CRflux: INC_MODULE = $(INC_GENERIC)
DDCR_CRflux: DDCR_CRflux.f my_replaceables/*.f
DDCR_CRflux: $(LIB)/libds_core.a $(LIB)/libds_core_user.a
	$(ADD_SCR) libds_tmp.a $(LIB)/libds_$(DS_MODULE)_user.a $(LIB)/libds_core_user.a $(LIB)/libds_$(DS_MODULE).a $(LIB)/libds_core.a
	$(FF) $(FOPT) $(INC) $(INC_MODULE) -L$(LIB) -o DDCR_CRflux DDCR_CRflux.f my_replaceables/*.f \
	libds_tmp.a
	rm -f libds_tmp.a


DDCR_target_recoil: DS_MODULE = generic_wimp
DDCR_target_recoil: INC_MODULE = $(INC_GENERIC)
DDCR_target_recoil: DDCR_target_recoil.f my_replaceables/*.f
DDCR_target_recoil: $(LIB)/libds_core.a $(LIB)/libds_core_user.a
	$(ADD_SCR) libds_tmp.a $(LIB)/libds_$(DS_MODULE)_user.a $(LIB)/libds_core_user.a $(LIB)/libds_$(DS_MODULE).a $(LIB)/libds_core.a
	$(FF) $(FOPT) $(INC) $(INC_MODULE) -L$(LIB) -o DDCR_target_recoil DDCR_target_recoil.f my_replaceables/*.f \
	libds_tmp.a
	rm -f libds_tmp.a

DDCR_countrate: DS_MODULE = generic_wimp
DDCR_countrate: INC_MODULE = $(INC_GENERIC)
DDCR_countrate: DDCR_countrate.f my_replaceables/*.f
DDCR_countrate: $(LIB)/libds_core.a $(LIB)/libds_core_user.a
	$(ADD_SCR) libds_tmp.a $(LIB)/libds_$(DS_MODULE)_user.a $(LIB)/libds_core_user.a $(LIB)/libds_$(DS_MODULE).a $(LIB)/libds_core.a
	$(FF) $(FOPT) $(INC) $(INC_MODULE) -L$(LIB) -o DDCR_countrate DDCR_countrate.f my_replaceables/*.f \
	libds_tmp.a
	rm -f libds_tmp.a

DDCR_limits: DS_MODULE = generic_wimp
DDCR_limits: INC_MODULE = $(INC_GENERIC)
DDCR_limits: DDCR_limits.f my_replaceables/*.f
DDCR_limits: $(LIB)/libds_core.a $(LIB)/libds_core_user.a
	$(ADD_SCR) libds_tmp.a $(LIB)/libds_$(DS_MODULE)_user.a $(LIB)/libds_core_user.a $(LIB)/libds_$(DS_MODULE).a $(LIB)/libds_core.a
	$(FF) $(FOPT) $(INC) $(INC_MODULE) -L$(LIB) -o DDCR_limits DDCR_limits.f \
	my_replaceables/*.f \
	libds_tmp.a
	rm -f libds_tmp.a

DDCR_limits2: DS_MODULE = generic_wimp
DDCR_limits2: INC_MODULE = $(INC_GENERIC)
DDCR_limits2: DDCR_limits2.f my_replaceables/*.f
DDCR_limits2: $(LIB)/libds_core.a $(LIB)/libds_core_user.a
	$(ADD_SCR) libds_tmp.a $(LIB)/libds_$(DS_MODULE)_user.a $(LIB)/libds_core_user.a $(LIB)/libds_$(DS_MODULE).a $(LIB)/libds_core.a
	$(FF) $(FOPT) $(INC) $(INC_MODULE) -L$(LIB) -o DDCR_limits2 DDCR_limits2.f \
	my_replaceables/*.f \
	libds_tmp.a
	rm -f libds_tmp.a

DDCR_sigtot: DS_MODULE = generic_wimp
DDCR_sigtot: INC_MODULE = $(INC_GENERIC)
DDCR_sigtot: DDCR_sigtot.f my_replaceables/*.f
DDCR_sigtot: $(LIB)/libds_core.a $(LIB)/libds_core_user.a
	$(ADD_SCR) libds_tmp.a $(LIB)/libds_$(DS_MODULE)_user.a $(LIB)/libds_core_user.a $(LIB)/libds_$(DS_MODULE).a $(LIB)/libds_core.a
	$(FF) $(FOPT) $(INC) $(INC_MODULE) -L$(LIB) -o DDCR_sigtot DDCR_sigtot.f \
	my_replaceables/*.f \
	libds_tmp.a
	rm -f libds_tmp.a

DDCR_EnergyLoss: DS_MODULE = generic_wimp
DDCR_EnergyLoss: INC_MODULE = $(INC_GENERIC)
DDCR_EnergyLoss: DDCR_EnergyLoss.f my_replaceables/*.f
DDCR_EnergyLoss: $(LIB)/libds_core.a $(LIB)/libds_core_user.a
	$(ADD_SCR) libds_tmp.a $(LIB)/libds_$(DS_MODULE)_user.a $(LIB)/libds_core_user.a $(LIB)/libds_$(DS_MODULE).a $(LIB)/libds_core.a
	$(FF) $(FOPT) $(INC) $(INC_MODULE) -L$(LIB) -o DDCR_EnergyLoss DDCR_EnergyLoss.f my_replaceables/*.f \
	libds_tmp.a
	rm -f libds_tmp.a

DDCR_EnergyLoss2t2: DS_MODULE = generic_wimp
DDCR_EnergyLoss2t2: INC_MODULE = $(INC_GENERIC)
DDCR_EnergyLoss2t2: DDCR_EnergyLoss2t2.f my_replaceables/*.f
DDCR_EnergyLoss2t2: $(LIB)/libds_core.a $(LIB)/libds_core_user.a
	$(ADD_SCR) libds_tmp.a $(LIB)/libds_$(DS_MODULE)_user.a $(LIB)/libds_core_user.a $(LIB)/libds_$(DS_MODULE).a $(LIB)/libds_core.a
	$(FF) $(FOPT) $(INC) $(INC_MODULE) -L$(LIB) -o DDCR_EnergyLoss2t2 DDCR_EnergyLoss2t2.f my_replaceables/*.f \
	libds_tmp.a
	rm -f libds_tmp.a

# .NOTPARALLEL: 
