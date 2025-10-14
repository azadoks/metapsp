# System-dependent makefile options for ONCVPSP / METAPSP
# This must be carefully edited before executing "make" in src
#
# Copyright (c) 1989-2024 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
# University
 
##### Edit the following lines to correspond to your compilers ####

F77        = gfortran
F90        = gfortran
CC         = gcc
FCCPP      = cpp

FLINKER     = $(F90)

#FCCPPFLAGS = -ansi -DLIBXC_VERSION=301  #Use this for versions <301
#FCCPPFLAGS = -ansi -DLIBXC_VERSION=400  #Use this for versions >400

##### Edit the following optimization flags for your system ####

#FFLAGS     = -O3 -ffast-math -funroll-loops

FFLAGS     = -O3 -ffast-math -funroll-loops  -g -Wall -Wextra -Wno-unused -Wsurprising -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=snan

#FFLAGS  = -ffpe-trap=invalid,zero,overflow -g -static

CFLAGS     = -O3		

##### Edit the following LAPACK and BLAS library paths for your system ####

LIBS = -L/usr/local/lapack/lib -llapack -lrefblas

##### Edit the following for to use libxc if available #####
##### Perhaps for a future release of metapsp #####

#OBJS_LIBXC =	exc_libxc_stub.o

# oncvpsp is compatible with libxc
# To build oncvpsp with libxc, uncomment the e following lines and edit
# the paths to point to your libxc library and include directories
# make clean in src before rebuilding after changing this

#LIBS += -L/usr/local/libxc/lib -lxcf90 -lxc
#FFLAGS += -I/usr/local/libxc/include

#LIBS += -L/home/drh/abinit/fallbacks/exports/lib -lxcf90 -lxc
#FFLAGS += -I/home/drh/abinit/fallbacks/exports/include

#OBJS_LIBXC =	functionals.o exc_libxc.o


LIBS += -L/usr/local/libxml/lib -lxmlf90
FFLAGS += -I/usr/local/libxml/include


