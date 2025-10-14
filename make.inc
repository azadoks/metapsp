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


##### Edit the following LAPACK and BLAS library paths for your system ####

LIBS = -L/usr/local/lapack/lib -llapack -lrefblas


##### [Optional] Uncomment and edit the following libxmlf90 paths ####
# Instructions for building libxmlf90 are in psml_build.txt

# LIBS += -L/usr/local/libxml/lib -lxmlf90
# FFLAGS += -I/usr/local/libxml/include


##### Edit the following optimization flags for your system ####

FFLAGS     = -O3 -ffast-math -funroll-loops -g -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=snan
CFLAGS     = -O3
