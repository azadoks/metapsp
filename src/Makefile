# Makefile for ONCVPSP
#
# Copyright (c) 1989-2024 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
# University
# Copyright (c) 2015-2017 by Alberto Garcia (PSML part)
#

include ../make.inc

OBJS_ALLM =	globalmath_loc.o r2scanmod.o derivs.o \
		exc_libxc_stub.o fcrossover.o fphsft_m.o gnu_script_m.o \
		lschfb_m.o lschfs_m.o lschkb_m.o lschps_m.o \
		lschvkbbe_m.o lschvkbb_m.o lschvkbs_m.o polysmooth.o \
                gauss_smooth.o \
		pspot_m.o rtloc_m.o run_diag_m.o excr2scan01_m.o upfout_m.o \
		run_optimize_m.o run_phsft_m.o run_plot_m.o run_vkb_m.o \
		sratom_m.o vkboutwf_m.o vkbphsft_m.o vout_m.o wellstate_m.o \
		run_ghosts_m.o lschpsbar_m.o

OBJS_ALL =	aeo.o check_data.o const_basis.o dpnint.o \
		eresid.o exchdl.o exc_off_pbe.f excpzca.o excwig.o \
		fpovlp.o linout.o lschfb.o \
		optimize.o \
		qroots.o sbf8.o sbf_basis_con.o \
		sbf_basis.o sbf_rc_der.o sratom.o tfapot.o \
		vout.o vpinteg.o vploc.o wf_rc_der.o

OBJS_M =	oncvpsp_m.o 

OBJS_MX =	oncvpsp_mx.o 

OBJS_PSML = 	m_psmlout_m.o m_getopts.o m_uuid.o


OBJS_OPSPM =	$(OBJS_PSML) $(OBJS_ALL) $(OBJS_ALLM) $(OBJS_M)

OBJS_OPSPMX =	$(OBJS_ALL) $(OBJS_ALLM) $(OBJS_MX)

oncvpspm:	$(OBJS_OPSPM)
		$(FLINKER) -o oncvpspm.x $(OBJS_OPSPM) $(LIBS)

oncvpspmx:	$(OBJS_OPSPMX)
		$(FLINKER) -o oncvpspmx.x $(OBJS_OPSPMX) $(LIBS)


m_psmlout_m.o: m_uuid.o m_getopts.o

#
# This file might need to be compiled without optimization
#
exc_off_pbe.o: FFLAGS=$(FFLAGS_NOCHECKS)

.SUFFIXES: .c .o .f .f90 .F90

.F90.o:
	$(FCCPP) $(FCCPPFLAGS) $*.F90  > $*_cpp.f90
	$(F90) $(FFLAGS) $(INC) -o $*.o -c $*_cpp.f90

.f90.o:
	$(F90) $(FFLAGS) $(INC) -c $*.f90

.f.o:
	$(F90) $(FFLAGS) $(INC) -c $*.f

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $*.c

clean:
	 /bin/rm -f *.o  *.d  *.mod *.log *_cpp.f90

