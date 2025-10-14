!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This module contains the following active subroutines:
!      Init_GlobalConstants, sphbes, bessandneumfunc, shift4,
!          linsol, minverse, filter, conthomas, thomas, jbessel, solvbes,
!          shapebes, CALERF, SolveAXeqB, SolveAXeqBM, SolveAXeqB_eig,
!          SVDINVERSE, dirachwfn
!  This module contains the following active functions:
!      ddlog, ddexp, ranx, factorial, hwfn, intjl, kummer,
!      DERF, DERFC, DERFCX, ASSOCIATEDLAGUERRE, THRJ2, GAMMAFUNC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

MODULE GlobalMath_loc

!DRH_edit
!  USE io_tools

  IMPLICIT NONE

  REAL(8), parameter :: inverse_fine_structure=137.035999139d0
  REAL(8), parameter :: ifsalpha2=inverse_fine_structure**2
  REAL(8), parameter :: fsalpha2=1.d0/inverse_fine_structure**2
  REAL(8) :: pi , machine_precision , machine_zero , machine_infinity

  REAL(8), PRIVATE :: minlog,maxlog,minexp,maxexp
  REAL(8), PRIVATE :: minlogarg,maxlogarg,minexparg,maxexparg

CONTAINS

  !****************************************************
  SUBROUTINE Init_GlobalConstants()

    INTEGER :: i
    REAL(8)    :: tmp , a1,a2,a3

    ! Calculate machine accuracy
    machine_precision = 0
    a1 = 4.d0/3.d0
    DO WHILE (machine_precision == 0.d0)
       a2 = a1 - 1.d0
       a3 = a2 + a2 + a2
       machine_precision = ABS(a3 - 1.d0)
    ENDDO

    !machine_zero= machine_precision**4
    machine_zero= machine_precision**5    ! suggested by Marc and Francois
    machine_infinity = 1.d0/machine_zero

    pi = ACOS(-1.d0)

    minlogarg=machine_precision; minlog=LOG(minlogarg)
    maxlogarg=1.d0/machine_precision; maxlog=LOG(maxlogarg)
    minexparg=LOG(machine_precision);  minexp=0.d0
    maxexparg=-LOG(machine_precision);  maxexp=EXP(maxexparg)

!DRH_edits 3 lines
!   write(std_out,*) 'Init_GlobalConstants: machine_precision ', machine_precision
!   write(std_out,*) 'Init_GlobalConstants: machine_zero ', machine_zero
!   write(std_out,*) 'Init_GlobalConstants: machine_infinity ', machine_infinity

    RETURN
  END SUBROUTINE Init_GlobalConstants

  !****************************************************
  FUNCTION ddlog(arg)
    REAL(8) :: arg,ddlog

    IF (arg>maxlogarg) THEN
       ddlog=maxlog
    ELSE IF (arg<minlogarg) THEN
       ddlog=minlog
    ELSE
       ddlog=LOG(arg)
    ENDIF

    RETURN
  END FUNCTION ddlog

  !****************************************************
  FUNCTION ddexp(arg)
    REAL(8) :: arg,ddexp

    IF (arg>maxexparg) THEN
       ddexp=maxexp
    ELSE IF (arg<minexparg) THEN
       ddexp=minexp
    ELSE
       ddexp=EXP(arg)
    ENDIF

    RETURN
  END FUNCTION ddexp

  !**********************************************

END MODULE GlobalMath_loc

