*DECK DLGAMS
      SUBROUTINE DLGAMS (X, DLGAM, SGNGAM)
C***BEGIN PROLOGUE  DLGAMS
C***PURPOSE  Compute the logarithm of the absolute value of the Gamma
C            function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7A
C***TYPE      REAL*8 (ALGAMS-S, DLGAMS-D)
C***KEYWORDS  ABSOLUTE VALUE OF THE LOGARITHM OF THE GAMMA FUNCTION,
C             FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DLGAMS(X,DLGAM,SGNGAM) calculates the double precision natural
C logarithm of the absolute value of the Gamma function for
C double precision argument X and stores the result in double
C precision argument DLGAM.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DLNGAM
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DLGAMS
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      integer INT
      REAL*8 X, DLGAM, SGNGAM, DLNGAM
C***FIRST EXECUTABLE STATEMENT  DLGAMS
      DLGAM = DLNGAM(X)
      SGNGAM = 1.0_r8
      IF (X.GT.0._r8) RETURN
C
      INT = MOD (-AINT(X), 2.0_r8) + 0.1_r8
      IF (INT.EQ.0) SGNGAM = -1.0_r8
C
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
