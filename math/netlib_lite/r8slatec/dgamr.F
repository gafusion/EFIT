*DECK DGAMR
      REAL*8 FUNCTION DGAMR (X)
C***BEGIN PROLOGUE  DGAMR
C***PURPOSE  Compute the reciprocal of the Gamma function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7A
C***TYPE      REAL*8 (GAMR-S, DGAMR-D, CGAMR-C)
C***KEYWORDS  FNLIB, RECIPROCAL GAMMA FUNCTION, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DGAMR(X) calculates the double precision reciprocal of the
C complete Gamma function for double precision argument X.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DGAMMA, DLGAMS, XERCLR, XGETF, XSETF
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900727  Added EXTERNAL statement.  (WRB)
C***END PROLOGUE  DGAMR
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER irold
C============
      REAL*8 X, ALNGX, SGNGX, DGAMMA
      EXTERNAL DGAMMA
C***FIRST EXECUTABLE STATEMENT  DGAMR
      DGAMR = 0.0_r8
      IF (X.LE.0.0_r8.AND. AINT(X).EQ.X) RETURN
C
      CALL XGETF (IROLD)
      CALL XSETF (1)
      IF (ABS(X).GT.10.0_r8) GO TO 10
      DGAMR = 1.0_r8/DGAMMA(X)
      CALL XERCLR
      CALL XSETF (IROLD)
      RETURN
C
 10   CALL DLGAMS (X, ALNGX, SGNGX)
      CALL XERCLR
      CALL XSETF (IROLD)
      DGAMR = SGNGX * EXP(-ALNGX)
      RETURN
C
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
