*DECK DACOSH
      REAL*8 FUNCTION DACOSH (X)
C***BEGIN PROLOGUE  DACOSH
C***PURPOSE  Compute the arc hyperbolic cosine.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4C
C***TYPE      REAL*8 (ACOSH-S, DACOSH-D, CACOSH-C)
C***KEYWORDS  ACOSH, ARC HYPERBOLIC COSINE, ELEMENTARY FUNCTIONS, FNLIB,
C             INVERSE HYPERBOLIC COSINE
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DACOSH(X) calculates the double precision arc hyperbolic cosine for
C double precision argument X.  The result is returned on the
C positive branch.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DACOSH
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      REAL*8 X, DLN2, XMAX,  D1MACH
      SAVE DLN2, XMAX
      DATA DLN2 / 0.69314718055994530941723212145818_r8/
      DATA XMAX / 0._r8/
C***FIRST EXECUTABLE STATEMENT  DACOSH
      IF (XMAX.EQ.0._r8) XMAX = 1.0_r8/SQRT(D1MACH(3))
C
      IF (X .LT. 1._r8) CALL XERMSG ('SLATEC', 'DACOSH',
     +   'X LESS THAN 1', 1, 2)
C
      IF (X.LT.XMAX) DACOSH = LOG (X+SQRT(X*X-1.0_r8))
      IF (X.GE.XMAX) DACOSH = DLN2 + LOG(X)
C
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
