*DECK DSINDG
      REAL*8 FUNCTION DSINDG (X)
C***BEGIN PROLOGUE  DSINDG
C***PURPOSE  Compute the sine of an argument in degrees.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4A
C***TYPE      REAL*8 (SINDG-S, DSINDG-D)
C***KEYWORDS  DEGREES, ELEMENTARY FUNCTIONS, FNLIB, SINE, TRIGONOMETRIC
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DSINDG(X) calculates the double precision sine for double
C precision argument X where X is in degrees.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DSINDG
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER n
C============
      REAL*8 X, RADDEG
      SAVE RADDEG
      DATA RADDEG / 0.017453292519943295769236907684886_r8/
C***FIRST EXECUTABLE STATEMENT  DSINDG
      DSINDG = SIN (RADDEG*X)
C
      IF (MOD(X,90._r8).NE.0._r8) RETURN
      N = ABS(X)/90._r8+ 0.5_r8
      N = MOD (N, 2)
      IF (N.EQ.0) DSINDG = 0._r8
      IF (N.EQ.1) DSINDG = SIGN (1.0_r8, DSINDG)
C
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
