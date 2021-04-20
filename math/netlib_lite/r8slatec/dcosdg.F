*DECK DCOSDG
      REAL*8 FUNCTION DCOSDG (X)
C***BEGIN PROLOGUE  DCOSDG
C***PURPOSE  Compute the cosine of an argument in degrees.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4A
C***TYPE      REAL*8 (COSDG-S, DCOSDG-D)
C***KEYWORDS  COSINE, DEGREES, ELEMENTARY FUNCTIONS, FNLIB,
C             TRIGONOMETRIC
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DCOSDG(X) calculates the double precision trigonometric cosine
C for double precision argument X in units of degrees.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DCOSDG
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER n
C============
      REAL*8 X, RADDEG
      SAVE RADDEG
      DATA RADDEG / 0.017453292519943295769236907684886_r8/
C***FIRST EXECUTABLE STATEMENT  DCOSDG
      DCOSDG = COS (RADDEG*X)
C
      IF (MOD(X,90._r8).NE.0._r8) RETURN
      N = ABS(X)/90._r8+ 0.5_r8
      N = MOD (N, 2)
      IF (N.EQ.0) DCOSDG = SIGN (1.0_r8, DCOSDG)
      IF (N.EQ.1) DCOSDG = 0.0_r8
C
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
