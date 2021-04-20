*DECK DEI
      REAL*8 FUNCTION DEI (X)
C***BEGIN PROLOGUE  DEI
C***PURPOSE  Compute the exponential integral Ei(X).
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C5
C***TYPE      REAL*8 (EI-S, DEI-D)
C***KEYWORDS  EI FUNCTION, EXPONENTIAL INTEGRAL, FNLIB,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DEI calculates the double precision exponential integral, Ei(X), for
C positive double precision argument X and the Cauchy principal value
C for negative X.  If principal values are used everywhere, then, for
C all X,
C
C    Ei(X) = -E1(-X)
C or
C    E1(X) = -Ei(-X).
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DE1
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   891115  Modified prologue description.  (WRB)
C   891115  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DEI
      IMPLICIT NONE
      REAL*8 X, DE1
C***FIRST EXECUTABLE STATEMENT  DEI
      DEI = -DE1(-X)
C
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
