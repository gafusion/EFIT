*DECK D9LGIC
      REAL*8 FUNCTION D9LGIC (A, X, ALX)
C***BEGIN PROLOGUE  D9LGIC
C***SUBSIDIARY
C***PURPOSE  Compute the log complementary incomplete Gamma function
C            for large X and for A .LE. X.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7E
C***TYPE      REAL*8 (R9LGIC-S, D9LGIC-D)
C***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, LARGE X,
C             LOGARITHM, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Compute the log complementary incomplete gamma function for large X
C and for A .LE. X.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900720  Routine changed from user-callable to subsidiary.  (WRB)
C***END PROLOGUE  D9LGIC
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER k
C============
      REAL*8 A, X, ALX, EPS, FK, P, R, S, T, XMA, XPA, D1MACH
      SAVE EPS
      DATA EPS / 0._r8/
C***FIRST EXECUTABLE STATEMENT  D9LGIC
      IF (EPS.EQ.0._r8) EPS = 0.5_r8*D1MACH(3)
C
      XPA = X + 1.0_r8- A
      XMA = X - 1._r8- A
C
      R = 0._r8
      P = 1._r8
      S = P
      DO 10 K=1,300
        FK = K
        T = FK*(A-FK)*(1._r8+R)
        R = -T/((XMA+2._r8*FK)*(XPA+2._r8*FK)+T)
        P = R*P
        S = S + P
        IF (ABS(P).LT.EPS*S) GO TO 20
 10   CONTINUE
      CALL XERMSG ('SLATEC', 'D9LGIC',
     +   'NO CONVERGENCE IN 300 TERMS OF CONTINUED FRACTION', 1, 2)
C
 20   D9LGIC = A*ALX - X + LOG(S/XPA)
C
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
