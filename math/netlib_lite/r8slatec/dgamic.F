*DECK DGAMIC
      REAL*8 FUNCTION DGAMIC (A, X)
C***BEGIN PROLOGUE  DGAMIC
C***PURPOSE  Calculate the complementary incomplete Gamma function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7E
C***TYPE      REAL*8 (GAMIC-S, DGAMIC-D)
C***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C   Evaluate the complementary incomplete Gamma function
C
C   DGAMIC = integral from X to infinity of EXP(-T) * T**(A-1.)  .
C
C   DGAMIC is evaluated for arbitrary real values of A and for non-
C   negative values of X (even though DGAMIC is defined for X .LT.
C   0.0), except that for X = 0 and A .LE. 0.0, DGAMIC is undefined.
C
C   DGAMIC, A, and X are REAL*8.
C
C   A slight deterioration of 2 or 3 digits accuracy will occur when
C   DGAMIC is very large or very small in absolute value, because log-
C   arithmic variables are used.  Also, if the parameter A is very close
C   to a negative INTEGER (but not a negative integer), there is a loss
C   of accuracy, which is reported if the result is less than half
C   machine precision.
C
C***REFERENCES  W. Gautschi, A computational procedure for incomplete
C                 gamma functions, ACM Transactions on Mathematical
C                 Software 5, 4 (December 1979), pp. 466-481.
C               W. Gautschi, Incomplete gamma functions, Algorithm 542,
C                 ACM Transactions on Mathematical Software 5, 4
C                 (December 1979), pp. 482-489.
C***ROUTINES CALLED  D1MACH, D9GMIC, D9GMIT, D9LGIC, D9LGIT, DLGAMS,
C                    DLNGAM, XERCLR, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
C***END PROLOGUE  DGAMIC
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER izero
C============
      REAL*8 A, X, AEPS, AINTA, ALGAP1, ALNEPS, ALNGS, ALX,
     1  BOT, E, EPS, GSTAR, H, SGA, SGNG, SGNGAM, SGNGS, SQEPS, T,
     2  D1MACH, DLNGAM, D9GMIC, D9GMIT, D9LGIC, D9LGIT
      LOGICAL FIRST
      SAVE EPS, SQEPS, ALNEPS, BOT, FIRST
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DGAMIC
      IF (FIRST) THEN
         EPS = 0.5_r8*D1MACH(3)
         SQEPS = SQRT(D1MACH(4))
         ALNEPS = -LOG (D1MACH(3))
         BOT = LOG (D1MACH(1))
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LT. 0._r8)
     +   CALL XERMSG ('SLATEC', 'DGAMIC', 'X IS NEGATIVE'
     +   , 2, 2)
C
      IF (X.GT.0._r8) GO TO 20
      IF (A .LE. 0._r8) CALL XERMSG ('SLATEC', 'DGAMIC',
     +   'X = 0 AND A LE 0 SO DGAMIC IS UNDEFINED', 3, 2)
C
      DGAMIC = EXP (DLNGAM(A+1._r8) - LOG(A))
      RETURN
C
 20   ALX = LOG (X)
      SGA = 1.0_r8
      IF (A.NE.0._r8) SGA = SIGN (1.0_r8, A)
      AINTA = AINT (A + 0.5_r8*SGA)
      AEPS = A - AINTA
C
      IZERO = 0
      IF (X.GE.1.0_r8) GO TO 40
C
      IF (A.GT.0.5_r8.OR. ABS(AEPS).GT.0.001_r8) GO TO 30
      E = 2.0_r8
      IF (-AINTA.GT.1._r8) E = 2._r8*
     >   (-AINTA+2._r8)/(AINTA*AINTA-1.0_r8)
      E = E - ALX * X**(-0.001_r8)
      IF (E*ABS(AEPS).GT.EPS) GO TO 30
C
      DGAMIC = D9GMIC (A, X, ALX)
      RETURN
C
 30   CALL DLGAMS (A+1.0_r8, ALGAP1, SGNGAM)
      GSTAR = D9GMIT (A, X, ALGAP1, SGNGAM, ALX)
      IF (GSTAR.EQ.0._r8) IZERO = 1
      IF (GSTAR.NE.0._r8) ALNGS = LOG (ABS(GSTAR))
      IF (GSTAR.NE.0._r8) SGNGS = SIGN (1.0_r8, GSTAR)
      GO TO 50
C
 40   IF (A.LT.X) DGAMIC = EXP (D9LGIC(A, X, ALX))
      IF (A.LT.X) RETURN
C
      SGNGAM = 1.0_r8
      ALGAP1 = DLNGAM (A+1.0_r8)
      SGNGS = 1.0_r8
      ALNGS = D9LGIT (A, X, ALGAP1)
C
C EVALUATION OF DGAMIC(A,X) IN TERMS OF TRICOMI-S INCOMPLETE GAMMA FN.
C
 50   H = 1._r8
      IF (IZERO.EQ.1) GO TO 60
C
      T = A*ALX + ALNGS
      IF (T.GT.ALNEPS) GO TO 70
      IF (T.GT.(-ALNEPS)) H = 1.0_r8- SGNGS*EXP(T)
C
      IF (ABS(H).LT.SQEPS) CALL XERCLR
      IF (ABS(H) .LT. SQEPS) CALL XERMSG ('SLATEC', 'DGAMIC',
     +   'RESULT LT HALF PRECISION', 1, 1)
C
 60   SGNG = SIGN (1.0_r8, H) * SGA * SGNGAM
      T = LOG(ABS(H)) + ALGAP1 - LOG(ABS(A))
      IF (T.LT.BOT) CALL XERCLR
      DGAMIC = SGNG * EXP(T)
      RETURN
C
 70   SGNG = -SGNGS * SGA * SGNGAM
      T = T + ALGAP1 - LOG(ABS(A))
      IF (T.LT.BOT) CALL XERCLR
      DGAMIC = SGNG * EXP(T)
      RETURN
C
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
