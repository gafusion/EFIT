 
*DECK DCHU
      REAL*8 FUNCTION DCHU (A, B, X)
C***BEGIN PROLOGUE  DCHU
C***PURPOSE  Compute the logarithmic confluent hypergeometric function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C11
C***TYPE      REAL*8 (CHU-S, DCHU-D)
C***KEYWORDS  FNLIB, LOGARITHMIC CONFLUENT HYPERGEOMETRIC FUNCTION,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DCHU(A,B,X) calculates the double precision logarithmic confluent
C hypergeometric function U(A,B,X) for double precision arguments
C A, B, and X.
C
C This routine is not valid when 1+A-B is close to zero if X is small.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, D9CHU, DEXPRL, DGAMMA, DGAMR, DPOCH,
C                    DPOCH1, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900727  Added EXTERNAL statement.  (WRB)
C***END PROLOGUE  DCHU
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER n,m,i,istrt
C============
      REAL*8 A, B, X, AINTB, ALNX, A0, BEPS, B0, C0, EPS,
     1  FACTOR, GAMRI1, GAMRNI, PCH1AI, PCH1I, PI, POCHAI, SUM, T,
     2  XEPS1, XI, XI1, XN, XTOEPS,  D1MACH, DPOCH, DGAMMA, DGAMR,
     3  DPOCH1, DEXPRL, D9CHU
      EXTERNAL DGAMMA
      SAVE PI, EPS
      DATA PI / 3.141592653589793238462643383279503_r8/
      DATA EPS / 0.0_r8/
C***FIRST EXECUTABLE STATEMENT  DCHU
      IF (EPS.EQ.0.0_r8) EPS = D1MACH(3)
C
      IF (X .EQ. 0.0_r8) CALL XERMSG ('SLATEC', 'DCHU',
     +   'X IS ZERO SO DCHU IS INFINITE', 1, 2)
      IF (X .LT. 0.0_r8) CALL XERMSG ('SLATEC', 'DCHU',
     +   'X IS NEGATIVE, USE CCHU', 2, 2)
C
      IF (MAX(ABS(A),1.0_r8)*
     1   MAX(ABS(1.0_r8+A-B),1.0_r8).LT.
     1   0.99_r8*ABS(X)) GO TO 120
C
C THE ASCENDING SERIES WILL BE USED, BECAUSE THE DESCENDING RATIONAL
C APPROXIMATION (WHICH IS BASED ON THE ASYMPTOTIC SERIES) IS UNSTABLE.
C
      IF (ABS(1.0_r8+A-B) .LT. SQRT(EPS))
     +   CALL XERMSG ('SLATEC', 'DCHU',
     +   'ALGORITHMIS BAD WHEN 1+A-B IS NEAR ZERO FOR SMALL X', 10, 2)
C
      IF (B.GE.0.0_r8) AINTB = AINT(B+0.5_r8)
      IF (B.LT.0.0_r8) AINTB = AINT(B-0.5_r8)
      BEPS = B - AINTB
      N = AINTB
C
      ALNX = LOG(X)
      XTOEPS = EXP (-BEPS*ALNX)
C
C EVALUATE THE FINITE SUM.     -----------------------------------------
C
      IF (N.GE.1) GO TO 40
C
C CONSIDER THE CASE B .LT. 1.0 FIRST.
C
      SUM = 1.0_r8
      IF (N.EQ.0) GO TO 30
C
      T = 1.0_r8
      M = -N
      DO 20 I=1,M
        XI1 = I - 1
        T = T*(A+XI1)*X/((B+XI1)*(XI1+1.0_r8))
        SUM = SUM + T
 20   CONTINUE
C
 30   SUM = DPOCH(1.0_r8+A-B, -A)*SUM
      GO TO 70
C
C NOW CONSIDER THE CASE B .GE. 1.0.
C
 40   SUM = 0.0_r8
      M = N - 2
      IF (M.LT.0) GO TO 70
      T = 1.0_r8
      SUM = 1.0_r8
      IF (M.EQ.0) GO TO 60
C
      DO 50 I=1,M
        XI = I
        T = T * (A-B+XI)*X/((1.0_r8-B+XI)*XI)
        SUM = SUM + T
 50   CONTINUE
C
 60   SUM = DGAMMA(B-1.0_r8) * DGAMR(A) * X**(1-N) * XTOEPS * SUM
C
C NEXT EVALUATE THE INFINITE SUM.     ----------------------------------
C
 70   ISTRT = 0
      IF (N.LT.1) ISTRT = 1 - N
      XI = ISTRT
C
      FACTOR = (-1.0_r8)**N * DGAMR(1.0_r8+A-B) * X**ISTRT
      IF (BEPS.NE.0.0_r8) FACTOR = FACTOR * BEPS*PI/SIN(BEPS*PI)
C
      POCHAI = DPOCH (A, XI)
      GAMRI1 = DGAMR (XI+1.0_r8)
      GAMRNI = DGAMR (AINTB+XI)
      B0 = FACTOR * DPOCH(A,XI-BEPS) * GAMRNI *
     >   DGAMR(XI+1.0_r8-BEPS)
C
      IF (ABS(XTOEPS-1.0_r8).GT.0.5_r8) GO TO 90
C
C X**(-BEPS) IS CLOSE TO 1.0D0, SO WE MUST BE CAREFUL IN EVALUATING THE
C DIFFERENCES.
C
      PCH1AI = DPOCH1 (A+XI, -BEPS)
      PCH1I = DPOCH1 (XI+1.0_r8-BEPS, BEPS)
      C0 = FACTOR * POCHAI * GAMRNI * GAMRI1 * (
     1  -DPOCH1(B+XI,-BEPS) + PCH1AI - PCH1I + BEPS*PCH1AI*PCH1I)
C
C XEPS1 = (1.0 - X**(-BEPS))/BEPS = (X**(-BEPS) - 1.0)/(-BEPS)
      XEPS1 = ALNX*DEXPRL(-BEPS*ALNX)
C
      DCHU = SUM + C0 + XEPS1*B0
      XN = N
      DO 80 I=1,1000
        XI = ISTRT + I
        XI1 = ISTRT + I - 1
        B0 = (A+XI1-BEPS)*B0*X/((XN+XI1)*(XI-BEPS))
        C0 = (A+XI1)*C0*X/((B+XI1)*XI)
     1    - ((A-1.0_r8)*(XN+2._r8*XI-1.0_r8) +
     2     XI*(XI-BEPS)) * B0
     3     / (XI*(B+XI1)*(A+XI1-BEPS))
        T = C0 + XEPS1*B0
        DCHU = DCHU + T
        IF (ABS(T).LT.EPS*ABS(DCHU)) GO TO 130
 80   CONTINUE
      CALL XERMSG ('SLATEC', 'DCHU',
     +   'NO CONVERGENCE IN 1000 TERMS OF THE ASCENDING SERIES', 3, 2)
C
C X**(-BEPS) IS VERY DIFFERENT FROM 1.0, SO THE STRAIGHTFORWARD
C FORMULATION IS STABLE.
C
 90   A0 = FACTOR * POCHAI * DGAMR(B+XI) * GAMRI1 / BEPS
      B0 = XTOEPS * B0 / BEPS
C
      DCHU = SUM + A0 - B0
      DO 100 I=1,1000
        XI = ISTRT + I
        XI1 = ISTRT + I - 1
        A0 = (A+XI1)*A0*X/((B+XI1)*XI)
        B0 = (A+XI1-BEPS)*B0*X/((AINTB+XI1)*(XI-BEPS))
        T = A0 - B0
        DCHU = DCHU + T
        IF (ABS(T).LT.EPS*ABS(DCHU)) GO TO 130
 100  CONTINUE
      CALL XERMSG ('SLATEC', 'DCHU',
     +   'NO CONVERGENCE IN 1000 TERMS OF THE ASCENDING SERIES', 3, 2)
C
C USE LUKE-S RATIONAL APPROXIMATION IN THE ASYMPTOTIC REGION.
C
 120  DCHU = X**(-A) * D9CHU(A,B,X)
C
 130  RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
