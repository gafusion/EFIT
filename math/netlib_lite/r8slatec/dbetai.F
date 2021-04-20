*DECK DBETAI
      REAL*8 FUNCTION DBETAI (X, PIN, QIN)
C***BEGIN PROLOGUE  DBETAI
C***PURPOSE  Calculate the incomplete Beta function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7F
C***TYPE      REAL*8 (BETAI-S, DBETAI-D)
C***KEYWORDS  FNLIB, INCOMPLETE BETA FUNCTION, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C   DBETAI calculates the REAL*8 incomplete beta function.
C
C   The incomplete beta function ratio is the probability that a
C   random variable from a beta distribution having parameters PIN and
C   QIN will be less than or equal to X.
C
C     -- Input Arguments -- All arguments are REAL*8.
C   X      upper limit of integration.  X must be in (0,1) inclusive.
C   PIN    first beta distribution parameter.  PIN must be .GT. 0.0.
C   QIN    second beta distribution parameter.  QIN must be .GT. 0.0.
C
C***REFERENCES  Nancy E. Bosten and E. L. Battiste, Remark on Algorithm
C                 179, Communications of the ACM 17, 3 (March 1974),
C                 pp. 156.
C***ROUTINES CALLED  D1MACH, DLBETA, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
C***END PROLOGUE  DBETAI
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER n,i,ib
C============
      REAL*8 X, PIN, QIN, ALNEPS, ALNSML, C, EPS, FINSUM, P,
     1  PS, Q, SML, TERM, XB, XI, Y, D1MACH, DLBETA, P1
      LOGICAL FIRST
      SAVE EPS, ALNEPS, SML, ALNSML, FIRST
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DBETAI
      IF (FIRST) THEN
         EPS = D1MACH(3)
         ALNEPS = LOG (EPS)
         SML = D1MACH(1)
         ALNSML = LOG (SML)
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LT. 0._r8.OR. X .GT. 1._r8)
     +   CALL XERMSG ('SLATEC', 'DBETAI',
     +   'X IS NOT IN THE RANGE (0,1)', 1, 2)
      IF (PIN .LE. 0._r8.OR. QIN .LE. 0._r8)
     +   CALL XERMSG ('SLATEC',
     +   'DBETAI', 'P AND/OR Q IS LE ZERO', 2, 2)
C
      Y = X
      P = PIN
      Q = QIN
      IF (Q.LE.P .AND. X.LT.0.8_r8) GO TO 20
      IF (X.LT.0.2_r8) GO TO 20
      Y = 1.0_r8- Y
      P = QIN
      Q = PIN
C
 20   IF ((P+Q)*Y/(P+1._r8).LT.EPS) GO TO 80
C
C EVALUATE THE INFINITE SUM FIRST.  TERM WILL EQUAL
C Y**P/BETA(PS,P) * (1.-PS)-SUB-I * Y**I / FAC(I) .
C
      PS = Q - AINT(Q)
      IF (PS.EQ.0._r8) PS = 1.0_r8
      XB = P*LOG(Y) - DLBETA(PS,P) - LOG(P)
      DBETAI = 0.0_r8
      IF (XB.LT.ALNSML) GO TO 40
C
      DBETAI = EXP (XB)
      TERM = DBETAI*P
      IF (PS.EQ.1.0_r8) GO TO 40
      N = MAX (ALNEPS/LOG(Y), 4.0_r8)
      DO 30 I=1,N
        XI = I
        TERM = TERM * (XI-PS)*Y/XI
        DBETAI = DBETAI + TERM/(P+XI)
 30   CONTINUE
C
C NOW EVALUATE THE FINITE SUM, MAYBE.
C
 40   IF (Q.LE.1.0_r8) GO TO 70
C
      XB = P*LOG(Y) + Q*LOG(1.0_r8-Y) - DLBETA(P,Q) - LOG(Q)
      IB = MAX (XB/ALNSML, 0.0_r8)
      TERM = EXP(XB - IB*ALNSML)
      C = 1.0_r8/(1._r8-Y)
      P1 = Q*C/(P+Q-1._r8)
C
      FINSUM = 0.0_r8
      N = Q
      IF (Q.EQ.real(N,r8)) N = N - 1
      DO 50 I=1,N
        IF (P1.LE.1.0_r8.AND. TERM/EPS.LE.FINSUM) GO TO 60
        XI = I
        TERM = (Q-XI+1.0_r8)*C*TERM/(P+Q-XI)
C
        IF (TERM.GT.1.0_r8) IB = IB - 1
        IF (TERM.GT.1.0_r8) TERM = TERM*SML
C
        IF (IB.EQ.0) FINSUM = FINSUM + TERM
 50   CONTINUE
C
 60   DBETAI = DBETAI + FINSUM
 70   IF (Y.NE.X .OR. P.NE.PIN) DBETAI = 1.0_r8- DBETAI
      DBETAI = MAX (MIN (DBETAI, 1.0_r8), 0.0_r8)
      RETURN
C
 80   DBETAI = 0.0_r8
      XB = P*LOG(MAX(Y,SML)) - LOG(P) - DLBETA(P,Q)
      IF (XB.GT.ALNSML .AND. Y.NE.0.0_r8) DBETAI = EXP(XB)
      IF (Y.NE.X .OR. P.NE.PIN) DBETAI = 1.0_r8- DBETAI
C
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
! 22Jun2000 fgtok
