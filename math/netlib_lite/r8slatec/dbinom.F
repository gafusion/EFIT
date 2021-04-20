*DECK DBINOM
      REAL*8 FUNCTION DBINOM (N, M)
C***BEGIN PROLOGUE  DBINOM
C***PURPOSE  Compute the binomial coefficients.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C1
C***TYPE      REAL*8 (BINOM-S, DBINOM-D)
C***KEYWORDS  BINOMIAL COEFFICIENTS, FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBINOM(N,M) calculates the double precision binomial coefficient
C for integer arguments N and M.  The result is (N!)/((M!)(N-M)!).
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, D9LGMC, DLNREL, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DBINOM
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER m,n,k,i
C============
      REAL*8 CORR, FINTMX, SQ2PIL, XK, XN, XNK, D9LGMC,
     1  DLNREL, D1MACH, BILNMX
      LOGICAL FIRST
      SAVE SQ2PIL, BILNMX, FINTMX, FIRST
      DATA SQ2PIL / 0.91893853320467274178032973640562_r8/
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DBINOM
      IF (FIRST) THEN
         BILNMX = LOG(D1MACH(2)) - 0.0001_r8
         FINTMX = 0.9_r8/D1MACH(3)
      ENDIF
      FIRST = .FALSE.
C
      IF (N .LT. 0 .OR. M .LT. 0) CALL XERMSG ('SLATEC', 'DBINOM',
     +   'N OR M LT ZERO', 1, 2)
      IF (N .LT. M) CALL XERMSG ('SLATEC', 'DBINOM', 'N LT M', 2, 2)
C
      K = MIN (M, N-M)
      IF (K.GT.20) GO TO 30
      IF (K*LOG(1.0_r8*max(N,1)).GT.BILNMX) GO TO 30
C
      DBINOM = 1.0_r8
      IF (K.EQ.0) RETURN
      DO 20 I=1,K
        XN = N - I + 1
        XK = I
        DBINOM = DBINOM * (XN/XK)
 20   CONTINUE
C
      IF (DBINOM.LT.FINTMX) DBINOM = AINT (DBINOM+0.5_r8)
      RETURN
C
C IF K.LT.9, APPROX IS NOT VALID AND ANSWER IS CLOSE TO THE OVERFLOW LIM
 30   IF (K .LT. 9) CALL XERMSG ('SLATEC', 'DBINOM',
     +   'RESULT OVERFLOWS BECAUSE N AND/OR M TOO BIG', 3, 2)
C
      XN = N + 1
      XK = K + 1
      XNK = N - K + 1
C
      CORR = D9LGMC(XN) - D9LGMC(XK) - D9LGMC(XNK)
      DBINOM = XK*LOG(XNK/XK) - XN*DLNREL(-(XK-1.0_r8)/XN)
     1  -0.5_r8*LOG(XN*XNK/XK) + 1.0_r8- SQ2PIL + CORR
C
      IF (DBINOM .GT. BILNMX) CALL XERMSG ('SLATEC', 'DBINOM',
     +   'RESULT OVERFLOWS BECAUSE N AND/OR M TOO BIG', 3, 2)
C
      DBINOM = EXP (DBINOM)
      IF (DBINOM.LT.FINTMX) DBINOM = AINT (DBINOM+0.5_r8)
C
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
