*DECK D9GMIC
      REAL*8 FUNCTION D9GMIC (A, X, ALX)
C***BEGIN PROLOGUE  D9GMIC
C***SUBSIDIARY
C***PURPOSE  Compute the complementary incomplete Gamma function for A
C            near a negative integer and X small.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7E
C***TYPE      REAL*8 (R9GMIC-S, D9GMIC-D)
C***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, SMALL X,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Compute the complementary incomplete gamma function for A near
C a negative integer and for small X.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DLNGAM, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900720  Routine changed from user-callable to subsidiary.  (WRB)
C***END PROLOGUE  D9GMIC
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER m,k,mm1
C============
      REAL*8 A, X, ALX, ALNG, BOT, EPS, EULER, FK, FKP1, FM,
     1  S, SGNG, T, TE, D1MACH, DLNGAM
      LOGICAL FIRST
      SAVE EULER, EPS, BOT, FIRST
      DATA EULER / 0.57721566490153286060651209008240_r8/
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  D9GMIC
      IF (FIRST) THEN
         EPS = 0.5_r8*D1MACH(3)
         BOT = LOG (D1MACH(1))
      ENDIF
      FIRST = .FALSE.
C
      IF (A .GT. 0._r8) CALL XERMSG ('SLATEC', 'D9GMIC',
     +   'A MUST BE NEAR A NEGATIVE INTEGER', 2, 2)
      IF (X .LE. 0._r8) CALL XERMSG ('SLATEC', 'D9GMIC',
     +   'X MUST BE GT ZERO', 3, 2)
C
      M = -(A - 0.5_r8)
      FM = M
C
      TE = 1.0_r8
      T = 1.0_r8
      S = T
      DO 20 K=1,200
        FKP1 = K + 1
        TE = -X*TE/(FM+FKP1)
        T = TE/FKP1
        S = S + T
        IF (ABS(T).LT.EPS*S) GO TO 30
 20   CONTINUE
      CALL XERMSG ('SLATEC', 'D9GMIC',
     +   'NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION', 4, 2)
C
 30   D9GMIC = -ALX - EULER + X*S/(FM+1.0_r8)
      IF (M.EQ.0) RETURN
C
      IF (M.EQ.1) D9GMIC = -D9GMIC - 1._r8+ 1._r8/X
      IF (M.EQ.1) RETURN
C
      TE = FM
      T = 1._r8
      S = T
      MM1 = M - 1
      DO 40 K=1,MM1
        FK = K
        TE = -X*TE/FK
        T = TE/(FM-FK)
        S = S + T
        IF (ABS(T).LT.EPS*ABS(S)) GO TO 50
 40   CONTINUE
C
 50   DO 60 K=1,M
        D9GMIC = D9GMIC + 1.0_r8/K
 60   CONTINUE
C
      SGNG = 1.0_r8
      IF (MOD(M,2).EQ.1) SGNG = -1.0_r8
      ALNG = LOG(D9GMIC) - DLNGAM(FM+1._r8)
C
      D9GMIC = 0._r8
      IF (ALNG.GT.BOT) D9GMIC = SGNG * EXP(ALNG)
      IF (S.NE.0._r8) D9GMIC = D9GMIC +
     1  SIGN (EXP(-FM*ALX+LOG(ABS(S)/FM)), S)
C
      IF (D9GMIC .EQ. 0._r8.AND. S .EQ. 0._r8)
     +   CALL XERMSG ('SLATEC',
     +   'D9GMIC', 'RESULT UNDERFLOWS', 1, 1)
      RETURN
C
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
