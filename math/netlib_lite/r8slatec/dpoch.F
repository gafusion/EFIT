*DECK DPOCH
      REAL*8 FUNCTION DPOCH (A, X)
C***BEGIN PROLOGUE  DPOCH
C***PURPOSE  Evaluate a generalization of Pochhammer's symbol.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C1, C7A
C***TYPE      REAL*8 (POCH-S, DPOCH-D)
C***KEYWORDS  FNLIB, POCHHAMMER, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Evaluate a double precision generalization of Pochhammer's symbol
C (A)-sub-X = GAMMA(A+X)/GAMMA(A) for double precision A and X.
C For X a non-negative integer, POCH(A,X) is just Pochhammer's symbol.
C This is a preliminary version that does not handle wrong arguments
C properly and may not properly handle the case when the result is
C computed to less than half of double precision.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D9LGMC, DFAC, DGAMMA, DGAMR, DLGAMS, DLNREL, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900727  Added EXTERNAL statement.  (WRB)
C***END PROLOGUE  DPOCH
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER n,ia,i
C============
      REAL*8 A, X, ABSA, ABSAX, ALNGA, ALNGAX, AX, B, PI,
     1  SGNGA, SGNGAX, DFAC, DLNREL, D9LGMC, DGAMMA, DGAMR, DCOT
      EXTERNAL DGAMMA,DCOT
      SAVE PI
      DATA PI / 3.141592653589793238462643383279503_r8/
C***FIRST EXECUTABLE STATEMENT  DPOCH
      AX = A + X
      IF (AX.GT.0.0_r8) GO TO 30
      IF (AINT(AX).NE.AX) GO TO 30
C
      IF (A .GT. 0.0_r8.OR. AINT(A) .NE. A) CALL XERMSG ('SLATEC',
     +   'DPOCH', 'A+X IS NON-POSITIVE INTEGER BUT A IS NOT', 2, 2)
C
C WE KNOW HERE THAT BOTH A+X AND A ARE NON-POSITIVE INTEGERS.
C
      DPOCH = 1.0_r8
      IF (X.EQ.0._r8) RETURN
C
      N = X
      IF (MIN(A+X,A).LT.(-20.0_r8)) GO TO 20
C
      IA = A
      DPOCH = (-1.0_r8)**N * DFAC(-IA)/DFAC(-IA-N)
      RETURN
C
 20   DPOCH = (-1.0_r8)**N *
     1   EXP ((A-0.5_r8)*DLNREL(X/(A-1.0_r8))
     1   + X*LOG(-A+1.0_r8-X) - X + D9LGMC(-A+1.0_r8) -
     1   D9LGMC(-A-X+1._r8))
      RETURN
C
C A+X IS NOT ZERO OR A NEGATIVE INTEGER.
C
 30   DPOCH = 0.0_r8
      IF (A.LE.0.0_r8.AND. AINT(A).EQ.A) RETURN
C
      N = ABS(X)
      IF (real(N,r8).NE.X .OR. N.GT.20) GO TO 50
C
C X IS A SMALL NON-POSITIVE INTEGER, PRESUMMABLY A COMMON CASE.
C
      DPOCH = 1.0_r8
      IF (N.EQ.0) RETURN
      DO 40 I=1,N
        DPOCH = DPOCH * (A+I-1)
 40   CONTINUE
      RETURN
C
 50   ABSAX = ABS(A+X)
      ABSA = ABS(A)
      IF (MAX(ABSAX,ABSA).GT.20.0_r8) GO TO 60
      DPOCH = DGAMMA(A+X) * DGAMR(A)
      RETURN
C
 60   IF (ABS(X).GT.0.5_r8*ABSA) GO TO 70
C
C ABS(X) IS SMALL AND BOTH ABS(A+X) AND ABS(A) ARE LARGE.  THUS,
C A+X AND A MUST HAVE THE SAME SIGN.  FOR NEGATIVE A, WE USE
C GAMMA(A+X)/GAMMA(A) = GAMMA(-A+1)/GAMMA(-A-X+1) *
C SIN(PI*A)/SIN(PI*(A+X))
C
      B = A
      IF (B.LT.0.0_r8) B = -A - X + 1.0_r8
      DPOCH = EXP ((B-0.5_r8)*DLNREL(X/B) + X*LOG(B+X) - X
     1  + D9LGMC(B+X) - D9LGMC(B) )
      IF (A.LT.0.0_r8.AND. DPOCH.NE.0.0_r8) DPOCH =
     1  DPOCH/(COS(PI*X) + DCOT(PI*A)*SIN(PI*X) )
      RETURN
C
 70   CALL DLGAMS (A+X, ALNGAX, SGNGAX)
      CALL DLGAMS (A, ALNGA, SGNGA)
      DPOCH = SGNGAX * SGNGA * EXP(ALNGAX-ALNGA)
C
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
