*DECK DBKISR
      SUBROUTINE DBKISR (X, N, SUM, IERR)
C***BEGIN PROLOGUE  DBKISR
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBSKIN
C***LIBRARY   SLATEC
C***TYPE      REAL*8 (BKISR-S, DBKISR-D)
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     DBKISR computes repeated integrals of the K0 Bessel function
C     by the series for N=0,1, and 2.
C
C***SEE ALSO  DBSKIN
C***ROUTINES CALLED  D1MACH, DPSIXN
C***REVISION HISTORY  (YYMMDD)
C   820601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DBKISR
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER I, IERR, K, KK, KKN, K1, N, NP
      REAL*8 AK, ATOL, BK, C, FK, FN, HX, HXS, POL, PR, SUM,
     * TKP, TOL, TRM, X, XLN
      REAL*8 DPSIXN, D1MACH
      DIMENSION C(2)
      SAVE C
C
      DATA C(1), C(2) /1.57079632679489662E00_r8,1.0_r8/
C***FIRST EXECUTABLE STATEMENT  DBKISR
      IERR=0
      TOL = MAX(D1MACH(4),1.0E-18_r8)
      IF (X.LT.TOL) GO TO 50
      PR = 1.0_r8
      POL = 0.0_r8
      IF (N.EQ.0) GO TO 20
      DO 10 I=1,N
        POL = -POL*X + C(I)
        PR = PR*X/I
   10 CONTINUE
   20 CONTINUE
      HX = X*0.5_r8
      HXS = HX*HX
      XLN = LOG(HX)
      NP = N + 1
      TKP = 3.0_r8
      FK = 2.0_r8
      FN = N
      BK = 4.0_r8
      AK = 2.0_r8/((FN+1.0_r8)*(FN+2.0_r8))
      SUM = AK*(DPSIXN(N+3)-DPSIXN(3)+DPSIXN(2)-XLN)
      ATOL = SUM*TOL*0.75_r8
      DO 30 K=2,20
        AK = AK*(HXS/BK)*((TKP+1.0_r8)/(TKP+FN+1.0_r8))*
     +      (TKP/(TKP+FN))
        K1 = K + 1
        KK = K1 + K
        KKN = KK + N
        TRM = (DPSIXN(K1)+DPSIXN(KKN)-DPSIXN(KK)-XLN)*AK
        SUM = SUM + TRM
        IF (ABS(TRM).LE.ATOL) GO TO 40
        TKP = TKP + 2.0_r8
        BK = BK + TKP
        FK = FK + 1.0_r8
   30 CONTINUE
      GO TO 80
   40 CONTINUE
      SUM = (SUM*HXS+DPSIXN(NP)-XLN)*PR
      IF (N.EQ.1) SUM = -SUM
      SUM = POL + SUM
      RETURN
C-----------------------------------------------------------------------
C     SMALL X CASE, X.LT.WORD TOLERANCE
C-----------------------------------------------------------------------
   50 CONTINUE
      IF (N.GT.0) GO TO 60
      HX = X*0.5_r8
      SUM = DPSIXN(1) - LOG(HX)
      RETURN
   60 CONTINUE
      SUM = C(N)
      RETURN
   80 CONTINUE
      IERR=2
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
