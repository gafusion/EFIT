*DECK DHKSEQ
      SUBROUTINE DHKSEQ (X, M, H, IERR)
C***BEGIN PROLOGUE  DHKSEQ
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBSKIN
C***LIBRARY   SLATEC
C***TYPE      REAL*8 (HKSEQ-S, DHKSEQ-D)
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C   DHKSEQ is an adaptation of subroutine DPSIFN described in the
C   reference below.  DHKSEQ generates the sequence
C   H(K,X) = (-X)**(K+1)*(PSI(K,X) PSI(K,X+0.5))/GAMMA(K+1), for
C            K=0,...,M.
C
C***SEE ALSO  DBSKIN
C***REFERENCES  D. E. Amos, A portable Fortran subroutine for
C                 derivatives of the Psi function, Algorithm 610, ACM
C                 Transactions on Mathematical Software 9, 4 (1983),
C                 pp. 494-502.
C***ROUTINES CALLED  D1MACH, J1MACH
C***REVISION HISTORY  (YYMMDD)
C   820601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
C***END PROLOGUE  DHKSEQ
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER I, IERR, J, K, M, MX, NX
      INTEGER J1MACH
      REAL*8 B, FK, FLN, FN, FNP, H, HRX, RLN, RXSQ, R1M5, S,
     * SLOPE, T, TK, TRM, TRMH, TRMR, TST, U, V, WDTOL, X, XDMY, XH,
     * XINC, XM, XMIN, YINT
      REAL*8 D1MACH
      DIMENSION B(22), TRM(22), TRMR(25), TRMH(25), U(25), V(25), H(*)
      SAVE B
C-----------------------------------------------------------------------
C             SCALED BERNOULLI NUMBERS 2.0*B(2K)*(1-2**(-2K))
C-----------------------------------------------------------------------
      DATA B(1), B(2), B(3), B(4), B(5), B(6), B(7), B(8), B(9), B(10),
     * B(11), B(12), B(13), B(14), B(15), B(16), B(17), B(18), B(19),
     * B(20), B(21), B(22) /1.00000000000000000E00_r8,
     * -5.00000000000000000E-01_r8,2.50000000000000000E-01_r8,
     * -6.25000000000000000E-02_r8,4.68750000000000000E-02_r8,
     * -6.64062500000000000E-02_r8,1.51367187500000000E-01_r8,
     * -5.06103515625000000E-01_r8,2.33319091796875000E00_r8,
     * -1.41840972900390625E01_r8,1.09941936492919922E02_r8,
     * -1.05824747562408447E03_r8,1.23842434241771698E04_r8,
     * -1.73160495905935764E05_r8,2.85103429084961116E06_r8,
     * -5.45964619322445132E07_r8,1.20316174668075304E09_r8,
     * -3.02326315271452307E10_r8,8.59229286072319606E11_r8,
     * -2.74233104097776039E13_r8,9.76664637943633248E14_r8,
     * -3.85931586838450360E16_r8/
C
C***FIRST EXECUTABLE STATEMENT  DHKSEQ
      IERR=0
      WDTOL = MAX(D1MACH(4),1.0E-18_r8)
      FN = M - 1
      FNP = FN + 1.0_r8
C-----------------------------------------------------------------------
C     COMPUTE XMIN
C-----------------------------------------------------------------------
      R1M5 = D1MACH(5)
      RLN = R1M5*J1MACH(14)
      RLN = MIN(RLN,18.06_r8)
      FLN = MAX(RLN,3.0_r8) - 3.0_r8
      YINT = 3.50_r8+ 0.40_r8*FLN
      SLOPE = 0.21_r8+
     >   FLN*(0.0006038_r8*FLN+0.008677_r8)
      XM = YINT + SLOPE*FN
      MX = INT(XM) + 1
      XMIN = MX
C-----------------------------------------------------------------------
C     GENERATE H(M-1,XDMY)*XDMY**(M) BY THE ASYMPTOTIC EXPANSION
C-----------------------------------------------------------------------
      XDMY = X
      XINC = 0.0_r8
      IF (X.GE.XMIN) GO TO 10
      NX = INT(X)
      XINC = XMIN - NX
      XDMY = X + XINC
   10 CONTINUE
      RXSQ = 1.0_r8/(XDMY*XDMY)
      HRX = 0.5_r8/XDMY
      TST = 0.5_r8*WDTOL
      T = FNP*HRX
C-----------------------------------------------------------------------
C     INITIALIZE COEFFICIENT ARRAY
C-----------------------------------------------------------------------
      S = T*B(3)
      IF (ABS(S).LT.TST) GO TO 30
      TK = 2.0_r8
      DO 20 K=4,22
        T = T*((TK+FN+1.0_r8)/
     >      (TK+1.0_r8))*((TK+FN)/(TK+2.0_r8))*RXSQ
        TRM(K) = T*B(K)
        IF (ABS(TRM(K)).LT.TST) GO TO 30
        S = S + TRM(K)
        TK = TK + 2.0_r8
   20 CONTINUE
      GO TO 110
   30 CONTINUE
      H(M) = S + 0.5_r8
      IF (M.EQ.1) GO TO 70
C-----------------------------------------------------------------------
C     GENERATE LOWER DERIVATIVES, I.LT.M-1
C-----------------------------------------------------------------------
      DO 60 I=2,M
        FNP = FN
        FN = FN - 1.0_r8
        S = FNP*HRX*B(3)
        IF (ABS(S).LT.TST) GO TO 50
        FK = FNP + 3.0_r8
        DO 40 K=4,22
          TRM(K) = TRM(K)*FNP/FK
          IF (ABS(TRM(K)).LT.TST) GO TO 50
          S = S + TRM(K)
          FK = FK + 2.0_r8
   40   CONTINUE
        GO TO 110
   50   CONTINUE
        MX = M - I + 1
        H(MX) = S + 0.5_r8
   60 CONTINUE
   70 CONTINUE
      IF (XINC.EQ.0.0_r8) RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD FROM XDMY TO X
C-----------------------------------------------------------------------
      XH = X + 0.5_r8
      S = 0.0_r8
      NX = INT(XINC)
      DO 80 I=1,NX
        TRMR(I) = X/(X+NX-I)
        U(I) = TRMR(I)
        TRMH(I) = X/(XH+NX-I)
        V(I) = TRMH(I)
        S = S + U(I) - V(I)
   80 CONTINUE
      MX = NX + 1
      TRMR(MX) = X/XDMY
      U(MX) = TRMR(MX)
      H(1) = H(1)*TRMR(MX) + S
      IF (M.EQ.1) RETURN
      DO 100 J=2,M
        S = 0.0_r8
        DO 90 I=1,NX
          TRMR(I) = TRMR(I)*U(I)
          TRMH(I) = TRMH(I)*V(I)
          S = S + TRMR(I) - TRMH(I)
   90   CONTINUE
        TRMR(MX) = TRMR(MX)*U(MX)
        H(J) = H(J)*TRMR(MX) + S
  100 CONTINUE
      RETURN
  110 CONTINUE
      IERR=2
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
