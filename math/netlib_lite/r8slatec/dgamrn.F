*DECK DGAMRN
      REAL*8 FUNCTION DGAMRN (X)
C***BEGIN PROLOGUE  DGAMRN
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBSKIN
C***LIBRARY   SLATEC
C***TYPE      REAL*8 (GAMRN-S, DGAMRN-D)
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Abstract   * A Double Precision Routine *
C         DGAMRN computes the GAMMA function ratio GAMMA(X)/GAMMA(X+0.5)
C         for real X.gt.0. If X.ge.XMIN, an asymptotic expansion is
C         evaluated. If X.lt.XMIN, an integer is added to X to form a
C         new value of X.ge.XMIN and the asymptotic expansion is eval-
C         uated for this new value of X. Successive application of the
C         recurrence relation
C
C                      W(X)=W(X+1)*(1+0.5/X)
C
C         reduces the argument to its original value. XMIN and comp-
C         utational tolerances are computed as a function of the number
C         of digits carried in a word by calls to J1MACH and D1MACH.
C         However, the computational accuracy is limited to the max-
C         imum of unit roundoff (=D1MACH(4)) and 1.0D-18 since critical
C         constants are given to only 18 digits.
C
C         Input      X is Double Precision
C           X      - Argument, X.gt.0.0D0
C
C         Output      DGAMRN is REAL*8
C           DGAMRN  - Ratio  GAMMA(X)/GAMMA(X+0.5)
C
C***SEE ALSO  DBSKIN
C***REFERENCES  Y. L. Luke, The Special Functions and Their
C                 Approximations, Vol. 1, Math In Sci. And
C                 Eng. Series 53, Academic Press, New York, 1969,
C                 pp. 34-35.
C***ROUTINES CALLED  D1MACH, J1MACH
C***REVISION HISTORY  (YYMMDD)
C   820601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C   920520  Added REFERENCES section.  (WRB)
C***END PROLOGUE  DGAMRN
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER I, I1M11, K, MX, NX
      INTEGER J1MACH
      REAL*8 FLN, GR, RLN, S, TOL, TRM, X, XDMY, XINC, XM,
     * XMIN, XP, XSQ
      REAL*8 D1MACH
      DIMENSION GR(12)
      SAVE GR
C
      DATA GR(1), GR(2), GR(3), GR(4), GR(5), GR(6), GR(7), GR(8),
     * GR(9), GR(10), GR(11), GR(12) /1.00000000000000000E00_r8,
     * -1.56250000000000000E-02_r8,2.56347656250000000E-03_r8,
     * -1.27983093261718750E-03_r8,1.34351104497909546E-03_r8,
     * -2.43289663922041655E-03_r8,6.75423753364157164E-03_r8,
     * -2.66369606131178216E-02_r8,1.41527455519564332E-01_r8,
     * -9.74384543032201613E-01_r8,8.43686251229783675E00_r8,
     * -8.97258321640552515E01_r8/
C
C***FIRST EXECUTABLE STATEMENT  DGAMRN
      NX = INT(X)
      TOL = MAX(D1MACH(4),1.0E-18_r8)
      I1M11 = J1MACH(14)
      RLN = D1MACH(5)*I1M11
      FLN = MIN(RLN,20.0_r8)
      FLN = MAX(FLN,3.0_r8)
      FLN = FLN - 3.0_r8
      XM = 2.0_r8+ FLN*(0.2366_r8+0.01723_r8*FLN)
      MX = INT(XM) + 1
      XMIN = MX
      XDMY = X - 0.25_r8
      XINC = 0.0_r8
      IF (X.GE.XMIN) GO TO 10
      XINC = XMIN - NX
      XDMY = XDMY + XINC
   10 CONTINUE
      S = 1.0_r8
      IF (XDMY*TOL.GT.1.0_r8) GO TO 30
      XSQ = 1.0_r8/(XDMY*XDMY)
      XP = XSQ
      DO 20 K=2,12
        TRM = GR(K)*XP
        IF (ABS(TRM).LT.TOL) GO TO 30
        S = S + TRM
        XP = XP*XSQ
   20 CONTINUE
   30 CONTINUE
      S = S/SQRT(XDMY)
      IF (XINC.NE.0.0_r8) GO TO 40
      DGAMRN = S
      RETURN
   40 CONTINUE
      NX = INT(XINC)
      XP = 0.0_r8
      DO 50 I=1,NX
        S = S*(1.0_r8+0.5_r8/(X+XP))
        XP = XP + 1.0_r8
   50 CONTINUE
      DGAMRN = S
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
