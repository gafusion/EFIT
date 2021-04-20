*DECK DPOCH1
      REAL*8 FUNCTION DPOCH1 (A, X)
C***BEGIN PROLOGUE  DPOCH1
C***PURPOSE  Calculate a generalization of Pochhammer's symbol starting
C            from first order.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C1, C7A
C***TYPE      REAL*8 (POCH1-S, DPOCH1-D)
C***KEYWORDS  FIRST ORDER, FNLIB, POCHHAMMER, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Evaluate a double precision generalization of Pochhammer's symbol
C for double precision A and X for special situations that require
C especially accurate values when X is small in
C        POCH1(A,X) = (POCH(A,X)-1)/X
C                   = (GAMMA(A+X)/GAMMA(A) - 1.0)/X .
C This specification is particularly suited for stably computing
C expressions such as
C        (GAMMA(A+X)/GAMMA(A) - GAMMA(B+X)/GAMMA(B))/X
C             = POCH1(A,X) - POCH1(B,X)
C Note that POCH1(A,0.0) = PSI(A)
C
C When ABS(X) is so small that substantial cancellation will occur if
C the straightforward formula is used, we use an expansion due
C to Fields and discussed by Y. L. Luke, The Special Functions and Their
C Approximations, Vol. 1, Academic Press, 1969, page 34.
C
C The ratio POCH(A,X) = GAMMA(A+X)/GAMMA(A) is written by Luke as
C        (A+(X-1)/2)**X * polynomial in (A+(X-1)/2)**(-2) .
C In order to maintain significance in POCH1, we write for positive a
C        (A+(X-1)/2)**X = EXP(X*LOG(A+(X-1)/2)) = EXP(Q)
C                       = 1.0 + Q*EXPREL(Q) .
C Likewise the polynomial is written
C        POLY = 1.0 + X*POLY1(A,X) .
C Thus,
C        POCH1(A,X) = (POCH(A,X) - 1) / X
C                   = EXPREL(Q)*(Q/X + Q*POLY1(A,X)) + POLY1(A,X)
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DCOT, DEXPRL, DPOCH, DPSI, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900727  Added EXTERNAL statement.  (WRB)
C***END PROLOGUE  DPOCH1
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER incr,nterms,k,j,ndx,ii,i
C============
      REAL*8 A, X, ABSA, ABSX, ALNEPS, ALNVAR, B, BERN(20),
     1  BINV, BP, GBERN(21), GBK, PI, POLY1, Q, RHO, SINPXX, SINPX2,
     2  SQTBIG, TERM, TRIG, VAR, VAR2, D1MACH, DPSI, DEXPRL, DCOT, DPOCH
      LOGICAL FIRST
      EXTERNAL DCOT
      SAVE BERN, PI, SQTBIG, ALNEPS, FIRST
      DATA BERN  (  1) / +.833333333333333333333333333333333E-1_r8/
      DATA BERN  (  2) / -.138888888888888888888888888888888E-2_r8/
      DATA BERN  (  3) / +.330687830687830687830687830687830E-4_r8/
      DATA BERN  (  4) / -.826719576719576719576719576719576E-6_r8/
      DATA BERN  (  5) / +.208767569878680989792100903212014E-7_r8/
      DATA BERN  (  6) / -.528419013868749318484768220217955E-9_r8/
      DATA BERN  (  7) / +.133825365306846788328269809751291E-10_r8/
      DATA BERN  (  8) / -.338968029632258286683019539124944E-12_r8/
      DATA BERN  (  9) / +.858606205627784456413590545042562E-14_r8/
      DATA BERN  ( 10) / -.217486869855806187304151642386591E-15_r8/
      DATA BERN  ( 11) / +.550900282836022951520265260890225E-17_r8/
      DATA BERN  ( 12) / -.139544646858125233407076862640635E-18_r8/
      DATA BERN  ( 13) / +.353470703962946747169322997780379E-20_r8/
      DATA BERN  ( 14) / -.895351742703754685040261131811274E-22_r8/
      DATA BERN  ( 15) / +.226795245233768306031095073886816E-23_r8/
      DATA BERN  ( 16) / -.574472439520264523834847971943400E-24_r8/
      DATA BERN  ( 17) / +.145517247561486490186626486727132E-26_r8/
      DATA BERN  ( 18) / -.368599494066531017818178247990866E-28_r8/
      DATA BERN  ( 19) / +.933673425709504467203255515278562E-30_r8/
      DATA BERN  ( 20) / -.236502241570062993455963519636983E-31_r8/
      DATA PI / 3.141592653589793238462643383279503_r8/
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DPOCH1
      IF (FIRST) THEN
         SQTBIG = 1.0_r8/SQRT(24.0_r8*D1MACH(1))
         ALNEPS = LOG(D1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      IF (X.EQ.0.0_r8) DPOCH1 = DPSI(A)
      IF (X.EQ.0.0_r8) RETURN
C
      ABSX = ABS(X)
      ABSA = ABS(A)
      IF (ABSX.GT.0.1_r8*ABSA) GO TO 70
      IF (ABSX*LOG(MAX(ABSA,2.0_r8)).GT.0.1_r8) GO TO 70
C
      BP = A
      IF (A.LT.(-0.5_r8)) BP = 1.0_r8- A - X
      INCR = 0
      IF (BP.LT.10.0_r8) INCR = 11.0_r8- BP
      B = BP + INCR
C
      VAR = B + 0.5_r8*(X-1.0_r8)
      ALNVAR = LOG(VAR)
      Q = X*ALNVAR
C
      POLY1 = 0.0_r8
      IF (VAR.GE.SQTBIG) GO TO 40
      VAR2 = (1.0_r8/VAR)**2
C
      RHO = 0.5_r8*(X+1.0_r8)
      GBERN(1) = 1.0_r8
      GBERN(2) = -RHO/12.0_r8
      TERM = VAR2
      POLY1 = GBERN(2)*TERM
C
      NTERMS = -0.5_r8*ALNEPS/ALNVAR + 1.0_r8
      IF (NTERMS .GT. 20) CALL XERMSG ('SLATEC', 'DPOCH1',
     +   'NTERMS IS TOO BIG, MAYBE D1MACH(3) IS BAD', 1, 2)
      IF (NTERMS.LT.2) GO TO 40
C
      DO 30 K=2,NTERMS
        GBK = 0.0_r8
        DO 20 J=1,K
          NDX = K - J + 1
          GBK = GBK + BERN(NDX)*GBERN(J)
 20     CONTINUE
        GBERN(K+1) = -RHO*GBK/K
C
        TERM = TERM * (2*K-2-X)*(2*K-1-X)*VAR2
        POLY1 = POLY1 + GBERN(K+1)*TERM
 30   CONTINUE
C
 40   POLY1 = (X-1.0_r8)*POLY1
      DPOCH1 = DEXPRL(Q)*(ALNVAR+Q*POLY1) + POLY1
C
      IF (INCR.EQ.0) GO TO 60
C
C WE HAVE DPOCH1(B,X), BUT BP IS SMALL, SO WE USE BACKWARDS RECURSION
C TO OBTAIN DPOCH1(BP,X).
C
      DO 50 II=1,INCR
        I = INCR - II
        BINV = 1.0_r8/(BP+I)
        DPOCH1 = (DPOCH1 - BINV) / (1.0_r8+ X*BINV)
 50   CONTINUE
C
 60   IF (BP.EQ.A) RETURN
C
C WE HAVE DPOCH1(BP,X), BUT A IS LT -0.5.  WE THEREFORE USE A REFLECTION
C FORMULA TO OBTAIN DPOCH1(A,X).
C
      SINPXX = SIN(PI*X)/X
      SINPX2 = SIN(0.5_r8*PI*X)
      TRIG = SINPXX*DCOT(PI*B) - 2.0_r8*SINPX2*(SINPX2/X)
C
      DPOCH1 = TRIG + (1.0_r8+ X*TRIG)*DPOCH1
      RETURN
C
 70   DPOCH1 = (DPOCH(A,X) - 1.0_r8) / X
      RETURN
C
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
