*DECK DGAMLN
      REAL*8 FUNCTION DGAMLN (Z, IERR)
C***BEGIN PROLOGUE  DGAMLN
C***SUBSIDIARY
C***PURPOSE  Compute the logarithm of the Gamma function
C***LIBRARY   SLATEC
C***CATEGORY  C7A
C***TYPE      REAL*8 (GAMLN-S, DGAMLN-D)
C***KEYWORDS  LOGARITHM OF GAMMA FUNCTION
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C               **** A REAL*8 ROUTINE ****
C         DGAMLN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION FOR
C         Z.GT.0.  THE ASYMPTOTIC EXPANSION IS USED TO GENERATE VALUES
C         GREATER THAN ZMIN WHICH ARE ADJUSTED BY THE RECURSION
C         G(Z+1)=Z*G(Z) FOR Z.LE.ZMIN.  THE FUNCTION WAS MADE AS
C         PORTABLE AS POSSIBLE BY COMPUTING ZMIN FROM THE NUMBER OF BASE
C         10 DIGITS IN A WORD, RLN=MAX(-ALOG10(R1MACH(4)),0.5E-18)
C         LIMITED TO 18 DIGITS OF (RELATIVE) ACCURACY.
C
C         SINCE INTEGER ARGUMENTS ARE COMMON, A TABLE LOOK UP ON 100
C         VALUES IS USED FOR SPEED OF EXECUTION.
C
C     DESCRIPTION OF ARGUMENTS
C
C         INPUT      Z IS D0UBLE PRECISION
C           Z      - ARGUMENT, Z.GT.0.0D0
C
C         OUTPUT      DGAMLN IS REAL*8
C           DGAMLN  - NATURAL LOG OF THE GAMMA FUNCTION AT Z.NE.0.0D0
C           IERR    - ERROR FLAG
C                     IERR=0, NORMAL RETURN, COMPUTATION COMPLETED
C                     IERR=1, Z.LE.0.0D0,    NO COMPUTATION
C
C
C***REFERENCES  COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C***ROUTINES CALLED  D1MACH, J1MACH
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   830501  REVISION DATE from Version 3.2
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   920128  Category corrected.  (WRB)
C   921215  DGAMLN defined for Z negative.  (WRB)
C***END PROLOGUE  DGAMLN
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      REAL*8 CF, CON, FLN, FZ, GLN, RLN, S, TLG, TRM, TST,
     * T1, WDTOL, Z, ZDMY, ZINC, ZM, ZMIN, ZP, ZSQ, D1MACH
      INTEGER I, IERR, I1M, K, MZ, NZ, J1MACH
      DIMENSION CF(22), GLN(100)
C           LNGAMMA(N), N=1,100
      DATA GLN(1), GLN(2), GLN(3), GLN(4), GLN(5), GLN(6), GLN(7),
     1     GLN(8), GLN(9), GLN(10), GLN(11), GLN(12), GLN(13), GLN(14),
     2     GLN(15), GLN(16), GLN(17), GLN(18), GLN(19), GLN(20),
     3     GLN(21), GLN(22)/
     4  0.00000000000000000E00_r8,  0.00000000000000000E00_r8,
     5  6.93147180559945309E-01_r8,  1.79175946922805500E00_r8,
     6  3.17805383034794562E00_r8,  4.78749174278204599E00_r8,
     7  6.57925121201010100E00_r8,  8.52516136106541430E00_r8,
     8  1.06046029027452502E01_r8,  1.28018274800814696E01_r8,
     9  1.51044125730755153E01_r8,  1.75023078458738858E01_r8,
     A  1.99872144956618861E01_r8,  2.25521638531234229E01_r8,
     B  2.51912211827386815E01_r8,  2.78992713838408916E01_r8,
     C  3.06718601060806728E01_r8,  3.35050734501368889E01_r8,
     D  3.63954452080330536E01_r8,  3.93398841871994940E01_r8,
     E  4.23356164607534850E01_r8,  4.53801388984769080E01_r8/
      DATA GLN(23), GLN(24), GLN(25), GLN(26), GLN(27), GLN(28),
     1     GLN(29), GLN(30), GLN(31), GLN(32), GLN(33), GLN(34),
     2     GLN(35), GLN(36), GLN(37), GLN(38), GLN(39), GLN(40),
     3     GLN(41), GLN(42), GLN(43), GLN(44)/
     4  4.84711813518352239E01_r8,  5.16066755677643736E01_r8,
     5  5.47847293981123192E01_r8,  5.80036052229805199E01_r8,
     6  6.12617017610020020E01_r8,  6.45575386270063311E01_r8,
     7  6.78897431371815350E01_r8,  7.12570389671680090E01_r8,
     8  7.46582363488301644E01_r8,  7.80922235533153106E01_r8,
     9  8.15579594561150372E01_r8,  8.50544670175815174E01_r8,
     A  8.85808275421976788E01_r8,  9.21361756036870925E01_r8,
     B  9.57196945421432025E01_r8,  9.93306124547874269E01_r8,
     C  1.02968198614513813E02_r8,  1.06631760260643459E02_r8,
     D  1.10320639714757395E02_r8,  1.14034211781461703E02_r8,
     E  1.17771881399745072E02_r8,  1.21533081515438634E02_r8/
      DATA GLN(45), GLN(46), GLN(47), GLN(48), GLN(49), GLN(50),
     1     GLN(51), GLN(52), GLN(53), GLN(54), GLN(55), GLN(56),
     2     GLN(57), GLN(58), GLN(59), GLN(60), GLN(61), GLN(62),
     3     GLN(63), GLN(64), GLN(65), GLN(66)/
     4  1.25317271149356895E02_r8,  1.29123933639127215E02_r8,
     5  1.32952575035616310E02_r8,  1.36802722637326368E02_r8,
     6  1.40673923648234259E02_r8,  1.44565743946344886E02_r8,
     7  1.48477766951773032E02_r8,  1.52409592584497358E02_r8,
     8  1.56360836303078785E02_r8,  1.60331128216630907E02_r8,
     9  1.64320112263195181E02_r8,  1.68327445448427652E02_r8,
     A  1.72352797139162802E02_r8,  1.76395848406997352E02_r8,
     B  1.80456291417543771E02_r8,  1.84533828861449491E02_r8,
     C  1.88628173423671591E02_r8,  1.92739047287844902E02_r8,
     D  1.96866181672889994E02_r8,  2.01009316399281527E02_r8,
     E  2.05168199482641199E02_r8,  2.09342586752536836E02_r8/
      DATA GLN(67), GLN(68), GLN(69), GLN(70), GLN(71), GLN(72),
     1     GLN(73), GLN(74), GLN(75), GLN(76), GLN(77), GLN(78),
     2     GLN(79), GLN(80), GLN(81), GLN(82), GLN(83), GLN(84),
     3     GLN(85), GLN(86), GLN(87), GLN(88)/
     4  2.13532241494563261E02_r8,  2.17736934113954227E02_r8,
     5  2.21956441819130334E02_r8,  2.26190548323727593E02_r8,
     6  2.30439043565776952E02_r8,  2.34701723442818268E02_r8,
     7  2.38978389561834323E02_r8,  2.43268849002982714E02_r8,
     8  2.47572914096186884E02_r8,  2.51890402209723194E02_r8,
     9  2.56221135550009525E02_r8,  2.60564940971863209E02_r8,
     A  2.64921649798552801E02_r8,  2.69291097651019823E02_r8,
     B  2.73673124285693704E02_r8,  2.78067573440366143E02_r8,
     C  2.82474292687630396E02_r8,  2.86893133295426994E02_r8,
     D  2.91323950094270308E02_r8,  2.95766601350760624E02_r8,
     E  3.00220948647014132E02_r8,  3.04686856765668715E02_r8/
      DATA GLN(89), GLN(90), GLN(91), GLN(92), GLN(93), GLN(94),
     1     GLN(95), GLN(96), GLN(97), GLN(98), GLN(99), GLN(100)/
     2  3.09164193580146922E02_r8,  3.13652829949879062E02_r8,
     3  3.18152639620209327E02_r8,  3.22663499126726177E02_r8,
     4  3.27185287703775217E02_r8,  3.31717887196928473E02_r8,
     5  3.36261181979198477E02_r8,  3.40815058870799018E02_r8,
     6  3.45379407062266854E02_r8,  3.49954118040770237E02_r8,
     7  3.54539085519440809E02_r8,  3.59134205369575399E02_r8/
C             COEFFICIENTS OF ASYMPTOTIC EXPANSION
      DATA CF(1), CF(2), CF(3), CF(4), CF(5), CF(6), CF(7), CF(8),
     1     CF(9), CF(10), CF(11), CF(12), CF(13), CF(14), CF(15),
     2     CF(16), CF(17), CF(18), CF(19), CF(20), CF(21), CF(22)/
     3  8.33333333333333333E-02_r8, -2.77777777777777778E-03_r8,
     4  7.93650793650793651E-04_r8, -5.95238095238095238E-04_r8,
     5  8.41750841750841751E-04_r8, -1.91752691752691753E-03_r8,
     6  6.41025641025641026E-03_r8, -2.95506535947712418E-02_r8,
     7  1.79644372368830573E-01_r8, -1.39243221690590112E00_r8,
     8  1.34028640441683920E01_r8, -1.56848284626002017E02_r8,
     9  2.19310333333333333E03_r8, -3.61087712537249894E04_r8,
     A  6.91472268851313067E05_r8, -1.52382215394074162E07_r8,
     B  3.82900751391414141E08_r8, -1.08822660357843911E10_r8,
     C  3.47320283765002252E11_r8, -1.23696021422692745E13_r8,
     D  4.88788064793079335E14_r8, -2.13203339609193739E16_r8/
C
C             LN(2*PI)
      DATA CON                    /  1.83787706640934548E00_r8/
C
C***FIRST EXECUTABLE STATEMENT  DGAMLN
      IERR=0
      IF (Z.LE.0.0_r8) GO TO 70
      IF (Z.GT.101.0_r8) GO TO 10
      NZ = Z
      FZ = Z - NZ
      IF (FZ.GT.0.0_r8) GO TO 10
      IF (NZ.GT.100) GO TO 10
      DGAMLN = GLN(NZ)
      RETURN
   10 CONTINUE
      WDTOL = D1MACH(4)
      WDTOL = MAX(WDTOL,0.5E-18_r8)
      I1M = J1MACH(14)
      RLN = D1MACH(5)*I1M
      FLN = MIN(RLN,20.0_r8)
      FLN = MAX(FLN,3.0_r8)
      FLN = FLN - 3.0_r8
      ZM = 1.8000_r8+ 0.3875_r8*FLN
      MZ = ZM + 1
      ZMIN = MZ
      ZDMY = Z
      ZINC = 0.0_r8
      IF (Z.GE.ZMIN) GO TO 20
      ZINC = ZMIN - NZ
      ZDMY = Z + ZINC
   20 CONTINUE
      ZP = 1.0_r8/ZDMY
      T1 = CF(1)*ZP
      S = T1
      IF (ZP.LT.WDTOL) GO TO 40
      ZSQ = ZP*ZP
      TST = T1*WDTOL
      DO 30 K=2,22
        ZP = ZP*ZSQ
        TRM = CF(K)*ZP
        IF (ABS(TRM).LT.TST) GO TO 40
        S = S + TRM
   30 CONTINUE
   40 CONTINUE
      IF (ZINC.NE.0.0_r8) GO TO 50
      TLG = LOG(Z)
      DGAMLN = Z*(TLG-1.0_r8) + 0.5_r8*(CON-TLG) + S
      RETURN
   50 CONTINUE
      ZP = 1.0_r8
      NZ = ZINC
      DO 60 I=1,NZ
        ZP = ZP*(Z+(I-1))
   60 CONTINUE
      TLG = LOG(ZDMY)
      DGAMLN = ZDMY*(TLG-1.0_r8) -
     >   LOG(ZP) + 0.5_r8*(CON-TLG) + S
      RETURN
C
C
   70 CONTINUE
      DGAMLN = D1MACH(2)
      IERR=1
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
