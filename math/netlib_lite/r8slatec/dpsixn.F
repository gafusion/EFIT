*DECK DPSIXN
      REAL*8 FUNCTION DPSIXN (N)
C***BEGIN PROLOGUE  DPSIXN
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DEXINT
C***LIBRARY   SLATEC
C***TYPE      REAL*8 (PSIXN-S, DPSIXN-D)
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     This subroutine returns values of PSI(X)=derivative of log
C     GAMMA(X), X.GT.0.0 at integer arguments. A table look-up is
C     performed for N .LE. 100, and the asymptotic expansion is
C     evaluated for N.GT.100.
C
C***SEE ALSO  DEXINT
C***ROUTINES CALLED  D1MACH
C***REVISION HISTORY  (YYMMDD)
C   800501  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DPSIXN
C
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER N, K
      REAL*8 AX, B, C, FN, RFN2, TRM, S, WDTOL
      REAL*8 D1MACH
      DIMENSION B(6), C(100)
C
C             DPSIXN(N), N = 1,100
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2     C(19), C(20), C(21), C(22), C(23), C(24)/
     3 -5.77215664901532861E-01_r8,  4.22784335098467139E-01_r8,
     4  9.22784335098467139E-01_r8,  1.25611766843180047E00_r8,
     5  1.50611766843180047E00_r8,  1.70611766843180047E00_r8,
     6  1.87278433509846714E00_r8,  2.01564147795561000E00_r8,
     7  2.14064147795561000E00_r8,  2.25175258906672111E00_r8,
     8  2.35175258906672111E00_r8,  2.44266167997581202E00_r8,
     9  2.52599501330914535E00_r8,  2.60291809023222227E00_r8,
     1  2.67434666166079370E00_r8,  2.74101332832746037E00_r8,
     2  2.80351332832746037E00_r8,  2.86233685773922507E00_r8,
     3  2.91789241329478063E00_r8,  2.97052399224214905E00_r8,
     4  3.02052399224214905E00_r8,  3.06814303986119667E00_r8,
     5  3.11359758531574212E00_r8,  3.15707584618530734E00_r8/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3  3.19874251285197401E00_r8,  3.23874251285197401E00_r8,
     4  3.27720405131351247E00_r8,  3.31424108835054951E00_r8,
     5  3.34995537406483522E00_r8,  3.38443813268552488E00_r8,
     6  3.41777146601885821E00_r8,  3.45002953053498724E00_r8,
     7  3.48127953053498724E00_r8,  3.51158256083801755E00_r8,
     8  3.54099432554389990E00_r8,  3.56956575411532847E00_r8,
     9  3.59734353189310625E00_r8,  3.62437055892013327E00_r8,
     1  3.65068634839381748E00_r8,  3.67632737403484313E00_r8,
     2  3.70132737403484313E00_r8,  3.72571761793728215E00_r8,
     3  3.74952714174680596E00_r8,  3.77278295570029433E00_r8,
     4  3.79551022842756706E00_r8,  3.81773245064978928E00_r8,
     5  3.83947158108457189E00_r8,  3.86074817682925274E00_r8/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2     C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/
     3  3.88158151016258607E00_r8,  3.90198967342789220E00_r8,
     4  3.92198967342789220E00_r8,  3.94159751656514710E00_r8,
     5  3.96082828579591633E00_r8,  3.97969621032421822E00_r8,
     6  3.99821472884273674E00_r8,  4.01639654702455492E00_r8,
     7  4.03425368988169777E00_r8,  4.05179754953082058E00_r8,
     8  4.06903892884116541E00_r8,  4.08598808138353829E00_r8,
     9  4.10265474805020496E00_r8,  4.11904819067315578E00_r8,
     1  4.13517722293122029E00_r8,  4.15105023880423617E00_r8,
     2  4.16667523880423617E00_r8,  4.18205985418885155E00_r8,
     3  4.19721136934036670E00_r8,  4.21213674247469506E00_r8,
     4  4.22684262482763624E00_r8,  4.24133537845082464E00_r8,
     5  4.25562109273653893E00_r8,  4.26970559977879245E00_r8/
      DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80),
     1     C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88),
     2     C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/
     3  4.28359448866768134E00_r8,  4.29729311880466764E00_r8,
     4  4.31080663231818115E00_r8,  4.32413996565151449E00_r8,
     5  4.33729786038835659E00_r8,  4.35028487337536958E00_r8,
     6  4.36310538619588240E00_r8,  4.37576361404398366E00_r8,
     7  4.38826361404398366E00_r8,  4.40060929305632934E00_r8,
     8  4.41280441500754886E00_r8,  4.42485260777863319E00_r8,
     9  4.43675736968339510E00_r8,  4.44852207556574804E00_r8,
     1  4.46014998254249223E00_r8,  4.47164423541605544E00_r8,
     2  4.48300787177969181E00_r8,  4.49424382683587158E00_r8,
     3  4.50535493794698269E00_r8,  4.51634394893599368E00_r8,
     4  4.52721351415338499E00_r8,  4.53796620232542800E00_r8,
     5  4.54860450019776842E00_r8,  4.55913081598724211E00_r8/
      DATA C(97), C(98), C(99), C(100)/
     1  4.56954748265390877E00_r8,  4.57985676100442424E00_r8,
     2  4.59006084263707730E00_r8,  4.60016185273808740E00_r8/
C             COEFFICIENTS OF ASYMPTOTIC EXPANSION
      DATA B(1), B(2), B(3), B(4), B(5), B(6)/
     1  8.33333333333333333E-02_r8, -8.33333333333333333E-03_r8,
     2  3.96825396825396825E-03_r8, -4.16666666666666666E-03_r8,
     3  7.57575757575757576E-03_r8, -2.10927960927960928E-02_r8/
C
C***FIRST EXECUTABLE STATEMENT  DPSIXN
      IF (N.GT.100) GO TO 10
      DPSIXN = C(N)
      RETURN
   10 CONTINUE
      WDTOL = MAX(D1MACH(4),1.0E-18_r8)
      FN = N
      AX = 1.0_r8
      S = -0.5_r8/FN
      IF (ABS(S).LE.WDTOL) GO TO 30
      RFN2 = 1.0_r8/(FN*FN)
      DO 20 K=1,6
        AX = AX*RFN2
        TRM = -B(K)*AX
        IF (ABS(TRM).LT.WDTOL) GO TO 30
        S = S + TRM
   20 CONTINUE
   30 CONTINUE
      DPSIXN = S + LOG(FN)
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
