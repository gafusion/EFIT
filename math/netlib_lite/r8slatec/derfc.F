*DECK DERFC
      REAL*8 FUNCTION DERFC (X)
C***BEGIN PROLOGUE  DERFC
C***PURPOSE  Compute the complementary error function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C8A, L5A1E
C***TYPE      REAL*8 (ERFC-S, DERFC-D)
C***KEYWORDS  COMPLEMENTARY ERROR FUNCTION, ERFC, FNLIB,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DERFC(X) calculates the double precision complementary error function
C for double precision argument X.
C
C Series for ERF        on the interval  0.          to  1.00000E+00
C                                        with weighted Error   1.28E-32
C                                         log weighted Error  31.89
C                               significant figures required  31.05
C                                    decimal places required  32.55
C
C Series for ERC2       on the interval  2.50000E-01 to  1.00000E+00
C                                        with weighted Error   2.67E-32
C                                         log weighted Error  31.57
C                               significant figures required  30.31
C                                    decimal places required  32.42
C
C Series for ERFC       on the interval  0.          to  2.50000E-01
C                                        with weighted error   1.53E-31
C                                         log weighted error  30.82
C                               significant figures required  29.47
C                                    decimal places required  31.70
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920618  Removed space from variable names.  (RWC, WRB)
C***END PROLOGUE  DERFC
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER nterf,nterfc,nterc2,initds
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 eta
C============
      REAL*8 X, ERFCS(21), ERFCCS(59), ERC2CS(49), SQEPS,
     1  SQRTPI, XMAX, TXMAX, XSML, Y, D1MACH, DCSEVL
      LOGICAL FIRST
      SAVE ERFCS, ERC2CS, ERFCCS, SQRTPI, NTERF,
     1 NTERFC, NTERC2, XSML, XMAX, SQEPS, FIRST
      DATA ERFCS(  1) / -.49046121234691808039984544033376E-1_r8/
      DATA ERFCS(  2) / -.14226120510371364237824741899631E+0_r8/
      DATA ERFCS(  3) / +.10035582187599795575754676712933E-1_r8/
      DATA ERFCS(  4) / -.57687646997674847650827025509167E-3_r8/
      DATA ERFCS(  5) / +.27419931252196061034422160791471E-4_r8/
      DATA ERFCS(  6) / -.11043175507344507604135381295905E-5_r8/
      DATA ERFCS(  7) / +.38488755420345036949961311498174E-7_r8/
      DATA ERFCS(  8) / -.11808582533875466969631751801581E-8_r8/
      DATA ERFCS(  9) / +.32334215826050909646402930953354E-10_r8/
      DATA ERFCS( 10) / -.79910159470045487581607374708595E-12_r8/
      DATA ERFCS( 11) / +.17990725113961455611967245486634E-13_r8/
      DATA ERFCS( 12) / -.37186354878186926382316828209493E-15_r8/
      DATA ERFCS( 13) / +.71035990037142529711689908394666E-17_r8/
      DATA ERFCS( 14) / -.12612455119155225832495424853333E-18_r8/
      DATA ERFCS( 15) / +.20916406941769294369170500266666E-20_r8/
      DATA ERFCS( 16) / -.32539731029314072982364160000000E-22_r8/
      DATA ERFCS( 17) / +.47668672097976748332373333333333E-24_r8/
      DATA ERFCS( 18) / -.65980120782851343155199999999999E-26_r8/
      DATA ERFCS( 19) / +.86550114699637626197333333333333E-28_r8/
      DATA ERFCS( 20) / -.10788925177498064213333333333333E-29_r8/
      DATA ERFCS( 21) / +.12811883993017002666666666666666E-31_r8/
      DATA ERC2CS(  1) / -.6960134660230950112739150826197E-1_r8/
      DATA ERC2CS(  2) / -.4110133936262089348982212084666E-1_r8/
      DATA ERC2CS(  3) / +.3914495866689626881561143705244E-2_r8/
      DATA ERC2CS(  4) / -.4906395650548979161280935450774E-3_r8/
      DATA ERC2CS(  5) / +.7157479001377036380760894141825E-4_r8/
      DATA ERC2CS(  6) / -.1153071634131232833808232847912E-4_r8/
      DATA ERC2CS(  7) / +.1994670590201997635052314867709E-5_r8/
      DATA ERC2CS(  8) / -.3642666471599222873936118430711E-6_r8/
      DATA ERC2CS(  9) / +.6944372610005012589931277214633E-7_r8/
      DATA ERC2CS( 10) / -.1371220902104366019534605141210E-7_r8/
      DATA ERC2CS( 11) / +.2788389661007137131963860348087E-8_r8/
      DATA ERC2CS( 12) / -.5814164724331161551864791050316E-9_r8/
      DATA ERC2CS( 13) / +.1238920491752753181180168817950E-9_r8/
      DATA ERC2CS( 14) / -.2690639145306743432390424937889E-10_r8/
      DATA ERC2CS( 15) / +.5942614350847910982444709683840E-11_r8/
      DATA ERC2CS( 16) / -.1332386735758119579287754420570E-11_r8/
      DATA ERC2CS( 17) / +.3028046806177132017173697243304E-12_r8/
      DATA ERC2CS( 18) / -.6966648814941032588795867588954E-13_r8/
      DATA ERC2CS( 19) / +.1620854541053922969812893227628E-13_r8/
      DATA ERC2CS( 20) / -.3809934465250491999876913057729E-14_r8/
      DATA ERC2CS( 21) / +.9040487815978831149368971012975E-15_r8/
      DATA ERC2CS( 22) / -.2164006195089607347809812047003E-15_r8/
      DATA ERC2CS( 23) / +.5222102233995854984607980244172E-16_r8/
      DATA ERC2CS( 24) / -.1269729602364555336372415527780E-16_r8/
      DATA ERC2CS( 25) / +.3109145504276197583836227412951E-17_r8/
      DATA ERC2CS( 26) / -.7663762920320385524009566714811E-18_r8/
      DATA ERC2CS( 27) / +.1900819251362745202536929733290E-18_r8/
      DATA ERC2CS( 28) / -.4742207279069039545225655999965E-19_r8/
      DATA ERC2CS( 29) / +.1189649200076528382880683078451E-19_r8/
      DATA ERC2CS( 30) / -.3000035590325780256845271313066E-20_r8/
      DATA ERC2CS( 31) / +.7602993453043246173019385277098E-21_r8/
      DATA ERC2CS( 32) / -.1935909447606872881569811049130E-21_r8/
      DATA ERC2CS( 33) / +.4951399124773337881000042386773E-22_r8/
      DATA ERC2CS( 34) / -.1271807481336371879608621989888E-22_r8/
      DATA ERC2CS( 35) / +.3280049600469513043315841652053E-23_r8/
      DATA ERC2CS( 36) / -.8492320176822896568924792422399E-24_r8/
      DATA ERC2CS( 37) / +.2206917892807560223519879987199E-24_r8/
      DATA ERC2CS( 38) / -.5755617245696528498312819507199E-25_r8/
      DATA ERC2CS( 39) / +.1506191533639234250354144051199E-25_r8/
      DATA ERC2CS( 40) / -.3954502959018796953104285695999E-26_r8/
      DATA ERC2CS( 41) / +.1041529704151500979984645051733E-26_r8/
      DATA ERC2CS( 42) / -.2751487795278765079450178901333E-27_r8/
      DATA ERC2CS( 43) / +.7290058205497557408997703680000E-28_r8/
      DATA ERC2CS( 44) / -.1936939645915947804077501098666E-28_r8/
      DATA ERC2CS( 45) / +.5160357112051487298370054826666E-29_r8/
      DATA ERC2CS( 46) / -.1378419322193094099389644800000E-29_r8/
      DATA ERC2CS( 47) / +.3691326793107069042251093333333E-30_r8/
      DATA ERC2CS( 48) / -.9909389590624365420653226666666E-31_r8/
      DATA ERC2CS( 49) / +.2666491705195388413323946666666E-31_r8/
      DATA ERFCCS(  1) / +.715179310202924774503697709496E-1_r8/
      DATA ERFCCS(  2) / -.265324343376067157558893386681E-1_r8/
      DATA ERFCCS(  3) / +.171115397792085588332699194606E-2_r8/
      DATA ERFCCS(  4) / -.163751663458517884163746404749E-3_r8/
      DATA ERFCCS(  5) / +.198712935005520364995974806758E-4_r8/
      DATA ERFCCS(  6) / -.284371241276655508750175183152E-5_r8/
      DATA ERFCCS(  7) / +.460616130896313036969379968464E-6_r8/
      DATA ERFCCS(  8) / -.822775302587920842057766536366E-7_r8/
      DATA ERFCCS(  9) / +.159214187277090112989358340826E-7_r8/
      DATA ERFCCS( 10) / -.329507136225284321486631665072E-8_r8/
      DATA ERFCCS( 11) / +.722343976040055546581261153890E-9_r8/
      DATA ERFCCS( 12) / -.166485581339872959344695966886E-9_r8/
      DATA ERFCCS( 13) / +.401039258823766482077671768814E-10_r8/
      DATA ERFCCS( 14) / -.100481621442573113272170176283E-10_r8/
      DATA ERFCCS( 15) / +.260827591330033380859341009439E-11_r8/
      DATA ERFCCS( 16) / -.699111056040402486557697812476E-12_r8/
      DATA ERFCCS( 17) / +.192949233326170708624205749803E-12_r8/
      DATA ERFCCS( 18) / -.547013118875433106490125085271E-13_r8/
      DATA ERFCCS( 19) / +.158966330976269744839084032762E-13_r8/
      DATA ERFCCS( 20) / -.472689398019755483920369584290E-14_r8/
      DATA ERFCCS( 21) / +.143587337678498478672873997840E-14_r8/
      DATA ERFCCS( 22) / -.444951056181735839417250062829E-15_r8/
      DATA ERFCCS( 23) / +.140481088476823343737305537466E-15_r8/
      DATA ERFCCS( 24) / -.451381838776421089625963281623E-16_r8/
      DATA ERFCCS( 25) / +.147452154104513307787018713262E-16_r8/
      DATA ERFCCS( 26) / -.489262140694577615436841552532E-17_r8/
      DATA ERFCCS( 27) / +.164761214141064673895301522827E-17_r8/
      DATA ERFCCS( 28) / -.562681717632940809299928521323E-18_r8/
      DATA ERFCCS( 29) / +.194744338223207851429197867821E-18_r8/
      DATA ERFCCS( 30) / -.682630564294842072956664144723E-19_r8/
      DATA ERFCCS( 31) / +.242198888729864924018301125438E-19_r8/
      DATA ERFCCS( 32) / -.869341413350307042563800861857E-20_r8/
      DATA ERFCCS( 33) / +.315518034622808557122363401262E-20_r8/
      DATA ERFCCS( 34) / -.115737232404960874261239486742E-20_r8/
      DATA ERFCCS( 35) / +.428894716160565394623737097442E-21_r8/
      DATA ERFCCS( 36) / -.160503074205761685005737770964E-21_r8/
      DATA ERFCCS( 37) / +.606329875745380264495069923027E-22_r8/
      DATA ERFCCS( 38) / -.231140425169795849098840801367E-22_r8/
      DATA ERFCCS( 39) / +.888877854066188552554702955697E-23_r8/
      DATA ERFCCS( 40) / -.344726057665137652230718495566E-23_r8/
      DATA ERFCCS( 41) / +.134786546020696506827582774181E-23_r8/
      DATA ERFCCS( 42) / -.531179407112502173645873201807E-24_r8/
      DATA ERFCCS( 43) / +.210934105861978316828954734537E-24_r8/
      DATA ERFCCS( 44) / -.843836558792378911598133256738E-25_r8/
      DATA ERFCCS( 45) / +.339998252494520890627359576337E-25_r8/
      DATA ERFCCS( 46) / -.137945238807324209002238377110E-25_r8/
      DATA ERFCCS( 47) / +.563449031183325261513392634811E-26_r8/
      DATA ERFCCS( 48) / -.231649043447706544823427752700E-26_r8/
      DATA ERFCCS( 49) / +.958446284460181015263158381226E-27_r8/
      DATA ERFCCS( 50) / -.399072288033010972624224850193E-27_r8/
      DATA ERFCCS( 51) / +.167212922594447736017228709669E-27_r8/
      DATA ERFCCS( 52) / -.704599152276601385638803782587E-28_r8/
      DATA ERFCCS( 53) / +.297976840286420635412357989444E-28_r8/
      DATA ERFCCS( 54) / -.126252246646061929722422632994E-28_r8/
      DATA ERFCCS( 55) / +.539543870454248793985299653154E-29_r8/
      DATA ERFCCS( 56) / -.238099288253145918675346190062E-29_r8/
      DATA ERFCCS( 57) / +.109905283010276157359726683750E-29_r8/
      DATA ERFCCS( 58) / -.486771374164496572732518677435E-30_r8/
      DATA ERFCCS( 59) / +.152587726411035756763200828211E-30_r8/
      DATA SQRTPI / 1.77245385090551602729816748334115_r8/
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DERFC
      IF (FIRST) THEN
         ETA = 0.1_r8*real(D1MACH(3),r8)
         NTERF = INITDS (ERFCS, 21, ETA)
         NTERFC = INITDS (ERFCCS, 59, ETA)
         NTERC2 = INITDS (ERC2CS, 49, ETA)
C
         XSML = -SQRT(-LOG(SQRTPI*D1MACH(3)))
         TXMAX = SQRT(-LOG(SQRTPI*D1MACH(1)))
         XMAX = TXMAX - 0.5_r8*LOG(TXMAX)/TXMAX - 0.01_r8
         SQEPS = SQRT(2.0_r8*D1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      IF (X.GT.XSML) GO TO 20
C
C ERFC(X) = 1.0 - ERF(X)  FOR  X .LT. XSML
C
      DERFC = 2.0_r8
      RETURN
C
 20   IF (X.GT.XMAX) GO TO 40
      Y = ABS(X)
      IF (Y.GT.1.0_r8) GO TO 30
C
C ERFC(X) = 1.0 - ERF(X)  FOR ABS(X) .LE. 1.0
C
      IF (Y.LT.SQEPS) DERFC = 1.0_r8- 2.0_r8*X/SQRTPI
      IF (Y.GE.SQEPS) DERFC = 1.0_r8- X*(1.0_r8+
     1   DCSEVL (2._r8*X*X-1._r8,
     1   ERFCS, NTERF))
      RETURN
C
C ERFC(X) = 1.0 - ERF(X)  FOR  1.0 .LT. ABS(X) .LE. XMAX
C
 30   Y = Y*Y
      IF (Y.LE.4._r8) DERFC = EXP(-Y)/ABS(X) * (0.5_r8+
     1   DCSEVL (
     1   (8._r8/Y-5._r8)/3._r8, ERC2CS, NTERC2) )
      IF (Y.GT.4._r8) DERFC = EXP(-Y)/ABS(X) * (0.5_r8+
     1   DCSEVL (
     1   8._r8/Y-1._r8, ERFCCS, NTERFC) )
      IF (X.LT.0._r8) DERFC = 2.0_r8- DERFC
      RETURN
C
 40   CALL XERMSG ('SLATEC', 'DERFC', 'X SO BIG ERFC UNDERFLOWS', 1, 1)
      DERFC = 0._r8
      RETURN
C
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
! 22Jun2000 fgtok
