*DECK DBSI1E
      REAL*8 FUNCTION DBSI1E (X)
C***BEGIN PROLOGUE  DBSI1E
C***PURPOSE  Compute the exponentially scaled modified (hyperbolic)
C            Bessel function of the first kind of order one.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B1
C***TYPE      REAL*8 (BESI1E-S, DBSI1E-D)
C***KEYWORDS  EXPONENTIALLY SCALED, FIRST KIND, FNLIB,
C             HYPERBOLIC BESSEL FUNCTION, MODIFIED BESSEL FUNCTION,
C             ORDER ONE, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBSI1E(X) calculates the double precision exponentially scaled
C modified (hyperbolic) Bessel function of the first kind of order
C one for double precision argument X.  The result is I1(X)
C multiplied by EXP(-ABS(X)).
C
C Series for BI1        on the interval  0.          to  9.00000E+00
C                                        with weighted error   1.44E-32
C                                         log weighted error  31.84
C                               significant figures required  31.45
C                                    decimal places required  32.46
C
C Series for AI1        on the interval  1.25000E-01 to  3.33333E-01
C                                        with weighted error   2.81E-32
C                                         log weighted error  31.55
C                               significant figures required  29.93
C                                    decimal places required  32.38
C
C Series for AI12       on the interval  0.          to  1.25000E-01
C                                        with weighted error   1.83E-32
C                                         log weighted error  31.74
C                               significant figures required  29.97
C                                    decimal places required  32.66
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DBSI1E
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER nti1,ntai1,ntai12,initds
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 eta
C============
      REAL*8 X, BI1CS(17), AI1CS(46), AI12CS(69), XMIN,
     1  XSML, Y, D1MACH, DCSEVL
      LOGICAL FIRST
      SAVE BI1CS, AI1CS, AI12CS, NTI1, NTAI1, NTAI12, XMIN, XSML,
     1  FIRST
      DATA BI1CS(  1) / -.19717132610998597316138503218149E-2_r8/
      DATA BI1CS(  2) / +.40734887667546480608155393652014E+0_r8/
      DATA BI1CS(  3) / +.34838994299959455866245037783787E-1_r8/
      DATA BI1CS(  4) / +.15453945563001236038598401058489E-2_r8/
      DATA BI1CS(  5) / +.41888521098377784129458832004120E-4_r8/
      DATA BI1CS(  6) / +.76490267648362114741959703966069E-6_r8/
      DATA BI1CS(  7) / +.10042493924741178689179808037238E-7_r8/
      DATA BI1CS(  8) / +.99322077919238106481371298054863E-10_r8/
      DATA BI1CS(  9) / +.76638017918447637275200171681349E-12_r8/
      DATA BI1CS( 10) / +.47414189238167394980388091948160E-14_r8/
      DATA BI1CS( 11) / +.24041144040745181799863172032000E-16_r8/
      DATA BI1CS( 12) / +.10171505007093713649121100799999E-18_r8/
      DATA BI1CS( 13) / +.36450935657866949458491733333333E-21_r8/
      DATA BI1CS( 14) / +.11205749502562039344810666666666E-23_r8/
      DATA BI1CS( 15) / +.29875441934468088832000000000000E-26_r8/
      DATA BI1CS( 16) / +.69732310939194709333333333333333E-29_r8/
      DATA BI1CS( 17) / +.14367948220620800000000000000000E-31_r8/
      DATA AI1CS(  1) / -.2846744181881478674100372468307E-1_r8/
      DATA AI1CS(  2) / -.1922953231443220651044448774979E-1_r8/
      DATA AI1CS(  3) / -.6115185857943788982256249917785E-3_r8/
      DATA AI1CS(  4) / -.2069971253350227708882823777979E-4_r8/
      DATA AI1CS(  5) / +.8585619145810725565536944673138E-5_r8/
      DATA AI1CS(  6) / +.1049498246711590862517453997860E-5_r8/
      DATA AI1CS(  7) / -.2918338918447902202093432326697E-6_r8/
      DATA AI1CS(  8) / -.1559378146631739000160680969077E-7_r8/
      DATA AI1CS(  9) / +.1318012367144944705525302873909E-7_r8/
      DATA AI1CS( 10) / -.1448423418183078317639134467815E-8_r8/
      DATA AI1CS( 11) / -.2908512243993142094825040993010E-9_r8/
      DATA AI1CS( 12) / +.1266388917875382387311159690403E-9_r8/
      DATA AI1CS( 13) / -.1664947772919220670624178398580E-10_r8/
      DATA AI1CS( 14) / -.1666653644609432976095937154999E-11_r8/
      DATA AI1CS( 15) / +.1242602414290768265232168472017E-11_r8/
      DATA AI1CS( 16) / -.2731549379672432397251461428633E-12_r8/
      DATA AI1CS( 17) / +.2023947881645803780700262688981E-13_r8/
      DATA AI1CS( 18) / +.7307950018116883636198698126123E-14_r8/
      DATA AI1CS( 19) / -.3332905634404674943813778617133E-14_r8/
      DATA AI1CS( 20) / +.7175346558512953743542254665670E-15_r8/
      DATA AI1CS( 21) / -.6982530324796256355850629223656E-16_r8/
      DATA AI1CS( 22) / -.1299944201562760760060446080587E-16_r8/
      DATA AI1CS( 23) / +.8120942864242798892054678342860E-17_r8/
      DATA AI1CS( 24) / -.2194016207410736898156266643783E-17_r8/
      DATA AI1CS( 25) / +.3630516170029654848279860932334E-18_r8/
      DATA AI1CS( 26) / -.1695139772439104166306866790399E-19_r8/
      DATA AI1CS( 27) / -.1288184829897907807116882538222E-19_r8/
      DATA AI1CS( 28) / +.5694428604967052780109991073109E-20_r8/
      DATA AI1CS( 29) / -.1459597009090480056545509900287E-20_r8/
      DATA AI1CS( 30) / +.2514546010675717314084691334485E-21_r8/
      DATA AI1CS( 31) / -.1844758883139124818160400029013E-22_r8/
      DATA AI1CS( 32) / -.6339760596227948641928609791999E-23_r8/
      DATA AI1CS( 33) / +.3461441102031011111108146626560E-23_r8/
      DATA AI1CS( 34) / -.1017062335371393547596541023573E-23_r8/
      DATA AI1CS( 35) / +.2149877147090431445962500778666E-24_r8/
      DATA AI1CS( 36) / -.3045252425238676401746206173866E-25_r8/
      DATA AI1CS( 37) / +.5238082144721285982177634986666E-27_r8/
      DATA AI1CS( 38) / +.1443583107089382446416789503999E-26_r8/
      DATA AI1CS( 39) / -.6121302074890042733200670719999E-27_r8/
      DATA AI1CS( 40) / +.1700011117467818418349189802666E-27_r8/
      DATA AI1CS( 41) / -.3596589107984244158535215786666E-28_r8/
      DATA AI1CS( 42) / +.5448178578948418576650513066666E-29_r8/
      DATA AI1CS( 43) / -.2731831789689084989162564266666E-30_r8/
      DATA AI1CS( 44) / -.1858905021708600715771903999999E-30_r8/
      DATA AI1CS( 45) / +.9212682974513933441127765333333E-31_r8/
      DATA AI1CS( 46) / -.2813835155653561106370833066666E-31_r8/
      DATA AI12CS(  1) / +.2857623501828012047449845948469E-1_r8/
      DATA AI12CS(  2) / -.9761097491361468407765164457302E-2_r8/
      DATA AI12CS(  3) / -.1105889387626237162912569212775E-3_r8/
      DATA AI12CS(  4) / -.3882564808877690393456544776274E-5_r8/
      DATA AI12CS(  5) / -.2512236237870208925294520022121E-6_r8/
      DATA AI12CS(  6) / -.2631468846889519506837052365232E-7_r8/
      DATA AI12CS(  7) / -.3835380385964237022045006787968E-8_r8/
      DATA AI12CS(  8) / -.5589743462196583806868112522229E-9_r8/
      DATA AI12CS(  9) / -.1897495812350541234498925033238E-10_r8/
      DATA AI12CS( 10) / +.3252603583015488238555080679949E-10_r8/
      DATA AI12CS( 11) / +.1412580743661378133163366332846E-10_r8/
      DATA AI12CS( 12) / +.2035628544147089507224526136840E-11_r8/
      DATA AI12CS( 13) / -.7198551776245908512092589890446E-12_r8/
      DATA AI12CS( 14) / -.4083551111092197318228499639691E-12_r8/
      DATA AI12CS( 15) / -.2101541842772664313019845727462E-13_r8/
      DATA AI12CS( 16) / +.4272440016711951354297788336997E-13_r8/
      DATA AI12CS( 17) / +.1042027698412880276417414499948E-13_r8/
      DATA AI12CS( 18) / -.3814403072437007804767072535396E-14_r8/
      DATA AI12CS( 19) / -.1880354775510782448512734533963E-14_r8/
      DATA AI12CS( 20) / +.3308202310920928282731903352405E-15_r8/
      DATA AI12CS( 21) / +.2962628997645950139068546542052E-15_r8/
      DATA AI12CS( 22) / -.3209525921993423958778373532887E-16_r8/
      DATA AI12CS( 23) / -.4650305368489358325571282818979E-16_r8/
      DATA AI12CS( 24) / +.4414348323071707949946113759641E-17_r8/
      DATA AI12CS( 25) / +.7517296310842104805425458080295E-17_r8/
      DATA AI12CS( 26) / -.9314178867326883375684847845157E-18_r8/
      DATA AI12CS( 27) / -.1242193275194890956116784488697E-17_r8/
      DATA AI12CS( 28) / +.2414276719454848469005153902176E-18_r8/
      DATA AI12CS( 29) / +.2026944384053285178971922860692E-18_r8/
      DATA AI12CS( 30) / -.6394267188269097787043919886811E-19_r8/
      DATA AI12CS( 31) / -.3049812452373095896084884503571E-19_r8/
      DATA AI12CS( 32) / +.1612841851651480225134622307691E-19_r8/
      DATA AI12CS( 33) / +.3560913964309925054510270904620E-20_r8/
      DATA AI12CS( 34) / -.3752017947936439079666828003246E-20_r8/
      DATA AI12CS( 35) / -.5787037427074799345951982310741E-22_r8/
      DATA AI12CS( 36) / +.7759997511648161961982369632092E-21_r8/
      DATA AI12CS( 37) / -.1452790897202233394064459874085E-21_r8/
      DATA AI12CS( 38) / -.1318225286739036702121922753374E-21_r8/
      DATA AI12CS( 39) / +.6116654862903070701879991331717E-22_r8/
      DATA AI12CS( 40) / +.1376279762427126427730243383634E-22_r8/
      DATA AI12CS( 41) / -.1690837689959347884919839382306E-22_r8/
      DATA AI12CS( 42) / +.1430596088595433153987201085385E-23_r8/
      DATA AI12CS( 43) / +.3409557828090594020405367729902E-23_r8/
      DATA AI12CS( 44) / -.1309457666270760227845738726424E-23_r8/
      DATA AI12CS( 45) / -.3940706411240257436093521417557E-24_r8/
      DATA AI12CS( 46) / +.4277137426980876580806166797352E-24_r8/
      DATA AI12CS( 47) / -.4424634830982606881900283123029E-25_r8/
      DATA AI12CS( 48) / -.8734113196230714972115309788747E-25_r8/
      DATA AI12CS( 49) / +.4045401335683533392143404142428E-25_r8/
      DATA AI12CS( 50) / +.7067100658094689465651607717806E-26_r8/
      DATA AI12CS( 51) / -.1249463344565105223002864518605E-25_r8/
      DATA AI12CS( 52) / +.2867392244403437032979483391426E-26_r8/
      DATA AI12CS( 53) / +.2044292892504292670281779574210E-26_r8/
      DATA AI12CS( 54) / -.1518636633820462568371346802911E-26_r8/
      DATA AI12CS( 55) / +.8110181098187575886132279107037E-28_r8/
      DATA AI12CS( 56) / +.3580379354773586091127173703270E-27_r8/
      DATA AI12CS( 57) / -.1692929018927902509593057175448E-27_r8/
      DATA AI12CS( 58) / -.2222902499702427639067758527774E-28_r8/
      DATA AI12CS( 59) / +.5424535127145969655048600401128E-28_r8/
      DATA AI12CS( 60) / -.1787068401578018688764912993304E-28_r8/
      DATA AI12CS( 61) / -.6565479068722814938823929437880E-29_r8/
      DATA AI12CS( 62) / +.7807013165061145280922067706839E-29_r8/
      DATA AI12CS( 63) / -.1816595260668979717379333152221E-29_r8/
      DATA AI12CS( 64) / -.1287704952660084820376875598959E-29_r8/
      DATA AI12CS( 65) / +.1114548172988164547413709273694E-29_r8/
      DATA AI12CS( 66) / -.1808343145039336939159368876687E-30_r8/
      DATA AI12CS( 67) / -.2231677718203771952232448228939E-30_r8/
      DATA AI12CS( 68) / +.1619029596080341510617909803614E-30_r8/
      DATA AI12CS( 69) / -.1834079908804941413901308439210E-31_r8/
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DBSI1E
      IF (FIRST) THEN
         ETA = 0.1_r8*real(D1MACH(3),r8)
         NTI1 = INITDS (BI1CS, 17, ETA)
         NTAI1 = INITDS (AI1CS, 46, ETA)
         NTAI12 = INITDS (AI12CS, 69, ETA)
C
         XMIN = 2.0_r8*D1MACH(1)
         XSML = SQRT(4.5_r8*D1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.3.0_r8) GO TO 20
C
      DBSI1E = 0.0_r8
      IF (Y.EQ.0._r8)  RETURN
C
      IF (Y .LE. XMIN) CALL XERMSG ('SLATEC', 'DBSI1E',
     +   'ABS(X) SO SMALL I1 UNDERFLOWS', 1, 1)
      IF (Y.GT.XMIN) DBSI1E = 0.5_r8*X
      IF (Y.GT.XSML) DBSI1E = X*(0.875_r8+
     1   DCSEVL (Y*Y/4.5_r8-1._r8,
     1   BI1CS, NTI1) )
      DBSI1E = EXP(-Y) * DBSI1E
      RETURN
C
 20   IF (Y.LE.8._r8) DBSI1E = (0.375_r8+
     1   DCSEVL ((48._r8/Y-11._r8)/5._r8,
     1   AI1CS, NTAI1))/SQRT(Y)
      IF (Y.GT.8._r8) DBSI1E = (0.375_r8+
     1   DCSEVL (16._r8/Y-1._r8, AI12CS,
     1   NTAI12))/SQRT(Y)
      DBSI1E = SIGN (DBSI1E, X)
C
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
! 22Jun2000 fgtok
