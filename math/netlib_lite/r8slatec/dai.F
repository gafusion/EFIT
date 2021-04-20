*DECK DAI
      REAL*8 FUNCTION DAI (X)
C***BEGIN PROLOGUE  DAI
C***PURPOSE  Evaluate the Airy function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10D
C***TYPE      REAL*8 (AI-S, DAI-D)
C***KEYWORDS  AIRY FUNCTION, FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DAI(X) calculates the double precision Airy function for double
C precision argument X.
C
C Series for AIF        on the interval -1.00000E+00 to  1.00000E+00
C                                        with weighted error   8.37E-33
C                                         log weighted error  32.08
C                               significant figures required  30.87
C                                    decimal places required  32.63
C
C Series for AIG        on the interval -1.00000E+00 to  1.00000E+00
C                                        with weighted error   7.47E-34
C                                         log weighted error  33.13
C                               significant figures required  31.50
C                                    decimal places required  33.68
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, D9AIMP, DAIE, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920618  Removed space from variable names.  (RWC, WRB)
C***END PROLOGUE  DAI
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER naif,naig,initds
C============
      REAL*8 X, AIFCS(13), AIGCS(13), THETA, XM, XMAX, X3SML,
     1  Z, D1MACH, DCSEVL, DAIE, XMAXT
      LOGICAL FIRST
      SAVE AIFCS, AIGCS, NAIF, NAIG, X3SML, XMAX, FIRST
      DATA AIFCS(  1) / -.37971358496669997496197089469414E-1_r8/
      DATA AIFCS(  2) / +.59191888537263638574319728013777E-1_r8/
      DATA AIFCS(  3) / +.98629280577279975365603891044060E-3_r8/
      DATA AIFCS(  4) / +.68488438190765667554854830182412E-5_r8/
      DATA AIFCS(  5) / +.25942025962194713019489279081403E-7_r8/
      DATA AIFCS(  6) / +.61766127740813750329445749697236E-10_r8/
      DATA AIFCS(  7) / +.10092454172466117901429556224601E-12_r8/
      DATA AIFCS(  8) / +.12014792511179938141288033225333E-15_r8/
      DATA AIFCS(  9) / +.10882945588716991878525295466666E-18_r8/
      DATA AIFCS( 10) / +.77513772196684887039238400000000E-22_r8/
      DATA AIFCS( 11) / +.44548112037175638391466666666666E-25_r8/
      DATA AIFCS( 12) / +.21092845231692343466666666666666E-28_r8/
      DATA AIFCS( 13) / +.83701735910741333333333333333333E-32_r8/
      DATA AIGCS(  1) / +.18152365581161273011556209957864E-1_r8/
      DATA AIGCS(  2) / +.21572563166010755534030638819968E-1_r8/
      DATA AIGCS(  3) / +.25678356987483249659052428090133E-3_r8/
      DATA AIGCS(  4) / +.14265214119792403898829496921721E-5_r8/
      DATA AIGCS(  5) / +.45721149200180426070434097558191E-8_r8/
      DATA AIGCS(  6) / +.95251708435647098607392278840592E-11_r8/
      DATA AIGCS(  7) / +.13925634605771399051150420686190E-13_r8/
      DATA AIGCS(  8) / +.15070999142762379592306991138666E-16_r8/
      DATA AIGCS(  9) / +.12559148312567778822703205333333E-19_r8/
      DATA AIGCS( 10) / +.83063073770821340343829333333333E-23_r8/
      DATA AIGCS( 11) / +.44657538493718567445333333333333E-26_r8/
      DATA AIGCS( 12) / +.19900855034518869333333333333333E-29_r8/
      DATA AIGCS( 13) / +.74702885256533333333333333333333E-33_r8/
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DAI
      IF (FIRST) THEN
         NAIF = INITDS (AIFCS, 13, 0.1_r8*real(D1MACH(3),r8))
         NAIG = INITDS (AIGCS, 13, 0.1_r8*real(D1MACH(3),r8))
C
         X3SML = D1MACH(3)**0.3334_r8
         XMAXT = (-1.5_r8*LOG(D1MACH(1)))**0.6667_r8
         XMAX = XMAXT - XMAXT*LOG(XMAXT)/
     *      (4.0_r8*SQRT(XMAXT)+1.0_r8)
     *      - 0.01_r8
      ENDIF
      FIRST = .FALSE.
C
      IF (X.GE.(-1._r8)) GO TO 20
      CALL D9AIMP (X, XM, THETA)
      DAI = XM * COS(THETA)
      RETURN
C
 20   IF (X.GT.1.0_r8) GO TO 30
      Z = 0.0_r8
      IF (ABS(X).GT.X3SML) Z = X**3
      DAI = 0.375_r8+ (DCSEVL (Z, AIFCS, NAIF) - X*(0.25_r8+
     1  DCSEVL (Z, AIGCS, NAIG)) )
      RETURN
C
 30   IF (X.GT.XMAX) GO TO 40
      DAI = DAIE(X) * EXP(-2.0_r8*X*SQRT(X)/3.0_r8)
      RETURN
C
 40   DAI = 0.0_r8
      CALL XERMSG ('SLATEC', 'DAI', 'X SO BIG AI UNDERFLOWS', 1, 1)
      RETURN
C
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
! 22Jun2000 fgtok
