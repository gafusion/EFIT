*DECK DBESY0
      REAL*8 FUNCTION DBESY0 (X)
C***BEGIN PROLOGUE  DBESY0
C***PURPOSE  Compute the Bessel function of the second kind of order
C            zero.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10A1
C***TYPE      REAL*8 (BESY0-S, DBESY0-D)
C***KEYWORDS  BESSEL FUNCTION, FNLIB, ORDER ZERO, SECOND KIND,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBESY0(X) calculates the double precision Bessel function of the
C second kind of order zero for double precision argument X.
C
C Series for BY0        on the interval  0.          to  1.60000E+01
C                                        with weighted error   8.14E-32
C                                         log weighted error  31.09
C                               significant figures required  30.31
C                                    decimal places required  31.73
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, D9B0MP, DBESJ0, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DBESY0
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER nty0,initds
C============
      REAL*8 X, BY0CS(19), AMPL, THETA, TWODPI, XSML,
     1  Y, D1MACH, DCSEVL, DBESJ0
      LOGICAL FIRST
      SAVE BY0CS, TWODPI, NTY0, XSML, FIRST
      DATA BY0CS(  1) / -.1127783939286557321793980546028E-1_r8/
      DATA BY0CS(  2) / -.1283452375604203460480884531838E+0_r8/
      DATA BY0CS(  3) / -.1043788479979424936581762276618E+0_r8/
      DATA BY0CS(  4) / +.2366274918396969540924159264613E-1_r8/
      DATA BY0CS(  5) / -.2090391647700486239196223950342E-2_r8/
      DATA BY0CS(  6) / +.1039754539390572520999246576381E-3_r8/
      DATA BY0CS(  7) / -.3369747162423972096718775345037E-5_r8/
      DATA BY0CS(  8) / +.7729384267670667158521367216371E-7_r8/
      DATA BY0CS(  9) / -.1324976772664259591443476068964E-8_r8/
      DATA BY0CS( 10) / +.1764823261540452792100389363158E-10_r8/
      DATA BY0CS( 11) / -.1881055071580196200602823012069E-12_r8/
      DATA BY0CS( 12) / +.1641865485366149502792237185749E-14_r8/
      DATA BY0CS( 13) / -.1195659438604606085745991006720E-16_r8/
      DATA BY0CS( 14) / +.7377296297440185842494112426666E-19_r8/
      DATA BY0CS( 15) / -.3906843476710437330740906666666E-21_r8/
      DATA BY0CS( 16) / +.1795503664436157949829120000000E-23_r8/
      DATA BY0CS( 17) / -.7229627125448010478933333333333E-26_r8/
      DATA BY0CS( 18) / +.2571727931635168597333333333333E-28_r8/
      DATA BY0CS( 19) / -.8141268814163694933333333333333E-31_r8/
      DATA TWODPI / 0.636619772367581343075535053490057_r8/
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DBESY0
      IF (FIRST) THEN
         NTY0 = INITDS (BY0CS, 19, 0.1_r8*real(D1MACH(3),r8))
         XSML = SQRT(4.0_r8*D1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LE. 0._r8) CALL XERMSG ('SLATEC', 'DBESY0',
     +   'X IS ZERO OR NEGATIVE', 1, 2)
      IF (X.GT.4.0_r8) GO TO 20
C
      Y = 0._r8
      IF (X.GT.XSML) Y = X*X
      DBESY0 = TWODPI*LOG(0.5_r8*X)*DBESJ0(X) +
     1   .375_r8+ DCSEVL (
     1   .125_r8*Y-1._r8, BY0CS, NTY0)
      RETURN
C
 20   CALL D9B0MP (X, AMPL, THETA)
      DBESY0 = AMPL * SIN(THETA)
      RETURN
C
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
! 22Jun2000 fgtok
