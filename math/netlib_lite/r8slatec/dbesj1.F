*DECK DBESJ1
      REAL*8 FUNCTION DBESJ1 (X)
C***BEGIN PROLOGUE  DBESJ1
C***PURPOSE  Compute the Bessel function of the first kind of order one.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10A1
C***TYPE      REAL*8 (BESJ1-S, DBESJ1-D)
C***KEYWORDS  BESSEL FUNCTION, FIRST KIND, FNLIB, ORDER ONE,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBESJ1(X) calculates the double precision Bessel function of the
C first kind of order one for double precision argument X.
C
C Series for BJ1        on the interval  0.          to  1.60000E+01
C                                        with weighted error   1.16E-33
C                                         log weighted error  32.93
C                               significant figures required  32.36
C                                    decimal places required  33.57
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, D9B1MP, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   780601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   910401  Corrected error in code which caused values to have the
C           wrong sign for arguments less than 4.0.  (WRB)
C***END PROLOGUE  DBESJ1
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER ntj1,initds
C============
      REAL*8 X, BJ1CS(19), AMPL, THETA, XSML, XMIN, Y,
     1  D1MACH, DCSEVL
      LOGICAL FIRST
      SAVE BJ1CS, NTJ1, XSML, XMIN, FIRST
      DATA BJ1CS(  1) / -.117261415133327865606240574524003E+0_r8/
      DATA BJ1CS(  2) / -.253615218307906395623030884554698E+0_r8/
      DATA BJ1CS(  3) / +.501270809844695685053656363203743E-1_r8/
      DATA BJ1CS(  4) / -.463151480962508191842619728789772E-2_r8/
      DATA BJ1CS(  5) / +.247996229415914024539124064592364E-3_r8/
      DATA BJ1CS(  6) / -.867894868627882584521246435176416E-5_r8/
      DATA BJ1CS(  7) / +.214293917143793691502766250991292E-6_r8/
      DATA BJ1CS(  8) / -.393609307918317979229322764073061E-8_r8/
      DATA BJ1CS(  9) / +.559118231794688004018248059864032E-10_r8/
      DATA BJ1CS( 10) / -.632761640466139302477695274014880E-12_r8/
      DATA BJ1CS( 11) / +.584099161085724700326945563268266E-14_r8/
      DATA BJ1CS( 12) / -.448253381870125819039135059199999E-16_r8/
      DATA BJ1CS( 13) / +.290538449262502466306018688000000E-18_r8/
      DATA BJ1CS( 14) / -.161173219784144165412118186666666E-20_r8/
      DATA BJ1CS( 15) / +.773947881939274637298346666666666E-23_r8/
      DATA BJ1CS( 16) / -.324869378211199841143466666666666E-25_r8/
      DATA BJ1CS( 17) / +.120223767722741022720000000000000E-27_r8/
      DATA BJ1CS( 18) / -.395201221265134933333333333333333E-30_r8/
      DATA BJ1CS( 19) / +.116167808226645333333333333333333E-32_r8/
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DBESJ1
      IF (FIRST) THEN
         NTJ1 = INITDS (BJ1CS, 19, 0.1_r8*real(D1MACH(3),r8))
C
         XSML = SQRT(8.0_r8*D1MACH(3))
         XMIN = 2.0_r8*D1MACH(1)
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.4.0_r8) GO TO 20
C
      DBESJ1 = 0.0_r8
      IF (Y.EQ.0.0_r8) RETURN
      IF (Y .LE. XMIN) CALL XERMSG ('SLATEC', 'DBESJ1',
     +   'ABS(X) SO SMALL J1 UNDERFLOWS', 1, 1)
      IF (Y.GT.XMIN) DBESJ1 = 0.5_r8*X
      IF (Y.GT.XSML) DBESJ1 =
     1   X*(.25_r8+ DCSEVL (.125_r8*Y*Y-1._r8,
     1   BJ1CS, NTJ1) )
      RETURN
C
 20   CALL D9B1MP (Y, AMPL, THETA)
      DBESJ1 = SIGN (AMPL, X) * COS(THETA)
C
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
! 22Jun2000 fgtok
