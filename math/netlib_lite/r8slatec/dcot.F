*DECK DCOT
      REAL*8 FUNCTION DCOT (X)
C***BEGIN PROLOGUE  DCOT
C***PURPOSE  Compute the cotangent.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4A
C***TYPE      REAL*8 (COT-S, DCOT-D, CCOT-C)
C***KEYWORDS  COTANGENT, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DCOT(X) calculates the double precision trigonometric cotangent
C for double precision argument X.  X is in units of radians.
C
C Series for COT        on the interval  0.          to  6.25000E-02
C                                        with weighted error   5.52E-34
C                                         log weighted error  33.26
C                               significant figures required  32.34
C                                    decimal places required  33.85
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920618  Removed space from variable names.  (RWC, WRB)
C***END PROLOGUE  DCOT
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER nterms,initds,ifn
C============
      REAL*8 X, COTCS(15), AINTY, AINTY2, PI2REC, SQEPS,
     1  XMAX, XMIN, XSML, Y, YREM, PRODBG, DCSEVL, D1MACH
      LOGICAL FIRST
      SAVE COTCS, PI2REC, NTERMS, XMAX, XSML, XMIN, SQEPS, FIRST
      DATA COTCS(  1) / +.240259160982956302509553617744970E+0_r8/
      DATA COTCS(  2) / -.165330316015002278454746025255758E-1_r8/
      DATA COTCS(  3) / -.429983919317240189356476228239895E-4_r8/
      DATA COTCS(  4) / -.159283223327541046023490851122445E-6_r8/
      DATA COTCS(  5) / -.619109313512934872588620579343187E-9_r8/
      DATA COTCS(  6) / -.243019741507264604331702590579575E-11_r8/
      DATA COTCS(  7) / -.956093675880008098427062083100000E-14_r8/
      DATA COTCS(  8) / -.376353798194580580416291539706666E-16_r8/
      DATA COTCS(  9) / -.148166574646746578852176794666666E-18_r8/
      DATA COTCS( 10) / -.583335658903666579477984000000000E-21_r8/
      DATA COTCS( 11) / -.229662646964645773928533333333333E-23_r8/
      DATA COTCS( 12) / -.904197057307483326719999999999999E-26_r8/
      DATA COTCS( 13) / -.355988551920600064000000000000000E-28_r8/
      DATA COTCS( 14) / -.140155139824298666666666666666666E-30_r8/
      DATA COTCS( 15) / -.551800436872533333333333333333333E-33_r8/
      DATA PI2REC / .011619772367581343075535053490057_r8/
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DCOT
      IF (FIRST) THEN
         NTERMS = INITDS (COTCS, 15, 0.1_r8*real(D1MACH(3),r8) )
         XMAX = 1.0_r8/D1MACH(4)
         XSML = SQRT(3.0_r8*D1MACH(3))
         XMIN = EXP (MAX(LOG(D1MACH(1)),
     >      -LOG(D1MACH(2))) + 0.01_r8)
         SQEPS = SQRT(D1MACH(4))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y .LT. XMIN) CALL XERMSG ('SLATEC', 'DCOT',
     +   'ABS(X) IS ZERO OR SO SMALL DCOT OVERFLOWS', 2, 2)
      IF (Y .GT. XMAX) CALL XERMSG ('SLATEC', 'DCOT',
     +   'NO PRECISION BECAUSE ABS(X) IS TOO BIG', 3, 2)
C
C CAREFULLY COMPUTE Y * (2/PI) = (AINT(Y) + REM(Y)) * (.625 + PI2REC)
C = AINT(.625*Y) + REM(.625*Y) + Y*PI2REC  =  AINT(.625*Y) + Z
C = AINT(.625*Y) + AINT(Z) + REM(Z)
C
      AINTY = AINT (Y)
      YREM = Y - AINTY
      PRODBG = 0.625_r8*AINTY
      AINTY = AINT (PRODBG)
      Y = (PRODBG-AINTY) + 0.625_r8*YREM + PI2REC*Y
      AINTY2 = AINT (Y)
      AINTY = AINTY + AINTY2
      Y = Y - AINTY2
C
      IFN = MOD (AINTY, 2.0_r8)
      IF (IFN.EQ.1) Y = 1.0_r8- Y
C
      IF (ABS(X) .GT. 0.5_r8.AND. Y .LT. ABS(X)*SQEPS) CALL XERMSG
     +   ('SLATEC', 'DCOT',
     +   'ANSWER LT HALF PRECISION, ABS(X) TOO BIG OR X NEAR N*PI ' //
     +   '(N.NE.0)', 1, 1)
C
      IF (Y.GT.0.25_r8) GO TO 20
      DCOT = 1.0_r8/X
      IF (Y.GT.XSML) DCOT = (0.5_r8+
     1   DCSEVL (32.0_r8*Y*Y-1._r8, COTCS,
     1   NTERMS)) / Y
      GO TO 40
C
 20   IF (Y.GT.0.5_r8) GO TO 30
      DCOT = (0.5_r8+ DCSEVL (8._r8*Y*Y-1._r8,
     >   COTCS, NTERMS))/(0.5_r8*Y)
      DCOT = (DCOT*DCOT-1._r8)*0.5_r8/DCOT
      GO TO 40
C
 30   DCOT = (0.5_r8+ DCSEVL (2._r8*Y*Y-1._r8,
     >   COTCS, NTERMS))/(.25_r8*Y)
      DCOT = (DCOT*DCOT-1._r8)*0.5_r8/DCOT
      DCOT = (DCOT*DCOT-1._r8)*0.5_r8/DCOT
C
 40   IF (X.NE.0._r8) DCOT = SIGN (DCOT, X)
      IF (IFN.EQ.1) DCOT = -DCOT
C
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
! 22Jun2000 fgtok
