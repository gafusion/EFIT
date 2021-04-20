*DECK DBESK1
      REAL*8 FUNCTION DBESK1 (X)
C***BEGIN PROLOGUE  DBESK1
C***PURPOSE  Compute the modified (hyperbolic) Bessel function of the
C            third kind of order one.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B1
C***TYPE      REAL*8 (BESK1-S, DBESK1-D)
C***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS,
C             THIRD KIND
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBESK1(X) calculates the double precision modified (hyperbolic)
C Bessel function of the third kind of order one for double precision
C argument X.  The argument must be large enough that the result does
C not overflow and small enough that the result does not underflow.
C
C Series for BK1        on the interval  0.          to  4.00000E+00
C                                        with weighted error   9.16E-32
C                                         log weighted error  31.04
C                               significant figures required  30.61
C                                    decimal places required  31.64
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DBESI1, DBSK1E, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DBESK1
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER ntk1,initds
C============
      REAL*8 X, BK1CS(16), XMAX, XMAXT, XMIN, XSML, Y,
     1  D1MACH, DCSEVL, DBESI1, DBSK1E
      LOGICAL FIRST
      SAVE BK1CS, NTK1, XMIN, XSML, XMAX, FIRST
      DATA BK1CS(  1) / +.25300227338947770532531120868533E-1_r8/
      DATA BK1CS(  2) / -.35315596077654487566723831691801E+0_r8/
      DATA BK1CS(  3) / -.12261118082265714823479067930042E+0_r8/
      DATA BK1CS(  4) / -.69757238596398643501812920296083E-2_r8/
      DATA BK1CS(  5) / -.17302889575130520630176507368979E-3_r8/
      DATA BK1CS(  6) / -.24334061415659682349600735030164E-5_r8/
      DATA BK1CS(  7) / -.22133876307347258558315252545126E-7_r8/
      DATA BK1CS(  8) / -.14114883926335277610958330212608E-9_r8/
      DATA BK1CS(  9) / -.66669016941993290060853751264373E-12_r8/
      DATA BK1CS( 10) / -.24274498505193659339263196864853E-14_r8/
      DATA BK1CS( 11) / -.70238634793862875971783797120000E-17_r8/
      DATA BK1CS( 12) / -.16543275155100994675491029333333E-19_r8/
      DATA BK1CS( 13) / -.32338347459944491991893333333333E-22_r8/
      DATA BK1CS( 14) / -.53312750529265274999466666666666E-25_r8/
      DATA BK1CS( 15) / -.75130407162157226666666666666666E-28_r8/
      DATA BK1CS( 16) / -.91550857176541866666666666666666E-31_r8/
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DBESK1
      IF (FIRST) THEN
         NTK1 = INITDS (BK1CS, 16, 0.1_r8*real(D1MACH(3),r8))
         XMIN = EXP(MAX(LOG(D1MACH(1)), -LOG(D1MACH(2))) +
     +      0.01_r8)
         XSML = SQRT(4.0_r8*D1MACH(3))
         XMAXT = -LOG(D1MACH(1))
         XMAX = XMAXT - 0.5_r8*XMAXT*LOG(XMAXT)/
     +      (XMAXT+0.5_r8)
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LE. 0._r8) CALL XERMSG ('SLATEC', 'DBESK1',
     +   'X IS ZERO OR NEGATIVE', 2, 2)
      IF (X.GT.2.0_r8) GO TO 20
C
      IF (X .LT. XMIN) CALL XERMSG ('SLATEC', 'DBESK1',
     +   'X SO SMALL K1 OVERFLOWS', 3, 2)
      Y = 0._r8
      IF (X.GT.XSML) Y = X*X
      DBESK1 = LOG(0.5_r8*X)*DBESI1(X) +
     1   (0.75_r8+ DCSEVL (.5_r8*Y-1._r8,
     1   BK1CS, NTK1))/X
      RETURN
C
 20   DBESK1 = 0._r8
      IF (X .GT. XMAX) CALL XERMSG ('SLATEC', 'DBESK1',
     +   'X SO BIG K1 UNDERFLOWS', 1, 1)
      IF (X.GT.XMAX) RETURN
C
      DBESK1 = EXP(-X) * DBSK1E(X)
C
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
! 22Jun2000 fgtok
