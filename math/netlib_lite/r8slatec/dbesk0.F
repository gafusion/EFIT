*DECK DBESK0
      REAL*8 FUNCTION DBESK0 (X)
C***BEGIN PROLOGUE  DBESK0
C***PURPOSE  Compute the modified (hyperbolic) Bessel function of the
C            third kind of order zero.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B1
C***TYPE      REAL*8 (BESK0-S, DBESK0-D)
C***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS,
C             THIRD KIND
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBESK0(X) calculates the double precision modified (hyperbolic)
C Bessel function of the third kind of order zero for double
C precision argument X.  The argument must be greater than zero
C but not so large that the result underflows.
C
C Series for BK0        on the interval  0.          to  4.00000E+00
C                                        with weighted error   3.08E-33
C                                         log weighted error  32.51
C                               significant figures required  32.05
C                                    decimal places required  33.11
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DBESI0, DBSK0E, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DBESK0
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER ntk0,initds
C============
      REAL*8 X, BK0CS(16), XMAX, XMAXT, XSML, Y,
     1  D1MACH, DCSEVL, DBESI0, DBSK0E
      LOGICAL FIRST
      SAVE BK0CS, NTK0, XSML, XMAX, FIRST
      DATA BK0CS(  1) / -.353273932339027687201140060063153E-1_r8/
      DATA BK0CS(  2) / +.344289899924628486886344927529213E+0_r8/
      DATA BK0CS(  3) / +.359799365153615016265721303687231E-1_r8/
      DATA BK0CS(  4) / +.126461541144692592338479508673447E-2_r8/
      DATA BK0CS(  5) / +.228621210311945178608269830297585E-4_r8/
      DATA BK0CS(  6) / +.253479107902614945730790013428354E-6_r8/
      DATA BK0CS(  7) / +.190451637722020885897214059381366E-8_r8/
      DATA BK0CS(  8) / +.103496952576336245851008317853089E-10_r8/
      DATA BK0CS(  9) / +.425981614279108257652445327170133E-13_r8/
      DATA BK0CS( 10) / +.137446543588075089694238325440000E-15_r8/
      DATA BK0CS( 11) / +.357089652850837359099688597333333E-18_r8/
      DATA BK0CS( 12) / +.763164366011643737667498666666666E-21_r8/
      DATA BK0CS( 13) / +.136542498844078185908053333333333E-23_r8/
      DATA BK0CS( 14) / +.207527526690666808319999999999999E-26_r8/
      DATA BK0CS( 15) / +.271281421807298560000000000000000E-29_r8/
      DATA BK0CS( 16) / +.308259388791466666666666666666666E-32_r8/
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DBESK0
      IF (FIRST) THEN
         NTK0 = INITDS (BK0CS, 16, 0.1_r8*real(D1MACH(3),r8))
         XSML = SQRT(4.0_r8*D1MACH(3))
         XMAXT = -LOG(D1MACH(1))
         XMAX = XMAXT -
     +      0.5_r8*XMAXT*LOG(XMAXT)/(XMAXT+0.5_r8)
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LE. 0._r8) CALL XERMSG ('SLATEC', 'DBESK0',
     +   'X IS ZERO OR NEGATIVE', 2, 2)
      IF (X.GT.2.0_r8) GO TO 20
C
      Y = 0._r8
      IF (X.GT.XSML) Y = X*X
      DBESK0 = -LOG(0.5_r8*X)*DBESI0(X) -
     1   0.25_r8+ DCSEVL (.5_r8*Y-1._r8,
     1   BK0CS, NTK0)
      RETURN
C
 20   DBESK0 = 0._r8
      IF (X .GT. XMAX) CALL XERMSG ('SLATEC', 'DBESK0',
     +   'X SO BIG K0 UNDERFLOWS', 1, 1)
      IF (X.GT.XMAX) RETURN
C
      DBESK0 = EXP(-X) * DBSK0E(X)
C
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
! 22Jun2000 fgtok
