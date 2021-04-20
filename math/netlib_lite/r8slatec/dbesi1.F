*DECK DBESI1
      REAL*8 FUNCTION DBESI1 (X)
C***BEGIN PROLOGUE  DBESI1
C***PURPOSE  Compute the modified (hyperbolic) Bessel function of the
C            first kind of order one.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B1
C***TYPE      REAL*8 (BESI1-S, DBESI1-D)
C***KEYWORDS  FIRST KIND, FNLIB, HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBESI1(X) calculates the double precision modified (hyperbolic)
C Bessel function of the first kind of order one and double precision
C argument X.
C
C Series for BI1        on the interval  0.          to  9.00000E+00
C                                        with weighted error   1.44E-32
C                                         log weighted error  31.84
C                               significant figures required  31.45
C                                    decimal places required  32.46
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DBSI1E, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DBESI1
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER nti1,initds
C============
      REAL*8 X, BI1CS(17), XMAX, XMIN, XSML, Y, D1MACH,
     1  DCSEVL, DBSI1E
      LOGICAL FIRST
      SAVE BI1CS, NTI1, XMIN, XSML, XMAX, FIRST
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
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DBESI1
      IF (FIRST) THEN
         NTI1 = INITDS (BI1CS, 17, 0.1_r8*real(D1MACH(3),r8))
         XMIN = 2.0_r8*D1MACH(1)
         XSML = SQRT(4.5_r8*D1MACH(3))
         XMAX = LOG (D1MACH(2))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.3.0_r8) GO TO 20
C
      DBESI1 = 0._r8
      IF (Y.EQ.0._r8)  RETURN
C
      IF (Y .LE. XMIN) CALL XERMSG ('SLATEC', 'DBESI1',
     +   'ABS(X) SO SMALL I1 UNDERFLOWS', 1, 1)
      IF (Y.GT.XMIN) DBESI1 = 0.5_r8*X
      IF (Y.GT.XSML) DBESI1 = X*(0.875_r8+
     1   DCSEVL (Y*Y/4.5_r8-1._r8,
     1   BI1CS, NTI1))
      RETURN
C
 20   IF (Y .GT. XMAX) CALL XERMSG ('SLATEC', 'DBESI1',
     +   'ABS(X) SO BIG I1 OVERFLOWS', 2, 2)
C
      DBESI1 = EXP(Y) * DBSI1E(X)
C
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
! 22Jun2000 fgtok
