*DECK D9LGMC
      REAL*8 FUNCTION D9LGMC (X)
C***BEGIN PROLOGUE  D9LGMC
C***SUBSIDIARY
C***PURPOSE  Compute the log Gamma correction factor so that
C            LOG(DGAMMA(X)) = LOG(SQRT(2*PI)) + (X-5.)*LOG(X) - X
C            + D9LGMC(X).
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7E
C***TYPE      REAL*8 (R9LGMC-S, D9LGMC-D, C9LGMC-C)
C***KEYWORDS  COMPLETE GAMMA FUNCTION, CORRECTION TERM, FNLIB,
C             LOG GAMMA, LOGARITHM, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Compute the log gamma correction factor for X .GE. 10. so that
C LOG (DGAMMA(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X + D9lGMC(X)
C
C Series for ALGM       on the interval  0.          to  1.00000E-02
C                                        with weighted error   1.28E-31
C                                         log weighted error  30.89
C                               significant figures required  29.81
C                                    decimal places required  31.48
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900720  Routine changed from user-callable to subsidiary.  (WRB)
C***END PROLOGUE  D9LGMC
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER nalgm,initds
C============
      REAL*8 X, ALGMCS(15), XBIG, XMAX, DCSEVL, D1MACH
      LOGICAL FIRST
      SAVE ALGMCS, NALGM, XBIG, XMAX, FIRST
      DATA ALGMCS(  1) / +.1666389480451863247205729650822E+0_r8/
      DATA ALGMCS(  2) / -.1384948176067563840732986059135E-4_r8/
      DATA ALGMCS(  3) / +.9810825646924729426157171547487E-8_r8/
      DATA ALGMCS(  4) / -.1809129475572494194263306266719E-10_r8/
      DATA ALGMCS(  5) / +.6221098041892605227126015543416E-13_r8/
      DATA ALGMCS(  6) / -.3399615005417721944303330599666E-15_r8/
      DATA ALGMCS(  7) / +.2683181998482698748957538846666E-17_r8/
      DATA ALGMCS(  8) / -.2868042435334643284144622399999E-19_r8/
      DATA ALGMCS(  9) / +.3962837061046434803679306666666E-21_r8/
      DATA ALGMCS( 10) / -.6831888753985766870111999999999E-23_r8/
      DATA ALGMCS( 11) / +.1429227355942498147573333333333E-24_r8/
      DATA ALGMCS( 12) / -.3547598158101070547199999999999E-26_r8/
      DATA ALGMCS( 13) / +.1025680058010470912000000000000E-27_r8/
      DATA ALGMCS( 14) / -.3401102254316748799999999999999E-29_r8/
      DATA ALGMCS( 15) / +.1276642195630062933333333333333E-30_r8/
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  D9LGMC
      IF (FIRST) THEN
         NALGM = INITDS (ALGMCS, 15, real(D1MACH(3),r8) )
         XBIG = 1.0_r8/SQRT(D1MACH(3))
         XMAX = EXP (MIN(LOG(D1MACH(2)/12._r8),
     >      -LOG(12._r8*D1MACH(1))))
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LT. 10._r8) CALL XERMSG ('SLATEC', 'D9LGMC',
     +   'X MUST BE GE 10', 1, 2)
      IF (X.GE.XMAX) GO TO 20
C
      D9LGMC = 1._r8/(12._r8*X)
      IF (X.LT.XBIG) D9LGMC =
     1   DCSEVL (2.0_r8*(10._r8/X)**2-1._r8, ALGMCS,
     1   NALGM) / X
      RETURN
C
 20   D9LGMC = 0._r8
      CALL XERMSG ('SLATEC', 'D9LGMC', 'X SO BIG D9LGMC UNDERFLOWS', 2,
     +   1)
      RETURN
C
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
! 22Jun2000 fgtok
