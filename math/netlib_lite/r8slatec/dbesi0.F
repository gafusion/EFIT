*DECK DBESI0
      REAL*8 FUNCTION DBESI0 (X)
C***BEGIN PROLOGUE  DBESI0
C***PURPOSE  Compute the hyperbolic Bessel function of the first kind
C            of order zero.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B1
C***TYPE      REAL*8 (BESI0-S, DBESI0-D)
C***KEYWORDS  FIRST KIND, FNLIB, HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBESI0(X) calculates the double precision modified (hyperbolic)
C Bessel function of the first kind of order zero and double
C precision argument X.
C
C Series for BI0        on the interval  0.          to  9.00000E+00
C                                        with weighted error   9.51E-34
C                                         log weighted error  33.02
C                               significant figures required  33.31
C                                    decimal places required  33.65
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DBSI0E, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DBESI0
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER nti0,initds
C============
      REAL*8 X, BI0CS(18), XMAX, XSML, Y, D1MACH,
     1  DCSEVL, DBSI0E
      LOGICAL FIRST
      SAVE BI0CS, NTI0, XSML, XMAX, FIRST
      DATA BI0CS(  1) / -.7660547252839144951081894976243285E-1_r8/
      DATA BI0CS(  2) / +.1927337953993808269952408750881196E+1_r8/
      DATA BI0CS(  3) / +.2282644586920301338937029292330415E+0_r8/
      DATA BI0CS(  4) / +.1304891466707290428079334210691888E-1_r8/
      DATA BI0CS(  5) / +.4344270900816487451378682681026107E-3_r8/
      DATA BI0CS(  6) / +.9422657686001934663923171744118766E-5_r8/
      DATA BI0CS(  7) / +.1434006289510691079962091878179957E-6_r8/
      DATA BI0CS(  8) / +.1613849069661749069915419719994611E-8_r8/
      DATA BI0CS(  9) / +.1396650044535669699495092708142522E-10_r8/
      DATA BI0CS( 10) / +.9579451725505445344627523171893333E-13_r8/
      DATA BI0CS( 11) / +.5333981859862502131015107744000000E-15_r8/
      DATA BI0CS( 12) / +.2458716088437470774696785919999999E-17_r8/
      DATA BI0CS( 13) / +.9535680890248770026944341333333333E-20_r8/
      DATA BI0CS( 14) / +.3154382039721427336789333333333333E-22_r8/
      DATA BI0CS( 15) / +.9004564101094637431466666666666666E-25_r8/
      DATA BI0CS( 16) / +.2240647369123670016000000000000000E-27_r8/
      DATA BI0CS( 17) / +.4903034603242837333333333333333333E-30_r8/
      DATA BI0CS( 18) / +.9508172606122666666666666666666666E-33_r8/
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DBESI0
      IF (FIRST) THEN
         NTI0 = INITDS (BI0CS, 18, 0.1_r8*real(D1MACH(3),r8))
         XSML = SQRT(4.5_r8*D1MACH(3))
         XMAX = LOG (D1MACH(2))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.3.0_r8) GO TO 20
C
      DBESI0 = 1.0_r8
      IF (Y.GT.XSML) DBESI0 = 2.75_r8+
     1   DCSEVL (Y*Y/4.5_r8-1._r8, BI0CS,
     1   NTI0)
      RETURN
C
 20   IF (Y .GT. XMAX) CALL XERMSG ('SLATEC', 'DBESI0',
     +   'ABS(X) SO BIG I0 OVERFLOWS', 2, 2)
C
      DBESI0 = EXP(Y) * DBSI0E(X)
C
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
! 22Jun2000 fgtok
