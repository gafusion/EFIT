*DECK DBESJ0
      REAL*8 FUNCTION DBESJ0 (X)
C***BEGIN PROLOGUE  DBESJ0
C***PURPOSE  Compute the Bessel function of the first kind of order
C            zero.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10A1
C***TYPE      REAL*8 (BESJ0-S, DBESJ0-D)
C***KEYWORDS  BESSEL FUNCTION, FIRST KIND, FNLIB, ORDER ZERO,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBESJ0(X) calculates the double precision Bessel function of
C the first kind of order zero for double precision argument X.
C
C Series for BJ0        on the interval  0.          to  1.60000E+01
C                                        with weighted error   4.39E-32
C                                         log weighted error  31.36
C                               significant figures required  31.21
C                                    decimal places required  32.00
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, D9B0MP, DCSEVL, INITDS
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DBESJ0
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER ntj0,initds
C============
      REAL*8 X, BJ0CS(19), AMPL, THETA, XSML, Y, D1MACH,
     1  DCSEVL
      LOGICAL FIRST
      SAVE BJ0CS, NTJ0, XSML, FIRST
      DATA BJ0CS(  1) / +.10025416196893913701073127264074E+0_r8/
      DATA BJ0CS(  2) / -.66522300776440513177678757831124E+0_r8/
      DATA BJ0CS(  3) / +.24898370349828131370460468726680E+0_r8/
      DATA BJ0CS(  4) / -.33252723170035769653884341503854E-1_r8/
      DATA BJ0CS(  5) / +.23114179304694015462904924117729E-2_r8/
      DATA BJ0CS(  6) / -.99112774199508092339048519336549E-4_r8/
      DATA BJ0CS(  7) / +.28916708643998808884733903747078E-5_r8/
      DATA BJ0CS(  8) / -.61210858663032635057818407481516E-7_r8/
      DATA BJ0CS(  9) / +.98386507938567841324768748636415E-9_r8/
      DATA BJ0CS( 10) / -.12423551597301765145515897006836E-10_r8/
      DATA BJ0CS( 11) / +.12654336302559045797915827210363E-12_r8/
      DATA BJ0CS( 12) / -.10619456495287244546914817512959E-14_r8/
      DATA BJ0CS( 13) / +.74706210758024567437098915584000E-17_r8/
      DATA BJ0CS( 14) / -.44697032274412780547627007999999E-19_r8/
      DATA BJ0CS( 15) / +.23024281584337436200523093333333E-21_r8/
      DATA BJ0CS( 16) / -.10319144794166698148522666666666E-23_r8/
      DATA BJ0CS( 17) / +.40608178274873322700800000000000E-26_r8/
      DATA BJ0CS( 18) / -.14143836005240913919999999999999E-28_r8/
      DATA BJ0CS( 19) / +.43910905496698880000000000000000E-31_r8/
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DBESJ0
      IF (FIRST) THEN
         NTJ0 = INITDS (BJ0CS, 19, 0.1_r8*real(D1MACH(3),r8))
         XSML = SQRT(8.0_r8*D1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.4.0_r8) GO TO 20
C
      DBESJ0 = 1.0_r8
      IF (Y.GT.XSML) DBESJ0 =
     >   DCSEVL (.125_r8*Y*Y-1._r8, BJ0CS, NTJ0)
      RETURN
C
 20   CALL D9B0MP (Y, AMPL, THETA)
      DBESJ0 = AMPL * COS(THETA)
C
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
! 22Jun2000 fgtok
