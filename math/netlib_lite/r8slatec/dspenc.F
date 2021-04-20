*DECK DSPENC
      REAL*8 FUNCTION DSPENC (X)
C***BEGIN PROLOGUE  DSPENC
C***PURPOSE  Compute a form of Spence's integral due to K. Mitchell.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C5
C***TYPE      REAL*8 (SPENC-S, DSPENC-D)
C***KEYWORDS  FNLIB, SPECIAL FUNCTIONS, SPENCE'S INTEGRAL
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DSPENC(X) calculates the double precision Spence's integral
C for double precision argument X.  Spence's function defined by
C        integral from 0 to X of  -LOG(1-Y)/Y  DY.
C For ABS(X) .LE. 1, the uniformly convergent expansion
C        DSPENC = sum K=1,infinity  X**K / K**2     is valid.
C This is a form of Spence's integral due to K. Mitchell which differs
C from the definition in the NBS Handbook of Mathematical Functions.
C
C Spence's function can be used to evaluate much more general integral
C forms.  For example,
C        integral from 0 to Z of  LOG(A*X+B)/(C*X+D)  DX  =
C             LOG(ABS(B-A*D/C))*LOG(ABS(A*(C*X+D)/(A*D-B*C)))/C
C             - DSPENC (A*(C*Z+D)/(A*D-B*C)) / C.
C
C Ref -- K. Mitchell, Philosophical Magazine, 40, p.351 (1949).
C        Stegun and Abromowitz, AMS 55, p.1004.
C
C
C Series for SPEN       on the interval  0.          to  5.00000E-01
C                                        with weighted error   4.74E-32
C                                         log weighted error  31.32
C                               significant figures required  30.37
C                                    decimal places required  32.11
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DCSEVL, INITDS
C***REVISION HISTORY  (YYMMDD)
C   780201  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891115  Corrected third argument in reference to INITDS.  (WRB)
C   891115  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DSPENC
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER nspenc,initds
C============
      REAL*8 X, SPENCS(38), ALN, PI26, XBIG, D1MACH, DCSEVL
      LOGICAL FIRST
      SAVE SPENCS, PI26, NSPENC, XBIG, FIRST
      DATA SPENCS(  1) / +.1527365598892405872946684910028E+0_r8/
      DATA SPENCS(  2) / +.8169658058051014403501838185271E-1_r8/
      DATA SPENCS(  3) / +.5814157140778730872977350641182E-2_r8/
      DATA SPENCS(  4) / +.5371619814541527542247889005319E-3_r8/
      DATA SPENCS(  5) / +.5724704675185826233210603054782E-4_r8/
      DATA SPENCS(  6) / +.6674546121649336343607835438589E-5_r8/
      DATA SPENCS(  7) / +.8276467339715676981584391689011E-6_r8/
      DATA SPENCS(  8) / +.1073315673030678951270005873354E-6_r8/
      DATA SPENCS(  9) / +.1440077294303239402334590331513E-7_r8/
      DATA SPENCS( 10) / +.1984442029965906367898877139608E-8_r8/
      DATA SPENCS( 11) / +.2794005822163638720201994821615E-9_r8/
      DATA SPENCS( 12) / +.4003991310883311823072580445908E-10_r8/
      DATA SPENCS( 13) / +.5823462892044638471368135835757E-11_r8/
      DATA SPENCS( 14) / +.8576708692638689278097914771224E-12_r8/
      DATA SPENCS( 15) / +.1276862586280193045989483033433E-12_r8/
      DATA SPENCS( 16) / +.1918826209042517081162380416062E-13_r8/
      DATA SPENCS( 17) / +.2907319206977138177795799719673E-14_r8/
      DATA SPENCS( 18) / +.4437112685276780462557473641745E-15_r8/
      DATA SPENCS( 19) / +.6815727787414599527867359135607E-16_r8/
      DATA SPENCS( 20) / +.1053017386015574429547019416644E-16_r8/
      DATA SPENCS( 21) / +.1635389806752377100051821734570E-17_r8/
      DATA SPENCS( 22) / +.2551852874940463932310901642581E-18_r8/
      DATA SPENCS( 23) / +.3999020621999360112770470379519E-19_r8/
      DATA SPENCS( 24) / +.6291501645216811876514149171199E-20_r8/
      DATA SPENCS( 25) / +.9933827435675677643803887752533E-21_r8/
      DATA SPENCS( 26) / +.1573679570749964816721763805866E-21_r8/
      DATA SPENCS( 27) / +.2500595316849476129369270954666E-22_r8/
      DATA SPENCS( 28) / +.3984740918383811139210663253333E-23_r8/
      DATA SPENCS( 29) / +.6366473210082843892691326293333E-24_r8/
      DATA SPENCS( 30) / +.1019674287239678367077061973333E-24_r8/
      DATA SPENCS( 31) / +.1636881058913518841111074133333E-25_r8/
      DATA SPENCS( 32) / +.2633310439417650117345279999999E-26_r8/
      DATA SPENCS( 33) / +.4244811560123976817224362666666E-27_r8/
      DATA SPENCS( 34) / +.6855411983680052916824746666666E-28_r8/
      DATA SPENCS( 35) / +.1109122433438056434018986666666E-28_r8/
      DATA SPENCS( 36) / +.1797431304999891457365333333333E-29_r8/
      DATA SPENCS( 37) / +.2917505845976095173290666666666E-30_r8/
      DATA SPENCS( 38) / +.4742646808928671061333333333333E-31_r8/
      DATA PI26 / +1.64493406684822643647241516664602518921_r8/
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DSPENC
      IF (FIRST) THEN
         NSPENC = INITDS (SPENCS, 38, 0.1_r8*real(D1MACH(3),r8))
         XBIG = 1.0_r8/D1MACH(3)
      ENDIF
      FIRST = .FALSE.
C
      IF (X.GT.2.0_r8) GO TO 60
      IF (X.GT.1.0_r8) GO TO 50
      IF (X.GT.0.5_r8) GO TO 40
      IF (X.GE.0.0_r8) GO TO 30
      IF (X.GT.(-1._r8)) GO TO 20
C
C HERE IF X .LE. -1.0
C
      ALN = LOG(1.0_r8-X)
      DSPENC = -PI26 - 0.5_r8*ALN*(2.0_r8*LOG(-X)-ALN)
      IF (X.GT.(-XBIG)) DSPENC = DSPENC
     1  + (1._r8+ DCSEVL (4._r8/(1._r8-X)-1._r8,
     1   SPENCS, NSPENC))/(1._r8-X)
      RETURN
C
C -1.0 .LT. X .LT. 0.0
C
 20   DSPENC = -0.5_r8*LOG(1.0_r8-X)**2
     1  - X*(1._r8+DCSEVL(4._r8*X/(X-1._r8)
     2   -1._r8, SPENCS, NSPENC))/(X-1._r8)
      RETURN
C
C 0.0 .LE. X .LE. 0.5
C
 30   DSPENC = X*(1._r8+ DCSEVL (4._r8*X-1._r8,
     1   SPENCS, NSPENC))
      RETURN
C
C 0.5 .LT. X .LE. 1.0
C
 40   DSPENC = PI26
      IF (X.NE.1._r8) DSPENC = PI26 - LOG(X)*LOG(1.0_r8-X)
     1  - (1._r8-X)*(1._r8+
     1   DCSEVL(4._r8*(1._r8-X)-1._r8, SPENCS,
     1   NSPENC))
      RETURN
C
C 1.0 .LT. X .LE. 2.0
C
 50   DSPENC = PI26 - 0.5_r8*LOG(X)*LOG((X-1._r8)**2/X)
     1  + (X-1._r8)*(1._r8+DCSEVL(4._r8*
     1   (X-1._r8)/X-1._r8, SPENCS, NSPENC))/X
      RETURN
C
C X .GT. 2.0
C
 60   DSPENC = 2.0_r8*PI26 - 0.5_r8*LOG(X)**2
      IF (X.LT.XBIG) DSPENC = DSPENC
     1  - (1._r8+
     2   DCSEVL (4._r8/X-1._r8, SPENCS, NSPENC))/X
      RETURN
C
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
! 22Jun2000 fgtok
