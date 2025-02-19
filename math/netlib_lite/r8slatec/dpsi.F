*DECK DPSI
      REAL*8 FUNCTION DPSI (X)
C***BEGIN PROLOGUE  DPSI
C***PURPOSE  Compute the Psi (or Digamma) function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7C
C***TYPE      REAL*8 (PSI-S, DPSI-D, CPSI-C)
C***KEYWORDS  DIGAMMA FUNCTION, FNLIB, PSI FUNCTION, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DPSI calculates the double precision Psi (or Digamma) function for
C double precision argument X.  PSI(X) is the logarithmic derivative
C of the Gamma function of X.
C
C Series for PSI        on the interval  0.          to  1.00000E+00
C                                        with weighted error   5.79E-32
C                                         log weighted error  31.24
C                               significant figures required  30.93
C                                    decimal places required  32.05
C
C
C Series for APSI       on the interval  0.          to  1.00000E-02
C                                        with weighted error   7.75E-33
C                                         log weighted error  32.11
C                               significant figures required  28.88
C                                    decimal places required  32.71
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DCOT, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900727  Added EXTERNAL statement.  (WRB)
C   920618  Removed space from variable name.  (RWC, WRB)
C***END PROLOGUE  DPSI
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER ntpsi,ntapsi,initds,n,i
C============
      REAL*8 X, PSICS(42), APSICS(16), AUX, DXREL, PI, XBIG,
     1  Y, DCOT, DCSEVL, D1MACH
      LOGICAL FIRST
      EXTERNAL DCOT
      SAVE PSICS, APSICS, PI, NTPSI, NTAPSI, XBIG, DXREL, FIRST
      DATA PSICS(  1) / -.38057080835217921520437677667039E-1_r8/
      DATA PSICS(  2) / +.49141539302938712748204699654277E+0_r8/
      DATA PSICS(  3) / -.56815747821244730242892064734081E-1_r8/
      DATA PSICS(  4) / +.83578212259143131362775650747862E-2_r8/
      DATA PSICS(  5) / -.13332328579943425998079274172393E-2_r8/
      DATA PSICS(  6) / +.22031328706930824892872397979521E-3_r8/
      DATA PSICS(  7) / -.37040238178456883592889086949229E-4_r8/
      DATA PSICS(  8) / +.62837936548549898933651418717690E-5_r8/
      DATA PSICS(  9) / -.10712639085061849855283541747074E-5_r8/
      DATA PSICS( 10) / +.18312839465484165805731589810378E-6_r8/
      DATA PSICS( 11) / -.31353509361808509869005779796885E-7_r8/
      DATA PSICS( 12) / +.53728087762007766260471919143615E-8_r8/
      DATA PSICS( 13) / -.92116814159784275717880632624730E-9_r8/
      DATA PSICS( 14) / +.15798126521481822782252884032823E-9_r8/
      DATA PSICS( 15) / -.27098646132380443065440589409707E-10_r8/
      DATA PSICS( 16) / +.46487228599096834872947319529549E-11_r8/
      DATA PSICS( 17) / -.79752725638303689726504797772737E-12_r8/
      DATA PSICS( 18) / +.13682723857476992249251053892838E-12_r8/
      DATA PSICS( 19) / -.23475156060658972717320677980719E-13_r8/
      DATA PSICS( 20) / +.40276307155603541107907925006281E-14_r8/
      DATA PSICS( 21) / -.69102518531179037846547422974771E-15_r8/
      DATA PSICS( 22) / +.11856047138863349552929139525768E-15_r8/
      DATA PSICS( 23) / -.20341689616261559308154210484223E-16_r8/
      DATA PSICS( 24) / +.34900749686463043850374232932351E-17_r8/
      DATA PSICS( 25) / -.59880146934976711003011081393493E-18_r8/
      DATA PSICS( 26) / +.10273801628080588258398005712213E-18_r8/
      DATA PSICS( 27) / -.17627049424561071368359260105386E-19_r8/
      DATA PSICS( 28) / +.30243228018156920457454035490133E-20_r8/
      DATA PSICS( 29) / -.51889168302092313774286088874666E-21_r8/
      DATA PSICS( 30) / +.89027730345845713905005887487999E-22_r8/
      DATA PSICS( 31) / -.15274742899426728392894971904000E-22_r8/
      DATA PSICS( 32) / +.26207314798962083136358318079999E-23_r8/
      DATA PSICS( 33) / -.44964642738220696772598388053333E-24_r8/
      DATA PSICS( 34) / +.77147129596345107028919364266666E-25_r8/
      DATA PSICS( 35) / -.13236354761887702968102638933333E-25_r8/
      DATA PSICS( 36) / +.22709994362408300091277311999999E-26_r8/
      DATA PSICS( 37) / -.38964190215374115954491391999999E-27_r8/
      DATA PSICS( 38) / +.66851981388855302310679893333333E-28_r8/
      DATA PSICS( 39) / -.11469986654920864872529919999999E-28_r8/
      DATA PSICS( 40) / +.19679385886541405920515413333333E-29_r8/
      DATA PSICS( 41) / -.33764488189750979801907200000000E-30_r8/
      DATA PSICS( 42) / +.57930703193214159246677333333333E-31_r8/
      DATA APSICS(  1) / -.832710791069290760174456932269E-3_r8/
      DATA APSICS(  2) / -.416251842192739352821627121990E-3_r8/
      DATA APSICS(  3) / +.103431560978741291174463193961E-6_r8/
      DATA APSICS(  4) / -.121468184135904152987299556365E-9_r8/
      DATA APSICS(  5) / +.311369431998356155521240278178E-12_r8/
      DATA APSICS(  6) / -.136461337193177041776516100945E-14_r8/
      DATA APSICS(  7) / +.902051751315416565130837974000E-17_r8/
      DATA APSICS(  8) / -.831542997421591464829933635466E-19_r8/
      DATA APSICS(  9) / +.101224257073907254188479482666E-20_r8/
      DATA APSICS( 10) / -.156270249435622507620478933333E-22_r8/
      DATA APSICS( 11) / +.296542716808903896133226666666E-24_r8/
      DATA APSICS( 12) / -.674686886765702163741866666666E-26_r8/
      DATA APSICS( 13) / +.180345311697189904213333333333E-27_r8/
      DATA APSICS( 14) / -.556901618245983607466666666666E-29_r8/
      DATA APSICS( 15) / +.195867922607736251733333333333E-30_r8/
      DATA APSICS( 16) / -.775195892523335680000000000000E-32_r8/
      DATA PI / 3.14159265358979323846264338327950_r8/
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DPSI
      IF (FIRST) THEN
         NTPSI = INITDS (PSICS, 42, 0.1_r8*real(D1MACH(3),r8) )
         NTAPSI = INITDS (APSICS, 16, 0.1_r8*real(D1MACH(3),r8) )
C
         XBIG = 1.0_r8/SQRT(D1MACH(3))
         DXREL = SQRT(D1MACH(4))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
C
      IF (Y.GT.10.0_r8) GO TO 50
C
C DPSI(X) FOR ABS(X) .LE. 2
C
      N = X
      IF (X.LT.0._r8) N = N - 1
      Y = X - N
      N = N - 1
      DPSI = DCSEVL (2._r8*Y-1._r8, PSICS, NTPSI)
      IF (N.EQ.0) RETURN
C
      IF (N.GT.0) GO TO 30
C
      N = -N
      IF (X .EQ. 0._r8)
     +   CALL XERMSG ('SLATEC', 'DPSI', 'X IS 0', 2, 2)
      IF (X .LT. 0._r8.AND. X+N-2 .EQ. 0._r8)
     +   CALL XERMSG ('SLATEC',
     +   'DPSI', 'X IS A NEGATIVE INTEGER', 3, 2)
      IF (X .LT. (-0.5_r8) .AND.
     +   ABS((X-AINT(X-0.5_r8))/X) .LT. DXREL)
     +   CALL XERMSG ('SLATEC', 'DPSI',
     +   'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER',
     +   1, 1)
C
      DO 20 I=1,N
        DPSI = DPSI - 1._r8/(X+I-1)
 20   CONTINUE
      RETURN
C
C DPSI(X) FOR X .GE. 2.0 AND X .LE. 10.0
C
 30   DO 40 I=1,N
        DPSI = DPSI + 1.0_r8/(Y+I)
 40   CONTINUE
      RETURN
C
C DPSI(X) FOR ABS(X) .GT. 10.0
C
 50   AUX = 0._r8
      IF (Y.LT.XBIG) AUX =
     1   DCSEVL (2._r8*(10._r8/Y)**2-1._r8, APSICS,
     1   NTAPSI)
C
      IF (X.LT.0._r8) DPSI = LOG(ABS(X)) - 0.5_r8/X + AUX
     1  - PI*DCOT(PI*X)
      IF (X.GT.0._r8) DPSI = LOG(X) - 0.5_r8/X + AUX
      RETURN
C
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
! 22Jun2000 fgtok
