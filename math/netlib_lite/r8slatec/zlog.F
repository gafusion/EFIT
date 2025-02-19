*DECK ZLOG
      SUBROUTINE ZLOG (AR, AI, BR, BI, IERR)
C***BEGIN PROLOGUE  ZLOG
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
C            ZBIRY
C***LIBRARY   SLATEC
C***TYPE      ALL (ZLOG-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     REAL*8 COMPLEX LOGARITHM B=CLOG(A)
C     IERR=0,NORMAL RETURN      IERR=1, Z=CMPLX(0.0,0.0)
C***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
C***ROUTINES CALLED  ZABS
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZLOG
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      REAL*8 AR, AI, BR, BI, ZM, DTHETA, DPI, DHPI
      REAL*8 ZABS
      INTEGER IERR
      EXTERNAL ZABS
      DATA DPI , DHPI  / 3.141592653589793238462643383E0_r8,
     1                   1.570796326794896619231321696E0_r8/
C***FIRST EXECUTABLE STATEMENT  ZLOG
      IERR=0
      IF (AR.EQ.0.0E0_r8) GO TO 10
      IF (AI.EQ.0.0E0_r8) GO TO 20
      DTHETA = atan(AI/AR)
      IF (DTHETA.LE.0.0E0_r8) GO TO 40
      IF (AR.LT.0.0E0_r8) DTHETA = DTHETA - DPI
      GO TO 50
   10 IF (AI.EQ.0.0E0_r8) GO TO 60
      BI = DHPI
      BR = LOG(ABS(AI))
      IF (AI.LT.0.0E0_r8) BI = -BI
      RETURN
   20 IF (AR.GT.0.0E0_r8) GO TO 30
      BR = LOG(ABS(AR))
      BI = DPI
      RETURN
   30 BR = LOG(AR)
      BI = 0.0E0_r8
      RETURN
   40 IF (AR.LT.0.0E0_r8) DTHETA = DTHETA + DPI
   50 ZM = ZABS(AR,AI)
      BR = LOG(ZM)
      BI = DTHETA
      RETURN
   60 CONTINUE
      IERR=1
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
