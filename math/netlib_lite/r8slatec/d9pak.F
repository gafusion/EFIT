*DECK D9PAK
      REAL*8 FUNCTION D9PAK (Y, N)
C***BEGIN PROLOGUE  D9PAK
C***PURPOSE  Pack a base 2 exponent into a floating point number.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  A6B
C***TYPE      REAL*8 (R9PAK-S, D9PAK-D)
C***KEYWORDS  FNLIB, PACK
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Pack a base 2 exponent into floating point number X.  This routine is
C almost the inverse of D9UPAK.  It is not exactly the inverse, because
C ABS(X) need not be between 0.5 and 1.0.  If both D9PAK and 2.d0**N
C were known to be in range we could compute
C               D9PAK = X *2.0d0**N
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, D9UPAK, J1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891009  Corrected error when XERROR called.  (WRB)
C   891009  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   901009  Routine used J1MACH(7) where it should use J1MACH(10),
C           Corrected (RWC)
C***END PROLOGUE  D9PAK
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER n,nmin,nmax,j1mach,ny,nsum
C============
      REAL*8 Y, A1N2B,A1N210,D1MACH
      LOGICAL FIRST
      SAVE NMIN, NMAX, A1N210, FIRST
      DATA A1N210 / 3.321928094887362347870319429489_r8/
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  D9PAK
      IF (FIRST) THEN
         A1N2B = 1.0_r8
         IF(J1MACH(10).NE.2) A1N2B=D1MACH(5)*A1N210
         NMIN = A1N2B*J1MACH(15)
         NMAX = A1N2B*J1MACH(16)
      ENDIF
      FIRST = .FALSE.
C
      CALL D9UPAK(Y,D9PAK,NY)
C
      NSUM=N+NY
      IF(NSUM.LT.NMIN)GO TO 40
      IF (NSUM .GT. NMAX) CALL XERMSG ('SLATEC', 'D9PAK',
     +   'PACKED NUMBER OVERFLOWS', 1, 2)
C
      IF (NSUM.EQ.0) RETURN
      IF(NSUM.GT.0) GO TO 30
C
 20   D9PAK = 0.5_r8*D9PAK
      NSUM=NSUM+1
      IF(NSUM.NE.0) GO TO 20
      RETURN
C
 30   D9PAK = 2.0_r8*D9PAK
      NSUM=NSUM - 1
      IF (NSUM.NE.0) GO TO 30
      RETURN
C
 40   CALL XERMSG ('SLATEC', 'D9PAK', 'PACKED NUMBER UNDERFLOWS', 1, 1)
      D9PAK = 0.0_r8
      RETURN
C
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
