*DECK DBSKES
      SUBROUTINE DBSKES (XNU, X, NIN, BKE)
C***BEGIN PROLOGUE  DBSKES
C***PURPOSE  Compute a sequence of exponentially scaled modified Bessel
C            functions of the third kind of fractional order.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B3
C***TYPE      REAL*8 (BESKES-S, DBSKES-D)
C***KEYWORDS  EXPONENTIALLY SCALED, FNLIB, FRACTIONAL ORDER,
C             MODIFIED BESSEL FUNCTION, SEQUENCE OF BESSEL FUNCTIONS,
C             SPECIAL FUNCTIONS, THIRD KIND
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBSKES(XNU,X,NIN,BKE) computes a double precision sequence
C of exponentially scaled modified Bessel functions
C of the third kind of order XNU + I at X, where X .GT. 0,
C XNU lies in (-1,1), and I = 0, 1, ... , NIN - 1, if NIN is positive
C and I = 0, -1, ... , NIN + 1, if NIN is negative.  On return, the
C vector BKE(.) contains the results at X for order starting at XNU.
C XNU, X, and BKE are double precision.  NIN is integer.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, D9KNUS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DBSKES
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER nin,n,iswtch,i
C============
      REAL*8 XNU, X, BKE(*), BKNU1, V, VINCR, VEND, ALNBIG,
     1  D1MACH, DIRECT
      SAVE ALNBIG
      DATA ALNBIG / 0._r8/
C***FIRST EXECUTABLE STATEMENT  DBSKES
      IF (ALNBIG.EQ.0._r8) ALNBIG = LOG (D1MACH(2))
C
      V = ABS(XNU)
      N = ABS(NIN)
C
      IF (V .GE. 1._r8) CALL XERMSG ('SLATEC', 'DBSKES',
     +   'ABS(XNU) MUST BE LT 1', 2, 2)
      IF (X .LE. 0._r8) CALL XERMSG ('SLATEC', 'DBSKES',
     +   'X IS LE 0', 3, 2)
      IF (N .EQ. 0) CALL XERMSG ('SLATEC', 'DBSKES',
     +   'N THE NUMBER IN THE SEQUENCE IS 0', 4, 2)
C
      CALL D9KNUS (V, X, BKE(1), BKNU1, ISWTCH)
      IF (N.EQ.1) RETURN
C
      VINCR = SIGN (1.0_r8, real(NIN,r8))
      DIRECT = VINCR
      IF (XNU.NE.0._r8) DIRECT = VINCR*SIGN(1._r8, XNU)
      IF (ISWTCH .EQ. 1 .AND. DIRECT .GT. 0._r8)
     +   CALL XERMSG ('SLATEC',
     +   'DBSKES', 'X SO SMALL BESSEL K-SUB-XNU+1 OVERFLOWS', 5, 2)
      BKE(2) = BKNU1
C
      IF (DIRECT.LT.0._r8)
     1   CALL D9KNUS (ABS(XNU+VINCR), X, BKE(2), BKNU1,
     1   ISWTCH)
      IF (N.EQ.2) RETURN
C
      VEND = ABS (XNU+NIN) - 1.0_r8
      IF ((VEND-.5_r8)*LOG(VEND)+0.27_r8-
     +   VEND*(LOG(X)-.694_r8) .GT.
     +   ALNBIG) CALL XERMSG ('SLATEC', 'DBSKES',
     +      'X SO SMALL OR ABS(NU) SO BIG THAT BESSEL K-SUB-NU ' //
     +      'OVERFLOWS', 5, 2)
C
      V = XNU
      DO 10 I=3,N
        V = V + VINCR
        BKE(I) = 2.0_r8*V*BKE(I-1)/X + BKE(I-2)
 10   CONTINUE
C
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
! 22Jun2000 fgtok
