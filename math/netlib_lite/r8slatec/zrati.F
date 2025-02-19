*DECK ZRATI
      SUBROUTINE ZRATI (ZR, ZI, FNU, N, CYR, CYI, TOL)
C***BEGIN PROLOGUE  ZRATI
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH, ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CRATI-A, ZRATI-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZRATI COMPUTES RATIOS OF I BESSEL FUNCTIONS BY BACKWARD
C     RECURRENCE.  THE STARTING INDEX IS DETERMINED BY FORWARD
C     RECURRENCE AS DESCRIBED IN J. RES. OF NAT. BUR. OF STANDARDS-B,
C     MATHEMATICAL SCIENCES, VOL 77B, P111-114, SEPTEMBER, 1973,
C     BESSEL FUNCTIONS I AND J OF COMPLEX ARGUMENT AND INTEGER ORDER,
C     BY D. J. SOOKNE.
C
C***SEE ALSO  ZBESH, ZBESI, ZBESK
C***ROUTINES CALLED  ZABS, ZDIV
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZRATI
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      REAL*8 AK, AMAGZ, AP1, AP2, ARG, AZ, CDFNUI, CDFNUR,
     * CONEI, CONER, CYI, CYR, CZEROI, CZEROR, DFNU, FDNU, FLAM, FNU,
     * FNUP, PTI, PTR, P1I, P1R, P2I, P2R, RAK, RAP1, RHO, RT2, RZI,
     * RZR, TEST, TEST1, TOL, TTI, TTR, T1I, T1R, ZI, ZR, ZABS
      INTEGER I, ID, IDNU, INU, ITIME, K, KK, MAGZ, N
      DIMENSION CYR(N), CYI(N)
      EXTERNAL ZABS
      DATA CZEROR,CZEROI,CONER,CONEI,RT2/
     1   0.0_r8, 0.0_r8, 1.0_r8, 0.0_r8,
     2   1.41421356237309505_r8/
C***FIRST EXECUTABLE STATEMENT  ZRATI
      AZ = ZABS(ZR,ZI)
      INU = FNU
      IDNU = INU + N - 1
      MAGZ = AZ
      AMAGZ = MAGZ+1
      FDNU = IDNU
      FNUP = MAX(AMAGZ,FDNU)
      ID = IDNU - MAGZ - 1
      ITIME = 1
      K = 1
      PTR = 1.0_r8/AZ
      RZR = PTR*(ZR+ZR)*PTR
      RZI = -PTR*(ZI+ZI)*PTR
      T1R = RZR*FNUP
      T1I = RZI*FNUP
      P2R = -T1R
      P2I = -T1I
      P1R = CONER
      P1I = CONEI
      T1R = T1R + RZR
      T1I = T1I + RZI
      IF (ID.GT.0) ID = 0
      AP2 = ZABS(P2R,P2I)
      AP1 = ZABS(P1R,P1I)
C-----------------------------------------------------------------------
C     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNU
C     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT
C     P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR
C     PREMATURELY.
C-----------------------------------------------------------------------
      ARG = (AP2+AP2)/(AP1*TOL)
      TEST1 = SQRT(ARG)
      TEST = TEST1
      RAP1 = 1.0_r8/AP1
      P1R = P1R*RAP1
      P1I = P1I*RAP1
      P2R = P2R*RAP1
      P2I = P2I*RAP1
      AP2 = AP2*RAP1
   10 CONTINUE
      K = K + 1
      AP1 = AP2
      PTR = P2R
      PTI = P2I
      P2R = P1R - (T1R*PTR-T1I*PTI)
      P2I = P1I - (T1R*PTI+T1I*PTR)
      P1R = PTR
      P1I = PTI
      T1R = T1R + RZR
      T1I = T1I + RZI
      AP2 = ZABS(P2R,P2I)
      IF (AP1.LE.TEST) GO TO 10
      IF (ITIME.EQ.2) GO TO 20
      AK = ZABS(T1R,T1I)*0.5_r8
      FLAM = AK + SQRT(AK*AK-1.0_r8)
      RHO = MIN(AP2/AP1,FLAM)
      TEST = TEST1*SQRT(RHO/(RHO*RHO-1.0_r8))
      ITIME = 2
      GO TO 10
   20 CONTINUE
      KK = K + 1 - ID
      AK = KK
      T1R = AK
      T1I = CZEROI
      DFNU = FNU + (N-1)
      P1R = 1.0_r8/AP2
      P1I = CZEROI
      P2R = CZEROR
      P2I = CZEROI
      DO 30 I=1,KK
        PTR = P1R
        PTI = P1I
        RAP1 = DFNU + T1R
        TTR = RZR*RAP1
        TTI = RZI*RAP1
        P1R = (PTR*TTR-PTI*TTI) + P2R
        P1I = (PTR*TTI+PTI*TTR) + P2I
        P2R = PTR
        P2I = PTI
        T1R = T1R - CONER
   30 CONTINUE
      IF (P1R.NE.CZEROR .OR. P1I.NE.CZEROI) GO TO 40
      P1R = TOL
      P1I = TOL
   40 CONTINUE
      CALL ZDIV(P2R, P2I, P1R, P1I, CYR(N), CYI(N))
      IF (N.EQ.1) RETURN
      K = N - 1
      AK = K
      T1R = AK
      T1I = CZEROI
      CDFNUR = FNU*RZR
      CDFNUI = FNU*RZI
      DO 60 I=2,N
        PTR = CDFNUR + (T1R*RZR-T1I*RZI) + CYR(K+1)
        PTI = CDFNUI + (T1R*RZI+T1I*RZR) + CYI(K+1)
        AK = ZABS(PTR,PTI)
        IF (AK.NE.CZEROR) GO TO 50
        PTR = TOL
        PTI = TOL
        AK = TOL*RT2
   50   CONTINUE
        RAK = CONER/AK
        CYR(K) = RAK*PTR*RAK
        CYI(K) = -RAK*PTI*RAK
        T1R = T1R - CONER
        K = K - 1
   60 CONTINUE
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
