*DECK DXPSI
      REAL*8 FUNCTION DXPSI (A, IPSIK, IPSIX)
C***BEGIN PROLOGUE  DXPSI
C***SUBSIDIARY
C***PURPOSE  To compute values of the Psi function for DXLEGF.
C***LIBRARY   SLATEC
C***CATEGORY  C7C
C***TYPE      REAL*8 (XPSI-S, DXPSI-D)
C***KEYWORDS  PSI FUNCTION
C***AUTHOR  Smith, John M., (NBS and George Mason University)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   820728  DATE WRITTEN
C   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
C   901019  Revisions to prologue.  (DWL and WRB)
C   901106  Changed all specific intrinsics to generic.  (WRB)
C           Corrected order of sections in prologue and added TYPE
C           section.  (WRB)
C   920127  Revised PURPOSE section of prologue.  (DWL)
C***END PROLOGUE  DXPSI
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER ipsix,ipsik,n,k1,i,k,m
C============
      REAL*8 A,B,C,CNUM,CDENOM
      DIMENSION CNUM(12),CDENOM(12)
      SAVE CNUM, CDENOM
C
C        CNUM(I) AND CDENOM(I) ARE THE ( REDUCED ) NUMERATOR
C        AND 2*I*DENOMINATOR RESPECTIVELY OF THE 2*I TH BERNOULLI
C        NUMBER.
C
      DATA CNUM(1),CNUM(2),CNUM(3),CNUM(4),CNUM(5),CNUM(6),CNUM(7),
     1CNUM(8),CNUM(9),CNUM(10),CNUM(11),CNUM(12)
     2 / 1._r8,-1._r8,   1._r8,-1._r8,
     3   1._r8,-691._r8, 1._r8,-3617._r8,
     4   43867._r8, -174611._r8, 77683._r8,
     5   -236364091._r8/
      DATA CDENOM(1),CDENOM(2),CDENOM(3),CDENOM(4),CDENOM(5),CDENOM(6),
     1   CDENOM(7),CDENOM(8),CDENOM(9),CDENOM(10),CDENOM(11),CDENOM(12)
     2  /12._r8,120._r8,   252._r8,
     3   240._r8,132._r8,
     4   32760._r8, 12._r8,  8160._r8,
     5   14364._r8, 6600._r8, 276._r8,
     6   65520._r8/
C***FIRST EXECUTABLE STATEMENT  DXPSI
      N=MAX(0,IPSIX-INT(A))
      B=N+A
      K1=IPSIK-1
C
C        SERIES EXPANSION FOR A .GT. IPSIX USING IPSIK-1 TERMS.
C
      C=0._r8
      DO 12 I=1,K1
      K=IPSIK-I
   12 C=(C+CNUM(K)/CDENOM(K))/B**2
      DXPSI=LOG(B)-(C+.5_r8/B)
      IF(N.EQ.0) GO TO 20
      B=0._r8
C
C        RECURRENCE FOR A .LE. IPSIX.
C
      DO 15 M=1,N
   15 B=B+1._r8/(N-M+A)
      DXPSI=DXPSI-B
   20 RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
