*DECK DASYIK
      SUBROUTINE DASYIK (X, FNU, KODE, FLGIK, RA, ARG, IN, Y)
C***BEGIN PROLOGUE  DASYIK
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBESI and DBESK
C***LIBRARY   SLATEC
C***TYPE      REAL*8 (ASYIK-S, DASYIK-D)
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C                    DASYIK computes Bessel functions I and K
C                  for arguments X.GT.0.0 and orders FNU.GE.35
C                  on FLGIK = 1 and FLGIK = -1 respectively.
C
C                                    INPUT
C
C      X    - Argument, X.GT.0.0D0
C      FNU  - Order of first Bessel function
C      KODE - A parameter to indicate the scaling option
C             KODE=1 returns Y(I)=        I/SUB(FNU+I-1)/(X), I=1,IN
C                    or      Y(I)=        K/SUB(FNU+I-1)/(X), I=1,IN
C                    on FLGIK = 1.0D0 or FLGIK = -1.0D0
C             KODE=2 returns Y(I)=EXP(-X)*I/SUB(FNU+I-1)/(X), I=1,IN
C                    or      Y(I)=EXP( X)*K/SUB(FNU+I-1)/(X), I=1,IN
C                    on FLGIK = 1.0D0 or FLGIK = -1.0D0
C     FLGIK - Selection parameter for I or K FUNCTION
C             FLGIK =  1.0D0 gives the I function
C             FLGIK = -1.0D0 gives the K function
C        RA - SQRT(1.+Z*Z), Z=X/FNU
C       ARG - Argument of the leading exponential
C        IN - Number of functions desired, IN=1 or 2
C
C                                    OUTPUT
C
C         Y - A vector whose first IN components contain the sequence
C
C     Abstract  **** A double precision routine ****
C         DASYIK implements the uniform asymptotic expansion of
C         the I and K Bessel functions for FNU.GE.35 and real
C         X.GT.0.0D0. The forms are identical except for a change
C         in sign of some of the terms. This change in sign is
C         accomplished by means of the FLAG FLGIK = 1 or -1.
C
C***SEE ALSO  DBESI, DBESK
C***ROUTINES CALLED  D1MACH
C***REVISION HISTORY  (YYMMDD)
C   750101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910408  Updated the AUTHOR section.  (WRB)
C***END PROLOGUE  DASYIK
C
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER IN, J, JN, K, KK, KODE, L
      REAL*8 AK,AP,ARG,C,COEF,CON,ETX,FLGIK,FN,FNU,GLN,RA,
     1 S1, S2, T, TOL, T2, X, Y, Z
      REAL*8 D1MACH
      DIMENSION Y(*), C(65), CON(2)
      SAVE CON, C
      DATA CON(1), CON(2)  /
     1   3.98942280401432678E-01_r8,1.25331413731550025E00_r8/
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2     C(19), C(20), C(21), C(22), C(23), C(24)/
     3    -2.08333333333333E-01_r8,     1.25000000000000E-01_r8,
     4     3.34201388888889E-01_r8,    -4.01041666666667E-01_r8,
     5     7.03125000000000E-02_r8,    -1.02581259645062E00_r8,
     6     1.84646267361111E00_r8,    -8.91210937500000E-01_r8,
     7     7.32421875000000E-02_r8,     4.66958442342625E00_r8,
     8    -1.12070026162230E01_r8,     8.78912353515625E00_r8,
     9    -2.36408691406250E00_r8,     1.12152099609375E-01_r8,
     1    -2.82120725582002E01_r8,     8.46362176746007E01_r8,
     2    -9.18182415432400E01_r8,     4.25349987453885E01_r8,
     3    -7.36879435947963E00_r8,     2.27108001708984E-01_r8,
     4     2.12570130039217E02_r8,    -7.65252468141182E02_r8,
     5     1.05999045252800E03_r8,    -6.99579627376133E02_r8/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3     2.18190511744212E02_r8,    -2.64914304869516E01_r8,
     4     5.72501420974731E-01_r8,    -1.91945766231841E03_r8,
     5     8.06172218173731E03_r8,    -1.35865500064341E04_r8,
     6     1.16553933368645E04_r8,    -5.30564697861340E03_r8,
     7     1.20090291321635E03_r8,    -1.08090919788395E02_r8,
     8     1.72772750258446E00_r8,     2.02042913309661E04_r8,
     9    -9.69805983886375E04_r8,     1.92547001232532E05_r8,
     1    -2.03400177280416E05_r8,     1.22200464983017E05_r8,
     2    -4.11926549688976E04_r8,     7.10951430248936E03_r8,
     3    -4.93915304773088E02_r8,     6.07404200127348E00_r8,
     4    -2.42919187900551E05_r8,     1.31176361466298E06_r8,
     5    -2.99801591853811E06_r8,     3.76327129765640E06_r8/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2     C(65)/
     3    -2.81356322658653E06_r8,     1.26836527332162E06_r8,
     4    -3.31645172484564E05_r8,     4.52187689813627E04_r8,
     5    -2.49983048181121E03_r8,     2.43805296995561E01_r8,
     6     3.28446985307204E06_r8,    -1.97068191184322E07_r8,
     7     5.09526024926646E07_r8,    -7.41051482115327E07_r8,
     8     6.63445122747290E07_r8,    -3.75671766607634E07_r8,
     9     1.32887671664218E07_r8,    -2.78561812808645E06_r8,
     1     3.08186404612662E05_r8,    -1.38860897537170E04_r8,
     2     1.10017140269247E02_r8/
C***FIRST EXECUTABLE STATEMENT  DASYIK
      TOL = D1MACH(3)
      TOL = MAX(TOL,1.0E-15_r8)
      FN = FNU
      Z  = (3.0_r8-FLGIK)/2.0_r8
      KK = INT(Z)
      DO 50 JN=1,IN
        IF (JN.EQ.1) GO TO 10
        FN = FN - FLGIK
        Z = X/FN
        RA = SQRT(1.0_r8+Z*Z)
        GLN = LOG((1.0_r8+RA)/Z)
        ETX = KODE - 1
        T = RA*(1.0_r8-ETX) + ETX/(Z+RA)
        ARG = FN*(T-GLN)*FLGIK
   10   COEF = EXP(ARG)
        T = 1.0_r8/RA
        T2 = T*T
        T = T/FN
        T = SIGN(T,FLGIK)
        S2 = 1.0_r8
        AP = 1.0_r8
        L = 0
        DO 30 K=2,11
          L = L + 1
          S1 = C(L)
          DO 20 J=2,K
            L = L + 1
            S1 = S1*T2 + C(L)
   20     CONTINUE
          AP = AP*T
          AK = AP*S1
          S2 = S2 + AK
          IF (MAX(ABS(AK),ABS(AP)) .LT.TOL) GO TO 40
   30   CONTINUE
   40   CONTINUE
      T = ABS(T)
      Y(JN) = S2*COEF*SQRT(T)*CON(KK)
   50 CONTINUE
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
