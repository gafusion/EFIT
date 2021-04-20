      SUBROUTINE r8cosgen (N,IJUMP,FNUM,FDEN,A)
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER ijump,n,k3,k4,k,k1,k5,i,k2,np1
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 fnum,fden,a,pi,r8pimach,dum,pibyn,x,y
!============
      DIMENSION       A(1)
C
C
C     THIS SUBROUTINE COMPUTES REQUIRED COSINE VALUES IN ASCENDING
C     ORDER.  WHEN IJUMP .GT. 1 THE ROUTINE COMPUTES VALUES
C
C        2*COS(J*PI/L) , J=1,2,...,L AND J .NE. 0(MOD N/IJUMP+1)
C
C     WHERE L = IJUMP*(N/IJUMP+1).
C
C
C     WHEN IJUMP = 1 IT COMPUTES
C
C            2*COS((J-FNUM)*PI/(N+FDEN)) ,  J=1, 2, ... ,N
C
C     WHERE
C        FNUM = 0.5, FDEN = 0.0,  FOR REGULAR REDUCTION VALUES
C        FNUM = 0.0, FDEN = 1.0, FOR B-R AND C-R WHEN ISTAG = 1
C        FNUM = 0.0, FDEN = 0.5, FOR B-R AND C-R WHEN ISTAG = 2
C        FNUM = 0.5, FDEN = 0.5, FOR B-R AND C-R WHEN ISTAG = 2
C                                IN POISN2 ONLY.
C
C
      PI = r8pimach(DUM)
      IF (N .EQ. 0) GO TO 105
      IF (IJUMP .EQ. 1) GO TO 103
      K3 = N/IJUMP+1
      K4 = K3-1
      PIBYN = PI/real(N+IJUMP,r8)
      DO 102 K=1,IJUMP
         K1 = (K-1)*K3
         K5 = (K-1)*K4
         DO 101 I=1,K4
            X = K1+I
            K2 = K5+I
            A(K2) = -2._r8*COS(X*PIBYN)
  101    CONTINUE
  102 CONTINUE
      GO TO 105
  103 CONTINUE
      NP1 = N+1
      Y = PI/(real(N,r8)+FDEN)
      DO 104 I=1,N
         X = real(NP1-I,r8)-FNUM
         A(I) = 2._r8*COS(X*Y)
  104 CONTINUE
  105 CONTINUE
      RETURN
      END
! 06Jul2005 fgtok -s r8_precision.sub fish_sub.txt "r8con.csh conversion"
! 06Jul2005 fgtok
