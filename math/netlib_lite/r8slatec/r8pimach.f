      FUNCTION r8pimach (DUM)
C
C     THIS SUBPROGRAM SUPPLIES THE VALUE OF THE CONSTANT PI CORRECT TO
C     MACHINE PRECISION WHERE
C
C     PI=3.1415926535897932384626433832795028841971693993751058209749446
C
!============
! idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 dum,r8pimach
!============
C      r8pimach = 3.14159265358979_r8
      r8pimach = 3.1415926535897932384_r8
      RETURN
      END
! 06Jul2005 fgtok -s r8_precision.sub fish_sub.txt "r8con.csh conversion"
