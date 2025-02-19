*DECK ZDIV
      SUBROUTINE ZDIV (AR, AI, BR, BI, CR, CI)
C***BEGIN PROLOGUE  ZDIV
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
C            ZBIRY
C***LIBRARY   SLATEC
C***TYPE      ALL (ZDIV-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     REAL*8 COMPLEX DIVIDE C=A/B.
C
C***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
C***ROUTINES CALLED  ZABS
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZDIV
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      REAL*8 AR, AI, BR, BI, CR, CI, BM, CA, CB, CC, CD
      REAL*8 ZABS
      EXTERNAL ZABS
C***FIRST EXECUTABLE STATEMENT  ZDIV
      BM = 1.0_r8/ZABS(BR,BI)
      CC = BR*BM
      CD = BI*BM
      CA = (AR*CC+AI*CD)*BM
      CB = (AI*CC-AR*CD)*BM
      CR = CA
      CI = CB
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
