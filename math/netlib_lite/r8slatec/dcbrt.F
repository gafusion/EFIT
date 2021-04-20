*DECK DCBRT
      REAL*8 FUNCTION DCBRT (X)
C***BEGIN PROLOGUE  DCBRT
C***PURPOSE  Compute the cube root.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C2
C***TYPE      REAL*8 (CBRT-S, DCBRT-D, CCBRT-C)
C***KEYWORDS  CUBE ROOT, ELEMENTARY FUNCTIONS, FNLIB, ROOTS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DCBRT(X) calculates the double precision cube root for
C double precision argument X.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, D9PAK, D9UPAK
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DCBRT
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER niter,n,ixpnt,irem,iter
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 z
C============
      REAL*8 X, CBRT2(5), Y, CBRTSQ,  D9PAK, D1MACH
      SAVE CBRT2, NITER
      DATA CBRT2(1) / 0.62996052494743658238360530363911_r8/
      DATA CBRT2(2) / 0.79370052598409973737585281963615_r8/
      DATA CBRT2(3) / 1.0_r8/
      DATA CBRT2(4) / 1.25992104989487316476721060727823_r8/
      DATA CBRT2(5) / 1.58740105196819947475170563927231_r8/
      DATA NITER / 0 /
C***FIRST EXECUTABLE STATEMENT  DCBRT
      IF (NITER.EQ.0) NITER =
     1   1.443_r8*
     1   LOG(-.106_r8*LOG(0.1_r8*real(D1MACH(3),r8)))
     1   + 1.0_r8
C
      DCBRT = 0._r8
      IF (X.EQ.0._r8) RETURN
C
      CALL D9UPAK (ABS(X), Y, N)
      IXPNT = N/3
      IREM = N - 3*IXPNT + 3
C
C THE APPROXIMATION BELOW IS A GENERALIZED CHEBYSHEV SERIES CONVERTED
C TO POLYNOMIAL FORM.  THE APPROX IS NEARLY BEST IN THE SENSE OF
C RELATIVE ERROR WITH 4.085 DIGITS ACCURACY.
C
      Z = Y
      DCBRT = .439581_r8+ Z*(.928549_r8+
     >   Z*(-.512653_r8+ Z*.144586_r8))
C
      DO 10 ITER=1,NITER
        CBRTSQ = DCBRT*DCBRT
        DCBRT = DCBRT + (Y-DCBRT*CBRTSQ)/(3._r8*CBRTSQ)
 10   CONTINUE
C
      DCBRT = D9PAK (CBRT2(IREM)*SIGN(DCBRT,X), IXPNT)
      RETURN
C
      END
c------------------------------------------
      real*8 function r8_dcbrt(x)
      real*8 x,dcbrt
      external dcbrt
 
      r8_dcbrt=dcbrt(x)
 
      return
      end
 
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
! 22Jun2000 fgtok
! 22Jun2000 fgtok
