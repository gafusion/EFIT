c-----------------------------------------------------------------------
c     file lsode_module.f
c     This contains the explicit interface for lsode.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0.  lsode_mod.
c-----------------------------------------------------------------------
c     subprogram 0. lsode_mod
c-----------------------------------------------------------------------
      MODULE lsode_mod
      IMPLICIT NONE

      INTERFACE exch_lagr
      SUBROUTINE LSODE (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1                  ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
        INTEGER, OPTIONAL :: MF
        INTEGER, OPTIONAL :: JAC
        EXTERNAL F
        !EXTERNAL JAC
        INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW
        REAL*8 Y, T, TOUT, RTOL, ATOL, RWORK
        DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW)
      END  SUBROUTINE 
      END INTERFACE  
c-----------------------------------------------------------------------
      END MODULE lsode_mod
