*DECK XSETUNP
      SUBROUTINE XSETUNP (LUN)
C-----------------------------------------------------------------------
C This routine resets the logical unit number for messages.
C
C Subroutine called by XSETUNP.. parmsetget
C Function routines called by XSETUNP.. None
C-----------------------------------------------------------------------
      INTEGER LUN
C
      IF (LUN .GT. 0) call parmsetget (1,LUN,.true.)
      RETURN
C----------------------- End of Subroutine XSETUNP ---------------------
      END
