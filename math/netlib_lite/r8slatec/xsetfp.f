*DECK XSETFP
      SUBROUTINE XSETFP (MFLAG)
C-----------------------------------------------------------------------
C This routine resets the print control flag MFLAG.
C
C Subroutine called by XSETFP.. parmsetget
C Function routines called by XSETFP.. None
C-----------------------------------------------------------------------
      INTEGER MFLAG
C
      IF (MFLAG .EQ. 0 .OR. MFLAG .EQ. 1)
     .   call parmsetget (2,MFLAG,.true.)
      RETURN
C----------------------- End of Subroutine XSETFP ----------------------
      END
