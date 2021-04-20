*DECK parmsetget
      SUBROUTINE parmsetget (IPAR, IVALUE, ISET)
      LOGICAL ISET
      INTEGER IPAR, IVALUE
C-----------------------------------------------------------------------
C parmsetget saves and recalls one of two error message parameters:
C   LUNIT, the logical unit number to which messages are printed, and
C   MESFLG, the message print flag.
C This is a modification of the SLATEC library routine J4SAVE.
C
C Saved local variables..
C  LUNIT  = Logical unit number for messages.
C           The default is 6 (machine-dependent).
C  MESFLG = Print control flag..
C           1 means print all messages (the default).
C           0 means no printing.
C
C On input..
C   IPAR   = Parameter indicator (1 for LUNIT, 2 for MESFLG).
C   IVALUE = The value to be set for the parameter, if ISET = .true.
c            or the existing value of the parameter, if ISET = .false.
C   ISET   = Logical flag to indicate whether to read or write.
C            If ISET = .true., the parameter will be given
C            the value IVALUE.  If ISET = .false., the parameter
C            will be unchanged, and IVALUE contains the parameter
c            value on exit.
C
C Subroutines/functions called by parmsetget.. None
C-----------------------------------------------------------------------
      INTEGER LUNIT, MESFLG
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this routine.
C-----------------------------------------------------------------------
      SAVE LUNIT, MESFLG
      DATA LUNIT/6/, MESFLG/1/
C
      IF (IPAR .EQ. 1) THEN
        IF (ISET) THEN
           LUNIT = IVALUE
        ELSE
           IVALUE = LUNIT
        ENDIF
      ENDIF
C
      IF (IPAR .EQ. 2) THEN
        IF (ISET) THEN
           MESFLG = IVALUE
        ELSE
           IVALUE = MESFLG
        ENDIF
      ENDIF
C
      RETURN
C----------------------- End of Subroutine parmsetget -------------------------
      END
