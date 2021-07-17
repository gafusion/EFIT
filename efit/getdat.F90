!**********************************************************************
! 
!  PRIMARY CALLING SEQUENCE USED BY REVIEW AND OTHER CODES.
!
!     GETDAT CALLS PTDATA AND RETURNS SHOT DATA IN A MORE USEFUL WAY.
!
!     NSHOT IS THE SHOT NUMBER
!
!     NAME IS THE POINT NAME(10 ASCII CHARACTERS)
!
!     ICAL=0 RETURNS DATA IN DIGITIZER COUNTS
!          1 RETURNS CALIBRATED DATA IN PHYSICAL UNITS
!          2 RETURNS THE VOLTAGE INPUT TO THE DIGITIZER
!          3 RETURNS THE VOLTAGE AT THE INTEGRATOR OUTPUT
!          4 RETURNS THE INTEGRATED SIGNAL IN V-SEC
!          10-19 ARE THE SAME AS 0-9 EXCEPT FOR BASELINE ALGORITHM
!
!     ICAL=0 HAS NO BASELINE SUBTRACTION
!           1-9 USES BASELINE FROM PTDATA (OBTAINED FROM EARLY SAMPLES)
!           10-19 USE DIGITIZER MIDPOINT AS ZERO
!
!     IER IS AN ERROR RETURN FLAG.
!     IER = POSITIVE NUMBER IS THE ERROR RETURN FROM PTDATA.
!     IER = -1 IS FOR OVERFLOW OF THE DATA IN THE DIGITIZER
!     IER = -2 IS FOR UNDERFLOW OF THE DATA
!     IER = -3 IS FOR BASELINE(ZERO) OUTSIDE THE RANGE OF DIGITIZER
!     IER = -6 DATA ARE RETURNED WITH SOME LOSS OF TIME RESOLUTION 
!  (TOO MANY SAMPLES FOR BUFFER SIZE)
!     IER = -7 DATA MAY HAVE SOME LOSS OF TIME RESOLUTION 
!  (TOO MANY SAMPLES FOR BUFFER SIZE)
!  STATUS UNKNOWN DUE TO NECESSITY OF CALLINBG PTDATA BY TIME
!
!     T IS THE ARRAY IN WHICH THE TIMES OF THE DATA SAMPLES WILL BE RETURNED
!
!     DATA IS THE ARRAY IN WHICH THE DATA WILL BE RETURNED
!
!     NP IS THE MAXIMUM NUMBER OF DATA POINTS YOU WANT BACK
!       ACTUAL NUMBER IS RETURNED
!
!     TMIN, TMAX (IN SECONDS) DEFINE THE TIME INTERVAL ON WHICH YOU WANT
!       DATA.  ACTUAL VALUES ARE RETURNED.
!
!     THE POINTNAME DATA IS ALWAYS MULTIPLIED BY SCALE BEFORE USE OF
!
!     MOP--WHICH IS THE OP CODE.
!     IF MOP=0 THE POINT IS WRITTEN INTO DATA, OVERWRITING THE
!        PREVIOUS CONTENTS.
!     IF MOP=1,THE POINT IS ADDED INTO THE ARRAY DATA.
!     IF MOP=2,THE POINT IS SUBTRACTED FROM THE CURRENT CONTENTS OF DATA.
!     IF MOP=3 MULTIPLY CURRENT CONTENTS OF DATA BY POINTNAME
!     IF MOP=4 DIVIDE CURRENT CONTENTS OF DATA BY POINTNAME
!
!     IBFLAG = 0 PTDATA'S "ZERO OFFSET" VALUE IS USED FOR THE BASELINE.
!              1 SAMPLES 6 THROUGH 35 ARE AVERAGED TO PROVIDE THE BASELINE.
!                    IN THIS CASE, THE FIRST 40 SAMPLES ARE NOT RETURNED.
!----- 
! 1% tolerance for signal clipping added 8/14/89
! at present the clipping errors (ier = -1, -2) are used only by the
! yoka routine in NEWSPEC.
!-----
!
      SUBROUTINE GETDAT(NSHOT,NAME,ICAL,IER,T,DATA,NP, &
      TMIN,TMAX,MOP,SCALE,jWAIT)

      use set_kinds
      parameter (npmax=262144)

      character*4 char4

! 20140905 tbt Dimension of 1 bombs when array bounds checked!
!     DIMENSION  DATA(1)
!     DIMENSION  FPDATA(1)

      DIMENSION  DATA(*)
      DIMENSION  FPDATA(npmax)

      dimension  int16(2),int32(2),iarr(50)
      dimension  iascii(12)
      CHARACTER*3  MONTH(12)
      CHARACTER  NAME*(*)
      character*12  pcst
      CHARACTER*10 PCSTIME
      CHARACTER  PHASE*4
      dimension  rarr(20+npmax)
      dimension  real32(100)
      REAL*8  REAL64(100)
      CHARACTER*10  SDATE, STIME

! 20140905 tbt Dimension of 1 bombs when array bounds checked!
!     DIMENSION  T(1)
      DIMENSION  T(*)

      REAL*8  TDUM(2+NPMAX)
      INTEGER*4  TEMP(npmax)
      REAL*8  TIME64(NPMAX)

! 20140905 tbt Dimension of 1 bombs when array bounds checked!
!     DIMENSION  TT(1)
      DIMENSION  TT(npmax)

      COMMON /TRAPPER/ IT_CHAN, ITRAPA

      equivalence (ichar,char4)
      equivalence (pcst,pcstime)
      EQUIVALENCE (TT,RARR(21)), (FPDATA(1),TEMP(1))

      DATA MONTH/'JAN','FEB','MAR','APR','MAY','JUN', &
         'JUL','AUG','SEP','OCT','NOV','DEC' /
      DATA PHASE /'.PLA'/
      data tolclip /.01/ !-- maximum fraction of clipped signals
      data kmax/10/
      data wait_time/30./

      LOGICAL INTCAL,efit
      character*10 cname, cname1, filnam
!----------------------------------------------------------------------

      INTCAL = .FALSE.
      efit = .false.
      rcfact = 1.0
      iwait = jwait

      GO TO 11111

!---------------------------------------------------------------
!
!  ALTERNATE ENTRY POINT FOR USE BY INTEGRATOR CALIBRATION CODE INTCAL.
!  THE ONLY DIFFERENCE IS THE RETURN OF ADDITIONAL PARAMETERS.
!
!---------------------------------------------------------------
      ENTRY GETDAT_I(NSHOT,NAME,ICAL,IER,T,DATA,NP, &
        TMIN,TMAX,MOP,SCALE,jWAIT, &
        RC_I,RCG_I,VPBIT_I,MZERO_I,NBITS_I,SDATE,STIME)

      INTCAL = .TRUE.
      efit = .false.
      rcfact = 1.0
      iwait = jwait

      GO TO 11111

!----------------------------------------------------------------------
!
! Alternate entry point for use by the EFIT code.
!
!----------------------------------------------------------------------

      ENTRY GETDAT_E(NSHOT,NAME,ICAL,IER,T,DATA,NP, &
      TMIN,TMAX,MOP,SCALE,bitone,ircfact,RC_E,RCG_E,VPBIT_E, &
       ZINHNO_E,T0_E)

      INTCAL = .FALSE.
      efit = .true.
      iwait = 0

        rcfact = 1.0

        IF (ircfact .eq. 1) THEN
                filnam='rcfact.dat'
                cname=to_upper(name)
                open(unit=21,file=filnam,status='old',err=1000)
1001            read (21,9,end=1000,err=1000) cname1, fact
9               format(x,a10,f15.5)
                if(cname .ne. cname1)   go to 1001
                if(fact .gt. 0.) rcfact = fact
1000            close(unit=21)
        ENDIF

!----------------------------------------------------------------------
11111 kount=0

1      IASCII(1)=5
      INT16(1)=0
      INT32(1)=0
      REAL32(1)=50

! ---------------------- SET UP FOR THE PTDATA CALL ----------------------

      itype = 12  ! point type call, variable timing
      iarr(1) = npmax  ! number of points requested
      iarr(3) = 1  ! start point
      nint = 1
      iarr(4) = nint  ! point interval
      pzero = 0.0
!
!  PTDATA CALL ... GET ALL THE DATA!!
!  Request to get data by point sampling with PTDATA returning time array
!
!vas f90 modifi
!10 call ptdata(itype,nshot,%ref(phase),%ref(name),temp,ier,iarr,rarr,
!vas 1 iascii,int16,int32,real32)

10    call ptdata(itype,nshot,%ref(phase),%ref(name),temp,ier,iarr,rarr,&
        iascii,int16,int32,real32)

! write(13,91300)(rarr(klg),klg=540,640)
91300 format(x,1pe17.10)

      if (ier .eq. 4 .or. ier.eq.2) ier = 0
      VPBIT  = real32(51)
      RC     = real32(52)

!
!  CHECK DATA FORMAT, TRY AGAIN IF NECESSARY
!  Attempt to get data based on time sampling, using a computed
!  start and delta time.
!

      IF (ier .ne. 35) go to 30 ! ier=35 ... wrong call type for dfi

      itype = 2
!vas f90 modifi
!vas20 call ptdata(itype,nshot,%ref(phase),%ref(name),temp,ier,iarr,rarr,
!vas 1 iascii,int16,int32,real32)
20    call ptdata(itype,nshot,%ref(phase),%ref(name),temp,ier,iarr,rarr,&
        iascii,int16,int32,real32)
      if (ier .eq. 4 .or. ier.eq.2) ier = 0

      VPBIT  = real32(13)
      RC     = 1.0

30    continue

!
!  WAIT FOR THE DATA TO COME IN
!

!vas IF (iwait .eq. 1 .and. (ier .eq. 3 .or. ier .eq. 1)
!vas 1   .and. kount .lt. kmax .and. itrapa .eq. 0) THEN
      IF (iwait .eq. 1 .and. (ier .eq. 3 .or. ier .eq. 1) &
         .and. kount .lt. kmax .and. itrapa .eq. 0) THEN
! MPI >>>
        !call rev_wait(iwait, nshot)
! MPI <<<
        if (itrapa .eq. 0) go to 1
      ENDIF

      !-- AT THIS POINT ANY REMAINING POSITIVE ERROR VALUE IS FATAL.
      !-- NON-FATAL POSITIVE ERRORS HAVE BEEN HANDLED AS FOLLOWS:
      !  ier = 4  reset ier to 0
      !  ier = 2  non-fatal warning:  pointname has been revised
      !  ier = 1 or 3 looped back to start if iwait=1
      !  ier = 35  re-called ptdata w/o time array for this dfi

      IF (IER .gt. 0 .AND. IER.NE.2) THEN
70000   np = 2
        t(1)=tmin
        t(2)=tmax
        data(1)=0.
        data(2)=0.
        RETURN
      ENDIF

      !
      !  At this point, data is available.  Save some required information.
      !

      idfi = iarr(31)
      nret = iarr(2)
      nsampl = iarr(32)   ! number of 16-bit words
      nbits = iarr(33)   ! number of bits per word
      ZINHNO = rarr(4)                        ! inherent number
      RCG    = rarr(5)
      if(nbits .le.  8) nsampl = nsampl*2   ! correct number of
      if(nbits .gt. 16) nsampl = nsampl/2 !       available samples
      mzero = iarr(30)
      yzero = mzero
      IF (idfi.eq.158 .OR. IDFI.EQ.2158) THEN !-- NEW DFI WITH MORE ACCURATE OFFSET
        yzero = real32(14)
        mzero = anint(yzero)
      ENDIF

!--- WAVEFORMS HAVE AN ADDITIONAL OFFSET FOR CONVERSION TO PHYSICAL UNITS
! --- this offset makes sense only with ical=1 or 11,
! --- by convention we will add offset in ical=1 but not 11
!vas f90 modifi
!vas if ((idfi .eq. 125 .or. idfi .eq. 158 .OR. IDFI.EQ.2158)
!vas 1 .and. ical .eq. 1)  pzero = real32(8)
      if ((idfi .eq. 125 .or. idfi .eq. 158 .OR. IDFI.EQ.2158) &
        .and. ical .eq. 1)  pzero = real32(8)

        !--- Rely on PTDATA's time-type call if the
        !    pointname is too big for the buffer.
        !--- This approach may return less than optimum resolution if the requested
        !    interval crosses a clock domain boundary.
        !
        !  Must handle data from new (4-95) Plasma Control System separately.
        !  Assume that PCS data will be less than max qty asked for already.
        !  (Apr 5, 1994:  John Ferron said some waveforms could have as much
        !  as 90K samples.  Max data requested in 1st PTDATA call 128K.)
        !
        !  If data is NOT PCS, then compute appropriate start and delta times
        !  in order to do 2nd PTDATA call.

        if (idfi.eq.2201 .or. idfi.eq.2202 .or. idfi.eq.2203) then
          if (idfi.eq.2201) then  ! digitizer data, int
            vpbit = real32(5)
            yzero = real32(6)
            if (real32(2).ge.5) rc = real32(7)
          else if (idfi.eq.2202) then  ! processed data, int
            vpbit = -1.0
            pzero = real32(3)
            yzero = 0.0
          else if (idfi.eq.2203) then  ! processed data, real
            vpbit = -1.0
            pzero = real32(3)
            yzero = 0.0
        end if
        ichar = iascii(3)
        pcst(1:4) = char4
        ichar = iascii(4)
        pcst(5:8) = char4
        ichar = iascii(5)
        pcst(9:12) = char4
        iarr(1) = npmax
        iarr(3) = 1
        iarr(4) = nint
        iascii(1) = 0
        int16(1) = 0
        int32(1) = 0
        real32(1) = 0
        real64(1) = 0
        call ptdata64(64,nshot,%ref(phase),%ref(pcstime),time64,ker, &
          iarr,rarr,iascii,int16,int32,real32,real64,tdum)
        if (ker.ne.0 .and. ker.ne.2 .and. ker.ne.4) then
          ier = 1
          go to 70000
        end if
        nrettim = iarr(2)
        do k = 1,nrettim
          tt(k) = time64(k)
        end do
        go to 85000
      else IF (itype .eq. 12 .or. itype .eq. 2) THEN
      !
      !  at this point, is there more digitizer data available?
      !
      IF (nret .lt. nsampl) THEN
        itype = itype - 1  ! call ptdata by time
        rarr(1) = tmin   ! start time
        rarr(2) = (tmax-tmin)/(npmax-1) ! minimum delta-time
        iarr(1) = npmax   ! guaranteed to cover tmax-tmin
        IF (itype .eq. 11) THEN
          go to 10
        ELSE
          go to 20
        ENDIF    ! itype = 11
      endif     ! nret < nsample
      endif     ! itype = 12

      !--- ier = -7 for time type call ... possible loss of time resolution

      if ((itype .eq. 1 .or. itype .eq. 11)  .and. ier .le. 0) ier = -7

      !  FILL IN TIME ARRAY IF NECESSARY

      IF (itype .eq. 2) THEN
        dt = rarr(9)
        t0 = rarr(8) - dt
        DO i=1, nret
          tt(i) = t0 + i * dt
        ENDDO
      ENDIF

85000 continue

      !  FIND START AND STOP TIME OF DATA TO BE RETURNED
      !  --- THE OVERLAP OF THE REQUESTED AND DIGITIZED TIME INTERVALS

      tmind = tt(1)
      tmaxd = tt(nret)

      IF (tmax .le. tmin) THEN
        tmin = tmind
        tmax = tmaxd
      ENDIF

      tmin = amax1(tmin,tmind)
      tmax = amin1(tmax,tmaxd)

      IF (tmax .lt. tmin) THEN
        np = 0
        return
      ENDIF

      !  FIND THE FIRST AND LAST DATA POINTS TO BE RETURNED
      ! Traditionally, the efit code has used a different limit test here.  So that
      ! efit will continue to get the same answer, one of two different limit tests
      ! is used here.
      !
      nfirst = 1
      nlast  = 0
      ITOTALPTS = NRET
      IF (NSAMPL.LT.ITOTALPTS) ITOTALPTS = NSAMPL

      if (efit) then
        DO i=1,ITOTALPTS
          ttime = tt(i)
          if (ttime .lt. tmin) nfirst=i+1
          if (ttime .lt. tmax) nlast =i
        ENDDO
      else
        DO i=1,ITOTALPTS
          ttime = tt(i)
          if (ttime .lt. tmin) nfirst=i+1
          if (ttime .le. tmax) nlast =i
        ENDDO
      endif

      ndata = nlast - nfirst + 1

      !  DECIDE ON DATA POINT INTERVAL

      if (np .eq. 0) then
        np = 16384
      else
        np = amin0(np,npmax)
      endif

      nint = (ndata-1) / np + 1
      ! if (nint .gt. 1 .and. ier .le. 0) ier = -6
      if (nint.gt.1) ier = -6

      ! SOME OF THE ICH DATA  (ICHPWR) IS
      ! JUST LIKE DIGITIZER DATA EXCEPT IN
      ! REAL FORMAT

      IF(IDFI .EQ. 119 .or. idfi.eq.2119) THEN      !ICH STUFF EXISTS IN FLOATING FORMAT
        II=0
        DO I = NFIRST, NLAST, NINT
          II=II+1
          T(II) = TT(I)
          DATA(II) = FPDATA(I)
        ENDDO
        NP = II
        RETURN
      ENDIF

      !  CHECK FOR ZERO OFFSET OUT OF RANGE
      !  2013/04/03 Revised for new digitzers signed binary numbers  LL/ES

      ! min = 0
      !
      !  THIS NEXT LINE HAS BEEN MODIFIED TO (NBITS - 2) FOR TWO REASONS:
      !  (1)  THE FIRST BIT IS ACTUALLY BIT #0 (THEREFORE FOR A 32-BIT VARIABLE
      !       THE LARGEST UNSIGNED VALUE THAT THE VARIABLE CAN CONTAIN IS
      !       2**(NBITS) - 1
      !  (2)  FOR A "SIGNED" DATA VALUE, THE LAST BIT IS USED TO INDICATE THE
      !       SIGN.  THEREFORE, THE LARGEST VALUE A SIGNED VARIABLE CAN CONTAIN
      !       IS     2**(NBITS - 1) - 1
      !  THE CALCULATION IN #2 ABOVE WILL RESULT IN AN INTEGER OVERFLOW BECAUSE
      !  THE MULTIPLICATION IS DONE BEFORE THE SUBTRACTION.  AND IT IS THE
      !  RESULT OF THE MULTIPLICATION THAT WILL EXCEED THE ABILITY OF THE
      !  INTEGER TO HOLD IT.
      !
      !  THE NEW LINE IS A MODIFIED FORM TO MAKE SURE THAT INTEGER OVERFLOW
      !  DOES NOT OCCUR
      !
      IF (NBITS.GE.16) THEN  !  NECESSARY ONLY FOR 16-BIT DATA
        min =-(2 ** (nbits - 2) - 1) * 2 + 1
        max = (2 ** (nbits - 2) - 1) * 2 + 1
      ELSE    !  12-BIT DIG.S NEED THIS ONE
        min = 0
        MAX = 2**NBITS - 1
      END IF


      !--- WAVEFORMS ARE ALLOWED TO HAVE ZERO OFFSET OUTSIDE THE GENERATOR RANGE
      !vas f90 modifi
      !vas IF ((mzero .lt. min .or. mzero .gt. max) .AND.
      !vas 1 (idfi .ne. 125 .and. idfi .ne. 158 .AND. &
      IF ((mzero .lt. min .or. mzero .gt. max) .AND. &
        (idfi .ne. 125 .and. idfi .ne. 158 .AND. &
        IDFI.NE.2158)) THEN
        IF (ier .eq. 0 ) ier = -3
      ENDIF

      !  CALIBRATION CODES TO RETURN ABSOLUTE SIGNAL LEVEL, NO BASELINE SUBTRACTION
      !  MUST ASSUME DIGITIZER ZERO IS IN THE MIDDLE OF ITS RANGE, SINCE THIS
      !  INFORMATION IS NOT STORED IN THE DATA BASE!
      !  WAVEFORMS (idfi = 125, 158, 2158) DO NOT HAVE A BASELINE, SINCE WE ARCHIVE
      !  ONLY THE PROGRAMMED VALUE AND NOT THE ACTUAL VALUE.
      ical0 = ical
      IF (ical .ge. 10) THEN
        if (idfi .ne. 125 .and. idfi .ne. 158 &
          .AND. IDFI.NE.2158) mzero = 2**(nbits-1)
        yzero = mzero
        ical0 = ical-10
      ENDIF

      !   FILL UP PLOT ARRAYS

      nclip = 0
      ii = 0

      DO i=nfirst, nlast, nint

        ii = ii + 1
        !
        ! --- time array
        !
        t(ii) = tt(i)
        ! write(14,91300)t(ii)
        !
        ! --- data calibration
        !     some waveforms for PCS have data that is REAL*4
        !
        if (idfi.eq.2203) then
          y = fpdata(i)
        else
          y = temp(i)
        endif

        ! --- check for data out of digitizer range
        !      IF (TEMP(I) .GE. MAX .AND. IER.EQ.0)  IER = -1
        !      IF (TEMP(I) .LE. MIN .AND. IER.EQ.0)  IER = -2
        if(y .ge. max .or. y .le. min) nclip = nclip + 1

        ! --- offset for all but calibration code 0
        IF (ical .ne. 0) THEN
          y = y - yzero
        !  y = y - yzero - deltab   ! ICAL > 0
        ENDIF

        ! --  scale factors for the different calibration codes

        IF (ical0 .eq. 0) THEN
          y = y    ! ICAL = 0
        ELSE IF (ical0 .eq. 2) THEN
          y = -y * vpbit   ! ICAL = 2
        ELSE IF (ical0 .eq. 3) THEN
          y = -y * rcg / rc  ! ICAL = 3
        ELSE IF (ical0 .eq. 4) THEN
          y = -y * rcg   ! ICAL = 4
        ELSE
          y = -y * rcg * zinhno + pzero ! ICAL = 1 or other
        ENDIF

        y = scale * y * rcfact    ! user's scale factor

        ! --- arithmetic combinations as determined by mop code
        IF (mop .eq. 0) THEN
          data(ii) = y
        ELSE IF (mop .eq. 1) THEN
          data(ii) = data(ii) + y
        ELSE IF (mop .eq. 2) THEN
          data(ii) = data(ii) - y
        ELSE IF (mop .eq. 3) THEN
          data(ii) = data(ii) * y
        ELSE IF (mop .eq. 4) THEN
          if (y.ne.0) then
            data(ii) = data(ii) / y
          else
            data(ii) = 0.
          end if
        ELSE
          data(ii) = y
        ENDIF

      ENDDO



      ! --- actual number of points returned
      np = ii

      ! --- clipping tolerance:
      !     set error flag only if clipped fraction is outside tolerance
      if(np .gt. 0) then
        fclip = real(nclip,dp)/np
        if (fclip .gt. tolclip .and.  &
          ier .le. 0) ier = -1
      endif
      IF (INTCAL) then
        !
        !   ADDITIONAL ARGUMENTS FOR INTEGRATOR CALIBRATION CALL.
        ! THESE DUMMY ARGUMENTS MUST ONLY BE REFERENCED WHEN USING THE
        ! ENTRY POINT GETDAT_I.  OTHERWISE ACCESS VIOLATION ERRORS OCCUR.
        !
        RC_I = RC
        RCG_I = RCG
        VPBIT_I = VPBIT
        MZERO_I = MZERO
        NBITS_I = NBITS
        WRITE(SDATE,880) IARR(19),MONTH(IARR(18)), &
          IARR(20)-IARR(20)/100*100
        WRITE(STIME,881) IARR(15),IARR(16),IARR(17)
880     FORMAT(I2,'-',A3,'-',I2)
881     FORMAT(I2,':',I2,':',I2)
      !
      endif

      if(efit) then
        if (ical0 .eq. 0) then
          bitone = 1                              ! ical = 0
        else if (ical0 .eq. 2) then
          bitone = vpbit                          ! ical = 2
        else if (ical0 .eq. 3) then
          bitone = rcg / rc                       ! ical = 3
        else
          bitone = rcg * zinhno                   ! ical = 1 or other
        endif

        bitone = scale * abs(bitone)
        !--------------------------------------------------------------------
        !-- New returned values for new error matrix                       --
        !--------------------------------------------------------------------
        RC_E = RC
        RCG_E = RCG
        VPBIT_E = VPBIT
        ZINHNO_E = ZINHNO
        T0_E = TT(1)
      endif

      RETURN

      CONTAINS
          Function to_upper (str)

          !   ==============================
          !   Changes a string to upper case
          !   ==============================

              Implicit None
              Character(*), Intent(In) :: str
              Character(LEN(str))      :: to_upper

              Integer :: ic, i

              Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
              Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

          !   Capitalize each letter if it is lowecase
              to_upper = str
              do i = 1, LEN_TRIM(str)
                  ic = INDEX(low, str(i:i))
                  if (ic > 0) to_upper(i:i) = cap(ic:ic)
              end do

           End Function to_upper
      END SUBROUTINE

