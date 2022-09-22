#include "config.f"
!**********************************************************************
!>
!!    PRIMARY CALLING SEQUENCE USED BY REVIEW AND OTHER CODES.\n
!!
!!    1% tolerance for signal clipping added 8/14/89
!!    at present the clipping errors (ier = -1, -2) are used only by the
!!    yoka routine in NEWSPEC.
!!
!!    WARNING: this subroutine uses both REAL*4 and REAL*8 variables 
!!             conversions must be handled carefully
!!             REAL*4 variables are both used for PTDATA calls and as
!!             inputs for non-EFIT entries to conform with use by other
!!             libraries (e.g. mselib)
!!
!!
!!    @param NSHOT :  THE SHOT NUMBER
!!
!!    @param NAME : THE POINT NAME(10 ASCII CHARACTERS)
!!
!!    @param ICAL : ICAL=0 RETURNS DATA IN DIGITIZER COUNTS\n
!!                       1 RETURNS CALIBRATED DATA IN PHYSICAL UNITS\n
!!                       2 RETURNS THE VOLTAGE INPUT TO THE DIGITIZER\n
!!                       3 RETURNS THE VOLTAGE AT THE INTEGRATOR OUTPUT\n
!!                       4 RETURNS THE INTEGRATED SIGNAL IN V-SEC\n
!!                       10-19 ARE THE SAME AS 0-9 EXCEPT FOR BASELINE ALGORITHM\n
!!                  ICAL=0 HAS NO BASELINE SUBTRACTION\n
!!                       1-9 USES BASELINE FROM PTDATA (OBTAINED FROM EARLY SAMPLES)\n
!!                       10-19 USE DIGITIZER MIDPOINT AS ZERO\n
!!
!!    @param IER : POSITIVE NUMBER IS THE ERROR RETURN FROM PTDATA\n
!!           IER = -1 IS FOR OVERFLOW OF THE DATA IN THE DIGITIZER\n
!!           IER = -2 IS FOR UNDERFLOW OF THE DATA\n
!!           IER = -3 IS FOR BASELINE(ZERO) OUTSIDE THE RANGE OF DIGITIZER\n
!!           IER = -6 DATA ARE RETURNED WITH SOME LOSS OF TIME RESOLUTION
!!                    (TOO MANY SAMPLES FOR BUFFER SIZE)\n
!!           IER = -7 DATA MAY HAVE SOME LOSS OF TIME RESOLUTION\n
!!                    (TOO MANY SAMPLES FOR BUFFER SIZE)\n
!!                    STATUS UNKNOWN DUE TO NECESSITY OF CALLING PTDATA BY TIME\n
!!
!!    @param T_R4 :  THE ARRAY IN WHICH THE TIMES OF THE DATA SAMPLES WILL BE RETURNED
!!
!!    @param DATA_R4 : IS THE ARRAY IN WHICH THE DATA WILL BE RETURNED
!!
!!    @param NP : THE MAXIMUM NUMBER OF DATA POINTS YOU WANT BACK
!!                ACTUAL NUMBER IS RETURNED
!!
!!    @param &TMIN_R4 : DEFINE THE TIME INTERVAL ON WHICH YOU WANT
!!                   DATA.  ACTUAL VALUES ARE RETURNED.
!!
!!    @param TMAX_R4 : DEFINE THE TIME INTERVAL ON WHICH YOU WANT
!!                  DATA.  ACTUAL VALUES ARE RETURNED.
!!
!!    @param MOP : WHICH IS THE OP CODE.\n
!!           MOP=0 THE POINT IS WRITTEN INTO DATA, OVERWRITING THE
!!                 PREVIOUS CONTENTS.\n
!!           MOP=1,THE POINT IS ADDED INTO THE ARRAY DATA.\n
!!           MOP=2,THE POINT IS SUBTRACTED FROM THE CURRENT CONTENTS OF DATA.\n
!!           MOP=3 MULTIPLY CURRENT CONTENTS OF DATA BY POINTNAME\n
!!           MOP=4 DIVIDE CURRENT CONTENTS OF DATA BY POINTNAME\n
!!
!!    @param SCALE_R4 : FACTOR TO SCALE RESULTS BY
!!
!!    @param jWAIT : UNUSED VARIABLE
!!
!**********************************************************************
      SUBROUTINE GETDAT(NSHOT,NAME,ICAL,IER,T_R4,DATA_R4,NP, &
                        TMIN_R4,TMAX_R4,MOP,SCALE_R4,jWAIT)

      use set_kinds
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      parameter (npmax=262144)

      character*4 char4
      real*4 t_r4, data_r4, tmin_r4, tmax_r4, scale_r4, rarr4, real32

      DIMENSION DATA_R4(*)
      DIMENSION T_R4(*)
      DIMENSION DATA(npmax)
      DIMENSION T(npmax)
      DIMENSION FPDATA(npmax)
      DIMENSION TT(npmax)

      dimension int16(2),int32(2),iarr(50)
      dimension iascii(12)
      CHARACTER*3 MONTH(12)
      CHARACTER NAME*(*)
      character*12 pcst
      CHARACTER*10 PCSTIME
      CHARACTER PHASE*4
      real*8 rarr, realdp, tt, dt, t0
      dimension rarr(20+npmax),rarr4(20+npmax)
      dimension real32(100),realdp(100)
      REAL*8 REAL64(100)
      CHARACTER*10 SDATE, STIME

      REAL*8 TDUM(2+NPMAX)
      INTEGER*4 TEMP(npmax)
      REAL*8 TIME64(NPMAX)

      equivalence (ichar,char4)
      equivalence (pcst,pcstime)
      equivalence (tt,rarr(21)), (fpdata(1),temp(1))

      DATA MONTH/'JAN','FEB','MAR','APR','MAY','JUN', &
         'JUL','AUG','SEP','OCT','NOV','DEC' /
      DATA PHASE /'.PLA'/
      data tolclip /.01/ !-- maximum fraction of clipped signals
      data kmax/10/

      LOGICAL INTCAL,efit
      character*10 cname, cname1, filnam
!----------------------------------------------------------------------

      INTCAL = .FALSE.
      efit = .false.
      rcfact = 1.0

      tmin = real(tmin_r4,dp)
      tmax = real(tmax_r4,dp)
      scale = real(scale_r4,dp)

      GO TO 11111

!---------------------------------------------------------------
!
!  ALTERNATE ENTRY POINT FOR USE BY INTEGRATOR CALIBRATION CODE INTCAL.
!  THE ONLY DIFFERENCE IS THE RETURN OF ADDITIONAL PARAMETERS.
!
!---------------------------------------------------------------
!      ENTRY GETDAT_I(NSHOT,NAME,ICAL,IER,T,DATA,NP, &
!                     TMIN,TMAX,MOP,SCALE,jWAIT, &
!                     RC_I,RCG_I,VPBIT_I,MZERO_I,NBITS_I,SDATE,STIME)

!      INTCAL = .TRUE.
!      efit = .false.
!      rcfact = 1.0

!      GO TO 11111

!----------------------------------------------------------------------
!
!     Alternate entry point for use by the EFIT code.
!
!----------------------------------------------------------------------

      ENTRY GETDAT_E(NSHOT,NAME,ICAL,IER,T,DATA,NP, &
                     TMIN_R8,TMAX_R8,MOP,SCALE_R8,bitone,ircfact,RC_E,RCG_E, &
                     VPBIT_E,ZINHNO_E,T0_E)

      INTCAL = .FALSE.
      efit = .true.

      rcfact = 1.0

      tmin = tmin_r8
      tmax = tmax_r8
      scale = scale_r8

      IF (ircfact .eq. 1) THEN
        filnam='rcfact.dat'
        cname=to_upper(name)
        open(unit=21,file=filnam,status='old',iostat=ioerr)
        if (ioerr.ne.0) then
1001      read (21,9,iostat=ioerr) cname1, fact
          if (ioerr.ne.0) then
9           format(x,a10,f15.5)
            if(cname .ne. cname1) go to 1001
            if(fact .gt. 0.) rcfact = fact
          else
            close(unit=21)
          endif
        else
          close(unit=21)
        endif
      ENDIF

!----------------------------------------------------------------------
11111 kount=0

1     IASCII(1)=5
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
10    continue
      rarr4=real(rarr,r4)
#ifdef __cray
      call ptdata(itype,nshot,phase,name,temp,ier,iarr,rarr4,&
                  iascii,int16,int32,real32)
#else
      call ptdata(itype,nshot,%ref(phase),%ref(name),temp,ier,iarr,&
                  rarr4,iascii,int16,int32,real32)
#endif
      rarr=real(rarr4,dp)
      realdp=real(real32,dp)

! write(13,91300)(rarr(klg),klg=540,640)
91300 format(x,1pe17.10)

      if(ier .eq. 4 .or. ier.eq.2) ier = 0
      VPBIT = realdp(51)
      RC    = realdp(52)
      !
      !  CHECK DATA FORMAT, TRY AGAIN IF NECESSARY
      !  Attempt to get data based on time sampling, using a computed
      !  start and delta time.
      !
      if(ier .ne. 35) go to 666 ! ier=35 ... wrong call type for dfi

        itype = 2
20      continue
        rarr4=real(rarr,r4)
#ifdef __cray
        call ptdata(itype,nshot,phase,name,temp,ier,iarr,rarr4,&
                    iascii,int16,int32,real32)
#else
        call ptdata(itype,nshot,%ref(phase),%ref(name),temp,ier,iarr,rarr4,&
                    iascii,int16,int32,real32)
#endif
        rarr=real(rarr4,dp)
        realdp=real(real32,dp)
        if(ier .eq. 4 .or. ier.eq.2) ier = 0

        VPBIT = realdp(13)
        RC    = 1.0

666   continue

      !-- AT THIS POINT ANY REMAINING POSITIVE ERROR VALUE IS FATAL.
      !-- NON-FATAL POSITIVE ERRORS HAVE BEEN HANDLED AS FOLLOWS:
      !  ier = 4  reset ier to 0
      !  ier = 2  non-fatal warning:  pointname has been revised
      !  ier = 35  re-called ptdata w/o time array for this dfi

      IF (IER .gt. 0 .AND. IER.NE.2) THEN
        np = 2
        if (efit) then
          t(1)=tmin
          t(2)=tmax
          data(1)=0.
          data(2)=0.
        else
          t_r4(1)=tmin_r4
          t_r4(2)=tmax_r4
          data_r4(1)=0.
          data_r4(2)=0.
        endif
        RETURN
      ENDIF
      !
      !  At this point, data is available.  Save some required information.
      !
      idfi = iarr(31)
      nret = iarr(2)
      nsampl = iarr(32)   ! number of 16-bit words
      nbits = iarr(33)   ! number of bits per word
      ZINHNO = rarr(4)             ! inherent number
      RCG    = rarr(5)
      if(nbits .le.  8) nsampl = nsampl*2   ! correct number of
      if(nbits .gt. 16) nsampl = nsampl/2 !       available samples
      mzero = iarr(30)
      yzero = mzero
      if (idfi.eq.158 .or. idfi.eq.2158) then !-- NEW DFI WITH MORE ACCURATE OFFSET
        yzero = realdp(14)
        mzero = anint(yzero)
      endif

!--- WAVEFORMS HAVE AN ADDITIONAL OFFSET FOR CONVERSION TO PHYSICAL UNITS
! --- this offset makes sense only with ical=1 or 11,
! --- by convention we will add offset in ical=1 but not 11
      if((idfi .eq. 125 .or. idfi .eq. 158 .OR. idfi.eq.2158) &
        .and. ical .eq. 1)  pzero = realdp(8)

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

      pt_again: if (idfi.eq.2201 .or. idfi.eq.2202 .or. idfi.eq.2203) then
        select case (idfi)
        case (2201) ! digitizer data, int
          vpbit = realdp(5)
          yzero = realdp(6)
          if (realdp(2).ge.5) rc = realdp(7)
        case (2202)  ! processed data, int
          vpbit = -1.0
          pzero = realdp(3)
          yzero = 0.0
        case (2203) ! processed data, real
          vpbit = -1.0
          pzero = realdp(3)
          yzero = 0.0
        end select
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
        rarr4=real(rarr,r4)
#ifdef __cray
        call ptdata64(64,nshot,phase,pcstime,time64,ker, &
                      iarr,rarr4,iascii,int16,int32, &
                      real32,real64,tdum)
#else
        call ptdata64(64,nshot,%ref(phase),%ref(pcstime),time64,ker, &
                      iarr,rarr4,iascii,int16,int32, &
                      real32,real64,tdum)
#endif
        rarr=real(rarr4,dp)
        if (ker.ne.0 .and. ker.ne.2 .and. ker.ne.4) then
          ier = 1
          np = 2
          if (efit) then
            t(1)=tmin
            t(2)=tmax
            data(1)=0.
            data(2)=0.
          else
            t_r4(1)=tmin_r4
            t_r4(2)=tmax_r4
            data_r4(1)=0.
            data_r4(2)=0.
          endif
          RETURN
        endif
        nrettim = iarr(2)
        do k = 1,nrettim
          tt(k) = time64(k)
        enddo
      elseif (itype .eq. 12 .or. itype .eq. 2) then pt_again
        !
        !  at this point, is there more digitizer data available?
        !
        if (nret .lt. nsampl) then
          itype = itype - 1  ! call ptdata by time
          rarr(1) = tmin   ! start time
          rarr(2) = (tmax-tmin)/(npmax-1) ! minimum delta-time
          iarr(1) = npmax   ! guaranteed to cover tmax-tmin
          IF (itype .eq. 11) THEN
            go to 10
          ELSE
            go to 20
          ENDIF
        endif
      endif pt_again

      if (idfi.ne.2201 .and. idfi.ne.2202 .and. idfi.ne.2203) then

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

      endif

      !  FIND START AND STOP TIME OF DATA TO BE RETURNED
      !  --- THE OVERLAP OF THE REQUESTED AND DIGITIZED TIME INTERVALS

      tmind = tt(1)
      tmaxd = tt(nret)

      IF (tmax .le. tmin) THEN
        tmin = tmind
        tmax = tmaxd
      ENDIF

      tmin = max(tmin,tmind)
      tmax = min(tmax,tmaxd)

      if (.not. efit) then
        tmin_r4 = real(tmin,r4)
        tmax_r4 = real(tmax,r4)
      endif

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
      IF(NSAMPL.LT.ITOTALPTS) ITOTALPTS = NSAMPL

      DO i=1,ITOTALPTS
        ttime = tt(i)
        if (ttime .lt. tmin) nfirst=i+1
        if (ttime .lt. tmax) nlast =i
      ENDDO

      ndata = nlast - nfirst + 1

      !  DECIDE ON DATA POINT INTERVAL

      if (np .eq. 0) then
        np = 16384
      else
        np = amin0(np,npmax)
      endif

      nint = (ndata-1) / np + 1
      ! if (nint .gt. 1 .and. ier .le. 0) ier = -6
      if(nint.gt.1) ier = -6

      ! SOME OF THE ICH DATA  (ICHPWR) IS
      ! JUST LIKE DIGITIZER DATA EXCEPT IN
      ! REAL FORMAT

      if (idfi .eq. 119 .or. idfi.eq.2119) then      !ICH STUFF EXISTS IN FLOATING FORMAT
        ii=0
        if (efit) then
          do i = nfirst, nlast, nint
            ii=ii+1
            t(ii) = tt(i)
            data(ii) = fpdata(i)
          enddo
        else
          do i = nfirst, nlast, nint
            ii=ii+1
            t_r4(ii) = real(tt(i),r4)
            data_r4(ii) = real(fpdata(i),r4)
          enddo
        endif
        np = ii
        return
      endif

      !  CHECK FOR ZERO OFFSET OUT OF RANGE
      !  2013/04/03 Revised for new digitzers signed binary numbers  LL/ES

      ! amin = 0
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
      if (nbits.ge.16) then  !  NECESSARY ONLY FOR 16-BIT DATA
        amin =-(2 ** (nbits - 2) - 1) * 2 + 1
        amax = (2 ** (nbits - 2) - 1) * 2 + 1
      else    !  12-BIT DIG.S NEED THIS ONE
        amin = 0
        amax = 2**nbits - 1
      endif


      !--- WAVEFORMS ARE ALLOWED TO HAVE ZERO OFFSET OUTSIDE THE GENERATOR RANGE
      if ((mzero .lt. amin .or. mzero .gt. amax) .and. &
        (idfi .ne. 125 .and. idfi .ne. 158 .and. &
        idfi.ne.2158) .and. (ier .eq. 0 )) ier = -3

      !  CALIBRATION CODES TO RETURN ABSOLUTE SIGNAL LEVEL, NO BASELINE SUBTRACTION
      !  MUST ASSUME DIGITIZER ZERO IS IN THE MIDDLE OF ITS RANGE, SINCE THIS
      !  INFORMATION IS NOT STORED IN THE DATA BASE!
      !  WAVEFORMS (idfi = 125, 158, 2158) DO NOT HAVE A BASELINE, SINCE WE ARCHIVE
      !  ONLY THE PROGRAMMED VALUE AND NOT THE ACTUAL VALUE.
      ical0 = ical
      if (ical .ge. 10) then
        if(idfi .ne. 125 .and. idfi .ne. 158 &
          .and. idfi.ne.2158) mzero = 2**(nbits-1)
        yzero = mzero
        ical0 = ical-10
      endif

      ! FILL UP PLOT ARRAYS

      nclip = 0
      ii = 0

      do i=nfirst, nlast, nint

        ii = ii + 1
        !
        ! --- time array
        !
        if (efit) then
          t(ii) = tt(i)
        else
          t_r4(ii) = real(tt(i),r4)
        endif
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
        if(y .ge. amax .or. y .le. amin) nclip = nclip + 1

        ! --- offset for all but calibration code 0
        if (ical .ne. 0) then
          y = y - yzero
        !  y = y - yzero - deltab   ! ICAL > 0
        endif

        ! --  scale factors for the different calibration codes

        select case (ical0)
        case (0)
          y = y    ! ICAL = 0
        case (2)
          y = -y * vpbit   ! ICAL = 2
        case (3)
          y = -y * rcg / rc  ! ICAL = 3
        case (4)
          y = -y * rcg   ! ICAL = 4
        case default
          y = -y * rcg * zinhno + pzero ! ICAL = 1 or other
        end select

        y = scale * y * rcfact    ! user's scale factor

        ! --- arithmetic combinations as determined by mop code
        select case (mop)
        case (1)
          if (efit) then
            data(ii) = data(ii) + y
          else
            data_r4(ii) = data_r4(ii) + real(y,r4)
          endif
        case (2)
          if (efit) then
            data(ii) = data(ii) - y
          else
            data_r4(ii) = data_r4(ii) - real(y,r4)
          endif
        case (3)
          if (efit) then
            data(ii) = data(ii) * y
          else
            data_r4(ii) = data_r4(ii) * real(y,r4)
          endif
        case (4)
          if (y.ne.0) then
            if (efit) then
              data(ii) = data(ii) / y
            else
              data_r4(ii) = data_r4(ii) / real(y,r4)
            endif
          else
            if (efit) then
              data(ii) = 0.
            else
              data_r4(ii) = 0.
            endif
          endif
        case default ! 0, ...
          if (efit) then
            data(ii) = y
          else
            data_r4(ii) = real(y,r4)
          endif
        end select

      enddo

      ! --- actual number of points returned
      np = ii

      ! --- clipping tolerance:
      !     set error flag only if clipped fraction is outside tolerance
      if (np .gt. 0) then
        fclip = real(nclip,dp)/np
        if(fclip .gt. tolclip .and.  &
          ier .le. 0) ier = -1
      endif
      if (intcal) then
        !
        !   ADDITIONAL ARGUMENTS FOR INTEGRATOR CALIBRATION CALL.
        ! THESE DUMMY ARGUMENTS MUST ONLY BE REFERENCED WHEN USING THE
        ! ENTRY POINT GETDAT_I.  OTHERWISE ACCESS VIOLATION ERRORS OCCUR.
        !
        rc_i = rc
        rcg_i = rcg
        vpbit_i = vpbit
        mzero_i = mzero
        nbits_i = nbits
        write(sdate,880) iarr(19),month(iarr(18)), &
          iarr(20)-iarr(20)/100*100
        write(stime,881) iarr(15),iarr(16),iarr(17)
880     format(I2,'-',A3,'-',I2)
881     format(I2,':',I2,':',I2)
      !
      endif

      if (efit) then
        select case (ical0)
        case (0)
          bitone = 1                              ! ical = 0
        case (2)
          bitone = vpbit                          ! ical = 2
        case (3)
          bitone = rcg / rc                       ! ical = 3
        case default
          bitone = rcg * zinhno                   ! ical = 1 or other
        end select

        bitone = scale * abs(bitone)
        !--------------------------------------------------------------------
        !-- New returned values for new error matrix                       --
        !--------------------------------------------------------------------
        rc_e = rc
        rcg_e = rcg
        vpbit_e = vpbit
        zinhno_e = zinhno
        t0_e = tt(1)
      endif

      return

      CONTAINS
          !   ==============================
          !   Changes a string to upper case
          !   ==============================
          Function to_upper (str)

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
            if(ic > 0) to_upper(i:i) = cap(ic:ic)
          end do

          End Function to_upper
      END SUBROUTINE GETDAT

