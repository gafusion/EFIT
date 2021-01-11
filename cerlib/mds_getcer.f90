
      ! Subroutine to return CER data from MDS+.
      ! IER = 0 means that the ions tree was successfully opened for       
      ! the shot.  If the tree was successfully, we presume that some
      ! valid data are being stored in record cerdat.  
      ! IER = anything else implies that there was a problem in opening
      ! the tree and that we do not have the possibility for valid data.
      ! Written: 1/27/99 - rjg - patterned after mdstanh_getts
      subroutine mds_getcer(ishot,code,ier)
      use cermod
      use cer_results
      use standard_model
      implicit none
      integer*4 ishot           ! Requested shot number
      character*(*) code        ! Analysis code (CERFIT, CERNEUR, CERQUICK)
      Integer, intent(out) :: ier      ! Error flag

      Integer, intent(out) :: ntimes   ! Number of distinct slices of CER data
      
      Integer, intent(in) :: islice    ! Number of timeframe
      Real*4, intent(out) :: timeout   ! The corresponding time

      ! The local data.  Declare an object to hold data for one CER
      ! chord.
      record /cer_idl/ cerdat(nchord)
      integer :: stored_ishot = -1
      character*20 :: stored_code = ' '
      integer :: stored_ier = -1
      ! Array of distinct times
      Real*4, dimension(0:max_nfiles*maxnumtokfrms) :: time_array
      integer :: stored_ntimes = 0
      Save

      ! Functions to get mds data
      logical Mds_Start
      logical mds_open_tree
      character*10 tree         ! Name of tree
      character*50 node_arr      ! Name of node_ for array data
      character*50 node_cal      ! Name of node for calibration data
      integer idum
      integer ii
      integer lencstr           ! Function to get length of a string
      character*2 cc
      logical ok                ! Function to check if status return
                                ! from MDS function calls is okay

      Integer :: lier      ! Local error parameter
                                         
      ! PASSED VARIABLES
      Real time               ! IN: Requested time
      REAL time_win_upp       ! IN: Upper window on time
      Real time_win_low       ! IN: Lower window on time
      Record /CER_SLICE_DATA/ Cer_slice  ! OUT: Record with data for one time

      Integer ichord        ! Chord counter

      ! timing variables
      Logical Might_find       ! True if we could still find the frame
                               ! of interest
      Integer endfrm
      Integer goodframe
      Integer iframe
      Real    minmiss
      Real    miss
      Real    start            ! Start time of CER integration interval
      Real    avet

      Integer :: ich                   ! Chord number
      
      Integer, parameter :: ibozo   = -1000000
      Real, parameter ::    rbozo   = -1.0E30

      Character schord*6, fchord*10
      Integer :: ident

      Character*1 CTYPE             ! T or V
      Character*6 CNUM            ! The last characters of ascii chord name
      Integer Inum                ! Integer in a chord name

   
      call str$upcase(code,code)

      ! See if data are already loaded
      if (ishot .eq. stored_ishot .and.                                 &
     &    code  .eq. stored_code  .and.                                 &
     &    stored_ier .eq. 0) then
         print *,' Mds_GetCer: Data already loaded'
         ier = 0
         return 
      endif

      if (code .ne. 'CERFIT' .and. code .ne. 'CERQUICK' .and.           &
     &    code .ne. 'CERNEUR') then
          ok = .false.
          print *,' Mds_GetCer: Do not recognize code ', code
          ier = -101
          return
      endif  

      ! Connect to mds
      ! We need to have an argument for successful linking
      idum = 0  
      ok = Mds_Start (idum)
      if (.not. ok) then
          print *,' error from Mds_Start' 
          ier = -201
          return
      endif
      
      ! Open the ions tree
      tree = 'IONS'          
      ok = mds_open_tree (tree, ishot)
      if (.not. ok) then
          print *,' error from Mds_open_tree'
          ier = -301
          return
      endif

      
      ! If we arrive here, we assume that there has been some valid 
      ! CER analysis for this shot stored in MDS+.  Since there are three
      ! codes which can put data into MDS+, we do not yet know if we have
      ! the desired data.  We set ier to 0 if any one of the chords returns
      ! valid data.
      print *,' Mds_GetCer: attempting to load ', code,                &
     &        ' data from MDS+ for shot ', ishot


      Do ichord = 1, nchord
         ! Convert chord number to space chord           
         CALL ichord_to_schord (ichord, SCHORD) 
         ! Get the parts of schord
!         print *,' schord = ', schord 

         ctype = 'R'
         cnum = 'Z'
         Call GET_CHORD_NAME_PARTS (SCHORD,CTYPE,CNUM,LIER) 
         
         ! Write the number in a chord's ascii name into a character string 
         write (cc,fmt='(i2.2)') inum

         If (ctype .eq. 'T') then
            node_arr = '\TOP.CER.'//code(1:lencstr(code))//             &
     &             '.TANGENTIAL.CHANNEL'//cc//':'
            node_cal = '\TOP.CER.CALIBRATION'//                         &
     &             '.TANGENTIAL.CHANNEL'//cc//':'
         Else
            node_arr = '\TOP.CER.'//code(1:lencstr(code))//             &
     &             '.VERTICAL.CHANNEL'//cc//':'
            node_cal = '\TOP.CER.CALIBRATION'//                         & 
     &             '.VERTICAL.CHANNEL'//cc//':'
         Endif

         Call Ld_Cer_Idl_From_Mds (node_arr, node_cal,                  &
     &                             cerdat(ichord), ok) 

         ! Load up things not loaded at this timme in Ld_Cer_Idl_From_Mds
         ! Compute ascii chord name
         cerdat(ichord).schord = schord
         

      enddo


!      print *,'ntimes=',cerdat(33).ntimes
!      print *,'schord=',cerdat(33).schord
!      print *,'lamda=',cerdat.lamda
!      print *,'order=',cerdat.iord
!      print *,'fid=',cerdat.fid
!      print *,'sfid=',cerdat.sfid
!      print *,'gain',cerdat.gain
!      print *,'vpbit=',cerdat.vpbit
!      print *,'linid=',cerdat.linid
!      print *,'lens_r=',cerdat.lens_r
!      print *,'lens_z=',cerdat.lens_z
!      print *,'lens_phi=',cerdat.lens_phi
!      print *,'geom_030lt=',cerdat.geom_030lt
!      print *,'geom_030rt=',cerdat.geom_030rt
!      print *,'geom_330lt=',cerdat.geom_330lt
!      print *,'geom_330rt=',cerdat.geom_330rt
!      print *,'r=',cerdat.r
      ! Do we close anything? - no - mds+ does some nice buffering for us.

      ! We need to store the present shot and code, even though it is
      ! possible that we did not return any valid data.
      stored_ishot = ishot
      stored_code  = code

      ! Now setup the time basis
      Call Mds_Setup_cer_time_basis (nchord, cerdat,                    &
     &                               stored_ntimes, time_array)
 
      ! If we have some valid times, then everything we have a good
      ! error return
      If ( stored_ntimes .gt. 0) then
         ier = 0
      Else
         ier = -501
      Endif

      ! Remember the error
      stored_ier   = ier
     
      return



      ! Return number of stored frames
      ! Written: 3/4/99 by rjg
      Entry Mds_Fetch_Ntimes (ntimes)
         ntimes = stored_ntimes
      return



      ! Return a time from the time array
      ! Modified: 2/25/99 by rjg            
      Entry Mds_Fetch_time (islice,timeout,ier)
      
        ier=0

        If ((islice .lt. 0) .or. (islice .gt. max_ntimes)) Then
            ier=1
            print *,                                                      &
     &        ' MDS_FETCH_TIME:  call for time after end of time basis.'
                Write(*,*)'     Time index requested:',islice
                Write(*,*)'     Max_ntimes:',max_ntimes
                timeout=0
        Else
                timeout=time_array(islice)
        End if

        Return



      ! A routine to fetch data for a particular timeframe from the
      ! record with all of the CER data.  No caching of data is done
      ! here.  That should be done by the calling routine.
      ! This routine loads as much of CER_SLICE as possible from MDS
      ! data but some parts of the data will be missing because MDS
      ! does not store everything. 
      ! Modified: 2/25/99 by rjg - 
      Entry MDS_Fetch_Cer_Slice  (IShot, Code,                           &
     &                       Time, time_win_upp, time_win_low,           &
     &                       Cer_Slice, ier)


      ! Initialize record for an error return
      Cer_slice.inshot = -1
      Cer_slice.time = -1.0
      Do ichord=1, nchord
         Cer_slice.chord(ichord)=Cer_slice.chord(0)
         Cer_slice.chord(ichord).lhave_data = .false.
      Enddo


      ! Make sure that record CERDAT contains desired data 
      If (Ishot .ne. stored_ishot .or.                                  &
     &    Code  .ne. stored_code) then
         Print *,' MDS_Fetch_Cer_Slice: Correct data not loaded '
         ier = - 101
         return
      Endif

      ! Initialize
      ier = -99

      ! Find a frame we like for each chord, if possible
      Do ichord=1, nchord

         ! Make sure the file covers the right time range
         ! Handle the possibility that nframes = 0.
         endfrm=Cerdat(ichord).ntimes

         If (endfrm .gt. 0) then       
            If (Cerdat(ichord).time(1)  .le.                              & 
     &         (time_win_upp + time) .and.                                &
     &         (Cerdat(ichord).time(endfrm)                               &
     &          + Cerdat(ichord).stime(endfrm)) .ge.                      &
     &                           time_win_low + time) Then 

               !+
                  ! Find the frame.  
               ! Now check to see if time range spanned by TIME to
               ! START+AVET is in time window specified by
               ! TIME+TIME_WIN_LOW to TIME+TIME_WIN_UPP.
               ! This can happen if the CER frame straddles the lower
               ! boundary of the user window, straddles the upper boundary
               ! of the user window or lies entirely inside the user window.
               ! If so, increment counter and fill in arrays.  If START
               ! is larger than TIME+TIME_WIN_UPP, then CER file
               ! has no more possibilities.
               ! The following algorithm assumes that times monotonically
               ! increase.
               !-
               minmiss=100000
               goodframe=0
               iframe = 0
               might_find = .true. 

               Do while (iframe .lt. endfrm .and. might_find)       
                  iframe = iframe + 1   


                  ! Calculate start time of CER frame and integration time.
                  Start= Cerdat(ichord).time(iframe)
                  Avet = Cerdat(ichord).stime(iframe)

                  IF ((START .LE. (TIME+TIME_WIN_LOW) .AND.              &
     &               (START+AVET) .GE. (TIME+TIME_WIN_LOW))              &
     &                .OR.                                               &
     &               (START .GE. (TIME+TIME_WIN_LOW) .AND.               &
     &               (START+AVET) .LE. (TIME+TIME_WIN_UPP))              &
     &               .OR.                                                &
     &               (START .LE. (TIME+TIME_WIN_UPP) .AND.               &
     &               (START+AVET) .GE. (TIME+TIME_WIN_UPP)))             &
     &                                                      THEN
    
                     ! See if this time is closer to desired time
                     ! than any other time
                     Miss = start + 0.5* avet - time
                     If (ABS(miss) .lt. ABS(minmiss)) Then
                        minmiss=miss
                        goodframe=iframe
                     Endif
                  ! If frame start time has passed the window of
                  ! interest, we no longer have the chance to find the
                  ! desired time.  This keeps us from looping through
                  ! all frames every time into the routine.
                  Else IF (START .GT. (TIME+TIME_WIN_UPP)) THEN
                     might_find = .false.
                  End if                                ! Test of window
               End do                                ! do iframe


!               print *,' for chord ', ichord, ' goodframe = ',goodframe

               ! if miss is acceptable, store good frame in Cer_slice
               If (goodframe .gt. 0) Then
                   
                   ier = 0
                  
!                  print *,' Cer_slice.chord(ichord).frame.time = ',
!     &                    Cerdat(ichord).time(goodframe)
                  Cer_slice.chord(ichord).frame.time =                   &
     &                Cerdat(ichord).time(goodframe)
                  Cer_slice.chord(ichord).frame.avet =                   &
     &                Cerdat(ichord).stime(goodframe)
                  Cer_slice.chord(ichord).frame.nave =                   &
     &                Cerdat(ichord).nave(goodframe)
                  Cer_slice.chord(ichord).frame.ttsub =                  &
     &                Cerdat(ichord).ttsub(goodframe)
                  Cer_slice.chord(ichord).frame.ttsub_stime =            &
     &                Cerdat(ichord).ttsub_stime(goodframe)
                  Cer_slice.chord(ichord).frame.sub_frame = ibozo

                  Cer_slice.chord(ichord).frame.cpe = rbozo

                  Cer_slice.chord(ichord).frame.amp = rbozo
                  Cer_slice.chord(ichord).frame.samp = rbozo
                  Cer_slice.chord(ichord).frame.bright =                 &
     &                Cerdat(ichord).amp(goodframe)
                  Cer_slice.chord(ichord).frame.sbright =                &
     &                Cerdat(ichord).samp(goodframe)

                  Cer_slice.chord(ichord).frame.temp =                   &
     &                Cerdat(ichord).temp(goodframe)
                  Cer_slice.chord(ichord).frame.stemp =                  &
     &                Cerdat(ichord).stemp(goodframe)

                  Cer_slice.chord(ichord).frame.pix = ibozo
                  Cer_slice.chord(ichord).frame.spix = ibozo

                  Cer_slice.chord(ichord).frame.losvel =                 &
     &                Cerdat(ichord).rot(goodframe)
                  Cer_slice.chord(ichord).frame.slosvel =                &
     &                Cerdat(ichord).srot(goodframe)

                  Cer_slice.chord(ichord).lhave_data = .true.


                  ! pass in the header and vector structures
                  Cer_slice.chord(ichord).header.id_format = ibozo
                  Cer_slice.chord(ichord).header.inshot = stored_ishot
                                 
                  Call DECODE_SPACECHORD (                               &
     &                   Cerdat(ichord).schord,                          &
     &                  Cer_slice.chord(ichord).header.ichord,           & 
     &                  SCHORD, FCHORD, IER)
!              print *,' Cerdat(ichord).linid = ', Cerdat(ichord).linid 
                  Call Get_Ident_From_Linid (Cerdat(ichord).linid,ident)
!              print *,' ident = ', ident
                  Cer_slice.chord(ichord).header.Z = atomic_num(ident)
                  Cer_slice.chord(ichord).header.lamda =                 &
     &                  Cerdat(ichord).lamda
                  Cer_slice.chord(ichord).header.iord =                  &
     &                  Cerdat(ichord).iord
                  Cer_slice.chord(ichord).header.disp = rbozo
                  Cer_slice.chord(ichord).header.schord=                 &
     &                  Cerdat(ichord).schord
                  Cer_slice.chord(ichord).header.disp_err=  rbozo
                  Cer_slice.chord(ichord).header.radius =                &
     &                  Cerdat(ichord).r(goodframe)
                  Cer_slice.chord(ichord).header.sig_radius_low = rbozo 
                  Cer_slice.chord(ichord).header.sig_radius_high= rbozo 
                  Cer_slice.chord(ichord).header.realet =   rbozo
                  Cer_slice.chord(ichord).header.gain =                  &
     &                  Cerdat(ichord).gain
                  Cer_slice.chord(ichord).header.v_bit=                  &
     &                  Cerdat(ichord).vpbit
                  Cer_slice.chord(ichord).header.fid=                    &
     &                  Cerdat(ichord).fid
                  Cer_slice.chord(ichord).header.sfid=                   &
     &                  Cerdat(ichord).sfid
                  Cer_slice.chord(ichord).header.idoppler = ibozo 
                  Cer_slice.chord(ichord).header.amp_norm = rbozo
                  Cer_slice.chord(ichord).header.br = rbozo
                  Cer_slice.chord(ichord).header.bz = rbozo
                  Cer_slice.chord(ichord).header.btor = rbozo
                  Cer_slice.chord(ichord).header.bpol = rbozo
                  Cer_slice.chord(ichord).header.psi = rbozo
                  Cer_slice.chord(ichord).header.rho = rbozo
                  Cer_slice.chord(ichord).header.exeid = 'unknown'
                  Cer_slice.chord(ichord).header.linid =                 &
     &                  Cerdat(ichord).linid
                  Cer_slice.chord(ichord).header.cchord =                &
     &                  Cerdat(ichord).schord(1:1)
                  Cer_slice.chord(ichord).header.cbeam = 'HELP!'

                  Cer_slice.chord(ichord).beams.geom_030lt =             &
     &                  Cerdat(ichord).geom_030lt
                  Cer_slice.chord(ichord).beams.geom_030rt =             &
     &                  Cerdat(ichord).geom_030rt
                  Cer_slice.chord(ichord).beams.geom_330lt =             &
     &                  Cerdat(ichord).geom_330lt
                  Cer_slice.chord(ichord).beams.geom_330rt =             &
     &                  Cerdat(ichord).geom_330rt

                  Cer_slice.chord(ichord).vector.lens_r =                &
     &                  Cerdat(ichord).lens_r
                  Cer_slice.chord(ichord).vector.lens_z =                &
     &                  Cerdat(ichord).lens_z
                  Cer_slice.chord(ichord).vector.lens_phi =              &
     &                  Cerdat(ichord).lens_phi
                  Cer_slice.chord(ichord).vector.view_r =                &
     &                  Cerdat(ichord).r(goodframe)
                  Cer_slice.chord(ichord).vector.view_z =                &
     &                  Cerdat(ichord).z(goodframe)
                  Cer_slice.chord(ichord).vector.view_phi =              &
     &                  Cerdat(ichord).view_phi(goodframe)


                  ! Load the beam id with default info.  This needs to be fixed
                  If (ichord .ge. 1 .and. ichord .le. 8) then
                     Cer_slice.chord(ichord).header.cbeam = '330B'
                  Else If (ichord .ge. 9 .and. ichord .le. 15) then
                     Cer_slice.chord(ichord).header.cbeam = '030L'
                  Else If (ichord .ge. 16 .and. ichord .le. 32) then
                     Cer_slice.chord(ichord).header.cbeam = '330B'
                  Else If (ichord .ge. 41 .and. ichord .le. 46) then
                     Cer_slice.chord(ichord).header.cbeam = '030L'
                  Else If (ichord .ge. 47 .and. ichord .le. 48) then
                     Cer_slice.chord(ichord).header.cbeam = '330B'
                  Else
                     Cer_slice.chord(ichord).header.cbeam = 'what?'
                  Endif

                End if                                ! if (miss acceptable)

              End if                                ! if (right time range)
         Endif                                  ! For endfrm .gt. 0
      End do                                         ! do ichord


      ! Store slice info in Cer_slice
      If (ier .eq.0) then
         Cer_slice.inshot = ishot
         Cer_slice.time = time
         Cer_slice.lblessed = .true.
         Cer_slice.time_win_upp = time_win_upp
         Cer_slice.time_win_low = time_win_low
      Endif


      Return   ! Subroutine MDS_Fetch_Cer_Slice



      End  ! Subroutine Mds_Fetch_time




      ! Get the array of unique times from the CER data
      ! Copied from Setup_cer_time_basis
      ! Modified: 2/25/99 by rjg
      Subroutine Mds_Setup_cer_time_basis (nchd, cerdat, ntimes, time_array)

      use cermod
      use cer_results
      implicit none

      integer, intent(in) :: nchd             ! Number of chords    
      record /cer_idl/ cerdat(nchd)  ! The CER data
      integer, intent(out) :: ntimes            ! Number of times  
      Real*4, dimension(0:max_nfiles*maxnumtokfrms),                     &
     &       intent(out) :: time_array          ! Array of times


      Real avet_array(max_nfiles*maxnumtokfrms)
      Integer :: tindex                ! Time index
      Integer :: tindex2               ! Time index
      Integer :: ich                   ! Chord number
      Integer :: iframe                ! Frame number
      Real, parameter :: junk_value = 100000.0

       
        ! This fills up the time_array with the times at which
        !  any data exist, in increasing order and returns
        !  an integer telling how many distinct times there are

        ! If the 0the time isn't 1, we need to do a time basis;
        !  otherwise, no new data has been loaded in, and we
        !  can use the old time basis
        tindex=0

        ! load the numbers into the array
        Do ich=1, nchd
           Do iframe=1, cerdat(ich).ntimes
              tindex=tindex+1
              time_array(tindex)=                                          &
     &           cerdat(ich).time(iframe) + 0.5*cerdat(ich).stime(iframe)
              avet_array(tindex) = cerdat(ich).stime(iframe)
           End do
        End do

         ! The number of times is what tindex is left at
        ntimes=tindex

        ! fill up rest of array with junk information
        Do tindex2=ntimes+1,maxnumtokfrms*max_nfiles
           time_array(tindex2)=junk_value
        End do

        ! Sort them
        call uniq_sort(time_array, ntimes, junk_value)
        if (time_array(ntimes) .eq. junk_value) then
           ntimes = ntimes - 1
        end if

        Write (*,*) 'MDS_SETUP_CER_TIME_BASIS:  ntimes=',ntimes
        time_array(0)=ntimes

      End 





      ! Routine to load up a structure of type CER_IDL with data
      ! for one chord from MDS+
      ! If OK is true, all data were obtained successfully.
      ! If OK is false, some or all data may be missing.
      ! Written: Feb 18, 1999 by rjg
      Subroutine Ld_Cer_Idl_From_Mds (node_arr, node_cal, cerdat, ok) 
      use cermod
      use cer_results
      implicit none
!      include 'parameters.inc'
!      include 'cer_basic_data.inc'
!      include 'cer_chord_data.inc'

      character*(*) node_arr      ! IN: MDS node_ for array data
      character*(*) node_cal      ! IN: MDS node for calibration data
      record /cer_idl/ cerdat     ! OUT: Object with CER data
      logical :: ok               ! OUT: TRUE means success 

      ! Things for call to mds_real_array
      integer :: maxarr                                ! Size of rarray
      real*4, dimension(size(cerdat.time)) :: rarray   ! Array to hold goodies
      character*10, dimension(size(cerdat.time)) :: cargarr ! Array to hold goodies
      integer :: narr                                  ! Number of values ret
      character*20 :: signal                           ! MDS signal
              
      integer :: ii 
      real :: rarg
      integer :: iarg
      character*20 :: carg

      maxarr = size (rarray)
 
      ! Initialize the structure
      Call Init_for_Idl (cerdat)

      ! Load array data
      ok = .true.
  
      ! Make sure that we have some data
      signal = 'TIME'  ! (ms)
      call mds_real_array (node_arr, signal, maxarr, rarray, narr, ok)
      if (.not. ok .or. narr .le. 0) then
!         print *, ' No data for signal ', signal
         ok = .false.
         return
      else
         cerdat.ntimes = narr
         cerdat.time = rarray    
      endif

      ! Load up averaging time (ms)
      signal = 'STIME'
      call mds_real_array (node_arr, signal, maxarr, rarray, narr, ok)
      if (ok) cerdat.stime = rarray    

      ! Load up intensity (photons/sec/m**2/sr)
      signal = 'INTENSITY'
      call mds_real_array (node_arr, signal, maxarr, rarray, narr, ok)
      if (ok) cerdat.amp = rarray    

      ! Load up intensity error (photons/sec/m**2/sr)
      signal = 'INTENSITYERR'
      call mds_real_array (node_arr, signal, maxarr, rarray, narr, ok)
      if (ok) cerdat.samp = rarray    

      ! Load up radius (m)
      signal = 'R'
      call mds_real_array (node_arr, signal, maxarr, rarray, narr, ok)
      if (ok) cerdat.r = rarray

      ! Load up line-of-sight velocity (km/sec)
      signal = 'ROT'
      call mds_real_array (node_arr, signal, maxarr, rarray, narr, ok)
      if (ok) cerdat.rot = rarray    

      ! Load up error in line-of-sight velocity (km/sec)
      signal = 'ROT_ERR'
      call mds_real_array (node_arr, signal, maxarr, rarray, narr, ok)
      if (ok) cerdat.srot = rarray    

      ! Load up temperature and convert from eV to keV
      signal = 'TEMP'
      call mds_real_array (node_arr, signal, maxarr, rarray, narr, ok)
      if (ok) cerdat.temp = 0.001*rarray    

      ! Load up temperature error and convert from eV to keV
      signal = 'TEMP_ERR'
      call mds_real_array (node_arr, signal, maxarr, rarray, narr, ok)
      if (ok) cerdat.stemp = 0.001*rarray

      ! Load up subtraction slice
      signal = 'TTSUB'
      call mds_real_array (node_arr, signal, maxarr, rarray, narr, ok)
      if (ok) cerdat.ttsub = rarray    

      ! Load up averaging time of subtraction slice (ms)
      signal = 'TTSUB_STIME'
      call mds_real_array (node_arr, signal, maxarr, rarray, narr, ok)
      if (ok) cerdat.ttsub_stime = rarray    

      ! Load up toroidal angle at viewing position (degrees)
      signal = 'VIEW_PHI'
      call mds_real_array (node_arr, signal, maxarr, rarray, narr, ok)
      if (ok) cerdat.view_phi = rarray

      ! Load up elevation (m)
      signal = 'Z'
      call mds_real_array (node_arr, signal, maxarr, rarray, narr, ok)
      if (ok) cerdat.z = rarray

      ! Get some 'calibration' data

      ! Load up beam pathlength info
      signal = 'BEAMGEOMETRY'
      call mds_real_array (node_cal, signal, maxarr, rarray, narr, ok)
      if (ok) then
         cerdat.geom_030lt = rarray(1)
         cerdat.geom_030rt = rarray(2)
         cerdat.geom_330lt = rarray(3)
         cerdat.geom_330rt = rarray(4)
      endif  

      ! Load up MCP control voltage (V)
      signal = 'GAIN'
      call mds_float (node_cal, signal, rarg, ok)
      if (ok) cerdat.gain = rarg

      ! Load up fiducial error (pixel)
      signal = 'FIDUCIAL_ERR'
      call mds_float (node_cal, signal, rarg, ok)
      if (ok) cerdat.sfid = rarg

      ! Load up fiducial (pixel)
      signal = 'FIDUCUAL'           ! Mis-spelled in MDS tree
      call mds_float (node_cal, signal, rarg, ok)
      if (ok) cerdat.fid = rarg

      ! Load up lens radius (m)
      signal = 'LENS_R'
      call mds_float (node_cal, signal, rarg, ok)
      if (ok) cerdat.lens_r = rarg

      ! Load up lens radius (m)
      signal = 'LENS_Z'
      call mds_float (node_cal, signal, rarg, ok)
      if (ok) cerdat.lens_z = rarg

      ! Load up lens radius (m)
      signal = 'LENS_PHI'
      call mds_float (node_cal, signal, rarg, ok)
      if (ok) cerdat.lens_phi = rarg

      ! Load up order!! (stored as real it appears!)
      signal = 'ORDER'
      call mds_int (node_cal, signal, iarg, ok)
      if (ok) cerdat.iord = iarg

      ! Load up digitizer volt per bit conversion
      signal = 'VOLTSPERBIT'
      call mds_float (node_cal, signal, rarg, ok)
      if (ok) cerdat.vpbit = rarg

      ! Load up wavelength (Angstroms)
      signal = 'WAVELENGTH'
      call mds_float (node_cal, signal, rarg, ok)
      if (ok) cerdat.lamda = rarg

      ! Load up line identification
      signal = 'LINEID'
!      call mds_char (node_cal, signal, carg, ok)
!      print *,' carg = ', carg(1:20)
!      if (ok) cerdat.linid = carg
!      print *,' ok = ', ok
!      print *,' cerdat.linid = ', cerdat.linid


      ! Load up beam identification.
      ! Don't bother calling MDS+ - we are not loading anything there
      ! so far as I know.
!      signal = 'BEAMID'
!      call mds_char_arr (node_cal, signal, cargarr, ok)
!      if (ok) cerdat.linid = carg
!      print *,' ok = ', ok
!      print *,' cerdat.linid = ', cerdat.linid

      End subroutine ! Ld_Cer_Idl_From_Mds


      ! Compute ascii name of a chord
      Subroutine Make_Schord (Carg,inum,schord)
      implicit none
      character*4 carg              ! VERT or TANG
      integer inum                  ! The chord number
      character*6 schord 

      if (inum .lt. 10) then
         write (schord, fmt='(A4,I1)') carg,inum
      else
         write (schord, fmt='(A4,I2)') carg,inum
      endif

      end


  
      ! Compare two structures
      subroutine print_things (ic, slice_files, slice)
      use cermod
      use cer_results
      implicit none
      
!      include 'parameters.inc'
!      include 'cer_basic_data.inc'
!      include 'cer_chord_data.inc'
!      include 'cer_file_data.inc'
!      include 'cer_slice_data.inc'

      integer*4 ic
      Record /CER_SLICE_DATA/ slice  
      Record /CER_SLICE_DATA/ slice_files


      print *,' printing results for chord ',  ic

      print *,'                 From files            From MDS+'

      print *,' inshot = ', slice_files.inshot, slice.inshot
      print *,' lblessed = ', slice_files.lblessed, slice.lblessed
      print *,' time = ', slice_files.time, slice.time
      print *,' TIME_WIN_LOW=', slice_files.TIME_WIN_LOW,               &
     &                          slice.TIME_WIN_LOW
      print *,' TIME_WIN_UPP=', slice_files.TIME_WIN_UPP,               &
     &                          slice.TIME_WIN_UPP
      print *,' nchord = ', slice_files.nchords, slice.nchords


      print *,' lhave_data = ', slice_files.chord(ic).lhave_data,       &
     &                          slice.chord(ic).lhave_data
      print *,' id_format = ', slice_files.chord(ic).header.id_format,  &
     &                         slice.chord(ic).header.id_format

      print *,' inshot = ', slice_files.chord(ic).header.inshot,        &
     &                         slice.chord(ic).header.inshot

      print *,' ichord = ', slice_files.chord(ic).header.ichord,        &
     &                         slice.chord(ic).header.ichord

      print *,' z = ', slice_files.chord(ic).header.z,                  & 
     &                         slice.chord(ic).header.z

      print *,' lamda = ', slice_files.chord(ic).header.lamda,          &
     &                         slice.chord(ic).header.lamda

      print *,' iord = ', slice_files.chord(ic).header.iord,            &
     &                         slice.chord(ic).header.iord

      print *,' disp = ', slice_files.chord(ic).header.disp,            &
     &                         slice.chord(ic).header.disp

      print *,' schord = ', slice_files.chord(ic).header.schord,        &
     &                         slice.chord(ic).header.schord

      print *,' disp_err = ', slice_files.chord(ic).header.disp_err,    &
     &                         slice.chord(ic).header.disp_err

      print *,' radius = ', slice_files.chord(ic).header.radius,        &
     &                         slice.chord(ic).header.radius

      print *,' sig_radius_low = ',                                     &
     &                    slice_files.chord(ic).header.sig_radius_low,  &
     &                    slice.chord(ic).header.sig_radius_low

      print *,' sig_radius_high = ',                                    &
     &                    slice_files.chord(ic).header.sig_radius_high, &
     &                    slice.chord(ic).header.sig_radius_high    

      print *,' realet = ', slice_files.chord(ic).header.realet,        &
     &                         slice.chord(ic).header.realet

      print *,' gain = ', slice_files.chord(ic).header.gain,            &
     &                         slice.chord(ic).header.gain

      print *,' v_bit = ', slice_files.chord(ic).header.v_bit,          &
     &                         slice.chord(ic).header.v_bit

      print *,' fid = ', slice_files.chord(ic).header.fid,              &
     &                         slice.chord(ic).header.fid
 
      print *,' sfid = ', slice_files.chord(ic).header.sfid,            &
     &                         slice.chord(ic).header.sfid

      print *,' idoppler = ', slice_files.chord(ic).header.idoppler,    &
     &                         slice.chord(ic).header.idoppler

      print *,' amp_norm = ', slice_files.chord(ic).header.amp_norm,    &
     &                         slice.chord(ic).header.amp_norm

      print *,' br = ', slice_files.chord(ic).header.br,                &
     &                         slice.chord(ic).header.br

      print *,' bz = ', slice_files.chord(ic).header.bz,                &
     &                         slice.chord(ic).header.bz

      print *,' btor = ', slice_files.chord(ic).header.btor,            &
     &                         slice.chord(ic).header.btor

      print *,' bpol = ', slice_files.chord(ic).header.bpol,            &
     &                         slice.chord(ic).header.bpol

      print *,' psi = ', slice_files.chord(ic).header.psi,              & 
     &                         slice.chord(ic).header.psi

      print *,' rho = ', slice_files.chord(ic).header.rho,              &
     &                         slice.chord(ic).header.rho

      print *,' exeid = ', slice_files.chord(ic).header.exeid,          &
     &                         slice.chord(ic).header.exeid

      print *,' linid = ', slice_files.chord(ic).header.linid,          &
     &                         slice.chord(ic).header.linid

      print *,' cchord = ', slice_files.chord(ic).header.cchord,        &
     &                         slice.chord(ic).header.cchord

      print *,' cbeam = ', slice_files.chord(ic).header.cbeam,          &
     &                         slice.chord(ic).header.cbeam

      print *,' geom_030lt = ', slice_files.chord(ic).beams.geom_030lt, &
     &                         slice.chord(ic).beams.geom_030lt

      print *,' geom_030rt = ', slice_files.chord(ic).beams.geom_030rt, &
     &                         slice.chord(ic).beams.geom_030rt

      print *,' geom_330rt = ', slice_files.chord(ic).beams.geom_330rt, &
     &                         slice.chord(ic).beams.geom_330rt

      print *,' geom_330lt = ', slice_files.chord(ic).beams.geom_330lt, &
     &                         slice.chord(ic).beams.geom_330lt

      print *,' lens_r = ', slice_files.chord(ic).vector.lens_r,        &
     &                         slice.chord(ic).vector.lens_r

      print *,' lens_z = ', slice_files.chord(ic).vector.lens_z,        &
     &                         slice.chord(ic).vector.lens_z

      print *,' lens_phi = ', slice_files.chord(ic).vector.lens_phi,    &
     &                         slice.chord(ic).vector.lens_phi

      print *,' view_r = ', slice_files.chord(ic).vector.view_r,        &
     &                         slice.chord(ic).vector.view_r

      print *,' view_z = ', slice_files.chord(ic).vector.view_z,        &
     &                         slice.chord(ic).vector.view_z

      print *,' view_phi = ', slice_files.chord(ic).vector.view_phi,    &
     &                         slice.chord(ic).vector.view_phi
 
      print *,' TIME = ', slice_files.chord(ic).frame.time,             &
     &                         slice.chord(ic).frame.time

      print *,' avet = ', slice_files.chord(ic).frame.avet,             &
     &                         slice.chord(ic).frame.avet

      print *,' nave = ', slice_files.chord(ic).frame.nave,             &
     &                         slice.chord(ic).frame.nave

      print *,' ttsub = ', slice_files.chord(ic).frame.ttsub,           &
     &                         slice.chord(ic).frame.ttsub

      print *,' ttsub_stime = ',                                        & 
     &                      slice_files.chord(ic).frame.ttsub_stime,    &
     &                         slice.chord(ic).frame.ttsub_stime

      print *,' sub_frame = ', slice_files.chord(ic).frame.sub_frame,   &
     &                         slice.chord(ic).frame.sub_frame

      print *,' cpe = ', slice_files.chord(ic).frame.cpe,               &
     &                         slice.chord(ic).frame.cpe

      print *,' amp = ', slice_files.chord(ic).frame.amp,               &
     &                         slice.chord(ic).frame.amp

      print *,' samp = ', slice_files.chord(ic).frame.samp,             &
     &                         slice.chord(ic).frame.samp

      print *,' temp = ', slice_files.chord(ic).frame.temp,             &
     &                         slice.chord(ic).frame.temp

      print *,' stemp = ', slice_files.chord(ic).frame.stemp,           &
     &                         slice.chord(ic).frame.stemp

      print *,' bright = ', slice_files.chord(ic).frame.bright,         &
     &                         slice.chord(ic).frame.bright

      print *,' sbright = ', slice_files.chord(ic).frame.sbright,       &
     &                         slice.chord(ic).frame.sbright

      print *,' pix = ', slice_files.chord(ic).frame.pix,               &
     &                         slice.chord(ic).frame.pix

      print *,' spix = ', slice_files.chord(ic).frame.spix,             &
     &                         slice.chord(ic).frame.spix

      print *,' losvel = ', slice_files.chord(ic).frame.losvel,         &
     &                         slice.chord(ic).frame.losvel

      print *,' slosvel = ', slice_files.chord(ic).frame.slosvel,       &
     &                         slice.chord(ic).frame.slosvel

      end



      subroutine blah !program mdstanh
      use cermod
      use cer_results
      implicit none
!      include 'parameters.inc'
!      include 'cer_basic_data.inc'
!      include 'cer_chord_data.inc'
!      include 'cer_file_data.inc'
!      include 'cer_slice_data.inc'


      integer*4 ishot
      integer*4 ntimes
      integer*4 islice
      integer*4 ic
      real*4 time, ztime
      integer*4 ierr
      logical ok
      character*10 code
      real time_win_upp, time_win_low
      Record /CER_SLICE_DATA/ slice  

      ! Things to get data from files, via the old access mechanism  
      Logical lblessed
      Logical LAMP_NORM  
      Integer IUNIT_CER
      Record /CER_SLICE_DATA/ slice_files

         ! Get data from mds+
  10     CONTINUE
         print *,' enter the SHOT '
         read (*,*) ISHOT

         code = 'CERQUICK' 
!         code = 'CERFIT'  
         call mds_getcer (ishot,code,ierr)
         Print *,'  ierr from MDS_getcer = ', ierr
         call Mds_Fetch_ntimes (ntimes)


         ! Get the same data from cerbl files
         lblessed = .true.
         lamp_norm = .false.
         iunit_cer = 55      
         CALL LD_SHOT_DATA(ISHOT, LBLESSED, IUNIT_CER, LAMP_NORM, IERR)
         Print *,'  error from LD_SHOT_DATA = ', ierr
         call Setup_cer_time_basis (ntimes)

         islice = 11           
         call Mds_Fetch_time (islice,ztime,ierr)
         print *,' for islice = ', islice, 'Mds_Fetch_time gives time',  &
     &                              ztime
         call Fetch_time (islice,ztime,ierr)
         print *,' for islice = ', islice, 'Fetch_time gives time',      &
     &                              ztime

         print *,' enter the time '
         read (*,*) time 

         print *,' time = ', time

         time_win_upp = 30.0
         time_win_low = -45.0

         call MDS_Fetch_Cer_Slice  (ishot, code, Time,                  &
     &                        time_win_upp, time_win_low,               &
     &                        Slice, ierr)
         Print *,'  error from MDS_Fetch_Cer_Slice = ', ierr

         call Fetch_Slice  (ishot, lblessed, time,                      &
     &           time_win_low, time_win_upp,                            &
     &           slice_files, ierr)
         Print *,'  error from Fetch_Slice = ', ierr
  
         print *,' enter ichord number '
         read (*,*) ic 

         call print_things (ic, slice_files, slice)

         goto 10      

      end

!
