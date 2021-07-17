!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          msels_hist returns synthetic or experimental MSE-LS     **
!**          data to EFIT                                            **
!**          For a given shot, returns complete time history         **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**          ishot: shot number                                      **
!**          synmls: synthetic 'SYN' or experimental 'EXP' data      **
!**                  or 'MSE' for EFIT02 MSE locations               **
!**          icmls: MSE-LS channel number                            **
!**          time: time array
!**          bbmls: MSE-LS magnetic field E/V_BEAM in Tesla          **
!**          sigbmls: MSE-LS magnetic field uncertainty in Tesla     **
!**          rrmls, zzmls: (R,Z) of MSE-LS channels in m             **
!**          L1mls, L2mls, L4mls: geometric cofficients of MSE-LS    **
!**                               channels                           **
!**          epotpmls: flux derivative of electrostatic potential    **
!**                    in 1/sec                                      **
!**          sigepmls: uncertainty in electrostatic potential flux   **
!**                    derivative in 1/sec                           **
!**          iermls: error flags of MSE-LS channels, 0=normal,       **
!**                  1=error (ignore data)                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          2015/02/12..........first created                       **
!**                                                                  **
!**********************************************************************
      SUBROUTINE  MSELS_HIST(ishot,synmls,icmls,ntimes,l_time,        &
           l_bbmls,l_sigbmls,l_rrmls,l_zzmls,l_L1mls,l_L2mls,         &
           l_L4mls,l_epotpmls,l_sigepmls,iermls)

      INTEGER*4 ishot,icmls,iermls,i
      REAL*4 time,bbmls,sigbmls,rrmls,zzmls, &
           L1mls,L2mls,L4mls, &
           epotmls,sigepmls
      CHARACTER*3 synmls
      
      ! Local variables
      LOGICAL file_exists
      CHARACTER(len=100) filename
      INTEGER*4 iostat,file_shot,ntimes,issyn,l_ntimes
      REAL*4 l_time(ntimes),l_bbmls(ntimes),l_sigbmls(ntimes),&
           l_rrmls(ntimes),l_zzmls(ntimes),l_L1mls(ntimes),l_L2mls(ntimes),&
           l_L4mls(ntimes),l_epotpmls(ntimes),l_sigepmls(ntimes)

      ! Define the filename that has the time history beginning
      IF (synmls == 'SYN') THEN
         WRITE(filename,"(A23,I6,A15,I2.2,A4)") "/u/grierson/efit/msels/", &
              ishot,"/msels_syn_chan",icmls,".dat"
      END IF

      IF (synmls == 'MSE') THEN
         WRITE(filename,"(A23,I6,A15,I2.2,A4)") "/u/grierson/efit/msels/", &
              ishot,"/msels_mse_chan",icmls,".dat"
      END IF

      IF (synmls == 'EXP') THEN
         WRITE(filename,"(A23,I6,A11,I2.2,A4)") "/u/grierson/efit/msels/", &
              ishot,"/msels_chan",icmls,".dat"
      END IF

      WRITE(*,*) "Using file ",filename

      ! See if file exists.  If not, return error code.
      INQUIRE(FILE=filename,EXIST=file_exists)
      IF (.NOT. file_exists) THEN
         WRITE(*,*) 'File Does not exist'
         iermls=1
         RETURN
      END IF

      ! Read the data
      OPEN(UNIT=1,FILE=filename,IOSTAT=iostat)
      WRITE(*,101) 'iostat: ',iostat
      ! Read the shot number in the file
      READ(1,*) file_shot
      WRITE(*,101) 'file shot:', file_shot

      ! Read the number of timeslices
      READ(1,*) l_ntimes
      WRITE(*,101) '# slices:', l_ntimes

      ! Read if synthetic(1) or measured(0)
      READ(1,*) issyn
      WRITE(*,101) 'issyn: ', issyn

      !      Assign logical to map 1->true and 0->false
      !      This doesn't work.
      !      IF (issyn .EQ. 1) THEN 
      !         synmls = .TRUE.
      !      END IF
      
      ! Read the channel number (1,2,etc...)
      READ(1,*) icmls
      WRITE(*,101) 'icmls: ', icmls

      ! Read the error code
      READ(1,*) iermls
      WRITE(*,101) 'iermls: ', iermls

      WRITE(*,100) 'Starting time loop'
      DO i=1,l_ntimes

         READ(1,*) l_time(i)
!         WRITE(*,102) 'Time:', l_time(i)

         READ(1,*) l_bbmls(i)
!         WRITE(*,102) 'B_LS:', l_bbmls(i)
         
         READ(1,*) l_sigbmls(i)
!         WRITE(*,102) 'sig B_LS:', l_sigbmls(i)
         
         READ(1,*) l_rrmls(i)
!         WRITE(*,102) 'R:', l_rrmls(i)
         
         READ(1,*) l_zzmls(i)
!         WRITE(*,102) 'Z:', l_zzmls(i)
         
         READ(1,*) l_L1mls(i)
!         WRITE(*,102) 'L1:', l_L1mls(i)
         
         READ(1,*) l_L2mls(i)
!         WRITE(*,102) 'L2:', l_L2mls(i)
         
!!         L3mls(i) = 1.0 - L1mls(i)
         
         READ(1,*) l_L4mls(i)
!         WRITE(*,102) 'L4:', l_L4mls(i)
         
         READ(1,*) l_epotpmls(i)
!         WRITE(*,102) 'epot:', l_epotpmls(i)
         
         READ(1,*) l_sigepmls(i)
!         WRITE(*,102) 'sigep:', l_sigepmls(i)
         
      END DO
      CLOSE(1)

100   FORMAT (A)
101   FORMAT (A, I6)
102   FORMAT (A, g15.4)
      RETURN

      END SUBROUTINE MSELS_HIST
