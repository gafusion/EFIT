!**********************************************************************
!>
!!    msels_hist returns synthetic or experimental MSE-LS
!!    data to EFIT
!!    For a given shot, returns complete time history
!!
!!    Note: this subroutine only uses REAL*4 for floating point
!!          variables to be consistent with mse files
!!
!!
!!    @param ishot : shot number   
!!
!!    @param synmls : synthetic 'SYN' or experimental 'EXP' data 
!!                  or 'MSE' for EFIT02 MSE locations  
!!
!!    @param icmls : MSE-LS channel number 
!!
!!    @param ntimes :  number of times
!!
!!    @param l_time : time array
!!
!!    @param l_bbmls : MSE-LS magnetic field E/V_BEAM in Tesla 
!!
!!    @param l_sigbmls :  MSE-LS magnetic field uncertainty in Tesla
!!
!!    @param l_rrmls : R of MSE-LS channels in m  
!!
!!    @param l_zzmls : Z of MSE-LS channels in m  
!!
!!    @param l_L1mls :  geometric cofficients of MSE-LS channels 
!!
!!    @param l_L2mls :  geometric cofficients of MSE-LS channels 
!!
!!    @param l_L4mls :  geometric cofficients of MSE-LS channels 
!!
!!    @param l_epotpmls : flux derivative of electrostatic potential
!!                    in 1/sec
!!
!!    @param l_sigepmls : uncertainty in electrostatic potential flux
!!                   derivative in 1/sec 
!!
!!    @param iermls :error flags of MSE-LS channels, 0=normal,  
!!                 1=error (ignore data)
!!
!**********************************************************************
      subroutine msels_hist(ishot,synmls,icmls,ntimes,l_time, &
           l_bbmls,l_sigbmls,l_rrmls,l_zzmls,l_L1mls,l_L2mls, &
           l_L4mls,l_epotpmls,l_sigepmls,iermls)

      integer*4 ishot,icmls,iermls,i
      real*4 time,bbmls,sigbmls,rrmls,zzmls, &
             L1mls,L2mls,L4mls, &
             epotmls,sigepmls
      character*3 synmls
      
      ! Local variables
      logical file_exists
      character(len=100) filename
      integer*4 iostat,file_shot,ntimes,issyn,l_ntimes
      real*4 l_time(ntimes),l_bbmls(ntimes),l_sigbmls(ntimes), &
             l_rrmls(ntimes),l_zzmls(ntimes),l_L1mls(ntimes), &
             l_L2mls(ntimes),l_L4mls(ntimes),l_epotpmls(ntimes), &
             l_sigepmls(ntimes)

      ! Define the filename that has the time history beginning
      if(synmls == 'SYN') &
         write(filename,"(A23,I6,A15,I2.2,A4)") &
          "/u/grierson/efit/msels/",ishot,"/msels_syn_chan",icmls,".dat"

      if(synmls == 'MSE') &
         write(filename,"(A23,I6,A15,I2.2,A4)") &
          "/u/grierson/efit/msels/",ishot,"/msels_mse_chan",icmls,".dat"

      if(synmls == 'EXP') &
         write(filename,"(A23,I6,A11,I2.2,A4)") &
          "/u/grierson/efit/msels/",ishot,"/msels_chan",icmls,".dat"

      write(*,*) "Using file ",filename

      ! See if file exists.  If not, return error code.
      inquire(file=filename,exist=file_exists)
      if (.not. file_exists) then
         write(*,*) 'File Does not exist'
         iermls=1
         return
      endif

      ! Read the data
      open(unit=1,file=filename,iostat=iostat)
      write(*,101) 'iostat: ',iostat
      ! Read the shot number in the file
      read(1,*) file_shot
      write(*,101) 'file shot:', file_shot

      ! Read the number of timeslices
      read(1,*) l_ntimes
      write(*,101) '# slices:', l_ntimes

      ! Read if synthetic(1) or measured(0)
      read(1,*) issyn
      write(*,101) 'issyn: ', issyn

      !      Assign logical to map 1->true and 0->false
      !      This doesn't work.
      !      IF (issyn .EQ. 1) THEN 
      !         synmls = .TRUE.
      !      END IF
      
      ! Read the channel number (1,2,etc...)
      read(1,*) icmls
      write(*,101) 'icmls: ', icmls

      ! Read the error code
      read(1,*) iermls
      write(*,101) 'iermls: ', iermls

      write(*,100) 'Starting time loop'
      do i=1,l_ntimes

         read(1,*) l_time(i)
!         write(*,102) 'Time:', l_time(i)

         read(1,*) l_bbmls(i)
!         write(*,102) 'B_LS:', l_bbmls(i)
         
         read(1,*) l_rrmls(i)
!         write(*,102) 'sig B_LS:', l_sigbmls(i)
         
         read(1,*) l_zzmls(i)
!         write(*,102) 'Z:', l_zzmls(i)
         
         read(1,*) l_L1mls(i)
!         write(*,102) 'L1:', l_L1mls(i)
         
         read(1,*) l_L2mls(i)
!         write(*,102) 'L2:', l_L2mls(i)
         
!!         L3mls(i) = 1.0 - L1mls(i)
         
         read(1,*) l_L4mls(i)
!         write(*,102) 'L4:', l_L4mls(i)
         
         read(1,*) l_epotpmls(i)
!         write(*,102) 'epot:', l_epotpmls(i)
         
         read(1,*) l_sigepmls(i)
!         write(*,102) 'sigep:', l_sigepmls(i)
         
      enddo
      close(1)

100   format (a)
101   format (a, i6)
102   format (a, g15.4)
      return

      end subroutine msels_hist
