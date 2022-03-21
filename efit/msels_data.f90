!**********************************************************************
!>
!!    msels_data returns synthetic or experimental MSE-LS
!!    data to EFIT
!!    For list of k-file times returs the available MSE-LS
!!    data for a given chord, averaged over avgtim.
!!
!!    Note: this subroutine only uses REAL*4 for floating point
!!          variables to be consistent with mse files
!!
!!
!!    @param ishot : shot number
!!
!!    @param atime : time array for MSE-LS data (k-file times) 
!!
!!    @param ktime : number of timeslices (number of k-files) 
!!
!!    @param avem : +- averaging time for MSE-LS data  
!!
!!    @param synmls : synthetic 'SYN' or experimental 'EXP' data
!!                  or 'MSE' for EFIT02 MSE locations 
!!
!!    @param icmls : MSE-LS channel number  
!!
!!    @param bbmls : MSE-LS magnetic field E/V_BEAM in Tesla
!!
!!    @param sigbmls : MSE-LS magnetic field uncertainty in Tesla
!!
!!    @param rrmls : R of MSE-LS channels in m 
!!
!!    @param zzmls : Z of MSE-LS channels in m 
!!
!!    @param L1mls : geometric cofficients of MSE-LS channels
!!
!!    @param L2mls : geometric cofficients of MSE-LS channels
!!
!!    @param L4mls : geometric cofficients of MSE-LS channels
!!
!!    @param epotpmls : flux derivative of electrostatic potential
!!                   in 1/sec 
!!
!!    @param sigepmls : uncertainty in electrostatic potential flux
!!                   derivative in 1/sec 
!!
!!    @param iermls : error flags of MSE-LS channels, 0=normal,
!!                 1=error (ignore data) 
!!
!**********************************************************************
      subroutine msels_data(ishot,atime,ktime,avem,synmls,icmls, &
             bbmls,sigbmls,rrmls,zzmls,L1mls,L2mls,L4mls,epotpmls, &
             sigepmls,iermls) 
      integer*4, parameter :: ntimes=2000
      character*3 synmls
      integer*4 i,j,count
      integer*4 ishot,ktime,icmls,iermls(ktime),l_iermls
      real*4 avem,atime(ktime),bbmls(ktime),sigbmls(ktime), &
             rrmls(ktime),zzmls(ktime), &
             L1mls(ktime),L2mls(ktime),L4mls(ktime), &
             epotpmls(ktime),sigepmls(ktime)

      ! local variables
      character(len=100) filename
      integer*4 file_shot
      real*4 :: l_time(ntimes),l_bbmls(ntimes),l_sigbmls(ntimes), &
                l_rrmls(ntimes),l_zzmls(ntimes),l_L1mls(ntimes), &
                l_L2mls(ntimes),l_L4mls(ntimes),l_epotpmls(ntimes), &
                l_sigepmls(ntimes)

      l_time=0.0
      l_bbmls=0.0
      l_sigbmls=0.0
      l_rrmls=0.0
      l_zzmls=0.0
      l_L1mls=0.0
      l_L2mls=0.0
      l_L4mls=0.0
      l_epotpmls=0.0
      l_sigepmls=0.0
      l_iermls=0

      ! Read complete time histories from
      ! /u/grierson/efit/msels/[shot]/msels_chan[xx].dat
      write(*,*)
      write(*,101) 'Getting data for channel ',icmls
      call msels_hist(ishot,synmls,icmls,ntimes,l_time,l_bbmls,l_sigbmls, &
                      l_rrmls,l_zzmls,l_L1mls,l_L2mls,l_L4mls,l_epotpmls, &
                      l_sigepmls,l_iermls)
      
      ! For each k-file time average the data for this chord
      write(*,*) ' Filling ktime arrays'
      do i=1,ktime
         ! Initialize counter
         count=0
         ! Zero arrays
         bbmls(i) = 0.0
         sigbmls(i) = 0.0
         rrmls(i) = 0.0
         zzmls(i) = 0.0
         L1mls(i) = 0.0
         L2mls(i) = 0.0
         L4mls(i) = 0.0
         epotpmls(i) = 0.0
         sigepmls(i) = 0.0
         iermls(i)=1
         if (l_iermls .eq. 0) then 
            do j=1,ntimes
               if (l_time(j) .ge. atime(i)-avem .and. l_time(j) .le. atime(i)+avem) then
                  bbmls(i) = bbmls(i) + l_bbmls(j)
                  sigbmls(i) = sigbmls(i) + l_sigbmls(j)
                  rrmls(i) = rrmls(i) + l_rrmls(j)
                  zzmls(i) = zzmls(i) + l_zzmls(j)
                  L1mls(i) = L1mls(i) + l_L1mls(j)
                  L2mls(i) = L2mls(i) + l_L2mls(j)
                  L4mls(i) = L4mls(i) + l_L4mls(j)
                  epotpmls(i) = epotpmls(i) + l_epotpmls(j)
                  sigepmls(i) = sigepmls(i) + l_sigepmls(j)
                  
                  ! Increase couter for mean
                  count = count+1
               endif
            enddo
            if (count .gt. 0) then
               bbmls(i) = bbmls(i)/count
               sigbmls(i) = sigbmls(i)/count
               rrmls(i) = rrmls(i)/count
               zzmls(i) = zzmls(i)/count
               L1mls(i) = L1mls(i)/count
               L2mls(i) = L2mls(i)/count
               L4mls(i) = L4mls(i)/count
               epotpmls(i) = epotpmls(i)/count
               sigepmls(i) = sigepmls(i)/count
               iermls(i)=0
            endif
         endif
      enddo
      
101   format (a, i6)
102   format (a, g15.4)
      return
      end subroutine msels_data
