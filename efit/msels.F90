#include "config.f"
#if defined(USE_SNAP)
!**********************************************************************
!>
!!    getmsels obtains MSE-LS data
!!    
!!
!!    @param ktime : number of time slices
!!
!**********************************************************************
      subroutine getmsels(ktime)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      character*3 synmlt
      integer*4 ktime,icmls,iermls(ntime)
      real*8 avemlt,atime(ntime),bbmls(ntime), &
             rrmls(ntime),zzmls(ntime), &
             l1mls(ntime),l2mls(ntime),l4mls(ntime), &
             epotpmls(ntime),sigepmls(ntime)
!
      avemlt=avemsels
      synmlt=synmsels
      atime=time*1000.0
      do i=1,nmsels
        icmls=i
        call msels_data(ishot,atime,ktime,avemlt,synmlt,icmls, &
             bbmls,rrmls,zzmls,l1mls,l2mls,l4mls,epotpmls, &
             sigepmls,iermls)
#ifdef DEBUG_LEVEL2
        write (6,*) 'GETMSELS bbmls= ',bbmls(1)
#endif
        bmselt(1:ktime,i)=bbmls(1:ktime)
        sbmselt(1:ktime,i)=0.0 ! unset
        rrmselt(1:ktime,i)=rrmls(1:ktime)
        zzmselt(1:ktime,i)=zzmls(1:ktime)
        l1mselt(1:ktime,i)=l1mls(1:ktime)
        l2mselt(1:ktime,i)=l2mls(1:ktime)
        l3mselt(1:ktime,i)=1.0-l1mls(1:ktime)
        l4mselt(1:ktime,i)=l4mls(1:ktime)
        emselt(1:ktime,i)=epotpmls(1:ktime)
        semselt(1:ktime,i)=sigepmls(1:ktime)
        iermselt(1:ktime,i)=iermls(1:ktime)
        fwtbmselt(1:ktime,i)=swtbmsels(1:ktime)
        fwtemselt(1:ktime,i)=swtemsels(1:ktime)
        do j=1,ktime
          if (iermselt(j,i).ne.0) then
            fwtbmselt(j,i)=0.0
            fwtemselt(j,i)=0.0
          endif 
        enddo
      enddo
!
      return
      end subroutine getmsels

!**********************************************************************
!>
!!    msels_data returns synthetic or experimental MSE-LS
!!    data to EFIT
!!    For list of k-file times returs the available MSE-LS
!!    data for a given chord, averaged over avgtim.
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
             bbmls,rrmls,zzmls,L1mls,L2mls,L4mls,epotpmls, &
             sigepmls,iermls) 
      integer*4, parameter :: ntimes=2000
      character*3 synmls
      integer*4 i,j,count
      integer*4 ishot,ktime,icmls,iermls(ktime),l_iermls
      real*8 avem,atime(ktime),bbmls(ktime), &
             rrmls(ktime),zzmls(ktime), &
             L1mls(ktime),L2mls(ktime),L4mls(ktime), &
             epotpmls(ktime),sigepmls(ktime)

      ! local variables
      character(len=100) filename
      integer*4 file_shot
      real*8 :: l_time(ntimes),l_bbmls(ntimes), &
                l_rrmls(ntimes),l_zzmls(ntimes),l_L1mls(ntimes), &
                l_L2mls(ntimes),l_L4mls(ntimes),l_epotpmls(ntimes), &
                l_sigepmls(ntimes)

      l_time=0.0
      l_bbmls=0.0
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
#ifdef DEBUG_LEVEL1
      write(*,*)
      write(*,101) 'Getting data for channel ',icmls
#endif
      call msels_hist(ishot,synmls,icmls,ntimes,l_time,l_bbmls, &
                      l_rrmls,l_zzmls,l_L1mls,l_L2mls,l_L4mls,l_epotpmls, &
                      l_sigepmls,l_iermls)
      
      ! For each k-file time average the data for this chord
#ifdef DEBUG_LEVEL1
      write(*,*) ' Filling ktime arrays'
#endif
      do i=1,ktime
         ! Initialize counter
         count=0
         ! Zero arrays
         bbmls(i) = 0.0
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

!**********************************************************************
!>
!!    msels_hist returns synthetic or experimental MSE-LS
!!    data to EFIT
!!    For a given shot, returns complete time history
!!
!!    WARNING: this subroutine uses both REAL*4 (used in mse files) and
!!             REAL*8 variables conversions must be handled carefully
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
!!    @param time : time array
!!
!!    @param bbmls : MSE-LS magnetic field E/V_BEAM in Tesla 
!!
!!    @param rrmls : R of MSE-LS channels in m  
!!
!!    @param zzmls : Z of MSE-LS channels in m  
!!
!!    @param L1mls :  geometric cofficients of MSE-LS channels 
!!
!!    @param L2mls :  geometric cofficients of MSE-LS channels 
!!
!!    @param L4mls :  geometric cofficients of MSE-LS channels 
!!
!!    @param epotpmls : flux derivative of electrostatic potential
!!                    in 1/sec
!!
!!    @param sigepmls : uncertainty in electrostatic potential flux
!!                   derivative in 1/sec 
!!
!!    @param iermls :error flags of MSE-LS channels, 0=normal,  
!!                 1=error (ignore data)
!!
!**********************************************************************
      subroutine msels_hist(ishot,synmls,icmls,ntimes,time, &
                            bbmls,rrmls,zzmls,L1mls,L2mls, &
                            L4mls,epotpmls,sigepmls,iermls)

      use set_kinds
      integer*4 ishot,icmls,iermls,i
      real*8 time(ntimes),bbmls(ntimes), &
             rrmls(ntimes),zzmls(ntimes), &
             L1mls(ntimes),L2mls(ntimes),L4mls(ntimes), &
             epotpmls(ntimes),sigepmls(ntimes)
      character*3 synmls
      
      ! Local variables
      logical file_exists
      character(len=100) filename
      integer*4 ioerr,file_shot,ntimes,issyn,l_ntimes
      real*4 l_time(ntimes),l_bbmls(ntimes), &
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

#ifdef DEBUG_LEVEL1
      write(*,*) "Using file ",filename
#endif

      ! See if file exists.  If not, return error code.
      inquire(file=filename,exist=file_exists)
      if (.not. file_exists) then
         write(*,*) synmls,' file does not exist'
         iermls=1
         return
      endif

      ! Read the data
      open(unit=1,file=filename,iostat=ioerr)
#ifdef DEBUG_LEVEL1
      write(*,101) 'iostat: ',ioerr
#endif
      ! Read the shot number in the file
      read(1,*) file_shot
#ifdef DEBUG_LEVEL1
      write(*,101) 'file shot:', file_shot
#endif

      ! Read the number of timeslices
      read(1,*) l_ntimes
#ifdef DEBUG_LEVEL1
      write(*,101) '# slices:', l_ntimes
#endif

      ! Read if synthetic(1) or measured(0)
      read(1,*) issyn
#ifdef DEBUG_LEVEL1
      write(*,101) 'issyn: ', issyn
#endif

      ! Read the channel number (1,2,etc...)
      read(1,*) icmls
#ifdef DEBUG_LEVEL1
      write(*,101) 'icmls: ', icmls
#endif

      ! Read the error code
      read(1,*) iermls
#ifdef DEBUG_LEVEL1
      write(*,101) 'iermls: ', iermls
#endif

#ifdef DEBUG_LEVEL1
      write(*,100) 'Starting time loop'
#endif
      do i=1,l_ntimes

         read(1,*) l_time(i)
#ifdef DEBUG_LEVEL3
         write(*,102) 'Time:', l_time(i)
#endif

         read(1,*) l_bbmls(i)
#ifdef DEBUG_LEVEL3
         write(*,102) 'B_LS:', l_bbmls(i)
#endif
         
         read(1,*) l_rrmls(i)
#ifdef DEBUG_LEVEL3
         write(*,102) 'sig B_LS:', l_rrmls(i)
#endif
         
         read(1,*) l_zzmls(i)
#ifdef DEBUG_LEVEL3
         write(*,102) 'Z:', l_zzmls(i)
#endif
         
         read(1,*) l_L1mls(i)
#ifdef DEBUG_LEVEL3
         write(*,102) 'L1:', l_L1mls(i)
#endif
         
         read(1,*) l_L2mls(i)
#ifdef DEBUG_LEVEL3
         write(*,102) 'L2:', l_L2mls(i)
#endif
         
         read(1,*) l_L4mls(i)
#ifdef DEBUG_LEVEL3
         write(*,102) 'L4:', l_L4mls(i)
#endif
         
         read(1,*) l_epotpmls(i)
#ifdef DEBUG_LEVEL3
         write(*,102) 'epot:', l_epotpmls(i)
#endif
         
         read(1,*) l_sigepmls(i)
#ifdef DEBUG_LEVEL3
         write(*,102) 'sigep:', l_sigepmls(i)
#endif
         
      enddo
      close(1)

      time = real(l_time,dp)
      bbmls = real(l_bbmls,dp)
      rrmls = real(l_rrmls,dp)
      zzmls = real(l_zzmls,dp)
      L1mls = real(l_L1mls,dp)
      L2mls = real(l_L2mls,dp)
      L4mls = real(l_L4mls,dp)
      epotpmls = real(l_epotpmls,dp)
      sigepmls = real(l_sigepmls,dp)

100   format (a)
101   format (a, i6)
102   format (a, g15.4)
      return

      end subroutine msels_hist
#endif

!**********************************************************************
!>
!!    erpote computes the stream function for the
!!    radial electric field.
!!    
!!
!!    @param ypsi :
!!
!!    @param nnn :
!!
!**********************************************************************
      real*8 function erpote(ypsi,nnn)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit none
      real*8 erppote
      integer*4, intent(in) :: nnn
      real*8, intent(in) :: ypsi
      real*8 xpsii(nercur)
!
      if (abs(ypsi).gt.1.0) then
        erpote=0.0
        return
      endif
      call seter(ypsi,xpsii)
      erpote=sum(cerer(1:nnn)*xpsii(1:nnn))
      return
!
      entry erppote(ypsi,nnn)
      if (abs(ypsi).gt.1.0) then
        erppote=0.0
        return
      endif
      call seterp(ypsi,xpsii)
      erppote=-sum(cerer(1:nnn)*xpsii(1:nnn))/sidif
      return
      end function erpote

!**********************************************************************
!>
!!    eradial computes the radial electric field.
!!    
!!
!!    @param ypsi :
!!
!!    @param nnn :
!!
!!    @param reee :
!!
!!    @param zeee :
!!
!**********************************************************************
      real*8 function eradial(ypsi,nnn,reee,zeee)
      use commonblocks,only: c,bkx,bky
      include 'eparm.inc'
      include 'modules1.inc'
      implicit none
      real*8 erpote,esradial,erppote,seval
      integer*4, intent(in) :: nnn
      real*8, intent(in) :: ypsi,reee,zeee
      integer*4 ier
      real*8 pds(6),fnow,bbnow
!
      if (abs(ypsi).ge.1.0) then
        eradial=0.0
        return
      endif
      call seva2d(bkx,lkx,bky,lky,c,reee,zeee,pds,ier,n333)
      eradial=erpote(ypsi,nnn)
      eradial=-pds(2)*eradial
      return
!
      entry esradial(ypsi,nnn,reee,zeee)
      if (abs(ypsi).ge.1.0) then
        esradial=0.0
        return
      endif
      call seva2d(bkx,lkx,bky,lky,c,reee,zeee,pds,ier,n333)
      esradial=erppote(ypsi,nnn)
      esradial=-(reee*pds(2))**2*esradial
      fnow=seval(nw,ypsi,sigrid,fpol,bbfpol,ccfpol,ddfpol)
      bbnow=sqrt(fnow**2+pds(2)**2+pds(3)**2)/reee
      esradial=esradial/bbnow
      return
      end function eradial
