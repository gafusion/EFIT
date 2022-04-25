!
!  Module to control reading of thomson mtanh data from mds+ database
!
   MODULE mds_mtanh

!!!!!!   use hpalias
   implicit none

   private
   public :: mtanh_ts, Mtanh_Init, GetNtime, GetTime,            &
             GetAlp, GetOff, GetPed, GetSym, GetWid, GetHwd, Shot,       &
             Id
    
   ! Define a structure to hold the Thomson mtanh data
   type mtanh_ts
      private
    ! DO NOT INITIALIZE THESE GUYS - THIS IS an F95 FEATURE WHICH THE
    ! HP COMPILER DOES NOT RECOGNIZE
      character(len=4) :: id  ! ='NONE'
      integer :: shot         !=  -1
      integer :: num          != 0
      real, pointer, dimension(:) :: time
      real, pointer, dimension(:) :: alp
      real, pointer, dimension(:) :: off
      real, pointer, dimension(:) :: ped
      real, pointer, dimension(:) :: sym
      real, pointer, dimension(:) :: wid
   end type mtanh_ts

   interface mtanh_init
      module procedure mtanh_ts_init
   end interface 

   interface GetTime
      module procedure mtanh_ts_GetTime
   end interface 

   interface GetNtime
      module procedure mtanh_ts_GetNtime
   end interface 

   interface GetAlp
      module procedure mtanh_ts_GetAlp
   end interface 

   interface GetOff
      module procedure mtanh_ts_GetOff
   end interface 

   interface GetPed
      module procedure mtanh_ts_GetPed
   end interface 

   interface GetSym
      module procedure mtanh_ts_GetSym
   end interface 

   interface GetWid
      module procedure mtanh_ts_GetWid
   end interface 

   interface GetHwd
      module procedure mtanh_ts_GetHwd
   end interface 

   interface Shot
      module procedure mtanh_ts_Shot
   end interface 

   interface Id
      module procedure mtanh_ts_Id
   end interface 


   contains
     

      ! Subroutine to read tanhfit goodies for electrons from mdsplus database
      ! Written: February 22, 2001 by rjg
      SUBROUTINE Mtanh_Ts_Init (self, ishot, cdata, err) 
                                                                        
      use mds_rg
      use mdslib
      IMPLICIT NONE 

                                                                        
      ! Data passed via argument list                                 
      type (mtanh_ts), intent(inout) :: self ! Object with mtanh fit data
      integer, intent(in) :: ishot           ! Desired shot number
      character(len=*), intent(in) :: cdata  ! Desired data, i.e., peped
      logical, intent(out) :: err            ! TRUE means an error occurred
                                                                        

      integer :: ntime          
      integer*4 :: istat
      integer*4 :: isz
      character(len=100) :: mds_node
      character(len=100) :: node_arr
      real, allocatable, dimension(:) :: rarray
      integer :: ilen
      logical :: ok
      character(len=20) :: signal
      integer length 
      character(len=9) :: tree = 'ELECTRONS'
      character(len=min(2,len(cdata))) :: cid


      err = .false.

      ! Connect to mds
      ! Need to connect to atlas as of 6/7/99 for mds+ on unix            
      ! status = MdsConnect('atlas.gat.com'//CHAR(0)) 
      ! WAKE_MDS knows if we are already connected and will only connect
      ! if that has not been done
      call wake_mds ()

      ! Want to work with upper case only but avoid str$upcase
      cid = cdata
      if (cid(1:1) == 'n') then
          cid(1:1) = 'N'
      else if (cid(1:1) == 'p') then
          cid(1:1) = 'P'
      else if (cid(1:1) == 't') then
          cid(1:1) = 'T'
      endif

      if (cid(2:2) == 'e') cid(2:2) = 'E'

!      print *,' cid = ', cid

      ! Open ELECTRONS tree                                           
      if (cid(2:2) .eq. 'E') then
         istat = MdsOpen(trim(tree)//CHAR(0), ishot) 
         err = MdsError (istat, 'MdsOpen', tree)
         if (err) return
      else
         print *,' do not recognize ', cdata
         err = .true.
         return
      endif
                                                                        
      ! Get to the tanhfit node for electron density
      if (cid(1:2) == 'NE') then 
         mds_node =  "\TOP.PROFILE_FITS.TANHFIT.DENSITY"//CHAR(0)
      ! Get to the tanhfit node for electron temperature
      else if (cid(1:2) == 'TE') then 
         mds_node =  "\TOP.PROFILE_FITS.TANHFIT.TEMPERATURE"//CHAR(0)
      ! Get to the tanhfit node for electron pressure
      else if (cid(1:2) == 'PE') then 
         mds_node =  "\TOP.PROFILE_FITS.TANHFIT.PRESSURE"//CHAR(0)
      ! Error condition
      else 
         print *,' Do not recognize id = ', cid
         err= .true.
         return
      endif

      istat = MdssetDefault(trim(mds_node))
      err = MdsError (istat, 'MdssetDefault', mds_node)
      if (err) return

      ntime = 0
!      istat = MdsValue("SIZE(\TSTIME_CORE)"//CHAR(0),               &
      istat = MdsValue("SIZE(TIME)"//CHAR(0),               &
     &                   descr(IDTYPE_LONG,isz,0),0,length)      
                                                                        
      ! Make sure to check istat from return
      if (good_status(istat)) then 
         ntime = isz
         if (allocated(rarray)) deallocate(rarray,stat=istat)
         allocate(rarray(ntime))
         call construct_mtanh_ts (self, cdata, ishot, ntime)         
      else
         ntime = 0
      endif

      ! Load array data
      signal = 'TIME'  
      istat  = MdsValue(trim(signal)//CHAR(0),                             &
     &                descr(IDTYPE_FLOAT,rarray,ntime,0),0,length)
      err = MdsError (istat, 'MdsValue', signal)
      if (err) return
      self%time = rarray

      signal = 'ALP'  
      istat  = MdsValue(trim(signal)//CHAR(0),                             &
     &                descr(IDTYPE_FLOAT,rarray,ntime,0),0,length)
      err = MdsError (istat, 'MdsValue', signal)
      if (err) return
      self%alp = rarray 

      signal = 'OFF'  
      istat  = MdsValue(trim(signal)//CHAR(0),                             &
     &                descr(IDTYPE_FLOAT,rarray,ntime,0),0,length)
      err = MdsError (istat, 'MdsValue', signal)
      if (err) return
      self%off = rarray

      signal = 'PED'  
      istat  = MdsValue(trim(signal)//CHAR(0),                             &
     &                descr(IDTYPE_FLOAT,rarray,ntime,0),0,length)
      err = MdsError (istat, 'MdsValue', signal)
      if (err) return
      self%ped = rarray

      signal = 'SYM'  
      istat  = MdsValue(trim(signal)//CHAR(0),                             &
     &                descr(IDTYPE_FLOAT,rarray,ntime,0),0,length)
      err = MdsError (istat, 'MdsValue', signal)
      if (err) return
      self%sym = rarray

      signal = 'WID'  
      istat  = MdsValue(trim(signal)//CHAR(0),                             &
     &                descr(IDTYPE_FLOAT,rarray,ntime,0),0,length)
      err = MdsError (istat, 'MdsValue', signal)
      if (err) return
      self%wid = rarray

      ! Close this tree from mds
      ! Note that we do not completely disconnect here.  It is possible
      ! that other parts of the code will run a bit faster if we remain
      ! connected.
      istat = MdsClose(trim(tree)// Char(0), ishot)
      err = MdsError (istat, 'MdsClose', tree)

      end subroutine Mtanh_Ts_Init


      ! Routine to initialize an mtanh_ts object to size ntime
      ! Written: 2/26/01 by rjg
      subroutine construct_mtanh_ts (self, cdata, ishot, ntime)         
      implicit none
      type (mtanh_ts), intent(out) :: self
      character(len=*), intent(in) :: cdata  ! Data type
      integer, intent(in) :: ishot          ! Shot number
      integer, intent(in) :: ntime         ! # of elements in each array
      integer :: istat
      
      self%id = cdata
      self%shot = ishot
      self%num = ntime


      ! Rather than going ahead and deallocating, I think I can first
      ! test to see if the pointer is associated.  If it isn't, it should
      ! not be allocated?

      ! Do we need to do any testing or deallocating?
      deallocate (self%time, stat=istat)
!      if (istat .ne. 0) print *,' error deallocating time'
      allocate (self%time(ntime))

      deallocate (self%alp, stat=istat)
!      if (istat .ne. 0) print *,' error deallocating alp '
      allocate (self%alp(ntime))

      deallocate (self%off, stat=istat)
!      if (istat .ne. 0) print *,' error deallocating off'
      allocate (self%off(ntime))

      deallocate (self%sym, stat=istat)
!      if (istat .ne. 0) print *,' error deallocating sym '
      allocate (self%sym(ntime))

      deallocate (self%ped, stat=istat)
!      if (istat .ne. 0) print *,' error deallocating ped '
      allocate (self%ped(ntime))

      deallocate (self%wid, stat=istat)
!      if (istat .ne. 0) print *,' error deallocating wid '
      allocate (self%wid(ntime))

 
      end subroutine construct_mtanh_ts
 
 
      ! Subroutine to extract number of times from an mtanh_ts object
      ! Written: 2/26/01 by rjg
      Subroutine Mtanh_Ts_GetNtime (self, num)
      implicit none
      type (mtanh_ts), intent(in) :: self
      integer, intent(out) :: num
           
      num = self%num

      End Subroutine Mtanh_Ts_GetNtime 


      ! Function to return shot number from an mtanh_ts object
      ! Written: 3/01/01 by rjg
      Integer function Mtanh_Ts_Shot (self)
      implicit none
      type (mtanh_ts), intent(in) :: self
           
      Mtanh_Ts_Shot = self%shot

      End function Mtanh_Ts_Shot


      ! Function to extract id from an mtanh_ts object
      ! Written: 3/01/01 by rjg
      Character(len=4) function Mtanh_Ts_Id (self)
      implicit none
      type (mtanh_ts), intent(in) :: self
      Mtanh_Ts_Id = self%id

      End function Mtanh_Ts_Id

 
      ! Subroutine to extract time from an mtanh_ts object
      ! Written: 2/26/01 by rjg
      Subroutine Mtanh_Ts_GetTime (self, time)
      implicit none
      type (mtanh_ts), intent(in) :: self
      real, pointer, dimension(:) :: time
      integer :: istat
           
      deallocate (time, stat=istat)
      allocate (time(size(self%time)))
      time = self%time

      End Subroutine Mtanh_Ts_GetTime 

 
      ! Subroutine to extract Alp from an mtanh_ts object
      ! Written: 2/26/01 by rjg
      Subroutine Mtanh_Ts_GetAlp (self, Alp)
      implicit none
      type (mtanh_ts), intent(in) :: self
      real, pointer, dimension(:) :: Alp
      integer :: istat
           
      deallocate (Alp, stat=istat)
      allocate (Alp(size(self%Alp)))
      Alp = self%Alp

      End Subroutine Mtanh_Ts_GetAlp 
 

      ! Subroutine to extract Off from an mtanh_ts object
      ! Written: 2/26/01 by rjg
      Subroutine Mtanh_Ts_GetOff (self, Off)
      implicit none
      type (mtanh_ts), intent(in) :: self
      real, pointer, dimension(:) :: Off
      integer :: istat
           
      deallocate (Off, stat=istat)
      allocate (Off(size(self%Off)))
      Off = self%Off

      End Subroutine Mtanh_Ts_GetOff 
 
      ! Subroutine to extract Ped from an mtanh_ts object
      ! Written: 2/26/01 by rjg
      Subroutine Mtanh_Ts_GetPed (self, Ped)
      implicit none
      type (mtanh_ts), intent(in) :: self
      real, pointer, dimension(:) :: Ped
      integer :: istat
           
      deallocate (Ped, stat=istat)
      allocate (Ped(size(self%Ped)))
      Ped = self%Ped

      End Subroutine Mtanh_Ts_GetPed 
 

      ! Subroutine to extract Sym from an mtanh_ts object
      ! Written: 2/26/01 by rjg
      Subroutine Mtanh_Ts_GetSym (self, Sym)
      implicit none
      type (mtanh_ts), intent(in) :: self
      real, pointer, dimension(:) :: Sym
      integer :: istat
           
      deallocate (Sym, stat=istat)
      allocate (Sym(size(self%Sym)))
      Sym = self%Sym

      End Subroutine Mtanh_Ts_GetSym 
 

      ! Subroutine to extract Wid from an mtanh_ts object
      ! Written: 2/26/01 by rjg
      Subroutine Mtanh_Ts_GetWid (self, Wid)
      implicit none
      type (mtanh_ts), intent(in) :: self
      real, pointer, dimension(:) :: Wid
      integer :: istat
           
      deallocate (Wid, stat=istat)
      allocate (Wid(size(self%Wid)))
      Wid = self%Wid

      End Subroutine Mtanh_Ts_GetWid 


      ! Subroutine to extract Half-wid from an mtanh_ts object
      ! Written: 4/30/01 by rjg
      Subroutine Mtanh_Ts_GetHwd (self, Hwd)
      implicit none
      type (mtanh_ts), intent(in) :: self
      real, pointer, dimension(:) :: Hwd
      integer :: istat
           
      deallocate (Hwd, stat=istat)
      allocate (Hwd(size(self%Wid)))
      Hwd = 0.5*self%Wid

      End Subroutine Mtanh_Ts_GetHwd 
       

   end module mds_mtanh
!
!
