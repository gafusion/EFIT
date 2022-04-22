!
!  Module to provide interface between efit and mtanh_mds module
!
   MODULE mds_mtanh_rdts
   
   use mds_mtanh
   implicit none

   private

   public :: mds_mtanh_ts_ld, time_range_ok, interp_mtanh, make_time_index 

   ! Declare module data
   type (mtanh_ts), save :: cached_mtanh_ts
   real, pointer, save, dimension(:) :: cached_time
   real, pointer, save, dimension(:) :: cached_wid
   real, pointer, save, dimension(:) :: cached_sym
   real, pointer, save, dimension(:) :: cached_ped
   integer, dimension(:,:), save, allocatable :: index 
   character(len=3), public, parameter :: cped = 'PED'
   character(len=3), public, parameter :: csym = 'SYM'
   character(len=3), public, parameter :: cwid = 'WID'

      contains

      ! Routine to make sure that the cached_mtanh_ts object has the proper data
      logical function mds_mtanh_ts_ld (ishot, cdata)
      integer, intent(in) :: ishot   ! Shot for which we want mtanh data
      character(len=*), intent(in) :: cdata ! Type of data for wc we want 
                                            ! mtanh data

      ! See if data are already loaded
      if (id(cached_mtanh_ts) == cdata .and.                     &
          shot(cached_mtanh_ts) == ishot)                       then
         mds_mtanh_ts_ld = .true.
      ! Else load the data
      else
         call Mtanh_Init (cached_mtanh_ts, ishot, cdata, mds_mtanh_ts_ld)
         ! If we got the object, load up the pointers with the goodies
         if (.not. mds_mtanh_ts_ld) then
            call GetTime (cached_mtanh_ts, cached_time)
            call GetSym  (cached_mtanh_ts, cached_sym)
            call GetWid  (cached_mtanh_ts, cached_wid)
            call GetPed  (cached_mtanh_ts, cached_ped)
         else
            print *,' Error loading mtanh data for shot ', ishot
         endif
      endif

      end function mds_mtanh_ts_ld


      ! Routine to make sure that all times in time array are encompassed
      ! by times in p_time array
      ! Written: 3/2/01 by rjg
      logical function time_range_ok (ntime, time)
      integer, intent(in) :: ntime
      real, dimension(ntime), intent(in) :: time
      if (time(1) .lt. cached_time(1)) then
         time_range_ok = .false.
         print *,' First requested time is less than 1st mtanh time ',   &
                   cached_time(1)
      else if (time(ntime) .gt. cached_time(size(cached_time)) ) then
         time_range_ok = .false.
         print *,' Last requested time is greater than last mtanh time ',   &
                   cached_time(size(cached_time))
      else
         time_range_ok = .true.
      endif
      end function time_range_ok



      ! Function to compute indices of locations in cached_time array to 
      ! use for interpolation
      ! We assume that time and cached_time arrays are both ordered in
      ! increasing value. 
      ! We assume that extremes of cached_time overlap extremes of time array
      ! Written: 3/2/01 by rjg
      subroutine make_time_index (ntime, time)
      implicit none
      integer, intent(in) :: ntime
      real, dimension(ntime), intent(in) :: time
      integer :: it    ! Index in time array
      integer :: ict   ! Index in cached_time array
      integer :: ict_strt   ! Location in cached_time array to start search
      integer :: isz   ! Size of cached_time array

      if (allocated(index)) deallocate(index)
      allocate (index(2,ntime))
  
      ict = 1
      isz = size (cached_time)
      loop_out: do it = 1, ntime
         ict_strt = ict
         loop_in: do ict = ict_strt , isz-1
            if (cached_time(ict) .le. time(it) .and.                 &
                time(it) .le. cached_time(ict+1))          then
                index(1,it) = ict
                index(2,it) = ict+1
!        write (6,fmt='(a,i5)') ' index(1,it) = ', index(1,it)
!        write (6,fmt='(a,i5)') ' index(2,it) = ', index(2,it)
                cycle loop_out
            endif
         enddo loop_in
      enddo loop_out

      end subroutine make_time_index


      
      ! Function to return desired values interpolated onto time array.
      ! Written: 3/2/01 by rjg
      ! Revised: 3/23/04 by rjg to fix syntax error
      function interp_mtanh (ntime,time,goodie) result(mtanh_out)
!      function interp_mtanh (ntime,time,goodie) 
      implicit none
      integer, intent(in) :: ntime
      real, dimension(ntime), intent(in) :: time
      character(len=*) :: goodie
      real, dimension(ntime) :: mtanh_out
!      real, dimension(ntime) :: interp_mtanh
      real, dimension(:), pointer :: p_goodie
      real :: delta
      integer :: it

      ! Find out where to point to
      if (associated(p_goodie)) nullify(p_goodie)
      if (goodie == cped) then
         p_goodie => cached_ped
      else if (goodie == csym) then
         p_goodie => cached_sym
      else if (goodie == cwid) then
         p_goodie => cached_wid
      endif

      ! Finally, we do the interpolation
      if (.not. associated(p_goodie)) then
         print *,' Coding error: do not recognize ', goodie
         print *,' Returning zeros for interpolation'
         mtanh_out = 0
!         interp_mtanh = 0
         return
      endif          
           
      do it = 1, ntime
         delta = (p_goodie(index(2,it)) - p_goodie(index(1,it)))  / &
                 (cached_time(index(2,it)) - cached_time(index(1,it))) 
!         interp_mtanh(it) = p_goodie(index(1,it)) +                 &
         mtanh_out(it) = p_goodie(index(1,it)) +                 &
                 delta * (time(it) - cached_time(index(1,it)))
      enddo 
   
      end function interp_mtanh


   end module mds_mtanh_rdts
     

      ! Subroutine to interface between efit and mtanh_mds module
      ! Written: March 01, 2001 by rjg
      ! Revised: March 20, 2001 by rjg to make stat an array
      subroutine get_mtanh_ts (ishot, cdata, ntime, time, sym, wid, ped, stat) 
      use mds_mtanh_rdts
      
      IMPLICIT NONE 

      integer, intent(in) :: ishot           ! Desired shot number
      character(len=*), intent(in) :: cdata  ! Desired data, i.e., 
                                             ! 'NE', 'TE' or 'PE'
                                                  ! which data are desired
      integer, intent(in) :: ntime           ! Number of times
      real, dimension(ntime), intent(in) :: time  ! Array of times (ms) for 
      real, dimension(ntime), intent(out) :: sym  ! Symmetry points (m) 
                                                  ! evaluated on time array
      real, dimension(ntime), intent(out) :: wid  ! Full widths (m) 
                                                  ! evaluated on time array
      real, dimension(ntime), intent(out) :: ped  ! Pedestal values (units 
                                                  ! vary on time array
      logical, dimension(ntime), intent(out) :: stat  ! True means no error for
                                                      ! the corresponding time
      logical :: stat_loc

      ! Initialize the status array for return 
      stat = .false.
      ! Initialize other arrays
      sym = 0
      wid = 0
      ped = 0

      ! Ask the efit_mtanh module to load the proper mtanh data.
      stat_loc = .not. mds_mtanh_ts_ld (ishot, cdata)
      if (stat_loc) then
         ! Make sure that the mtahn timebase covers the requested time range.
         stat_loc = time_range_ok (ntime, time)
         if (.not. stat_loc) then
            print *,' MTANH data do not cover requested time range '
            return
         endif
         ! Load indices array for interpolation
         call make_time_index (ntime, time)
         ped = interp_mtanh (ntime,time,cped)
         sym = interp_mtanh (ntime,time,csym)
         wid = interp_mtanh (ntime,time,cwid)
         stat = .true.
      endif

 

      end subroutine get_mtanh_ts


!       
!
!
