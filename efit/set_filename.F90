#include "config.f"
!**********************************************************************
!>
!!    This subroutine sets the file names for all standard output files
!!    
!!    @param prefix : filename prefix (letter or name)
!!    @param ishot : shot number
!!    @param itime : time in milliseconds
!!    @param itimeu : time in microseconds 
!!    @param suffix : filename suffix
!!    @param fname : file name 
!!
!**********************************************************************
      subroutine setfnm(prefix,ishot,itime,itimeu,suffix,fname)
      use var_nio, only: nttyo
      use var_input, only: istore
      use extvars, only: store_dir
      implicit none
      integer*4, intent(in) :: itimeu,ishot,itime
      character*(*), intent(in) :: prefix,suffix
      character*(*), intent(inout) :: fname
      integer*4 i,iitime,fnlen,length
      character*4 timeu
      character*9 shot,time

      ! Warn if a filename length risks overflow (store_dir is 256 chars)
      fnlen=len(fname)
      if(fnlen<300) &
        write(nttyo,'(A)') &
          'Warning: The filename length is short for: '//prefix

      ! Warn if the prefix and suffix still risk overflowing the filename
      length=len(prefix)+len(suffix)
      if(length>20) &
        write(nttyo,'(A)') &
          'Warning: The prefix and suffix are too long for: ' &
            //prefix//"+"//suffix

      ! Initialize empty name
      do i=1,fnlen
        fname(i:i)=' '
      enddo

      ! Handle negative times (unclear why this would be wanted...)
      iitime=itime
      if(iitime.lt.0) iitime=90000-iitime
#ifdef DEBUG_LEVEL2
      write (6,*) &
        'setfnm/prefix/ishot/iitime/itimeu/istore = ', &
        prefix,ishot,iitime,itimeu,istore
#endif

      ! Convert the shot number to a string
      if (ishot<1e6) then
        write(shot,600) ishot
      elseif (ishot<1e7) then
        write(shot,700) ishot
      elseif (ishot<1e8) then
        write(shot,800) ishot
      elseif (ishot<1e9) then
        write(shot,900) ishot
      else
        write(nttyo,*) &
          'The shot number is too large for existing file formats: ',ishot
        stop
      endif

      ! Convert times to a string
      if (iitime<1e5) then
        write(time,500) iitime
      elseif (iitime<1e6) then
        write(time,600) iitime
      elseif (iitime<1e7) then
        write(time,700) iitime
      elseif (iitime<1e8) then
        write(time,800) iitime
      elseif (iitime<1e9) then
        write(time,900) iitime
      else
        write(nttyo,*) &
          'The run time is too large for existing file formats: ',iitime
        stop
      endif

      ! Convert the micro sec time to a string if non-zero
      timeu='    '
      if (itimeu.gt.0) then
        write(timeu,300) itimeu
        timeu='_'//trim(timeu)
      endif

      ! Put the filename together
      fname = prefix//trim(shot)//'.'//trim(time)//trim(timeu)//suffix

      ! If istore = 1 then collect EFIT results in store_dir directory
      if (istore .eq. 1) then
        fname = trim(store_dir)//trim(fname)
      endif
#ifdef DEBUG_LEVEL2
      write (6,*) &
        'setfnm/shot/time/timeu/istore/fname = ', &
        shot,time,timeu,store_dir,fname
#endif
  300 format(i3.3)
  500 format(i5.5)
  600 format(i6.6)
  700 format(i7.7)
  800 format(i8.8)
  900 format(i9.9)
      return
      end subroutine setfnm

!**********************************************************************
!>
!!    This subroutine sets the omas output file name
!!    
!!    @param ishot : shot number
!!
!!    @param fname : file name 
!!
!**********************************************************************
      subroutine setfnmomas(ishot,fname)
      implicit none
      integer*4, intent(in) :: ishot
      character*(*), intent(inout) ::  fname
      integer*4 i,length

      ! Warn if a filename length risks overflow (store_dir is 256 chars)
      fnlen=len(fname)
      if(fnlen<275) &
        write(nttyo,'(A)') &
          'Warning: The filename length is short for: '//prefix

      ! Initialize empty name
      do i=1,fnlen
        fname(i:i)=' '
      enddo

      ! Convert the shot number to filename
      if (ishot<1e6) then
        write(shot,1600) ishot
      elseif (ishot<1e7) then
        write(shot,1700) ishot
      elseif (ishot<1e8) then
        write(shot,1800) ishot
      elseif (ishot<1e9) then
        write(shot,1900) ishot
      else
        write(nttyo,*) &
          'The shot number is too large for existing file formats: ',ishot
        stop
      endif

      ! If istore = 1 then put file in store_dir directory
      if (istore .eq. 1) then
        fname = trim(store_dir)//trim(fname)
      endif
 1600 format(i6.6,'.h5')
 1700 format(i7.7,'.h5')
 1800 format(i8.8,'.h5')
 1900 format(i9.9,'.h5')
      return
      end subroutine setfnmomas
