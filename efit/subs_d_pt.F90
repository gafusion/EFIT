!**********************************************************************
!>
!!    ---------------------------------------------------------------------
!!    --   Subroutines for non-portable functions
!!    Most have been eliminated, but if there are any compiler-dependent
!!    functions they should be done via ifdef
!!    ---------------------------------------------------------------------
!!    SEK:  Not sure why we should keep these at all.
!!    Since they do nothing, why not remove?  need to understand what
!!    DIII-D EFIT is doing with these.
!!    TODO
!!    ----------------------------------------------------------------------
!!
!!    @param ishot :
!!
!!    @param itime :
!!
!!    @param header :
!!
!**********************************************************************
     subroutine db_header(ishot,itime,header)
     character*(*) header
     header = ' '
     return
     end subroutine db_header
     subroutine getzeff()
     return
     end subroutine getzeff
     subroutine donepl()
     return
     end subroutine donepl


!**********************************************************************
!>
!!    This is to replace the lib$wait calls
!!    See: https://stackoverflow.com/questions/6931846/sleep-in-fortran
!!    
!!
!**********************************************************************
module Fortran_Sleep

   use, intrinsic :: iso_c_binding, only: c_int

   implicit none

   interface

      !  should be unsigned int ... not available in Fortran
      !  OK until highest bit gets set.
      function FortSleep (seconds)  bind ( C, name="sleep" )
          import
          integer (c_int) :: FortSleep
          integer (c_int), intent (in), VALUE :: seconds
      end function FortSleep

   end interface

end module Fortran_Sleep
