      module nio

      implicit none

      integer*4, parameter :: nin=11,nout=10,ntty=5
      integer*4, parameter :: nrsppc=25,nrspfc=26,ncontr=35
      character*4 :: ch1, ch2

      contains
!**********************************************************************
!>
!!    This subroutine set the defines file name suffixes to look for
!!    (e.g. if running nw=129 ,nw=129, efit will look for Green's function
!!     tables with ec129129, re129129, etc.)
!!    
!!
!!    @param nw : number of grid point (width)
!!
!!    @param nh : number of grid point (height)
!!
!!    @param ch1 : character of number of grid point (width)
!!
!!    @param ch2 : character of number of grid point (height)
!!
!**********************************************************************

      subroutine inp_file_ch(nw,nh,ch1,ch2)
        
        integer*4,   intent(in)  :: nw, nh 
        character*4, intent(out) :: ch1, ch2
        
        write(ch1,'(i4)') nw
        ch1 = adjustl(ch1)
        write(ch2,'(i4)') nh
        ch2 = adjustl(ch2)
        
        return
        
      end subroutine inp_file_ch

      end module nio
