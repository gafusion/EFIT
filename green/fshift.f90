      module fshift
!module for fshift
!Revisions:
!$Log: fshift.f90,v $
!Revision 1.1  2007/06/01 05:16:20  renq
!*** empty log message ***
!
!
      use exparm,only:nfcoil,magpr2
      implicit none
      public

      integer*4, dimension (:), allocatable :: nshiftrz
      real*8, dimension (:), allocatable :: rshift,zshift,pshift
      real*8, dimension (:), allocatable :: pmprobe

      end module fshift
