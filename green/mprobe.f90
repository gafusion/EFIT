      module mprobe
!module for mprobe
!Revisions:
!$Log: mprobe.f90,v $
!Revision 1.1  2007/06/01 05:16:26  renq
!*** empty log message ***
!
!
      use exparm,only:magpr2
      implicit none
      public

      real*8,dimension (:), allocatable :: xmp2,ymp2,amp2,smp2
      integer*4 :: nsmp2

      end module mprobe
