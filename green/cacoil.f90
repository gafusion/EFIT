      module cacoil
!module for cacoil
!Revisions:
!$Log: cacoil.f90,v $
!Revision 1.1  2007/06/01 05:15:54  renq
!*** empty log message ***
!
!
      use exparm,only:nacoil
      implicit none
      public

      real*8, dimension (:), allocatable :: racoil,zacoil,wacoil,hacoil
      integer*4 :: iacoil

      end module cacoil
