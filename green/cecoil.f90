      module cecoil
!module for cecoil
!Revisions:
!$Log: cecoil.f90,v $
!Revision 1.1  2007/06/01 05:15:57  renq
!*** empty log message ***
!
!
      use exparm,only:necoil
      implicit none
      public

      real*8, dimension (:), allocatable :: re,ze,we,he,ecid,ecturn

      end module cecoil
