      module siloop
!module for siloop
!Revisions:
!$Log: siloop.f90,v $
!Revision 1.1  2007/06/01 05:16:38  renq
!*** empty log message ***
!
!
      use exparm,only:nsilop
      implicit none
      public

      real*8, dimension (:), allocatable :: rsi,zsi,wsi,hsi,&
                                as,as2

      end module siloop
