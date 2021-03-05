      module fcoil
!module for fcoil
!Revisions:
!$Log: fcoil.f90,v $
!Revision 1.1  2007/06/01 05:16:17  renq
!*** empty log message ***
!
!
      use exparm,only:nfcoil
      implicit none
      public

      real*8, dimension (:), allocatable :: rf,zf,wf,hf,af,af2,turnfc,&
                                fcid,fcturn

      end module fcoil
