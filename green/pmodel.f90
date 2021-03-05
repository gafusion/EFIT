      module pmodel
!module for pmodel
!Revisions:
!$Log: pmodel.f90,v $
!Revision 1.1  2007/06/01 05:16:32  renq
!*** empty log message ***
!
!
      implicit none
      public

      real*8,dimension(:),allocatable,save :: rgrid
      real*8,dimension(:),allocatable,save :: zgrid
      real*8 :: dr,dz

      end module pmodel
