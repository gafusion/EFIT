      module bfgrid
!module for bfgrid
!Revisions:
!$Log: bfgrid.f90,v $
!Revision 1.1  2007/06/01 05:15:51  renq
!*** empty log message ***
!
!
      implicit none
      public

      real*8,dimension(:,:),allocatable,save :: brgridfc,bzgridfc

      end module bfgrid
