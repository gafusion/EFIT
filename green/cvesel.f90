      module cvesel
!module for cvesel
!Revisions:
!$Log: cvesel.f90,v $
!Revision 1.1  2007/06/01 05:16:07  renq
!*** empty log message ***
!
!
      use exparm,only:nvesel
      implicit none
      public

      real*8, dimension (:), allocatable :: rvs,zvs,wvs,hvs,avs,avs2,&
                                rsisvs,vsid

      end module cvesel
