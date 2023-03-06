      module cvesel
      use exparm, only: nvesel
      implicit none
      public

      real*8, dimension (:), allocatable :: rvs,zvs,wvs,hvs,avs,avs2, &
                                rsisvs
      integer*4, dimension(:), allocatable :: vsid

      end module cvesel
