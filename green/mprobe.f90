      module mprobe

      use exparm, only: magpr2
      implicit none
      public

      integer*4 :: nsmp2
      real*8,dimension (:), allocatable :: xmp2,ymp2,amp2,smp2

      end module mprobe
