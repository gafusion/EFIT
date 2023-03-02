      module fcoil
      use exparm, only: nfcoil
      implicit none
      public

      real*8, dimension(:), allocatable :: rf,zf,wf,hf,af,af2,turnfc,fcturn
      integer*4, dimension(:), allocatable :: fcid

      end module fcoil
