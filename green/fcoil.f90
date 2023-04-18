      module fcoil

      use machine, only: nfcoil,nfsum
      implicit none
      public

      integer*4 ifcoil
      real*8, dimension(:), allocatable :: rf,zf,wf,hf,af,af2,turnfc,fcturn
      integer*4, dimension(:), allocatable :: fcid

      end module fcoil
