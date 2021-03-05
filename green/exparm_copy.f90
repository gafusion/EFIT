      module exparm
! module for experimental parameters
! Revisions:
! $Log: exparm.f90,v $
! Revision 1.1.2.1  2008/11/18 22:24:35  radhakri
! *** empty log message ***
!
! Revision 1.1  2007/06/01 05:16:13  renq
! *** empty log message ***
!nesum changed from 18 to 6 by srini to match with public domain ver 09/25/2008
!
      implicit none
      public
      
      integer*4,parameter :: nfcoil = 18,nsilop = 44,magpr2 = 76,&
                           nrogow = 1,necoil = 122,nesum = 6,&
                           nfsum = 18,nvsum = 24,nvesel = 24,&
                           nacoil = 1,mgaus1 = 8,mgaus2 = 10
      integer*4 :: nw,nh,nwnh
      end module exparm
