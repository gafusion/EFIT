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
      
      integer*4:: nfcoil,nsilop,magpr2,&
                           nrogow,necoil,nesum,&
                           nfsum,nvsum,nvesel,&
                           nacoil,mgaus1,mgaus2

      integer*4 :: nw,nh,nwnh
      character*10:: device


      integer*4::magprirdp,magudom,maglds,&
                 nnece,nnecein,neceo,mse315,mse45,mse15,&
                 mse1h,mse315_2,mse210,mpress,libim,nmsels,&
                 nnnte,ngam_vars,ngam_u,ngam_w,nlimit,nlimbd,nangle,&
                 ntangle,nfbcoil,mccoil,micoil,ndata,nwwcur,&
                 nffcur,nppcur,nercur,ntime,ndim,kxiter,mqwant,& 
                 mbdry,mbdry1,nxtram,nxtlim,nco2v,nco2r,modef,&
                 modep,modew,kubics,icycred_loopmax,nfourier,nsilds
      end module exparm
