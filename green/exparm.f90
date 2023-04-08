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


      integer*4::magpri67,magpri322,magprirdp,magudom,maglds,&
                 nnece,nnecein,neceo,mse315,mse45,mse15,&
                 mse1h,mse315_2,mse210,mpress,libim,nmsels,&
                 nnnte,ngam_vars,ngam_u,ngam_w,nlimit,nlimbd,nangle,&
                 ntangle,nfbcoil,mccoil,micoil,ndata,nwwcur,&
                 nffcur,nppcur,nercur,ntime,ndim,kxiter,mqwant,& 
                 mbdry,mbdry1,nxtram,nxtlim,nco2v,nco2r,modef,&
                 modep,modew,nfourier,nsilds,nsilol
      end module exparm
!errlims
      module errlims
        ! Experiment dependant checks on the solution error
        real*8 li_max,li_min,betap_max,plasma_diff, &
               aminor_max,aminor_min,elong_max,elong_min, &
               rout_max,rout_min,zout_max,zout_min, &
               rcurrt_max,rcurrt_min,zcurrt_max,zcurrt_min, &
               qstar_max,qstar_min,betat_max, &
               gapin_min,gapout_min,gaptop_min, &
               sepin_check,qout_max,qout_min, &
               dbpli_diff,delbp_diff
      end module errlims
