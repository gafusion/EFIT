!**********************************************************************
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          write dprobe.dat file for EFIT                          **
!**                                                                  **
!**********************************************************************
      SUBROUTINE dprobe(mpnam2,lpname,patmp2)
     
      USE exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                      nfsum,nvsum,nvesel,nacoil,nw,nh,nwnh
      USE siloop
      USE cecoil
      USE fcoil
      USE pmodel
      USE consta
      USE input
      USE cacoil
      USE nio
      USE mprobe
      USE cvesel
      USE fshift
      USE var_filech
      
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
      REAL*8 :: patmp2(magpr2)
      CHARACTER*10:: mpnam2(magpr2),lpname(nsilop)
      
      NAMELIST/in3/mpnam2,xmp2,ymp2,amp2,smp2,rsi,zsi,wsi, &
                   as,as2,lpname,rsisvs,turnfc,patmp2, &
                   zacoil,wacoil,hacoil, &
                   rf,zf,fcid,wf,hf,wvs,hvs,avs,avs2,af,af2, &
                   re,ze,ecid,rvs,zvs,we,he,fcturn

      OPEN(unit=nout,status='unknown',file='dprobe.dat', &
           position='append',delim='quote')

      WRITE (nout,in3)
      CLOSE (nout)
      RETURN
      END SUBROUTINE dprobe


!**********************************************************************
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          write general machinein for dprobe.dat                  **
!**          Assumes dprobes to one toroidal angle                   **
!**          exparm variables passed to handle conflicting names     **
!**          with efit                                               **
!**                                                                  **
!**********************************************************************
      SUBROUTINE dprobe_machinein(nfcoil_efund,nsilop_efund,magpr2_efund,&
                                  nrogow_efund,necoil_efund,nesum_efund,&
                                  nfsum_efund,nvsum_efund,nvesel_efund,nacoil_efund)
     
      USE nio
      USE exparm, only : magpri67,magpri322,magprirdp,magudom,maglds,& 
                 nnece,nnecein,neceo,mse315,mse45,mse15, &
                 mse1h,mse315_2,mse210,mpress,libim,nmsels, &
                 nnnte,ngam_vars,ngam_u,ngam_w,nlimit,nlimbd,nangle, &
                 ntangle,nfbcoil,mccoil,micoil,ndata,nwwcur, &
                 nffcur,nppcur,nercur,ntime,ndim,kxiter,mqwant, & 
                 mbdry,mbdry1,nxtram,nxtlim,nco2v,nco2r,modef, &
                 modep,modew,kubics,icycred_loopmax,nfourier,nsilds

      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
      
      NAMELIST/machinein/nsilds,nsilol,nfcoil,nrogow,nacoil,mfcoil,necoil,nvesel, &
          mpress,nesum,magpri67,magpri322,magprirdp,magudom,maglds,mse315,mse45, &
          mse15,mse1h,mse315_2,mse210,libim,nmsels,nnece,nnecein,neceo,nnnte, &
          ngam_vars,ngam_u,ngam_w,nlimit,nlimbd,nangle,ntangle,nfbcoil,mccoil, &
          micoil,ndata,nwwcur,nffcur,nppcur,nercur,ntime,ndim,kxiter,mqwant, &
          mbdry,mbdry1,nxtram,nxtlim,nco2v,nco2r,modef,modep,modew,kubics, &
          icycred_loopmax,nfourier
           

      nsilol = nsilop_efund - nsilds
      nfcoil = nfsum_efund
      nrogow = nrogow_efund
      nacoil = nacoil_efund
      mfcoil = nfcoil_efund
      necoil = necoil_efund
      nvesel = nvesel_efund
      nesum = nesum_efund
          
      OPEN(unit=nout,status='unknown',file='dprobe.dat')

      WRITE (nout,machinein)
      CLOSE (nout)
      RETURN
      END SUBROUTINE dprobe_machinein
