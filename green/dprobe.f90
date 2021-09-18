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
      
      NAMELIST/in3/ mpnam2,xmp2,ymp2,amp2,smp2,rsi,zsi,wsi, &
           as,as2,lpname,rsisvs,turnfc,patmp2, &
           zacoil,wacoil,hacoil, &
           rf,zf,fcid,wf,hf,wvs,hvs,avs,af,af2, &
           re,ze,ecid,rvs,zvs,we,he

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
      
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
      
      NAMELIST/machinein/nsilds,nsilol,nfcoil,nrogow,nacoil,mfcoil,necoil,nvesel, &
          mpress,nesum,magpri67,magpri322,magprirdp,magudom,maglds,mse315,mse45, &
          mse15,mse1h,mse315_2,mse210,libim,nmsels,nnece,nnecein,neceo,nnnte, &
          ngam_vars,ngam_u,ngam_w,nlimit,nlimbd,nangle,ntangle,nfbcoil,mccoil, &
          micoil,ndata,nwwcur,nffcur,nppcur,nercur,ntime,ndim,kxiter,mqwant, &
          mbdry,mbdry1,nxtram,nxtlim,nco2v,nco2r,modef,modep,modew,kubics, &
          icycred_loopmax,nfourier
           
      nsilds = 1 ! default all psi loops to one side except 1
      nsilol = nsilop_efund - nsilds
      nfcoil = nfsum_efund
      nrogow=nrogow_efund
      nacoil=nacoil_efund
      mfcoil=nfcoil_efund
      necoil = necoil_efund
      nvesel = nvesel_efund
      nesum= nesum_efund
      magpri67=1
      magprirdp=1
      magudom=1
      maglds=1
      magpri322 = magpr2_efund - magpri67 - magprirdp - magudom - maglds

      ! Remaining default to DIII-D values
      ! efund knows nothing about mse and ece array lengths

      nnece=40
      nnecein=80
      neceo=1
      mse315=40
      mse45=0
      mse15=0
      mse1h=0
      mse315_2=0
      mse210=0
      mpress=201
      libim=32
      nmsels=16
      nnnte=801
      ngam_vars=9
      ngam_u=5
      ngam_w=3
      nlimit=160
      nlimbd=6
      nangle=64
      ntangle=12
      nfbcoil=12
      mccoil=6
      micoil=12     
      ndata=61
      nwwcur=32
      nffcur=32
      nppcur=32
      nercur=32
      npcurn=nffcur+nppcur
      nwcurn=nwwcur+npcurn
      nwcur2=nwcurn*2
      ntime=1001
      ndim=3200
      kxiter=515
      mqwant=30
      mbdry=300
      mbdry1=110
      nxtram=10
      nxtlim=9
      nco2v=3
      nco2r=2
      modef=4
      modep=4
      modew=4
      kubics=4
      icycred_loopmax=1290
      nfourier=5
          
      OPEN(unit=nout,status='new',file='dprobe.dat')

      WRITE (nout,machinein)
      CLOSE (nout)
      RETURN
      END SUBROUTINE dprobe_machinein
      
      
!**********************************************************************
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          write general machinein for dprobe.dat                  **
!**          Assumes dprobes to one toroidal angle                   **
!**          exparm variables passed to handle conflicting names     **
!**          with efit                                               **
!**                                                                  **
!**********************************************************************
      SUBROUTINE dprobe_machinein_d3d
     
      USE nio
      
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
      
      NAMELIST/machinein/nsilds,nsilol,nfcoil,nrogow,nacoil,mfcoil,necoil,nvesel, &
          mpress,nesum,magpri67,magpri322,magprirdp,magudom,maglds,mse315,mse45, &
          mse15,mse1h,mse315_2,mse210,libim,nmsels,nnece,nnecein,neceo,nnnte, &
          ngam_vars,ngam_u,ngam_w,nlimit,nlimbd,nangle,ntangle,nfbcoil,mccoil, &
          micoil,ndata,nwwcur,nffcur,nppcur,nercur,ntime,ndim,kxiter,mqwant, &
          mbdry,mbdry1,nxtram,nxtlim,nco2v,nco2r,modef,modep,modew,kubics, &
          icycred_loopmax,nfourier
           
      nsilds=3
      nsilol=41
      nfcoil=18
      nrogow=1
      nacoil=1
      mfcoil=18
      necoil=122
      nvesel=24
      mpress=201
      nesum=6
      magpri67=29
      magpri322=31
      magprirdp=8
      magudom=5
      maglds=3
      mse315=11
      mse45=15
      mse15=10
      mse1h=4
      mse315_2=5
      mse210=24
      libim=32
      nmsels=16
      nnece=40
      nnecein=80
      neceo=1
      nnnte=801
      ngam_vars=9
      ngam_u=5
      ngam_w=3
      nlimit=160
      nlimbd=6
      nangle=64
      ntangle=12
      nfbcoil=12
      mccoil=6
      micoil=12

      ndata=61
      nwwcur=32
      nffcur=32
      nppcur=32
      nercur=32

      npcurn=nffcur+nppcur
      nwcurn=nwwcur+npcurn
      nwcur2=nwcurn*2
      ntime=1001
      ndim=3200
      kxiter=515
      mqwant=30
      mbdry=300
      mbdry1=110
      nxtram=10
      nxtlim=9
      nco2v=3
      nco2r=2
      modef=4
      modep=4
      modew=4
      kubics=4
      icycred_loopmax=1290
      nfourier=5
          
      OPEN(unit=nout,status='new',file='dprobe.dat')

      WRITE (nout,machinein)
      CLOSE (nout)
      RETURN
      END SUBROUTINE dprobe_machinein_d3d
      
