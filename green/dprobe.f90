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
           re,ze,ecid,ecturn,vsid,rvs,zvs,we,he

      OPEN(unit=nout,status='new',file='dprobe.dat')

      WRITE (nout,in3)
      CLOSE (nout)
      RETURN
      END SUBROUTINE dprobe
