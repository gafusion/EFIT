      subroutine chkerr(mtime)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          chkerr checks for mhd fitting errors.                   **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          04/10/83..........first created                         **
!**          03/05/84..........increased error band on chisq to 25.  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use set_kinds
      include 'eparmdud129.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!      include 'ecomdu1.f90'
      data ercmin/0.01_dp/
!
      m=mtime
      do 20 k=1,30
        erflag(m,k)=0
   20 continue
      if (imagsigma.eq.0) then
        chisqerr=80.0
      elseif (imagsigma.gt.0) then
        chisqerr=300.0
      endif
      if (tsaisq(m).ge.chisqerr) erflag(m,1)=1
      if (ali(m).ge.2.5_dp.or.ali(m).le.0.05_dp) erflag(m,2)=2
      if (betap(m).ge.6.0.or.betap(m).le.0.) erflag(m,3)=3
      if (abs((cpasma(m)-pasmat(m))/cpasma(m)).ge.0.08_dp) erflag(m,4)=4
      if ((aout(m).ge.75.0).or.(aout(m).le.30.)) erflag(m,5)=5
      if (eout(m).le.0.8_dp.or.eout(m).ge.4.0) erflag(m,6)=6
      if (rout(m).gt.240..or.rout(m).lt.90.0) erflag(m,7)=7
      if (rcurrt(m).gt.240..or.rcurrt(m).lt.90.0) erflag(m,8)=8
      if (zout(m).gt.100..or.zout(m).lt.-100.) erflag(m,9)=9
      if (zcurrt(m).gt.100..or.zcurrt(m).lt.-100.) erflag(m,10)=10
      if (qsta(m).gt.200..or.qsta(m).lt.1.) erflag(m,13)=13
      if (betat(m).lt.0..or.betat(m).gt.25.) erflag(m,14)=14
      if (oleft(m).lt.-0.2_dp .or. oright(m).lt.-0.2_dp .or. otop(m).lt.-0.2_dp) &
                erflag(m,15)=15
      if (olefs(m).lt.-90.0) go to 32
      if (qout(m).gt.200..or.qout(m).lt.1.) erflag(m,18)=18
      go to 34
   32 continue
      if (qout(m).lt.1.) erflag(m,18)=18
   34 continue
      if (terror(m).ge.ercmin) erflag(m,19)=19
      if (dbpli(m).ge.0.05_dp) erflag(m,20)=20
      if (delbp(m).ge.0.08_dp) erflag(m,21)=21
      if ((eout(m).gt.elomin).or.(fwtdlc.gt.0.0)) go to 50
      betap(m)=0.0
      betat(m)=0.0
      ali(m)=0.0
      wplasm(m)=0.0
      terror(m)=0.0
      erflag(m,3)=0
      erflag(m,2)=0
      erflag(m,14)=0
      erflag(m,19)=0
   50 continue
      do 55 j=1,30
        kflag(j)=0
   55 continue
      lflag=0
!
      do 65 k = 1,30
      if(erflag(m,k).ne.0)go to 68
   65 continue
      return
   68 continue
!----------------------------------------------------------------------
!-- now write out errors to the erminal and error file ...
!----------------------------------------------------------------------
      open(unit=40,file='errfil.out',status='unknown',access='append' &
                                 )
      ktimeo=0
      ictr=mtime+ktimeo
      write (nttyo,1000)ishot,time(ictr)
      write(40,1000)ishot,time(ictr)
!
      do 60 k=1,30
      if (erflag(m,k).gt.0) kflag(k)=erflag(m,k)
      if (erflag(m,k).gt.0) lflag=kflag(k)
      if(kflag(k).eq.1)write (nttyo,1010) chisqerr
      if(kflag(k).eq.1)write(40,1010) chisqerr
      if(kflag(k).eq.2)write (nttyo,1020)
      if(kflag(k).eq.2)write(40,1020)
      if(kflag(k).eq.3)write (nttyo,1025)
      if(kflag(k).eq.3)write(40,1025)
      if(kflag(k).eq.4)write (nttyo,1030)
      if(kflag(k).eq.4)write(40,1030)
      if(kflag(k).eq.5)write (nttyo,1040)
      if(kflag(k).eq.5)write(40,1040)
      if(kflag(k).eq.6)write (nttyo,1050)
      if(kflag(k).eq.6)write(40,1050)
      if(kflag(k).eq.7)write (nttyo,1060)
      if(kflag(k).eq.7)write(40,1060)
      if(kflag(k).eq.8)write (nttyo,1070)
      if(kflag(k).eq.8)write(40,1070)
      if(kflag(k).eq.9)write (nttyo,1080)
      if(kflag(k).eq.9)write(40,1080)
      if(kflag(k).eq.10)write (nttyo,1090)
      if(kflag(k).eq.10)write(40,1090)
      if(kflag(k).eq.13)write (nttyo,1100)
      if(kflag(k).eq.13)write(40,1100)
      if(kflag(k).eq.14)write (nttyo,1110)
      if(kflag(k).eq.14)write(40,1110)
      if(kflag(k).eq.15)write (nttyo,1120)
      if(kflag(k).eq.15)write(40,1120)
      if(kflag(k).eq.16)write (nttyo,1130)
      if(kflag(k).eq.16)write(40,1130)
      if(kflag(k).eq.18)write (nttyo,1150)
      if(kflag(k).eq.18)write(40,1150)
      if(kflag(k).eq.19)write (nttyo,1170) errmin
      if(kflag(k).eq.19)write(40,1170) errmin
      if(kflag(k).eq.20)write (nttyo,1180) dbpli(m)
      if(kflag(k).eq.20)write(40,1180) dbpli(m)
      if(kflag(k).eq.21)write (nttyo,1190) delbp(m)
      if(kflag(k).eq.21)write(40,1190) delbp(m)
   60 continue
      close(unit=40)
!
      return
 1000 format(2x,/,'  fit errors, shot ',i7,2x,f5.0, &
             ' msec. no eqdsks will be written')
 1010 format(5x,'Error #1, Chisq > ',f6.0)
 1020 format(5x,'Error #2, Li > 2.5 or < 0.05')
 1025 format(5x,'Error #3, Betap > 6.0 or < 0.')
 1030 format(5x,'Error #4, (MHD Ip-Exp Ip)/MHD Ip > 8%')
 1040 format(5x,'Error #5, a large, small ')
 1050 format(5x,'Error #6, b/a < 0.8 or > 2.5')
 1060 format(5x,'Error #7 Rout large, small    ')
 1070 format(5x,'Error #8, Rcurrt > large, small  ')
 1080 format(5x,'Error #9, Zout > large, small ')
 1090 format(5x,'Error #10, Zcurrt > large, small ')
 1100 format(5x,'Error #13, Q* > 200. or < 1.')
 1110 format(5x,'Error #14, Betat > 25. or < 0.')
 1120 format(5x,'Error #15, Oleft<-.2 or Oright<-.2 or Otop<-.2')
 1130 format(5x,'Error #16, CO2 chord lengths large, small  ')
 1150 format(5x,'Error #18, Qout > 200. or < 1.')
 1170 format(5x,'Error #19, error > ',e10.3)
 1180 format(5x,'Error #20, Bp+li/2 not consistent , error = ',e10.3)
 1190 format(5x,'Error #21, Bp not consistent , error = ',e10.3)
      end
      subroutine lenco2(xplt,yplt,nplt,jges)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          lenco2 calculates the co2 path lengths.                 **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          05/08/83..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use set_kinds
      include 'eparmdud129.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!      include 'ecomdu1.f90'
      dimension xplt(*),yplt(*)
      dimension zuper(nco2v),zlower(nco2v),rco2(nco2r),rco2in(nco2r)
!
      do 45 i=1,nco2r
        rco2(i)=100.
        rco2in(i)=0.
   45 continue
      do 50 i=1,nco2v
        zuper(i)=100.
        zlower(i)=0.
   50 continue
      zuperts(jges)=100.
      zlowerts=0.
      do 300 i=1,nplt-1
        do 200 k=1,nco2r
          yr1=chordr(k)-yplt(i)
          yr2=chordr(k)-yplt(i+1)
          if (yr1*yr2.gt.0.0) go to 200
          rpath=xplt(i)+yr1*(xplt(i+1)-xplt(i))/(yplt(i+1)-yplt(i))
          if (rpath.ge.rcentr) rco2(k)=rpath
          if (rpath.le.rcentr) rco2in(k)=rpath
  200   continue
          yr1=zlibim   -yplt(i)
          yr2=zlibim   -yplt(i+1)
          if (yr1*yr2.gt.0.0) go to 202
          rpath=xplt(i)+yr1*(xplt(i+1)-xplt(i))/(yplt(i+1)-yplt(i))
          if (rpath.ge.rcentr) rlibim(jges)=rpath*100.0
  202     continue
        do 100 k=1,nco2v
          yr1=chordv(k)-xplt(i)
          yr2=chordv(k)-xplt(i+1)
          if (yr1*yr2.gt.0.0) go to 100
          zpath=yplt(i)+yr1*(yplt(i+1)-yplt(i))/(xplt(i+1)-xplt(i))
          if (zpath.ge.zcentr) zuper(k)=zpath
          if (zpath.le.zcentr) zlower(k)=zpath
  100   continue
          yr1=rmajts-xplt(i)
          yr2=rmajts-xplt(i+1)
          if (yr1*yr2.gt.0.0) go to 280
          zpath=yplt(i)+yr1*(yplt(i+1)-yplt(i))/(xplt(i+1)-xplt(i))
          if (zpath.ge.zcentr) zuperts(jges)=zpath*100.
          if (zpath.le.zcentr) zlowerts=zpath*100.
  280   continue
  300 continue
      do 320 k=1,nco2v
        rco2v(k,jges)=100.0*(zuper(k)-zlower(k))
        dco2v(jges,k)=denvt(jges,k)/rco2v(k,jges)
  320 continue
      do 330 k=1,nco2r
        rco2r(k,jges)=100.0*(rco2(k)-rco2in(k))
        dco2r(jges,k)=denrt(jges,k)/rco2r(k,jges)
  330 continue
!
      return
      end
