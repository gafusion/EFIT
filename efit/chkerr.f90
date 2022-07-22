!**********************************************************************
!! 
!>   chkerr checks for mhd fitting errors.
!!
!!   @param : mtime is time index
!!
!*********************************************************************
      subroutine chkerr(mtime)
      use set_kinds
      include 'eparm.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: mtime
      integer*4 m,k
      real*8 ercmin,chisqerr
      integer*4 kflag(nflag)
      data ercmin/0.01_dp/,chisqerr/80.0_dp/
!
      m=mtime
      erflag(m,:)=0
      if(tsaisq(m).ge.chisqerr) erflag(m,1)=1
      if(ali(m).ge.2.5_dp.or.ali(m).le.0.05_dp) erflag(m,2)=2
      if(betap(m).ge.6.0.or.betap(m).le.0.) erflag(m,3)=3
      if(abs((cpasma(m)-pasmat(m))/cpasma(m)).ge.0.08_dp) erflag(m,4)=4
      if((aout(m).ge.75.0).or.(aout(m).le.30.)) erflag(m,5)=5
      if(eout(m).le.0.8_dp.or.eout(m).ge.4.0) erflag(m,6)=6
      if(rout(m).gt.240..or.rout(m).lt.90.0) erflag(m,7)=7
      if(rcurrt(m).gt.240..or.rcurrt(m).lt.90.0) erflag(m,8)=8
      if(zout(m).gt.100..or.zout(m).lt.-100.) erflag(m,9)=9
      if(zcurrt(m).gt.100..or.zcurrt(m).lt.-100.) erflag(m,10)=10
      if(qsta(m).gt.200..or.qsta(m).lt.1.) erflag(m,13)=13
      if(betat(m).lt.0..or.betat(m).gt.25.) erflag(m,14)=14
      if(oleft(m).lt.-0.2_dp .or. oright(m).lt.-0.2_dp .or. otop(m).lt.-0.2_dp) &
        erflag(m,15)=15
      if (olefs(m).lt.-90.0) then
        if(qout(m).lt.1.) erflag(m,18)=18
      else
        if(qout(m).gt.200..or.qout(m).lt.1.) erflag(m,18)=18
      endif
      if(terror(m).ge.ercmin) erflag(m,19)=19
      if(dbpli(m).ge.0.05_dp) erflag(m,20)=20
      if(delbp(m).ge.0.08_dp) erflag(m,21)=21
      if ((eout(m).le.elomin).and.(fwtdlc.le.0.0)) then
        betap(m)=0.0
        betat(m)=0.0
        ali(m)=0.0
        wplasm(m)=0.0
        terror(m)=0.0
        erflag(m,3)=0
        erflag(m,2)=0
        erflag(m,14)=0
        erflag(m,19)=0
      endif
      kflag=0
      lflag=0
!
      if (sum(erflag(m,:)).eq.0) return
!----------------------------------------------------------------------
!--   write out errors to the terminal and error file
!----------------------------------------------------------------------
      open(unit=40,file='errfil.out',status='unknown',position='append')
      select case (ierchk)
      case (3)
        write(nttyo,980) ishot,time(mtime)
        write(40,980) ishot,time(mtime)
      case (2)
        write(nttyo,990) ishot,time(mtime)
        write(40,990) ishot,time(mtime)
      case default
        write(nttyo,1000) ishot,time(mtime)
        write(40,1000) ishot,time(mtime)
      end select
!
      do k=1,nflag
        if (erflag(m,k).gt.0) kflag(k)=erflag(m,k)
        if (erflag(m,k).gt.0) lflag=kflag(k)
        if (kflag(k).eq.1) write(nttyo,1010) chisqerr
        if (kflag(k).eq.1) write(40,1010) chisqerr
        if (kflag(k).eq.2) write(nttyo,1020)
        if (kflag(k).eq.2) write(40,1020)
        if (kflag(k).eq.3) write(nttyo,1025)
        if (kflag(k).eq.3) write(40,1025)
        if (kflag(k).eq.4) write(nttyo,1030)
        if (kflag(k).eq.4) write(40,1030)
        if (kflag(k).eq.5) write(nttyo,1040)
        if (kflag(k).eq.5) write(40,1040)
        if (kflag(k).eq.6) write(nttyo,1050)
        if (kflag(k).eq.6) write(40,1050)
        if (kflag(k).eq.7) write(nttyo,1060)
        if (kflag(k).eq.7) write(40,1060)
        if (kflag(k).eq.8) write(nttyo,1070)
        if (kflag(k).eq.8) write(40,1070)
        if (kflag(k).eq.9) write(nttyo,1080)
        if (kflag(k).eq.9) write(40,1080)
        if (kflag(k).eq.10) write(nttyo,1090)
        if (kflag(k).eq.10) write(40,1090)
        if (kflag(k).eq.13) write(nttyo,1100)
        if (kflag(k).eq.13) write(40,1100)
        if (kflag(k).eq.14) write(nttyo,1110)
        if (kflag(k).eq.14) write(40,1110)
        if (kflag(k).eq.15) write(nttyo,1120)
        if (kflag(k).eq.15) write(40,1120)
        if (kflag(k).eq.16) write(nttyo,1130)
        if (kflag(k).eq.16) write(40,1130)
        if (kflag(k).eq.18) write(nttyo,1150)
        if (kflag(k).eq.18) write(40,1150)
        if (kflag(k).eq.19) write(nttyo,1170) errmin
        if (kflag(k).eq.19) write(40,1170) errmin
        if (kflag(k).eq.20) write(nttyo,1180) dbpli(m)
        if (kflag(k).eq.20) write(40,1180) dbpli(m)
        if (kflag(k).eq.21) write(nttyo,1190) delbp(m)
        if (kflag(k).eq.21) write(40,1190) delbp(m)
      enddo
      close(unit=40)
!
      return
  980 format(2x,/,'  fit errors, shot ',i7,2x,f5.0, &
             ' msec. no eqdsks will be written')
  990 format(2x,/,'  fit errors, shot ',i7,2x,f5.0, &
             ' msec. only a-eqdsk will be written')
 1000 format(2x,/,'  fit errors, shot ',i7,2x,f5.0, &
             ' msec. all eqdsks are still written')
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
      end subroutine chkerr
