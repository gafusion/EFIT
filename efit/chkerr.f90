!**********************************************************************
!! 
!>   chkerr checks for mhd fitting errors.
!!
!!   @param mtime : time index
!!
!*********************************************************************
      subroutine chkerr(mtime)
      use errlims
      include 'eparm.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: mtime
      integer*4 m
      real*8 error_max,chisq_max
!
      if (iconvr.eq.2) then
        error_max=errmin
        chisq_max=saicon
      else
        error_max=error
        chisq_max=saimin
      endif
      m=mtime
      erflag(m,:)=0
      if(chisq(m).ge.chisq_max) erflag(m,1)=1
      if(li(m).ge.li_max .or. li(m).le.li_min) erflag(m,2)=2
      if(betap(m).ge.betap_max .or. betap(m).le.0.) erflag(m,3)=3
      if(abs((ipmhd(m)-ipmeas(m))/ipmhd(m)).ge.plasma_diff) erflag(m,4)=4
      if(aminor(m).ge.aminor_max .or. aminor(m).le.aminor_min) erflag(m,5)=5
      if(elong(m).ge.elong_max .or. elong(m).le.elong_min) erflag(m,6)=6
      if(rout(m).gt.rout_max .or. rout(m).lt.rout_min) erflag(m,7)=7
      if(zout(m).gt.zout_max .or. zout(m).lt.zout_min) erflag(m,9)=9
      if(rcurrt(m).gt.rcurrt_max .or. rcurrt(m).lt.rcurrt_min) erflag(m,8)=8
      if(zcurrt(m).gt.zcurrt_max .or. zcurrt(m).lt.zcurrt_min) erflag(m,10)=10
      if(qstar(m).gt.qstar_max .or. qstar(m).lt.qstar_min) erflag(m,13)=13
      if(betat(m).gt.betat_max .or. betat(m).lt.0.) erflag(m,14)=14
      if(gapin(m).lt.gapin_min .or. gapout(m).lt.gapout_min  &
                               .or. gaptop(m).lt.gaptop_min) &
        erflag(m,15)=15
      if (sepin(m).lt.sepin_check) then
        if(qout(m).lt.qout_min) erflag(m,18)=18
      else
        if(qout(m).gt.qout_max .or. qout(m).lt.qout_min) erflag(m,18)=18
      endif
      if(terror(m).ge.error_max) erflag(m,19)=19
      if(dbpli(m).ge.dbpli_diff) erflag(m,20)=20
      if(delbp(m).ge.delbp_diff) erflag(m,21)=21
      if ((elong(m).le.elomin).and.(fwtdlc.le.0.0)) then
        betap(m)=0.0
        betat(m)=0.0
        li(m)=0.0
        wmhd(m)=0.0
        !terror(m)=0.0
        erflag(m,3)=0
        erflag(m,2)=0
        erflag(m,14)=0
        erflag(m,19)=0
      endif
      lflag=0
!
      if(sum(erflag(m,:)).eq.0) return
!----------------------------------------------------------------------
!--   write out errors to the terminal and error file
!----------------------------------------------------------------------
      open(unit=40,file='errfil.out',status='unknown',position='append')
      select case (ierchk)
      case (2,3)
        write(nttyo,980) ishot,time(mtime)
        write(40,980) ishot,time(mtime)
      case (-1,-2,-3)
        write(nttyo,990) ishot,time(mtime)
        write(40,990) ishot,time(mtime)
      case default
        write(nttyo,1000) ishot,time(mtime)
        write(40,1000) ishot,time(mtime)
      end select
!
      if (erflag(m,1).eq.1) then
        write(nttyo,1010) chisq_max
        write(40,1010) chisq_max
        if (abs(ierchk).ne.2) lflag=erflag(m,1)
      endif
      if (erflag(m,2).eq.2) then
        write(nttyo,1020) li_max,li_min
        write(40,1020) li_max,li_min
        lflag=erflag(m,2)
      endif
      if (erflag(m,3).eq.3) then
        write(nttyo,1025) betap_max
        write(40,1025) betap_max
        lflag=erflag(m,3)
      endif
      if (erflag(m,4).eq.4) then
        write(nttyo,1030) plasma_diff
        write(40,1030) plasma_diff
        lflag=erflag(m,4)
      endif
      if (erflag(m,5).eq.5) then
        write(nttyo,1040) aminor_max,aminor_min
        write(40,1040) aminor_max,aminor_min
        lflag=erflag(m,5)
      endif
      if (erflag(m,6).eq.6) then
        write(nttyo,1050) elong_max,elong_min
        write(40,1050) elong_max,elong_min
        lflag=erflag(m,6)
      endif
      if (erflag(m,7).eq.7) then
        write(nttyo,1060) rout_max,rout_min
        write(40,1060) rout_max,rout_min
        lflag=erflag(m,7)
      endif
      if (erflag(m,8).eq.8) then
        write(nttyo,1070) rcurrt_max,rcurrt_min
        write(40,1070) rcurrt_max,rcurrt_min
        lflag=erflag(m,8)
      endif
      if (erflag(m,9).eq.9) then
        write(nttyo,1080) zout_max,zout_min
        write(40,1080) zout_max,zout_min
        lflag=erflag(m,9)
      endif
      if (erflag(m,10).eq.10) then
        write(nttyo,1090) zcurrt_max,zcurrt_min
        write(40,1090) zcurrt_max,zcurrt_min
        lflag=erflag(m,10)
      endif
      if (erflag(m,13).eq.13) then
        write(nttyo,1100) qstar_max,qstar_min
        write(40,1100) qstar_max,qstar_min
        lflag=erflag(m,13)
      endif
      if (erflag(m,14).eq.14) then
        write(nttyo,1110) betat_max
        write(40,1110) betat_max
        lflag=erflag(m,14)
      endif
      if (erflag(m,15).eq.15) then
        write(nttyo,1120) gapin_min,gapout_min,gaptop_min
        write(40,1120) gapin_min,gapout_min,gaptop_min
        lflag=erflag(m,15)
      endif
      if (erflag(m,18).eq.18) then 
        if (sepin(m).lt.sepin_check) then
          write(nttyo,1150) qout_min
          write(40,1150) qout_min
        else
          write(nttyo,1160) qout_max,qout_min
          write(40,1160) qout_max,qout_min
        endif
        lflag=erflag(m,18)
      endif
      if (erflag(m,19).eq.19) then 
        write(nttyo,1170) error_max
        write(40,1170) error_max
        if (abs(ierchk).ne.2) lflag=erflag(m,19)
      endif
      if (erflag(m,20).eq.20) then
        write(nttyo,1180) dbpli(m)
        write(40,1180) dbpli(m)
        lflag=erflag(m,20)
      endif
      if (erflag(m,21).eq.21) then
        write(nttyo,1190) delbp(m)
        write(40,1190) delbp(m)
        lflag=erflag(m,21)
      endif
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
 1020 format(5x,'Error #2, Li > ',f6.0,' or < ',f6.0)
 1025 format(5x,'Error #3, Betap > ',f6.0,' or < 0')
 1030 format(5x,'Error #4, (MHD Ip-Exp Ip)/MHD Ip > ',f6.0)
 1040 format(5x,'Error #5, a > ',f6.0,' or < ',f6.0)
 1050 format(5x,'Error #6, b/a > ',f6.0,' or < ',f6.0)
 1060 format(5x,'Error #7, Rout > ',f6.0,' or < ',f6.0)
 1070 format(5x,'Error #8, Rcurrt > ',f6.0,' or < ',f6.0)
 1080 format(5x,'Error #9, Zout > ',f6.0,' or < ',f6.0)
 1090 format(5x,'Error #10, Zcurrt > ',f6.0,' or < ',f6.0)
 1100 format(5x,'Error #13, Q* > ',f6.0,' or < ',f6.0)
 1110 format(5x,'Error #14, Betat > ',f6.0,' or < 0')
 1120 format(5x,'Error #15, Oleft<',f6.0,' or Oright<',f6.0' or Otop<',f6.0)
 1150 format(5x,'Error #18, Qout < ',f6.0)
 1160 format(5x,'Error #18, Qout > ',f6.0,' or < ',f6.0)
 1170 format(5x,'Error #19, error > ',e10.3)
 1180 format(5x,'Error #20, Bp+li/2 not consistent , error = ',e10.3)
 1190 format(5x,'Error #21, Bp not consistent , error = ',e10.3)
      end subroutine chkerr
