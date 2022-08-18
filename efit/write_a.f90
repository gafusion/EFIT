!**********************************************************************
!>
!!    writes out the plasma shape and quality parameters
!!    
!!
!!    @param ktime : Number of time slices
!!
!!    @param jtime : Time index
!!
!**********************************************************************
      subroutine write_a(ktime,jtime)
      use set_kinds
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      real*8,dimension(:),allocatable :: coils,expmp2
      namelist/in1/ishot,itime,qvfit,plasma,expmp2,coils,btor,ierchk, &
           fwtsi,fwtcur,limitr,xlim,ylim,fwtmp2,brsp,kffcur,kppcur, &
           pcurbd,fcurbd,kzeroj,rzeroj,vbit
      namelist/basis/kppfnc,kppknt,ppknt,pptens, &
                   kfffnc,kffknt,ffknt,fftens, &
                   kwwfnc,kwwknt,wwknt,wwtens, &
                   ppbdry,pp2bdry,kppbdry,kpp2bdry, &
                   ffbdry,ff2bdry,kffbdry,kff2bdry, &
                   wwbdry,ww2bdry,kwwbdry,kww2bdry, &
                   keefnc,keeknt,eeknt,eetens, &
                   eebdry,ee2bdry,keebdry,kee2bdry
      namelist/inwant/psiwant,vzeroj
      namelist/invt/kwwcur,kvtor,rvtor,wcurbd
      character(14) :: sfile
      character eqdsk*72,header*42,qmflag*3,fit_type*3
      character wform*20,let
      data nlold/40/,nlnew/41/
      save xdum

      allocate(coils(nsilop),expmp2(magpri))
      xdum = 0
      cprof=icprof
      iyes=0
      idup=0
      tavem=2*iavem
      jflag=1
      rcencm=rcentr*100.
      ktime1=1
      mco2v=nco2v
      mco2r=nco2r
      ! zero unused variables before write (for consistency)
      if (keecur.le.0) then
        eeknt=0.0
        keefnc=0
        keeknt=0
        slantl=0.0
        slantu=0.0
        vsurfa=0.0
        wpdot=0.0
        wbdot=0.0
      endif
!----------------------------------------------------------------------
!--   set fit type                                                   --
!----------------------------------------------------------------------
      if ((kprfit.gt.0).and.(kstark.gt.0)) then
        fit_type = 'KIM'
      elseif (kprfit.gt.0) then
        fit_type = 'KIN'
      elseif (kstark.gt.0) then
        fit_type = 'MSE'
      else
        fit_type = 'MAG'
      endif
!----------------------------------------------------------------------
!--   set fit type equilibrium mode                                  --
!----------------------------------------------------------------------
      if (iconvr.eq.3) then
        fit_type = 'EQU'
        if(npsi_ext.gt.0) fit_type='EQG'
      endif
      ijtime=time(jtime)
      jj=jtime
!-----------------------------------------------------------------------
!--   check error, on return lflag > 0 for error and the type of error
!-----------------------------------------------------------------------
      lflag=0
      if(ierchk.gt.0) call chkerr(jj)
      if (lflag.gt.0) then
        if(ierchk.eq.3) return
        jflag=0
      endif
      if (keqdsk.ge.1) then
        wform='formatted'
      else
        wform='unformatted'
      endif
!
      let = 'a'
      call setfnmeq(itimeu,let,ishot,ijtime,eqdsk)
!----------------------------------------------------------------------
!--   If (ISTORE = 1) Then                                           --
!--   Central directory to collect EFIT results is store_dir         --
!----------------------------------------------------------------------
      if (istore .eq. 1) then
        eqdsk = store_dir(1:lstdir)//eqdsk
      endif
      header = ' '
      open(unit=neqdsk,file=eqdsk,status='old', &
           form=wform,iostat=ioerr)
      if(ioerr.eq.0) close(unit=neqdsk,status='delete')
      open(unit=neqdsk,file=eqdsk,status='new', &
           form=wform,delim='quote')
!
      xbetapr=kbetapr
      a_eqdsk_format: if (keqdsk.ge.1) then
      write (neqdsk,1055) efitdatealt,efitdate
      if (ishot.le.99999) then
        write (neqdsk,1050) ishot,ktime1
      else
        write (neqdsk,1053) ishot,ktime1
      endif
      write (neqdsk,1040) (time(jj))
!
      if(fwtqa.gt.0.0.and.qvfit.gt.0.0) qmflag='FIX'
      if(fwtqa.le.0.0) qmflag='CLC'
!
      write (neqdsk,1060) time(jj),jflag,lflag,limloc(jj), &
                          mco2v,mco2r,qmflag,nlold,nlnew
      write (neqdsk,1040) tsaisq(jj),rcencm,bcentr(jj),pasmat(jj)
      write (neqdsk,1040) cpasma(jj),rout(jj),zout(jj),aout(jj)
      write (neqdsk,1040) eout(jj),doutu(jj),doutl(jj),vout(jj)
      write (neqdsk,1040) rcurrt(jj),zcurrt(jj),qsta(jj),betat(jj)
      write (neqdsk,1040) betap(jj),ali(jj),oleft(jj),oright(jj)
      write (neqdsk,1040) otop(jj),obott(jj),qpsib(jj),vertn(jj)
      write (neqdsk,1040) (rco2v(k,jj),k=1,mco2v)
      write (neqdsk,1040) (dco2v(jj,k),k=1,mco2v)
      write (neqdsk,1040) (rco2r(k,jj),k=1,mco2r)
      write (neqdsk,1040) (dco2r(jj,k),k=1,mco2r)
      write (neqdsk,1040) shearb(jj),bpolav(jj),s1(jj),s2(jj)
      write (neqdsk,1040) s3(jj),qout(jj),olefs(jj),orighs(jj)
      write (neqdsk,1040) otops(jj),sibdry(jj),areao(jj),wplasm(jj)
      write (neqdsk,1040) terror(jj),elongm(jj),qqmagx(jj),cdflux(jj)
      write (neqdsk,1040) alpha(jj),rttt(jj),psiref(jj),xndnt(jj)
      write (neqdsk,1040) rseps(1,jj),zseps(1,jj),rseps(2,jj) &
                          ,zseps(2,jj)
      write (neqdsk,1040) sepexp(jj),obots(jj),btaxp(jj),btaxv(jj)
      write (neqdsk,1040) aaq1(jj),aaq2(jj),aaq3(jj),seplim(jj)
      write (neqdsk,1040) rmagx(jj),zmagx(jj),simagx(jj),taumhd(jj)
      fluxx=diamag(jj)*1.0e-03_dp
      write (neqdsk,1040) betapd(jj),betatd(jj),wplasmd(jj),fluxx
      write (neqdsk,1040) vloopt(jj),taudia(jj),qmerci(jj),tavem
!
      if (ishot.lt.91000) then
        write (neqdsk,1041) nsilop,magpri67+magpri322,nfcoil,nesum
        write (neqdsk,1040) (csilop(k,jj),k=1,nsilop), &
                            (cmpr2(k,jj),k=1,magpri67+magpri322)
      else
        write (neqdsk,1041) nsilop,magpri,nfcoil,nesum
        write (neqdsk,1040) (csilop(k,jj),k=1,nsilop), &
                            (cmpr2(k,jj),k=1,magpri)
      endif
      write (neqdsk,1040) (ccbrsp(k,jj),k=1,nfcoil)
      write (neqdsk,1040) (eccurt(jj,k),k=1,nesum)
!
      write (neqdsk,1040) pbinj(jj),rvsin(jj),zvsin(jj),rvsout(jj)
      write (neqdsk,1040) zvsout(jj),vsurfa(jj),wpdot(jj),wbdot(jj)
      write (neqdsk,1040) slantu(jj),slantl(jj),zuperts(jj),chipre
      write (neqdsk,1040) cjor95(jj),pp95(jj),ssep(jj),yyy2(jj)
      write (neqdsk,1040) xnnc(jj),cprof,oring(jj),cjor0(jj)
      write (neqdsk,1040) fexpan,qqmin,chigamt,ssi01
      write (neqdsk,1040) fexpvs,sepnose,ssi95(jj),rqqmin
      write (neqdsk,1040) cjor99(jj),cj1ave(jj),rmidin(jj),rmidout(jj)
      write (neqdsk,1040) psurfa(jj),peak(jj),dminux(jj),dminlx(jj)
      write (neqdsk,1040) dolubaf(jj),dolubafm(jj),diludom(jj), &
                          diludomm(jj)
      write (neqdsk,1040) ratsol(jj),rvsiu(jj),zvsiu(jj),rvsid(jj)
      write (neqdsk,1040) zvsid(jj),rvsou(jj),zvsou(jj),rvsod(jj)
      write (neqdsk,1040) zvsod(jj),condno,psin32(jj),psin21(jj)
      write (neqdsk,1040) rq32in(jj),rq21top(jj),chilibt,ali3(jj)
      write (neqdsk,1040) xbetapr,tflux(jj),tchimls,twagap(jj)
!-----------------------------------------------------------------------
!--   one xxx replaced with ring gap 8/15/91                        --
!--   cjor0=flux surface average current density normalized to      --
!--         I/A at magnetic axis                                    --
!--   cjor95=jor95/jor0, cjor99=jor99/(I/A)                         --
!--   xbetapr=kbetapr                        2014/06/25             --
!-----------------------------------------------------------------------
      write (neqdsk,1042) header,fit_type
!-------------------------------------------------------------------
!--   binary format                                             --
!-------------------------------------------------------------------
      else a_eqdsk_format
      write (neqdsk) efitdatealt,efitdate
      if (ishot.le.99999) then
        write (neqdsk) ishot,ktime1
      else
        write (neqdsk) ishot,ktime1
      endif
      write (neqdsk) (time(jj))
!
      if(fwtqa.gt.0.0.and.qvfit.gt.0.0) qmflag='FIX'
      if(fwtqa.le.0.0) qmflag='CLC'
!
      write (neqdsk) real(time(jj),r4),jflag,lflag,limloc(jj), &
                          mco2v,mco2r,qmflag,nlold,nlnew
      write (neqdsk) real(tsaisq(jj),r4),real(rcencm,r4), &
                     real(bcentr(jj),r4),real(pasmat(jj),r4)
      write (neqdsk) real(cpasma(jj),r4),real(rout(jj),r4), &
                     real(zout(jj),r4),real(aout(jj),r4)
      write (neqdsk) real(eout(jj),r4),real(doutu(jj),r4), &
                     real(doutl(jj),r4),real(vout(jj),r4)
      write (neqdsk) real(rcurrt(jj),r4),real(zcurrt(jj),r4), &
                     real(qsta(jj),r4),real(betat(jj),r4)
      write (neqdsk) real(betap(jj),r4),real(ali(jj),r4), &
                     real(oleft(jj),r4),real(oright(jj),r4)
      write (neqdsk) real(otop(jj),r4),real(obott(jj),r4), &
                     real(qpsib(jj),r4),real(vertn(jj),r4)
      write (neqdsk) (real(rco2v(k,jj),r4),k=1,mco2v)
      write (neqdsk) (real(dco2v(jj,k),r4),k=1,mco2v)
      write (neqdsk) (real(rco2r(k,jj),r4),k=1,mco2r)
      write (neqdsk) (real(dco2r(jj,k),r4),k=1,mco2r)
      write (neqdsk) real(shearb(jj),r4),real(bpolav(jj),r4), &
                     real(s1(jj),r4),real(s2(jj),r4)
      write (neqdsk) real(s3(jj),r4),real(qout(jj),r4), &
                     real(olefs(jj),r4),real(orighs(jj),r4)
      write (neqdsk) real(otops(jj),r4),real(sibdry(jj),r4), &
                     real(areao(jj),r4),real(wplasm(jj),r4)
      write (neqdsk) real(terror(jj),r4),real(elongm(jj),r4), &
                     real(qqmagx(jj),r4),real(cdflux(jj),r4)
      write (neqdsk) real(alpha(jj),r4),real(rttt(jj),r4), &
                     real(psiref(jj),r4),real(xndnt(jj),r4)
      write (neqdsk) real(rseps(1,jj),r4),real(zseps(1,jj),r4), &
                     real(rseps(2,jj),r4),real(zseps(2,jj),r4)
      write (neqdsk) real(sepexp(jj),r4),real(obots(jj),r4), &
                     real(btaxp(jj),r4),real(btaxv(jj),r4)
      write (neqdsk) real(aaq1(jj),r4),real(aaq2(jj),r4), &
                     real(aaq3(jj),r4),real(seplim(jj),r4)
      write (neqdsk) real(rmagx(jj),r4),real(zmagx(jj),r4), &
                     real(simagx(jj),r4),real(taumhd(jj),r4)
      fluxx=diamag(jj)*1.0e-03_dp
      write (neqdsk) real(betapd(jj),r4),real(betatd(jj),r4), &
                     real(wplasmd(jj),r4),real(fluxx,r4)
      write (neqdsk) real(vloopt(jj),r4),real(taudia(jj),r4), &
                     real(qmerci(jj),r4),real(tavem,r4)
!
      if (ishot.lt.91000) then
        write (neqdsk) nsilop,magpri67+magpri322,nfcoil,nesum
        write (neqdsk) (real(csilop(k,jj),r4),k=1,nsilop), &
                       (real(cmpr2(k,jj),r4),k=1,magpri67+magpri322)
      else
        write (neqdsk) nsilop,magpri,nfcoil,nesum
        write (neqdsk) (real(csilop(k,jj),r4),k=1,nsilop), &
                       (real(cmpr2(k,jj),r4),k=1,magpri)
      endif
      write (neqdsk) (real(ccbrsp(k,jj),r4),k=1,nfcoil)
      write (neqdsk) (real(eccurt(jj,k),r4),k=1,nesum)
!
      write (neqdsk) real(pbinj(jj),r4),real(rvsin(jj),r4), &
                     real(zvsin(jj),r4),real(rvsout(jj),r4)
      write (neqdsk) real(zvsout(jj),r4),real(vsurfa(jj),r4), &
                     real(wpdot(jj),r4),real(wbdot(jj),r4)
      write (neqdsk) real(slantu(jj),r4),real(slantl(jj),r4), &
                     real(zuperts(jj),r4),real(chipre,r4)
      write (neqdsk) real(cjor95(jj),r4),real(pp95(jj),r4), &
                     real(ssep(jj),r4),real(yyy2(jj),r4)
      write (neqdsk) real(xnnc(jj),r4),real(cprof,r4), &
                     real(oring(jj),r4),real(cjor0(jj),r4)
      write (neqdsk) real(fexpan,r4),real(qqmin,r4), &
                     real(chigamt,r4),real(ssi01,r4)
      write (neqdsk) real(fexpvs,r4),real(sepnose,r4), &
                     real(ssi95(jj),r4),real(rqqmin,r4)
      write (neqdsk) real(cjor99(jj),r4),real(cj1ave(jj),r4), &
                     real(rmidin(jj),r4),real(rmidout(jj),r4)
      write (neqdsk) real(psurfa(jj),r4),real(peak(jj),r4), &
                     real(dminux(jj),r4),real(dminlx(jj),r4)
      write (neqdsk) real(dolubaf(jj),r4),real(dolubafm(jj),r4), &
                     real(diludom(jj),r4),real(diludomm(jj),r4)
      write (neqdsk) real(ratsol(jj),r4),real(rvsiu(jj),r4), &
                     real(zvsiu(jj),r4),real(rvsid(jj),r4)
      write (neqdsk) real(zvsid(jj),r4),real(rvsou(jj),r4), &
                     real(zvsou(jj),r4),real(rvsod(jj),r4)
      write (neqdsk) real(zvsod(jj),r4),real(condno,r4), &
                     real(psin32(jj),r4),real(psin21(jj),r4)
      write (neqdsk) real(rq32in(jj),r4),real(rq21top(jj),r4), &
                     real(chilibt,r4),real(ali3(jj),r4)
      write (neqdsk) real(xbetapr,r4),real(tflux(jj),r4), &
                     real(tchimls,r4),real(twagap(jj),r4)
!
      write (neqdsk) header,fit_type
      endif a_eqdsk_format
!
      close(unit=neqdsk)
!
! --- add flag for writting out esave.dat. IOUT contains 16.
!
      if ((jtime.eq.ktime).and.(iand(iout,16).ne.0)) then
        if (nproc.gt.1) then
          WRITE(sfile,fmt='(i5.5)') rank
          sfile='esave'//TRIM(sfile)//'.dat'
        else
          sfile='esave.dat'
        endif
        open(unit=nsave,status='old',form='unformatted',file=sfile, &
             iostat=ioerr)
        if(ioerr.eq.0) close(unit=nsave,status='delete')
        open(unit=nsave,status='new',form='unformatted',file=sfile)
        write (nsave) nw,nh
        write (nsave) xpsi
        write (nsave) brsp
        write (nsave) www
        write (nsave) emaxis, rmaxis, fcentr
        close(unit=nsave)
      endif
!
      if (kinput.gt.0) then
        plasma=cpasma(jtime)
        btor=bcentr(jtime)
        coils(1:nsilop)=csilop(1:nsilop,jtime)
        expmp2(1:magpri)=cmpr2(1:magpri,jtime)
        ierold=ierchk
        ierchk=0
        fwtcur=1.
        fwtmp2(1:magpri)=1.
        fwtsi(1:nsilop)=1.
        iitime=time(jj)+1
        limitr=limitr-1
        qvfit=qenp
        let = 'm'
        call setfnmeq(itimeu,let,ishot,iitime,eqdsk)
        open(unit=neqdsk,file=eqdsk,status='old',iostat=ioerr)
        if(ioerr.eq.0) close(unit=neqdsk,status='delete')
        open(unit=neqdsk,file=eqdsk,status='new',delim='quote')
        itsave=itime
        itime=iitime
        write (neqdsk,in1)
        call ppstore
        call ffstore
        call wwstore
        call eestore
        write (neqdsk,basis)
        write (neqdsk,inwant)
        if(kvtor.gt.0) &
          write (neqdsk,invt)
        itime=itsave
        close(unit=neqdsk)
        limitr=limitr+1
        ierchk=ierold
      endif
      return
!
 1020 format ('p',i5,'.',i3)
 1040 format (1x,4e16.9)
 1041 format (1x,4i5)
 1042 format (1x,a42,1x,a3)
 1050 format (1x,i5,11x,i5)
 1053 format (1x,i6,11x,i5)
 1055 format (1x,2a10)
 1060 format ('*',f8.3,9x,i5,11x,i5,1x,a3,1x,i3,1x,i3,1x,a3,1x,2i5)
 1070 format(16a5)
      end subroutine write_a


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
      integer*4 kflag(nflag)
      ! TODO: all limits were choosen for DIII-D but likely need to
      !       change for other experiments and should be moved to
      !       dprobe.dat
      real*8, parameter :: ercmin=0.01_dp,chisqerr=80.0_dp
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
