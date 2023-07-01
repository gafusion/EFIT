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
      use set_kinds, only: r4
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
      integer*4, parameter :: nlold=40,nlnew=41
      integer*4, parameter :: magpri67=29,magpri322=31
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
        if(abs(ierchk).gt.1) return
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
      write (neqdsk,1040) chisq(jj),rcencm,bcentr(jj),ipmeas(jj)
      write (neqdsk,1040) ipmhd(jj),rout(jj),zout(jj),aminor(jj)
      write (neqdsk,1040) elong(jj),utri(jj),ltri(jj),volume(jj)
      write (neqdsk,1040) rcurrt(jj),zcurrt(jj),qstar(jj),betat(jj)
      write (neqdsk,1040) betap(jj),li(jj),gapin(jj),gapout(jj)
      write (neqdsk,1040) gaptop(jj),gapbot(jj),q95(jj),vertn(jj)
      write (neqdsk,1040) (rco2v(k,jj),k=1,mco2v)
      write (neqdsk,1040) (dco2v(jj,k),k=1,mco2v)
      write (neqdsk,1040) (rco2r(k,jj),k=1,mco2r)
      write (neqdsk,1040) (dco2r(jj,k),k=1,mco2r)
      write (neqdsk,1040) shearb(jj),bpolav(jj),s1(jj),s2(jj)
      write (neqdsk,1040) s3(jj),qout(jj),sepin(jj),sepout(jj)
      write (neqdsk,1040) septop(jj),sibdry(jj),area(jj),wmhd(jj)
      write (neqdsk,1040) terror(jj),elongm(jj),qm(jj),cdflux(jj)
      write (neqdsk,1040) alpha(jj),rttt(jj),psiref(jj),indent(jj)
      write (neqdsk,1040) rseps(1,jj),zseps(1,jj),rseps(2,jj), &
                          zseps(2,jj)
      write (neqdsk,1040) sepexp(jj),sepbot(jj),btaxp(jj),btaxv(jj)
      write (neqdsk,1040) aaq1(jj),aaq2(jj),aaq3(jj),dsep(jj)
      write (neqdsk,1040) rm(jj),zm(jj),psim(jj),taumhd(jj)
      fluxx=diamag(jj)*1.0e-03_dp
      write (neqdsk,1040) betapd(jj),betatd(jj),wdia(jj),fluxx
      write (neqdsk,1040) vloopt(jj),taudia(jj),qmerci(jj),tavem
!
      if (ishot.lt.91000) then
        write (neqdsk,1041) nsilop,magpri67+magpri322,nfsum,nesum
        write (neqdsk,1040) (csilop(k,jj),k=1,nsilop), &
                            (cmpr2(k,jj),k=1,magpri67+magpri322)
      else
        write (neqdsk,1041) nsilop,magpri,nfsum,nesum
        write (neqdsk,1040) (csilop(k,jj),k=1,nsilop), &
                            (cmpr2(k,jj),k=1,magpri)
      endif
      write (neqdsk,1040) (ccbrsp(k,jj),k=1,nfsum)
      write (neqdsk,1040) (eccurt(jj,k),k=1,nesum)
!
      write (neqdsk,1040) pbinj(jj),rvsin(jj),zvsin(jj),rvsout(jj)
      write (neqdsk,1040) zvsout(jj),vsurfa(jj),wpdot(jj),wbdot(jj)
      write (neqdsk,1040) slantu(jj),slantl(jj),zuperts(jj),chipre
      write (neqdsk,1040) cjor95(jj),pp95(jj),drsep(jj),yyy2(jj)
      write (neqdsk,1040) xnnc(jj),cprof,oring(jj),cjor0(jj)
      write (neqdsk,1040) fexpan,qmin,chigamt,ssi01
      write (neqdsk,1040) fexpvs,sepnose,ssi95(jj),rhoqmin
      write (neqdsk,1040) cjor99(jj),cj1ave(jj),rmidin(jj),rmidout(jj)
      write (neqdsk,1040) psurfa(jj),peak(jj),dminux(jj),dminlx(jj)
      write (neqdsk,1040) dolubaf(jj),dolubafm(jj),diludom(jj), &
                          diludomm(jj)
      write (neqdsk,1040) ratsol(jj),rvsiu(jj),zvsiu(jj),rvsid(jj)
      write (neqdsk,1040) zvsid(jj),rvsou(jj),zvsou(jj),rvsod(jj)
      write (neqdsk,1040) zvsod(jj),condno,psin32(jj),psin21(jj)
      write (neqdsk,1040) rq32in(jj),rq21top(jj),chilibt,li3(jj)
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
      write (neqdsk) real(chisq(jj),r4),real(rcencm,r4), &
                     real(bcentr(jj),r4),real(ipmeas(jj),r4)
      write (neqdsk) real(ipmhd(jj),r4),real(rout(jj),r4), &
                     real(zout(jj),r4),real(aminor(jj),r4)
      write (neqdsk) real(elong(jj),r4),real(utri(jj),r4), &
                     real(ltri(jj),r4),real(volume(jj),r4)
      write (neqdsk) real(rcurrt(jj),r4),real(zcurrt(jj),r4), &
                     real(qstar(jj),r4),real(betat(jj),r4)
      write (neqdsk) real(betap(jj),r4),real(li(jj),r4), &
                     real(gapin(jj),r4),real(gapout(jj),r4)
      write (neqdsk) real(gaptop(jj),r4),real(gapbot(jj),r4), &
                     real(q95(jj),r4),real(vertn(jj),r4)
      write (neqdsk) (real(rco2v(k,jj),r4),k=1,mco2v)
      write (neqdsk) (real(dco2v(jj,k),r4),k=1,mco2v)
      write (neqdsk) (real(rco2r(k,jj),r4),k=1,mco2r)
      write (neqdsk) (real(dco2r(jj,k),r4),k=1,mco2r)
      write (neqdsk) real(shearb(jj),r4),real(bpolav(jj),r4), &
                     real(s1(jj),r4),real(s2(jj),r4)
      write (neqdsk) real(s3(jj),r4),real(qout(jj),r4), &
                     real(sepin(jj),r4),real(sepout(jj),r4)
      write (neqdsk) real(septop(jj),r4),real(sibdry(jj),r4), &
                     real(area(jj),r4),real(wmhd(jj),r4)
      write (neqdsk) real(terror(jj),r4),real(elongm(jj),r4), &
                     real(qm(jj),r4),real(cdflux(jj),r4)
      write (neqdsk) real(alpha(jj),r4),real(rttt(jj),r4), &
                     real(psiref(jj),r4),real(indent(jj),r4)
      write (neqdsk) real(rseps(1,jj),r4),real(zseps(1,jj),r4), &
                     real(rseps(2,jj),r4),real(zseps(2,jj),r4)
      write (neqdsk) real(sepexp(jj),r4),real(sepbot(jj),r4), &
                     real(btaxp(jj),r4),real(btaxv(jj),r4)
      write (neqdsk) real(aaq1(jj),r4),real(aaq2(jj),r4), &
                     real(aaq3(jj),r4),real(dsep(jj),r4)
      write (neqdsk) real(rm(jj),r4),real(zm(jj),r4), &
                     real(psim(jj),r4),real(taumhd(jj),r4)
      fluxx=diamag(jj)*1.0e-03_dp
      write (neqdsk) real(betapd(jj),r4),real(betatd(jj),r4), &
                     real(wdia(jj),r4),real(fluxx,r4)
      write (neqdsk) real(vloopt(jj),r4),real(taudia(jj),r4), &
                     real(qmerci(jj),r4),real(tavem,r4)
!
      if (ishot.lt.91000) then
        write (neqdsk) nsilop,magpri67+magpri322,nfsum,nesum
        write (neqdsk) (real(csilop(k,jj),r4),k=1,nsilop), &
                       (real(cmpr2(k,jj),r4),k=1,magpri67+magpri322)
      else
        write (neqdsk) nsilop,magpri,nfsum,nesum
        write (neqdsk) (real(csilop(k,jj),r4),k=1,nsilop), &
                       (real(cmpr2(k,jj),r4),k=1,magpri)
      endif
      write (neqdsk) (real(ccbrsp(k,jj),r4),k=1,nfsum)
      write (neqdsk) (real(eccurt(jj,k),r4),k=1,nesum)
!
      write (neqdsk) real(pbinj(jj),r4),real(rvsin(jj),r4), &
                     real(zvsin(jj),r4),real(rvsout(jj),r4)
      write (neqdsk) real(zvsout(jj),r4),real(vsurfa(jj),r4), &
                     real(wpdot(jj),r4),real(wbdot(jj),r4)
      write (neqdsk) real(slantu(jj),r4),real(slantl(jj),r4), &
                     real(zuperts(jj),r4),real(chipre,r4)
      write (neqdsk) real(cjor95(jj),r4),real(pp95(jj),r4), &
                     real(drsep(jj),r4),real(yyy2(jj),r4)
      write (neqdsk) real(xnnc(jj),r4),real(cprof,r4), &
                     real(oring(jj),r4),real(cjor0(jj),r4)
      write (neqdsk) real(fexpan,r4),real(qmin,r4), &
                     real(chigamt,r4),real(ssi01,r4)
      write (neqdsk) real(fexpvs,r4),real(sepnose,r4), &
                     real(ssi95(jj),r4),real(rhoqmin,r4)
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
                     real(chilibt,r4),real(li3(jj),r4)
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
        plasma=ipmhd(jtime)
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
