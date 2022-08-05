!**********************************************************************
!>
!!    shipit writes out the plasma shape and quality parameters
!!    
!!
!!    @param ktime : Number of time slices
!!
!!    @param jtime : Time index
!!
!**********************************************************************
      subroutine shipit(ktime,jtime)
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
      if (lflag.eq.3) then
        return
      elseif (lflag.gt.0) then
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
      end subroutine shipit


!**********************************************************************
!>
!!    weqdsk writes out the GAQ type eqdsk.
!!    
!!
!!    @param jtime : time index
!!
!**********************************************************************
      subroutine weqdsk(jtime)
      use set_kinds
      use commonblocks,only: psirz
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8,dimension(:),allocatable :: coils,expmp2, prexp, &
                tgamma,sgamma,rrrgam, &
                zzzgam,aa1gam,aa2gam, aa3gam,aa4gam,aa5gam, &
                aa6gam,aa7gam,aa8gam,tgammauncor
      real*8,dimension(:),allocatable :: bmsels,sbmsels,fwtbmsels, &
                rrmsels,zzmsels,l1msels,l2msels, &
                l4msels,emsels,semsels,fwtemsels
      namelist/out1/ishot,itime,betap0,rzero,qenp,enp,emp,plasma, &
           expmp2,coils,btor,rcentr,brsp,icurrt,rbdry,zbdry, &
           nbdry,fwtsi,fwtcur,mxiter,nxiter,limitr,xlim,ylim,error, &
           iconvr,ibunmn,pressr,rpress,nqpsi,npress,sigpre
      namelist/mseout/rrrgam,zzzgam,aa1gam,aa2gam,aa3gam,aa4gam, &
                      tgamma,sgamma,aa5gam,aa6gam,aa7gam,aa8gam, &
                      msebkp,fwtgam,tgammauncor
      namelist/in_msels/bmsels,sbmsels,fwtbmsels,rrmsels,zzmsels, &
            l1msels,l2msels,l4msels,emsels,semsels
      namelist/vtout/npresw,presw,sigprw,rpresw
      namelist/in1/ishot,itime,betap0,rzero,qenp,enp,emp,plasma, &
           coils,btor,rcentr,icurrt,ierchk,kzeroj,rzeroj, &
           fwtsi,fwtcur,mxiter,nxiter,limitr,xlim,ylim,error, &
           iconvr,ibunmn,brsp,kffcur,kppcur,pcurbd,fcurbd,vbit
      namelist/basis/kppfnc,kppknt,ppknt,pptens, &
                   kfffnc,kffknt,ffknt,fftens, &
                   kwwfnc,kwwknt,wwknt,wwtens, &
                   ppbdry,pp2bdry,kppbdry,kpp2bdry, &
                   ffbdry,ff2bdry,kffbdry,kff2bdry, &
                   wwbdry,ww2bdry,kwwbdry,kww2bdry, &
                   keefnc,keeknt,eeknt,eetens, &
                   eebdry,ee2bdry,keebdry,kee2bdry
      namelist/inwant/psiwant,vzeroj
      namelist/chiout/saisil,saimpi,saiip
      namelist/eccd/kkstark,chigamt,chigam,bzmse,psiecn,dpsiecn, &
              saisq,cjeccd
      character eqdsk*72,header*42,wform*20,let,fit_type*3
      character*10 vers(6)
      real*8,dimension(:),allocatable :: workk,dmion,bworm,cworm,dworm
!
      allocate(workk(nw),dmion(nw),bworm(nmass),cworm(nmass), &
               dworm(nmass),coils(nsilop),expmp2(magpri),prexp(nrogow))

      allocate(tgamma(nstark),sgamma(nstark),rrrgam(nstark), &
               zzzgam(nstark),aa1gam(nstark),aa2gam(nstark), &
               aa3gam(nstark),aa4gam(nstark),aa5gam(nstark), &
               aa6gam(nstark),aa7gam(nstark),aa8gam(nstark), &
               tgammauncor(nstark))
      allocate(bmsels(nmsels),sbmsels(nmsels),fwtbmsels(nmsels), &
         rrmsels(nmsels),zzmsels(nmsels),l1msels(nmsels),l2msels(nmsels), &
         l4msels(nmsels),emsels(nmsels),semsels(nmsels))
!
      xdum = 0
      if (keqdsk.eq.-1.or.keqdsk.eq.2) return
      if (kdot.gt.0.and.jtime.ne.kdot+1) return
      ijtime=time(jtime)
      ! zero unused variables before write (for consistency)
      if (keecur.le.0) then
        keefnc=0
        keeknt=0
        eeknt=0.0
      endif
!
      if (keqdsk.eq.1) then
        wform='formatted'
      else
        wform='unformatted'
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
!
      plasma=cpasma(jtime)
      btor=bcentr(jtime)
      coils(1:nsilop)=csilop(1:nsilop,jtime)
      expmp2(1:magpri)=cmpr2(1:magpri,jtime)
      prexp(1:nrogow)=crogow(1:nrogow,jtime)
      nbsave=nbdry
      nbabs=nbbbs/min(mbdry,nbdrymx)+1
      jb=0
      do i=1,nbbbs,nbabs
        jb=jb+1
        rbdry(jb)=rbbbs(i)
        zbdry(jb)=zbbbs(i)
      enddo
      nbdry=jb
      if(kdopre.eq.0) pressr(1:npress)=precal(1:npress)
!----------------------------------------------------------------------
!--   set up the ion mass density profile if available               --
!----------------------------------------------------------------------
      if (nmass.gt.0) then
        call zpline(nmass,sibeam,dmass,bworm,cworm,dworm)
        do i=1,nw
          xn=real(i-1,dp)/(nw-1)
          dmion(i)=seval(nmass,xn,sibeam,dmass,bworm,cworm,dworm)
        enddo
      endif
!
      xdim=rgrid(nw)-rgrid(1)
      zdim=zgrid(nh)-zgrid(1)
      zmid=(zgrid(1)+zgrid(nh))/2.0
      do i=1,nw
        do j=1,nh
          kk=(i-1)*nh+j
          if (ivacum.eq.0) then
            if (pasmat(jtime).gt.0.0) then
              psirz(i,j)=-psi(kk)
            else
              psirz(i,j)=psi(kk)
            endif
          else
            psirz(i,j)=-psi(kk)
          endif
        enddo
      enddo
      write (vers(1),1040)
      write (vers(2),1050) efitdate(1:5)
      write (vers(3),1060) efitdate(6:10)
      if (ishot.le.99999) then
        write (vers(4),1070) ishot
      else
        write (vers(4),1073) ishot
      endif
      write (vers(5),1080) ijtime
      vers(6)=' '
      let = 'g'
      call setfnmeq(itimeu,let,ishot,ijtime,eqdsk)
!----------------------------------------------------------------------
!--   If (ISTORE = 1) Then                                           --
!--   Central directory to collect EFIT results is store_dir         --
!----------------------------------------------------------------------
      if (istore .eq. 1) then
        eqdsk = store_dir(1:lstdir)//eqdsk
      endif
      open(unit=neqdsk,file=eqdsk,status='old', &
           form=wform,iostat=ioerr)
      if (ioerr.eq.0) close(unit=neqdsk,status='delete')
      open(unit=neqdsk,file=eqdsk,status='new', &
           form=wform,delim='quote')
      if (pasmat(jtime).gt.0.0) then
        ssimag=-simag
        ssibry=-psibry
      else
        ssimag=simag
        ssibry=psibry
      endif
      eqdsk_format: if (keqdsk.eq.1) then
      write (neqdsk,2000) (vers(i),i=1,6),0,nw,nh
      write (neqdsk,2020) xdim,zdim,rzero,rgrid(1),zmid
      write (neqdsk,2020) rmaxis,zmaxis,ssimag,ssibry,bcentr(jtime)
      write (neqdsk,2020) cpasma(jtime),ssimag,xdum,rmaxis,xdum
      write (neqdsk,2020) zmaxis,xdum,ssibry,xdum,xdum
      write (neqdsk,2020) (fpol(i),i=1,nw)
      write (neqdsk,2020) (pres(i),i=1,nw)
      if (pasmat(jtime).gt.0.0) then
        workk=-ffprim
      else
        workk=ffprim
      endif
      write (neqdsk,2020) (workk(i),i=1,nw)
      if (pasmat(jtime).gt.0.0) then
        workk=-pprime
      else
        workk=pprime
      endif
      write (neqdsk,2020) (workk(i),i=1,nw)
      write (neqdsk,2020) ((psirz(i,j),i=1,nw),j=1,nh)
      write (neqdsk,2020) (qpsi(i),i=1,nw)
      write (neqdsk,2022) nbbbs,limitr
      write (neqdsk,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)
      write (neqdsk,2020) (xlim(i),ylim(i),i=1,limitr)
!----------------------------------------------------------------------
!--   write out rotation information                                 --
!----------------------------------------------------------------------
      write (neqdsk,2024) kvtor,rvtor,nmass
      if (kvtor.gt.0) then
        write (neqdsk,2020) (pressw(i),i=1,nw)
        if (pasmat(jtime).gt.0.0) then
          workk=-pwprim
        else
          workk=pwprim
        endif
        write (neqdsk,2020) (workk(i),i=1,nw)
      endif
!----------------------------------------------------------------------
!--   write out ion mass density profile if available                --
!----------------------------------------------------------------------
      if (nmass.gt.0) then
        write (neqdsk,2020) (dmion(i),i=1,nw)
      endif
      write (neqdsk,2020) (rhovn(i),i=1,nw)
      write (neqdsk,2026) keecur
      if (keecur.gt.0) then
        if (pasmat(jtime).gt.0.0) then
          workk=-epoten
        else
          workk=epoten
        endif
        write (neqdsk,2020) (workk(i),i=1,nw)
      endif
!
!jal 4/23/2004
      if (iplcout.gt.0) then
        if (ishot.le.99999) then
        write (neqdsk,3000) nw,nh,ishot,itime
        else
          write (neqdsk,3003) nw,nh,ishot,itime
        endif
        write (neqdsk,2020) rgrid(1),rgrid(nw),zgrid(1),zgrid(nh)
        write (neqdsk,2020) (brsp(i),i=1,nfcoil)
        write (neqdsk,2020) (ecurrt(i),i=1,nesum)
        write (neqdsk,2020) (pcurrt(i),i=1,nwnh)
      endif
!---------------------------------------------------------------------
!--   Append SNAP file                                              --
!---------------------------------------------------------------------
      if (appendsnap.eq.'G'.or.appendsnap.eq.'KG') then
        if (snapfile/='none') then
          open(unit=nsnapf,status='old', &
               file=snapfile,iostat=ioerr)
          if (ioerr.eq.0) then
            do i=1,1000000
              read (nsnapf,9991,iostat=ioerr) tmpdata
              if (ioerr.ne.0) exit
              if (INDEX(tmpdata,'&efitin')/=0) exit
            enddo
          endif
          if (ioerr.eq.0) then
            do i=1,1000000
              write (neqdsk,9991) tmpdata
              read (nsnapf,9991,iostat=ioerr) tmpdata
              if (ioerr.ne.0) exit
              if (INDEX(tmpdata,'/')/=0) then
                write (neqdsk,9991) tmpdata
                exit
              endif
            enddo
            close (unit=nsnapf)
          endif
 9991     format (a)
        endif
      endif
!
      nqpsi=nw
      limitr=limitr-1
      write (neqdsk,out1)

      call ppstore
      call ffstore
      call wwstore
      call eestore
      write (neqdsk,basis)
      limitr=limitr+1
      MSE: if (kdomse.gt.0.or.kstark.gt.0) then
        if (kdomse.gt.0) then
          tgamma(1:nstark)=cmgam(1:nstark,jtime)
          tgammauncor(1:nstark)=tgamma(1:nstark)
          sgamma(1:nstark)=tgamma(1:nstark)*0.05_dp
        else
          tgamma(1:nstark)=tangam(jtime,1:nstark)
          tgammauncor(1:nstark)=tangam_uncor(jtime,1:nstark)
          sgamma(1:nstark)=siggam(jtime,1:nstark)
        endif
        rrrgam(1:nstark)=rrgam(jtime,1:nstark)
        zzzgam(1:nstark)=zzgam(jtime,1:nstark)
        aa1gam(1:nstark)=a1gam(jtime,1:nstark)
        aa2gam(1:nstark)=a2gam(jtime,1:nstark)
        aa3gam(1:nstark)=a3gam(jtime,1:nstark)
        aa4gam(1:nstark)=a4gam(jtime,1:nstark)
        aa5gam(1:nstark)=a5gam(jtime,1:nstark)
        aa6gam(1:nstark)=a6gam(jtime,1:nstark)
        aa7gam(1:nstark)=a7gam(jtime,1:nstark)
        aa8gam(1:nstark)=a8gam(jtime,1:nstark)
        write (neqdsk,mseout)
        kkstark=nstark
        saisq=tsaisq(jtime)
        write (neqdsk,eccd)
      endif MSE
!-----------------------------------------------------------------------
!--   Write out MSE-LS namelist                                       --
!-----------------------------------------------------------------------
      if (mmbmsels.gt.0.or.kdomsels.gt.0) then         
        if (kdomsels.gt.0) then
          bmsels(1:nmsels)=cmmls(jtime,1:nmsels)
          sbmsels(1:nmsels)=0.05_dp*bmsels(1:nmsels)
        else
          bmsels(1:nmsels)=bmselt(jtime,1:nmsels)
          sbmsels(1:nmsels)=sbmselt(jtime,1:nmsels)
        endif
        rrmsels(1:nmsels)=rrmselt(jtime,1:nmsels)
        zzmsels(1:nmsels)=zzmselt(jtime,1:nmsels)
        l1msels(1:nmsels)=l1mselt(jtime,1:nmsels)
        l2msels(1:nmsels)=l2mselt(jtime,1:nmsels)
        l4msels(1:nmsels)=l4mselt(jtime,1:nmsels)
        emsels(1:nmsels)=emselt(jtime,1:nmsels)
        semsels(1:nmsels)=semselt(jtime,1:nmsels)
        write (neqdsk,in_msels)
      endif
!
      if (kdovt.gt.0) then
        npresw=0
        ndel=nw/26
        do i=1,nw,ndel
          npresw=npresw+1
          presw(npresw)=pressw(i)
          sigprw(npresw)=0.1_dp*presw(npresw)
          rpresw(npresw)=-real(i-1,dp)/(nw-1)
        enddo
        write (neqdsk,vtout)
      endif
      write (neqdsk,chiout)
      header = ' '
      write (neqdsk,1042) header,fit_type
!-----------------------------------------------------------------------
!--   binary format                                                   --
!-----------------------------------------------------------------------
      else eqdsk_format
      write (neqdsk) (vers(i),i=1,6),0,nw,nh
      write (neqdsk) real(xdim,r4),real(zdim,r4),real(rzero,r4), &
                     real(rgrid(1),r4),real(zmid,r4)
      write (neqdsk) real(rmaxis,r4),real(zmaxis,r4),real(ssimag,r4), &
                     real(ssibry,r4),real(bcentr(jtime),r4)
      write (neqdsk) real(cpasma(jtime),r4),real(ssimag,r4),real(xdum,r4), &
                     real(rmaxis,r4),real(xdum,r4)
      write (neqdsk) real(zmaxis,r4),real(xdum,r4),real(ssibry,r4), &
                     real(xdum,r4),real(xdum,r4)
      write (neqdsk) (real(fpol(i),r4),i=1,nw)
      write (neqdsk) (real(pres(i),r4),i=1,nw)
      if (pasmat(jtime).gt.0.0) then
        workk=-ffprim
      else
        workk=ffprim
      endif
      write (neqdsk) (real(workk(i),r4),i=1,nw)
      if (pasmat(jtime).gt.0.0) then
        workk=-pprime
      else
        workk=pprime
      endif
      write (neqdsk) (real(workk(i),r4),i=1,nw)
      write (neqdsk) ((real(psirz(i,j),r4),i=1,nw),j=1,nh)
      write (neqdsk) (real(qpsi(i),r4),i=1,nw)
      write (neqdsk) nbbbs,limitr
      write (neqdsk) (real(rbbbs(i),r4),real(zbbbs(i),r4),i=1,nbbbs)
      write (neqdsk) (real(xlim(i),r4),real(ylim(i),r4),i=1,limitr)
!----------------------------------------------------------------------
!--   write out rotation information                                 --
!----------------------------------------------------------------------
      write (neqdsk) kvtor,real(rvtor,r4),nmass
      if (nmass.gt.0) then
        write (neqdsk) (real(pressw(i),r4),i=1,nw)
        write (neqdsk) (real(pwprim(i),r4),i=1,nw)
      endif
!----------------------------------------------------------------------
!--   write out ion mass density profile if available                --
!----------------------------------------------------------------------
      if (nmass.gt.0) then
        write (neqdsk) (real(dmion(i),r4),i=1,nw)
      endif
      write (neqdsk) (real(rhovn(i),r4),i=1,nw)
      write (neqdsk) keecur
      if (keecur.gt.0) then
        if (pasmat(jtime).gt.0.0) then
          workk=-epoten
        else
          workk=epoten
        endif
        write (neqdsk) (real(workk(i),r4),i=1,nw)
      endif
!jal 4/23/2004
      if (iplcout.gt.0) then
        if (ishot.le.99999) then
          write (neqdsk) nw,nh,ishot,itime
        else
          write (neqdsk) nw,nh,ishot,itime
        endif
        write (neqdsk) real(rgrid(1),r4),real(rgrid(nw),r4), &
                       real(zgrid(1),r4),real(zgrid(nh),r4)
        write (neqdsk) (real(brsp(i),r4),i=1,nfcoil)
        write (neqdsk) (real(ecurrt(i),r4),i=1,nesum)
        write (neqdsk) (real(pcurrt(i),r4),i=1,nwnh)
      endif
!
      header = ' '
      write (neqdsk) header,fit_type
      endif eqdsk_format
!
      close(unit=neqdsk)
      nbdry=nbsave
      fixed_bdry: if (nbdry.gt.2) then
        ierold=ierchk
        ierchk=0
        ibunmn=3
        fwtcur=1.
        mxiter=1
        nxiter=99
        error=1.0e-04_dp
        enps=enp
        enp=0.5_dp
        fwtsi(1:nfcoil)=1.
        fwtsi((nfcoil+1):nsilop)=0.0
        itime=itime+1
        limitr=limitr-1
        let = 'x'
        call setfnmeq(itimeu,let,ishot,itime,eqdsk)
        open(unit=neqdsk,file=eqdsk,status='old',iostat=ioerr)
        if (ioerr.eq.0) close(unit=neqdsk,status='delete')
        open(unit=neqdsk,file=eqdsk,status='new',delim='quote')
        write (neqdsk,in1)
        call ppstore
        call ffstore
        call wwstore
        call eestore
        write (neqdsk,basis)
        write (neqdsk,inwant)
        close(unit=neqdsk)
        enp=enps
        limitr=limitr+1
        itime=itime-1
        ierchk=ierold
      endif fixed_bdry
!
      deallocate(workk,dmion,bworm,cworm,dworm)
!
      return
 1020 format ('x',i5,'.',i3)
 1040 format ('EFIT-AI ')
 1042 format (1x,a42,1x,a3)
 1050 format ('   ',a5)
 1060 format (a5,'   ')
 1070 format (' # ',i5)
 1073 format (' #',i6)
 1080 format ('  ',i4,2hms)
 2000 format (6a8,3i4)
 2020 format (5e16.9)
 2022 format (2i5)
 2024 format (i5,e16.9,i5)
 2026 format (i5)
 3000 format (4i5)
 3003 format (2i5,i6,i5)
      end subroutine weqdsk
