      subroutine shipit(ktime,ifirsttime,ilasttime)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          shipit writes out the plasma shape parameters.          **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          29/06/83..........first created                         **
!**          24/07/85..........revised                               **
!**          24/07/96 revised  by Q.Peng to add the surface area of  **
!**                            the last closed flux surface (psurfa) **
!**                            boundary R when Z=0.0 (rmidin/out)    **
!**                            to write real*4 to unformatted A_eqdsk**
!**          15/05/97 Q.P. added nsilop,magpri,nfcoil,nesum to A     **
!**          23/04/04..JAL iplcout added instead of iecurr=2 print   **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
      include 'basiscomdu.f90'
      dimension jflag(ntime)
      dimension coils(nsilop),expmp2(magpri)
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
      character*10 uday
      character eqdsk*72,header*42,qmflag*3,fit_type*3
      character wform*20,let
!     data nlold/39/,nlnew/40/
      data nlold/40/,nlnew/41/
      save xdum
      xdum = 0
!
!     if (ivacum.gt.0) return
!
      cprof=icprof
      iyes=0
      idup=0
      tavem=2*iavem
      do 300 i=ifirsttime,ilasttime
        jflag(i)=1
  300 continue
      uday = ' '
      call date(uday)
      rcencm=rcentr*100.
      ktime1=1
      mco2v=nco2v
      mco2r=nco2r
!----------------------------------------------------------------------
!-- set fit type                                                     --
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
!-- set fit type equilibrium mode                                    --
!----------------------------------------------------------------------
        if (iconvr.eq.3) then
          fit_type = 'EQU'
          if (npsi_ext.gt.0) fit_type='EQG'
        endif
      do 500 i=ifirsttime,ilasttime
        ijtime=time(i)
        jj=i
!----------------------------------------------------------------------
!-- check error, on return LFLAG > 0 for error and the type of error --
!----------------------------------------------------------------------
        lflag=0
        jj=i
        if (ierchk.gt.0) call chkerr(jj)
        if (lflag.le.0) go to 9306
        jflag(jj)=0
 9306   continue
        if (keqdsk.ge.1) then
         wform='formatted'
        else
         wform='unformatted'
        endif
!
        let = 'a'
        call getfnmu(itimeu,let,ishot,ijtime,eqdsk)
!----------------------------------------------------------------------
!-- If (ISTORE = 1) Then       --
!-- Central directory to collect EFIT results is store_dir  --
!----------------------------------------------------------------------
     if (istore .eq. 1) then
        eqdsk = store_dir(1:lstdir)//eqdsk
     endif
        call db_header(ishot,ijtime,header)
  305   open(unit=neqdsk,file=eqdsk,status='old', &
                  form=wform,err=12929)
             close(unit=neqdsk,status='delete')
12929      continue
        open(unit=neqdsk,                       file=eqdsk,status='new', &
             form=wform)
!
  310   continue
        xbetapr=kbetapr
        if (keqdsk.ge.1) then
        write (neqdsk,1055) uday,(mfvers(j),j=1,2)
        if (ishot.le.99999) then
        write (neqdsk,1050) ishot,ktime1
        else
        write (neqdsk,1053) ishot,ktime1
        endif
        write (neqdsk,1040) (time(j),j=i,i)
        ltime=i
!
        jj=i
  320   continue
        if (fwtqa.gt.0.0.and.qvfit.gt.0.0) qmflag='FIX'
        if (fwtqa.le.0.0) qmflag='CLC'
!
        write (neqdsk,1060) time(jj),jflag(jj),lflag,limloc(jj), &
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
        fluxx=diamag(jj)*1.0e-03
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
!-- one xxx replaced with ring gap 8/15/91                        --
!--     cjor0=flux surface average current density normalized to      --
!--           I/A at magnetic axis                                    --
!--     cjor95=jor95/jor0, cjor99=jor99/(I/A)                         --
!--     xbetapr=kbetapr                        2014/06/25             --
!-----------------------------------------------------------------------
        write (neqdsk,1042) header,fit_type
!-------------------------------------------------------------------
!-- binary format                                                 --
!-------------------------------------------------------------------
        else
        write (neqdsk) uday,(mfvers(j),j=1,2)
        if (ishot.le.99999) then
        write (neqdsk) ishot,ktime1
        else
        write (neqdsk) ishot,ktime1
        endif
        write (neqdsk) (time(j),j=i,i)
        ltime=i
!
        jj=i
30320   continue
        if (fwtqa.gt.0.0.and.qvfit.gt.0.0) qmflag='FIX'
        if (fwtqa.le.0.0) qmflag='CLC'
!
        write (neqdsk) real(time(jj)),jflag(jj),lflag,limloc(jj), &
                            mco2v,mco2r,qmflag,nlold,nlnew
        write (neqdsk) real(tsaisq(jj)),real(rcencm), &
                       real(bcentr(jj)),real(pasmat(jj))
        write (neqdsk) real(cpasma(jj)),real(rout(jj)), &
                       real(zout(jj)),real(aout(jj))
        write (neqdsk) real(eout(jj)),real(doutu(jj)), &
                       real(doutl(jj)),real(vout(jj))
        write (neqdsk) real(rcurrt(jj)),real(zcurrt(jj)), &
                       real(qsta(jj)),real(betat(jj))
        write (neqdsk) real(betap(jj)),real(ali(jj)), &
                       real(oleft(jj)),real(oright(jj))
        write (neqdsk) real(otop(jj)),real(obott(jj)), &
                       real(qpsib(jj)),real(vertn(jj))
        write (neqdsk) (real(rco2v(k,jj)),k=1,mco2v)
        write (neqdsk) (real(dco2v(jj,k)),k=1,mco2v)
        write (neqdsk) (real(rco2r(k,jj)),k=1,mco2r)
        write (neqdsk) (real(dco2r(jj,k)),k=1,mco2r)
        write (neqdsk) real(shearb(jj)),real(bpolav(jj)), &
                       real(s1(jj)),real(s2(jj))
        write (neqdsk) real(s3(jj)),real(qout(jj)), &
                       real(olefs(jj)),real(orighs(jj))
        write (neqdsk) real(otops(jj)),real(sibdry(jj)), &
                       real(areao(jj)),real(wplasm(jj))
        write (neqdsk) real(terror(jj)),real(elongm(jj)), &
                       real(qqmagx(jj)),real(cdflux(jj))
        write (neqdsk) real(alpha(jj)),real(rttt(jj)), &
                       real(psiref(jj)),real(xndnt(jj))
        write (neqdsk) real(rseps(1,jj)),real(zseps(1,jj)), &
                       real(rseps(2,jj)),real(zseps(2,jj))
        write (neqdsk) real(sepexp(jj)),real(obots(jj)), &
                       real(btaxp(jj)),real(btaxv(jj))
        write (neqdsk) real(aaq1(jj)),real(aaq2(jj)), &
                       real(aaq3(jj)),real(seplim(jj))
        write (neqdsk) real(rmagx(jj)),real(zmagx(jj)), &
                       real(simagx(jj)),real(taumhd(jj))
        fluxx=diamag(jj)*1.0e-03
        write (neqdsk) real(betapd(jj)),real(betatd(jj)), &
                       real(wplasmd(jj)),real(fluxx)
        write (neqdsk) real(vloopt(jj)),real(taudia(jj)), &
                       real(qmerci(jj)),real(tavem)
!
        if (ishot.lt.91000) then
        write (neqdsk) nsilop,magpri67+magpri322,nfcoil,nesum
        write (neqdsk) (real(csilop(k,jj)),k=1,nsilop), &
                       (real(cmpr2(k,jj)),k=1,magpri67+magpri322)
        else
        write (neqdsk) nsilop,magpri,nfcoil,nesum
        write (neqdsk) (real(csilop(k,jj)),k=1,nsilop), &
                       (real(cmpr2(k,jj)),k=1,magpri)
        endif
        write (neqdsk) (real(ccbrsp(k,jj)),k=1,nfcoil)
        write (neqdsk) (real(eccurt(jj,k)),k=1,nesum)
!
        write (neqdsk) real(pbinj(jj)),real(rvsin(jj)), &
                       real(zvsin(jj)),real(rvsout(jj))
        write (neqdsk) real(zvsout(jj)),real(vsurfa(jj)), &
                       real(wpdot(jj)),real(wbdot(jj))
        write (neqdsk) real(slantu(jj)),real(slantl(jj)), &
                       real(zuperts(jj)),real(chipre)
        write (neqdsk) real(cjor95(jj)),real(pp95(jj)), &
                       real(ssep(jj)),real(yyy2(jj))
        write (neqdsk) real(xnnc(jj)),real(cprof), &
                       real(oring(jj)),real(cjor0(jj))
        write (neqdsk) real(fexpan),real(qqmin), &
                       real(chigamt),real(ssi01)
        write (neqdsk) real(fexpvs),real(sepnose), &
                       real(ssi95(jj)),real(rqqmin)
        write (neqdsk) real(cjor99(jj)),real(cj1ave(jj)), &
                       real(rmidin(jj)),real(rmidout(jj))
        write (neqdsk) real(psurfa(jj)),real(peak(jj)), &
                       real(dminux(jj)),real(dminlx(jj))
        write (neqdsk) real(dolubaf(jj)),real(dolubafm(jj)), &
                       real(diludom(jj)),real(diludomm(jj))
        write (neqdsk) real(ratsol(jj)),real(rvsiu(jj)), &
                       real(zvsiu(jj)),real(rvsid(jj))
        write (neqdsk) real(zvsid(jj)),real(rvsou(jj)), &
                       real(zvsou(jj)),real(rvsod(jj))
        write (neqdsk) real(zvsod(jj)),real(condno), &
                       real(psin32(jj)),real(psin21(jj))
        write (neqdsk) real(rq32in(jj)),real(rq21top(jj)), &
                       real(chilibt),real(ali3(jj))
        write (neqdsk) real(xbetapr),real(tflux(jj)),real(tchimls),real(twagap(jj))
!
        write (neqdsk) header,fit_type
        endif
!
        close(unit=neqdsk)
  500 continue
!
! --- add flag for writting out esave.dat. IOUT contains 16.
!
      if ((ifirsttime.eq.ktime).and.(iand(iout,16).ne.0)) then
      open(unit=nsave,status='old',form='unformatted',file='esave.dat', &
           err=12930)
      close(unit=nsave,status='delete')
12930  continue
      open(unit=nsave,status='new',form='unformatted',file='esave.dat')
      mw=nw
      mh=nh
      write (nsave) mw,mh
      write (nsave) xpsi
      write (nsave) brsp
      write (nsave) www
      write (nsave) emaxis, rmaxis, fcentr
      close(unit=nsave)
      endif
!
      if (kinput.le.0) go to 700
      do 695 ij=ifirsttime,ilasttime
      jtime=ij
      plasma=cpasma(jtime)
      btor=bcentr(jtime)
      do 600 i=1,nsilop
        coils(i)=csilop(i,jtime)
  600 continue
      do 630 i=1,magpri
        expmp2(i)=cmpr2(i,jtime)
  630 continue
      ierold=ierchk
      ierchk=0
      fwtcur=1.
      do 670 i=1,magpri
        fwtmp2(i)=1.
  670 continue
      do 690 i=1,nsilop
        fwtsi(i)=1.
  690 continue
      iitime=time(ij)
      iitime=iitime+1
      limitr=limitr-1
      qvfit=qenp
      let = 'm'
      call getfnmu(itimeu,let,ishot,iitime,eqdsk)
      open(unit=neqdsk,file=eqdsk,status='old',err=12931)
      close(unit=neqdsk,status='delete')
12931  continue
      open(unit=neqdsk,                       file=eqdsk,status='new')
      itsave=itime
      itime=iitime
      write (neqdsk,in1)
      call ppstore
      call ffstore
      call wwstore
      call eestore
      write (neqdsk,basis)
      write (neqdsk,inwant)
      if (kvtor.gt.0) then
        write (neqdsk,invt)
      endif
      itime=itsave
      close(unit=neqdsk)
      limitr=limitr+1
      ierchk=ierold
  695 continue
  700 continue
      return
!
 1020 format ('p',i5,'.',i3)
 1040 format (1x,4e16.9)
 1041 format (1x,4i5)
 1042 format (1x,a42,1x,a3)
 1050 format (1x,i5,11x,i5)
 1053 format (1x,i6,11x,i5)
 1055 format (1x,a10,2a5)
 1060 format ('*',f8.3,9x,i5,11x,i5,1x,a3,1x,i3,1x,i3,1x,a3,1x,2i5)
 1070 format(16a5)
      end
      subroutine weqdsk(jtime)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          weqdsk writes out the GAQ type eqdsk.                   **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          29/06/83..........first created                         **
!**          24/07/85..........revised                               **
!**          22/07/96 revised  by Q.Peng                             **
!**                            to write real*4 to unformatted G_eqdsk**
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: psirz
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      include 'basiscomdu.f90'
      dimension coils(nsilop),expmp2(magpri),prexp(nrogow)
      namelist/out1/ishot,itime,betap0,rzero,qenp,enp,emp,plasma, &
           expmp2,coils,btor,rcentr,brsp,icurrt,rbdry,zbdry, &
           nbdry,fwtsi,fwtcur,mxiter,nxiter,limitr,xlim,ylim,error, &
           iconvr,ibunmn,pressr,rpress,nqpsi, &
           npress,sigpre
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
      namelist/chiout/saisil,saimpi,saipr,saiip
      namelist/eccd/kkstark,chigamt,chigam,bzmse,psiecn,dpsiecn, &
              saisq,cjeccd
      character eqdsk*72,header*42,wform*20,let,fit_type*3
      character*10 case(6)
      integer :: parameter
      parameter (pltnw=1025)
      real*8,dimension(:),allocatable :: workk,dmion,bworm,cworm, &
             dworm
      real*8 tgamma(nstark),sgamma(nstark),rrrgam(nstark) &
               ,zzzgam(nstark),aa1gam(nstark),aa2gam(nstark) &
               ,aa3gam(nstark),aa4gam(nstark),aa5gam(nstark) &
               ,aa6gam(nstark),aa7gam(nstark),aa8gam(nstark) &
               ,tgammauncor(nstark)
      real*8  bmsels(nmsels),sbmsels(nmsels),fwtbmsels(nmsels), &
         rrmsels(nmsels),zzmsels(nmsels),l1msels(nmsels),l2msels(nmsels), &
         l4msels(nmsels),emsels(nmsels),semsels(nmsels)
!
      ALLOCATE(workk(pltnw),dmion(pltnw),bworm(pltnw),cworm(pltnw),dworm(pltnw))
!
      xdum = 0
      if (keqdsk.eq.-1.or.keqdsk.eq.2) return
      if (kdot.gt.0.and.jtime.ne.kdot+1) return
      ijtime=time(jtime)
      mw=nw
      mh=nh
  100 continue
!
      if (keqdsk.eq.1) then
        wform='formatted'
      else
        wform='unformatted'
      endif
!----------------------------------------------------------------------
!-- set fit type                                                     --
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
      do 500 i=1,nsilop
        coils(i)=csilop(i,jtime)
  500 continue
      do 520 i=1,magpri
        expmp2(i)=cmpr2(i,jtime)
  520 continue
      do 530 i=1,nrogow
        prexp(i)=crogow(i,jtime)
  530 continue
      nbsave=nbdry
      nbabs=nbbbs/min(mbdry,nbdrymx)+1
      jb=0
      do 540 i=1,nbbbs,nbabs
        jb=jb+1
        rbdry(jb)=rbbbs(i)
        zbdry(jb)=zbbbs(i)
  540 continue
      nbdry=jb
      if (kdopre.eq.0) then
      do 541 i=1,npress
        pressr(i)=precal(i)
  541 continue
      endif
!----------------------------------------------------------------------
!-- set up the ion mass density profile if available                 --
!----------------------------------------------------------------------
      if (nmass.gt.0) then
      call zpline(nmass,sibeam,dmass,bworm,cworm,dworm)
      do 93003 i=1,nw
        xn=float(i-1)/float(nw-1)
        dmion(i)=seval(nmass,xn,sibeam,dmass,bworm,cworm,dworm)
93003 continue
      endif
!
      xdim=rgrid(nw)-rgrid(1)
      zdim=zgrid(nh)-zgrid(1)
      zmid=(zgrid(1)+zgrid(nh))/2.0
      do 200 i=1,nw
      do 200 j=1,nh
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
  200 continue
      write (case(1),1040)
      write (case(2),1050) mfvers(1)
      write (case(3),1060) mfvers(2)
      if (ishot.le.99999) then
       write (case(4),1070) ishot
      else
       write (case(4),1073) ishot
      endif
      write (case(5),1080) ijtime
      case(6)=' '
      let = 'g'
      call getfnmu(itimeu,let,ishot,ijtime,eqdsk)
!----------------------------------------------------------------------
!-- If (ISTORE = 1) Then       --
!-- Central directory to collect EFIT results is store_dir  --
!----------------------------------------------------------------------
     if (istore .eq. 1) then
        eqdsk = store_dir(1:lstdir)//eqdsk
     endif
  205 open(unit=neqdsk,file=eqdsk,status='old', &
           form=wform,err=12932)
           close(unit=neqdsk,status='delete')
12932   continue
      open(unit=neqdsk,file=eqdsk,status='new', &
           form=wform)
  210 continue
      idum=3
      if (pasmat(jtime).gt.0.0) then
        ssimag=-simag
        ssibry=-psibry
      else
        ssimag=simag
        ssibry=psibry
      endif
      if (keqdsk.eq.1) then
      write (neqdsk,2000) (case(i),i=1,6),idum,mw,mh
      write (neqdsk,2020) xdim,zdim,rzero,rgrid(1),zmid
      write (neqdsk,2020) rmaxis,zmaxis,ssimag,ssibry,bcentr(jtime)
      write (neqdsk,2020) cpasma(jtime),ssimag,xdum,rmaxis,xdum
      write (neqdsk,2020) zmaxis,xdum,ssibry,xdum,xdum
      write (neqdsk,2020) (fpol(i),i=1,nw)
      write (neqdsk,2020) (pres(i),i=1,nw)
      do 310 i=1,nw
      if (pasmat(jtime).gt.0.0) then
        workk(i)=-ffprim(i)
      else
        workk(i)=ffprim(i)
      endif
  310 continue
      write (neqdsk,2020) (workk(i),i=1,nw)
      do 315 i=1,nw
      if (pasmat(jtime).gt.0.0) then
        workk(i)=-pprime(i)
      else
        workk(i)=pprime(i)
      endif
  315 continue
      write (neqdsk,2020) (workk(i),i=1,nw)
      write (neqdsk,2020) ((psirz(i,j),i=1,nw),j=1,nh)
      write (neqdsk,2020) (qpsi(i),i=1,nw)
      write (neqdsk,2022) nbbbs,limitr
      write (neqdsk,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)
      write (neqdsk,2020) (xlim(i),ylim(i),i=1,limitr)
!----------------------------------------------------------------------
!--  write out rotation information                                  --
!----------------------------------------------------------------------
      write (neqdsk,2024) kvtor,rvtor,nmass
      if (kvtor.gt.0) then
        write (neqdsk,2020) (pressw(i),i=1,nw)
        do i=1,nw
          if (pasmat(jtime).gt.0.0) then
            workk(i)=-pwprim(i)
          else
            workk(i)=pwprim(i)
          endif
        enddo
        write (neqdsk,2020) (workk(i),i=1,nw)
      endif
!----------------------------------------------------------------------
!--  write out ion mass density profile if available                 --
!----------------------------------------------------------------------
      if (nmass.gt.0) then
        write (neqdsk,2020) (dmion(i),i=1,nw)
      endif
      write (neqdsk,2020) (rhovn(i),i=1,nw)
      write (neqdsk,2026) keecur
      if (keecur.gt.0) then
        do i=1,nw
          if (pasmat(jtime).gt.0.0) then
            workk(i)=-epoten(i)
          else
            workk(i)=epoten(i)
          endif
        enddo
        write (neqdsk,2020) (workk(i),i=1,nw)
      endif
!
!jal 4/23/2004
      if (iplcout.le.0) go to 350
      if (ishot.le.99999) then
      write (neqdsk,3000) mw,mh,ishot,itime
      else
      write (neqdsk,3003) mw,mh,ishot,itime
      endif
      write (neqdsk,2020) rgrid(1),rgrid(nw),zgrid(1),zgrid(nh)
      write (neqdsk,2020) (brsp(i),i=1,nfcoil)
      write (neqdsk,2020) (ecurrt(i),i=1,nesum)
      write (neqdsk,2020) (pcurrt(i),i=1,nwnh)
  350 continue
!---------------------------------------------------------------------
!-- Append SNAP file                                                --
!---------------------------------------------------------------------
      if (appendsnap.eq.'G'.or.appendsnap.eq.'KG') then
      if (snapfile/='none') then
         open(unit=nsnapf,status='old', &
           file=snapfile,err=9981)
      do i=1,1000000
         read (nsnapf,9991,end=9981) tmpdata
         if (INDEX(tmpdata,'&efitin')/=0) go to 381
      enddo
 381  continue
      do i=1,1000000
         write (neqdsk,9991) tmpdata
         read (nsnapf,9991,end=9979) tmpdata
         if (INDEX(tmpdata,'/')/=0) then
            write (neqdsk,9991) tmpdata
            go to 9979
         endif
      enddo
 9979 close (unit=nsnapf)
 9981 continue
 9991 format (a)
      endif
      endif
!
      nqpsi=nw
      limitr=limitr-1
!      write (neqdsk,out1)
      write(neqdsk,*) '&OUT1'
!       write (neqdsk,*) 'QPSI =', qpsi(1),','
!       do i=2,SIZE(qpsi)
!        write (neqdsk,*) '      ', qpsi(i),','
!       enddo
       write (neqdsk,*) 'ISHOT =',ishot
       write (neqdsk,*) 'ITIME =',itime
       write (neqdsk,*) 'BETAP0 =',betap0
       write (neqdsk,*) 'RZERO =',rzero
       write (neqdsk,*) 'QENP =',qenp
       write (neqdsk,*) 'ENP =',enp
       write (neqdsk,*) 'EMP =',emp
       write (neqdsk,*) 'PLASMA =',plasma
       write (neqdsk,*) 'EXPMP2 =',expmp2
       write (neqdsk,*) 'COILS =',coils
       write (neqdsk,*) 'BTOR =',btor
       write (neqdsk,*) 'RCENTR =',rcentr
       write (neqdsk,*) 'BRSP =',brsp
       write (neqdsk,*) 'ICURRT =',icurrt
       write (neqdsk,*) 'RBDRY =',rbdry
       write (neqdsk,*) 'ZBDRY =',zbdry
       write (neqdsk,*) 'NBDRY =',nbdry
       write (neqdsk,*) 'FWTSI =',fwtsi
       write (neqdsk,*) 'FWTCUR =',fwtcur
       write (neqdsk,*) 'MXITER =',mxiter
       write (neqdsk,*) 'NXITER =',nxiter
       write (neqdsk,*) 'LIMITR =',limitr
       write (neqdsk,*) 'XLIM =',xlim
       write (neqdsk,*) 'YLIM =',ylim
       write (neqdsk,*) 'ERROR =',error
       write (neqdsk,*) 'ICONVR =',iconvr
       write (neqdsk,*) 'IBUNMN =',ibunmn
       write (neqdsk,*) 'PRESSR =',pressr
       write (neqdsk,*) 'RPRESS =',rpress
       write (neqdsk,*) 'QPSI =',qpsi
       write (neqdsk,*) 'PRESSW =',pressw
       write (neqdsk,*) 'PRES =',pres
       write (neqdsk,*) 'NQPSI =',nqpsi
       write (neqdsk,*) 'NPRESS =',npress
       write (neqdsk,*) 'SIGPRE =',sigpre
       write (neqdsk,*) '/'

       call ppstore
       call ffstore
       call wwstore
       call eestore
      write (neqdsk,basis)
      limitr=limitr+1
      if (kdomse.gt.0.or.kstark.gt.0) then
         do i=1,nstark
           if (kdomse.gt.0) then
             tgamma(i)=cmgam(i,jtime)
             tgammauncor(i) = tgamma(i)
             sgamma(i)=tgamma(i)*0.05
           else
             tgamma(i)=tangam(jtime,i)
             tgammauncor(i)=tangam_uncor(jtime,i)
             sgamma(i)=siggam(jtime,i)
           endif
           rrrgam(i)=rrgam(jtime,i)
           zzzgam(i)=zzgam(jtime,i)
           aa1gam(i)=a1gam(jtime,i)
           aa2gam(i)=a2gam(jtime,i)
           aa3gam(i)=a3gam(jtime,i)
           aa4gam(i)=a4gam(jtime,i)
           aa5gam(i)=a5gam(jtime,i)
           aa6gam(i)=a6gam(jtime,i)
           aa7gam(i)=a7gam(jtime,i)
           aa8gam(i)=a8gam(jtime,i)
         enddo
         write (neqdsk,mseout)
         kkstark=nstark
         saisq=tsaisq(jtime)
         write (neqdsk,eccd)
      endif
!-----------------------------------------------------------------------
!--  Write out MSE-LS namelist                                        --
!-----------------------------------------------------------------------
      if (mmbmsels.gt.0.or.kdomsels.gt.0) then         
         do i=1,nmsels
           if (kdomsels.gt.0) then
             bmsels(i)=cmmls(jtime,i)
             sbmsels(i)=0.05*bmsels(i)
           else
             bmsels(i)=bmselt(jtime,i)
             sbmsels(i)=sbmselt(jtime,i)
           endif
           rrmsels(i)=rrmselt(jtime,i)
           zzmsels(i)=zzmselt(jtime,i)
           l1msels(i)=l1mselt(jtime,i)
           l2msels(i)=l2mselt(jtime,i)
           l4msels(i)=l4mselt(jtime,i)
           emsels(i)=emselt(jtime,i)
           semsels(i)=semselt(jtime,i)
         enddo
         write (neqdsk,in_msels)
      endif
!
      if (kdovt.gt.0) then
         npresw=0
         ndel=nw/26
         do i=1,nw,ndel
           npresw=npresw+1
           presw(npresw)=pressw(i)
           sigprw(npresw)=0.1*presw(npresw)
           rpresw(npresw)=-float(i-1)/float(nw-1)
         enddo
         write (neqdsk,vtout)
      endif
      write (neqdsk,chiout)
      call db_header(ishot,itime,header)
      write (neqdsk,1042) header,fit_type
!-----------------------------------------------------------------------
!--   binary format                                                   --
!-----------------------------------------------------------------------
      else
      write (neqdsk) (case(i),i=1,6),idum,mw,mh
      write (neqdsk) real(xdim),real(zdim),real(rzero), &
                     real(rgrid(1)),real(zmid)
      write (neqdsk) real(rmaxis),real(zmaxis),real(ssimag), &
                     real(ssibry),real(bcentr(jtime))
      write (neqdsk) real(cpasma(jtime)),real(ssimag),real(xdum), &
                     real(rmaxis),real(xdum)
      write (neqdsk) real(zmaxis),real(xdum),real(ssibry), &
                     real(xdum),real(xdum)
      write (neqdsk) (real(fpol(i)),i=1,nw)
      write (neqdsk) (real(pres(i)),i=1,nw)
      do 30310 i=1,nw
      if (pasmat(jtime).gt.0.0) then
        workk(i)=-ffprim(i)
      else
        workk(i)=ffprim(i)
      endif
30310 continue
      write (neqdsk) (real(workk(i)),i=1,nw)
      do 30315 i=1,nw
      if (pasmat(jtime).gt.0.0) then
        workk(i)=-pprime(i)
      else
        workk(i)=pprime(i)
      endif
30315 continue
      write (neqdsk) (real(workk(i)),i=1,nw)
      write (neqdsk) ((real(psirz(i,j)),i=1,nw),j=1,nh)
      write (neqdsk) (real(qpsi(i)),i=1,nw)
      write (neqdsk) nbbbs,limitr
      write (neqdsk) (real(rbbbs(i)),real(zbbbs(i)),i=1,nbbbs)
      write (neqdsk) (real(xlim(i)),real(ylim(i)),i=1,limitr)
!----------------------------------------------------------------------
!--  write out rotation information                                  --
!----------------------------------------------------------------------
      write (neqdsk) kvtor,real(rvtor),nmass
      if (nmass.gt.0) then
        write (neqdsk) (real(pressw(i)),i=1,nw)
        write (neqdsk) (real(pwprim(i)),i=1,nw)
      endif
!----------------------------------------------------------------------
!--  write out ion mass density profile if available                 --
!----------------------------------------------------------------------
      if (nmass.gt.0) then
        write (neqdsk) (real(dmion(i)),i=1,nw)
      endif
      write (neqdsk) (real(rhovn(i)),i=1,nw)
      write (neqdsk) keecur
      if (keecur.gt.0) then
        do i=1,nw
          if (pasmat(jtime).gt.0.0) then
            workk(i)=-epoten(i)
          else
            workk(i)=epoten(i)
          endif
        enddo
        write (neqdsk) (real(workk(i)),i=1,nw)
      endif
!jal 4/23/2004
      if (iplcout.le.0) go to 30350
      if (ishot.le.99999) then
      write (neqdsk) mw,mh,ishot,itime
      else
      write (neqdsk) mw,mh,ishot,itime
      endif
      write (neqdsk) real(rgrid(1)),real(rgrid(nw)), &
                     real(zgrid(1)),real(zgrid(nh))
      write (neqdsk) (real(brsp(i)),i=1,nfcoil)
      write (neqdsk) (real(ecurrt(i)),i=1,nesum)
      write (neqdsk) (real(pcurrt(i)),i=1,nwnh)
30350 continue
!
      call db_header(ishot,itime,header)
      write (neqdsk) header,fit_type
      endif
!
      close(unit=neqdsk)
      nbdry=nbsave
      if (nbdry.le.2) go to 600
      ierold=ierchk
      ierchk=0
      ibunmn=3
      fwtcur=1.
      mxiter=1
      nxiter=99
      error=1.0e-04
      enps=enp
      enp=0.5
      do 590 i=1,nsilop
        fwtsi(i)=1.
        if (i.gt.nfcoil) fwtsi(i)=0.0
  590 continue
      itime=itime+1
      limitr=limitr-1
      let = 'x'
      call getfnmu(itimeu,let,ishot,itime,eqdsk)
      open(unit=neqdsk,file=eqdsk,status='old',err=12933)
      close(unit=neqdsk,status='delete')
12933  continue
      open(unit=neqdsk,                       file=eqdsk,status='new')
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
  600 continue
!
      DEALLOCATE(workk,dmion,bworm,cworm,dworm)
!
      return
 1020 format ('x',i5,'.',i3)
 1040 format ('  EFITD ')
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
      end
      subroutine timdot(kdot,time,wplasm,sibdry,wbpol,vloopt, &
                        wpdot,wbdot,vsurfa,pbinj,taumhd,cpasma, &
                        wplasmd,taudia)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          shipit writes out the plasma shape parameters.          **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          12/03/87..........first created                         **
!**                                                                  **
!**********************************************************************
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension time(1),wplasm(1),sibdry(1),wbpol(1),vloopt(1), &
                wpdot(1),wbdot(1),vsurfa(1),pbinj(1),taumhd(1), &
                cpasma(1),wplasmd(1),taudia(1)
      data mychoice/2/,twopi/6.283185/
!
      if (kdot.eq.0) return
      if (mychoice.eq.2) go to 10000
      return
!----------------------------------------------------------------------
!-- derivative by finite differencing                                --
!----------------------------------------------------------------------
10000 continue
      kkkk=kdot+1
      wpdot(kkkk)=0.
      wbdot(kkkk)=0.
      vsurfa(kkkk)=0.
      wpddot=0.0
      do 10050 i=1,kdot
       ip=2*kdot+2-i
       deltt=(time(ip)-time(i))/1000./float(kdot)
       wpdot(kkkk)=wpdot(kkkk)+(wplasm(ip)-wplasm(i))/deltt
       wpddot=wpddot+(wplasmd(ip)-wplasmd(i))/deltt
       wbdot(kkkk)=wbdot(kkkk)+(wbpol(ip)-wbpol(i))/deltt
       vsurfa(kkkk)=vsurfa(kkkk)+(sibdry(ip)-sibdry(i))/deltt
10050 continue
      vsurfa(kkkk)=vloopt(kkkk)-twopi*vsurfa(kkkk)
      ptotal=vsurfa(kkkk)*cpasma(kkkk)-wbdot(kkkk)+pbinj(kkkk)
      taumhd(kkkk)=wplasm(kkkk)/(ptotal-wpdot(kkkk))*1000.
      taudia(kkkk)=wplasmd(kkkk)/(ptotal-wpddot)*1000.

      return
      end

      subroutine wmeasure(ktime,ifirst,ilast,itype)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          Subroutine wmeasure writes measurements (MSE, magnetic) **
!**          and some of their corresponding calculated values in    **
!**          a netCDF file of name m0sssss.nc/m0sssss.0tttt.         **
!**                                                                  **
!**          The input variable IOUT decides the mode,               **
!**          IOUT.and.2.ne.0 - one file for all slices m0sssss.nc    **
!**          IOUT.and.4.ne.0 - one file for each slice m0sssss.0tttt **
!**          IOUT.and.6.ne.0 - one file for all m0sssss.0tttt_nnnn   **
!**               tttt as the first time,nnnn as the number of slices**
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**          ktime - total number of time slices                     **
!**          ifirst- index of starting time                          **
!**          ilast - index of ending time                            **
!**          itype - 1 called from main routine with time loop       **
!**                  2 called from main routine out of time loop,    **
!**                    write all slices at one time                  **
!**                                                                  **
!**     REFERENCES:                                                  **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          12-10-97 created by Q.Peng                              **
!**          15-06-98 Q.Peng added measured pressure as premea.      **
!**          01-11-98 Q.Peng added mseport (vport in stark_multi)    **
!**          01-20-2000 Q.P. following the convention of A & G files **
!**                          use mssssss.0tttt_nnn for sub-millisec. **
!**                                                                  **
!**********************************************************************
!
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
      include 'netcdf.inc'   ! from the netCDF package..
!                            ..this must be symlinked to local directory
! --- common block ids is local to wmeasure
!
      common/ids/nceq,idim_time,idim_1,idim_nstark,idim_nsilop &
      ,idim_magpri,idim_nfcoil,idim_nesum,idim_npress,idim_npcurn &
      ,idim_nitera &
      ,id_shot,id_time,id_tangam,id_cmgam,id_siggam,id_rrgam,id_zzgam &
      ,id_a1gam,id_a2gam,id_a3gam,id_a4gam,id_a5gam,id_a6gam,id_a7gam &
      ,id_a8gam,id_fgam, id_tangam_uncor &
      ,id_fwtgam,id_msebkp,id_mseport &
      ,id_silopt,id_csilop,id_fwtsi,id_expmpi,id_fwtmp2,id_cmpr2 &
      ,id_fccurt,id_ccbrsp,id_fwtfc,id_eccurt,id_cecurr,id_fwtec &
      ,id_diamag,id_cdflux,id_sigdia,id_plasma,id_cpasma &
      ,id_pressr,id_cpress,id_rpress,id_zpress,id_sigpre,id_saipre &
      ,id_chigam,id_saimpi,id_saisil,id_czmaxi,id_cchisq,id_cerror &
      ,id_fungam,id_gaingam,id_slopegam,id_scalegam,id_offsetgam
      dimension xrsp(npcurn)
!sri  to make it compatiable to mpi version
!      character eqdsk*20,let,title*80,last*3
      character eqdsk*80,let,title*80,last*3
      character*30 sfname
      integer dim2(2),c11(2),cnn(2),imap(2),stride(2)
! --- temporay variables to convert double to single
      real*4 zwork(ntime+nsilop+nstark+nfcoil+nesum+magpri+npress) &
      ,zcmgam(nstark,ntime),zcsilop(nsilop,ntime) &
      ,zcmpr2(magpri,ntime),zccbrsp(nfcoil,ntime),zstark(ntime,nstark) &
      ,zsilopt(ntime,nsilop),zexpmpi(ntime,magpri) &
      ,zfccurt(ntime,nfcoil),zeccurt(ntime,nesum) 
!vas f90 modifi
      character*85 vastext
!-----------------------------------------------------------------------
!--  write out t(shot).(time)_X files                                    --
!-----------------------------------------------------------------------
      if ((kwripre.eq.11).and.(itype.eq.2)) then
        call getfnmd('t',ishot,itime,sfname)
        sfname=sfname(1:13)//'_chi2'
          open(unit=74,status='old',file=sfname,err=19920)
          close(unit=74,status='delete')
19920     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirst,ilast
          write (74,*) time(i),tsaisq(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_error'
          open(unit=74,status='old',file=sfname,err=19922)
          close(unit=74,status='delete')
19922     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirst,ilast
          write (74,*) time(i),terror(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_j1ave'
          open(unit=74,status='old',file=sfname,err=19930)
          close(unit=74,status='delete')
19930     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirst,ilast
          write (74,*) time(i),cj1ave(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_li'
          open(unit=74,status='old',file=sfname,err=19940)
          close(unit=74,status='delete')
19940     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirst,ilast
          write (74,*) time(i),ali(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_betat'
          open(unit=74,status='old',file=sfname,err=19950)
          close(unit=74,status='delete')
19950     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirst,ilast
          write (74,*) time(i),betat(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_q95'
          open(unit=74,status='old',file=sfname,err=19960)
          close(unit=74,status='delete')
19960     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirst,ilast
          write (74,*) time(i),qpsib(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_q0'
          open(unit=74,status='old',file=sfname,err=19970)
          close(unit=74,status='delete')
19970     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirst,ilast
          write (74,*) time(i),qqmagx(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_q0'
          open(unit=74,status='old',file=sfname,err=19980)
          close(unit=74,status='delete')
19980     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirst,ilast
          write (74,*) time(i),qqmagx(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_eout'
          open(unit=74,status='old',file=sfname,err=19990)
          close(unit=74,status='delete')
19990     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirst,ilast
          write (74,*) time(i),eout(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_vout'
          open(unit=74,status='old',file=sfname,err=19992)
          close(unit=74,status='delete')
19992     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirst,ilast
          vm3=vout(i)/1.e6
          write (74,*) time(i),vm3,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_betan'
          open(unit=74,status='old',file=sfname,err=19995)
          close(unit=74,status='delete')
19995     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirst,ilast
          pasman=cpasma(i)/1.e4/aout(i)/abs(bcentr(i))
          pasman=pasman*rout(i)/100./rcentr
          betatnx=betat(i)/pasman
          write (74,*) time(i),betatnx,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_zts'
          open(unit=74,status='old',file=sfname,err=19998)
          close(unit=74,status='delete')
19998     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirst,ilast
          write (74,*) time(i),zuperts(i),xdum,xdum
        enddo
        close(unit=74)
      endif
!
      if ((iand(iout,2).eq.0).and.(iand(iout,4).eq.0)) &
           return
!
      call NCPOPT (NCVERBOS)              ! get nc error message but not fatal
!
! --- time-dependent and the first slice and first time called, 
! --- create file and enter define mode
! --- if iout.and.2.ne.0 = create m0sssss.nc
! --- if iout.and.6.ne.0 = create m0sssss.0tttt_nnnn, 
! ---    while tttt is the first time and nnn is the number of slices.

      if ((iand(iout,2).ne.0).and.(iand(iout,4).ne.0).and. &
           ((ifirst.eq.1).and.(itype.eq.1))) then
         let = 'm'
         iitime = time(ifirst)
         call getfnmd(let,ishot,iitime,eqdsk)
         if (ktime.le.9) encode (4,1030,last) ktime
         if ((ktime.gt.9).and.(ktime.le.99))  &
              encode (4,1040,last) ktime
         if ((ktime.gt.99).and.(ktime.le.999)) &
              encode (4,1050,last) ktime
         if (ktime.gt.999) encode (4,1060,last) ktime
 1030    format ('000',i1)
 1040    format ('00',i2)
 1050    format ('0',i3)
 1060    format (i4)
         eqdsk = eqdsk(1:13)//'_'//last
         nceq = NCCRE(eqdsk,NCCLOB,ierr)              
!
      elseif ((iand(iout,2).ne.0).and.(iand(iout,4).eq.0).and. &
              (ifirst.eq.1).and.(itype.eq.1)) then
         if (ishot.le.99999) encode (10,1010,eqdsk) ishot
         if (ishot.gt.99999) encode (10,1020,eqdsk) ishot
 1010    format ('m0',i5,'.nc')
 1020    format ('m',i6,'.nc')
         nceq = NCCRE(eqdsk,NCCLOB,ierr) ! create file, overwrite if exists
!
! --- creates one file for each slice
!
      elseif ((iand(iout,4).ne.0).and.(iand(iout,2).eq.0).and. &
              (itype.eq.1)) then
         let = 'm'
         iitime = time(ifirst)
         call getfnmu(itimeu,let,ishot,iitime,eqdsk)
         nceq = NCCRE(eqdsk,NCCLOB,ierr)  
!
! --- time-dependent but NOT the first slice NOR the first time,
! --- or single slice but NOT the first time, 
! --- skip define mode and go to write directly
!
      else
         go to 100
      endif
      
      if (ierr.ne.0) return
!
! --- nc file has been created successfully, now define variables, etc.
!
      title = 'EFIT measurement file Mssssss.ttttt/Mssssss.nc'
      call NCAPTC (nceq,NC_GLOBAL,'title',NCCHAR,46,title,ierr)
! 
! --- define unlimited time dimension and scalar and array dimensions
!
      idim_time = NCDDEF (nceq,'dim_time',NCUNLIM,ierr)
      idim_1 = NCDDEF (nceq,'dim_scalar',1,ierr)
      idim_nstark = NCDDEF (nceq,'dim_nstark',nstark,ierr)
      idim_nsilop = NCDDEF (nceq,'dim_nsilop',nsilop,ierr)
      idim_magpri = NCDDEF (nceq,'dim_magpri',magpri,ierr)
      idim_nfcoil = NCDDEF (nceq,'dim_nfcoil',nfcoil,ierr)
      idim_nesum  = NCDDEF (nceq,'dim_nesum' ,nesum ,ierr)
      npress1 = npress
      if (npress.eq.0) npress1 = 1
      idim_npress = NCDDEF (nceq,'dim_npress',npress1,ierr)
      idim_npcurn = NCDDEF (nceq,'dim_npcurn',npcurn,ierr)
      idim_nitera = NCDDEF (nceq,'dim_nitera',nitera,ierr)
      dim2(2) = idim_time
!
! --- define variables
!
      id_shot = NCVDEF (nceq,'shot',NCLONG,1,idim_1,ierr)
      call NCAPTC (nceq,id_shot,'long_name',NCCHAR,11, &
                  'shot number',ierr)
!
      id_time = NCVDEF (nceq,'time',NCFLOAT,1,idim_time,ierr)
      call NCAPTC (nceq,id_time,'units',NCCHAR,4,'msec',ierr)
!
! --- MSE measurements
!
      id_fungam = NCVDEF (nceq,'msefitfun',NCLONG,1,idim_1,ierr)
      call NCAPTC (nceq,id_fungam,'long_name',NCCHAR,16, &
                  'MSE fit function',ierr)

      dim2(1) = idim_nstark
      id_tangam = NCVDEF (nceq,'tangam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_tangam,'long_name',NCCHAR,35, &
      'tangent of measured MSE pitch angle',ierr)
!
      id_tangam_uncor = NCVDEF (nceq,'tangam_uncor',NCFLOAT,&
                             2,dim2,ierr)
      call NCAPTC (nceq,id_tangam_uncor,'long_name',NCCHAR,54, &
      'tangent of measured MSE pitch angle w/o cer correction',ierr)
!
      id_fgam = NCVDEF (nceq,'fixgam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_fgam,'long_name',NCCHAR,58, &
      'radians correction of tangam for spatial averaging effects',ierr)
!
      id_cmgam = NCVDEF (nceq,'cmgam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_cmgam,'long_name',NCCHAR,30, &
      'calculated polarimetry signals',ierr)
!
      id_siggam = NCVDEF (nceq,'siggam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_siggam,'long_name',NCCHAR,21, &
      'uncertainty of tangam',ierr)
!
      id_rrgam = NCVDEF (nceq,'rrgam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_rrgam,'long_name',NCCHAR,22, &
      'radius of MSE channels',ierr)
!
      id_zzgam = NCVDEF (nceq,'zzgam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_zzgam,'long_name',NCCHAR,26, &
      'Z position of MSE channels',ierr)
!
      id_a1gam = NCVDEF (nceq,'a1gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_a1gam,'long_name',NCCHAR,45, &
      'viewing geometry coefficients of MSE channels',ierr)
!
      id_a2gam = NCVDEF (nceq,'a2gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_a2gam,'long_name',NCCHAR,45, &
      'viewing geometry coefficients of MSE channels',ierr)
!
      id_a3gam = NCVDEF (nceq,'a3gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_a3gam,'long_name',NCCHAR,45, &
      'viewing geometry coefficients of MSE channels',ierr)
!
      id_a4gam = NCVDEF (nceq,'a4gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_a4gam,'long_name',NCCHAR,45, &
      'viewing geometry coefficients of MSE channels',ierr)
!
      id_a5gam = NCVDEF (nceq,'a5gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_a5gam,'long_name',NCCHAR,45, &
      'viewing geometry coefficients of MSE channels',ierr)
!
      id_a6gam = NCVDEF (nceq,'a6gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_a6gam,'long_name',NCCHAR,45, &
      'viewing geometry coefficients of MSE channels',ierr)
!
      id_a7gam = NCVDEF (nceq,'a7gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_a7gam,'long_name',NCCHAR,45, &
      'viewing geometry coefficients of MSE channels',ierr)
!
      id_a8gam = NCVDEF (nceq,'a8gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_a8gam,'long_name',NCCHAR,45, &
      'viewing geometry coefficients of MSE channels',ierr)
!
      id_fwtgam = NCVDEF (nceq,'fwtgam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_fwtgam,'long_name',NCCHAR,31, &
      'fitting weight for MSE channels',ierr)
!
      id_msebkp = NCVDEF (nceq,'msebkp',NCFLOAT,1,idim_time,ierr)
      call NCAPTC (nceq,id_msebkp,'long_name',NCCHAR,30, &
      'background substraction switch',ierr)
!     
      id_mseport = NCVDEF (nceq,'vport',NCLONG,1,idim_nstark,ierr)
      call NCAPTC (nceq,id_mseport,'long_name',NCCHAR,33, &
      'view port locations of MSE system',ierr)
!
      if( msefitfun .eq. 1)then

          id_gaingam = NCVDEF (nceq,'mcal_gain',NCFLOAT,2,dim2,ierr)
          call NCAPTC (nceq,id_gaingam,'long_name',NCCHAR,38, &
      'gain param for tangent offset function',ierr)

          id_slopegam = NCVDEF (nceq,'mcal_slope',NCFLOAT, &
                              2,dim2,ierr)
          call NCAPTC (nceq,id_slopegam,'long_name',NCCHAR,39, &
      'slope param for tangent offset function',ierr)

          id_scalegam = NCVDEF (nceq,'mcal_scale',NCFLOAT, &
                              2,dim2,ierr)
          call NCAPTC (nceq,id_scalegam,'long_name',NCCHAR,39, &
      'scale param for tangent offset function',ierr)

          id_offsetgam = NCVDEF (nceq,'mcal_offset',NCFLOAT, &
                              2,dim2,ierr)
          call NCAPTC (nceq,id_offsetgam,'long_name',NCCHAR,40, &
      'offset param for tangent offset function',ierr)


      else if( msefitfun .eq. 3)then

          id_gaingam = NCVDEF (nceq,'mcal3_gain',NCFLOAT,2,dim2,ierr)
          call NCAPTC (nceq,id_gaingam,'long_name',NCCHAR,37, &
      'gain param for tangent slope function',ierr)

          id_slopegam = NCVDEF (nceq,'mcal3_phase',NCFLOAT,2,dim2,ierr)
          call NCAPTC (nceq,id_slopegam,'long_name',NCCHAR,38, &
      'phase param for tangent slope function',ierr)

          id_scalegam = NCVDEF (nceq,'mcal3_btscale',NCFLOAT, &
                              2,dim2,ierr)
          call NCAPTC (nceq,id_scalegam,'long_name',NCCHAR,40, &
      'btscale param for tangent slope function',ierr)

          id_offsetgam = NCVDEF (nceq,'mcal3_dc_offset',NCFLOAT, &
                              2,dim2,ierr)
          call NCAPTC (nceq,id_offsetgam,'long_name',NCCHAR,42, &
      'dc_offset param for tangent slope function',ierr)

      else if( msefitfun .eq. 4)then

          id_gaingam = NCVDEF (nceq,'mcal4_gain',NCFLOAT,2,dim2,ierr)
          call NCAPTC (nceq,id_gaingam,'long_name',NCCHAR,37, &
      'gain param for tangent slope function',ierr)

          id_slopegam = NCVDEF (nceq,'mcal4_phase',NCFLOAT,2,dim2,ierr)
          call NCAPTC (nceq,id_slopegam,'long_name',NCCHAR,38, &
      'phase param for tangent slope function',ierr)

          id_scalegam = NCVDEF (nceq,'mcal4_btscale',NCFLOAT, &
                              2,dim2,ierr)
          call NCAPTC (nceq,id_scalegam,'long_name',NCCHAR,40, &
      'btscale param for tangent slope function',ierr)

          id_offsetgam = NCVDEF (nceq,'mcal4_dc_offset',NCFLOAT, &
                              2,dim2,ierr)
          call NCAPTC (nceq,id_offsetgam,'long_name',NCCHAR,42, &
      'dc_offset param for tangent slope function',ierr)

      endif
!
! --- Magnetic measurements
!
      dim2(1) = idim_nsilop
      id_silopt = NCVDEF (nceq,'silopt',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_silopt,'long_name',NCCHAR,19, &
      'measured flux loops',ierr)
!
      id_csilop = NCVDEF (nceq,'csilop',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_csilop,'long_name',NCCHAR,21, &
      'calculated flux loops',ierr)
!
      id_fwtsi = NCVDEF (nceq,'fwtsi',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_fwtsi,'long_name',NCCHAR,21, &
      'weight for flux loops',ierr)
!
      dim2(1) = idim_magpri
      id_expmpi = NCVDEF (nceq,'expmpi',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_expmpi,'long_name',NCCHAR,24, &
      'measured magnetic probes',ierr)
!
      id_fwtmp2 = NCVDEF (nceq,'fwtmp2',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_fwtmp2,'long_name',NCCHAR,26, &
      'weight for magnetic probes',ierr)
!
      id_cmpr2 = NCVDEF (nceq,'cmpr2',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_cmpr2,'long_name',NCCHAR,26, &
      'calculated magnetic probes',ierr)
!
      dim2(1) = idim_nfcoil
      id_fccurt = NCVDEF (nceq,'fccurt',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_fccurt,'long_name',NCCHAR,30, &
      'measured F-coil currents (Amp)',ierr)
!
      id_ccbrsp = NCVDEF (nceq,'ccbrsp',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_ccbrsp,'long_name',NCCHAR,32, &
      'calculated F-coil currents (Amp)',ierr)
!
      id_fwtfc = NCVDEF (nceq,'fwtfc',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_fwtfc,'long_name',NCCHAR,26, &
      'weight for F-coil currents',ierr)
!
      dim2(1) = idim_nesum
      id_eccurt = NCVDEF (nceq,'eccurt',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_eccurt,'long_name',NCCHAR,30, &
      'measured E-coil currents (Amp)',ierr)
!
      id_cecurr = NCVDEF (nceq,'cecurr',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_cecurr,'long_name',NCCHAR,32, &
      'calculated E-coil currents (Amp)',ierr)
!
      id_fwtec = NCVDEF (nceq,'fwtec',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_fwtec,'long_name',NCCHAR,26, &
      'weight for E-coil currents',ierr)
!
      id_diamag = NCVDEF (nceq,'diamag',NCFLOAT,1,idim_time,ierr)
      call NCAPTC (nceq,id_diamag,'long_name',NCCHAR,25, &
      'measured diamagnetic flux',ierr)
!
      id_cdflux = NCVDEF (nceq,'cdflux',NCFLOAT,1,idim_time,ierr)
      call NCAPTC (nceq,id_cdflux,'long_name',NCCHAR,27, &
      'calculated diamagnetic flux',ierr)
!
      id_sigdia = NCVDEF (nceq,'sigdia',NCFLOAT,1,idim_time,ierr)
      call NCAPTC (nceq,id_sigdia,'long_name',NCCHAR,31, &
      'uncertainty of diamagnetic flux',ierr)
!
! --- plasma current
!
      id_plasma = NCVDEF (nceq,'plasma',NCFLOAT,1,idim_time,ierr)
      call NCAPTC (nceq,id_plasma,'long_name',NCCHAR,29, &
      'measured plasma current (Amp)',ierr)
!
      id_cpasma = NCVDEF (nceq,'cpasma',NCFLOAT,1,idim_time,ierr)
      call NCAPTC (nceq,id_cpasma,'long_name',NCCHAR,31, &
      'calculated plasma current (Amp)',ierr)
!
      id_darea = NCVDEF (nceq,'darea',NCFLOAT,1,idim_time,ierr)
      call NCAPTC (nceq,id_darea,'long_name',NCCHAR,33, &
      'plasma coefficients normalization',ierr)
!
      dim2(1) = idim_npcurn
      id_xrsp = NCVDEF (nceq,'xrsp',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_xrsp,'long_name',NCCHAR,19, &
      'plasma coefficients',ierr)
!
! --- pressure
!
      dim2(1) = idim_npress
      id_pressr = NCVDEF (nceq,'pressr',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_pressr,'long_name',NCCHAR,57, &
      'measured pressure vs. normalized flux (kinetic fits only)',ierr)      
!
      id_cpress = NCVDEF (nceq,'cpress',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_cpress,'long_name',NCCHAR,59, &
      'calculated pressure vs. normalized flux (kinetic fits only)', &
      ierr)      
!
      id_rpress = NCVDEF (nceq,'rpress',NCFLOAT,2,dim2,ierr)
!vas f90 modifi.
vastext(1:59) = '<0 - input pressure profile vs. flux; >0 - R coordinates of'
vastext(60:85)= 'input pressure profile (m)'

      call NCAPTC (nceq,id_rpress,'long_name',NCCHAR,86, &
!vas      '<0 - input pressure profile vs. flux; >0 - R coordinates of input &
!vas       pressure profile (m)',ierr)
      vastext,ierr)
!
      id_zpress = NCVDEF (nceq,'zpress',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_zpress,'long_name',NCCHAR,43, &
      'Z coordinates of input pressure profile (m)',ierr)
!
      id_sigpre = NCVDEF (nceq,'sigpre',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_sigpre,'long_name',NCCHAR,24, &
      'uncertainty for pressure',ierr)
!
      id_saipre = NCVDEF (nceq,'saipre',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_saipre,'long_name',NCCHAR,17, &
      'chisq of pressure',ierr)

!
! --- fitting parameters
!
      dim2(1) = idim_nstark
      id_chigam = NCVDEF (nceq,'chigam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_chigam,'long_name',NCCHAR,23, &
      'chisq vs. polarimetries',ierr)
!
      dim2(1) = idim_magpri
      id_saimpi = NCVDEF (nceq,'saimpi',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_saimpi,'long_name',NCCHAR,25, &
      'chisq vs. magnetic probes',ierr)
!
      dim2(1) = idim_nsilop
      id_saisil = NCVDEF (nceq,'saisil',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_saisil,'long_name',NCCHAR,19, &
      'chisq vs. PSI loops',ierr)
!
      dim2(1) = idim_nitera
      id_czmaxi = NCVDEF (nceq,'czmaxi',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_czmaxi,'long_name',NCCHAR,21, &
      'Zm (cm) vs. iteration',ierr)
!
      id_cchisq = NCVDEF (nceq,'cchisq',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_cchisq,'long_name',NCCHAR,19, &
      'chisq vs. iteration',ierr)
!
      id_cerror = NCVDEF (nceq,'cerror',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_cerror,'long_name',NCCHAR,19, &
      'error vs. iteration',ierr)
!
      call NCENDF (nceq,ierr)             ! leave define mode
!
! --- write values of variables
!
      call NCVPT1 (nceq,id_shot,1,ishot,ierr)
 100  continue
      m = 1
      if (iand(iout,2).ne.0) &
      m = ifirst
      n = ilast-ifirst+1
      c11(1) = 1
      c11(2) = m
      cnn(2) = n
!
      if (itype.eq.1) then
!
! --- itype = 1, called by main routine within time loop,
! --- writes variables that are time-dependent but do not have time dimension
!
      cnn(1) = nstark
      do j=1, nstark
         zwork(j) = real(fwtgam(j))
      enddo
      call NCVPT (nceq,id_fwtgam,c11,cnn,zwork,ierr)

      call NCVPT (nceq,id_msebkp,m,n,real(msebkp),ierr)
!
      cnn(1) = nsilop
      do j=1, nsilop
         zwork(j) = real(fwtsi(j))
      enddo
      call NCVPT (nceq,id_fwtsi,c11,cnn,zwork,ierr)
      cnn(1) = magpri
      do j=1, magpri
         zwork(j) = real(fwtmp2(j))
      enddo
      call NCVPT (nceq,id_fwtmp2,c11,cnn,zwork,ierr)
      cnn(1) = nfcoil
      do j=1, nfcoil
         zwork(j) = real(fwtfc(j))
      enddo
      call NCVPT (nceq,id_fwtfc,c11,cnn,zwork,ierr)
      cnn(1) = nesum
      do j=1, nesum
         zwork(j) = real(fwtec(j))
      enddo
      call NCVPT (nceq,id_fwtec,c11,cnn,zwork,ierr)
!
!     Note that pressr has been changed to calculated pressure from the 
!     measured one. It is switched back in measurement file.
!
      cnn(1) = npress1
      do j=1, npress1
         zwork(j) = real(premea(j))
      enddo
      call NCVPT (nceq,id_pressr,c11,cnn,zwork,ierr)
      do j=1, npress1
         zwork(j) = real(pressr(j))
      enddo
      call NCVPT (nceq,id_cpress,c11,cnn,zwork,ierr)
      do j=1, npress1
         zwork(j) = real(rpress(j))
      enddo
      call NCVPT (nceq,id_rpress,c11,cnn,zwork,ierr)
      do j=1, npress1
         zwork(j) = real(zpress(j))
      enddo
      call NCVPT (nceq,id_zpress,c11,cnn,zwork,ierr)
      do j=1, npress1
         zwork(j) = real(sigpre(j))
      enddo
      call NCVPT (nceq,id_sigpre,c11,cnn,zwork,ierr)
      do j=1, npress1
         zwork(j) = real(saipre2(j))
      enddo
      call NCVPT (nceq,id_saipre,c11,cnn,zwork,ierr)
!
      cnn(1) = nesum
      do j=1, nesum
         zwork(j) = real(cecurr(j))
      enddo
      call NCVPT (nceq,id_cecurr,c11,cnn,zwork,ierr)
      cnn(1) = kwcurn
      do i=1,kwcurn
         zwork(i)=brsp(nfcoil+i)/darea
      enddo
      call NCVPT (nceq,id_darea,m,n,real(darea),ierr)
      call NCVPT (nceq,id_xrsp,c11,cnn,zwork,ierr)
!
      cnn(1) = nstark
      do i=1, nstark
         zwork(i) = real(chigam(i))
      enddo
      call NCVPT (nceq,id_chigam,c11,cnn,zwork,ierr)
      cnn(1) = magpri
      do i=1, magpri
         zwork(i) = real(saimpi(i))
      enddo
      call NCVPT (nceq,id_saimpi,c11,cnn,zwork,ierr)
      cnn(1) = nsilop
      do i=1, nsilop
         zwork(i) = real(saisil(i))
      enddo
      call NCVPT (nceq,id_saisil,c11,cnn,zwork,ierr)
!
      cnn(1) = nitera
      do i=1, nitera
         zwork(i) = real(czmaxi(i))
      enddo
      call NCVPT (nceq,id_czmaxi,c11,cnn,zwork,ierr)
      do i=1, nitera
         zwork(i) = real(cchisq(i))
      enddo
      call NCVPT (nceq,id_cchisq,c11,cnn,zwork,ierr)
      do i=1, nitera
         zwork(i) = real(cerror(i))
      enddo
      call NCVPT (nceq,id_cerror,c11,cnn,zwork,ierr)
!
! --- following variables do NOT have time dimension.
!
      call NCVPT (nceq,id_mseport,1,nstark,mseport,ierr)
      endif
!
      if (((itype.eq.1).and.(iand(iout,2).eq.0)).or. &
          ((itype.eq.2).and.(iand(iout,2).ne.0))) then
!
! --- following variables have time dimension.
! --- if called by main routine within time loop and in individual mode,
! --- writes a slice of the variables according the time specified;
! --- if called by main routine out of time loop (itype = 2) and in 
! --- accumulative mode, writes variables in block.
!
      do i=1, ntime
         zwork(i) = real(time(i))
      enddo
      call NCVPT (nceq,id_time,m,n,zwork(ifirst),ierr)
      cnn(1) = nstark
      do i=1, nstark
      do j=1, ntime
         zcmgam(i,j) = real(cmgam(i,j))
      enddo
      enddo
      call NCVPT (nceq,id_cmgam,c11,cnn,zcmgam(1,ifirst),ierr)
      cnn(1) = nsilop
      do i=1, nsilop
      do j=1, ntime
         zcsilop(i,j) = real(csilop(i,j))
      enddo
      enddo
      call NCVPT (nceq,id_csilop,c11,cnn,zcsilop(1,ifirst),ierr)
      cnn(1) = magpri
      do i=1, magpri
      do j=1, ntime
         zcmpr2(i,j) = real(cmpr2(i,j))
      enddo
      enddo
      call NCVPT (nceq,id_cmpr2,c11,cnn,zcmpr2(1,ifirst),ierr)
      cnn(1) = nfcoil
      do i=1, nfcoil
      do j=1, ntime
         zccbrsp(i,j) = real(ccbrsp(i,j))
      enddo
      enddo
      call NCVPT (nceq,id_ccbrsp,c11,cnn,zccbrsp(1,ifirst),ierr)
      cnn(1) = nesum
      do i=1, ntime
         zwork(i) = real(diamag(i))
      enddo
      call NCVPT (nceq,id_diamag,m,n,zwork(ifirst),ierr)
      do i=1, ntime
         zwork(i) = real(sigdia(i))
      enddo
      call NCVPT (nceq,id_sigdia,m,n,zwork(ifirst),ierr)
      do i=1, ntime
         zwork(i) = real(cdflux(i))
      enddo
      call NCVPT (nceq,id_cdflux,m,n,zwork(ifirst),ierr)
!
      do i=1, ntime
         zwork(i) = real(pasmat(i))
      enddo
      call NCVPT (nceq,id_plasma,m,n,zwork(ifirst),ierr)
      do i=1, ntime
         zwork(i) = real(cpasma(i))
      enddo
      call NCVPT (nceq,id_cpasma,m,n,zwork(ifirst),ierr)
!
! --- following data and their corresponding ncdf variables have 
! --- reversed dimensions
!
      call NCVPT1 (nceq,id_fungam,1,msefitfun,ierr)

      stride(1) = 1
      stride(2) = 1
      cnn(1) = nstark
      imap(2) = 4                         ! number of bytes in float
!     imap(2) = imap(2)*2                 ! number of bytes in double
      imap(1) = imap(2)*ntime
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(tangam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_tangam,c11,cnn,stride,imap, &
           zstark(ifirst,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(tangam_uncor(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_tangam_uncor,c11,cnn,stride,imap, &
           zstark(ifirst,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(spatial_fix(j,i))
      enddo
      enddo
      call NCVPTG (nceq,id_fgam,c11,cnn,stride,imap, &
           zstark(ifirst,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(siggam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_siggam,c11,cnn,stride,imap, &
           zstark(ifirst,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(rrgam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_rrgam,c11,cnn,stride,imap, &
           zstark(ifirst,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(zzgam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_zzgam,c11,cnn,stride,imap, &
           zstark(ifirst,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(a1gam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_a1gam,c11,cnn,stride,imap, &
           zstark(ifirst,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(a2gam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_a2gam,c11,cnn,stride,imap, &
           zstark(ifirst,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(a3gam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_a3gam,c11,cnn,stride,imap, &
           zstark(ifirst,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(a4gam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_a4gam,c11,cnn,stride,imap, &
           zstark(ifirst,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(a5gam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_a5gam,c11,cnn,stride,imap, &
           zstark(ifirst,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(a6gam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_a6gam,c11,cnn,stride,imap, &
           zstark(ifirst,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(a7gam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_a7gam,c11,cnn,stride,imap, &
           zstark(ifirst,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(a8gam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_a8gam,c11,cnn,stride,imap, &
           zstark(ifirst,1),ierr)

      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(rmse_gain(j))
      enddo
      enddo
      call NCVPTG (nceq,id_gaingam,c11,cnn,stride,imap, &
           zstark(ifirst,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(rmse_slope(j))
      enddo
      enddo
      call NCVPTG (nceq,id_slopegam,c11,cnn,stride,imap, &
           zstark(ifirst,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(rmse_scale(j))
      enddo
      enddo
      call NCVPTG (nceq,id_scalegam,c11,cnn,stride,imap, &
           zstark(ifirst,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(rmse_offset(j))
      enddo
      enddo
      call NCVPTG (nceq,id_offsetgam,c11,cnn,stride,imap, &
           zstark(ifirst,1),ierr)
!
      cnn(1) = nsilop
      do i=1, ntime
      do j=1, nsilop
         zsilopt(i,j) = real(silopt(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_silopt,c11,cnn,stride,imap, &
           zsilopt(ifirst,1),ierr)
      cnn(1) = magpri
      do i=1, ntime
      do j=1, magpri
         zexpmpi(i,j) = real(expmpi(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_expmpi,c11,cnn,stride,imap, &
           zexpmpi(ifirst,1),ierr)
      cnn(1) = nfcoil
      do i=1, ntime
      do j=1, nfcoil
         zfccurt(i,j) = real(fccurt(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_fccurt,c11,cnn,stride,imap, &
           zfccurt(ifirst,1),ierr)
      cnn(1) = nesum
      do i=1, ntime
      do j=1, nesum
         zeccurt(i,j) = real(eccurt(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_eccurt,c11,cnn,stride,imap, &
           zeccurt(ifirst,1),ierr)
!
      call NCCLOS (nceq,ierr)             ! close the file
      endif
!
      return
      end
