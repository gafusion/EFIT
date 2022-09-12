!**********************************************************************
!>
!!    writes out the GAQ type eqdsk.
!!    
!!
!!    @param jtime : time index
!!
!**********************************************************************
      subroutine write_g(jtime)
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
        ibunmn=1
        fwtcur=1.
        mxiter=1
        nxiter=99
        error=1.0e-04_dp ! TODO: why is this set at a fixed value?
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
      end subroutine write_g
