!**********************************************************************
!>
!!    writes out the GAQ type eqdsk.
!!    
!!
!!    @param jtime : time index
!!
!**********************************************************************
      subroutine write_g(jtime)
      use set_kinds, only: i4,r4
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      real*8 seval
      integer*4, intent(in) :: jtime
      integer*4 ijtime,i,j,jb,kk,ierold,ioerr,kkstark,nbsave,nbabs,ndel, &
                nqpsi,nsw,nsh
      real*8 btor,enps,plasma,saisq,ssibry,ssimag,xdiff,xdim,xdum,xn, &
             zdiff,zdim,zmid
      real*8 psirz(nw,nh),pcurrz(nw,nh),workk(nw),dmion(nw), &
             bworm(nmass),cworm(nmass),dworm(nmass),coils(nsilop), &
             expmp2(magpri),prexp(nrogow),tgamma(nstark),sgamma(nstark), &
             rrrgam(nstark),zzzgam(nstark),aa1gam(nstark),aa2gam(nstark), &
             aa3gam(nstark),aa4gam(nstark),aa5gam(nstark),aa6gam(nstark), &
             aa7gam(nstark),aa8gam(nstark),tgammauncor(nstark), &
             bmsels(nmsels),sbmsels(nmsels),fwtbmsels(nmsels), &
             rrmsels(nmsels),zzmsels(nmsels),l1msels(nmsels), &
             l2msels(nmsels),l4msels(nmsels),emsels(nmsels),semsels(nmsels)
      character eqdsk*72,header*42,wform*20,let,fit_type*3
      character*10 vers(6)
      namelist/out1/ishot,itime,betap0,rzero,qenp,enp,emp,plasma, &
           expmp2,coils,btor,rcentr,brsp,icurrt,rbdry,zbdry, &
           nbdry,fwtsi,fwtcur,mxiter,nxiter,limitr,xlim,ylim,error, &
           iconvr,ibunmn,pressr,rpress,nqpsi,npress,sigpre,qpsi, &
           pressr,pressw
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
      namelist/chiout/saisil,saimpi,chipasma,saipre
      namelist/eccd/kkstark,chigamt,chigam,bzmse,psiecn,dpsiecn, &
              saisq,cjeccd
!
      xdum=0.0
      psirz=0.0
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
        fit_type='KIM'
      elseif (kprfit.gt.0) then
        fit_type='KIN'
      elseif (kstark.gt.0) then
        fit_type='MSE'
      else
        fit_type='MAG'
      endif
!
      plasma=ipmhd(jtime)
      btor=bcentr(jtime)
      coils=csilop(:,jtime)
      expmp2=cmpr2(:,jtime)
      prexp=crogow(:,jtime)
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
!----------------------------------------------------------------------
!--   set grid dimensions, flux sign, and grid subsampling step-size
!----------------------------------------------------------------------
      xdim=rgrid(nw)-rgrid(1)
      zdim=zgrid(nh)-zgrid(1)
      zmid=(zgrid(1)+zgrid(nh))/2.0
      do i=1,nw
        do j=1,nh
          kk=(i-1)*nh+j
          if (ivacum.eq.0) then
            if (ipmeas(jtime).gt.0.0) then
              psirz(i,j)=-psi(kk)
            else
              psirz(i,j)=psi(kk)
            endif
          else
            psirz(i,j)=-psi(kk)
          endif
        enddo
      enddo
      nsw=(nw-1)/(nw_sub-1)
      nsh=(nh-1)/(nh_sub-1)
!----------------------------------------------------------------------------
!--   convert the pointwise current density to spatial density and match the
!--   orientation of psi (if requested)
!----------------------------------------------------------------------------
      if (iplcout.eq.2) then
        pcurrz=0.0
        if (ivacum.eq.0) then
          xdiff=xdim/(nw-1)
          zdiff=zdim/(nh-1)
          do i=1,nw
            do j=1,nh
              kk=(i-1)*nh+j
                  pcurrz(i,j)=pcurrt(kk)/xdiff/zdiff
            enddo
          enddo
        endif
      endif
!----------------------------------------------------------------------
!--   setup the file header and name
!----------------------------------------------------------------------
      write(vers(1),1040)
      write(vers(2),1050) efitdate(1:5)
      write(vers(3),1060) efitdate(6:10)
      if (ishot.le.99999) then
        write(vers(4),1070) ishot
      else
        write(vers(4),1073) ishot
      endif
      write(vers(5),1080) ijtime
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
      if (ipmeas(jtime).gt.0.0) then
        ssimag=-simag
        ssibry=-psibry
      else
        ssimag=simag
        ssibry=psibry
      endif
      eqdsk_format: if (keqdsk.eq.1) then
      write(neqdsk,2000) (vers(i),i=1,6),0,abs(nw_sub),abs(nh_sub)
      write(neqdsk,2020) xdim,zdim,rzero,rgrid(1),zmid
      write(neqdsk,2020) rmaxis,zmaxis,ssimag,ssibry,bcentr(jtime)
      write(neqdsk,2020) ipmhd(jtime),ssimag,xdum,rmaxis,xdum
      write(neqdsk,2020) zmaxis,xdum,ssibry,xdum,xdum
      write(neqdsk,2020) (fpol(i),i=1,nw,nsw)
      write(neqdsk,2020) (pres(i),i=1,nw,nsw)
      if (ipmeas(jtime).gt.0.0) then
        workk=-ffprim
      else
        workk=ffprim
      endif
      write(neqdsk,2020) (workk(i),i=1,nw,nsw)
      if (ipmeas(jtime).gt.0.0) then
        workk=-pprime
      else
        workk=pprime
      endif
      write(neqdsk,2020) (workk(i),i=1,nw,nsw)
      write(neqdsk,2020) ((psirz(i,j),i=1,nw,nsw),j=1,nh,nsh)
      write(neqdsk,2020) (qpsi(i),i=1,nw,nsw)
      write(neqdsk,2022) nbbbs,limitr
      write(neqdsk,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)
      write(neqdsk,2020) (xlim(i),ylim(i),i=1,limitr)
!----------------------------------------------------------------------
!--   write out rotation information                                 --
!----------------------------------------------------------------------
      write(neqdsk,2024) kvtor,rvtor,nmass
      if (kvtor.gt.0) then
        write(neqdsk,2020) (pressw(i),i=1,nw,nsw)
        if (ipmeas(jtime).gt.0.0) then
          workk=-pwprim
        else
          workk=pwprim
        endif
        write(neqdsk,2020) (workk(i),i=1,nw,nsw)
      endif
!----------------------------------------------------------------------
!--   write out ion mass density profile if available                --
!----------------------------------------------------------------------
      if (nmass.gt.0) then
        write(neqdsk,2020) (dmion(i),i=1,nw,nsw)
      endif
      write(neqdsk,2020) (rhovn(i),i=1,nw,nsw)
      write(neqdsk,2026) keecur
      if (keecur.gt.0) then
        if (ipmeas(jtime).gt.0.0) then
          workk=-epoten
        else
          workk=epoten
        endif
        write(neqdsk,2020) (workk(i),i=1,nw,nsw)
      endif
! note: unlike the rest of the file, these optional extras have
!       never been described completely with variables available here
!       (since being added initially in 2004)
      if (iplcout.gt.0) then
        if (iplcout.eq.1) then
          ! grid sub-sampling not setup here
          if (ishot.le.99999) then
            write(neqdsk,3000) nw,nh,ishot,itime
          else
            write(neqdsk,3003) nw,nh,ishot,itime
          endif
          write(neqdsk,2020) rgrid(1),rgrid(nw),zgrid(1),zgrid(nh)
          write(neqdsk,2020) (brsp(i),i=1,nfsum)  ! also in m and a-files
          write(neqdsk,2020) (ecurrt(i),i=1,nesum) ! also in m and a-files
          write(neqdsk,2020) (pcurrt(i),i=1,nwnh)
        elseif (iplcout.eq.2) then
          write(neqdsk,2020) ((pcurrz(i,j),i=1,nw,nsw),j=1,nh,nsh)
        endif
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
              write(neqdsk,9991) tmpdata
              read (nsnapf,9991,iostat=ioerr) tmpdata
              if (ioerr.ne.0) exit
              if (INDEX(tmpdata,'/')/=0) then
                write(neqdsk,9991) tmpdata
                exit
              endif
            enddo
            close (unit=nsnapf)
          endif
 9991     format (a)
        endif
      endif
!
      ! grid subsampling not setup here
      nqpsi=nw
      limitr=limitr-1
      write(neqdsk,out1)

      call ppstore
      call ffstore
      call wwstore
      call eestore
      write(neqdsk,basis)
      limitr=limitr+1
      MSE: if (kdomse.gt.0.or.kstark.gt.0) then
        if (kdomse.gt.0) then
          tgamma=cmgam(:,jtime)
          tgammauncor=tgamma
          sgamma=tgamma*0.05_dp
        else
          tgamma=tangam(jtime,:)
          tgammauncor=tangam_uncor(jtime,:)
          sgamma=siggam(jtime,:)
        endif
        rrrgam=rrgam(jtime,:)
        zzzgam=zzgam(jtime,:)
        aa1gam=a1gam(jtime,:)
        aa2gam=a2gam(jtime,:)
        aa3gam=a3gam(jtime,:)
        aa4gam=a4gam(jtime,:)
        aa5gam=a5gam(jtime,:)
        aa6gam=a6gam(jtime,:)
        aa7gam=a7gam(jtime,:)
        aa8gam=a8gam(jtime,:)
        write(neqdsk,mseout)
        kkstark=nstark
        saisq=chisq(jtime)
        write(neqdsk,eccd)
      endif MSE
!-----------------------------------------------------------------------
!--   Write out MSE-LS namelist                                       --
!-----------------------------------------------------------------------
      if (mmbmsels.gt.0.or.kdomsels.gt.0) then         
        if (kdomsels.gt.0) then
          bmsels=cmmls(jtime,:)
          sbmsels=0.05_dp*bmsels
        else
          bmsels=bmselt(jtime,:)
          sbmsels=sbmselt(jtime,:)
        endif
        rrmsels=rrmselt(jtime,:)
        zzmsels=zzmselt(jtime,:)
        l1msels=l1mselt(jtime,:)
        l2msels=l2mselt(jtime,:)
        l4msels=l4mselt(jtime,:)
        emsels=emselt(jtime,:)
        semsels=semselt(jtime,:)
        write(neqdsk,in_msels)
      endif
!
      if (kdovt.gt.0) then
        ! grid subsampling not setup here
        npresw=0
        ndel=nw/26
        do i=1,nw,ndel
          npresw=npresw+1
          presw(npresw)=pressw(i)
          sigprw(npresw)=0.1_dp*presw(npresw)
          rpresw(npresw)=-real(i-1,dp)/(nw-1)
        enddo
        write(neqdsk,vtout)
      endif
      write(neqdsk,chiout)
      header = ' '
      write(neqdsk,1042) header,fit_type
!-----------------------------------------------------------------------
!--   binary format                                                   --
!-----------------------------------------------------------------------
      else eqdsk_format
      write(neqdsk) (vers(i),i=1,6),0,abs(nw_sub),abs(nh_sub)
      write(neqdsk) real(xdim,r4),real(zdim,r4),real(rzero,r4), &
                     real(rgrid(1),r4),real(zmid,r4)
      write(neqdsk) real(rmaxis,r4),real(zmaxis,r4),real(ssimag,r4), &
                     real(ssibry,r4),real(bcentr(jtime),r4)
      write(neqdsk) real(ipmhd(jtime),r4),real(ssimag,r4),real(xdum,r4), &
                     real(rmaxis,r4),real(xdum,r4)
      write(neqdsk) real(zmaxis,r4),real(xdum,r4),real(ssibry,r4), &
                     real(xdum,r4),real(xdum,r4)
      write(neqdsk) (real(fpol(i),r4),i=1,nw,nsw)
      write(neqdsk) (real(pres(i),r4),i=1,nw,nsw)
      if (ipmeas(jtime).gt.0.0) then
        workk=-ffprim
      else
        workk=ffprim
      endif
      write(neqdsk) (real(workk(i),r4),i=1,nw,nsw)
      if (ipmeas(jtime).gt.0.0) then
        workk=-pprime
      else
        workk=pprime
      endif
      write(neqdsk) (real(workk(i),r4),i=1,nw,nsw)
      write(neqdsk) ((real(psirz(i,j),r4),i=1,nw,nsw),j=1,nh,nsh)
      write(neqdsk) (real(qpsi(i),r4),i=1,nw,nsw)
      write(neqdsk) nbbbs,limitr
      write(neqdsk) (real(rbbbs(i),r4),real(zbbbs(i),r4),i=1,nbbbs)
      write(neqdsk) (real(xlim(i),r4),real(ylim(i),r4),i=1,limitr)
!----------------------------------------------------------------------
!--   write out rotation information                                 --
!----------------------------------------------------------------------
      write(neqdsk) kvtor,real(rvtor,r4),nmass
      if (nmass.gt.0) then
        write(neqdsk) (real(pressw(i),r4),i=1,nw,nsw)
        write(neqdsk) (real(pwprim(i),r4),i=1,nw,nsw)
      endif
!----------------------------------------------------------------------
!--   write out ion mass density profile if available                --
!----------------------------------------------------------------------
      if (nmass.gt.0) then
        write(neqdsk) (real(dmion(i),r4),i=1,nw,nsw)
      endif
      write(neqdsk) (real(rhovn(i),r4),i=1,nw,nsw)
      write(neqdsk) keecur
      if (keecur.gt.0) then
        if (ipmeas(jtime).gt.0.0) then
          workk=-epoten
        else
          workk=epoten
        endif
        write(neqdsk) (real(workk(i),r4),i=1,nw,nsw)
      endif
! note: unlike the rest of the file, these optional extras have
!       never been described completely with variables available here
!       (since being added initially in 2004)
      if (iplcout.gt.0) then
        if (iplcout.eq.1) then
          ! grid subsampling not setup here
          if (ishot.le.99999) then
            write(neqdsk) nw,nh,ishot,itime
          else
            write(neqdsk) nw,nh,ishot,itime
          endif
          write(neqdsk) real(rgrid(1),r4),real(rgrid(nw),r4), &
                         real(zgrid(1),r4),real(zgrid(nh),r4)
          write(neqdsk) (real(brsp(i),r4),i=1,nfsum)  ! also in m and a-files
          write(neqdsk) (real(ecurrt(i),r4),i=1,nesum) ! also in m and a-files
          write(neqdsk) (real(pcurrt(i),r4),i=1,nwnh)
        elseif (iplcout.eq.2) then
          write(neqdsk) ((real(pcurrz(i,j),r4),i=1,nw,nsw),j=1,nh,nsh)
        endif
      endif
!
      header = ' '
      write(neqdsk) header,fit_type
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
        fwtsi(1:nfsum)=1.
        fwtsi((nfsum+1):nsilop)=0.0
        itime=itime+1
        limitr=limitr-1
        let = 'x'
        call setfnmeq(itimeu,let,ishot,itime,eqdsk)
        open(unit=neqdsk,file=eqdsk,status='old',iostat=ioerr)
        if (ioerr.eq.0) close(unit=neqdsk,status='delete')
        open(unit=neqdsk,file=eqdsk,status='new',delim='quote')
        write(neqdsk,in1)
        call ppstore
        call ffstore
        call wwstore
        call eestore
        write(neqdsk,basis)
        write(neqdsk,inwant)
        close(unit=neqdsk)
        enp=enps
        limitr=limitr+1
        itime=itime-1
        ierchk=ierold
      endif fixed_bdry
!
      return
 1020 format ('x',i5,'.',i3)
 1040 format ('  EFIT  ')
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
