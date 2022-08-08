#include "config.f"
!**********************************************************************
!>
!!    Subroutine for writing K file as the exclusive output (from SNAP)
!!    
!!    
!!    @param ktime : number of time slices requested
!!
!!    @param kerror : Error flag
!!
!**********************************************************************
      subroutine write_k(ktime,kerror)
      use set_kinds
      use opt_input, only: cmdfile_in,shotfile_in
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
#if defined(USEMPI)
      include 'mpif.h'
#endif

      integer*4, intent(inout) :: ktime
      integer*4, intent(out) :: kerror
      integer*4 i,jtime,idtime,ioerr,itimeb,mtime,irshot,istat,istop, &
                mmemsels,mmstark,ktear
      real*8 plasma,btor,dflux,xltype,xltype_180,vloop,siref, &
             pnbeam,timeus,timems,currn1,currc79,currc139,currc199, &
             curriu30,curriu90,curriu150,curril30,curril90,curril150
      real*8 delt,times
      character ishotime*12,news*72, &
                eqdsk*20,prefix1*1,header*42,fit_type*3
      real*8,dimension(:),allocatable :: coils,expmp2, &
                denr,denv, tgamma,sgamma,rrrgam, &
                zzzgam,aa1gam,aa2gam, aa3gam,aa4gam,aa5gam, &
                aa6gam,aa7gam,tgammauncor
      real*8,dimension(:),allocatable :: bmsels,sbmsels,fwtbmsels, &
                rrmsels,zzmsels,l1msels,l2msels, &
                l4msels,emsels,semsels,fwtemsels
      real*8 :: dnmin
      character*86 snap_ext
      character(1000) :: line
      character(256) table_save
      namelist/in1/ishot,itime,itimeu,qvfit,plasma,expmp2,coils,btor, &
           fwtsi,fwtcur,limitr,fwtmp2,kffcur,kppcur,fwtqa,ierchk, &
           fwtbp,serror,nextra,scrape,itrace,itek,xltype,rcentr,bitip, &
           psibit,bitmpi,denr,denv,siref,fwtfc,brsp,bitfc,iecurr,iplim, &
           ecurrt,ifitvs,vloop,dflux,ifcurr,iavem,icprof,currn1,n1coil, &
           pnbeam,error,errmin,mxiter,xltype_180,icutfp,keqdsk,ibtcomp, &
           fcurbd,pcurbd,kbound,alphafp,kskipvs,vsdamp,kframe,zelip, &
           fwtdlc,sigdlc,elomin,kcalpa,kcgama,calpa,cgama,xalpa,xgama, &
           kzeroj,rzeroj,iaveus,relax,fwtec,bitec,fitsiref, &
           kppfnc,kppknt,ppknt,pptens,kfffnc,kffknt,ffknt,fftens, &
           kwwfnc,kwwknt,wwknt,wwtens,nbdry,rbdry,zbdry,vbit, nbdrymx, &
           ppbdry,kppbdry,pp2bdry,kpp2bdry, &
           ffbdry,kffbdry,ff2bdry,kff2bdry, &
           wwbdry,kwwbdry,ww2bdry,kww2bdry, &
           ktear,kersil,iout,ixray,table_dir,input_dir,store_dir, &
           kpphord,kffhord,keehord,psiecn,dpsiecn,fitzts,isolve, &
           iplcout,imagsigma,errmag,saimin,errmagb,fitfcsum,fwtfcsum,efitversion, &
           kwripre,ifindopt,tolbndpsi,siloplim,use_previous
      namelist/inwant/psiwant,vzeroj,nccoil,currc79,currc139,rexpan, &
           znose,sizeroj,fitdelz,relaxdz,errdelz,oldccomp,nicoil, &
           oldcomp,currc199,curriu30,curriu90, &
           curriu150,curril30,curril90,curril150,ifitdelz,scaledz
      namelist/ink/isetfb,ioffr,ioffz,ishiftz,gain,gainp,idplace, &
           symmetrize,backaverage
      namelist/ins/tgamma,sgamma,fwtgam,rrrgam,zzzgam,aa1gam,aa2gam, &
           aa3gam,aa4gam,aa5gam,aa6gam,aa7gam,msebkp, &
           msefitfun,mse_quiet,mse_spave_on,kwaitmse, &
           dtmsefull,mse_strict,t_max_beam_off,ok_30rt,ok_210lt,&
           mse_usecer,mse_certree,mse_use_cer330,mse_use_cer210,&
           tgammauncor,v30lt,v30rt,v210lt,v210rt
      namelist/in_msels/bmsels,sbmsels,fwtbmsels,rrmsels,zzmsels, &
           l1msels,l2msels,l4msels,emsels,semsels,fwtemsels,kdomsels, &
           fmlscut,synmsels,avemsels
      namelist/ina/spatial_avg_gam
      namelist/inece/necein,teecein0,feece0,errorece0,fwtece0,fwtecebz0, &
           ecefit,ecebzfit,kfitece,kinputece,kcallece,nharm, &
           kfixro,rteo,zteo,kfixrece,rtep,rtem,rpbit,rmbit,robit, &
           nfit,kcmin,fwtnow,mtxece
      namelist/iner/keecur,ecurbd,keefnc,eetens,keebdry,kee2bdry, &
           eebdry,ee2bdry,eeknt,keeknt,keehord
      namelist/efitin/ishot,istore,timeb,dtime,mtime,scrape,nextra, &
           iexcal,itrace,xltype,ivesel,fwtsi,fwtmp2,fwtcur,iprobe, &
           itek,limid,qvfit,fwtbp,kffcur,kppcur,fwtqa,mxiter,  &
           serror,ibatch,ifitvs,fwtfc,iecurr,itimeb,idtime,znose, &
           iavem,iaved,iavev,idite,ifcurr,imerci,iacoil,iaveus, &
           cutip,lookfw,error,errmin,xltype_180,icprof,condin, &
           icutfp,keqdsk,kcaldia,fcurbd,pcurbd,ircfact,zelip, &
           kbound,alphafp,kskipvs,vsdamp,kframe,dnmin,vzeroj, &
           fwtdlc,elomin,fwtgam,saicon,fwacoil,itimeu,nccoil, &
           kcalpa,kcgama,calpa,cgama,xalpa,xgama,n1coil,rexpan, &
           psiwant,ibtcomp,icinit,iplim,kwripre,relax,rzeroj,kzeroj, &
           kppfnc,kppknt,ppknt,pptens,kfffnc,kffknt,ffknt,fftens, &
           kwwfnc,kwwknt,wwknt,wwtens,sizeroj,fwtec, &
           ppbdry,kppbdry,pp2bdry,kpp2bdry,nicoil,oldcomp, &
           ffbdry,kffbdry,ff2bdry,kff2bdry,msefitfun, &
           wwbdry,kwwbdry,ww2bdry,kww2bdry,fwtec,fitdelz,fitsiref, &
           nbdry,rbdry,zbdry,sigrbd,sigzbd,nbskip,msebkp, &
           ktear,keecur,ecurbd,keefnc,eetens,eeknt,keeknt, &
           keebdry,kee2bdry,eebdry,ee2bdry,kersil,iout,ixray, &
           use_alternate_pointnames, alternate_pointname_file, &
           do_spline_fit,table_dir,input_dir,store_dir, &
           kedgep,pedge,pe_psin,pe_width,kedgef,f2edge,fe_psin,fe_width, &
           psiecn,dpsiecn,relaxdz,fitzts,isolve,stabdz, &
           iplcout,errdelz,imagsigma,errmag,ksigma,saimin,errmagb, &
           write_Kfile,fitfcsum,fwtfcsum,appendsnap, &
           mse_quiet,mse_spave_on,kwaitmse,dtmsefull, &
           mse_strict,t_max_beam_off,ifitdelz,scaledz, &
           mse_usecer,mse_certree,mse_use_cer330,mse_use_cer210, &
           ok_30rt,ok_210lt,vbit,nbdrymx,fwtbmsels,fwtemsels, &
           idebug,jdebug,synmsels,avemsels,kwritime, &
           v30lt,v30rt,v210lt,v210rt,ifindopt,tolbndpsi, &
           siloplim,use_previous,ierchk
      namelist/efitink/isetfb,ioffr,ioffz,ishiftz,gain,gainp,idplace, &
           symmetrize,backaverage,lring
      data currn1/0.0/,currc79/0.0/,currc139/0.0/,currc199/0.0/, &
                       curriu30/0.0/,curriu90/0.0/,curriu150/0.0/, &
                       curril30/0.0/,curril90/0.0/,curril150/0.0/

      ALLOCATE(coils(nsilop),expmp2(magpri), &
               denr(nco2r),denv(nco2v), &
               tgamma(nmtark),sgamma(nmtark),rrrgam(nmtark), &
               zzzgam(nmtark),aa1gam(nmtark),aa2gam(nmtark), &
               aa3gam(nmtark),aa4gam(nmtark),aa5gam(nmtark), &
               aa6gam(nmtark),aa7gam(nmtark),tgammauncor(nmtark))
      ALLOCATE(bmsels(nmsels),sbmsels(nmsels),fwtbmsels(nmsels), &
               rrmsels(nmsels),zzmsels(nmsels),l1msels(nmsels), &
               l2msels(nmsels),l4msels(nmsels),emsels(nmsels), &
               semsels(nmsels),fwtemsels(nmsels))

      kerror=0
!---------------------------------------------------------------------
!--   generate input files and command file for running EFIT-AI     --
!---------------------------------------------------------------------
      if (cmdfile_in.ne.'0' .and. cmdfile_in.ne.'') then
        open(unit=nffile,status='old',file=cmdfile_in,iostat=ioerr)
        if (ioerr.ne.0) close(unit=nffile,status='delete')
        open(unit=nffile,status='new',file=cmdfile_in)
        write (nffile,4958)
      endif
      if (shotfile_in.ne.'0' .and. shotfile_in.ne.'') then
        open(unit=nout,status='old',file=shotfile_in,iostat=ioerr)
        if (ioerr.ne.0) open(unit=nout,status='old', &
          file='phys_data:[d3phys.diiid.gsl]'//shotfile_in,iostat=ioerr)
        if (ioerr.ne.0) then
          call errctrl_msg('write_K', &
                           'could not open good shot list file')
          stop
        endif
      endif
      if (snapextin.eq.'none') then ! could come from efit.input
        snap_file = 'efit_snap.dat'
      else
        snap_ext = adjustl(snapextin)
        snap_file = 'efit_snap.dat_'//snap_ext
      endif
      open(unit=neqdsk,status='old',file=snap_file,iostat=ioerr)
      if (ioerr.eq.0) then
        snapfile=snap_file
      else
        open(unit=neqdsk,status='old', &
           file=input_dir(1:lindir)//snapfile)
        snapfile=input_dir(1:lindir)//snapfile
      endif
!
      fwtbmsels=0.0
      ktear=0
      mtime=ktime
      table_save=table_dir
      xltype=0.
      xltype_180=0.0
!
      call set_defaults
!
      read (neqdsk,efitin,iostat=istat)
      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist efitin: '//trim(line)
        stop
      endif
      read (neqdsk,efitink,iostat=istat)
      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist efitink: '//trim(line)
        stop
      endif
      if(ioerr.ne.0) close(unit=neqdsk)
      if (link_efit(1:1).ne.'') then
        ! use directories specified in efit.setup if available
        table_dir=trim(link_efit)//'green/'
        input_dir=trim(link_efit)
      elseif (table_dir.ne.table_save .and. link_efit.eq.'') then
        ! error if table_dir changes (global_allocs already set)
        call errctrl_msg('write_K', &
          'changing experiment during run is not supported (table_dir)')
        kerror=1
        return
      endif
      if (link_store(1:1).ne.'')  store_dir=trim(link_store)
!--   warn that idebug, jdebug, and ktear inputs are deprecated
      if (idebug.ne.0) write(*,*) &
      "idebug input variable is deprecated, set cmake variable instead"
      if (jdebug.ne."NONE") write(*,*) &
      "jdebug input variable is deprecated, set cmake variable instead"
      if(ktear.ne.0) write(*,*) &
      "tearing calculations don't exist, ktear is deprecated"
#ifdef DEBUG_LEVEL2
      write (6,*) 'WRITE_K fwtbmsels= ',(fwtbmsels(i),i=1,nmsels)
#endif
!---------------------------------------------------------------------
!--   specific choice of current profile                           --
!--       ICPROF=1  no edge current density allowed                 --
!--       ICPROF=2  free edge current density                       --
!--       ICPROF=3  weak edge current density constraint            --
!---------------------------------------------------------------------
      select case (icprof)
      case (1)
        kffcur=2
        kppcur=2
        fcurbd=1.
        pcurbd=1.
        fwtbp=1.
        fwtqa=0.
        qvfit=0.
      case (2)
        kffcur=2
        kppcur=2
        fcurbd=0.
        pcurbd=0.
        fwtbp=0.
        fwtqa=0.
        qvfit=0.
      case (3)
        kffcur=3
        kppcur=2
        fcurbd=0.
        pcurbd=0.
        fwtbp=0.
        fwtqa=0.
        qvfit=0.
        kcalpa=1
        calpa(1,1)=0.1_dp
        calpa(2,1)=0.1_dp
        calpa(3,1)=0.1_dp
        xalpa(1)=0.0
        kcgama=1
        cgama(1,1)=0.1_dp
        cgama(2,1)=0.1_dp
        cgama(3,1)=0.1_dp
        xgama(1)=0.0
      end select
      if(mse_usecer .eq. 1) keecur = 0
      if (mse_usecer .eq. 2 .and. keecur .eq. 0) then
        keecur = 2
        keefnc = 0
        itek = 5
      endif
!
      limitr=-limid
      swtfc(1:nfcoil)=fwtfc(1:nfcoil)
      swtec(1:nesum)=fwtec(1:nesum)
      swtmp2(1:magpri)=fwtmp2(1:magpri)
      swtsi(1:nsilop)=fwtsi(1:nsilop)
      swtcur=fwtcur
      swtgam(1:nstark)=fwtgam(1:nstark)
      swtbmsels(1:nmsels)=fwtbmsels(1:nmsels)
      swtemsels(1:nmsels)=fwtemsels(1:nmsels)
 3046 continue
#ifdef DEBUG_LEVEL2
      write (6,*) 'WRITE_K fwtbmsels= ',(fwtbmsels(i),i=1,nmsels)
      write (6,*) 'WRITE_K swtbmsels= ',(swtbmsels(i),i=1,nmsels)
#endif
!
      if (shotfile_in.ne.'0' .and. shotfile_in.ne.'') then
        read (nout,4970,iostat=ioerr) ishotime
        if (ioerr.eq.0) then
          read (ishotime,fmt='(i6,1x,i5)',err=3046) ishot,itime
          times=itime/1000.
          delt=0.002_dp
          ktime=1
          timeb=itime
          dtime=0.
        endif
      else
        times=timeb/1000.
        delt=dtime/1000.
        ioerr=0
      endif
      no_err: if (ioerr.eq.0) then
      if (ifitvs.gt.0) then
        istop=-1
      else
        istop=0
      endif
      if(ishot.gt.108281) n1coil = 0
!----------------------------------------------------------------------
!--   Fetch data                                                     --
!----------------------------------------------------------------------
#if defined(USEMPI)
      if (nproc == 1) then
        call get_constraints(ishot,times,delt,ktime,istop)
      else
        call get_constraints_mpi(ishot,times,delt,ktime,istop)
      endif
#else
      call get_constraints(ishot,times,delt,ktime,istop)
#endif
      if (istop.gt.0) then
        write (6,20000)
        if (shotfile_in.ne.'0' .and. shotfile_in.ne.'') then
          go to 3046
        else
          kerror = 1
          call errctrl_msg('write_K', 'shot data not on disk')
          return
        endif
      endif
!
      mmstark=0
      do i=1,nstark
        if (swtgam(i).gt.1.e-06_dp) mmstark=mmstark+1
      enddo
      if (mmstark.gt.0) then
#if defined(USE_MSE)
#if defined(USEMPI)
        if (nproc == 1) then
          call getstark(ktime)
        else
          call getstark_mpi(ktime)
        endif
#else
        call getstark(ktime)
#endif
#else
        kerror = 1
        call errctrl_msg('write_K','MSE library not available')
        return
#endif
      endif
!
      mmbmsels=0
      mmemsels=0
      do i=1,nmsels
        if (swtbmsels(i).gt.1.e-06_dp) mmbmsels=mmbmsels+1
        if (swtemsels(i).gt.1.e-06_dp) mmemsels=mmemsels+1
      enddo
#ifdef DEBUG_LEVEL2
      write (6,*) 'WRITE_K mmbmsels= ', mmbmsels
#endif
      if (mmbmsels.gt.0) then
        call getmsels(ktime)
      endif
!----------------------------------------------------------------------
!--   get xltype, xltype_180 without reading limiters                --
!----------------------------------------------------------------------
      call getlim(0,xltype,xltype_180,.false.)
!----------------------------------------------------------------------
!--   adjust fitting weights based on FITWEIGHT.DAT                  --
!----------------------------------------------------------------------
      if (lookfw.ge.0) then
        rwtmp2(1:magpri)=0.0
        rwtsi(1:nsilop)=0.0
        open(unit=neqdsk,status='old', &
             file=input_dir(1:lindir)//'fitweight.dat')
 3050   read (neqdsk,*,iostat=ioerr) irshot
        if ((ioerr.eq.0).and.(irshot.le.ishot)) then
          if (irshot.lt.124985) then
            read (neqdsk,*) (rwtsi(i),i=1,nsilol)
          else
            read (neqdsk,*) (rwtsi(i),i=1,nsilop)
          endif
          if (irshot.lt.59350) then
            read (neqdsk,*) (rwtmp2(i),i=1,magpri67)
          elseif (irshot.lt.91000) then
            read (neqdsk,*) (rwtmp2(i),i=1,magpri67+magpri322)
          elseif (irshot.lt.100771) then
            read (neqdsk,*) (rwtmp2(i),i=1,magpri67+magpri322 &
                                                   +magprirdp)
          elseif (irshot.lt.124985) then
            read (neqdsk,*) (rwtmp2(i),i=1,magpri67+magpri322 &
                                         +magprirdp+magudom)
          else
            read (neqdsk,*) (rwtmp2(i),i=1,magpri)
          endif
          go to 3050
        endif
        close(unit=neqdsk)
      endif
!
      time(1:ktime)=time(1:ktime)*1000.
      if (lookfw.gt.0) then
        do i=1,magpri
           if (fwtmp2(i).gt.0.0) fwtmp2(i)=rwtmp2(i)
        enddo
        do i=1,nsilop
           if (fwtsi(i).gt.0.0) fwtsi(i)=rwtsi(i)
        enddo
      endif
      swtmp2(1:magpri)=fwtmp2(1:magpri)
      swtsi(1:nsilop)=fwtsi(1:nsilop)
!-----------------------------------------------------------------------
!--   Get edge pedestal tanh paramters                                --
!-----------------------------------------------------------------------
      if(fitzts.eq.'te') then
#if defined(USE_MDS)
        call gettanh(ishot,fitzts,ktime,time,ztssym,ztswid, &
                     ptssym,ztserr)
#else
        kerror = 1
        call errctrl_msg('write_K','Cannot fit pedestal without MDS+')
        return
#endif
      endif
!
      do jtime=1,ktime
        itime=time(jtime)
        timems=itime
        timeus=(time(jtime)-timems)*1000.
        itimeu=timeus
!-----------------------------------------------------------------------
!--     correction for truncation                                     --
!-----------------------------------------------------------------------
        if (itimeu.ge.990) then
          itime=itime+1
          itimeu=0
          time(jtime)=itime
        endif
        coils(1:nsilop)=silopt(jtime,1:nsilop)
        fwtsi(1:nsilop)=swtsi(1:nsilop)
        do i=1,nsilop
          if (lookfw.gt.0.and.fwtsi(i).gt.0.0) fwtsi(i)=rwtsi(i)
          if (ierpsi(i).ne.0) fwtsi(i)=0.0
        enddo
        expmp2(1:magpri)=expmpi(jtime,1:magpri)
        fwtmp2(1:magpri)=swtmp2(1:magpri)
        do i=1,magpri
          if (lookfw.gt.0.and.fwtmp2(i).gt.0.0) fwtmp2(i)=rwtmp2(i)
          if (iermpi(i).ne.0) fwtmp2(i)=0.0
        enddo
        brsp(1:nfcoil)=fccurt(jtime,1:nfcoil)
        fwtfc(1:nfcoil)=swtfc(1:nfcoil)
        do i=1,nfcoil
          if (ierfc(i).ne.0) fwtfc(i)=0.0
        enddo
        ecurrt=eccurt(jtime,:)
        fwtec(1:nesum)=swtec(1:nesum)
        do i=1,nesum
          if (ierec(i).ne.0) fwtec(i)=0.0
        enddo
        denr(1:nco2r)=denrt(jtime,1:nco2r)
        denv(1:nco2v)=denvt(jtime,1:nco2v)
        if (mmstark.gt.0) then
          tgamma(1:nmtark)=tangam(jtime,1:nmtark)
          tgammauncor(1:nmtark)=tangam_uncor(jtime,1:nmtark)
          sgamma(1:nmtark)=siggam(jtime,1:nmtark)
          rrrgam(1:nmtark)=rrgam(jtime,1:nmtark)
          zzzgam(1:nmtark)=zzgam(jtime,1:nmtark)
          aa1gam(1:nmtark)=a1gam(jtime,1:nmtark)
          aa2gam(1:nmtark)=a2gam(jtime,1:nmtark)
          aa3gam(1:nmtark)=a3gam(jtime,1:nmtark)
          aa4gam(1:nmtark)=a4gam(jtime,1:nmtark)
          aa5gam(1:nmtark)=a5gam(jtime,1:nmtark)
          aa6gam(1:nmtark)=a6gam(jtime,1:nmtark)
          aa7gam(1:nmtark)=a7gam(jtime,1:nmtark)
          fwtgam(1:nmtark)=swtgam(1:nmtark)
          do i=1,nmtark
            if (iergam(i).ne.0) fwtgam(i)=0.0
          enddo
        endif
!
        if (mmbmsels.gt.0) then
          bmsels(1:nmsels)=bmselt(jtime,1:nmsels)
          sbmsels(1:nmsels)=sbmselt(jtime,1:nmsels)
          fwtbmsels(1:nmsels)=swtbmsels(1:nmsels)
          rrmsels(1:nmsels)=rrmselt(jtime,1:nmsels)
          zzmsels(1:nmsels)=zzmselt(jtime,1:nmsels)
          l1msels(1:nmsels)=l1mselt(jtime,1:nmsels)
          l2msels(1:nmsels)=l2mselt(jtime,1:nmsels)
          l4msels(1:nmsels)=l4mselt(jtime,1:nmsels)
          emsels(1:nmsels)=emselt(jtime,1:nmsels)
          semsels(1:nmsels)=semselt(jtime,1:nmsels)
          fwtemsels(1:nmsels)=swtemsels(1:nmsels)
          do i=1,nmsels
            if (iermselt(jtime,i).ne.0) then
              fwtbmsels(i)= 0.0
              fwtemsels(i)= 0.0
            endif
          enddo
        endif
!
        fwtcur=swtcur
        if(ierpla.ne.0) fwtcur=0.0
        btor=bcentr(jtime)
        plasma=pasmat(jtime)
        siref=psiref(jtime)
        vloop=vloopt(jtime)
        dflux=1.0e+03_dp*diamag(jtime)
        sigdlc=1.0e+03_dp*sigdia(jtime)
        pnbeam=pbinj(jtime)
        if(n1coil.gt.0) currn1=curtn1(jtime)
        if (nccoil.gt.0) then
          currc79=curc79(jtime)
          currc139=curc139(jtime)
          currc199=curc199(jtime)
        endif
        curriu30=curiu30(jtime)
        curriu90=curiu90(jtime)
        curriu150=curiu150(jtime)
        curril30=curil30(jtime)
        curril90=curil90(jtime)
        curril150=curil150(jtime)
!-----------------------------------------------------------------------
!--     Set edge pedestal tanh paramters                              --
!-----------------------------------------------------------------------
        if (fitzts.eq.'te'.and.ztserr(jtime)) then
          nbdry=1
          rbdry(1)=1.94_dp
          zbdry(1)=ztssym(jtime)+0.5_dp*ztswid(jtime)
        endif
        call setfnmeq(itimeu,'k',ishot,itime,eqdsk)
        open(unit=neqdsk,status='old',file=eqdsk,iostat=ioerr)
        if(ioerr.eq.0) close(unit=neqdsk,status='delete')
!-----------------------------------------------------------------------
!--     Set bit noise for ishot > 152000                              --
!-----------------------------------------------------------------------
!
!05-03-2016 Bob Johnson - problem seen on iris but not sun
!The next line gets rid of valgrind complaints on the
!following if statement checking ishot
        if(ishot.gt.152000) vbit = 80
        open(unit=neqdsk,file=eqdsk,status='new', &
             delim='quote')
!             delim='APOSTROPHE')
        write (neqdsk,in1)
        write (neqdsk,inwant)
        if(isetfb.ne.0) write (neqdsk,ink)
        if(mmstark.gt.0) write (neqdsk,ins)
        if(mmbmsels.gt.0) write (neqdsk,in_msels)
        if(kwaitmse.ne.0) write (neqdsk,ina)
        if(kfitece.gt.0) write (neqdsk,inece)
        if(keecur.gt.0) write (neqdsk,iner)
!-----------------------------------------------------------------------
!--     fitting type flag                                             --
!-----------------------------------------------------------------------
        if ((kprfit.gt.0).and.(mmstark.gt.0)) then
          fit_type = 'KIM'
        elseif (kprfit.gt.0) then
          fit_type = 'KIN'
        elseif (mmstark.gt.0) then
          fit_type = 'MSE'
          if (mmbmsels.gt.0) fit_type = 'MSL'
        elseif (mmbmsels.gt.0) then
          fit_type = 'MLS'
        else
          fit_type = 'MAG'
        endif
!
        header = ' '
        write (neqdsk,4042) header,fit_type
!---------------------------------------------------------------------
!--     Append SNAP file                                            --
!---------------------------------------------------------------------
        if (appendsnap.eq.'K'.or.appendsnap.eq.'KG') then
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
            endif
            if (ioerr.eq.0) close (unit=nsnapf)
 9991       format (a)
          endif
        endif
!
        close(unit=neqdsk)
        if (shotfile_in.eq.'0' .or. shotfile_in.eq.'') &
          read (eqdsk,6700) prefix1,ishotime
        if (cmdfile_in.ne.'0' .and. cmdfile_in.ne.'') &
          write (nffile,4960) ishotime
      enddo
      if (shotfile_in.ne.'0' .and. shotfile_in.ne.'') go to 3046
      endif no_err
      if (cmdfile_in.ne.'0' .and. cmdfile_in.ne.'') then
        write (nffile,4962)
        close(unit=nffile)
      endif
      if (shotfile_in.ne.'0' .and. shotfile_in.ne.'') close(unit=nout)
 4042 format (1x,a42,1x,a3)
 4958 format ('#!/bin/csh -f')
 4960 format ('      runefit.sc k',a12)
 4962 format ('#',/,'exit')
 4970 format (2x,a,1x,a,1x,a,1x,a)
 6700 format (a1,a12)
20000 format (/,1x,'shot data not on disk')
      end subroutine write_k

!**********************************************************************
!>
!!    Subroutine for writing K file during execution (with equilibrium)
!!    
!!
!!    @param jtime : time index
!!
!!    @param kerror : error flag
!!
!**********************************************************************
      subroutine write_k2(jtime,kerror)
      use set_kinds
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
#if defined(USEMPI)
      include 'mpif.h'
#endif

      integer*4, intent(in) :: jtime
      integer*4, intent(out) :: kerror
      integer*4 i,ioerr,limitrss,mmstark,nbdryss
      real*8 plasma,btor,dflux,xltype,xltype_180,vloop,siref,ktear, &
             pnbeam,timeus,timems,currn1,currc79,currc139,currc199, &
             curriu30,curriu90,curriu150,curril30,curril90,curril150
      real*8 kffcurt,kppcurt
      character eqdsk*20,header*42,fit_type*3
      integer*4, dimension(npcurn) :: kffbdryss,kff2bdryss, &
                                      kppbdryss,kpp2bdryss, &
                                      kwwbdryss,kww2bdryss
      real*8, dimension(npcurn) :: ffbdryss,ff2bdryss, &
                                   ppbdryss,pp2bdryss, &
                                   wwbdryss,ww2bdryss
      real*8, dimension(mbdry) :: rbdryss,zbdryss
      real*8 :: brspss(nrsmat),coils(nsilop),expmp2(magpri), &
                denr(nco2r),denv(nco2v)
      real*8, dimension(nmtark) :: tgamma,sgamma,rrrgam,zzzgam, &
                                   aa1gam,aa2gam,aa3gam,aa4gam,aa5gam, &
                                   aa6gam,aa7gam,tgammauncor

      character*82 snap_ext
      namelist/in1/ishot,itime,itimeu,qvfit,plasma,expmp2,coils,btor, &
           fwtsi,fwtcur,limitr,fwtmp2,kffcur,kppcur,fwtqa,ierchk, &
           fwtbp,serror,nextra,scrape,itrace,itek,xltype,rcentr,bitip, &
           psibit,bitmpi,denr,denv,siref,fwtfc,brsp,bitfc,iecurr,iplim, &
           ecurrt,ifitvs,vloop,dflux,ifcurr,iavem,icprof,currn1,n1coil, &
           pnbeam,error,errmin,mxiter,xltype_180,icutfp,keqdsk,ibtcomp, &
           fcurbd,pcurbd,kbound,alphafp,kskipvs,vsdamp,kframe,zelip, &
           fwtdlc,sigdlc,elomin,kcalpa,kcgama,calpa,cgama,xalpa,xgama, &
           kzeroj,rzeroj,iaveus,relax,fwtec,bitec,fitsiref, &
           kppfnc,kppknt,ppknt,pptens,kfffnc,kffknt,ffknt,fftens, &
           kwwfnc,kwwknt,wwknt,wwtens,nbdry,rbdry,zbdry,vbit, nbdrymx, &
           ppbdry,kppbdry,pp2bdry,kpp2bdry, &
           ffbdry,kffbdry,ff2bdry,kff2bdry, &
           wwbdry,kwwbdry,ww2bdry,kww2bdry, &
           ktear,kersil,iout,ixray,table_dir,input_dir,store_dir, &
           kpphord,kffhord,keehord,psiecn,dpsiecn,fitzts,isolve, &
           iplcout,imagsigma,errmag,saimin,errmagb,fitfcsum,fwtfcsum, &
           kwripre,ifindopt,tolbndpsi,siloplim,efitversion,use_previous
      namelist/inwant/psiwant,vzeroj,nccoil,currc79,currc139,rexpan, &
           znose,sizeroj,fitdelz,relaxdz,errdelz,oldccomp,nicoil, &
           oldcomp,currc199,curriu30,curriu90, &
           curriu150,curril30,curril90,curril150,ifitdelz,scaledz
      namelist/ink/isetfb,ioffr,ioffz,ishiftz,gain,gainp,idplace, &
           symmetrize,backaverage
      namelist/ins/tgamma,sgamma,fwtgam,rrrgam,zzzgam,aa1gam,aa2gam, &
           aa3gam,aa4gam,aa5gam,aa6gam,aa7gam,msebkp, &
           msefitfun,mse_quiet,mse_spave_on,kwaitmse, &
           dtmsefull,mse_strict,t_max_beam_off,ok_30rt,ok_210lt,&
           mse_usecer,mse_certree,mse_use_cer330,mse_use_cer210,&
           tgammauncor,v30lt,v30rt,v210lt,v210rt
      namelist/ina/spatial_avg_gam
      namelist/inece/necein,teecein0,feece0,errorece0,fwtece0,fwtecebz0, &
           ecefit,ecebzfit,kfitece,kinputece,kcallece,nharm, &
           kfixro,rteo,zteo,kfixrece,rtep,rtem,rpbit,rmbit,robit, &
           nfit,kcmin,fwtnow,mtxece
      namelist/iner/keecur,ecurbd,keefnc,eetens,keebdry,kee2bdry, &
           eebdry,ee2bdry,eeknt,keeknt,keehord
      data currn1/0.0/,currc79/0.0/,currc139/0.0/,currc199/0.0/, &
                       curriu30/0.0/,curriu90/0.0/,curriu150/0.0/, &
                       curril30/0.0/,curril90/0.0/,curril150/0.0/

      kerror = 0
!
      itime=time(jtime)
      timems=itime
      timeus=(time(jtime)-timems)*1000.
      itimeu=timeus
      siref=psiref(jtime)
      coils(1:nsilop)=silopt(jtime,1:nsilop) !-siref
      fwtsi(1:nsilop)=swtsi(1:nsilop)
      expmp2(1:magpri)=expmpi(jtime,1:magpri)
      fwtmp2(1:magpri)=swtmp2(1:magpri)
      brspss=brsp
      brsp=0.0
      brsp(1:nfcoil)=fccurt(jtime,1:nfcoil)
      fwtfc(1:nfcoil)=swtfc(1:nfcoil)
      ecurrt=eccurt(jtime,:)
      fwtec(1:nesum)=swtec(1:nesum)
      denr(1:nco2r)=denrt(jtime,1:nco2r)
      denv(1:nco2v)=denvt(jtime,1:nco2v)
      mmstark=0
      do i=1,nstark
        if (swtgam(i).gt.1.e-06_dp) mmstark=mmstark+1
      enddo
      if (mmstark.gt.0) then
        tgamma(1:nmtark)=tangam(jtime,1:nmtark)
        tgammauncor(1:nmtark)=tangam_uncor(jtime,1:nmtark)
        sgamma(1:nmtark)=siggam(jtime,1:nmtark)
        rrrgam(1:nmtark)=rrgam(jtime,1:nmtark)
        zzzgam(1:nmtark)=zzgam(jtime,1:nmtark)
        aa1gam(1:nmtark)=a1gam(jtime,1:nmtark)
        aa2gam(1:nmtark)=a2gam(jtime,1:nmtark)
        aa3gam(1:nmtark)=a3gam(jtime,1:nmtark)
        aa4gam(1:nmtark)=a4gam(jtime,1:nmtark)
        aa5gam(1:nmtark)=a5gam(jtime,1:nmtark)
        aa6gam(1:nmtark)=a6gam(jtime,1:nmtark)
        aa7gam(1:nmtark)=a7gam(jtime,1:nmtark)
        fwtgam(1:nmtark)=swtgam(1:nmtark)
      endif
      fwtcur=swtcur
      btor=bcentr(jtime)
      plasma=pasmat(jtime)
      nbdryss=nbdry
      rbdryss=rbdry
      zbdryss=zbdry
      rbdry=0.0
      zbdry=0.0
      if (fitzts.eq.'te'.and.ztserr(jtime)) then
        nbdry=1
        rbdry(1)=1.94_dp
        zbdry(1)=ztssym(jtime)+0.5_dp*ztswid(jtime)
      endif
      ppbdryss=ppbdry
      pp2bdryss=pp2bdry
      ffbdryss=ffbdry
      ff2bdryss=ff2bdry
      wwbdryss=wwbdry
      ww2bdryss=ww2bdry
      kppbdryss=kppbdry
      kpp2bdryss=kpp2bdry
      kffbdryss=kffbdry
      kff2bdryss=kff2bdry
      kwwbdryss=kwwbdry
      kww2bdryss=kww2bdry
      limitrss=limitr
!
      ppbdry=0.0
      pp2bdry=0.0
      ffbdry=0.0
      ff2bdry=0.0
      wwbdry=0.0
      ww2bdry=0.0
      kppbdry=0.0
      kpp2bdry=0.0
      kffbdry=0.0
      kff2bdry=0.0
      kwwbdry=0.0
      kww2bdry=0.0
      limitr=-limid
      vloop=vloopt(jtime)
      dflux=1.0e+03_dp*diamag(jtime)
      sigdlc=1.0e+03_dp*sigdia(jtime)
      pnbeam=pbinj(jtime)
      currn1=curtn1(jtime)
      currc79=curc79(jtime)
      currc139=curc139(jtime)
      currc199=curc199(jtime)
      curriu30=curiu30(jtime)
      curriu90=curiu90(jtime)
      curriu150=curiu150(jtime)
      curril30=curil30(jtime)
      curril90=curil90(jtime)
      curril150=curil150(jtime)
      itekt=itek
      mxitert=mxiter
      n1coilt=n1coil
      zeliptt=zelip
      kffcurt=kffcur
      kppcurt=kppcur
      itek=iteks
      mxiter=mxiters
      n1coil=n1coils
      zelip=zelipss
      kffcur=kffcurs
      kppcur=kppcurs
      if(ishot.gt.108281) n1coil = 0
!----------------------------------------------------------------------
!--   get xltype, xltype_180 without reading limiters                --
!----------------------------------------------------------------------
      call getlim(0,xltype,xltype_180,.false.)
!-----------------------------------------------------------------------
!--   Write K file                                                    --
!-----------------------------------------------------------------------
      call setfnmeq(itimeu,'k',ishot,itime,eqdsk)
      open(unit=neqdsk,status='old',file=eqdsk,iostat=ioerr)
      if(ioerr.eq.0) close(unit=neqdsk,status='delete')
      open(unit=neqdsk,file=eqdsk,status='new', &
           delim='quote')
!           delim='APOSTROPHE')
      write (neqdsk,in1)
      write (neqdsk,inwant)
      if(isetfb.ne.0) write (neqdsk,ink)
      if(mmstark.gt.0) write (neqdsk,ins)
      if(kwaitmse.ne.0) write (neqdsk,ina)
      if(kfitece.gt.0) write (neqdsk,inece)
      if(keecur.gt.0) write (neqdsk,iner)
!-----------------------------------------------------------------------
!--   fitting type flag                                               --
!-----------------------------------------------------------------------
      if ((kprfit.gt.0).and.(mmstark.gt.0)) then
        fit_type = 'KIM'
      elseif (kprfit.gt.0) then
        fit_type = 'KIN'
      elseif (mmstark.gt.0) then
        fit_type = 'MSE'
      else
        fit_type = 'MAG'
      endif
!
      header = ' '
      write (neqdsk,4042) header,fit_type
!----------------------------------------------------------------------
!--   Restore variables                                              --
!----------------------------------------------------------------------
      limitr=limitrss
      rbdry=rbdryss
      zbdry=zbdryss
      nbdry=nbdryss
      brsp=brspss
      ppbdry=ppbdryss
      pp2bdry=pp2bdryss
      ffbdry=ffbdryss
      ff2bdry=ff2bdryss
      wwbdry=wwbdryss
      ww2bdry=ww2bdryss
      kppbdry=kppbdryss
      kpp2bdry=kpp2bdryss
      kffbdry=kffbdryss
      kff2bdry=kff2bdryss
      kwwbdry=kwwbdryss
      kww2bdry=kww2bdryss
      itek=itekt
      mxiter=mxitert
      n1coil=n1coilt
      zelip=zeliptt
      kffcur=kffcurt
      kppcur=kppcurt
!---------------------------------------------------------------------
!--   Append SNAP file                                              --
!---------------------------------------------------------------------
      open(unit=nsnapf,status='old', &
           file=snap_file,iostat=ioerr)
      if (ioerr.eq.0) then
        snapfile=snap_file
      else
        open(unit=nsnapf,status='old', &
             file=input_dir(1:lindir)//snap_file,iostat=ioerr)
        if (ioerr.eq.0) then
          snapfile=input_dir(1:lindir)//snap_file
        elseif (snapextin.ne.'none') then
          open(unit=nsnapf,status='old', &
               file=snapextin)
          snapfile=snapextin
        endif
      endif
      if (appendsnap.eq.'K'.or.appendsnap.eq.'KG') then
        do i=1,1000000
          read (nsnapf,9991,iostat=ioerr) tmpdata
          if(ioerr.ne.0) exit
          if (INDEX(tmpdata,'&efitin')/=0) exit
        enddo
        if (ioerr.eq.0) then
          do i=1,1000000
            write (neqdsk,9991) tmpdata
            read (nsnapf,9991,iostat=ioerr) tmpdata
            if(ioerr.ne.0) exit
            if (INDEX(tmpdata,'/')/=0) then
              write (neqdsk,9991) tmpdata
              exit
            endif
          enddo
        endif
        if(ioerr.eq.0) close (unit=nsnapf)
 9991   format (a)
      endif
!
      close(unit=neqdsk)
 4042 format (1x,a42,1x,a3)
      end subroutine write_k2
