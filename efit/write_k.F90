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
      real*8 coils(nsilop),expmp2(magpri),denr(nco2r),denv(nco2v)
      real*8,dimension(nmselp) :: tgamma,sgamma,rrrgam,zzzgam, &
                                  aa1gam,aa2gam,aa3gam,aa4gam, &
                                  aa5gam,aa6gam,aa7gam,tgammauncor
      real*8,dimension(nmsels) :: bmsels,sbmsels,fwtbmsels, &
                                  rrmsels,zzmsels,l1msels,l2msels, &
                                  l4msels,emsels,semsels,fwtemsels
      real*8 :: dnmin
      character*86 snap_ext
      character(1000) :: line
      character(256) table_save
      ! DIIID specific parameters for reading fitweight.dat file
      integer*4, parameter :: nsilds=3,nsilol=41, &
                              magpri67=29,magpri322=31,magprirdp=8, &
                              magudom=5,maglds=3
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
           efitversion,kwripre,ifindopt,tolbndpsi,siloplim,use_previous, &
           req_valid,nw_sub,nh_sub
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
           write_kfile,fitfcsum,fwtfcsum,appendsnap, &
           mse_quiet,mse_spave_on,kwaitmse,dtmsefull, &
           mse_strict,t_max_beam_off,ifitdelz,scaledz, &
           mse_usecer,mse_certree,mse_use_cer330,mse_use_cer210, &
           ok_30rt,ok_210lt,vbit,nbdrymx,fwtbmsels,fwtemsels, &
           idebug,jdebug,synmsels,avemsels,kwritime, &
           v30lt,v30rt,v210lt,v210rt,ifindopt,tolbndpsi, &
           siloplim,use_previous,ierchk,req_valid,ibunmn,nw_sub,nh_sub
      namelist/efitink/isetfb,ioffr,ioffz,ishiftz,gain,gainp,idplace, &
           symmetrize,backaverage,lring
      data currn1/0.0/,currc79/0.0/,currc139/0.0/,currc199/0.0/, &
                       curriu30/0.0/,curriu90/0.0/,curriu150/0.0/, &
                       curril30/0.0/,curril90/0.0/,curril150/0.0/

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
        ! the snap file wasn't found in the CWD, look in the 
        ! support files instead
        open(unit=neqdsk,status='old', &
             file=input_dir(1:lindir)//'snapfiles/'//snap_file,iostat=ioerr)
        if (ioerr.eq.0) then
          snapfile=input_dir(1:lindir)//'snapfiles/'//snap_file
        elseif (snapextin.ne.'none') then
          ! the snap file still wasn't found, check for any snapfile
          ! in the CWD without an extension
          snap_file = snap_ext
          open(unit=neqdsk,status='old',file=snap_file,iostat=ioerr)
          if (ioerr.ne.0) then
            call errctrl_msg('write_k','could not find snap file')
            stop
          endif
          snapfile=snap_file
        endif
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
      rewind(neqdsk)
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
          'changing machine during run is not supported (table_dir)')
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
 3046 continue
!
      if (shotfile_in.ne.'0' .and. shotfile_in.ne.'') then
        read(nout,4970,iostat=ioerr) ishotime
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
      if (ifitvs.eq.1 .or. ivesel.eq.3) then
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
        call get_measurements(ishot,times,delt,ktime,istop)
      else
        call get_measurements_mpi(ishot,times,delt,ktime,istop)
      endif
#else
      call get_measurements(ishot,times,delt,ktime,istop)
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
        if(fwtgam(i).gt.1.e-06_dp) mmstark=mmstark+1
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
      swtbmsels=fwtbmsels
      swtemsels=fwtemsels
#ifdef DEBUG_LEVEL2
      write (6,*) 'WRITE_K swtbmsels= ',(swtbmsels(i),i=1,nmsels)
#endif
      do i=1,nmsels
        if(fwtbmsels(i).gt.1.e-06_dp) mmbmsels=mmbmsels+1
        if(fwtemsels(i).gt.1.e-06_dp) mmemsels=mmemsels+1
      enddo
#ifdef DEBUG_LEVEL2
      write (6,*) 'WRITE_K mmbmsels= ', mmbmsels
#endif
      if(mmbmsels.gt.0) call getmsels(ktime)
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
 3050   read(neqdsk,*,iostat=ioerr) irshot
        if ((ioerr.eq.0).and.(irshot.le.ishot)) then
          ! Note: These DIIID specific probe choices
          ! could be removed if Green function tables are regenerated...
          if (irshot.lt.124985) then
            read(neqdsk,*) (rwtsi(i),i=1,nsilol)
          else
            read(neqdsk,*) (rwtsi(i),i=1,nsilop)
          endif
          if (irshot.lt.59350) then
            read(neqdsk,*) (rwtmp2(i),i=1,magpri67)
          elseif (irshot.lt.91000) then
            read(neqdsk,*) (rwtmp2(i),i=1,magpri67+magpri322)
          elseif (irshot.lt.100771) then
            read(neqdsk,*) (rwtmp2(i),i=1,magpri67+magpri322 &
                                                   +magprirdp)
          elseif (irshot.lt.124985) then
            read(neqdsk,*) (rwtmp2(i),i=1,magpri67+magpri322 &
                                         +magprirdp+magudom)
          else
            read(neqdsk,*) (rwtmp2(i),i=1,magpri)
          endif
          go to 3050
        endif
        close(unit=neqdsk)
      endif
!
      time(1:ktime)=time(1:ktime)*1000.
      do i=1,nsilop
        if (ierpsi(i).ne.0) then
          fwtsi(i)=0.0
        elseif (lookfw.gt.0.and.fwtsi(i).gt.0.0) then
          fwtsi(i)=rwtsi(i)
        endif
      enddo
      do i=1,magpri
        if (iermpi(i).ne.0) then
          fwtmp2(i)=0.0
        elseif (lookfw.gt.0 .and. fwtmp2(i).gt.0.0) then
          fwtmp2(i)=rwtmp2(i)
        endif
      enddo
      do i=1,nfsum
        if(ierfc(i).ne.0) fwtfc(i)=0.0
      enddo
      do i=1,nesum
        if(ierec(i).ne.0) fwtec(i)=0.0
      enddo
      if (mmstark.gt.0) then
        do i=1,nmselp
          if(iergam(i).ne.0) fwtgam(i)=0.0
        enddo
      endif
      if(ierpla.ne.0) fwtcur=0.0
!-----------------------------------------------------------------------
!--   Get edge pedestal tanh paramters                                --
!-----------------------------------------------------------------------
      if (fitzts.eq.'te') then
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
        call errctrl_setstate(rank,time(jtime))
        if (req_valid) then
          ! don't write times without mse
          if (sum(abs(fwtgam(1:nmselp))).gt.nmselp*1.e-6_dp .and. &
              sum(abs(siggam(jtime,1:nmselp))).lt.nmselp*1.e-10_dp) then
            call errctrl_msg('write_k', &
              'No MSE data found, not writing k-file')
            cycle
          endif
          ! don't write times without cer
          if (mse_usecer.ne.0 .and. &
              maxval(abs(tangam(jtime,1:nmselp) &
                        -tangam_uncor(jtime,1:nmselp))).le.1.e-10_dp) then
            call errctrl_msg('write_k', &
              'No CER correction used, not writing k-file')
            cycle
          endif
        endif
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
        coils=silopt(jtime,:)
        expmp2=expmpi(jtime,:)
        brsp(1:nfsum)=fccurt(jtime,:)
        ecurrt=eccurt(jtime,:)
        denr=denrt(jtime,:)
        denv=denvt(jtime,:)
        if (mmstark.gt.0) then
          tgamma=tangam(jtime,1:nmselp)
          tgammauncor=tangam_uncor(jtime,1:nmselp)
          sgamma=siggam(jtime,1:nmselp)
          rrrgam=rrgam(jtime,1:nmselp)
          zzzgam=zzgam(jtime,1:nmselp)
          aa1gam=a1gam(jtime,1:nmselp)
          aa2gam=a2gam(jtime,1:nmselp)
          aa3gam=a3gam(jtime,1:nmselp)
          aa4gam=a4gam(jtime,1:nmselp)
          aa5gam=a5gam(jtime,1:nmselp)
          aa6gam=a6gam(jtime,1:nmselp)
          aa7gam=a7gam(jtime,1:nmselp)
        endif
!
        if (mmbmsels.gt.0) then
          bmsels=bmselt(jtime,:)
          sbmsels=sbmselt(jtime,:)
          fwtbmsels=swtbmsels
          rrmsels=rrmselt(jtime,:)
          zzmsels=zzmselt(jtime,:)
          l1msels=l1mselt(jtime,:)
          l2msels=l2mselt(jtime,:)
          l4msels=l4mselt(jtime,:)
          emsels=emselt(jtime,:)
          semsels=semselt(jtime,:)
          fwtemsels=swtemsels
          do i=1,nmsels
            if (iermselt(jtime,i).ne.0) then
              fwtbmsels(i)= 0.0
              fwtemsels(i)= 0.0
            endif
          enddo
        endif
!
        btor=bcentr(jtime)
        plasma=ipmeas(jtime)
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
      integer*4 i,ioerr,limitrss,mmstark,nbdryss,kffcurt,kppcurt,ktear
      real*8 plasma,btor,dflux,xltype,xltype_180,vloop,siref, &
             pnbeam,timeus,timems,currn1,currc79,currc139,currc199, &
             curriu30,curriu90,curriu150,curril30,curril90,curril150
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
      real*8, dimension(nmselp) :: tgamma,sgamma,rrrgam,zzzgam, &
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
           kzeroj,rzeroj,iaveus,relax,fwtec,bitec,fitsiref,nw_sub,nh_sub, &
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
      coils=silopt(jtime,:) !-siref
      fwtsi=swtsi
      expmp2=expmpi(jtime,:)
      fwtmp2=swtmp2
      brspss=brsp
      brsp=0.0
      brsp(1:nfsum)=fccurt(jtime,:)
      fwtfc=swtfc
      ecurrt=eccurt(jtime,:)
      fwtec=swtec
      denr=denrt(jtime,:)
      denv=denvt(jtime,:)
      fwtgam=swtgam
      mmstark=0
      do i=1,nstark
        if (fwtgam(i).gt.1.e-06_dp) mmstark=mmstark+1
      enddo
      if (mmstark.gt.0) then
        tgamma=tangam(jtime,1:nmselp)
        tgammauncor=tangam_uncor(jtime,1:nmselp)
        sgamma=siggam(jtime,1:nmselp)
        rrrgam=rrgam(jtime,1:nmselp)
        zzzgam=zzgam(jtime,1:nmselp)
        aa1gam=a1gam(jtime,1:nmselp)
        aa2gam=a2gam(jtime,1:nmselp)
        aa3gam=a3gam(jtime,1:nmselp)
        aa4gam=a4gam(jtime,1:nmselp)
        aa5gam=a5gam(jtime,1:nmselp)
        aa6gam=a6gam(jtime,1:nmselp)
        aa7gam=a7gam(jtime,1:nmselp)
      endif
      fwtcur=swtcur
      btor=bcentr(jtime)
      plasma=ipmeas(jtime)
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
      kppbdry=0
      kpp2bdry=0
      kffbdry=0
      kff2bdry=0
      kwwbdry=0
      kww2bdry=0
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
