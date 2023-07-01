#include "config.f"
!**********************************************************************
!>
!!    setup_data_fetch performs setup and initialization, primarily
!!      for fetching data in snap mode
!!
!!                                                                  
!!    @param ktime : number of time slices requested
!!
!!    @param kerror : error flag
!!
!**********************************************************************
      subroutine setup_data_fetch(ktime,kerror)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
#if defined(USEMPI)
      include 'mpif.h'
#endif
      integer*4, intent(inout) :: ktime
      integer*4, intent(out) :: kerror
      integer*4 i,jtime,idtime,ioerr,istat,itimeb,iread,ireadold,irshot, &
                istop,mdoskip,mmstark,mmemsels,mtime,ktear
      real*8 xltype,xltype_180,dnmin,delt,times
      logical lopened,exists
      character filenm*15,ishotime*12,news*72,snap_ext*86, &
                eqdsk*20,comfile*15,prefix1*1,header*42,fit_type*3
      character(1000) line
      character(256) table_save
      real*8,dimension(:),allocatable :: fwtbmsels,fwtemsels
      ! DIIID specific parameters for reading fitweight.dat file
      integer*4, parameter :: nsilds=3,nsilol=41, &
                              magpri67=29,magpri322=31,magprirdp=8, &
                              magudom=5,maglds=3
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
           use_alternate_pointnames,alternate_pointname_file, &
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

      ALLOCATE(fwtbmsels(nmsels),fwtemsels(nmsels))
      fwtbmsels=0.0

      kerror=0
      mdoskip=0
      iout=1                 ! default - write fitout.dat
      appendsnap='KG'
      patmp2(1)=-1.
      ibatch=0 ! never used in code (just gets output)
      ilaser=0
      kffcurs=0
      kppcurs=0
      kwritime=0
      mtime=ktime
      table_save=table_dir

      ! initialize time dependent variables before loop
      rvsin=0.0
      zvsin=0.0
      rvsout=0.0
      zvsout=0.0
      tangam_uncor=0.0
      taumhd=0.0
      taudia=0.0
      chisq=0.0

! ONLY root process can check for existence of fitout.dat file
      mpi_rank0: if (rank == 0) then
        ! Delete fitout.dat if already exists
        open(unit=nout,status='old',file='fitout.dat',iostat=ioerr)
        if (ioerr.eq.0) close(unit=nout,status='delete')
      endif mpi_rank0
      if (iand(iout,1).ne.0) then
#if defined(USEMPI)
        if (nproc > 1) then
          if (rank == 0) then
            open(unit=nout,status='new',file='fitout.dat',delim='quote')
          endif
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          if (rank > 0) then
            open(unit=nout,status='old',file='fitout.dat',delim='quote')
          endif
        else
          open(unit=nout,status='new',file='fitout.dat',delim='quote')
        endif
#else
        open(unit=nout,status='new',file='fitout.dat',delim='quote')
#endif
      endif

#if defined(USE_SNAP)
      snap: if ((kdata.ne.1).and.(kdata.ne.2)) then
!----------------------------------------------------------------------
!--  look up magnetic data directly                                  --
!----------------------------------------------------------------------
      ktear=0
      xltype=0.0
      xltype_180=0.
!
      call set_defaults
!
      select case (kdata)
      case (3,7)
        if (snapextin.eq.'none') then
          snap_file = 'efit_snap.dat'
        else
          snap_ext = adjustl(snapextin)
          snap_file = 'efit_snap.dat_'//snap_ext
        endif
        open(unit=neqdsk,status='old', &
             file=snap_file,iostat=ioerr)
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
            snap_file = 'efit_snap.dat'
            open(unit=neqdsk,status='old',file=snap_file,iostat=ioerr)
            if (ioerr.ne.0) then
              call errctrl_msg('setup_data_fetch', &
                               'could not find snap file')
              stop
            endif
            snapfile=snap_file
          endif
        endif
      case (4)
        open(unit=neqdsk,status='old', &
             file='efit_time.dat',iostat=ioerr)
        if(ioerr.ne.0) &
          open(unit=neqdsk,status='old', &
               file=input_dir(1:lindir)//'efit_time.dat')
      end select

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
      close(unit=neqdsk)
      if (link_efit(1:1).ne.'') then
        ! use directories specified in efit.setup if available
        table_dir=trim(link_efit)//'green/'
        input_dir=trim(link_efit)
      elseif (table_dir.ne.table_save .and. link_efit.eq.'') then
        ! error if table_dir changes (global_allocs already set)
        call errctrl_msg('setup_data_fetch', &
          'changing machine during run is not supported (table_dir)')
        kerror=1
        return
      endif
      if(link_store(1:1).ne.'') store_dir=trim(link_store)
!--   warn that idebug, jdebug, and ktear inputs are deprecated
      if(idebug.ne.0) write(*,*) &
      "idebug input variable is deprecated, set cmake variable instead"
      if(jdebug.ne."NONE") write(*,*) &
      "jdebug input variable is deprecated, set cmake variable instead"
      if(ktear.ne.0) write(*,*) &
      "tearing calculations don't exist, ktear is deprecated"
!----------------------------------------------------------------------
!--   writes out the efitin namelist. Flag iout = 32.                --
!----------------------------------------------------------------------
      if (iand(iout,32).ne.0) then
        open(unit=nin,status='unknown',file='efit_snap.dat_out', &
             delim='quote',iostat=ioerr)
        if(ioerr.eq.0) write(nin,efitin)
        close(unit=nin)
      endif

      iteks=itek
      mxiters=mxiter
      zelipss=zelip
      kffcurs=kffcur
      kppcurs=kppcur
      n1coils=n1coil
      iexcals=iexcal
      ibunmns=ibunmn
      ierchks=ierchk
!
      ktime=mtime ! don't know why this would be wanted...
      if(qvfit.gt.0.0) qenp=qvfit ! iconvr=2 since not in efitin
      if (itek.lt.0) then
        ixray=1
        itek=iabs(itek)
      endif
!--------------------------------------------------------------------------
!--   itek > 100, write out PLTOUT.OUT individually                      --
!--------------------------------------------------------------------------
      kgraph=0
      if (itek.gt.100) then
        itek=itek-100
        kgraph=1
      endif
      itell=0
      if (mxiter.lt.0) then
        mxiter=-mxiter
        itell=1
        if(fitdelz) itell=4
      endif
      call set_basis_params
!
! -- Qilong Ren
      iishot = ishot
      kktime = ktime
!-------------------------------------------------------------------------------
!--   Set bit noise for ishot > 152000                                        --
!-------------------------------------------------------------------------------
      if(ishot.gt.152000) vbit = 80
!-------------------------------------------------------------------------------
!--   read in limiter data                                                    --
!-------------------------------------------------------------------------------
      call getlim(1,xltype,xltype_180,.false.)
!
  100 continue
      if (lookfw.ge.0) then
        rwtmp2(1:magpri)=0.0
        rwtsi(1:nsilop)=0.0
        open(unit=neqdsk,status='old', &
             file=input_dir(1:lindir)//'fitweight.dat')
  105   read(neqdsk,*,iostat=ioerr) irshot
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
          go to 105
        endif
        close(unit=neqdsk)
      endif
!
      times=timeb/1000.
      delt=dtime/1000.
      if (ifitvs.eq.1) then
        istop=-1
      else
        istop=0
      endif

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
        kerror = 1
        call errctrl_msg('setup_data_fetch','shot data not found')
        return
      endif

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
        call errctrl_msg('setup_data_fetch','MSE library not available')
        return
#endif
      endif
!
      do jtime=1,ktime
        swtbmselt(jtime,:)=fwtbmsels(:)
        swtemselt(jtime,:)=fwtemsels(:)
      enddo
      mmbmsels=0
      mmemsels=0
      swtbmsels=fwtbmsels
      do i=1,nmsels
        if(fwtbmsels(i).gt.1.e-06_dp) mmbmsels=mmbmsels+1
        if(fwtemsels(i).gt.1.e-06_dp) mmemsels=mmemsels+1
      enddo
      if(mmbmsels.gt.0) call getmsels(ktime)
!
      time(1:ktime)=time(1:ktime)*1000.
!-----------------------------------------------------------------------
!--   Get edge pedestal tanh paramters                                --
!-----------------------------------------------------------------------
      if (fitzts.eq.'te') then
#if defined(USE_MDS)
        call gettanh(ishot,fitzts,ktime,time,ztssym,ztswid, &
                     ptssym,ztserr)
#else
        kerror = 1
        call errctrl_msg('setup_data_fetch', &
                         'Cannot fit pedestal without MDS+')
        return
#endif
      endif
!----------------------------------------------------------------------
!--   Zero or replace bad fitting weights
!----------------------------------------------------------------------
      do i=1,nsilop
        if (ierpsi(i).ne.0) then
          fwtsi(i)=0.0
        elseif (lookfw.gt.0 .and. fwtsi(i).ne.0.0) then
          fwtsi(i)=rwtsi(i)
        endif
      enddo
      do i=1,nfsum
        if(ierfc(i).ne.0) fwtfc(i)=0.0
      enddo
      if (iecurr.eq.2) then
        do i=1,nesum
          if(ierec(i).ne.0) fwtec(i)=0.0
        enddo
      endif
      do i=1,magpri
        if (iermpi(i).ne.0) then
          fwtmp2(i)=0.0
        elseif (lookfw.gt.0 .and. fwtmp2(i).ne.0.0) then
          fwtmp2(i)=rwtmp2(i)
        endif
      enddo
      do i=1,nstark
        if(iergam(i).ne.0) fwtgam(i)=0.0
      enddo
      fwtece0=swtece
      do i=1,nnece
        if(ierece(i).ne.0) fwtece0(i)=0.0
      enddo
      fwtecebz0=swtecebz
      if(ierecebz.ne.0) fwtecebz0=0.0
      if(fwtqa.ne.0.0) fwtqa=1.
      if(fwtbp.ne.0.0) fwtbp=1.
      if(ierpla.ne.0) fwtcur=0.0
      if(ierrdi.ne.0) fwtdlc=0.0
!----------------------------------------------------------------------
!--   Save fitting weights
!----------------------------------------------------------------------
      swtdlc=fwtdlc
      swtcur=fwtcur
      swtfc=fwtfc
      swtec=fwtec
      swtmp2=fwtmp2
      swtsi=fwtsi
      swtgam=fwtgam
!
      call zlim(zero,nw,nh,limitr,xlim,ylim,rgrid,zgrid,limfag)
      endif snap
#endif
!----------------------------------------------------------------------
!--   Set grid parameters
!----------------------------------------------------------------------
      drgrid=rgrid(2)-rgrid(1)
      dzgrid=zgrid(2)-zgrid(1)
      darea=drgrid*dzgrid
      tmu2=-pi*tmu*dzgrid/drgrid
!
      return
      end subroutine setup_data_fetch
