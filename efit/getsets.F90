#include "config.f"
!**********************************************************************
!>
!!    getsets performs inputing and initialization.
!!
!!                                                                  
!!    @param ktime : number of time slices requested
!!
!!    @param mtear :
!!
!!    @param kerror : error flag
!!
!**********************************************************************
      subroutine getsets(ktime,mtear,kerror)
      use set_kinds
      use Fortran_Sleep
      use, intrinsic :: iso_c_binding, only: c_int
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
#if defined(USEMPI)
      include 'mpif.h'
#endif
      integer (c_int) :: retwait

      logical lopened
      character filenm*15,ishotime*12,news*72, &
                eqdsk*20,comfile*15,prefix1*1,header*42,fit_type*3
      real*8,dimension(:),allocatable :: fwtbmsels,fwtemsels
      character*82 snap_ext
      character(256) table_save
      character(len=1000) :: line
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
           do_spline_fit,table_dir,input_dir,store_dir,kedgep, &
           pedge,pe_psin,pe_width,kedgef,f2edge,fe_psin,fe_width, &
           psiecn,dpsiecn,relaxdz,fitzts,isolve,stabdz, &
           iplcout,errdelz,imagsigma,errmag,ksigma,saimin,errmagb, &
           write_Kfile,fitfcsum,fwtfcsum,appendsnap, &
           mse_quiet,mse_spave_on,kwaitmse,dtmsefull, &
           mse_strict,t_max_beam_off,ifitdelz,scaledz, &
           mse_usecer,mse_certree,mse_use_cer330,mse_use_cer210, &
           ok_30rt,ok_210lt,vbit,nbdrymx,fwtbmsels,fwtemsels,idebug,jdebug, &
           synmsels,avemsels,kwritime,v30lt,v30rt,v210lt,v210rt,ifindopt,tolbndpsi
      namelist/efitink/isetfb,ioffr,ioffz,ishiftz,gain,gainp,idplace, &
           symmetrize,backaverage,lring
      data mcontr/35/,lfile/36/,ifpsi/0/
      data currn1/0.0/,currc79/0.0/,currc139/0.0/,currc199/0.0/, &
                       curriu30/0.0/,curriu90/0.0/,curriu150/0.0/, &
                       curril30/0.0/,curril90/0.0/,curril150/0.0/
      logical exists
      integer*4, intent(inout) :: kerror

      ALLOCATE(fwtbmsels(nmsels),fwtemsels(nmsels))
      fwtbmsels=0.0

      kerror=0
      mdoskip=0
      iout=1                 ! default - write fitout.dat
      appendsnap='KG'
      patmp2(1)=-1.
      tmu0=twopi*tmu
      tmu02=tmu0*2.0
      errorm=1.
      ibatch=0 ! never used in code (just gets output)
      ilaser=0
      kffcurs=0
      kppcurs=0
      mtime=ktime
      table_save=table_dir
!----------------------------------------------------------------------
!--   news and help information                                      --
!----------------------------------------------------------------------
      ! ONLY root process displays EFIT news and help information
      mpi_rank: if (rank == 0) then
        open(unit=80,status='old', &
             file=input_dir(1:lindir)//'efithelp.txt',err=83220)
        do i=1,100
          read (80,83210,end=83220,err=83220) news
          write (nttyo,83210) news
        enddo
        close(unit=80)
      endif mpi_rank
83210 format (a)
!
83220 continue
#if defined(USE_SNAP)
!----------------------------------------------------------------------
!--   K-file from snap mode                                          --
!----------------------------------------------------------------------
      if (kdata.eq.5.or.kdata.eq.6) then
        call write_K(ktime,kerror)
        !if (kerror.gt.0) return ! don't return here because we're stopping anyway
        call errctrl_msg('getsets','Done writing k-files',3)
#if defined(USEMPI)
        call mpi_finalize(ierr)
#endif
        stop
      endif
#endif
      if (kwake.eq.1.and.mdoskip.eq.0.and.(iand(iout,1).ne.0)) close(unit=nout)

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
      not_kwake: if ((kwake.ne.1).or.(mdoskip.ne.0)) then
!----------------------------------------------------------------------
!--  look up magnetic data directly                                  --
!----------------------------------------------------------------------
      psibit(1:nsilop)=0.0
      fwtsi(1:nsilop)=0.0
      bitmpi(1:magpri)=0.0
      fwtmp2(1:magpri)=0.0
      fwtgam(1:nstark)=0.0
      fwtece0(1:nnece)=0.0
      fwtecebz0=0.0
      backaverage=.false.
      bitip=0.0
      betap0=0.50_dp
      brsp(1)=-1.e+20_dp
      cfcoil=-1.
      cutip=80000.
      ecurrt(1:nesum)=0.0
      rsisec(1:nesum)=-1.
      emf=1.00
      emp=1.00
      enf=1.00
      enp=1.00
      error=1.0e-03_dp
      fbetap=0.0
      fbetat=0.0
      fcurbd=1.
      fcsum(1:nfcoil)=1.0
      fczero(1:nfcoil)=1.0
      fwtfc(1:nfcoil)=0.
      rsisfc(1:nfcoil)=-1.
      fwtec(1:nesum)=0.0
      fwtbdry(1:mbdry)=1.0
      fwtsol(1:mbdry)=1.0
      sigrbd(1:mbdry)=1.e10_dp
      sigzbd(1:mbdry)=1.e10_dp
      fli=0.0
      fwtbp=0.0
      fwtdlc=0.0
      fwtqa=0.0
      gammap=1.0e+10_dp
      iaved=5
      iavem=5
      iavev=10
      ibound=0
      ibunmn=3
      icinit=2
      icondn=-1
      iconsi=-1
      iconvr=2
      icprof=0
      icurrt=2
      icutfp=0
      idite=0
      iecoil=0
      ierchk=1
      iecurr=1
      iexcal=0
!jal 04/23/2004
      iplcout=0
      ifcurr=0
      ifitvs=0
      ifref=-1
      itimeu=0
      iplim=0
      iprobe=0
      iqplot=1
      isetfb=0
      idplace=0
      islve=0
      isumip=0
      itek=0
      itrace=1
      ivacum=0
      ivesel=0
      n1coil=0
      ibtcomp=1
      iweigh=0
      ixray=0
      ixstrt=1
      keqdsk=1
      kinput=0
      kffcur=1
      kppcur=3
      kprfit=0
      limfag=2
      limitr=-33
      lookfw=1
      mxiter=25
      nbdry=0
      ncstfp=1
      ncstpp=1
      nextra=1
      nxiter=1
      pcurbd=1.
      psibry=0.0
      qemp=0.0
      qenp=0.95_dp
      qvfit=0.95_dp
      scrape=0.030_dp
      serror=0.03_dp
      sidif=-1.0e+10_dp
      symmetrize=.false.
      xltype=0.0
      xltype_180=0.
      gammap=1./gammap
      gammaf=gammap
      rmaxis=rzero
      mtear=0
      ktear=0
      snapfile='none'
      nsnapf=66
! -- Qilong Ren
      write_Kfile=.false.
      fitfcsum=.false.
      ifindopt=2
      tolbndpsi=1.0e-12_dp
!----------------------------------------------------------------------
!--   Initialize istore = 0                                          --
!--   Central directory to collect EFIT results is the default       --
!--   directory. Otherwise in store_dir (default to /link/store/)    --
!----------------------------------------------------------------------
      istore = 0
!
      select case (kdata)
      case (3)
        open(unit=neqdsk,status='old', &
             file='efit_snap.dat',iostat=ioerr)
        if (ioerr.eq.0) then
          snapfile='efit_snap.dat'
        else
          open(unit=neqdsk,status='old',       &
               file= input_dir(1:lindir)//'efit_snap.dat')
          snapfile=input_dir(1:lindir)//'efit_snap.dat'
        endif
      case (4)
        open(unit=neqdsk,status='old', &
             file='efit_time.dat',iostat=ioerr)
        if(ioerr.ne.0) &
          open(unit=neqdsk,status='old', &
               file= input_dir(1:lindir)//'efit_time.dat')
      case (7)
!----------------------------------------------------------------------
!--     Snap-Extension mode                                          --
!----------------------------------------------------------------------
        snap_ext = adjustl(snapextin)
        snap_file = 'efit_snap.dat_'//snap_ext
        open(unit=neqdsk,status='old', &
             file=snap_file,iostat=ioerr)
        if (ioerr.eq.0) then
          snapfile=snap_file
        else
          open(unit=neqdsk,status='old', &
               file= input_dir(1:lindir)//snap_file,iostat=ioerr)
          if (ioerr.eq.0) then
            snapfile=input_dir(1:lindir)//snap_file
          else
            snap_file = snap_ext
            open(unit=neqdsk,status='old',file= snap_file)
            snapfile=snap_file
          endif
        endif
      end select

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
      close(unit=neqdsk)
!--warn that idebug and jdebug inputs are depreciated
      if(idebug.ne.0) write(*,*) &
      "idebug input variable is depreciated, set cmake variable instead"
      if(jdebug.ne."NONE") write(*,*) &
      "jdebug input variable is depreciated, set cmake variable instead"
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
!
      ktime=mtime ! don't know why this would be wanted...
      mtear=ktear
      if((iconvr.ne.3).and.(qvfit.gt.0.0)) qenp=qvfit 
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
      endif not_kwake
!---------------------------------------------------------------------
!--   wakeup mode (KDATA=16)                                        --
!---------------------------------------------------------------------
10999 continue
      if (kwake.eq.1) then
28000   inquire(file='wakeefit.dat',opened=lopened)
        if(lopened) close(unit=neqdsk)
        open(unit=neqdsk, file='wakeefit.dat', status='old', err=28002)
        go to 28005
28002   continue
        retwait=FortSleep(10)
        go to 28000
28005   iread=0
28010   read (neqdsk,*,end=28020,err=28000) ishot,timeb,dtime,ktime
        iread=iread+1
        if(iread.ge.ireadold+1) go to 28020
        go to 28010
28020   close(unit=neqdsk)
        if (iread.le.ireadold) then
          retwait=FortSleep(10)
          go to 28000
        endif
        if (ishot.lt.0) then
          kerror = 1
          call errctrl_msg('getsets','shot not found')
          return
        endif
        ireadold=iread
      endif

!
! -- Qilong Ren
      iishot = ishot
      kktime = ktime
!----------------------------------------------------------------------- 
!--   reset table dir if necessary                                    --
!-----------------------------------------------------------------------
      if (table_dir.ne.table_save) then
        call set_table_dir
        call read_eparmdud
        call get_eparmdud_dependents
        call efit_read_tables
      endif
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
  105   read (neqdsk,*,iostat=ioerr) irshot
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
          go to 105
        endif
        close(unit=neqdsk)
      endif
!
      times=timeb/1000.
      delt=dtime/1000.
      if (ifitvs.gt.0) then
        istop=-1
      else
        istop=0
      endif

#if defined(USEMPI)
      if (nproc == 1) then
        call getpts(ishot,times,delt,ktime,istop)
      else
        call getpts_mpi(ishot,times,delt,ktime,istop)
      endif
#else
      call getpts(ishot,times,delt,ktime,istop)
#endif
      if (istop.gt.0) then
        kerror = 1
        call errctrl_msg('getsets','shot data not found')
        return
      endif

      mmstark=0
      swtgam(1:nstark)=fwtgam(1:nstark)
      do i=1,nstark
        if (fwtgam(i).gt.1.e-06_dp) mmstark=mmstark+1
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
        call errctrl_msg('getsets','MSE library not available')
        return
#endif
      endif
!
      do jtime=1,ktime
        swtbmselt(jtime,1:nmsels)=fwtbmsels(1:nmsels)
        swtemselt(jtime,1:nmsels)=fwtemsels(1:nmsels)
      enddo
      mmbmsels=0
      mmemsels=0
      swtbmsels(1:nmsels)=fwtbmsels(1:nmsels)
      do i=1,nmsels
        if (fwtbmsels(i).gt.1.e-06_dp) mmbmsels=mmbmsels+1
        if (fwtemsels(i).gt.1.e-06_dp) mmemsels=mmemsels+1
      enddo
      if(mmbmsels.gt.0) &
        call getmsels(ktime)
!
      time(1:ktime)=time(1:ktime)*1000.
!-----------------------------------------------------------------------
!--   Get edge pedestal tanh paramters                                --
!-----------------------------------------------------------------------
      if(fitzts.eq.'te') then
#if defined(USE_MDS)
        call gettanh(ishot,fitzts,ktime,time,ztssym,ztswid, &
                     ptssym,ztserr)
#else
        kerror = 1
        call errctrl_msg('getpts','Cannot fit pedestal without MDS+')
        return
#endif
      endif
!----------------------------------------------------------------------
!--   save fitting weights for SNAP modes                            --
!----------------------------------------------------------------------
      swtdlc=fwtdlc
      swtcur=fwtcur
      swtfc(1:nfcoil)=fwtfc(1:nfcoil)
      swtec(1:nesum)=fwtec(1:nesum)
      if (lookfw.gt.0) then
        do i=1,magpri
           if(fwtmp2(i).gt.0.0) fwtmp2(i)=rwtmp2(i)
        enddo
      endif
      swtmp2(1:magpri)=fwtmp2(1:magpri)
      if (lookfw.gt.0) then
        do i=1,nsilop
           if(fwtsi(i).gt.0.0) fwtsi(i)=rwtsi(i)
        enddo
      endif
      swtsi(1:nsilop)=fwtsi(1:nsilop)
!
      call zlim(zero,nw,nh,limitr,xlim,ylim,rgrid,zgrid,limfag)
      endif snap
#endif
!----------------------------------------------------------------------
!--   read in the plasma response function                           --
!----------------------------------------------------------------------
!
      drgrid=rgrid(2)-rgrid(1)
      dzgrid=zgrid(2)-zgrid(1)
      darea=drgrid*dzgrid
      tmu2=-pi*tmu*dzgrid/drgrid
!
      return
      end subroutine getsets
