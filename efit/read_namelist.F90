#include "config.f"
!**********************************************************************
!>
!!    This subroutine read efit.input file
!!
!**********************************************************************
      subroutine read_optin()
      use commonblocks
      use set_kinds
      use opt_input
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
#if defined(USEMPI)
      include 'mpif.h'
#endif
      integer*4 maxinpfile,mode,shot,steps,istat
      real*8 starttime,deltatime
      logical input_flag
      character cmdfile*15,shotfile*15,snapext*86
      character(1000) :: line

      namelist/setup/link_efit,link_store,maxinpfile,ntims,kxiter
      namelist/optin/mode,cmdfile,shotfile,shot,starttime,deltatime, &
                     steps,snapext,inpfile

      ! Initialize variables
      input_flag = .false.
      link_efit = ''
      link_store = ''
      maxinpfile = 1001
      ntims = 8192 ! sufficient for ms data from 8s shot
      use_opt_input = .false.
      mode = 99
      shot = -1
      steps = -1
      starttime = -1
      deltatime = -1
      cmdfile = '0'
      shotfile = '0'
      snapext = 'none'
      kxiter = -1
      kxiter_save = -1

      ! Determine if input file exists
      if(rank == 0) inquire(file='efit.input',exist=input_flag)
#if defined(USEMPI)
      if(nproc > 1) &
        call MPI_BCAST(input_flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
#endif
      if(.not.input_flag) return

      if(rank == 0) write(*,*) 'Reading efit.input file'
      ! All processes open and read efit.input file
      open(unit=nin,status='old',file='efit.input')

      ! Read setup and initialize input variables
      read (nin,setup,iostat=istat)
      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist setup: '//trim(line)
        stop
      endif
      if(rank == 0) write(*,*) & 
        'Using setup (if present) from efit.input file'
      allocate(inpfile(maxinpfile))
      if(kxiter.ne.-1) kxiter_save=kxiter

      ! Read inputs
      rewind(nin)
      read (nin,optin,iostat=istat)
      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist optin: '//trim(line)
        stop
      endif
      close(nin)

      ! check that variables are set (namelist may not exist or be used)
      if(mode.eq.99) return

      if(rank == 0) write(*,*) 'Using inputs from efit.input file'
      use_opt_input = .true.
      mode_in = mode
      cmdfile_in = cmdfile
      shotfile_in = shotfile
      shot_in = shot
      starttime_in = starttime
      deltatime_in = deltatime
      steps_in = steps
      snapext_in = snapext
      return

      end subroutine read_optin

!**********************************************************************
!>
!!    This subroutine reads experiment dependent parameters
!!
!**********************************************************************
      subroutine read_machinein()
      use var_nio, only : nin
      use errlims
      include 'eparm.inc'
      implicit none
      integer*4 :: istat,kubics
      character(1000) :: line

      ! it would probably be more clear to call this namelist experimentin or eparmin
      namelist/machinein/nsilds,nsilol,nfcoil,nrogow,nacoil,mfcoil,necoil,nvesel, &
        mpress,nesum,magpri67,magpri322,magprirdp,magudom,maglds,mse315,mse45, &
        mse15,mse1h,mse315_2,mse210,libim,nmsels,nnece,nnecein,neceo,nnnte, &
        ngam_vars,ngam_u,ngam_w,nlimit,nlimbd,nangle,ntangle,nfbcoil,mccoil, &
        micoil,ndata,nwwcur,nffcur,nppcur,nercur,ntime,ndim,kxiter,mqwant, &
        mbdry,mbdry1,nxtram,nxtlim,nco2v,nco2r,modef,modep,modew,kubics, &
        icycred_loopmax,nfourier,device, &
        ali_upper,ali_lower,betap_lim,plasma_diff, &
        aout_upper,aout_lower,eout_upper,eout_lower, &
        rout_upper,rout_lower,zout_upper,zout_lower, &
        rcurrt_upper,rcurrt_lower,zcurrt_upper,zcurrt_lower, &
        qsta_upper,qsta_lower,betat_lim, &
        oleft_lim,oright_lim,otop_lim, &
        olefs_check,qout_upper,qout_lower, &
        dbpli_lim,delbp_lim

      ! it would probably be more clear to call this experiment.dat or eparm.dat
      open(unit=nin,status='old',file=table_di2(1:ltbdi2)//'dprobe.dat')
      read (nin,machinein,iostat=istat)
      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist machinein: '//trim(line)
        stop
      endif
      close(unit=nin)
      
      if(kxiter_save.ne.-1) kxiter=kxiter_save

      return
      end subroutine read_machinein

!**********************************************************************
!>
!!    this subroutine grabs necessary information from file to properly
!!    set up directory paths
!!
!**********************************************************************
      subroutine read_dirs_shot(filename)
      use var_exdata, only: ishot,ifitvs
      use var_cecoil, only: iecurr
      use var_vessel, only: ivesel
      use var_input, only: icutfp
      use extvars, only: table_dir,input_dir,store_dir,efitversion
      implicit none
      character (*), intent(in) :: filename
      integer*4 nin,istat,itime,itek,itrace,nxiter,kffcur,kppcur,mxiter, &
                ierchk,limitr,nbdry,nslref,ibunmn,icurrt,icinit,iweigh, &
                iconvr,icprof,nextra,ixstrt,itimeu,islve,icntour,iprobe, &
                ifref,isumip,n1coil,ifcurr,iecoil,iplim,kinput,limfag, &
                kprfit,npress,keqdsk,npteth,nption,npneth,nbeam,iexcal, &
                iconsi,kcalpa,kcgama,iacoil,limid,iqplot,nptionf, &
                iavem,ktear,ndokin,kskipvs,limvs,kpressb,kzeroj, &
                ibtcomp,klabel,nmass,iaveus,kwripre,kbound,kframe, &
                kppfnc,kppknt,kfffnc,kffknt,kwwfnc,kwwknt,nbskip, &
                kersil,iout,ixray,kedgef,kedgep,kautoknt,kakloop,kakiter, &
                kpphord,kffhord,keehord,isolve,iplcout,imagsigma,ksigma, &
                nbdrymx,nsol,kbetapr,nbdryp,idebug,ifindopt,npnef,nptef
      real*8 plasma,btor,fwtcur,fwtqa,qemp,error,serror,psibry, &
             bitip,qenp,fwtbp,relip,zelip,aelip,eelip,qvfit,fwtdlc, &
             betap0,emp,enp,scrape,errmin,rbound,fwacoil,rcentr,rzero, &
             gammap,cfcoil,salpha,srm,sbeta,co2cor,xltype,xltype_180, &
             vsgnemin,sgtemin,sgnethi,sgtethi,dnmin,dflux,sigdlc, &
             sigtii,sigprebi,fwtxx,zeffvs,vloop,siref, &
             currn1,relax,saimin,cutip,pnbeam,sgprmin,elomin,siloplim, &
             fcurbd,pcurbd,prbdry,zlowimp,errmag,errmagb, &
             pressbi,prespb,sigppb,rminvs,rmaxvs,zmaxvs,errbry,condin, &
             sgtimin,alphafp,zbound,vsdamp,zminvs,saicon,tolbndpsi, &
             pptens,fftens,wwtens,scalesir,scalea,errsil,vbit, &
             f2edge,fe_width,fe_psin,pedge,pe_width,pe_psin, &
             akchiwt,akerrwt,aktol,akgamwt,akprewt,psiecn,dpsiecn
      integer*4, dimension(2000) :: irfila,jzfila
      real*8, dimension(2000) :: expmp2,coils,fwtsi,fwtmp2,psibit,bitmpi, &
                                 denr,denv,fwtfc,acoilc,brsp,bitfc, &
                                 ecurrt,xalpa,xgama,rzeroj,fwtec,bitec, &
                                 ppknt,ffknt,wwknt,rbdry,zbdry, &
                                 ppbdry,kppbdry,pp2bdry,kpp2bdry, &
                                 ffbdry,kffbdry,ff2bdry,kff2bdry, &
                                 wwbdry,kwwbdry,ww2bdry,kww2bdry, &
                                 fwtfcsum,fczero,fcsum,fwtbdry,xlim,ylim, &
                                 rpress,zpress,pressr,sigpre, &
                                 fwtpre,tethom,rteth,zteth,sgteth,tionex, &
                                 rion,zion,dnethom,rneth,zneth, &
                                 sibeam,pbeam,dnbeam,dmass,vcurfb,vcurrt, &
                                 brsptu,sigti,sgneth,scalepr, &
                                 sigrbd,sigzbd,rsol,zsol,fwtsol, &
                                 chordv,chordr
      real*8, dimension(256,256) :: calpa,cgama
      character(1000) :: line, fitzts
      character appendsnap*2,jdebug*4
      logical :: fitsiref,fitfcsum,use_previous,req_valid
      namelist/in1/ishot,itime,plasma,itek,itrace,nxiter,fwtcur,kffcur, &
        coils,fwtsi,expmp2,fwtmp2,kppcur,mxiter,ierchk,fwtqa,qemp,error, &
        limitr,xlim,ylim,serror,nbdry,rbdry,zbdry,psibry,nslref,ibunmn, &
        btor,psibit,bitmpi,bitip,icurrt,icinit,brsp,iweigh,qenp,fwtbp, &
        relip,zelip,aelip,eelip,qvfit,fwtdlc,betap0,emp,enp,iconvr,icprof, &
        nextra,ixstrt,scrape,errmin,rbound,npnef,nptef,fwacoil,itimeu, &
        rcentr,rzero,gammap,cfcoil,fczero,fcsum,islve,icntour,iprobe, &
        salpha,srm,sbeta,ifref,isumip,n1coil,ifcurr,iecurr,ecurrt,iecoil, &
        co2cor,vsgnemin,sigtii,sgtemin,sgnethi,sgtethi,dnmin, &
        vcurrt,dflux,sigdlc,iplim,kinput,limfag,sigprebi,fwtxx, &
        kprfit,pressr,rpress,zpress,sigpre,npress,tethom,rteth,keqdsk, &
        zteth,sgteth,npteth,tionex,rion,zion,sigti,nption,dnethom,zeffvs, &
        rneth,zneth,sgneth,npneth,pbeam,sibeam,nbeam,rzeroj,xalpa,cgama, &
        ivesel,iexcal,iconsi,fwtfc,xltype,kcalpa,kcgama,calpa,iacoil, &
        limid,irfila,jzfila,vloop,iqplot,siref,denr,denv,xgama,&
        nptionf,currn1,ifitvs,bitfc,relax,saimin,icutfp,acoilc, &
        cutip,iavem,pnbeam,xltype_180,sgprmin,elomin,ktear, &
        fcurbd,pcurbd,prbdry,ndokin,zlowimp,kskipvs,limvs, &
        vcurfb,kpressb,pressbi,prespb,sigppb,kzeroj,rminvs,rmaxvs,errbry, &
        fwtpre,ibtcomp,klabel,zmaxvs,dnbeam,dmass,nmass,condin,iaveus, &
        sgtimin,kwripre,kbound,alphafp,kframe,zbound,vsdamp,zminvs,saicon, &
        kppfnc,kppknt,ppknt,pptens,kfffnc,kffknt,ffknt,fftens,fwtbdry, &
        kwwfnc,kwwknt,wwknt,wwtens,fwtec,fitsiref,bitec,scalepr,scalesir, &
        ppbdry,kppbdry,pp2bdry,kpp2bdry,scalea,sigrbd,sigzbd,nbskip, &
        ffbdry,kffbdry,ff2bdry,kff2bdry,errsil,vbit,kersil,iout,ixray, &
        wwbdry,kwwbdry,ww2bdry,kww2bdry,f2edge,fe_width,fe_psin,kedgef, &
        pedge,kedgep,pe_width,pe_psin,chordv,chordr, &
        kautoknt,akchiwt,akerrwt,kakloop,aktol,kakiter,akgamwt,akprewt, &
        kpphord,kffhord,keehord,psiecn,dpsiecn,fitzts,isolve,iplcout, &
        imagsigma,errmag,ksigma,errmagb,brsptu,fitfcsum,fwtfcsum,appendsnap, &
        nbdrymx,nsol,rsol,zsol,fwtsol,efitversion,kbetapr,nbdryp, &
        idebug,jdebug,ifindopt,tolbndpsi,siloplim,use_previous,req_valid, &
        table_dir,input_dir,store_dir
      parameter(nin=343)

      open(unit=nin,status='old',file=filename)
      read(nin,in1,iostat=istat)

      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist in1: '//trim(line)
        stop
      endif
      close(nin)
      end subroutine read_dirs_shot

!**********************************************************************
!>
!!    this subroutine grabs necessary information from file to properly
!!    set up directory paths from an IMAS equilibrium hdf5 file
!!
!!    timeslice 0 is used, but could be made variable
!!
!**********************************************************************
      subroutine read_dirs_shot_imas(filename)
      use var_exdata, only: ishot,ifitvs
      use var_cecoil, only: iecurr
      use var_vessel, only: ivesel
      use var_input, only: icutfp
      use var_nio
      use error_control
      use extvars, only: table_dir,input_dir,store_dir,efitversion
      implicit none
      character(*), intent(in) :: filename
      integer*4  :: istat
      integer*4 :: dims
      logical :: file_stat
 
#if defined(USE_HDF5)
      inquire(file=trim(filename),exist=file_stat)
      if (.not. file_stat) then
        call errctrl_msg('read_imas_in1',trim(filename)//' not found')
        stop
      endif
      call fch5init
      call open_oldh5file(trim(filename),fileid,rootgid,h5in,h5err)
      call test_group(rootgid,"equilibrium",file_stat,h5err)
      if (.not. file_stat) then
        call errctrl_msg('read_imas_in1','equilibrium group not found')
        stop
      endif
      call open_group(rootgid,"equilibrium",eqid,h5err)
      call test_group(eqid,"code",file_stat,h5err)
      if (.not. file_stat) then
        call errctrl_msg('read_imas_in1','code group not found')
        stop
      endif
      call open_group(eqid,"code",cid,h5err)
      call test_group(cid,"parameters",file_stat,h5err)
      if (.not. file_stat) then
        call errctrl_msg('read_imas_in1','parameters group not found')
        stop
      endif
      call open_group(cid,"parameters",pid,h5err)
      call test_group(pid,"time_slice",file_stat,h5err)
      if (.not. file_stat) then
        call errctrl_msg('read_imas_in1','time_slice group not found')
        stop
      endif
      call open_group(pid,"time_slice",tid,h5err)
      call test_group(tid,"0",file_stat,h5err)
      if (.not. file_stat) then
        call errctrl_msg('read_imas_in1','100 group not found')
        stop
      endif
      call open_group(tid,"0",sid,h5err)
      call test_group(sid,"in1",file_stat,h5err)
      if (file_stat) then
        call open_group(sid,"in1",nid,h5err)
        call read_h5_ex(nid,"ishot",ishot,h5in,h5err)
        call read_h5_ex(nid,"iecurr",iecurr,h5in,h5err)
        call read_h5_ex(nid,"ivesel",ivesel,h5in,h5err)
        call read_h5_ex(nid,"ifitvs",ifitvs,h5in,h5err)
        call read_h5_ex(nid,"icutfp",icutfp,h5in,h5err)
        call read_h5_ex(nid,"table_dir",table_dir,h5in,h5err)
        call read_h5_ex(nid,"input_dir",input_dir,h5in,h5err)
        call read_h5_ex(nid,"efitversion",efitversion,h5in,h5err)
        call close_group("in1",nid,h5err)
      endif
      call close_group("0",sid,h5err)
      call close_group("time_slice",tid,h5err)
      call close_group("parameters",pid,h5err)
      call close_group("code",cid,h5err)
      call close_group("equilibrium",eqid,h5err)
      call close_h5file(fileid,rootgid,h5err)

#else
      ! this code should not be reachable
      call errctrl_msg('read_imas_in1','HDF5 needs to be linked')
      stop
#endif
      end subroutine read_dirs_shot_imas
