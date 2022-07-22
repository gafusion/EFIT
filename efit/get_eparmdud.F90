#include "config.f"
!**********************************************************************
!>
!!    sets default DIII-D values for eparmdud module
!!
!**********************************************************************
subroutine get_eparmdud_defaults()
  use eparm
  use var_cecoil, only: iecurr
  nsilds=3
  nsilol=41
  nfcoil=18
  nrogow=1
  nacoil=1
  mfcoil=18
  necoil=122
  nvesel=24
  mpress=201
  nesum=6
  magpri67=29
  magpri322=31
  magprirdp=8
  magudom=5
  maglds=3
  mse315=11
  mse45=15
  mse15=10
  mse1h=4
  mse315_2=5
  mse210=24
  libim=32
  nmsels=16
  nnece=40
  nnecein=80
  neceo=1
  nnnte=801
  ngam_vars=9
  ngam_u=5
  ngam_w=3
  nlimit=160
  nlimbd=6
  nangle=64
  ntangle=12
  nfbcoil=12
  mccoil=6
  micoil=12

  ndata=61
  nwwcur=32
  nffcur=32
  nppcur=32
  nercur=32

  ntime=1001
  ndim=3200
  kxiter=515
  mqwant=30
  mbdry=300
  mbdry1=110
  nxtram=10
  nxtlim=9
  nco2v=3
  nco2r=2
  modef=4
  modep=4
  modew=4
  kubics=4
  icycred_loopmax=1290
  nfourier=5

  iecurr=1

end subroutine get_eparmdud_defaults


!**********************************************************************
!>
!!    this subroutine calculates eparmdud variables that are other 
!!    eparmdud variables
!!
!**********************************************************************
subroutine get_eparmdud_dependents()

  use eparm

  nmtark=mse315+mse45+mse15+mse1h+mse315_2+mse210
  nstark=nmtark+libim

  magpol=magpri67+magpri322+magprirdp+magudom
  magpri=magpol+maglds

  nsilop=nsilds+nsilol

  nbwork=nsilop
  msbdry=mbdry+nsilop+nfcoil+1
  msbdr2=2*msbdry

  npcurn=nffcur+nppcur
  necur2=nercur*2
  mfnpcr=nfcoil+npcurn+nvesel+nwwcur+nesum+nfcoil+nercur 
  npcur2=npcurn*2
  nrsmat=nsilop+magpri+nrogow+nffcur+1+npcurn+nwwcur+mpress+nfcoil+nstark+nnece+neceo
  nrsma2=2*nrsmat

  nwcurn=nwwcur+npcurn
  nwcur2=nwcurn*2

  ncurrt=nvesel+nesum+nfcoil
end subroutine get_eparmdud_dependents


!**********************************************************************
!>
!!    this subroutine grabs necessary information from file to properly
!!    set up directory paths
!!
!**********************************************************************
subroutine read_dirs_shot(filename)
  use var_exdata, only: ishot
  use var_cecoil, only: iecurr
  use var_vessel, only: ivesel
  use exvars, only: table_dir,input_dir,store_dir,efitversion
  implicit integer*4 (i-n), real*8 (a-h,o-z)
  real*8, dimension(2000) :: expmp2,coils,fwtsi,fwtmp2,psibit,bitmpi,denr,denv,fwtfc, &
                            acoilc,brsp,bitfc,ecurrt,xalpa,xgama,rzeroj,fwtec,bitec, &
                            ppknt,ffknt,wwknt,rbdry,zbdry,ppbdry,kppbdry,pp2bdry,kpp2bdry, &
                            ffbdry,kffbdry,ff2bdry,kff2bdry,wwbdry,kwwbdry,ww2bdry,kww2bdry,&
                            fwtfcsum,fczero,fcsum,fwtbdry,xlim,ylim,rpress,zpress,pressr,sigpre,&
                            fwtpre,tethom,rteth,zteth,sgteth,tionex,rion,zion,dnethom,rneth,zneth, &
                            sibeam,pbeam,dnbeam,dmass,vcurfb,vcurrt,brsptu,sigti,sgneth,scalepr, &
                            sigrbd,sigzbd,rsol,zsol,fwtsol
  integer*4, dimension(2000) :: irfila,jzfila
  real*8, dimension(256,256) :: calpa,cgama
  real*8 :: siloplim
  integer*4 :: istat
  character(len=1000) :: line, fitzts
  character (*) :: filename
  character appendsnap*2,jdebug*4
  logical :: fitsiref, fitfcsum
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
       pedge,kedgep,pe_width,pe_psin, &
       kautoknt,akchiwt,akerrwt,kakloop,aktol,kakiter,akgamwt,akprewt, &
       kpphord,kffhord,keehord,psiecn,dpsiecn,fitzts,isolve,iplcout, &
       imagsigma,errmag,ksigma,errmagb,brsptu,fitfcsum,fwtfcsum,appendsnap, &
       nbdrymx,nsol,rsol,zsol,fwtsol,efitversion,kbetapr,nbdryp, &
       idebug,jdebug,ifindopt,tolbndpsi,siloplim,use_previous, &
       table_dir,input_dir,store_dir

  nin=343
  open(unit=nin,status='old',file=filename)
  read (nin,in1,iostat=istat)

  if (istat>0) then
    read(nin,fmt='(A)') line
    write(*,'(A)') &
      'Invalid current line in namelist: '//trim(line)
    backspace(nin)
    backspace(nin)
    read(nin,fmt='(A)') line
    write(*,'(A)') &
      'Invalid line in namelist: '//trim(line)
    backspace(nin)
    backspace(nin)
    read(nin,fmt='(A)') line
    write(*,'(A)') &
      'Previous line in namelist: '//trim(line)
    backspace(nin)
    backspace(nin)
    read(nin,fmt='(A)') line
    write(*,'(A)') &
      'Previous line in namelist: '//trim(line)
    backspace(nin)
    backspace(nin)
    read(nin,fmt='(A)') line
    write(*,'(A)') &
      'Previous line in namelist: '//trim(line)
    backspace(nin)
    backspace(nin)
    read(nin,fmt='(A)') line
    write(*,'(A)') &
      'Previous line in namelist: '//trim(line)
    backspace(nin)
    backspace(nin)
    read(nin,fmt='(A)') line
    write(*,'(A)') &
      'Previous line in namelist: '//trim(line)
  endif
  close(nin)
end subroutine read_dirs_shot
!**********************************************************************
!>
!!    read the in1 namelist from an OMAS equilibrium hdf5 file
!!
!!    this only includes necessary information to set up directory paths
!!    
!!    timeslice 0 is used, but could be made variable
!!
!**********************************************************************
subroutine read_omas_in1(filename)
  use var_exdata, only: ishot
  use var_cecoil, only: iecurr
  use var_vessel, only: ivesel
  use var_nio
  use error_control
  use exvars, only: table_dir,input_dir,store_dir,efitversion
  implicit integer*4 (i-n), real*8 (a-h,o-z)
  integer*4  :: istat
  integer*4 :: dims
  character (*) :: filename
  logical :: file_stat
 
#if defined(USE_HDF5)
  inquire(file=trim(filename),exist=file_stat)
  if (.not. file_stat) then
    call errctrl_msg('read_omas_in1',trim(filename)//' not found')
    stop
  endif
  call fch5init
  call open_oldh5file(trim(filename),fileid,rootgid,h5in,h5err)
  call test_group(rootgid,"equilibrium",file_stat,h5err)
  if (.not. file_stat) then
    call errctrl_msg('read_omas_in1','equilibrium group not found')
    stop
  endif
  call open_group(rootgid,"equilibrium",eqid,h5err)
  call test_group(eqid,"code",file_stat,h5err)
  if (.not. file_stat) then
    call errctrl_msg('read_omas_in1','code group not found')
    stop
  endif
  call open_group(eqid,"code",cid,h5err)
  call test_group(cid,"parameters",file_stat,h5err)
  if (.not. file_stat) then
    call errctrl_msg('read_omas_in1','parameters group not found')
    stop
  endif
  call open_group(cid,"parameters",pid,h5err)
  call test_group(pid,"time_slice",file_stat,h5err)
  if (.not. file_stat) then
    call errctrl_msg('read_omas_in1','time_slice group not found')
    stop
  endif
  call open_group(pid,"time_slice",tid,h5err)
  call test_group(tid,"0",file_stat,h5err)
  if (.not. file_stat) then
    call errctrl_msg('read_omas_in1','100 group not found')
    stop
  endif
  call open_group(tid,"0",sid,h5err)
  call test_group(sid,"in1",file_stat,h5err)
  if (file_stat) then
    call open_group(sid,"in1",nid,h5err)
    call read_h5_ex(nid,"ishot",ishot,h5in,h5err)
    call read_h5_ex(nid,"iecurr",iecurr,h5in,h5err)
    call read_h5_ex(nid,"ivesel",ivesel,h5in,h5err)
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
  call errctrl_msg('read_omas_in1','HDF5 needs to be linked')
  stop
#endif
end subroutine read_omas_in1
