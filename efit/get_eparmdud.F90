#include "config.f"
!**********************************************************************
!>
!!    sets default DIII-D values for eparmdud module
!!
!**********************************************************************
subroutine get_eparmdud_defaults()
  use eparm
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

end subroutine


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
end subroutine


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
  use expath, only: table_dir,input_dir,store_dir
  implicit integer*4 (i-n), real*8 (a-h,o-z)
  real*8, dimension(2000):: expmp2,coils,fwtsi,fwtmp2,psibit,bitmpi,denr,denv,fwtfc, &
                            brsp,bitfc,ecurrt,xalpa,xgama,rzeroj,fwtec,bitec, &
                            ppknt,ffknt,wwknt,rbdry,zbdry,ppbdry,kppbdry,pp2bdry,kpp2bdry, &
                            ffbdry,kffbdry,ff2bdry,kff2bdry,wwbdry,kwwbdry,ww2bdry,kww2bdry,&
                            fwtfcsum,fczero,fcsum,fwtbdry,xlim,ylim,rpress,pressr,sigpre,&
                            fwtpre,sibeam,pbeam,dnbeam,dmass,vcurfb,vcurrt
  real*8, dimension(256,256)::calpa,cgama
  integer*4  :: istat
  character(len=1000) :: line, fitzts
  character (*) :: filename
  logical :: fitsiref, fitfcsum
 
    NAMELIST/in1/ishot,itime,plasma,itek,itrace,nxiter,fwtcur,kffcur &
      ,coils,fwtsi,expmp2,fwtmp2,kppcur,mxiter,ierchk,fwtqa,qemp,error &
      ,limitr,xlim,ylim,serror,nbdry,rbdry,zbdry,psibry,nslref,ibunmn &
      ,btor,psibit,bitmpi,bitip,icurrt,icinit,brsp,iweigh,qenp,fwtbp &
      ,relip,zelip,aelip,eelip,qvfit,fwtdlc,betap0,emp,enp,iconvr,icprof &
      ,nextra,ixstrt,scrape,errmin,rbound,npnef,nptef,fwacoil,itimeu &
      ,rcentr,rzero,gammap,cfcoil,fczero,fcsum,islve,icntour,iprobe &
      ,salpha,srm,sbeta,ifref,isumip,n1coil,ifcurr,iecurr,ecurrt,iecoil &
      ,co2cor,vcurrt,dflux,sigdlc,iplim,kinput,limfag,sigprebi,fwtxx &
      ,kprfit,pressr,rpress,zpress,sigpre,npress,tethom,rteth,keqdsk &
      ,zteth,sgteth,npteth,tionex,rion,zion,sigti,nption,dnethom,zeffvs &
      ,rneth,zneth,sgneth,npneth,pbeam,sibeam,nbeam,rzeroj,xalpa,cgama &
      ,ivesel,iexcal,iconsi,fwtfc,xltype,kcalpa,kcgama,calpa,iacoil &
      ,limid,irfila,jzfila,vloop,iqplot,siref,denr,denv,xgama,sgnemin &
      ,nptionf,currn1,ifitvs,bitfc,idfila,relax,saimin,icutfp,acoilc &
      ,sigtii,cutip,iavem,pnbeam,xltype_180,sgtemin,sgprmin,elomin,dnmin &
      ,sgnethi,fcurbd,pcurbd,prbdry,sgtethi,ndokin,zlowimp,kskipvs,limvs &
      ,vcurfb,kpressb,pressbi,prespb,sigppb,kzeroj,rminvs,rmaxvs,errbry &
      ,fwtpre,ibtcomp,klabel,zmaxvs,dnbeam,dmass,nmass,condin,iaveus &
      ,sgtimin,kwripre,kbound,alphafp,kframe,zbound,vsdamp,zminvs,saicon &
      ,kppfnc,kppknt,ppknt,pptens,kfffnc,kffknt,ffknt,fftens,fwtbdry &
      ,kwwfnc,kwwknt,wwknt,wwtens,fwtec,fitsiref,bitec,scalepr,scalesir &
      ,ppbdry,kppbdry,pp2bdry,kpp2bdry,scalea,sigrbd,sigzbd,nbskip &
      ,ffbdry,kffbdry,ff2bdry,kff2bdry,errsil,vbit &
      ,wwbdry,kwwbdry,ww2bdry,kww2bdry,f2edge,fe_width,fe_psin,kedgef &
      ,ktear,kersil,iout,ixray,pedge,kedgep,pe_width,pe_psin &
      ,table_dir,input_dir,store_dir,kautoknt,akchiwt,akerrwt &
      ,kakloop,aktol,kakiter,akgamwt,akprewt &
      ,kpphord,kffhord,keehord,psiecn,dpsiecn,fitzts,isolve,iplcout &
      ,imagsigma,errmag,ksigma,errmagb,brsptu,fitfcsum,fwtfcsum,appendsnap &
      ,idebug,nbdrymx,nsol,rsol,zsol,fwtsol,efitversion,kbetapr,nbdryp,jdebug &
      ,ifindopt,tolbndpsi

  nin=343
  open(unit=nin,status='old',file=filename)
  read (nin,in1,iostat=istat)

  if (istat>0) then
    backspace(nin)
    read(nin,fmt='(A)') line
    write(*,'(A)') &
      'Invalid line in namelist: '//trim(line)
  end if
  close(nin)
end subroutine
!**********************************************************************
!>
!!    read the in1 namelist from an OMAS equilibrium hdf5 file
!!
!!    this includes necessary information to set up directory paths and
!!      much more...
!!    
!!    every scalar is read, but vectors need to have the correct length
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
  use expath, only: table_dir,input_dir,store_dir
  implicit integer*4 (i-n), real*8 (a-h,o-z)
  integer*4  :: istat
  integer :: dims,fitsiref_int,scalea_int,fitfcsum_int,fitdelz_int
  integer :: writepc_int,oldccomp_int,oldcomp_int
  integer :: symmetrize_int,backaverage_int,shape_ext_int,fixpp_int
  character (*) :: filename
  logical :: fitsiref,scalea,fitfcsum
  logical :: file_stat
 
  fitsiref_int=0
  scalea_int=0
  fitfcsum_int=0
  writepc_int=0
  oldccomp_int=0
  oldcomp_int=0
  symmetrize_int=0
  backaverage_int=0
  shape_ext_int=0
  fixpp_int=0

#ifdef HAVE_HDF5
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
    if (obj_exists(nid,"ishot",h5err)) &
      call read_h5(nid,"ishot",ishot,h5in,h5err)
    if (obj_exists(nid,"itime",h5err)) &
      call read_h5(nid,"itime",itime,h5in,h5err)
    if (obj_exists(nid,"iplasma",h5err)) &
      call read_h5(nid,"iplasma",iplasma,h5in,h5err)
    if (obj_exists(nid,"itek",h5err)) &
      call read_h5(nid,"itek",itek,h5in,h5err)
    if (obj_exists(nid,"itrace",h5err)) &
      call read_h5(nid,"itrace",itrace,h5in,h5err)
    if (obj_exists(nid,"nxiter",h5err)) &
      call read_h5(nid,"nxiter",nxiter,h5in,h5err)
    if (obj_exists(nid,"fwtcur",h5err)) &
      call read_h5(nid,"fwtcur",fwtcur,h5in,h5err)
    if (obj_exists(nid,"kffcur",h5err)) &
      call read_h5(nid,"kffcur",kffcur,h5in,h5err)
    if (obj_exists(nid,"kppcur",h5err)) &
      call read_h5(nid,"kppcur",kppcur,h5in,h5err)
    if (obj_exists(nid,"mxiter",h5err)) &
      call read_h5(nid,"mxiter",mxiter,h5in,h5err)
    if (obj_exists(nid,"ierchk",h5err)) &
      call read_h5(nid,"ierchk",ierchk,h5in,h5err)
    if (obj_exists(nid,"fwtqa",h5err)) &
      call read_h5(nid,"fwtqa",fwtqa,h5in,h5err)
    if (obj_exists(nid,"qemp",h5err)) &
      call read_h5(nid,"qemp",qemp,h5in,h5err)
    if (obj_exists(nid,"error",h5err)) &
      call read_h5(nid,"error",error,h5in,h5err)
    if (obj_exists(nid,"limitr",h5err)) &
      call read_h5(nid,"limitr",limitr,h5in,h5err)
    if (obj_exists(nid,"serror",h5err)) &
      call read_h5(nid,"serror",serror,h5in,h5err)
    if (obj_exists(nid,"nbdry",h5err)) &
      call read_h5(nid,"nbdry",nbdry,h5in,h5err)
    if (obj_exists(nid,"psibry",h5err)) &
      call read_h5(nid,"psibry",psibry,h5in,h5err)
    if (obj_exists(nid,"nslref",h5err)) &
      call read_h5(nid,"nslref",nslref,h5in,h5err)
    if (obj_exists(nid,"ibunmn",h5err)) &
      call read_h5(nid,"ibunmn",ibunmn,h5in,h5err)
    if (obj_exists(nid,"btor",h5err)) &
      call read_h5(nid,"btor",btor,h5in,h5err)
    if (obj_exists(nid,"bitip",h5err)) &
      call read_h5(nid,"bitip",bitip,h5in,h5err)
    if (obj_exists(nid,"icurrt",h5err)) &
      call read_h5(nid,"icurrt",icurrt,h5in,h5err)
    if (obj_exists(nid,"icinit",h5err)) &
      call read_h5(nid,"icinit",icinit,h5in,h5err)
    if (obj_exists(nid,"iweigh",h5err)) &
      call read_h5(nid,"iweigh",iweigh,h5in,h5err)
    if (obj_exists(nid,"qenp",h5err)) &
      call read_h5(nid,"qenp",qenp,h5in,h5err)
    if (obj_exists(nid,"fwtbp",h5err)) &
      call read_h5(nid,"fwtbp",fwtbp,h5in,h5err)
    if (obj_exists(nid,"relip",h5err)) &
      call read_h5(nid,"relip",relip,h5in,h5err)
    if (obj_exists(nid,"zelip",h5err)) &
      call read_h5(nid,"zelip",zelip,h5in,h5err)
    if (obj_exists(nid,"aelip",h5err)) &
      call read_h5(nid,"aelip",aelip,h5in,h5err)
    if (obj_exists(nid,"eelip",h5err)) &
      call read_h5(nid,"eelip",eelip,h5in,h5err)
    if (obj_exists(nid,"qvfit",h5err)) &
      call read_h5(nid,"qvfit",qvfit,h5in,h5err)
    if (obj_exists(nid,"fwtdlc",h5err)) &
      call read_h5(nid,"fwtdlc",fwtdlc,h5in,h5err)
    if (obj_exists(nid,"betap0",h5err)) &
      call read_h5(nid,"betap0",betap0,h5in,h5err)
    if (obj_exists(nid,"emp",h5err)) &
      call read_h5(nid,"emp",emp,h5in,h5err)
    if (obj_exists(nid,"enp",h5err)) &
      call read_h5(nid,"enp",enp,h5in,h5err)
    if (obj_exists(nid,"iconvr",h5err)) &
      call read_h5(nid,"iconvr",iconvr,h5in,h5err)
    if (obj_exists(nid,"icprof",h5err)) &
      call read_h5(nid,"icprof",icprof,h5in,h5err)
    if (obj_exists(nid,"nextra",h5err)) &
      call read_h5(nid,"nextra",nextra,h5in,h5err)
    if (obj_exists(nid,"ixstrt",h5err)) &
      call read_h5(nid,"ixstrt",ixstrt,h5in,h5err)
    if (obj_exists(nid,"scrape",h5err)) &
      call read_h5(nid,"scrape",scrape,h5in,h5err)
    if (obj_exists(nid,"errmin",h5err)) &
      call read_h5(nid,"errmin",errmin,h5in,h5err)
    if (obj_exists(nid,"rbound",h5err)) &
      call read_h5(nid,"rbound",rbound,h5in,h5err)
    if (obj_exists(nid,"npnef",h5err)) &
      call read_h5(nid,"npnef",npnef,h5in,h5err)
    if (obj_exists(nid,"nptef",h5err)) &
      call read_h5(nid,"nptef",nptef,h5in,h5err)
    if (obj_exists(nid,"fwacoil",h5err)) &
      call read_h5(nid,"fwacoil",fwacoil,h5in,h5err)
    if (obj_exists(nid,"itimeu",h5err)) &
      call read_h5(nid,"itimeu",itimeu,h5in,h5err)
    if (obj_exists(nid,"rcentr",h5err)) &
      call read_h5(nid,"rcentr",rcentr,h5in,h5err)
    if (obj_exists(nid,"gammap",h5err)) &
      call read_h5(nid,"gammap",gammap,h5in,h5err)
    if (obj_exists(nid,"cfcoil",h5err)) &
      call read_h5(nid,"cfcoil",cfcoil,h5in,h5err)
    if (obj_exists(nid,"islve",h5err)) &
      call read_h5(nid,"islve",islve,h5in,h5err)
    if (obj_exists(nid,"icntour",h5err)) &
      call read_h5(nid,"icntour",icntour,h5in,h5err)
    if (obj_exists(nid,"iprobe",h5err)) &
      call read_h5(nid,"iprobe",iprobe,h5in,h5err)
    if (obj_exists(nid,"salpha",h5err)) &
      call read_h5(nid,"salpha",salpha,h5in,h5err)
    if (obj_exists(nid,"srm",h5err)) &
      call read_h5(nid,"srm",srm,h5in,h5err)
    if (obj_exists(nid,"sbeta",h5err)) &
      call read_h5(nid,"sbeta",sbeta,h5in,h5err)
    if (obj_exists(nid,"ifref",h5err)) &
      call read_h5(nid,"ifref",ifref,h5in,h5err)
    if (obj_exists(nid,"isumip",h5err)) &
      call read_h5(nid,"isumip",isumip,h5in,h5err)
    if (obj_exists(nid,"n1coil",h5err)) &
      call read_h5(nid,"n1coil",n1coil,h5in,h5err)
    if (obj_exists(nid,"ifcurr",h5err)) &
      call read_h5(nid,"ifcurr",ifcurr,h5in,h5err)
    if (obj_exists(nid,"iecurr",h5err)) &
      call read_h5(nid,"iecurr",iecurr,h5in,h5err)
    if (obj_exists(nid,"iecoil",h5err)) &
      call read_h5(nid,"iecoil",iecoil,h5in,h5err)
    if (obj_exists(nid,"co2cor",h5err)) &
      call read_h5(nid,"co2cor",co2cor,h5in,h5err)
    if (obj_exists(nid,"dflux",h5err)) &
      call read_h5(nid,"dflux",dflux,h5in,h5err)
    if (obj_exists(nid,"sigdlc",h5err)) &
      call read_h5(nid,"sigdlc",sigdlc,h5in,h5err)
    if (obj_exists(nid,"iplim",h5err)) &
      call read_h5(nid,"iplim",iplim,h5in,h5err)
    if (obj_exists(nid,"kinput",h5err)) &
      call read_h5(nid,"kinput",kinput,h5in,h5err)
    if (obj_exists(nid,"limfag",h5err)) &
      call read_h5(nid,"limfag",limfag,h5in,h5err)
    if (obj_exists(nid,"sigprebi",h5err)) &
      call read_h5(nid,"sigprebi",sigprebi,h5in,h5err)
    if (obj_exists(nid,"fwtxx",h5err)) &
      call read_h5(nid,"fwtxx",fwtxx,h5in,h5err)
    if (obj_exists(nid,"kprfit",h5err)) &
      call read_h5(nid,"kprfit",kprfit,h5in,h5err)
    if (obj_exists(nid,"zpress",h5err)) &
      call read_h5(nid,"zpress",zpress,h5in,h5err)
    if (obj_exists(nid,"npress",h5err)) &
      call read_h5(nid,"npress",npress,h5in,h5err)
    if (obj_exists(nid,"tethom",h5err)) &
      call read_h5(nid,"tethom",tethom,h5in,h5err)
    if (obj_exists(nid,"rteth",h5err)) &
      call read_h5(nid,"rteth",rteth,h5in,h5err)
    if (obj_exists(nid,"keqdsk",h5err)) &
      call read_h5(nid,"keqdsk",keqdsk,h5in,h5err)
    if (obj_exists(nid,"zteth",h5err)) &
      call read_h5(nid,"zteth",zteth,h5in,h5err)
    if (obj_exists(nid,"sgteth",h5err)) &
      call read_h5(nid,"sgteth",sgteth,h5in,h5err)
    if (obj_exists(nid,"npteth",h5err)) &
      call read_h5(nid,"npteth",npteth,h5in,h5err)
    if (obj_exists(nid,"tionex",h5err)) &
      call read_h5(nid,"tionex",tionex,h5in,h5err)
    if (obj_exists(nid,"rion",h5err)) &
      call read_h5(nid,"rion",rion,h5in,h5err)
    if (obj_exists(nid,"zion",h5err)) &
      call read_h5(nid,"zion",zion,h5in,h5err)
    if (obj_exists(nid,"sigti",h5err)) &
      call read_h5(nid,"sigti",sigti,h5in,h5err)
    if (obj_exists(nid,"nption",h5err)) &
      call read_h5(nid,"nption",nption,h5in,h5err)
    if (obj_exists(nid,"dnethom",h5err)) &
      call read_h5(nid,"dnethom",dnethom,h5in,h5err)
    if (obj_exists(nid,"zeffvs",h5err)) &
      call read_h5(nid,"zeffvs",zeffvs,h5in,h5err)
    if (obj_exists(nid,"rneth",h5err)) &
      call read_h5(nid,"rneth",rneth,h5in,h5err)
    if (obj_exists(nid,"zneth",h5err)) &
      call read_h5(nid,"zneth",zneth,h5in,h5err)
    if (obj_exists(nid,"sgneth",h5err)) &
      call read_h5(nid,"sgneth",sgneth,h5in,h5err)
    if (obj_exists(nid,"npneth",h5err)) &
      call read_h5(nid,"npneth",npneth,h5in,h5err)
    if (obj_exists(nid,"nbeam",h5err)) &
      call read_h5(nid,"nbeam",nbeam,h5in,h5err)
    if (obj_exists(nid,"rzeroj",h5err)) &
      call read_h5(nid,"rzeroj",rzeroj,h5in,h5err)
    if (obj_exists(nid,"ivesel",h5err)) &
      call read_h5(nid,"ivesel",ivesel,h5in,h5err)
    if (obj_exists(nid,"iexcal",h5err)) &
      call read_h5(nid,"iexcal",iexcal,h5in,h5err)
    if (obj_exists(nid,"iconsi",h5err)) &
      call read_h5(nid,"iconsi",iconsi,h5in,h5err)
    if (obj_exists(nid,"xltype",h5err)) &
      call read_h5(nid,"xltype",xltype,h5in,h5err)
    if (obj_exists(nid,"kcalpa",h5err)) &
      call read_h5(nid,"kcalpa",kcalpa,h5in,h5err)
    if (obj_exists(nid,"kcgama",h5err)) &
      call read_h5(nid,"kcgama",kcgama,h5in,h5err)
    if (obj_exists(nid,"iacoil",h5err)) &
      call read_h5(nid,"iacoil",iacoil,h5in,h5err)
    if (obj_exists(nid,"limid",h5err)) &
      call read_h5(nid,"limid",limid,h5in,h5err)
    if (obj_exists(nid,"irfila",h5err)) &
      call read_h5(nid,"irfila",irfila,h5in,h5err)
    if (obj_exists(nid,"jzfila",h5err)) &
      call read_h5(nid,"jzfila",jzfila,h5in,h5err)
    if (obj_exists(nid,"vloop",h5err)) &
      call read_h5(nid,"vloop",vloop,h5in,h5err)
    if (obj_exists(nid,"iqplot",h5err)) &
      call read_h5(nid,"iqplot",iqplot,h5in,h5err)
    if (obj_exists(nid,"siref",h5err)) &
      call read_h5(nid,"siref",siref,h5in,h5err)
    if (obj_exists(nid,"sgnemin",h5err)) &
      call read_h5(nid,"sgnemin",sgnemin,h5in,h5err)
    if (obj_exists(nid,"nptionf",h5err)) &
      call read_h5(nid,"nptionf",nptionf,h5in,h5err)
    if (obj_exists(nid,"currn1",h5err)) &
      call read_h5(nid,"currn1",currn1,h5in,h5err)
    if (obj_exists(nid,"ifitvs",h5err)) &
      call read_h5(nid,"ifitvs",ifitvs,h5in,h5err)
    if (obj_exists(nid,"idfila",h5err)) &
      call read_h5(nid,"idfila",idfila,h5in,h5err)
    if (obj_exists(nid,"relax",h5err)) &
      call read_h5(nid,"relax",relax,h5in,h5err)
    if (obj_exists(nid,"saimin",h5err)) &
      call read_h5(nid,"saimin",saimin,h5in,h5err)
    if (obj_exists(nid,"icutfp",h5err)) &
      call read_h5(nid,"icutfp",icutfp,h5in,h5err)
    if (obj_exists(nid,"acoilc",h5err)) &
      call read_h5(nid,"acoilc",acoilc,h5in,h5err)
    if (obj_exists(nid,"sigtii",h5err)) &
      call read_h5(nid,"sigtii",sigtii,h5in,h5err)
    if (obj_exists(nid,"cutip",h5err)) &
      call read_h5(nid,"cutip",cutip,h5in,h5err)
    if (obj_exists(nid,"iavem",h5err)) &
      call read_h5(nid,"iavem",iavem,h5in,h5err)
    if (obj_exists(nid,"pnbeam",h5err)) &
      call read_h5(nid,"pnbeam",pnbeam,h5in,h5err)
    if (obj_exists(nid,"xltype_180",h5err)) &
      call read_h5(nid,"xltype_180",xltype_180,h5in,h5err)
    if (obj_exists(nid,"sgtemin",h5err)) &
      call read_h5(nid,"sgtemin",sgtemin,h5in,h5err)
    if (obj_exists(nid,"sgprmin",h5err)) &
      call read_h5(nid,"sgprmin",sgprmin,h5in,h5err)
    if (obj_exists(nid,"elomin",h5err)) &
      call read_h5(nid,"elomin",elomin,h5in,h5err)
    if (obj_exists(nid,"dnmin",h5err)) &
      call read_h5(nid,"dnmin",dnmin,h5in,h5err)
    if (obj_exists(nid,"sgnethi",h5err)) &
      call read_h5(nid,"sgnethi",sgnethi,h5in,h5err)
    if (obj_exists(nid,"fcurbd",h5err)) &
      call read_h5(nid,"fcurbd",fcurbd,h5in,h5err)
    if (obj_exists(nid,"pcurbd",h5err)) &
      call read_h5(nid,"pcurbd",pcurbd,h5in,h5err)
    if (obj_exists(nid,"prbdry",h5err)) &
      call read_h5(nid,"prbdry",prbdry,h5in,h5err)
    if (obj_exists(nid,"sgtethi",h5err)) &
      call read_h5(nid,"sgtethi",sgtethi,h5in,h5err)
    if (obj_exists(nid,"ndokin",h5err)) &
      call read_h5(nid,"ndokin",ndokin,h5in,h5err)
    if (obj_exists(nid,"zlowimp",h5err)) &
      call read_h5(nid,"zlowimp",zlowimp,h5in,h5err)
    if (obj_exists(nid,"kskipvs",h5err)) &
      call read_h5(nid,"kskipvs",kskipvs,h5in,h5err)
    if (obj_exists(nid,"limvs",h5err)) &
      call read_h5(nid,"limvs",limvs,h5in,h5err)
    if (obj_exists(nid,"kpressb",h5err)) &
      call read_h5(nid,"kpressb",kpressb,h5in,h5err)
    if (obj_exists(nid,"pressbi",h5err)) &
      call read_h5(nid,"pressbi",pressbi,h5in,h5err)
    if (obj_exists(nid,"prespb",h5err)) &
      call read_h5(nid,"prespb",prespb,h5in,h5err)
    if (obj_exists(nid,"sigppb",h5err)) &
      call read_h5(nid,"sigppb",sigppb,h5in,h5err)
    if (obj_exists(nid,"kzeroj",h5err)) &
      call read_h5(nid,"kzeroj",kzeroj,h5in,h5err)
    if (obj_exists(nid,"rminvs",h5err)) &
      call read_h5(nid,"rminvs",rminvs,h5in,h5err)
    if (obj_exists(nid,"rmaxvs",h5err)) &
      call read_h5(nid,"rmaxvs",rmaxvs,h5in,h5err)
    if (obj_exists(nid,"errbry",h5err)) &
      call read_h5(nid,"errbry",errbry,h5in,h5err)
    if (obj_exists(nid,"ibtcomp",h5err)) &
      call read_h5(nid,"ibtcomp",ibtcomp,h5in,h5err)
    if (obj_exists(nid,"klabel",h5err)) &
      call read_h5(nid,"klabel",klabel,h5in,h5err)
    if (obj_exists(nid,"zmaxvs",h5err)) &
      call read_h5(nid,"zmaxvs",zmaxvs,h5in,h5err)
    if (obj_exists(nid,"nmass",h5err)) &
      call read_h5(nid,"nmass",nmass,h5in,h5err)
    if (obj_exists(nid,"condin",h5err)) &
      call read_h5(nid,"condin",condin,h5in,h5err)
    if (obj_exists(nid,"iaveus",h5err)) &
      call read_h5(nid,"iaveus",iaveus,h5in,h5err)
    if (obj_exists(nid,"sgtimin",h5err)) &
      call read_h5(nid,"sgtimin",sgtimin,h5in,h5err)
    if (obj_exists(nid,"kwripre",h5err)) &
      call read_h5(nid,"kwripre",kwripre,h5in,h5err)
    if (obj_exists(nid,"kbound",h5err)) &
      call read_h5(nid,"kbound",kbound,h5in,h5err)
    if (obj_exists(nid,"alphafp",h5err)) &
      call read_h5(nid,"alphafp",alphafp,h5in,h5err)
    if (obj_exists(nid,"kframe",h5err)) &
      call read_h5(nid,"kframe",kframe,h5in,h5err)
    if (obj_exists(nid,"zbound",h5err)) &
      call read_h5(nid,"zbound",zbound,h5in,h5err)
    if (obj_exists(nid,"vsdamp",h5err)) &
      call read_h5(nid,"vsdamp",vsdamp,h5in,h5err)
    if (obj_exists(nid,"zminvs",h5err)) &
      call read_h5(nid,"zminvs",zminvs,h5in,h5err)
    if (obj_exists(nid,"saicon",h5err)) &
      call read_h5(nid,"saicon",saicon,h5in,h5err)
    if (obj_exists(nid,"kppfnc",h5err)) &
      call read_h5(nid,"kppfnc",kppfnc,h5in,h5err)
    if (obj_exists(nid,"kppknt",h5err)) &
      call read_h5(nid,"kppknt",kppknt,h5in,h5err)
    if (obj_exists(nid,"pptens",h5err)) &
      call read_h5(nid,"pptens",pptens,h5in,h5err)
    if (obj_exists(nid,"kfffnc",h5err)) &
      call read_h5(nid,"kfffnc",kfffnc,h5in,h5err)
    if (obj_exists(nid,"kffknt",h5err)) &
      call read_h5(nid,"kffknt",kffknt,h5in,h5err)
    if (obj_exists(nid,"fftens",h5err)) &
      call read_h5(nid,"fftens",fftens,h5in,h5err)
    if (obj_exists(nid,"kwwfnc",h5err)) &
      call read_h5(nid,"kwwfnc",kwwfnc,h5in,h5err)
    if (obj_exists(nid,"kwwknt",h5err)) &
      call read_h5(nid,"kwwknt",kwwknt,h5in,h5err)
    if (obj_exists(nid,"wwtens",h5err)) &
      call read_h5(nid,"wwtens",wwtens,h5in,h5err)
    if (obj_exists(nid,"fitsiref",h5err)) &
      call read_h5(nid,"fitsiref",fitsiref_int,h5in,h5err)
    if (fitsiref_int.ne.0) fitsiref=.true.  ! default off
    if (obj_exists(nid,"scalepr",h5err)) &
      call read_h5(nid,"scalepr",scalepr,h5in,h5err)
    if (obj_exists(nid,"scalesir",h5err)) &
      call read_h5(nid,"scalesir",scalesir,h5in,h5err)
    if (obj_exists(nid,"scalea",h5err)) &
      call read_h5(nid,"scalea",scalea_int,h5in,h5err)
    if (scalea_int.ne.0) scalea=.true.  ! default off
    if (obj_exists(nid,"sigrbd",h5err)) &
      call read_h5(nid,"sigrbd",sigrbd,h5in,h5err)
    if (obj_exists(nid,"sigzbd",h5err)) &
      call read_h5(nid,"sigzbd",sigzbd,h5in,h5err)
    if (obj_exists(nid,"nbskip",h5err)) &
      call read_h5(nid,"nbskip",nbskip,h5in,h5err)
    if (obj_exists(nid,"errsil",h5err)) &
      call read_h5(nid,"errsil",errsil,h5in,h5err)
    if (obj_exists(nid,"vbit",h5err)) &
      call read_h5(nid,"vbit",vbit,h5in,h5err)
    if (obj_exists(nid,"f2edge",h5err)) &
      call read_h5(nid,"f2edge",f2edge,h5in,h5err)
    if (obj_exists(nid,"fe_width",h5err)) &
      call read_h5(nid,"fe_width",fe_width,h5in,h5err)
    if (obj_exists(nid,"fe_psin",h5err)) &
      call read_h5(nid,"fe_psin",fe_psin,h5in,h5err)
    if (obj_exists(nid,"kedgef",h5err)) &
      call read_h5(nid,"kedgef",kedgef,h5in,h5err)
    if (obj_exists(nid,"ktear",h5err)) &
      call read_h5(nid,"ktear",ktear,h5in,h5err)
    if (obj_exists(nid,"kersil",h5err)) &
      call read_h5(nid,"kersil",kersil,h5in,h5err)
    if (obj_exists(nid,"iout",h5err)) &
      call read_h5(nid,"iout",iout,h5in,h5err)
    if (obj_exists(nid,"ixray",h5err)) &
      call read_h5(nid,"ixray",ixray,h5in,h5err)
    if (obj_exists(nid,"pedge",h5err)) &
      call read_h5(nid,"pedge",pedge,h5in,h5err)
    if (obj_exists(nid,"kedgep",h5err)) &
      call read_h5(nid,"kedgep",kedgep,h5in,h5err)
    if (obj_exists(nid,"pe_width",h5err)) &
      call read_h5(nid,"pe_width",pe_width,h5in,h5err)
    if (obj_exists(nid,"pe_psin",h5err)) &
      call read_h5(nid,"pe_psin",pe_psin,h5in,h5err)
    if (obj_exists(nid,"table_dir",h5err)) &
      call read_h5(nid,"table_dir",table_dir,h5in,h5err)
    if (obj_exists(nid,"input_dir",h5err)) &
      call read_h5(nid,"input_dir",input_dir,h5in,h5err)
    if (obj_exists(nid,"store_dir",h5err)) &
      call read_h5(nid,"store_dir",store_dir,h5in,h5err)
    if (obj_exists(nid,"kautoknt",h5err)) &
      call read_h5(nid,"kautoknt",kautoknt,h5in,h5err)
    if (obj_exists(nid,"akchiwt",h5err)) &
      call read_h5(nid,"akchiwt",akchiwt,h5in,h5err)
    if (obj_exists(nid,"akerrwt",h5err)) &
      call read_h5(nid,"akerrwt",akerrwt,h5in,h5err)
    if (obj_exists(nid,"kakloop",h5err)) &
      call read_h5(nid,"kakloop",kakloop,h5in,h5err)
    if (obj_exists(nid,"aktol",h5err)) &
      call read_h5(nid,"aktol",aktol,h5in,h5err)
    if (obj_exists(nid,"kakiter",h5err)) &
      call read_h5(nid,"kakiter",kakiter,h5in,h5err)
    if (obj_exists(nid,"akgamwt",h5err)) &
      call read_h5(nid,"akgamwt",akgamwt,h5in,h5err)
    if (obj_exists(nid,"akprewt",h5err)) &
      call read_h5(nid,"akprewt",akprewt,h5in,h5err)
    if (obj_exists(nid,"kpphord",h5err)) &
      call read_h5(nid,"kpphord",kpphord,h5in,h5err)
    if (obj_exists(nid,"kffhord",h5err)) &
      call read_h5(nid,"kffhord",kffhord,h5in,h5err)
    if (obj_exists(nid,"keehord",h5err)) &
      call read_h5(nid,"keehord",keehord,h5in,h5err)
    if (obj_exists(nid,"psiecn",h5err)) &
      call read_h5(nid,"psiecn",psiecn,h5in,h5err)
    if (obj_exists(nid,"dpsiecn",h5err)) &
      call read_h5(nid,"dpsiecn",dpsiecn,h5in,h5err)
    if (obj_exists(nid,"isolve",h5err)) &
      call read_h5(nid,"isolve",isolve,h5in,h5err)
    if (obj_exists(nid,"iplcout",h5err)) &
      call read_h5(nid,"iplcout",iplcout,h5in,h5err)
    if (obj_exists(nid,"imagsigma",h5err)) &
      call read_h5(nid,"imagsigma",imagsigma,h5in,h5err)
    if (obj_exists(nid,"errmag",h5err)) &
      call read_h5(nid,"errmag",errmag,h5in,h5err)
    if (obj_exists(nid,"ksigma",h5err)) &
      call read_h5(nid,"ksigma",ksigma,h5in,h5err)
    if (obj_exists(nid,"errmagb",h5err)) &
      call read_h5(nid,"errmagb",errmagb,h5in,h5err)
    if (obj_exists(nid,"brsptu",h5err)) &
      call read_h5(nid,"brsptu",brsptu,h5in,h5err)
    if (obj_exists(nid,"fitfcsum",h5err)) &
      call read_h5(nid,"fitfcsum",fitfcsum_int,h5in,h5err)
    fitfcsum=.false.  ! default off?
    if (obj_exists(nid,"appendsnap",h5err)) &
      call read_h5(nid,"appendsnap",appendsnap,h5in,h5err)
    if (obj_exists(nid,"idebug",h5err)) &
      call read_h5(nid,"idebug",idebug,h5in,h5err)
    if (obj_exists(nid,"nbdrymx",h5err)) &
      call read_h5(nid,"nbdrymx",nbdrymx,h5in,h5err)
    if (obj_exists(nid,"nsol",h5err)) &
      call read_h5(nid,"nsol",nsol,h5in,h5err)
    if (obj_exists(nid,"rsol",h5err)) &
      call read_h5(nid,"rsol",rsol,h5in,h5err)
    if (obj_exists(nid,"zsol",h5err)) &
      call read_h5(nid,"zsol",zsol,h5in,h5err)
    if (obj_exists(nid,"fwtsol",h5err)) &
      call read_h5(nid,"fwtsol",fwtsol,h5in,h5err)
    if (obj_exists(nid,"efitversion",h5err)) &
      call read_h5(nid,"efitversion",efitversion,h5in,h5err)
    if (obj_exists(nid,"kbetapr",h5err)) &
      call read_h5(nid,"kbetapr",kbetapr,h5in,h5err)
    if (obj_exists(nid,"nbdryp",h5err)) &
      call read_h5(nid,"nbdryp",nbdryp,h5in,h5err)
    if (obj_exists(nid,"jdebug",h5err)) &
      call read_h5(nid,"jdebug",jdebug,h5in,h5err)
    if (obj_exists(nid,"ifindopt",h5err)) &
      call read_h5(nid,"ifindopt",ifindopt,h5in,h5err)
    if (obj_exists(nid,"tolbndpsi",h5err)) &
      call read_h5(nid,"tolbndpsi",tolbndpsi,h5in,h5err)
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
end subroutine
