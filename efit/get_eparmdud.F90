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
  use var_inputc, only: efitversion
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
  use var_inputc, only: efitversion
  use var_nio
  use error_control
  use expath, only: table_dir,input_dir,store_dir
  implicit integer*4 (i-n), real*8 (a-h,o-z)
  integer*4  :: istat
  integer :: dims
  character (*) :: filename
  logical :: fitsiref,fitfcsum
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
    call read_h5_ex(nid,"itime",itime,h5in,h5err)
    call read_h5_ex(nid,"plasma",plasma,h5in,h5err)
    call read_h5_ex(nid,"itek",itek,h5in,h5err)
    call read_h5_ex(nid,"itrace",itrace,h5in,h5err)
    call read_h5_ex(nid,"nxiter",nxiter,h5in,h5err)
    call read_h5_ex(nid,"fwtcur",fwtcur,h5in,h5err)
    call read_h5_ex(nid,"kffcur",kffcur,h5in,h5err)
    call read_h5_ex(nid,"kppcur",kppcur,h5in,h5err)
    call read_h5_ex(nid,"mxiter",mxiter,h5in,h5err)
    call read_h5_ex(nid,"ierchk",ierchk,h5in,h5err)
    call read_h5_ex(nid,"fwtqa",fwtqa,h5in,h5err)
    call read_h5_ex(nid,"qemp",qemp,h5in,h5err)
    call read_h5_ex(nid,"error",error,h5in,h5err)
    call read_h5_ex(nid,"limitr",limitr,h5in,h5err)
    call read_h5_ex(nid,"serror",serror,h5in,h5err)
    call read_h5_ex(nid,"nbdry",nbdry,h5in,h5err)
    call read_h5_ex(nid,"psibry",psibry,h5in,h5err)
    call read_h5_ex(nid,"nslref",nslref,h5in,h5err)
    call read_h5_ex(nid,"ibunmn",ibunmn,h5in,h5err)
    call read_h5_ex(nid,"btor",btor,h5in,h5err)
    call read_h5_ex(nid,"bitip",bitip,h5in,h5err)
    call read_h5_ex(nid,"icurrt",icurrt,h5in,h5err)
    call read_h5_ex(nid,"icinit",icinit,h5in,h5err)
    call read_h5_ex(nid,"iweigh",iweigh,h5in,h5err)
    call read_h5_ex(nid,"qenp",qenp,h5in,h5err)
    call read_h5_ex(nid,"fwtbp",fwtbp,h5in,h5err)
    call read_h5_ex(nid,"relip",relip,h5in,h5err)
    call read_h5_ex(nid,"zelip",zelip,h5in,h5err)
    call read_h5_ex(nid,"aelip",aelip,h5in,h5err)
    call read_h5_ex(nid,"eelip",eelip,h5in,h5err)
    call read_h5_ex(nid,"qvfit",qvfit,h5in,h5err)
    call read_h5_ex(nid,"fwtdlc",fwtdlc,h5in,h5err)
    call read_h5_ex(nid,"betap0",betap0,h5in,h5err)
    call read_h5_ex(nid,"emp",emp,h5in,h5err)
    call read_h5_ex(nid,"enp",enp,h5in,h5err)
    call read_h5_ex(nid,"iconvr",iconvr,h5in,h5err)
    call read_h5_ex(nid,"icprof",icprof,h5in,h5err)
    call read_h5_ex(nid,"nextra",nextra,h5in,h5err)
    call read_h5_ex(nid,"ixstrt",ixstrt,h5in,h5err)
    call read_h5_ex(nid,"scrape",scrape,h5in,h5err)
    call read_h5_ex(nid,"errmin",errmin,h5in,h5err)
    call read_h5_ex(nid,"rbound",rbound,h5in,h5err)
    call read_h5_ex(nid,"npnef",npnef,h5in,h5err)
    call read_h5_ex(nid,"nptef",nptef,h5in,h5err)
    call read_h5_ex(nid,"fwacoil",fwacoil,h5in,h5err)
    call read_h5_ex(nid,"itimeu",itimeu,h5in,h5err)
    call read_h5_ex(nid,"rcentr",rcentr,h5in,h5err)
    call read_h5_ex(nid,"gammap",gammap,h5in,h5err)
    call read_h5_ex(nid,"cfcoil",cfcoil,h5in,h5err)
    call read_h5_ex(nid,"islve",islve,h5in,h5err)
    call read_h5_ex(nid,"icntour",icntour,h5in,h5err)
    call read_h5_ex(nid,"iprobe",iprobe,h5in,h5err)
    call read_h5_ex(nid,"salpha",salpha,h5in,h5err)
    call read_h5_ex(nid,"srm",srm,h5in,h5err)
    call read_h5_ex(nid,"sbeta",sbeta,h5in,h5err)
    call read_h5_ex(nid,"ifref",ifref,h5in,h5err)
    call read_h5_ex(nid,"isumip",isumip,h5in,h5err)
    call read_h5_ex(nid,"n1coil",n1coil,h5in,h5err)
    call read_h5_ex(nid,"ifcurr",ifcurr,h5in,h5err)
    call read_h5_ex(nid,"iecurr",iecurr,h5in,h5err)
    call read_h5_ex(nid,"iecoil",iecoil,h5in,h5err)
    call read_h5_ex(nid,"co2cor",co2cor,h5in,h5err)
    call read_h5_ex(nid,"dflux",dflux,h5in,h5err)
    call read_h5_ex(nid,"sigdlc",sigdlc,h5in,h5err)
    call read_h5_ex(nid,"iplim",iplim,h5in,h5err)
    call read_h5_ex(nid,"kinput",kinput,h5in,h5err)
    call read_h5_ex(nid,"limfag",limfag,h5in,h5err)
    call read_h5_ex(nid,"sigprebi",sigprebi,h5in,h5err)
    call read_h5_ex(nid,"fwtxx",fwtxx,h5in,h5err)
    call read_h5_ex(nid,"kprfit",kprfit,h5in,h5err)
    call read_h5_ex(nid,"zpress",zpress,h5in,h5err)
    call read_h5_ex(nid,"npress",npress,h5in,h5err)
    call read_h5_ex(nid,"tethom",tethom,h5in,h5err)
    call read_h5_ex(nid,"rteth",rteth,h5in,h5err)
    call read_h5_ex(nid,"keqdsk",keqdsk,h5in,h5err)
    call read_h5_ex(nid,"zteth",zteth,h5in,h5err)
    call read_h5_ex(nid,"sgteth",sgteth,h5in,h5err)
    call read_h5_ex(nid,"npteth",npteth,h5in,h5err)
    call read_h5_ex(nid,"tionex",tionex,h5in,h5err)
    call read_h5_ex(nid,"rion",rion,h5in,h5err)
    call read_h5_ex(nid,"zion",zion,h5in,h5err)
    call read_h5_ex(nid,"sigti",sigti,h5in,h5err)
    call read_h5_ex(nid,"nption",nption,h5in,h5err)
    call read_h5_ex(nid,"dnethom",dnethom,h5in,h5err)
    call read_h5_ex(nid,"zeffvs",zeffvs,h5in,h5err)
    call read_h5_ex(nid,"rneth",rneth,h5in,h5err)
    call read_h5_ex(nid,"zneth",zneth,h5in,h5err)
    call read_h5_ex(nid,"sgneth",sgneth,h5in,h5err)
    call read_h5_ex(nid,"npneth",npneth,h5in,h5err)
    call read_h5_ex(nid,"nbeam",nbeam,h5in,h5err)
    call read_h5_ex(nid,"ivesel",ivesel,h5in,h5err)
    call read_h5_ex(nid,"iexcal",iexcal,h5in,h5err)
    call read_h5_ex(nid,"iconsi",iconsi,h5in,h5err)
    call read_h5_ex(nid,"xltype",xltype,h5in,h5err)
    call read_h5_ex(nid,"kcalpa",kcalpa,h5in,h5err)
    call read_h5_ex(nid,"kcgama",kcgama,h5in,h5err)
    call read_h5_ex(nid,"iacoil",iacoil,h5in,h5err)
    call read_h5_ex(nid,"limid",limid,h5in,h5err)
    call read_h5_ex(nid,"irfila",irfila,h5in,h5err)
    call read_h5_ex(nid,"jzfila",jzfila,h5in,h5err)
    call read_h5_ex(nid,"vloop",vloop,h5in,h5err)
    call read_h5_ex(nid,"iqplot",iqplot,h5in,h5err)
    call read_h5_ex(nid,"siref",siref,h5in,h5err)
    call read_h5_ex(nid,"sgnemin",sgnemin,h5in,h5err)
    call read_h5_ex(nid,"nptionf",nptionf,h5in,h5err)
    call read_h5_ex(nid,"currn1",currn1,h5in,h5err)
    call read_h5_ex(nid,"ifitvs",ifitvs,h5in,h5err)
    call read_h5_ex(nid,"idfila",idfila,h5in,h5err)
    call read_h5_ex(nid,"relax",relax,h5in,h5err)
    call read_h5_ex(nid,"saimin",saimin,h5in,h5err)
    call read_h5_ex(nid,"icutfp",icutfp,h5in,h5err)
    call read_h5_ex(nid,"acoilc",acoilc,h5in,h5err)
    call read_h5_ex(nid,"sigtii",sigtii,h5in,h5err)
    call read_h5_ex(nid,"cutip",cutip,h5in,h5err)
    call read_h5_ex(nid,"iavem",iavem,h5in,h5err)
    call read_h5_ex(nid,"pnbeam",pnbeam,h5in,h5err)
    call read_h5_ex(nid,"xltype_180",xltype_180,h5in,h5err)
    call read_h5_ex(nid,"sgtemin",sgtemin,h5in,h5err)
    call read_h5_ex(nid,"sgprmin",sgprmin,h5in,h5err)
    call read_h5_ex(nid,"elomin",elomin,h5in,h5err)
    call read_h5_ex(nid,"dnmin",dnmin,h5in,h5err)
    call read_h5_ex(nid,"sgnethi",sgnethi,h5in,h5err)
    call read_h5_ex(nid,"fcurbd",fcurbd,h5in,h5err)
    call read_h5_ex(nid,"pcurbd",pcurbd,h5in,h5err)
    call read_h5_ex(nid,"prbdry",prbdry,h5in,h5err)
    call read_h5_ex(nid,"sgtethi",sgtethi,h5in,h5err)
    call read_h5_ex(nid,"ndokin",ndokin,h5in,h5err)
    call read_h5_ex(nid,"zlowimp",zlowimp,h5in,h5err)
    call read_h5_ex(nid,"kskipvs",kskipvs,h5in,h5err)
    call read_h5_ex(nid,"limvs",limvs,h5in,h5err)
    call read_h5_ex(nid,"kpressb",kpressb,h5in,h5err)
    call read_h5_ex(nid,"pressbi",pressbi,h5in,h5err)
    call read_h5_ex(nid,"prespb",prespb,h5in,h5err)
    call read_h5_ex(nid,"sigppb",sigppb,h5in,h5err)
    call read_h5_ex(nid,"kzeroj",kzeroj,h5in,h5err)
    call read_h5_ex(nid,"rminvs",rminvs,h5in,h5err)
    call read_h5_ex(nid,"rmaxvs",rmaxvs,h5in,h5err)
    call read_h5_ex(nid,"errbry",errbry,h5in,h5err)
    call read_h5_ex(nid,"ibtcomp",ibtcomp,h5in,h5err)
    call read_h5_ex(nid,"klabel",klabel,h5in,h5err)
    call read_h5_ex(nid,"zmaxvs",zmaxvs,h5in,h5err)
    call read_h5_ex(nid,"nmass",nmass,h5in,h5err)
    call read_h5_ex(nid,"condin",condin,h5in,h5err)
    call read_h5_ex(nid,"iaveus",iaveus,h5in,h5err)
    call read_h5_ex(nid,"sgtimin",sgtimin,h5in,h5err)
    call read_h5_ex(nid,"kwripre",kwripre,h5in,h5err)
    call read_h5_ex(nid,"kbound",kbound,h5in,h5err)
    call read_h5_ex(nid,"alphafp",alphafp,h5in,h5err)
    call read_h5_ex(nid,"kframe",kframe,h5in,h5err)
    call read_h5_ex(nid,"zbound",zbound,h5in,h5err)
    call read_h5_ex(nid,"vsdamp",vsdamp,h5in,h5err)
    call read_h5_ex(nid,"zminvs",zminvs,h5in,h5err)
    call read_h5_ex(nid,"saicon",saicon,h5in,h5err)
    call read_h5_ex(nid,"kppfnc",kppfnc,h5in,h5err)
    call read_h5_ex(nid,"kppknt",kppknt,h5in,h5err)
    call read_h5_ex(nid,"pptens",pptens,h5in,h5err)
    call read_h5_ex(nid,"kfffnc",kfffnc,h5in,h5err)
    call read_h5_ex(nid,"kffknt",kffknt,h5in,h5err)
    call read_h5_ex(nid,"fftens",fftens,h5in,h5err)
    call read_h5_ex(nid,"kwwfnc",kwwfnc,h5in,h5err)
    call read_h5_ex(nid,"kwwknt",kwwknt,h5in,h5err)
    call read_h5_ex(nid,"wwtens",wwtens,h5in,h5err)
    call read_h5_ex(nid,"fitsiref",fitsiref,h5in,h5err)
    call read_h5_ex(nid,"scalepr",scalepr,h5in,h5err)
    call read_h5_ex(nid,"scalesir",scalesir,h5in,h5err)
    call read_h5_ex(nid,"scalea",scalea,h5in,h5err)
    call read_h5_ex(nid,"sigrbd",sigrbd,h5in,h5err)
    call read_h5_ex(nid,"sigzbd",sigzbd,h5in,h5err)
    call read_h5_ex(nid,"nbskip",nbskip,h5in,h5err)
    call read_h5_ex(nid,"errsil",errsil,h5in,h5err)
    call read_h5_ex(nid,"vbit",vbit,h5in,h5err)
    call read_h5_ex(nid,"f2edge",f2edge,h5in,h5err)
    call read_h5_ex(nid,"fe_width",fe_width,h5in,h5err)
    call read_h5_ex(nid,"fe_psin",fe_psin,h5in,h5err)
    call read_h5_ex(nid,"kedgef",kedgef,h5in,h5err)
    call read_h5_ex(nid,"ktear",ktear,h5in,h5err)
    call read_h5_ex(nid,"kersil",kersil,h5in,h5err)
    call read_h5_ex(nid,"iout",iout,h5in,h5err)
    call read_h5_ex(nid,"ixray",ixray,h5in,h5err)
    call read_h5_ex(nid,"pedge",pedge,h5in,h5err)
    call read_h5_ex(nid,"kedgep",kedgep,h5in,h5err)
    call read_h5_ex(nid,"pe_width",pe_width,h5in,h5err)
    call read_h5_ex(nid,"pe_psin",pe_psin,h5in,h5err)
    call read_h5_ex(nid,"table_dir",table_dir,h5in,h5err)
    call read_h5_ex(nid,"input_dir",input_dir,h5in,h5err)
    call read_h5_ex(nid,"store_dir",store_dir,h5in,h5err)
    call read_h5_ex(nid,"kautoknt",kautoknt,h5in,h5err)
    call read_h5_ex(nid,"akchiwt",akchiwt,h5in,h5err)
    call read_h5_ex(nid,"akerrwt",akerrwt,h5in,h5err)
    call read_h5_ex(nid,"kakloop",kakloop,h5in,h5err)
    call read_h5_ex(nid,"aktol",aktol,h5in,h5err)
    call read_h5_ex(nid,"kakiter",kakiter,h5in,h5err)
    call read_h5_ex(nid,"akgamwt",akgamwt,h5in,h5err)
    call read_h5_ex(nid,"akprewt",akprewt,h5in,h5err)
    call read_h5_ex(nid,"kpphord",kpphord,h5in,h5err)
    call read_h5_ex(nid,"kffhord",kffhord,h5in,h5err)
    call read_h5_ex(nid,"keehord",keehord,h5in,h5err)
    call read_h5_ex(nid,"psiecn",psiecn,h5in,h5err)
    call read_h5_ex(nid,"dpsiecn",dpsiecn,h5in,h5err)
    call read_h5_ex(nid,"isolve",isolve,h5in,h5err)
    call read_h5_ex(nid,"iplcout",iplcout,h5in,h5err)
    call read_h5_ex(nid,"imagsigma",imagsigma,h5in,h5err)
    call read_h5_ex(nid,"errmag",errmag,h5in,h5err)
    call read_h5_ex(nid,"ksigma",ksigma,h5in,h5err)
    call read_h5_ex(nid,"errmagb",errmagb,h5in,h5err)
    call read_h5_ex(nid,"brsptu",brsptu,h5in,h5err)
    call read_h5_ex(nid,"fitfcsum",fitfcsum,h5in,h5err)
    call read_h5_ex(nid,"appendsnap",appendsnap,h5in,h5err)
    call read_h5_ex(nid,"idebug",idebug,h5in,h5err)
    call read_h5_ex(nid,"nbdrymx",nbdrymx,h5in,h5err)
    call read_h5_ex(nid,"nsol",nsol,h5in,h5err)
    call read_h5_ex(nid,"rsol",rsol,h5in,h5err)
    call read_h5_ex(nid,"zsol",zsol,h5in,h5err)
    call read_h5_ex(nid,"fwtsol",fwtsol,h5in,h5err)
    ! H5 TODO: this could be a string or a vector so we need a type check...
!    call read_h5_ex(nid,"efitversion",efitversion,h5in,h5err)
    call read_h5_ex(nid,"kbetapr",kbetapr,h5in,h5err)
    call read_h5_ex(nid,"nbdryp",nbdryp,h5in,h5err)
    call read_h5_ex(nid,"jdebug",jdebug,h5in,h5err)
    call read_h5_ex(nid,"ifindopt",ifindopt,h5in,h5err)
    call read_h5_ex(nid,"tolbndpsi",tolbndpsi,h5in,h5err)
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
end subroutine
