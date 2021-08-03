!**********************************************************************
!>
!!    sets default DIII-D values for eparmdud129 module
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

  nwcur2=nwcurn*2
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
