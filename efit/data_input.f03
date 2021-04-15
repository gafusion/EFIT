      subroutine data_input(jtime,kconvr,ktime,mtear,kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          data sets up the magnetic data and weighting arrays.    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          29/06/83..........first created                         **
!**          24/07/85..........revised                               **
!**          23/04/04...JAL iplcout added to namelist                **
!**          01/08/07...DPB namelist for mag uncertainty added       **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky,wgridpc,rfcpc
      use set_kinds
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      parameter(mfila=10)
      parameter (m_ext=101)
      
      character(len=1000) :: line

      real*8,allocatable :: gridpf(:,:),gwork(:,:),rgrids(:),zgrids(:)
      real*8,dimension(:),allocatable ::  coils,expmp2,acoilc &
         ,tgamma,sgamma,rrrgam,zzzgam &
         ,aa1gam,aa2gam,aa3gam,aa4gam &
         ,aa5gam,aa6gam,aa7gam,tgammauncor
      real*8,dimension(:),allocatable ::  bmsels,sbmsels,fwtbmsels, &
         rrmsels,zzmsels,l1msels,l2msels, &
         l4msels,emsels,semsels,fwtemsels
      integer*4,dimension(:),allocatable ::  iemsels
      real*8,dimension(:),allocatable ::  tlibim,slibim,rrrlib &
         ,zzzlib,aa1lib,aa8lib,fwtlib
      real*8,dimension(:),allocatable ::  pds,denr,denv
      integer*8,dimension(:),allocatable ::  ilower
      real*8,dimension(:),allocatable ::  devxmpin,rnavxmpin &
               ,devpsiin,rnavpsiin,devfcin,rnavfcin &
               ,devein,rnavecin,brsptu

      namelist/in1/ishot,itime,plasma,itek,itrace,nxiter,fwtcur,kffcur &
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
      namelist/inwant/psiwant,vzeroj,fwtxxj,fbetap,fbetan,fli,fq95,fqsiw &
           ,jbeta,jli,alpax,gamax,jwantm,fwtxxq,fwtxxb,fwtxli,znose &
           ,fwtbdry,nqwant,siwantq,n_write,kccoils,ccoils,rexpan &
           ,xcoils,kcloops,cloops,xloops,currc79,currc139,nccoil,sizeroj &
           ,fitdelz,ndelzon,relaxdz,stabdz,writepc,table_dir,errdelz &
           ,oldccomp,nicoil,oldcomp,currc199,curriu30,curriu90 &
           ,curriu150,curril30,curril90,curril150,ifitdelz,scaledz
      namelist/inms/xmprcg,xmp_k,vresxmp,t0xmp,psircg,psi_k,vrespsi &
           ,t0psi,fcrcg,fc_k,vresfc,t0fc,ercg,e_k,vrese,t0e,bcrcg &
           ,bc_k,vresbc,t0bc,prcg,p_k,vresp,t0p,bti322in,curc79in &
           ,curc139in,curc199in,devxmpin,rnavxmpin,devpsiin,rnavpsiin &
           ,devfcin,rnavfcin,devein,rnavecin
      namelist/ink/isetfb,ioffr,ioffz,ishiftz,gain,gainp,idplace &
           ,symmetrize,backaverage,lring,cupdown
      namelist/ins/tgamma,sgamma,fwtgam,rrrgam,zzzgam,aa1gam,aa2gam, &
            aa3gam,aa4gam,aa5gam,aa6gam,aa7gam,iplots,kdomse, &
            msebkp,msefitfun,mse_quiet,mse_spave_on,kwaitmse, &
            dtmsefull,mse_strict,t_max_beam_off,ok_30rt,ok_210lt, &
            mse_usecer,mse_certree,mse_use_cer330,mse_use_cer210, &
            tgammauncor,v30lt,v30rt,v210lt,v210rt
      namelist/in_msels/bmsels,sbmsels,fwtbmsels,rrmsels,zzmsels, &
            l1msels,l2msels,l4msels,emsels,semsels,fwtemsels,kdomsels, &
            fmlscut,synmsels,avemsels
      namelist/ina/spatial_avg_gam
      namelist/inlibim/tlibim,slibim,fwtlib,rrrlib,zzzlib,aa1lib,aa8lib
      namelist/inece/necein,teecein0,feece0,errorece0,fwtece0,fwtecebz0 &
      ,ecefit,ecebzfit,kfitece,kinputece,kcallece,nharm,kfixro,rteo,zteo &
      ,kfixrece,rtep,rtem,rpbit,rmbit,robit,nfit,kcmin,fwtnow,kdoece &
      ,mtxece,nconstr,eceiter,eceerror
      namelist/iner/keecur,ecurbd,keefnc,eetens,keebdry,kee2bdry, &
                    eebdry,ee2bdry,eeknt,keeknt,keehord
      namelist/insxr/ksxr0,ksxr2,idosxr
      namelist/in3/mpnam2,xmp2,ymp2,amp2,smp2,rsi,zsi,wsi,hsi,as, &
        as2,lpname,rsisvs,vsname,turnfc,patmp2,racoil,zacoil, &
        hacoil,wacoil
      namelist/edgep/symmetrize, &
       rpress , pressr , sigpre, npress , kprfit, kpressb, ndokin, &
       kppfnc,  kfffnc,  kffcur, kppcur,  mxiter, error, errmin, keecur
      namelist/edat/nption,tionex,sigti,rion,zion, &
                    npneth,dnethom,sgneth,rneth,zneth, &
                    npteth,tethom,sgteth,rteth,zteth, &
                    nbrmcrd,bremin,bremsig,brmrtan,brmzelev,ivbcuse, &
                    sigtii,sgnethi,sgtethi,bremsigi, &
                    npress,rpress,zpress,pressr,sigpre
      namelist/invt/omegat,nomegat,enw,emw,betapw0,kvtor,kwwcur,rvtor, &
                    wcurbd,preswb,fwtprw,npresw,presw,sigprw,rpresw, &
                    zpresw,kplotp,sbetaw,nsplot,comega,kcomega,xomega &
                    ,kdovt,romegat,zomegat,sigome,scalepw &
                    ,kwwfnc,kwwknt,wwknt,wwtens
      namelist/efitin/istore,scrape,nextra, &
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
           ppbdry,kppbdry,pp2bdry,kpp2bdry, &
           ffbdry,kffbdry,ff2bdry,kff2bdry,errsil,nicoil, &
           wwbdry,kwwbdry,ww2bdry,kww2bdry,fwtec,fitdelz,fitsiref, &
           nbdry,rbdry,zbdry,sigrbd,sigzbd,nbskip,msebkp, &
           ktear,keecur,ecurbd,keefnc,eetens,eeknt,keeknt, &
           keebdry,kee2bdry,eebdry,ee2bdry,kersil,iout,ixray, &
           use_alternate_pointnames, alternate_pointname_file, &
           do_spline_fit,table_dir,input_dir,store_dir &
          ,kautoknt,akchiwt,akerrwt,akgamwt,akprewt, &
           kakloop,aktol,kakiter,psiecn,dpsiecn,relaxdz,isolve &
          ,iplcout,errdelz,oldcomp,imagsigma,errmag,errmagb &
          ,fitfcsum,fwtfcsum,fixpp &
          ,mse_quiet,mse_spave_on,kwaitmse,dtmsefull &
          ,mse_strict,t_max_beam_off,ifitdelz,scaledz &
          ,mse_usecer,mse_certree,mse_use_cer330,mse_use_cer210 &
          ,ok_30rt,ok_210lt,vbit,nbdrymx,fwtbmsels,fwtemsels,idebug,jdebug &
          ,synmsels,avemsels,kwritime,v30lt,v30rt,v210lt,v210rt &
          ,ifindopt,tolbndpsi
      namelist/profile_ext/npsi_ext,pprime_ext,ffprim_ext,psin_ext, &
      geqdsk_ext,sign_ext,scalepp_ext,scaleffp_ext,shape_ext,dr_ext, &
      dz_ext,rc_ext,zc_ext,a_ext,eup_ext,elow_ext,dup_ext,dlow_ext, &
      setlim_ext,reflect_ext,fixpp
!
      integer :: nw_ext, nh_ext
      real*8 :: c_ext, dr_ext, dz_ext,rc_ext,zc_ext, a_ext
      real*8 :: eup_ext, elow_ext, dup_ext, dlow_ext, setlim_ext
      real*8 :: r0min,r0max,z0min,z0max,zr0min,zr0max,rz0min,rz0max
      real*8 :: r0ave,z0ave,a0ave,e0top,e0bot,d0top,d0bot
      character*10 case_ext(6)
      character*50 edatname
      character*82 table_nam
      character*10 namedum
      character*2 :: reflect_ext
      logical :: shape_ext
      !real*4 spatial_avg_ham(nmtark,ngam_vars,ngam_u,ngam_w)
      data nsq/1/
      data ersil8/1.0e-03_dp/,currn1/0.0/
      data idodo/0/,idovs/0/,zetafc/2.5e-08_dp/
      data co2cor/1.0/,idoac/0/,fq95/0.0/
      data mcontr/35/
      data ten2m3/1.0e-03_dp/
      data idtime/0/,itimeb/0/
      save idodo, idovs, idoac
!
      ALLOCATE(gridpf(nwnh,mfila),gwork(nbwork,nwnh), &
         rgrids(nw),zgrids(nh))
      ALLOCATE(coils(nsilop),expmp2(magpri),acoilc(nacoil) &
         ,tgamma(nmtark),sgamma(nmtark),rrrgam(nmtark),zzzgam(nmtark) &
         ,aa1gam(nmtark),aa2gam(nmtark),aa3gam(nmtark),aa4gam(nmtark) &
         ,aa5gam(nmtark),aa6gam(nmtark),aa7gam(nmtark) &
         ,tgammauncor(nmtark))
      ALLOCATE(bmsels(nmsels),sbmsels(nmsels),fwtbmsels(nmsels), &
         rrmsels(nmsels),zzmsels(nmsels),l1msels(nmsels),l2msels(nmsels), &
         l4msels(nmsels),emsels(nmsels),semsels(nmsels),fwtemsels(nmsels))
      ALLOCATE(iemsels(nmsels))
      ALLOCATE( tlibim(libim),slibim(libim),rrrlib(libim) &
         ,zzzlib(libim),aa1lib(libim),aa8lib(libim),fwtlib(libim))
      ALLOCATE(pds(6),denr(nco2r),denv(nco2v))
      ALLOCATE(ilower(mbdry))
      ALLOCATE(devxmpin(magpri),rnavxmpin(magpri) &
               ,devpsiin(nsilop),rnavpsiin(nsilop) &
               ,devfcin(nfcoil),rnavfcin(nfcoil) &
               ,devein(nesum),rnavecin(nesum),brsptu(nfcoil))

      brsptu(1)=-1.e-20_dp
!
      kerror = 0
      idone=0
      sicont=tmu*drslop/aaslop
!
      if (kdata.eq.2) go to 180
      if (jtime.gt.1) go to 178
!----------------------------------------------------------------------
!-- normalize fitting weights, SNAP mode                             --
!----------------------------------------------------------------------
      do 150 i=1,nsilop
        if ((kdata.ge.3).and.(fwtsi(i).ne.0.0)) then
          fwtsi(i)=swtsi(i)
          if (lookfw.gt.0) fwtsi(i)=rwtsi(i)
        endif
        if (ierpsi(i).ne.0) fwtsi(i)=0.0
  150 continue
      do 152 i=1,nfcoil
        if ((kdata.ge.3).and.(fwtfc(i).ne.0.0)) fwtfc(i)=swtfc(i)
        if (ierfc(i).ne.0) fwtfc(i)=0.0
  152 continue
      if (iecurr.eq.2) then
      do i=1,nesum
        if ((kdata.ge.3).and.(fwtec(i).ne.0.0)) fwtec(i)=swtec(i)
        if (ierec(i).ne.0) fwtec(i)=0.0
      enddo
      endif
      do 160 i=1,magpri
        if ((kdata.ge.3).and.(fwtmp2(i).ne.0.0)) then
          fwtmp2(i)=swtmp2(i)
          if (lookfw.gt.0) fwtmp2(i)=rwtmp2(i)
        endif
        if (iermpi(i).ne.0) fwtmp2(i)=0.0
  160 continue
      do 170 i=1,nstark
        fwtgam(i)=swtgam(i)
        if (iergam(i).ne.0) fwtgam(i)=0.0
  170 continue
      do 172 i=1,nnece
        fwtece0(i)=swtece(i)
        if (ierece(i).ne.0) fwtece0(i)=0.0
  172 continue
        fwtecebz0=swtecebz
        if (ierecebz.ne.0) fwtecebz0=0.0
      if (fwtcur.ne.0.0) fwtcur=swtcur
      if (fwtqa.ne.0.0) fwtqa=1.
      if (fwtbp.ne.0.0) fwtbp=1.
      if (fwtdlc.ne.0.0) fwtdlc=swtdlc
      if (ierpla.ne.0) fwtcur=0.0
      if (ierrdi.ne.0) fwtdlc=0.0
!---------------------------------------------------------------------
!--  Save fitting weights                                           --
!---------------------------------------------------------------------
      swtdlc=fwtdlc
      swtcur=fwtcur
      do i=1,nfcoil
        swtfc(i)=fwtfc(i)
      enddo
      do i=1,nesum
        swtec(i)=fwtec(i)
      enddo
      do i=1,magpri
        swtmp2(i)=fwtmp2(i)
      enddo
      do i=1,nsilop
        swtsi(i)=fwtsi(i)
      enddo
      do i=1,nstark
        swtgam(i)=fwtgam(i)
      enddo
      go to 179
  178 continue
!---------------------------------------------------------------------
!--  Restore fitting weights for time slices > 1                    --
!---------------------------------------------------------------------
      fwtdlc=swtdlc
      fwtcur=swtcur
      if (fwtqa.ne.0.0) fwtqa=1.
      if (fwtbp.ne.0.0) fwtbp=1.
      do i=1,nfcoil
        fwtfc(i)=swtfc(i)
      enddo
      do i=1,nesum
        fwtec(i)=swtec(i)
      enddo
      do i=1,magpri
        fwtmp2(i)=swtmp2(i)
      enddo
      do i=1,nsilop
        fwtsi(i)=swtsi(i)
      enddo
      do i=1,nstark
        fwtgam(i)=swtgam(i)
      enddo
  179 continue
!-----------------------------------------------------------------------
!-- Set edge pedestal tanh paramters                                  --
!-----------------------------------------------------------------------
      if (fitzts.eq.'te'.and.ztserr(jtime)) then
        nbdry=1
        rbdry(1)=1.94_dp
        zbdry(1)=ztssym(jtime)+0.5_dp*ztswid(jtime)
      endif
      go to 280
!----------------------------------------------------------------------
!-- file mode                                                        --
!----------------------------------------------------------------------
  180 continue
      do 182 i=1,nsilop
        psibit(i)=0.0
  182 continue
      do 184 i=1,magpri
        bitmpi(i)=0.0
  184 continue
      alpax(jbeta)=1.e4_dp
      backaverage=.false.
      bitip=0.0
      betap0=0.50_dp
      brsp(1)=-1.e+20_dp
      cfcoil=-1.
      cutip=80000.
      do 188 i=1,nco2v
        denv(i)=0.
  188 continue
      do 190 i=1,nco2r
        denr(i)=0.
  190 continue
      do 194 i=1,nesum
        rsisec(i)=-1.
  194 continue
      emf=1.00
      emp=1.00
      enf=1.00
      enp=1.00
      error=1.0e-03_dp
      fbetap=0.0
      fbetat=0.0
      do 196 i=1,nfcoil
        fcsum(i)=1.0
        fczero(i)=1.0
        fwtfc(i)=0.
        rsisfc(i)=-1.
  196 continue
      do i=1,nesum
        fwtec(i)=0.0
      enddo
      do 197 i=1,mpress
        fwtpre(i)=1.
  197 continue
      fcurbd=1.0
      fli=0.0
      fwtbp=0.0
      fwtdlc=0.0
      do 11199 i=1,nstark
       fwtgam(i)=0.0
11199 continue
      do i=1,nmsels
        fwtbmsels(i)=0.0
        fwtemsels(i)=0.0
      enddo
      do 11299 i=1,nnece
       fwtece0(i)=0.0
11299 continue
      fwtecebz0=0.0
      do 12399 i=1,mbdry
       fwtbdry(i)=1.0
       sigrbd(i)=1.e10_dp
       sigzbd(i)=1.e10_dp
       fwtsol(i)=1.0
12399 continue
      akchiwt=1.0
      akprewt=0.0
      akgamwt=0.0
      akerrwt=0.0
      aktol=0.1_dp
      fwtqa=0.0
      gammap=1.0e+10_dp
      gamax(jli)=-1.e6_dp
      iavem=5
      ibatch=0
      ibound=0
      ibunmn=3
      icinit=2
      icondn=-1
      iconsi=-1
      iconvr=2
      icprof=0
      icurrt=2
      icutfp=0
      iecoil=0
      ierchk=1
      iecurr=1
      iexcal=0
!jal 04/23/2004
      iplcout=0
      ifcurr=0
      ifitvs=0
      ifref=-1
      iplim=0
      iprobe=0
      iqplot=1
      isetfb=0
      idplace=0
      islve=0
      isumip=0
      itek=0
      iout=1                 ! default - write fitout.dat
      itrace=1
      ivacum=0
      ivesel=0
      n1coil=0
      ibtcomp=1
      iweigh=0
      ixray=0
      ixstrt=1
      kautoknt=0
      kakloop=1
      kakiter=25
      kbetapr=0
      keqdsk=1
      kffcur=1
      kinput=0
      kplotpr=1
      kfcurb=0
      kpcurb=0
      kppcur=3
      kpressb=0
      kprfit=0
      kzeroj=0
      limfag=2
      limitr=-33
      mxiter=25
      nbdry=0
      ncstte=1
      ncstne=1
      ncstfp=1
      ncstpp=1
      nextra=1
      nsq=1
      nxiter=1
      pcurbd=1.0
      pnbeam=0.0
      prbdry=0.
      psibry=0.0
      qemp=0.0
      qenp=0.95_dp
      qvfit=0.95_dp
      rzeroj(1)=0.0
      salpha=1._dp/40._dp
      sbeta=1._dp/8._dp
      sbetaw=0.0
      scrape=0.030_dp
      serror=0.03_dp
      sgnemin=0.0
      sgprmin=0.0
      sgtemin=0.0
      sidif=-1.0e+10_dp
      srm=-3.5_dp
      symmetrize=.false.
      xltype=0.0
      xltype_180=0.0
      rmaxis=rzero
      siref=0.
      vcurfb(1)=0.0
      vcurfb(2)= 500.
      vcurfb(3)= 500.
      vloop=0.
      scalepr(1)=-1.
      scalepw(1)=-1.
      isolve=0
      ifindopt = 2
      tolbndpsi = 1.0e-12_dp
!----------------------------------------------------------------------
!--   Read input file for KDATA = 2                                  --
!----------------------------------------------------------------------
      open(unit=nin,status='old',file=ifname(jtime))
!
      xlim(1)=-1.0
      rbdry(1)=-1.0
      itimeu=0
      table_nam = table_dir
      nbdryp=-1
      ktear=0

      read (nin,in1,iostat=istat)

      if (istat/=0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') &
         'Invalid line in namelist: '//trim(line)
      end if

      if (nbdryp==-1) nbdryp=nbdry
      read (nin,ink,err=11111,end=101)
101    continue
11111 close(unit=nin)
      open(unit=nin,status='old',file=ifname(jtime))
      read (nin,ins,err=11113,end=103)
103    continue
11113 close(unit=nin)
      open(unit=nin,status='old',file=ifname(jtime))
      rrmsels(1)=-10.
      read (nin,in_msels,err=91113,end=91103)
91103    continue
         if (kdomsels.gt.0) then
!----------------------------------------------------------------------
!--   Read msels_all.dat if needed                                   --
!----------------------------------------------------------------------
     if (rrmsels(1).lt.-5.0) then
             close(unit=nin)
             open(unit=nin,status='old',file='msels_all.dat')
             do i=1,nmsels
         read (nin,91008,err=91110) bmsels(i),sbmsels(i),rrmsels(i), &
                    zzmsels(i), l1msels(i),l2msels(i),l4msels(i),emsels(i), &
                    semsels(i),iemsels(i)
               iermselt(jtime,i)=iemsels(i)
             enddo
             go to 91113
91008        format(9e12.5,i2)
91110        continue
             im1=i-1
             write (6,*) 'msels_all i=',i,' bmsels(i-1)= ',bmsels(im1)
             stop
           endif
         endif
91113 close(unit=nin)
      open(unit=nin,status='old',file=ifname(jtime) &
                                 )
      read (nin,ina,err=11331,end=1031)
1031    continue
11331 close(unit=nin)
      open(unit=nin,status='old',file=ifname(jtime) &
                                 )
      ecefit = 0.0
      ecebzfit = 0.0
      read (nin,inece,err=11112,end=102)
102    continue
11112 close(unit=nin)
      open(unit=nin,status='old',file=ifname(jtime) &
                                 )
      read (nin,edgep,err=11117,end=1021)
1021  continue
11117 close(unit=nin)
      open(unit=nin,status='old',file=ifname(jtime) &
                                 )
      read (nin,iner,err=11193,end=12193)
12193  continue
11193 close(unit=nin)
      open(unit=nin,status='old',file=ifname(jtime) &
                                 )
      read (nin,insxr,err=11114,end=104)
104    continue
11114 close(unit=nin)
      open(unit=nin,status='old',file=ifname(jtime) &
                                 )
      read (nin,inms,err=11115,end=105)
105    continue
11115 close(unit=nin)
      bti322(jtime)=bti322in
      curc79(jtime)=curc79in
      curc139(jtime)=curc139in
      curc199(jtime)=curc199in
      devxmp(jtime,:)=devxmpin(:)
      rnavxmp(jtime,:)=rnavxmpin(:)
      devpsi(jtime,:)=devpsiin(:)
      rnavpsi(jtime,:)=rnavpsiin(:)
      devfc(jtime,:)=devfcin(:)
      rnavfc(jtime,:)=rnavfcin(:)
      deve(jtime,:)=devein(:)
      rnavec(jtime,:)=rnavecin(:)
      devbc(jtime)=devbcin
      rnavbc(jtime)=rnavbcin
      devp(jtime)=devpin
      rnavp(jtime)=rnavpin
!
      open(unit=nin,status='old',file=ifname(jtime))
      read (nin,inwant,err=11219,end=106)
106    continue
11219 close(unit=nin)
      open(unit=nin,status='old',file=ifname(jtime))
      read (nin,invt,err=11230,end=11221)
11221 continue
11230 close(unit=nin)
!----------------------------------------------------------------------
!--   Input FF', P' arrays                                           --
!----------------------------------------------------------------------
      geqdsk_ext = 'none'
      psin_ext(1) = -1000.0
      sign_ext = 1.0
      cratio_ext = 1.0
      cratiop_ext = 1.0
      cratiof_ext = 1.0
      scalepp_ext=1.0
      scaleffp_ext=1.0
      dr_ext=0.0
      dz_ext=0.0
      shape_ext=.false.
      rc_ext=-10.
      zc_ext=-10.
      a_ext=-10.
      eup_ext=-10.
      elow_ext=-10.
      dup_ext=-10.
      dlow_ext=-10.
      setlim_ext=-10.
      reflect_ext='no'
!
      open(unit=nin,status='old',file=ifname(jtime))
      if (idebug /= 0) write (nttyo,*) 'DATA_INPUT geqdsk_ext=',geqdsk_ext
      read(nin,profile_ext,err=11777,end=777)
777   continue
      if (idebug /= 0) write (nttyo,*) 'DATA_INPUT geqdsk_ext=',geqdsk_ext
      if (geqdsk_ext.ne.'none') then
        open(unit=neqdsk,status='old',file=geqdsk_ext)
        read (neqdsk,11775) (case_ext(i),i=1,6),nh_ext,nw_ext,nh_ext
        npsi_ext=nw_ext
        if (idebug /= 0) write (nttyo,*) 'npsi_ext,nw_ext=',npsi_ext, &
           nw_ext
        do i = 1,2
          read (neqdsk,11773)
        enddo
        read (neqdsk,11776) plasma_ext,c_ext,c_ext,c_ext,c_ext
        if (plasma_ext > 0.0) then
          sign_ext = -1.0
        endif
        do i = 1,1
          read (neqdsk,11773)
        enddo
        read (neqdsk,11776) (ffprim_ext(i),i=1,nw_ext)
        read (neqdsk,11776) (pprime_ext(i),i=1,nw_ext)
        prbdry=pprime_ext(nw_ext)
        read (neqdsk,11776) (ffprim_ext(i),i=1,nw_ext)
        read (neqdsk,11776) (pprime_ext(i),i=1,nw_ext)
        read (neqdsk,11776) ((psirz_ext,i=1,nw_ext),j=1,nh_ext)
        read (neqdsk,11776,err=11777) (qpsi_ext(i),i=1,nw_ext)
        read (neqdsk,11774,err=11777) nbdry_ext,limitr_ext
        read (neqdsk,11776,err=11777) (rbdry_ext(i),zbdry_ext(i),i=1,nbdry_ext)
        read (neqdsk,11776,err=11777) (xlim_ext(i),ylim_ext(i),i=1,limitr_ext)
11773   format (a)
11774   format (2i5)
11775   format (6a8,3i4)
11776   format (5e16.9)
!
        if (nbdry.le.0) then
          nbabs_ext=nbdry_ext/mbdry1+1
          jb_ext=0
          do i=1,nbdry_ext,nbabs_ext
            jb_ext=jb_ext+1
            rbdry(jb_ext)=rbdry_ext(i)
            zbdry(jb_ext)=zbdry_ext(i)
          enddo
          nbdry=jb_ext
!
          rbmin_e=rbdry(1)
          rbmax_e=rbdry(1)
          zbmin_e=zbdry(1)
          zbmax_e=zbdry(1)
          nrmin_e=1
          nrmax_e=1
          nzmin_e=1
          nzmax_e=1
          do i=1,nbdry
            fwtbdry(i)=1.
            if (rbdry(i).lt.rbmin_e) then
              rbmin_e=rbdry(i)
              nrmin_e=i
            endif
            if (rbdry(i).gt.rbmax_e) then
              rbmax_e=rbdry(i)
              nrmax_e=i
            endif
            if (zbdry(i).lt.zbmin_e) then
              zbmin_e=zbdry(i)
              nzmin_e=i
            endif
            if (zbdry(i).gt.zbmax_e) then
              zbmax_e=zbdry(i)
              nzmax_e=i
            endif
          enddo
          fwtbdry(nrmin_e)=10.
          fwtbdry(nrmax_e)=10.
          fwtbdry(nzmin_e)=10.
          fwtbdry(nzmax_e)=10.
        endif
!
        if (limitr.le.0) then
          do i=1,limitr_ext
            xlim(i)=xlim_ext(i)
            ylim(i)=ylim_ext(i)
          enddo
        endif
      endif
11777 close(nin)
      if (idebug /= 0) write (nttyo,*) 'npsi_ext=',npsi_ext
      if (npsi_ext > 0) then
        if (idebug /= 0) then
          write (nttyo,*) 'scalepp_ext,pprime_ext= ',scalepp_ext, &
             pprime_ext(1)
          write (nttyo,*) 'scaleffp_ext,ffpprim_ext= ',scaleffp_ext, &
             ffprim_ext(1)
        endif
        pprime_ext = pprime_ext*darea*sign_ext*scalepp_ext
        ffprim_ext = ffprim_ext*darea/twopi/tmu*sign_ext*scaleffp_ext
        prbdry=prbdry*scalepp_ext*scalepp_ext
        if (idebug /= 0) then
          write (nttyo,*) 'scalepp_ext,pprime_ext= ',scalepp_ext, &
             pprime_ext(1)
          write (nttyo,*) 'scaleffp_ext,ffpprim_ext= ',scaleffp_ext, &
             ffprim_ext(1)
        endif

        if (psin_ext(1) < 0) then
          do i = 1, npsi_ext
            psin_ext(i) = real(i-1,dp)/real(npsi_ext-1,dp)
          enddo
        endif
        call zpline(npsi_ext,psin_ext,pprime_ext,bpp_ext,cpp_ext,dpp_ext)
        call zpline(npsi_ext,psin_ext,ffprim_ext,bfp_ext,cfp_ext,dfp_ext)
      endif
!----------------------------------------------------------------------
!-- Scale boundary points                                            --
!----------------------------------------------------------------------
      rbdry(1:nbdry)=rbdry(1:nbdry)+dr_ext
      zbdry(1:nbdry)=zbdry(1:nbdry)+dz_ext
      if (shape_ext) then
        rbdry0(1:nbdry)=rbdry(1:nbdry)
        zbdry0(1:nbdry)=zbdry(1:nbdry)
        r0min=rbdry0(1)
        r0max=r0min
        z0min=zbdry0(1)
        z0max=z0min
        do i=1,nbdry
          if (rbdry0(i).le.r0min) then
            r0min=rbdry0(i)
            zr0min=zbdry0(i)
          endif
          if (rbdry0(i).ge.r0max) then
            r0max=rbdry0(i)
            zr0max=zbdry0(i)
          endif
          if (zbdry0(i).le.z0min) then
            z0min=zbdry0(i)
            rz0min=rbdry0(i)
          endif
          if (zbdry0(i).ge.z0max) then
            z0max=zbdry0(i)
            rz0max=rbdry0(i)
          endif
        enddo
        r0ave=0.5*(r0min+r0max)
        z0ave=0.5*(zr0min+zr0max)
        a0ave=0.5*(r0max-r0min)
        e0top=(z0max-z0ave)/a0ave
        e0bot=(z0ave-z0min)/a0ave
        d0top=(r0ave-rz0max)/a0ave
        d0bot=(r0ave-rz0min)/a0ave
        if (rc_ext.le.-10.0) rc_ext=r0ave
        if (zc_ext.le.-10.0) zc_ext=z0ave
        if (a_ext.le.-10.0) a_ext=a0ave
        if (eup_ext.le.-10.0) eup_ext=e0top
        if (elow_ext.le.-10.0) elow_ext=e0bot
        if (dup_ext.le.-10.0) dup_ext=d0top
        if (dlow_ext.le.-10.0) dlow_ext=d0bot
        do i=1,nbdry
          if (zbdry0(i).gt.z0ave) then
            rbdry(i)=rc_ext+a_ext*(rbdry0(i)-r0ave)/a0ave               &
                 +a_ext*(d0top-dup_ext)*((zbdry0(i)-z0ave)/e0top/a0ave)**2
            zbdry(i)=zc_ext+eup_ext*a_ext*(zbdry0(i)-z0ave)/a0ave/e0top
          endif
          if (zbdry0(i).le.z0ave) then
            rbdry(i)=rc_ext+a_ext*(rbdry0(i)-r0ave)/a0ave               &
                +a_ext*(d0bot-dlow_ext)*((z0ave-zbdry0(i))/e0bot/a0ave)**2
            zbdry(i)=zc_ext+elow_ext*a_ext*(zbdry0(i)-z0ave)/a0ave/e0bot
          endif
        enddo
      endif
      if (reflect_ext.eq.'UL') then
        rbdry0(1:nbdry)=rbdry(1:nbdry)
        zbdry0(1:nbdry)=zbdry(1:nbdry)
        zbdry(1:nbdry)=-zbdry0(1:nbdry)
      endif
!----------------------------------------------------------------------
!-- Reflection, Lower = -Upper                                       --
!----------------------------------------------------------------------
      if (reflect_ext.eq.'UU') then
        rbdry0(1:nbdry)=rbdry(1:nbdry)
        zbdry0(1:nbdry)=zbdry(1:nbdry)
        nbdry0=nbdry
        if (zbdry0(1).le.0.0) then
        nbdry=0
        do i=1,nbdry0
          if (zbdry0(i).ge.0.0) then
            nbdry=nbdry+1
            rbdry(nbdry)=rbdry0(i)
            zbdry(nbdry)=zbdry0(i)
          endif
        enddo
        nbdry0=nbdry
        do i=1,nbdry0
          j=nbdry0-i+1
          if (zbdry(j).gt.0.0) then
            nbdry=nbdry+1
            rbdry(nbdry)=rbdry(j)
            zbdry(nbdry)=-zbdry(j)
          endif
        enddo
        endif
!
        if (zbdry0(1).gt.0.0) then
        nbdry=0
        nbranch=0
        do i=1,nbdry0
          if ((zbdry0(i).ge.0.0).and.(nbranch.eq.0)) then
            nbdry=nbdry+1
            rbdry(nbdry)=rbdry0(i)
            zbdry(nbdry)=zbdry0(i)
          else
            nbranch=1
          endif
        enddo
        nbdry1=nbdry
        do i=1,nbdry1
          if (zbdry(i).gt.0.0) then
            nbdry=nbdry+1
            rbdry(nbdry)=rbdry(i)
            zbdry(nbdry)=-zbdry(i)
          endif
        enddo
        nbdry2=nbdry
        nbranch=0
        do i=1,nbdry0
          j=nbdry0-i+1
          if ((zbdry0(j).ge.0.0).and.(nbranch.eq.0)) then
            nbdry=nbdry+1
            rbdry(nbdry)=rbdry0(j)
            zbdry(nbdry)=-zbdry0(j)
          else
            nbranch=1
          endif
        enddo
        nbdry3=nbdry
        do i=1,nbdry3-nbdry2
          j=nbdry3-i+1
          if (zbdry(j).lt.0.0) then
            nbdry=nbdry+1
            rbdry(nbdry)=rbdry(j)
            zbdry(nbdry)=-zbdry(j)
          endif
        enddo
        endif
      endif
!----------------------------------------------------------------------
!-- Reflection, Upper = - Lower                                      --
!----------------------------------------------------------------------
      if (reflect_ext.eq.'LL') then
        rbdry0(1:nbdry)=rbdry(1:nbdry)
        zbdry0(1:nbdry)=zbdry(1:nbdry)
        nbdry0=nbdry
        if (zbdry0(1).ge.0.0) then
        nbdry=0
        do i=1,nbdry0
          if (zbdry0(i).le.0.0) then
            nbdry=nbdry+1
            rbdry(nbdry)=rbdry0(i)
            zbdry(nbdry)=zbdry0(i)
          endif
        enddo
        nbdry0=nbdry
        do i=1,nbdry0
          j=nbdry0-i+1
          if (zbdry(j).lt.0.0) then
            nbdry=nbdry+1
            rbdry(nbdry)=rbdry(j)
            zbdry(nbdry)=-zbdry(j)
          endif
        enddo
        endif
!
        if (zbdry0(1).lt.0.0) then
        nbdry=0
        nbranch=0
        do i=1,nbdry0
          if ((zbdry0(i).le.0.0).and.(nbranch.eq.0)) then
            nbdry=nbdry+1
            rbdry(nbdry)=rbdry0(i)
            zbdry(nbdry)=zbdry0(i)
          else
            nbranch=1
          endif
        enddo
        nbdry1=nbdry
        do i=1,nbdry1
          if (zbdry(i).lt.0.0) then
            nbdry=nbdry+1
            rbdry(nbdry)=rbdry(i)
            zbdry(nbdry)=-zbdry(i)
          endif
        enddo
        nbdry2=nbdry
        nbranch=0
        do i=1,nbdry0
          j=nbdry0-i+1
          if ((zbdry0(j).le.0.0).and.(nbranch.eq.0)) then
            nbdry=nbdry+1
            rbdry(nbdry)=rbdry0(j)
            zbdry(nbdry)=-zbdry0(j)
          else
            nbranch=1
          endif
        enddo
        nbdry3=nbdry
        do i=1,nbdry3-nbdry2
          j=nbdry3-i+1
          if (zbdry(j).gt.0.0) then
            nbdry=nbdry+1
            rbdry(nbdry)=rbdry(j)
            zbdry(nbdry)=-zbdry(j)
          endif
        enddo
        endif
      endif
!----------------------------------------------------------------------
!-- Scale limiter                                                    --
!----------------------------------------------------------------------
      if (setlim_ext.gt.-10.0) then
      rbdry0(1:nbdry)=rbdry(1:nbdry)
      zbdry0(1:nbdry)=zbdry(1:nbdry)
      r0min=rbdry0(1)
      r0max=r0min
      z0min=zbdry0(1)
      z0max=z0min
      do i=1,nbdry
        if (rbdry0(i).le.r0min) then
           r0min=rbdry0(i)
           zr0min=zbdry0(i)
        endif
        if (rbdry0(i).ge.r0max) then
           r0max=rbdry0(i)
           zr0max=zbdry0(i)
        endif
        if (zbdry0(i).le.z0min) then
          z0min=zbdry0(i)
          rz0min=rbdry0(i)
        endif
        if (zbdry0(i).ge.z0max) then
          z0max=zbdry0(i)
          rz0max=rbdry0(i)
        endif
      enddo
        limitr=7
        xlim(1)=r0min-setlim_ext*(r0max-r0min)
        ylim(1)=0.5_dp*(z0min+z0max)
        xlim(2)=xlim(1)
        ylim(2)=z0max+setlim_ext*(z0max-z0min)
        xlim(3)=r0max+setlim_ext*(r0max-r0min)
        ylim(3)=ylim(2)
        xlim(4)=xlim(3)
        ylim(4)=ylim(1)
        xlim(5)=xlim(4)
        ylim(5)=z0min-setlim_ext*(z0max-z0min)
        xlim(6)=xlim(1)
        ylim(6)=ylim(5)
        xlim(7)=xlim(1)
        ylim(7)=ylim(1)
      endif
!----------------------------------------------------------------------
!--   Read Li beam data                                              --
!----------------------------------------------------------------------
      do i=1,libim
        fwtlib(i)=0.0
        rrrlib(i)=0.0
      enddo
      open(unit=nin,status='old',file=ifname(jtime))
      read (nin,inlibim,err=11237,end=11233)
11233 continue
11237 close(unit=nin)
!----------------------------------------------------------------------
!--   recalculate length of default directories in case any change   --
!----------------------------------------------------------------------
      call set_table_dir
      call efit_read_tables

11337 continue
!---------------------------------------------------------------------
!--  specific choice of current profile                             --
!--       ICPROF=1  no edge current density allowed                 --
!--       ICPROF=2  free edge current density                       --
!--       ICPROF=3  weak edge current density constraint            --
!---------------------------------------------------------------------
      if (icprof.eq.1) then
        kffcur=2
        kppcur=2
        fcurbd=1.
        pcurbd=1.
        fwtbp=1.
        fwtqa=0.
        qvfit=0.
      elseif (icprof.eq.2) then
        kffcur=2
        kppcur=2
        fcurbd=0.
        pcurbd=0.
        fwtbp=0.
        fwtqa=0.
        qvfit=0.
      elseif (icprof.eq.3) then
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
      endif
      if(mse_usecer .eq. 1)keecur = 0
      if(mse_usecer .eq. 2 .and. keecur .eq. 0) then
           keecur = 2
           keefnc = 0
           itek = 5
      endif
      mtear=ktear
      if (kedgep.gt.0) then
        s1edge=(1.0-pe_psin)/pe_width
        tpedge=tanh(s1edge)
        s1edge=(1.0-fe_psin)/fe_width
        tfedge=tanh(s1edge)
      endif
      if (imagsigma.gt.0) then
         do_spline_fit=.false.
         saimin=300.
      endif
      if (kzeroj.eq.1.and.sizeroj(1).lt.0.0) sizeroj(1)=psiwant
      print *, 'before save fitting weights'
!---------------------------------------------------------------------
!--  save fitting weights for FILE mode                             --
!---------------------------------------------------------------------
      swtdlc=fwtdlc
      swtcur=fwtcur
      do i=1,nfcoil
        swtfc(i)=fwtfc(i)
      enddo
      do i=1,nesum
        swtec(i)=fwtec(i)
      enddo
      do i=1,magpri
        swtmp2(i)=fwtmp2(i)
      enddo
      do i=1,nsilop
        swtsi(i)=fwtsi(i)
      enddo
      do i=1,nmtark
        swtgam(i)=fwtgam(i)
      enddo
      do i=1,nmsels
        swtbmsels(i)=fwtbmsels(i)
        swtemsels(i)=fwtemsels(i)
      enddo
      do i=1,nnece
        swtece(i)=fwtece0(i)
      enddo
      swtecebz=fwtecebz0
      print *, 'adjust fit parameters based on basis function selected '
!-----------------------------------------------------------------------
!-- adjust fit parameters based on basis function selected            --
!-----------------------------------------------------------------------
       if (kppfnc .eq. 3) then
          kppcur = 4 * (kppknt - 1)
       endif
       if (kppfnc .eq. 4) then
          kppcur = 4 * (kppknt - 1)
       endif
       if (kppfnc .eq. 5) then
          kppcur = kppcur * (kppknt - 1)
       endif
       if (kppfnc .eq. 6) then
          kppcur = kppknt * 2
       endif
       if (kfffnc .eq. 3) then
          kffcur = 4 * (kffknt - 1)
       endif
       if (kfffnc .eq. 4) then
          kffcur = 4 * (kffknt - 1)
       endif
       if (kfffnc .eq. 5) then
          kffcur = kffcur * (kffknt - 1)
       endif
       if (kfffnc .eq. 6) then
          kffcur = kffknt * 2
       endif
       if (kwwfnc .eq. 3) then
          kwwcur = 4 * (kwwknt - 1)
       endif
       if (kwwfnc .eq. 4) then
          kwwcur = 4 * (kwwknt - 1)
       endif
       if (kwwfnc .eq. 5) then
          kwwcur = kwwcur * (kwwknt - 1)
       endif
       if (kwwfnc .eq. 6) then
          kwwcur = kwwknt * 2
       endif
       if (keecur.gt.0) then
       if (keefnc .eq. 3) then
          keecur = 4 * (keeknt - 1)
       endif
       if (keefnc .eq. 4) then
          keecur = 4 * (keeknt - 1)
       endif
       if (keefnc .eq. 5) then
          keecur = keecur * (keeknt - 1)
       endif
       if (keefnc .eq. 6) then
          keecur = keeknt * 2
       endif
       endif
!
      if (fbetan.gt.0.0) brsp(nfcoil+jbeta)=alpax(jbeta)*darea
      if (fli.gt.0.0) brsp(nfcoil+kppcur+jli)=gamax(jli)*darea
      if (kedgep.gt.0) pedge=pedge*darea
      if (kedgef.gt.0) f2edge=f2edge*darea
      if (npress.lt.0) then
        kdopre=-npress
        npress=0
      endif
      if (itek.lt.0) then
        ixray=1
        itek=iabs(itek)
      endif
      if (psiwant.le.0.0) psiwant=1.e-5_dp

      print *, 'itek > 100, write out PLTOUT.OUT individually '
!--------------------------------------------------------------------------
!-- itek > 100, write out PLTOUT.OUT individually                        --
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
        if (fitdelz) itell=4
      endif
      if (nxiter.lt.0) then
        nxiter=-nxiter
        itell=1
        if (fbetan.gt.0.0) itell=2
        if (fli.gt.0.0)    itell=2
        if (nqwant.gt.0.0)  itell=2
        if ((symmetrize).and.(nbdry.gt.1)) itell=3
      endif
      if ((iconvr.ne.3).and.(qvfit.gt.0.0)) qenp=qvfit
      if (ishot.eq.-1) then
          kerror = 1
          call errctrl_msg('data_input','shot number not set')
          return
      end if
      if ((limitr.gt.0).and.(xlim(1).le.-1.0)) read (nin,5000) &
           (xlim(i),ylim(i),i=1,limitr)
      if ((nbdry.gt.0).and.(rbdry(1).le.-1.0)) read (nin,5020) &
          (rbdry(i),zbdry(i),i=1,nbdry)
      if (kprfit.gt.0) then
        do i=1,npress
          premea(i)=pressr(i)
        enddo
      endif
      if (kprfit.gt.0.and.sigpre(1).lt.0.0) then
        scalep=abs(sigpre(1))
        scalemin=abs(sigpre(2))
        do 40010 i=1,npress
          sigpre(i)=scalep*pressr(i)
          sigpre(i)=max(sigpre(i),scalemin)
40010   continue
      endif
      if (kprfit.gt.0.and.scalepr(1).gt.0.0) then
        do i=1,npress
          pressr(i)=pressr(i)*scalepr(i)
          sigpre(i)=sigpre(i)*scalepr(i)
        enddo
      endif
      if (kprfit.gt.0.and.scalepw(1).gt.0.0) then
        if (npresw.gt.0) then
          do i=1,npresw
            presw(i)=presw(i)*scalepw(i)
            sigprw(i)=sigprw(i)*scalepw(i)
          enddo
        elseif (nomegat.gt.0) then
          do i=1,nomegat
            omegat(i)=omegat(i)*scalepw(i)
            sigome(i)=sigome(i)*scalepw(i)
          enddo
        endif
      endif
      print *, 'option to symmetrize added 8/14/91 eal  '
!--------------------------------------------------------------
!-- option to symmetrize added 8/14/91 eal                   --
!--------------------------------------------------------------
      if (symmetrize) then  ! symmetrize the fixed boundery
        isetfb=0  ! be sure vertical feedback is off
        zelip=0 ! must be symmetric about midplane
        if(nbdry.gt.1)then  ! remove duplicate point
          nbryup=0
          delx2=(rbdry(1)-rbdry(nbdry))**2
          dely2=(zbdry(1)-zbdry(nbdry))**2
          if ((delx2+dely2).le.1.0e-08_dp) nbdry=nbdry-1
          rbar=0
          do i=1,nbdry
            rbar=rbar+rbdry(i)/nbdry
          enddo
          do i=1,nbdry
            ilower(i)=-1
            if(zbdry(i).gt.0.)then
              aup=atan2(zbdry(i),rbdry(i)-rbar)
              iup=i
              close=1e30
              do j=1,nbdry
                if(zbdry(j).lt.0.)then
                  adn=atan2(zbdry(j),rbdry(j)-rbar)
                  val=abs(aup+adn)
                  if(val.lt.close)then
                    close=val
                    idn=j
                  endif
                endif
              enddo
              rnow=(rbdry(iup)+rbdry(idn))/2
              znow=(zbdry(iup)-zbdry(idn))/2
              rbdry(iup)=rnow
              rbdry(idn)=rnow
              zbdry(iup)=znow
              zbdry(idn)=-znow
              nbryup=nbryup+1
              ilower(iup)=idn
            endif
          enddo
        endif
      endif ! end boundary symmertization
!---------------------------------------------------------------------
!--  Symmetrize the limiter positions for fixed boundary if request --
!--  set LIMITR=1000+points for this option                         --
!---------------------------------------------------------------------
      if ((symmetrize).and.(nbdry.gt.1)) then
       if (limitr.gt.1000) then
        limitr=limitr-1000
        limupper=limitr
        limitr=limupper+limupper-2
        do i=2,limupper-1
          lwant=limupper-i+1
          xlim(i-1+limupper)=xlim(lwant)
          ylim(i-1+limupper)=-ylim(lwant)
        enddo
       endif
      endif
!
      if (kpressb.eq.2) pcurbd=0.0
      if (kzeroj.gt.0) then
       pcurbd=0.0
       fcurbd=0.0
      endif
      close(unit=nin)
      if (kprfit.eq.1) then
        if (npress.lt.0) then
        call getfnmu(itimeu,'k',ishot,itime,edatname)
        edatname='edat_'//edatname(2:7)// &
                       '_'//edatname(9:13)//'.pressure'
        open(unit=nin,status='old',file=edatname &
                                 )
        read (nin,edat)
        close(unit=nin)
        endif
      endif
!
      if (kprfit.eq.2) then
        if (npteth.lt.0) then
        nptef=-npteth
        npnef=-npneth
        call getfnmu(itimeu,'k',ishot,itime,edatname)
        edatname='edat_'//edatname(2:7)// &
                       '_'//edatname(9:13)//'.thomson'
        open(unit=nin,status='old',file=edatname &
                                 )
        bfract=-1.
        if (tethom(1).lt.0.0) bfract=-tethom(1)
        read (nin,edat)
        close(unit=nin)
        endif
        if (nbeam.gt.0) then
          do 43901 i=1,nbeam
           dnbeam(i)=dnbeam(i)*1.e-19_dp
43901     continue
        endif

      print *, 'reorder TS data points '  
!---------------------------------------------------------------------
!--  reorder TS data points                                         --
!---------------------------------------------------------------------
        call tsorder(npteth,zteth,dnethom,tethom,sgneth,sgteth)
        if (sgnemin.lt.0.0) sgnemin=abs(sgnemin)*dnethom(1)*1.e-19_dp &
                                    *co2cor
        do 40020 i=1,npneth
          dnethom(i)=dnethom(i)*1.e-19_dp*co2cor
          sgneth(i)=sgneth(i)*1.e-19_dp*sgnethi*co2cor
          sgneth(i)=max(sgneth(i),sgnemin)
40020   continue
        if (sgtemin.lt.0.0) sgtemin=abs(sgtemin)*tethom(1)
        temax=tethom(1)
        demax=dnethom(1)
        do 40030 i=1,npteth
          sgteth(i)=sgteth(i)*sgtethi
          if (bfract.gt.0.0) then
            tethom(i)=tethom(i)*bfract
            sgteth(i)=sgteth(i)*bfract
          endif
          sgteth(i)=max(sgteth(i),sgtemin)
          temax=max(temax,tethom(i))
          demax=max(demax,dnethom(i))
40030   continue
        if (cstabte.lt.0.0) cstabte=abs(cstabte)*100./temax
        if (cstabne.lt.0.0) cstabne=abs(cstabne)*100./demax
        if (nption.lt.0) then
        nptionf=-nption
        if (nptionf.lt.100) then
        call getfnmu(itimeu,'k',ishot,itime,edatname)
        edatname='edat_'//edatname(2:7)//'_'//edatname(9:13)//'.cer'
        open(unit=nin,status='old',file=edatname &
                                 )
        bfract=-1.
        if (tionex(1).lt.0.0) bfract=-tionex(1)
        read (nin,edat)
        close(unit=nin)
        do 40040 i=1,nption
          sigti(i)=sigti(i)*sigtii
          if (bfract.gt.0.0) then
            sigti(i)=sigti(i)*bfract
            tionex(i)=tionex(i)*bfract
          endif
40040   continue
        endif
        if (nptionf.gt.100) then
          nptionf=nptionf-100
          nption=npteth
          nptionf=nptef
          bfract=tionex(1)
          do 40050 i=1,nption
            sigti(i)=sgteth(i)
            tionex(i)=tethom(i)*bfract
            rion(i)=rteth(i)
            zion(i)=zteth(i)
40050     continue
        endif
        endif
      endif

      print *, 'read in limiter data'
!-----------------------------------------------------------------------
!---- read in limiter data                                            --
!-----------------------------------------------------------------------
      call getlim(1,xltype,xltype_180)
!
      if (iconvr.ge.0) go to 214
      iecurr=1
      ivesel=1
  214 continue
      if (iand(iout,1).ne.0) then
      write (nout,in1)
      write (nout,inwant)
      write (nout,ink)
      write (nout,ins)
      write (nout,in_msels)
      if (kwaitmse.ne.0) write (neqdsk,ina)
      write (nout,inece)
      write (nout,insxr)
      write (nout,invt)
      write (nout,inlibim)
      endif
      if (islve.gt.0) nbdry=40
!
      diamag(jtime)=1.0e-03_dp*dflux
      sigdia(jtime)=1.0e-03_dp*abs(sigdlc)
      pbinj(jtime)=pnbeam
      do 220 i=1,nsilop
        silopt(jtime,i)=coils(i)
  220 continue
      do 230 i=1,nfcoil
        fccurt(jtime,i)=brsp(i)
  230 continue
      do 240 i=1,nmtark
        tangam(jtime,i)=tgamma(i)
        tangam_uncor(jtime,i)=tgammauncor(i)
        siggam(jtime,i)=sgamma(i)
        rrgam(jtime,i)=rrrgam(i)
        zzgam(jtime,i)=zzzgam(i)
        a1gam(jtime,i)=aa1gam(i)
        a2gam(jtime,i)=aa2gam(i)
        a3gam(jtime,i)=aa3gam(i)
        a4gam(jtime,i)=aa4gam(i)
        a5gam(jtime,i)=aa5gam(i)
        a6gam(jtime,i)=aa6gam(i)
        a7gam(jtime,i)=aa7gam(i)
        a8gam(jtime,i)=0.0
  240 continue
      do i=nmtark+1, nstark
        ii = i - nmtark
        tangam(jtime,i)=tlibim(ii)
        siggam(jtime,i)=slibim(ii)
        rrgam(jtime,i)=rrrlib(ii)
        zzgam(jtime,i)=zzzlib(ii)
        a1gam(jtime,i)=aa1lib(ii)
        a2gam(jtime,i)=1.0
        a3gam(jtime,i)=0.0
        a4gam(jtime,i)=0.0
        a5gam(jtime,i)=0.0
        a6gam(jtime,i)=0.0
        a7gam(jtime,i)=0.0
        a8gam(jtime,i)=aa8lib(ii)
        fwtgam(i)=fwtlib(ii)
        swtgam(i)=fwtlib(ii)
      enddo
!
      do 91240 i=1,nmsels
        bmselt(jtime,i)=bmsels(i)
        sbmselt(jtime,i)=sbmsels(i)
        fwtbmselt(jtime,i)=fwtbmsels(i)
        rrmselt(jtime,i)=rrmsels(i)
        zzmselt(jtime,i)=zzmsels(i)
        l1mselt(jtime,i)=l1msels(i)
        l2mselt(jtime,i)=l2msels(i)
        l3mselt(jtime,i)=1.-l1msels(i)
        l4mselt(jtime,i)=l4msels(i)
        emselt(jtime,i)=emsels(i)
        semselt(jtime,i)=semsels(i)
        fwtemselt(jtime,i)=fwtemsels(i)
91240 continue
!----------------------------------------------------------------------
!  give the constraint value for matrix routine
!----------------------------------------------------------------------
      do 250 i=1,nnece
        brspece(jtime,i)=ecefit(i)
  250 continue
      brspecebz(jtime)=ecebzfit
      do 260 i=1,magpri
        expmpi(jtime,i)=expmp2(i)
  260 continue
!------------------------------------------------------------------------
!--  New E-coil connections                   LLao, 95/07/11           --
!------------------------------------------------------------------------
      if (ecurrt(3).le.-1.e10_dp) ecurrt(3)=ecurrt(1)
      if (ecurrt(5).le.-1.e10_dp) ecurrt(5)=ecurrt(1)
      if (ecurrt(4).le.-1.e10_dp) ecurrt(4)=ecurrt(2)
      if (ecurrt(6).le.-1.e10_dp) ecurrt(6)=ecurrt(2)
      do 261 i=1,nesum
        eccurt(jtime,i)=ecurrt(i)
  261 continue
      pasmat(jtime)=plasma
      curtn1(jtime)=currn1
      curc79(jtime)=currc79
      curc139(jtime)=currc139
      curc199(jtime)=currc199
      curiu30(jtime)=curriu30
      curiu90(jtime)=curriu90
      curiu150(jtime)=curriu150
      curil30(jtime)=curril30
      curil90(jtime)=curril90
      curil150(jtime)=curril150
!-----------------------------------------------------------------------
!--  + 0.01 to take care of truncation problem      03/16/91          --
!-----------------------------------------------------------------------
      timeus=itimeu
      timems=itime
      time(jtime)=timems+timeus/1000.
      bcentr(jtime)=btor
      do 262 i=1,nco2v
        denvt(jtime,i)=denv(i)
  262 continue
      do 264 i=1,nco2r
        denrt(jtime,i)=denr(i)
  264 continue
      do 267 i=1,nacoil
        accurt(jtime,i)=acoilc(i)
        caccurt(jtime,i)=acoilc(i)
  267 continue
      vloopt(jtime)=vloop
      psiref(jtime)=siref
!
      gammap=1./gammap
      gammaf=gammap
      psibry0=psibry
      if (idebug >= 2) write (6,*) 'DATA_INPUT PSIBRY0= ', psibry0
!
      xlmint=xlmin
      xlmaxt=xlmax
!
      call zlim(zero,nw,nh,limitr,xlim,ylim,rgrid,zgrid,limfag)
!-----------------------------------------------------------------------
!--  set up parameters for all modes                                  --
!-----------------------------------------------------------------------
  280 continue
      if (kfffnc.eq.8) then
        rkec=pi/(2.0*dpsiecn)
      endif
      chigam=0.0
      tchimls=0.0
!-----------------------------------------------------------------------
!-- DIII-D shot > 112000 use a new set of table for magnetic probes   --
!--        shot > 156000 new 2014 set                                 --
!--        shot >= 168191 new 2017 set                                --
!-----------------------------------------------------------------------
      if (kdata.ne.2) then
        call set_table_dir
        call efit_read_tables
      endif
!
      if (pasmat(jtime).le.-1.e3_dp) then
        negcur=1
      else
        negcur=0
      endif
      if (abs(pasmat(jtime)).gt.cutip.or.iconvr.lt.0) go to 285
      if (iconsi.eq.-1) iconsi=55
      if (ivesel.gt.10) iconsi=0
      iexcal=1
      ivacum=1
      ibunmn=0
      ierchk=0
      nxiter=1
      iconvr=3
  285 continue
!---------------------------------------------------------------------
!--  correction to 322 degree probes due to N1 coil                 --
!---------------------------------------------------------------------
      if (oldccomp) then
      if (n1coil.eq.2.and.ishot.le.108281) then
      open(unit=60,file=input_dir(1:lindir)//'n1coil.ddd', &
           status='old'                                )
      j=jtime
      do 11368 i=30,magpri67+magpri322
       read(60,*)   namedum,xxxdum,signn1(i)
        expmpi(j,i)=expmpi(j,i)-signn1(i)*curtn1(j)
11368 continue
      close(unit=60)
      endif
      endif
!---------------------------------------------------------------------
!--  correction to 322 and 67 degree probes due to C coil           --
!---------------------------------------------------------------------
      if (nccoil.eq.1.and.oldccomp) then
      open(unit=60,file=input_dir(1:lindir)//'ccoil.ddd', &
           status='old'                                )
      j=jtime
      do i=30,magpri67+magpri322
       read(60,*)   namedum,signc139
        expmpi(j,i)=expmpi(j,i)-signc139*curc139(j)
      enddo
      do i=1,29
       read(60,*)   namedum,signc79
        expmpi(j,i)=expmpi(j,i)-signc79*curc79(j)
      enddo
      close(unit=60)
      endif
!
      if (ifitvs.gt.0) then
        ivesel=5
      endif
      if (.not.fitsiref) then
      if (iecurr.gt.0.or.nslref.lt.0) then
        do 287 m=1,nsilop
          silopt(jtime,m)=silopt(jtime,m)+psiref(jtime)
  287   continue
        psirefs(jtime)=psiref(jtime)
        psiref(jtime)=0.
      endif
      endif
      kconvr=iconvr
      do 290 kk=1,nwnh
        www(kk)=zero(kk)
  290 continue
!
      fwtref=fwtsi(iabs(nslref))
      if (kersil.eq.2) go to 325
      if (kersil.eq.3) go to 325
      do 300 m=1,nsilop
        tdata1=serror*abs(silopt(jtime,m))
        tdata2=abs(psibit(m))*vbit
        tdata=max(tdata1,tdata2)
        sigmafl0(m)=tdata
        if (tdata.gt.1.0e-10_dp) fwtsi(m)=fwtsi(m)/tdata**nsq
        if (tdata.le.1.0e-10_dp) fwtsi(m)=0.0
  300 continue
!----------------------------------------------------------------------
!--  signal at psi loop # NSLREF is used as reference                --
!----------------------------------------------------------------------
      if (abs(psibit(iabs(nslref))).gt.1.0e-10_dp) go to 340
      coilmx=abs(silopt(jtime,1))
      do 320 i=2,nsilop
        abcoil=abs(silopt(jtime,i))
        coilmx=max(abcoil,coilmx)
  320 continue
      fwtsi(iabs(nslref))=1.0/coilmx**nsq/serror**nsq*fwtref
      go to 340
!
  325 continue
 !     if (kdata.ne.2) then
 !     if (jtime.le.1) then
 !       open(unit=80,status='old',file=table_di2(1:ltbdi2)//'dprobe.dat')
 !       rsi(1)=-1.
 !       read (80,in3)
 !       read (80,10200) (rf(i),zf(i),wf(i),hf(i),af(i),af2(i), &
 !               i=1,mfcoil)
 !       if (rsi(1).lt.0.) &
 !       read (80,10200) (rsi(i),zsi(i),wsi(i),hsi(i),as(i),as2(i), &
 !               i=1,nsilop)
 !       read (80,10220) (re(i),ze(i),we(i),he(i),ecid(i), &
 !                                       i=1,necoil)
 !       if (ifitvs.gt.0.or.icutfp.eq.2) then
 !         read (80,10200) (rvs(i),zvs(i),wvs(i),hvs(i), &
 !                                       avs(i),avs2(i),i=1,nvesel)
 !       endif
 !       close(unit=80)
 !     endif
 !     endif
!-----------------------------------------------------------------------
!--  Fourier expansion of vessel sgments                              --
!-----------------------------------------------------------------------
      if (ifitvs.gt.0. .and. nfourier.gt.1) then
      do i=1,nvesel
      if(rvs(i).ge.1.75_dp.and.zvs(i).ge.0.) &
      thetav(i)=dasin(zvs(i)/sqrt ((rvs(i)-1.75_dp)**2+(zvs(i))**2))
      if(rvs(i).lt.1.75_dp.and.zvs(i).ge.0.) &
      thetav(i)=pi-dasin(zvs(i)/sqrt ((rvs(i)-1.75_dp)**2+(zvs(i))**2))
      if(rvs(i).lt.1.75_dp.and.zvs(i).lt.0.) &
      thetav(i)=pi-dasin(zvs(i)/sqrt ((rvs(i)-1.75_dp)**2+(zvs(i))**2))
      if(rvs(i).ge.1.75_dp.and.zvs(i).lt.0.) &
      thetav(i)=2*pi+dasin(zvs(i)/sqrt ((rvs(i)-1.75_dp)**2+(zvs(i))**2))
        do j=1,nfourier
          sinta(j,i)=sin(thetav(i)*j)
          costa(j,i)=cos(thetav(i)*j)
        enddo
      enddo
      do i=1,(2*nfourier+1)
        do j=1,nvesel
          if(i.eq.1) vecta(i,j)=1.0
          if(i.gt.1.and.i.le.(nfourier+1)) vecta(i,j)=costa(i,j)
          if(i.gt.(nfourier+1))  vecta(i,j)=sinta(i,j)
        enddo
      enddo
      endif
!
!     if (brsptu(1).le.-1.e-20_dp) &
!        brsp(1:nfcoil)=brsptu(1:nfcoil)*turnfc(1:nfcoil)
      if (brsptu(1).gt.-1.e-20_dp) &
         brsp(1:nfcoil)=brsptu(1:nfcoil)*turnfc(1:nfcoil)
      reflux=silopt(jtime,iabs(nslref))
      do 330 m=1,nsilop
        tdata1=errsil*abs(silopt(jtime,m)-reflux)
        tdata2=sicont*rsi(m)*abs(pasmat(jtime))
        tdata=max(tdata1,tdata2)
        tdata2=abs(psibit(m))*vbit
        tdata=max(tdata,tdata2)
        sigmafl0(m)=tdata
        if (tdata.gt.1.0e-10_dp) fwtsi(m)=fwtsi(m)/tdata**nsq
        if (tdata.le.1.0e-10_dp) fwtsi(m)=0.0
  330 continue
!----------------------------------------------------------------------
!--  signal at psi loop #8 is set to zero and used as reference      --
!----------------------------------------------------------------------
      if (kersil.ne.3) then
        fwtsi(iabs(nslref))=1.0/ersil8**nsq*fwtref
!----------------------------------------------------------------------
!-- New option for reference flux loop uncertainty                   --
!----------------------------------------------------------------------
      else
        m=iabs(nslref)
        tdata1=errsil*abs(silopt(jtime,m))
        tdata2=sicont*rsi(m)*abs(pasmat(jtime))
        tdata=max(tdata1,tdata2)
        tdata2=abs(psibit(m))*vbit
        tdata=max(tdata,tdata2)
!       sigmafl0(m)=tdata
        if (tdata.gt.1.0e-10_dp) fwtsi(m)=fwtref/tdata**nsq
        if (tdata.le.1.0e-10_dp) fwtsi(m)=0.0
      endif
      sigmafl0(m)=tdata
!
  340 continue
      fwtref=fwtsi(iabs(nslref))
      do 350 m=1,magpri
        tdata1=serror*abs(expmpi(jtime,m))
        tdata2=abs(bitmpi(m))*vbit
        tdata=max(tdata1,tdata2)
        sigmamp0(m)=tdata
        if (tdata.gt.1.0e-10_dp) fwtmp2(m)=fwtmp2(m)/tdata**nsq
        if (tdata.le.1.0e-10_dp) fwtmp2(m)=0.0
  350 continue
      do 400 m=1,nstark
        tdata=abs(siggam(jtime,m))
        if (tdata.gt.1.0e-10_dp) fwtgam(m)=fwtgam(m)/tdata**nsq
        if (tdata.le.1.0e-10_dp) fwtgam(m)=0.0
  400 continue
!
      if (idebug >= 2) then
         write (6,*) 'DATA fwtbmselt = ',(fwtbmselt(jtime,i),i=1,nmsels)
         write (6,*) 'DATA sbmselt = ',(sbmselt(jtime,i),i=1,nmsels)
      endif
      do i=1,nmsels
        tdata=abs(sbmselt(jtime,i))
        if (tdata.gt.1.0e-10_dp) fwtbmselt(jtime,i)=fwtbmselt(jtime,i)/tdata**nsq
        if (tdata.le.1.0e-10_dp) fwtbmselt(jtime,i)=0.0
      enddo
      do i=1,nmsels
        tdata=abs(semselt(jtime,i))
        if (tdata.gt.1.0e-10_dp) fwtemselt(jtime,i)=fwtemselt(jtime,i)/tdata**nsq
        if (tdata.le.1.0e-10_dp) fwtemselt(jtime,i)=0.0
      enddo
      if (idebug >= 2) write (6,*) 'DATA fwtbmselt = ', (fwtbmselt(jtime,i),i=1,nmsels)
!
      do 402 m=1,nfcoil
        tdata1=serror*abs(fccurt(jtime,m))
        tdata2=abs(bitfc(m))*vbit
        tdata=max(tdata1,tdata2)
        sigmaf0(m)=tdata
        if (tdata.gt.1.0e-10_dp) fwtfc(m)=fwtfc(m)/tdata**nsq
        if (tdata.le.1.0e-10_dp) fwtfc(m)=0.0
  402 continue
      if (iecurr.eq.2) then
      do m=1,nesum
        tdata1=serror*abs(ecurrt(m))
        tdata2=abs(bitec(m))*vbit
        tdata=max(tdata1,tdata2)
        sigmae0(m)=tdata
        if (tdata.gt.1.0e-10_dp) fwtec(m)=fwtec(m)/tdata**nsq
        if (tdata.le.1.0e-10_dp) fwtec(m)=0.0
      enddo
      endif
      tdata1=serror*abs(pasmat(jtime))
      tdata2=abs(bitip)*vbit
      tdata=max(tdata1,tdata2)
      sigmaip0=tdata
      if (tdata.gt.1.0e-10_dp) fwtcur=fwtcur/tdata**nsq
      if (tdata.le.1.0e-10_dp) fwtcur=0.0
!----------------------------------------------------------------------
!-- diamagetic flux                                                  --
!----------------------------------------------------------------------
      tdata=abs(sigdia(jtime))
      if (tdata.gt.1.0e-10_dp) fwtdlc=fwtdlc/tdata**nsq
!
      if (sidif.le.-1.0e+10_dp) then
        sidif=tmu*pasmat(jtime)*rcentr/2.0
      endif
      errcut=max(ten2m3,error*10.)
      fbrdy=bcentr(jtime)*rcentr/tmu
      constf2=darea*tmu/2.0/twopi
      fcentr=fbrdy
      rbetap=(1.-betap0)/betap0
      rbetaw=betapw0/betap0
      fconst=rzero**2*rbetap
      pbetap=betap0/(1.0-betap0)/rzero**2
      emf=emp
      enf=enp
      kpcurn=kppcur+kffcur
      nfnpcr=nfcoil+kpcurn
      nbase=nfcoil+kppcur
      nfnwcr=nfnpcr
      if (kvtor.gt.0) then
        nfnwcr=nfnwcr+kwwcur
        kwcurn=kpcurn+kwwcur
      else
        kwcurn=kpcurn
      endif
      nqaxis=0
      if (fwtqa.gt.1.0e-03_dp) nqaxis=1
      nparam=nfnwcr
      if (kprfit.gt.0) nparam=nparam+1
      if (fitdelz) nparam=nparam+1
      if (fitsiref) nparam=nparam+1
      if (kedgep.gt.0) nparam=nparam+1
      if (kedgef.gt.0) nparam=nparam+1
      if (fwtqa.gt.0.0) fwtqa=fwtqa/errorq
      if (fwtbp.gt.0.0) fwtbp=fwtbp/errorq
      if (fbetap.gt.0.0) betap0=fbetap
!
      ipsi(jtime)=0
      do 420 i=1,nsilop
        if (fwtsi(i).gt.0.0) ipsi(jtime)=ipsi(jtime)+1
  420 continue
      ifc(jtime)=0
      do 428 i=1,nfcoil
        if (fwtfc(i).gt.0.0) ifc(jtime)=ifc(jtime)+1
  428 continue
      iec(jtime)=0
      do i=1,nesum
        if (fwtec(i).gt.0.0) iec(jtime)=iec(jtime)+1
      enddo
      imag2(jtime)=0
      do 460 i=1,magpri
        if (fwtmp2(i).gt.0.0) imag2(jtime)=imag2(jtime)+1
  460 continue
      kmtark=0
      klibim=0
      do 463 i=1,nmtark
        if (fwtgam(i).gt.0.0) kmtark=kmtark+1
  463 continue
      do i=nmtark+1, nstark
        if (fwtgam(i).gt.0.0) klibim=klibim+1
      enddo
      kstark=kmtark+klibim
!
      mmbmsels=0
      mmemsels=0
      do i=1,nmsels
        if (fwtbmselt(jtime,i).gt.1.e-06_dp) mmbmsels=mmbmsels+1
        if (fwtemselt(jtime,i).gt.1.e-06_dp) mmemsels=mmemsels+1
      enddo
      if (jdebug.eq.'MSEL') then
         write (6,*) 'DATA mmbmsels = ', mmbmsels
         write (6,*) 'bmselt ',(bmselt(1,i),i=1,nmsels)
         write (6,*) 'iermselt ',(iermselt(1,i),i=1,nmsels)
      endif
!
      iplasm(jtime)=0
      if (fwtcur.gt.0.0) iplasm(jtime)=1
      idlopc(jtime)=0
      if (fwtdlc.gt.0.0) idlopc(jtime)=1
      cpasma(jtime)=pasmat(jtime)
      if (iconvr.eq.3.and.ivesel.eq.1) then
        do 469 i=1,nvesel
          cpasma(jtime)=cpasma(jtime)-vcurrt(i)
  469   continue
      endif
      do 500 n=1,nsilop
        csilop(n,jtime)=silopt(jtime,n)
  500 continue
      itime=time(jtime)
      timems=itime
      timeus=(time(jtime)-timems)*1000.
      timeus=timeus+0.4_dp
      itimeu=timeus
!-----------------------------------------------------------------------
!-- correction for truncation                                         --
!-----------------------------------------------------------------------
      if (itimeu.ge.990) then
        itime=itime+1
        itimeu=0
        time(jtime)=itime
      endif
      csiref=psiref(jtime)
      do 510 i=1,nesum
        ecurrt(i)=eccurt(jtime,i)
        cecurr(i)=ecurrt(i)
  510 continue
 !     if (kdata.ne.2) then
 !     if ((iecurr.le.0).or.(idodo.gt.0)) go to 520
 !     open(unit=nrsppc,status='old',form='unformatted', &
 !          file=table_dir(1:ltbdir)//'re'//trim(ch1)//trim(ch2)//'.ddd')
 !     read (nrsppc) rsilec
 !     read (nrsppc) rmp2ec
 !     read (nrsppc) gridec
 !     close(unit=nrsppc)
 !     idodo=1
 ! 520 continue
!
!      if ((ivesel.le.0).or.(idovs.gt.0)) go to 525
!      open(unit=nrsppc,status='old',form='unformatted', &
!           file=table_dir(1:ltbdir)//'rv'//trim(ch1)//trim(ch2)//'.ddd')
!      read (nrsppc) rsilvs
!      read (nrsppc) rmp2vs
!      read (nrsppc) gridvs
!      close(unit=nrsppc)
!      idovs=1
!      if (ivesel.le.10) go to  525
!      open(unit=nffile,status='old',form='unformatted', &
!           file=table_dir(1:ltbdir)//'fc'//trim(ch1)//trim(ch2)//'.ddd')
!      read (nffile) rfcfc
!      close(unit=nffile)
!      endif
!
      do 522 i=1,nfcoil
        if (rsisfc(i).le.-1.0) &
        rsisfc(i)=turnfc(i)**2*twopi*rf(i)/wf(i)/hf(i)*zetafc
  522 continue
  525 continue
!-----------------------------------------------------------------------
!-- read in the advance divertor coil response functions if needed    --
!-----------------------------------------------------------------------
      if ((iacoil.gt.0).and.(idoac.eq.0)) then
        open(unit=nrsppc,status='old',form='unformatted', &
             file=table_di2(1:ltbdi2)//'ra'//trim(ch1)//trim(ch2)//'.ddd')
        read (nrsppc) gridac
        read (nrsppc) rsilac
        read (nrsppc) rmp2ac
        close(unit=nrsppc)
        idoac=1
      endif
!--------------------------------------------------------------------
!--  optional vertical feedback                                    --
!--------------------------------------------------------------------
      if (isetfb.ne.0.and.idofb.le.0) then
        open(unit=mcontr,status='old',form='unformatted', &
             file=table_di2(1:ltbdi2)//'ef'//trim(ch1)//trim(ch2)//'.ddd')
        read (mcontr) mw,mh
        read (mcontr) rgrid,zgrid
        read (mcontr) grdfdb
!       read (mcontr) grddum
        close(unit=mcontr)
!---------------------------------------------------------------------
!--   DO DOUBLE PRECISION SUM OF GRDFDB IN M=1 CONFIGURATION        --
!---------------------------------------------------------------------
        do 40001 i=1,nw
        do 40001 j=1,nh
          kk=(i-1)*nh+j
          gsum=0
          do 39999 mmf=1,nfbcoil/2
            gsum=gsum+grdfdb(kk,mmf)
            gsum=gsum-grdfdb(kk,mmf+nfbcoil/2)
39999     continue
          grdfdb(kk,1)=gsum
40001   continue
        idofb=1
      endif
!---------------------------------------------------------------------
!--  polarimetry ?                                                  --
!---------------------------------------------------------------------
      if (kstark.gt.0.or.kdomse.gt.0) then
         call setstark(jtime)
      endif
      if (kdomse.gt.0.and.keecur.gt.0) then
        do i=1,keecur
         if (keefnc.le.2) then
          cerer(i)=eebdry(i)
         elseif (keefnc.eq.6) then
          cerer(2*i-1)=eebdry(i)
          cerer(2*i)=ee2bdry(i)
         endif
        enddo
      endif
!---------------------------------------------------------------------
!--  ECE --set kece, kecebz=0 before call setece                    --
!---------------------------------------------------------------------
      kece=0
      kecebz=0
      sigrid(1)=0.0
      sigrid(nw)=1.0
      do 80 i=2,nw-1
        sigrid(i)=1./real(nw-1,dp)*(i-1)
  80  continue
!--------------------------------------------------------------------
!-- kinputece=1, get Te, fe, error array from ECE data routine
!--     hecedata.for with get_hece.for,fftabl.11,fts.pst( copy from
!--     /u/austin/efit/hecefit/     (MAX AUSTIN))
!--     necein,teecein0(necein),feece0(necein),errorece0(necein)
!--     fe(GHz), Te(Kev)
!--( when kinputece>1,data from K-file )
!--    (data order: from low field to high field)
!--   change data order :( from high field to low field)
!--------------------------------------------------------------------
      do k=1,necein
         kk=necein-k+1
         if (kfitece.eq.3) kk = k
         feece(kk)=feece0(k)
         teecein(kk)=teecein0(k)
         errorece(kk)=errorece0(k)
      enddo
!---------------------------------------------------------------------
!--  toroidal rotation ? Then set up geometric parameters           --
!---------------------------------------------------------------------
      if (kvtor.gt.0) then
        do i=1,nw
         rgrvt(i)=(rgrid(i)/rvtor)**2
         rgrvt(i)=rgrvt(i)-1.
         rgsvt(i)=rgrid(i)*rgrvt(i)
        enddo
      endif
!----------------------------------------------------------------------
!-- make filement Green's tables only                                --
!----------------------------------------------------------------------
      if (iconvr.ge.0) go to 990
      if (iconvr.le.-20) go to 700
      mx=iabs(iconvr)
      do 690 k=1,mx
        if (aelip.le.0.0) then
        i=irfila(k)
        j=jzfila(k)
        else
        th=twopi*(k-1)/real(mx,dp)
        rmx(k)=relip-aelip*cos(th)
        zmx(k)=zelip+eelip*aelip*sin(th)
        ix=1
        if (k.gt.(mx/2+1)) ix=2
        i=(rmx(k)-rgrid(1))/drgrid+1
        j=(zmx(k)-zgrid(1))/dzgrid+ix
        zdif=zmx(k)-zgrid(j)
        if (abs(zdif).gt.0.6_dp*dzgrid) then
          if (zdif.gt.0.0) j=j+1
          if (zdif.lt.0.0) j=j-1
        endif
        endif
        irfila(k)=i
        jzfila(k)=j
        rmx(k)=rgrid(i)
        zmx(k)=zgrid(j)
        kk=(i-1)*nh+j
        do 600 m=1,nsilop
  600     rsilpf(m,k)=gsilpc(m,kk)
        do 605 m=1,magpri
  605     rmp2pf(m,k)=gmp2pc(m,kk)
        do 610 ii=1,nw
        do 610 jj=1,nh
          kkkk=(ii-1)*nh+jj
          mj=iabs(j-jj)+1
          mk=(i-1)*nh+mj
  610     gridpf(kkkk,k)=gridpc(mk,ii)
  690 continue
      mw=nw
      mh=nh
      open(unit=nffile,status='old',form='unformatted', &
           file='rpfxx.dat',err=12925)
      close(unit=nffile,status='delete')
12925 continue
      open(unit=nffile,status='new',form='unformatted', &
           file='rpfxx.dat')
      write (nffile) mx,rmx,zmx
      write (nffile) rsilpf
      write (nffile) rmp2pf
      write (nffile) mw,mh,rgrid,zgrid
      write (nffile) gridpf
      close(unit=nffile)
      write (nttyo,6550) (irfila(i),i=1,mx)
      write (nttyo,6555) (jzfila(i),i=1,mx)
      go to 2000
!
  700 continue
      if (aelip.gt.0.0) then
      do 710 i=1,nw
      do 710 j=1,nh
        kk=(i-1)*nh+j
        erho=sqrt((rgrid(i)-relip)**2+((zgrid(j)-zelip)/eelip)**2)
        xpsi(kk)=(erho/aelip)**2
  710 continue
      else
      endif
      do 720 m=1,nsilop
        wsilpc(m)=0.0
  720 continue
      do 722 m=1,magpri
        wmp2pc(m)=0.0
  722 continue
      do 724 m=1,nfcoil
        wfcpc(m)=0.0
  724 continue
      do 726 m=1,nesum
        wecpc(m)=0.0
  726 continue
      do 728 m=1,nvesel
        wvspc(m)=0.0
  728 continue
      do 730 m=1,nwnh
  730   wgridpc(m)=0.0
      wpcpc=0.0
      npc=0
!
          open(unit=nffile,status='old',form='unformatted', &
               file=table_di2(1:ltbdi2)//'fc'//trim(ch1)//trim(ch2)//'.ddd')
      read (nffile) rfcfc
      read (nffile) rfcpc
      close(unit=nffile)
!
      do 800 i=1,nw
      do 800 j=1,nh
        kk=(i-1)*nh+j
        if (xpsi(kk).lt.0.0.or.xpsi(kk).gt.1.0) go to 800
        npc=npc+1
        do 750 m=1,nsilop
  750     wsilpc(m)=wsilpc(m)+gsilpc(m,kk)
        do 755 m=1,magpri
  755     wmp2pc(m)=wmp2pc(m)+gmp2pc(m,kk)
        do 760 m=1,nfcoil
  760     wfcpc(m)=wfcpc(m)+rfcpc(m,kk)
        do 765 m=1,nesum
  765     wecpc(m)=wecpc(m)+gridec(kk,m)
        do 770 m=1,nvesel
  770     wvspc(m)=wvspc(m)+gridvs(kk,m)
        do 780 ii=1,nw
        do 780 jj=1,nh
          kkkk=(ii-1)*nh+jj
          mj=iabs(j-jj)+1
          mk=(i-1)*nh+mj
          wgridpc(kkkk)=wgridpc(kkkk)+gridpc(mk,ii)
          if (xpsi(kkkk).lt.0.0.or.xpsi(kkkk).gt.1.0) go to 780
          wpcpc=wpcpc+gridpc(mk,ii)
  780   continue
  800 continue
        xnpc=real(npc,dp)
        do 810 m=1,nsilop
  810     wsilpc(m)=wsilpc(m)/xnpc
        do 815 m=1,magpri
  815     wmp2pc(m)=wmp2pc(m)/xnpc
        do 820 m=1,nfcoil
  820     wfcpc(m)=wfcpc(m)/xnpc
        do 825 m=1,nesum
  825     wecpc(m)=wecpc(m)/xnpc
        do 830 m=1,nvesel
  830     wvspc(m)=wvspc(m)/xnpc
        do 835 m=1,nwnh
  835     wgridpc(m)=wgridpc(m)/xnpc
      wpcpc=wpcpc/xnpc**2
!
      open(unit=nffile,status='old',form='unformatted', &
           file='rpcxx.dat',err=12926)
      close(unit=nffile,status='delete')
12926 continue
      open(unit=nffile,status='new',form='unformatted', &
           file='rpcxx.dat')
      write (nffile) wsilpc
      write (nffile) wmp2pc
      write (nffile) wfcpc
      write (nffile) wecpc
      write (nffile) wvspc
      write (nffile) wpcpc
      write (nffile) wgridpc
      write (nffile) npc
      write (nffile) relip,zelip,aelip,eelip
      close(unit=nffile)
      write (nttyo,6557) npc
      go to 2000
!
  990 continue
      if (ivesel.gt.0) call vescur(jtime)
      if (nbdry.le.0) go to 2000
      if (islve.le.0) go to 1180
!------------------------------------------------------------------------------
!--  Solve equilibrium                                                      --
!------------------------------------------------------------------------------
      icurrt=1
      iecurr=0
      iconvr=3
      itrace=0
      ierchk=0
      nextra=0
      ssrm=srm
      srm=abs(srm)
      seee=1./srm/sqrt(salpha)
      rbetaw=0.0
!---------------------------------------------------------------------
!-- no rotation                                                     --
!---------------------------------------------------------------------
      if (kvtor.le.0) then
        scc1=sqrt(2./sbeta)/srm**2
        if (ssrm.lt.0.0) saaa=xlmint/srm/sqrt(1.-2.*scc1)
        if (ssrm.gt.0.0) saaa=xlmax/srm/sqrt(1.+2.*scc1)
        srma=srm*saaa
        dth=twopi/real(nbdry,dp)
        do 1120 i=1,nbdry
          th=(i-1)*dth
          rbdry(i)=srma*sqrt(1.-2.*scc1*cos(th))
          zbdry(i)=sin(th)
          zbdry(i)=saaa*zbdry(i)*seee
 1120   continue
      else
!----------------------------------------------------------------------
!--  toroidal rotation                                               --
!----------------------------------------------------------------------
        saaa=0.50_dp
        do i=1,nw
          rgrids(i)=rgrid(i)/saaa
        enddo
        do i=1,nh
          zgrids(i)=zgrid(i)/saaa
        enddo
        do i=1,nw
          xrm2=(rgrids(i)-srm)*(rgrids(i)+srm)
          xrm2=xrm2*xrm2
          xrvt=(rgrids(i)/srm)**2-1.
          do j=1,nh
            kk=(i-1)*nh+j
            psi(kk)=sbeta/8.*xrm2+(zgrids(j)/seee)**2 &
                    +sbetaw/24.*xrm2*xrvt
          enddo
        enddo
        siwant=1.0
        drgrids=rgrids(2)-rgrids(1)
        dzgrids=zgrids(2)-zgrids(1)
        xmin=rgrids(3)
        xmax=rgrids(nw-2)
        ymin=zgrids(3)
        ymax=zgrids(nh-2)
        npack=1
        rnow=0.5_dp*(rgrids(1)+rgrids(nw))
        znow=0.0
        call surfac(siwant,psi,nw,nh,rgrids,zgrids,xout,yout,nfound, &
                    npoint,drgrids,dzgrids,xmin,xmax,ymin,ymax,npack, &
                    rnow,znow,negcur,kerror)
        if (kerror.gt.0) return
        xmin=xout(1)
        xmax=xmin
        do i=2,nfound
          if (xout(i).lt.xmin) xmin=xout(i)
          if (xout(i).gt.xmax) xmax=xout(i)
        enddo
        if (ssrm.lt.0.0) saaa=xlmint/xmin
        if (ssrm.gt.0.0) saaa=xlmax/xmax
        nskip=nfound/mbdry+1
        j=0
        do i=1,nfound,nskip
          j=j+1
          rbdry(j)=xout(i)*saaa
          zbdry(j)=yout(i)*saaa
        enddo
        nbdry=j
        srma=srm*saaa
        rvtor=srma
        rbetaw=sbetaw/sbeta
      endif
 1180 continue
!-----------------------------------------------------------------------------
!--   set up plasma response                                                --
!-----------------------------------------------------------------------------
      do 72070 m=1,nbdry
        do 72060 i=1,nw
        rdif=rbdry(m)-rgrid(i)
        do 72060 j=1,nh
          k=(i-1)*nh+j
          zdif=zbdry(m)-zgrid(j)
          rsum=rdif**2+zdif**2
          if (rsum.gt.dselsum) go to 72054
!          mk=(i-1)*nh+1
!          rbdrpc(m,k)=gridpc(mk,i)
          zdif=dselsum
          rselsum=rgrid(i)-dselsum
          rbdrpc(m,k)=psical(rbdry(m),rselsum,zdif)*tmu
          go to 72056
72054     continue
          rbdrpc(m,k)=psical(rbdry(m),rgrid(i),zdif)*tmu
72056    continue
72060   continue
72070 continue
!-----------------------------------------------------------------------------
!--   SOL plasma response                                                   --
!-----------------------------------------------------------------------------
      if (nsol.gt.0) then
        do m=1,nsol
          do i=1,nw
            rdif=rsol(m)-rgrid(i)
            do j=1,nh
              k=(i-1)*nh+j
              zdif=zsol(m)-zgrid(j)
              rsum=rdif**2+zdif**2
              if (rsum.le.dselsum) then
!                mk=(i-1)*nh+1
!                rsolpc(m,k)=gridpc(mk,i)
                 zdif=dselsum
                 rselsum=rgrid(i)-dselsum
                 rsolpc(m,k)=psical(rsol(m),rselsum,zdif)*tmu
              else
                rsolpc(m,k)=psical(rsol(m),rgrid(i),zdif)*tmu
              endif
            enddo
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!--  set up parameters for fixed boundary calculations                --
!-----------------------------------------------------------------------
      if (ifref.eq.-1) ifref=1
      if (nbdry.gt.1) then
        delx2=(rbdry(1)-rbdry(nbdry))**2
        dely2=(zbdry(1)-zbdry(nbdry))**2
        if ((delx2+dely2).le.1.0e-08_dp) nbdry=nbdry-1
      endif
      if (nbdry.lt.10) go to 1210
      xmin=rbdry(1)
      xmax=xmin
      ymin=zbdry(1)
      ymax=ymin
      do 1200 i=2,nbdry
        xmin=min(xmin,rbdry(i))
        xmax=max(xmax,rbdry(i))
        ymin=min(ymin,zbdry(i))
        ymax=max(ymax,zbdry(i))
 1200 continue
      relip=(xmin+xmax)/2.
      zelip=(ymin+ymax)/2.
      aelip=(xmax-xmin)/2.
      eelip=(ymax-ymin)/(xmax-xmin)
 1210 continue
      if (cfcoil.lt.0.) cfcoil=100./pasmat(jtime)*abs(cfcoil)
      if (cupdown.lt.0.) cupdown=100./pasmat(jtime)*abs(cupdown)
!-----------------------------------------------------------------------
!--  symmetrize  F coil responses if needed                           --
!-----------------------------------------------------------------------
      if ((symmetrize).and.(nbdry.gt.1)) then
        do i=1,nw
         do j=1,nh
          kkl=(i-1)*nh+j
          kku=i*nh-j+1
          do m=nfcoil/2+1,nfcoil
            gridfc(kkl,m)=gridfc(kku,m-nfcoil/2)
          enddo
         enddo
        enddo
      endif
!-----------------------------------------------------------------------
!--  interpolate to get boundary response functions, first F coils    --
!-----------------------------------------------------------------------
      do 1500 n=1,nfcoil
      call sets2d(gridfc(1,n),c,rgrid,nw,bkx,lkx,zgrid,nh,bky, &
                       lky,wk,ier)
        do 1450 i=1,nbdry
        call seva2d(bkx,lkx,bky,lky,c,rbdry(i),zbdry(i),pds,ier,n111)
        rbdrfc(i,n)=pds(1)
 1450   continue
        if (nsol.gt.0) then
          do i=1,nsol
            call seva2d(bkx,lkx,bky,lky,c,rsol(i),zsol(i),pds,ier,n111)
            rsolfc(i,n)=pds(1)
          enddo
        endif
 1500 continue
!----------------------------------------------------------------------
!--  make sure interpolations are symmetrized                        --
!----------------------------------------------------------------------
      if ((symmetrize).and.(nbdry.gt.1)) then
        do i=1,nbryup
        if (ilower(i).ne.-1) then
        do j=nfcoil/2 +1, nfcoil
           jupper=j-nfcoil/2
           rbdrfc(i,j)=rbdrfc(ilower(i),jupper)
           rbdrfc(ilower(i),j)=rbdrfc(i,jupper)
        enddo
        endif
        enddo
      endif
!-----------------------------------------------------------------------
!--  advance divertor coil                                            --
!-----------------------------------------------------------------------
      if (iacoil.gt. 0) then
      do 1539 n=1,nacoil
      call sets2d(gridac(1,n),c,rgrid,nw,bkx,lkx,zgrid,nh,bky, &
                       lky,wk,ier)
        do 1517 i=1,nbdry
        call seva2d(bkx,lkx,bky,lky,c,rbdry(i),zbdry(i),pds,ier,n111)
        rbdrac(i,n)=pds(1)
 1517   continue
 1539 continue
      endif
!-----------------------------------------------------------------------
!-- Ohmic coils                                                       --
!-----------------------------------------------------------------------
      if (iecurr.le.0) go to 1710
      do 1700 n=1,nesum
      call sets2d(gridec(1,n),c,rgrid,nw,bkx,lkx,zgrid,nh,bky, &
                       lky,wk,ier)
        do 1650 i=1,nbdry
        call seva2d(bkx,lkx,bky,lky,c,rbdry(i),zbdry(i),pds,ier,n111)
        rbdrec(i,n)=pds(1)
 1650   continue
        if (nsol.gt.0) then
          do i=1,nsol
            call seva2d(bkx,lkx,bky,lky,c,rsol(i),zsol(i),pds,ier,n111)
            rsolec(i,n)=pds(1)
          enddo
        endif
 1700 continue
!-----------------------------------------------------------------------
!-- now vessel                                                        --
!-----------------------------------------------------------------------
 1710 if (ivesel.le.0) go to 2000
      do 1800 n=1,nvesel
      call sets2d(gridvs(1,n),c,rgrid,nw,bkx,lkx,zgrid,nh,bky, &
                       lky,wk,ier)
        do 1750 i=1,nbdry
        call seva2d(bkx,lkx,bky,lky,c,rbdry(i),zbdry(i),pds,ier,n111)
        rbdrvs(i,n)=pds(1)
 1750   continue
 1800 continue
!
 2000 continue
!
      DEALLOCATE(rgrids,zgrids,gridpf,gwork)
!
      return
 5000 format (2e12.6)
 5020 format (1x,4e16.9)
 6550 format (/,' ir = ',10i4)
 6555 format (' iz = ',10i4)
 6557 format (/,' npc = ',i4)
10200 format (6e12.6)
10220 format (5e10.4)
      end
