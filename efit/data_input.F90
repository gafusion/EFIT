#include "config.f"
!********************************************************************** 
!**                                                                  ** 
!**     SUBPROGRAM DESCRIPTION:                                      ** 
!**          data sets up the magnetic data and weighting arrays.    ** 
!**                                                                  ** 
!**     CALLING ARGUMENTS:                                           ** 
!**                                                                  ** 
!**     RECORD OF MODIFICATION:                                      ** 
!**          29/06/83..........first created                         ** 
!**          24/07/85..........revised                               ** 
!**          23/04/04...JAL iplcout added to namelist                ** 
!**          01/08/07...DPB namelist for mag uncertainty added       ** 
!**                                                                  ** 
!********************************************************************** 
      subroutine data_input(jtime,kconvr,ktime,mtear,kerror) 
      use commonblocks,only: c,wk,copy,bkx,bky,wgridpc,rfcpc 
      use set_kinds 
      include 'eparm.inc' 
      include 'modules2.inc' 
      include 'modules1.inc' 
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
      real*8,dimension(nppcur*nppcur) ::  calpa_flat
      real*8,dimension(nffcur*nffcur) ::  cgama_flat
      integer :: nw_ext,nh_ext
      integer :: dims,fitsiref_int,scalea_int,fitfcsum_int,fitdelz_int
      integer :: writepc_int,oldccomp_int,oldcomp_int
      integer :: symmetrize_int,backaverage_int,shape_ext_int,fixpp_int
      real*8 :: c_ext,dr_ext,dz_ext,rc_ext,zc_ext,a_ext 
      real*8 :: eup_ext,elow_ext,dup_ext,dlow_ext,setlim_ext 
      real*8 :: r0min,r0max,z0min,z0max,zr0min,zr0max,rz0min,rz0max 
      real*8 :: r0ave,z0ave,a0ave,e0top,e0bot,d0top,d0bot 
      character*10 case_ext(6) 
      character*50 edatname 
      character*82 table_nam 
      character*10 namedum 
      character*2 :: reflect_ext 
      logical :: shape_ext 
      logical :: file_stat
      !real*4 spatial_avg_ham(nmtark,ngam_vars,ngam_u,ngam_w) 
      data nsq/1/ 
      data ersil8/1.0e-03_dp/,currn1/0.0/ 
      data idodo/0/,idovs/0/,zetafc/2.5e-08_dp/ 
      data co2cor/1.0/,idoac/0/,fq95/0.0/ 
      data mcontr/35/ 
      data ten2m3/1.0e-03_dp/ 
      data idtime/0/,itimeb/0/ 
      save idodo,idovs,idoac 

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
      namelist/profile_ext/npsi_ext,pprime_ext,ffprim_ext,psin_ext, & 
      geqdsk_ext,sign_ext,scalepp_ext,scaleffp_ext,shape_ext,dr_ext, & 
      dz_ext,rc_ext,zc_ext,a_ext,eup_ext,elow_ext,dup_ext,dlow_ext, & 
      setlim_ext,reflect_ext,fixpp 
! 
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
      sbmsels=0.0
      tlibim=0.0
      slibim=0.0
      aa1lib=0.0
      aa8lib=0.0
      dflux=0.0
      coils=0.0
      acoilc=0.0
      expmp2=0.0
      sigdlc=0.0
      sigprebi=0.0
      zpress=0.0
      tethom=0.0
      rteth=0.0
      zteth=0.0
      sgteth=0.0
      tionex=0.0
      rion=0.0
      zion=0.0
      sigti=0.0
      dnethom=0.0
      zeffvs=0.0
      rneth=0.0
      zneth=0.0
      sgneth=0.0
      pbeam=0.0
      sibeam=0.0
      xalpa=0.0
      cgama=0.0
      calpa=0.0
      xgama=0.0
      bitfc=0.0
      sigtii=0.0
      sgnethi=0.0
      sgtethi=0.0
      zlowimp=0.0
      pressbi=0.0
      prespb=0.0
      sigppb=0.0
      dnbeam=0.0
      dmass=0.0
      sgtimin=0.0
      bitec=0.0
      scalepr=0.0
      brsptu=0.0
      fwtfcsum=0.0
      rsol=0.0
      zsol=0.0
      fbetan=0.0
      fqsiw=0.0
      alpax=0.0
      gamax=0.0
      siwantq=0.0
      ccoils=0.0
      xcoils=0.0
      cloops=0.0
      xloops=0.0
      currc79=0.0
      currc139=0.0
      sizeroj=0.0
      currc199=0.0
      curriu30=0.0
      curriu90=0.0
      curriu150=0.0
      curril30=0.0
      curril90=0.0
      curril150=0.0
      tgamma=0.0
      sgamma=0.0
      rrrgam=0.0
      zzzgam=0.0
      aa1gam=0.0
      aa2gam=0.0
      aa3gam=0.0
      aa4gam=0.0
      aa5gam=0.0
      aa6gam=0.0
      aa7gam=0.0
      tgammauncor=0.0
      bmsels=0.0
      rrmsels=0.0
      zzmsels=0.0
      l1msels=0.0
      l2msels=0.0
      l4msels=0.0
      emsels=0.0
      semsels=0.0
      teecein0=0.0
      feece0=0.0
      errorece0=0.0
      rteo=0.0
      zteo=0.0
      rtep=0.0
      rtem=0.0
      rpbit=0.0
      rmbit=0.0
      robit=0.0
      fwtnow=0.0
      omegat=0.0
      enw=0.0
      emw=0.0
      fwtprw=0.0
      presw=0.0
      sigprw=0.0
      rpresw=0.0
      zpresw=0.0
      comega=0.0
      xomega=0.0
      romegat=0.0
      zomegat=0.0
      sigome=0.0
      scalepw=0.0
! 
      kerror=0 
      idone=0 
      sicont=tmu*drslop/aaslop
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
! 
      if ((kdata.ne.2).and.(kdata.ne.9)) then
!---------------------------------------------------------------------- 
!--   normalize fitting weights, SNAP mode                           -- 
!---------------------------------------------------------------------- 
      if (jtime.le.1) then
        do i=1,nsilop 
          if ((kdata.ge.3).and.(fwtsi(i).ne.0.0)) then 
            fwtsi(i)=swtsi(i) 
            if (lookfw.gt.0) fwtsi(i)=rwtsi(i) 
          endif 
          if (ierpsi(i).ne.0) fwtsi(i)=0.0 
        enddo
        do i=1,nfcoil 
          if ((kdata.ge.3).and.(fwtfc(i).ne.0.0)) fwtfc(i)=swtfc(i) 
          if (ierfc(i).ne.0) fwtfc(i)=0.0 
        enddo
        if (iecurr.eq.2) then 
        do i=1,nesum 
          if ((kdata.ge.3).and.(fwtec(i).ne.0.0)) fwtec(i)=swtec(i) 
          if (ierec(i).ne.0) fwtec(i)=0.0 
        enddo 
        endif 
        do i=1,magpri 
          if ((kdata.ge.3).and.(fwtmp2(i).ne.0.0)) then 
            fwtmp2(i)=swtmp2(i) 
            if (lookfw.gt.0) fwtmp2(i)=rwtmp2(i) 
          endif 
          if (iermpi(i).ne.0) fwtmp2(i)=0.0 
        enddo
        do i=1,nstark 
          fwtgam(i)=swtgam(i) 
          if (iergam(i).ne.0) fwtgam(i)=0.0 
        enddo
        do i=1,nnece 
          fwtece0(i)=swtece(i) 
          if (ierece(i).ne.0) fwtece0(i)=0.0 
        enddo
        fwtecebz0=swtecebz 
        if (ierecebz.ne.0) fwtecebz0=0.0 
        if (fwtcur.ne.0.0) fwtcur=swtcur 
        if (fwtqa.ne.0.0) fwtqa=1. 
        if (fwtbp.ne.0.0) fwtbp=1. 
        if (fwtdlc.ne.0.0) fwtdlc=swtdlc 
        if (ierpla.ne.0) fwtcur=0.0 
        if (ierrdi.ne.0) fwtdlc=0.0 
        !--  Save fitting weights                                           -- 
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
      else
!--------------------------------------------------------------------- 
!--     Restore fitting weights for time slices > 1                 -- 
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
      endif
!----------------------------------------------------------------------- 
!--   Set edge pedestal tanh paramters                                -- 
!----------------------------------------------------------------------- 
      if (fitzts.eq.'te'.and.ztserr(jtime)) then 
        nbdry=1 
        rbdry(1)=1.94_dp 
        zbdry(1)=ztssym(jtime)+0.5_dp*ztswid(jtime) 
      endif 
      else ! (kdata.eq.2).or.(kdata.eq.9)
!---------------------------------------------------------------------- 
!--   file mode - initialize inputs                                  -- 
!---------------------------------------------------------------------- 
      do i=1,nsilop 
        psibit(i)=0.0 
      enddo
      do i=1,magpri 
        bitmpi(i)=0.0 
      enddo
      alpax(jbeta)=1.e4_dp 
      backaverage=.false. 
      bitip=0.0 
      betap0=0.50_dp 
      brsp(1)=-1.e+20_dp 
      cfcoil=-1. 
      cutip=80000. 
      do i=1,nco2v 
        denv(i)=0. 
      enddo
      do i=1,nco2r 
        denr(i)=0. 
      enddo
      do i=1,nesum 
        rsisec(i)=-1. 
      enddo
      emf=1.00 
      emp=1.00 
      enf=1.00 
      enp=1.00 
      error=1.0e-03_dp 
      fbetap=0.0 
      fbetat=0.0 
      do i=1,nfcoil 
        fcsum(i)=1.0 
        fczero(i)=1.0 
        fwtfc(i)=0. 
        rsisfc(i)=-1. 
      enddo
      do i=1,nesum 
        fwtec(i)=0.0 
      enddo 
      do i=1,mpress 
        fwtpre(i)=1. 
      enddo 
      fcurbd=1.0 
      fli=0.0 
      fwtbp=0.0 
      fwtdlc=0.0 
      do i=1,nstark 
        fwtgam(i)=0.0 
      enddo 
      do i=1,nmsels 
        fwtbmsels(i)=0.0 
        fwtemsels(i)=0.0 
      enddo 
      do i=1,nnece 
        fwtece0(i)=0.0 
      enddo 
      fwtecebz0=0.0 
      do i=1,mbdry 
        fwtbdry(i)=1.0 
        sigrbd(i)=1.e10_dp 
        sigzbd(i)=1.e10_dp 
        fwtsol(i)=1.0 
      enddo 
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

      xlim(1)=-1.0 
      rbdry(1)=-1.0 
      itimeu=0 
      nbdryp=-1 
      ktear=0
      rrmsels(1)=-10. 
      ecefit = 0.0 
      ecebzfit = 0.0 

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

      do i=1,libim
        fwtlib(i)=0.0
        rrrlib(i)=0.0
        zzzlib(i)=0.0
      enddo

      if (kdata.eq.2) then 
!---------------------------------------------------------------------- 
!--   Read ascii input files                                         -- 
!---------------------------------------------------------------------- 
      open(unit=nin,status='old',file=ifname(jtime)) 
      read (nin,in1,iostat=istat) 
      if (istat>0) then 
        backspace(nin) 
        read(nin,fmt='(A)') line 
        write(*,'(A)') 'Invalid line in namelist in1: '//trim(line) 
      endif

      read (nin,ink,err=11111,end=101) 
101   continue 
11111 close(unit=nin) 
      open(unit=nin,status='old',file=ifname(jtime)) 
      read (nin,ins,err=11113,end=103) 
103   continue 
11113 close(unit=nin) 
      open(unit=nin,status='old',file=ifname(jtime)) 

      read (nin,in_msels,iostat=istat) 
      if (istat>0) then 
        backspace(nin) 
        read(nin,fmt='(A)') line 
        !if (trim(line)/="/") write(*,'(A)') 'Invalid line in namelist in_msels: '//trim(line) 
      endif
      close(unit=nin) 

      open(unit=nin,status='old',file=ifname(jtime)) 
      read (nin,ina,iostat=istat)
      close(unit=nin) 

      open(unit=nin,status='old',file=ifname(jtime)) 
      read (nin,inece,iostat=istat) 
      close(unit=nin) 

      open(unit=nin,status='old',file=ifname(jtime)) 
      read (nin,edgep,iostat=istat) 
      close(unit=nin) 

      open(unit=nin,status='old',file=ifname(jtime)) 
      read (nin,iner,iostat=istat) 
      close(unit=nin) 

      open(unit=nin,status='old',file=ifname(jtime)) 
      read (nin,insxr,iostat=istat) 
      close(unit=nin) 

      open(unit=nin,status='old',file=ifname(jtime))
      read (nin,inms,iostat=istat) 
      close(unit=nin) 

      open(unit=nin,status='old',file=ifname(jtime)) 
      read (nin,inwant,iostat=istat) 
      close(unit=nin) 

      open(unit=nin,status='old',file=ifname(jtime)) 
      read (nin,invt,iostat=istat) 
      close(unit=nin) 

!--   Input FF', P' arrays
      open(unit=nin,status='old',file=ifname(jtime)) 
      if (idebug /= 0) write (nttyo,*) 'DATA_INPUT geqdsk_ext=',geqdsk_ext 
      read(nin,profile_ext,err=11777,iostat=istat) 
      if (idebug /= 0) write (nttyo,*) 'DATA_INPUT geqdsk_ext=',geqdsk_ext 
!      if (geqdsk_ext.ne.'none') then 
!        open(unit=neqdsk,status='old',file=geqdsk_ext) 
!        read (neqdsk,11775) (case_ext(i),i=1,6),nh_ext,nw_ext,nh_ext 
!        do i = 1,2 
!          read (neqdsk,11773) 
!        enddo 
!        read (neqdsk,11776) plasma_ext,c_ext,c_ext,c_ext,c_ext 
!        read (neqdsk,11773) 
!        read (neqdsk,11776) (ffprim_ext(i),i=1,nw_ext) 
!        read (neqdsk,11776) (pprime_ext(i),i=1,nw_ext) 
!        read (neqdsk,11776) (ffprim_ext(i),i=1,nw_ext) 
!        read (neqdsk,11776) (pprime_ext(i),i=1,nw_ext) 
!        read (neqdsk,11776) ((psirz_ext,i=1,nw_ext),j=1,nh_ext) 
!        read (neqdsk,11776,err=11777) (qpsi_ext(i),i=1,nw_ext) 
!        read (neqdsk,11774,err=11777) nbdry_ext,limitr_ext 
!        read (neqdsk,11776,err=11777) (rbdry_ext(i),zbdry_ext(i),i=1,nbdry_ext) 
!        read (neqdsk,11776,err=11777) (xlim_ext(i),ylim_ext(i),i=1,limitr_ext) 
!11773   format (a) 
!11774   format (2i5) 
!11775   format (6a8,3i4) 
!11776   format (5e16.9)
!      endif
11777 close(nin) 

!--   Read Li beam data 
      open(unit=nin,status='old',file=ifname(jtime))
      read (nin,inlibim,err=11237,end=11233) 
11233 continue 
11237 close(unit=nin) 
!---------------------------------------------------------------------- 
!--   HDF5 file mode                                                 -- 
!---------------------------------------------------------------------- 
      else ! kdata.eq.9
#ifdef HAVE_HDF5
        inquire(file=trim(ifname(1)),exist=file_stat)
        if (.not. file_stat) then
          call errctrl_msg('data_input',trim(line)//' not found')
          stop
        endif
        call fch5init
        call open_oldh5file(trim(ifname(1)),fileid,rootgid,h5in,h5err)
        call test_group(rootgid,"equilibrium",file_stat,h5err)
        if (.not. file_stat) then
          call errctrl_msg('data_input','equilibrium group not found')
          stop
        endif
        call open_group(rootgid,"equilibrium",eqid,h5err)
        call test_group(eqid,"code",file_stat,h5err)
        if (.not. file_stat) then
          call errctrl_msg('data_input','code group not found')
          stop
        endif
        call open_group(eqid,"code",cid,h5err)
        call test_group(cid,"parameters",file_stat,h5err)
        if (.not. file_stat) then
          call errctrl_msg('data_input','parameters group not found')
          stop
        endif
        call open_group(cid,"parameters",pid,h5err)
        call test_group(pid,"time_slice",file_stat,h5err)
        if (.not. file_stat) then
          call errctrl_msg('data_input','time_slice group not found')
          stop
        endif
        call open_group(pid,"time_slice",tid,h5err)
        write(line,"(I0)") jtime-1
        call test_group(tid,trim(line),file_stat,h5err)
        if (.not. file_stat) then
          call errctrl_msg('data_input',trim(line)//' group not found')
          stop
        endif
        call open_group(tid,trim(line),sid,h5err)

        call test_group(sid,"in1",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"in1",nid,h5err)
          if (obj_exists(nid,"ishot",h5err)) &
            call read_h5(nid,"ishot",ishot,h5in,h5err)
          if (obj_exists(nid,"itime",h5err)) &
            call read_h5(nid,"itime",itime,h5in,h5err)
          if (obj_exists(nid,"plasma",h5err)) &
            call read_h5(nid,"plasma",plasma,h5in,h5err)
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
          if (obj_exists(nid,"coils",h5err)) &
            call read_h5(nid,"coils",coils,h5in,h5err)
          if (obj_exists(nid,"fwtsi",h5err)) &
            call read_h5(nid,"fwtsi",fwtsi,h5in,h5err)
          if (obj_exists(nid,"expmp2",h5err)) &
            call read_h5(nid,"expmp2",expmp2,h5in,h5err)
          if (obj_exists(nid,"fwtmp2",h5err)) &
            call read_h5(nid,"fwtmp2",fwtmp2,h5in,h5err)
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
          if (obj_exists(nid,"xlim",h5err)) &
            call read_h5(nid,"xlim",xlim,h5in,h5err)
          if (obj_exists(nid,"ylim",h5err)) &
            call read_h5(nid,"ylim",ylim,h5in,h5err)
          if (obj_exists(nid,"serror",h5err)) &
            call read_h5(nid,"serror",serror,h5in,h5err)
          if (obj_exists(nid,"nbdry",h5err)) &
            call read_h5(nid,"nbdry",nbdry,h5in,h5err)
          ! Not currently in OMAS files!
          rbdry=0.0
          zbdry=0.0
!          if (obj_exists(nid,"rbdry",h5err)) &
!            call read_h5(nid,"rbdry",rbdry,h5in,h5err)
!          if (obj_exists(nid,"zbdry",h5err)) &
!            call read_h5(nid,"zbdry",zbdry,h5in,h5err)
          if (obj_exists(nid,"psibry",h5err)) &
            call read_h5(nid,"psibry",psibry,h5in,h5err)
          if (obj_exists(nid,"nslref",h5err)) &
            call read_h5(nid,"nslref",nslref,h5in,h5err)
          if (obj_exists(nid,"ibunmn",h5err)) &
            call read_h5(nid,"ibunmn",ibunmn,h5in,h5err)
          if (obj_exists(nid,"btor",h5err)) &
            call read_h5(nid,"btor",btor,h5in,h5err)
          if (obj_exists(nid,"psibit",h5err)) &
            call read_h5(nid,"psibit",psibit,h5in,h5err)
          if (obj_exists(nid,"bitmpi",h5err)) &
            call read_h5(nid,"bitmpi",bitmpi,h5in,h5err)
          if (obj_exists(nid,"bitip",h5err)) &
            call read_h5(nid,"bitip",bitip,h5in,h5err)
          if (obj_exists(nid,"icurrt",h5err)) &
            call read_h5(nid,"icurrt",icurrt,h5in,h5err)
          if (obj_exists(nid,"icinit",h5err)) &
            call read_h5(nid,"icinit",icinit,h5in,h5err)
          if (obj_exists(nid,"brsp",h5err)) &
            call read_h5(nid,"brsp",brsp,h5in,h5err)
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
          if (obj_exists(nid,"rzero",h5err)) &
            call read_h5(nid,"rzero",rzero,h5in,h5err)
          if (obj_exists(nid,"gammap",h5err)) &
            call read_h5(nid,"gammap",gammap,h5in,h5err)
          if (obj_exists(nid,"cfcoil",h5err)) &
            call read_h5(nid,"cfcoil",cfcoil,h5in,h5err)
          if (obj_exists(nid,"fczero",h5err)) &
            call read_h5(nid,"fczero",fczero,h5in,h5err)
          if (obj_exists(nid,"fcsum",h5err)) &
            call read_h5(nid,"fcsum",fcsum,h5in,h5err)
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
          if (obj_exists(nid,"ecurrt",h5err)) &
            call read_h5(nid,"ecurrt",ecurrt,h5in,h5err)
          if (obj_exists(nid,"iecoil",h5err)) &
            call read_h5(nid,"iecoil",iecoil,h5in,h5err)
          if (obj_exists(nid,"co2cor",h5err)) &
            call read_h5(nid,"co2cor",co2cor,h5in,h5err)
          if (obj_exists(nid,"vcurrt",h5err)) &
            call read_h5(nid,"vcurrt",vcurrt,h5in,h5err)
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
          if (obj_exists(nid,"pressr",h5err)) &
            call read_h5(nid,"pressr",pressr,h5in,h5err)
          if (obj_exists(nid,"rpress",h5err)) &
            call read_h5(nid,"rpress",rpress,h5in,h5err)
          if (obj_exists(nid,"zpress",h5err)) &
            call read_h5(nid,"zpress",zpress,h5in,h5err)
          if (obj_exists(nid,"sigpre",h5err)) &
            call read_h5(nid,"sigpre",sigpre,h5in,h5err)
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
          if (obj_exists(nid,"pbeam",h5err)) &
            call read_h5(nid,"pbeam",pbeam,h5in,h5err)
          if (obj_exists(nid,"sibeam",h5err)) &
            call read_h5(nid,"sibeam",sibeam,h5in,h5err)
          if (obj_exists(nid,"nbeam",h5err)) &
            call read_h5(nid,"nbeam",nbeam,h5in,h5err)
          if (obj_exists(nid,"rzeroj",h5err)) &
            call read_h5(nid,"rzeroj",rzeroj,h5in,h5err)
          if (obj_exists(nid,"xalpa",h5err)) &
            call read_h5(nid,"xalpa",xalpa,h5in,h5err)
!          DIIID files have 1d arrays here
          if (obj_exists(nid,"cgama",h5err)) then
            call read_ndims(nid,"cgama",dims,h5in,h5err)
            if (dims.eq.2) then
              call read_h5(nid,"cgama",cgama,h5in,h5err)
            elseif (dims.eq.1) then
              call read_h5(nid,"cgama",cgama_flat,h5in,h5err)
              cgama=reshape(cgama_flat,(/nffcur,nffcur/))
            else
              call errctrl_msg('data_input','cgama abnormal dimension')
              stop
            endif
          endif
          if (obj_exists(nid,"ivesel",h5err)) &
            call read_h5(nid,"ivesel",ivesel,h5in,h5err)
          if (obj_exists(nid,"iexcal",h5err)) &
            call read_h5(nid,"iexcal",iexcal,h5in,h5err)
          if (obj_exists(nid,"iconsi",h5err)) &
            call read_h5(nid,"iconsi",iconsi,h5in,h5err)
          if (obj_exists(nid,"fwtfc",h5err)) &
            call read_h5(nid,"fwtfc",fwtfc,h5in,h5err)
          if (obj_exists(nid,"xltype",h5err)) &
            call read_h5(nid,"xltype",xltype,h5in,h5err)
          if (obj_exists(nid,"kcalpa",h5err)) &
            call read_h5(nid,"kcalpa",kcalpa,h5in,h5err)
          if (obj_exists(nid,"kcgama",h5err)) &
            call read_h5(nid,"kcgama",kcgama,h5in,h5err)
!          DIIID files have 1d arrays here
          if (obj_exists(nid,"calpa",h5err)) then
            call read_ndims(nid,"calpa",dims,h5in,h5err)
            if (dims.eq.2) then
              call read_h5(nid,"calpa",calpa,h5in,h5err)
            elseif (dims.eq.1) then
              call read_h5(nid,"calpa",calpa_flat,h5in,h5err)
              calpa=reshape(calpa_flat,(/nppcur,nppcur/))
            else
              call errctrl_msg('data_input','calpa abnormal dimension')
              stop
            endif
          endif
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
          if (obj_exists(nid,"denr",h5err)) &
            call read_h5(nid,"denr",denr,h5in,h5err)
          if (obj_exists(nid,"denv",h5err)) &
            call read_h5(nid,"denv",denv,h5in,h5err)
          if (obj_exists(nid,"xgama",h5err)) &
            call read_h5(nid,"xgama",xgama,h5in,h5err)
          if (obj_exists(nid,"sgnemin",h5err)) &
            call read_h5(nid,"sgnemin",sgnemin,h5in,h5err)
          if (obj_exists(nid,"nptionf",h5err)) &
            call read_h5(nid,"nptionf",nptionf,h5in,h5err)
          if (obj_exists(nid,"currn1",h5err)) &
            call read_h5(nid,"currn1",currn1,h5in,h5err)
          if (obj_exists(nid,"ifitvs",h5err)) &
            call read_h5(nid,"ifitvs",ifitvs,h5in,h5err)
          if (obj_exists(nid,"bitfc",h5err)) &
            call read_h5(nid,"bitfc",bitfc,h5in,h5err)
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
          if (obj_exists(nid,"vcurfb",h5err)) &
            call read_h5(nid,"vcurfb",vcurfb,h5in,h5err)
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
          if (obj_exists(nid,"fwtpre",h5err)) &
            call read_h5(nid,"fwtpre",fwtpre,h5in,h5err)
          if (obj_exists(nid,"ibtcomp",h5err)) &
            call read_h5(nid,"ibtcomp",ibtcomp,h5in,h5err)
          if (obj_exists(nid,"klabel",h5err)) &
            call read_h5(nid,"klabel",klabel,h5in,h5err)
          if (obj_exists(nid,"zmaxvs",h5err)) &
            call read_h5(nid,"zmaxvs",zmaxvs,h5in,h5err)
          if (obj_exists(nid,"dnbeam",h5err)) &
            call read_h5(nid,"dnbeam",dnbeam,h5in,h5err)
          if (obj_exists(nid,"dmass",h5err)) &
            call read_h5(nid,"dmass",dmass,h5in,h5err)
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
          if (obj_exists(nid,"ppknt",h5err)) &
            call read_h5(nid,"ppknt",ppknt,h5in,h5err)
          if (obj_exists(nid,"pptens",h5err)) &
            call read_h5(nid,"pptens",pptens,h5in,h5err)
          if (obj_exists(nid,"kfffnc",h5err)) &
            call read_h5(nid,"kfffnc",kfffnc,h5in,h5err)
          if (obj_exists(nid,"kffknt",h5err)) &
            call read_h5(nid,"kffknt",kffknt,h5in,h5err)
          if (obj_exists(nid,"ffknt",h5err)) &
            call read_h5(nid,"ffknt",ffknt,h5in,h5err)
          if (obj_exists(nid,"fftens",h5err)) &
            call read_h5(nid,"fftens",fftens,h5in,h5err)
          if (obj_exists(nid,"fwtbdry",h5err)) &
            call read_h5(nid,"fwtbdry",fwtbdry,h5in,h5err)
          if (obj_exists(nid,"kwwfnc",h5err)) &
            call read_h5(nid,"kwwfnc",kwwfnc,h5in,h5err)
          if (obj_exists(nid,"kwwknt",h5err)) &
            call read_h5(nid,"kwwknt",kwwknt,h5in,h5err)
          if (obj_exists(nid,"wwknt",h5err)) &
            call read_h5(nid,"wwknt",wwknt,h5in,h5err)
          if (obj_exists(nid,"wwtens",h5err)) &
            call read_h5(nid,"wwtens",wwtens,h5in,h5err)
          if (obj_exists(nid,"fwtec",h5err)) &
            call read_h5(nid,"fwtec",fwtec,h5in,h5err)
          if (obj_exists(nid,"fitsiref",h5err)) &
            call read_h5(nid,"fitsiref",fitsiref_int,h5in,h5err)
          if (fitsiref_int.ne.0) fitsiref=.true.  ! default off
          if (obj_exists(nid,"bitec",h5err)) &
            call read_h5(nid,"bitec",bitec,h5in,h5err)
          if (obj_exists(nid,"scalepr",h5err)) &
            call read_h5(nid,"scalepr",scalepr,h5in,h5err)
          if (obj_exists(nid,"scalesir",h5err)) &
            call read_h5(nid,"scalesir",scalesir,h5in,h5err)
          if (obj_exists(nid,"ppbdry",h5err)) &
            call read_h5(nid,"ppbdry",ppbdry,h5in,h5err)
          if (obj_exists(nid,"kppbdry",h5err)) &
            call read_h5(nid,"kppbdry",kppbdry,h5in,h5err)
          if (obj_exists(nid,"pp2bdry",h5err)) &
            call read_h5(nid,"pp2bdry",pp2bdry,h5in,h5err)
          if (obj_exists(nid,"kpp2bdry",h5err)) &
            call read_h5(nid,"kpp2bdry",kpp2bdry,h5in,h5err)
          if (obj_exists(nid,"scalea",h5err)) &
            call read_h5(nid,"scalea",scalea_int,h5in,h5err)
          if (scalea_int.ne.0) scalea=.true.  ! default off
          if (obj_exists(nid,"sigrbd",h5err)) &
            call read_h5(nid,"sigrbd",sigrbd,h5in,h5err)
          if (obj_exists(nid,"sigzbd",h5err)) &
            call read_h5(nid,"sigzbd",sigzbd,h5in,h5err)
          if (obj_exists(nid,"nbskip",h5err)) &
            call read_h5(nid,"nbskip",nbskip,h5in,h5err)
          if (obj_exists(nid,"ffbdry",h5err)) &
            call read_h5(nid,"ffbdry",ffbdry,h5in,h5err)
          if (obj_exists(nid,"kffbdry",h5err)) &
            call read_h5(nid,"kffbdry",kffbdry,h5in,h5err)
          if (obj_exists(nid,"ff2bdry",h5err)) &
            call read_h5(nid,"ff2bdry",ff2bdry,h5in,h5err)
          if (obj_exists(nid,"kff2bdry",h5err)) &
            call read_h5(nid,"kff2bdry",kff2bdry,h5in,h5err)
          if (obj_exists(nid,"errsil",h5err)) &
            call read_h5(nid,"errsil",errsil,h5in,h5err)
          if (obj_exists(nid,"vbit",h5err)) &
            call read_h5(nid,"vbit",vbit,h5in,h5err)
          if (obj_exists(nid,"wwbdry",h5err)) &
            call read_h5(nid,"wwbdry",wwbdry,h5in,h5err)
          if (obj_exists(nid,"kwwbdry",h5err)) &
            call read_h5(nid,"kwwbdry",kwwbdry,h5in,h5err)
          if (obj_exists(nid,"ww2bdry",h5err)) &
            call read_h5(nid,"ww2bdry",ww2bdry,h5in,h5err)
          if (obj_exists(nid,"kww2bdry",h5err)) &
            call read_h5(nid,"kww2bdry",kww2bdry,h5in,h5err)
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
          if (obj_exists(nid,"fitzts",h5err)) &
            call read_h5(nid,"fitzts",fitzts,h5in,h5err)
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
          if (fitfcsum_int.ne.0) fitfcsum=.true.
          if (obj_exists(nid,"fwtfcsum",h5err)) &
            call read_h5(nid,"fwtfcsum",fwtfcsum,h5in,h5err)
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
          call close_group("in1",nid,h5err)
        endif
   
        call test_group(sid,"ink",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"ink",nid,h5err)
          if (obj_exists(nid,"isetfb",h5err)) &
            call read_h5(nid,"isetfb",isetfb,h5in,h5err)
          if (obj_exists(nid,"ioffr",h5err)) &
            call read_h5(nid,"ioffr",ioffr,h5in,h5err)
          if (obj_exists(nid,"ioffz",h5err)) &
            call read_h5(nid,"ioffz",ioffz,h5in,h5err)
          if (obj_exists(nid,"ishiftz",h5err)) &
            call read_h5(nid,"ishiftz",ishiftz,h5in,h5err)
          if (obj_exists(nid,"gain",h5err)) &
            call read_h5(nid,"gain",gain,h5in,h5err)
          if (obj_exists(nid,"gainp",h5err)) &
            call read_h5(nid,"gainp",gainp,h5in,h5err)
          if (obj_exists(nid,"idplace",h5err)) &
            call read_h5(nid,"idplace",idplace,h5in,h5err)
          if (obj_exists(nid,"symmetrize",h5err)) &
            call read_h5(nid,"symmetrize",symmetrize_int,h5in,h5err)
          if (symmetrize_int.ne.0) symmetrize=.true.  ! default off
          if (obj_exists(nid,"backaverage",h5err)) &
            call read_h5(nid,"backaverage",backaverage_int,h5in,h5err)
          if (backaverage_int.ne.0) backaverage=.true.  ! default off
          if (obj_exists(nid,"lring",h5err)) &
            call read_h5(nid,"lring",lring,h5in,h5err)
          if (obj_exists(nid,"cupdown",h5err)) &
            call read_h5(nid,"cupdown",cupdown,h5in,h5err)
          call close_group("ink",nid,h5err)
        endif
   
        call test_group(sid,"ins",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"ins",nid,h5err)
          if (obj_exists(nid,"tgamma",h5err)) &
            call read_h5(nid,"tgamma",tgamma,h5in,h5err)
          if (obj_exists(nid,"sgamma",h5err)) &
            call read_h5(nid,"sgamma",sgamma,h5in,h5err)
          if (obj_exists(nid,"fwtgam",h5err)) &
            call read_h5(nid,"fwtgam",fwtgam,h5in,h5err)
          if (obj_exists(nid,"rrrgam",h5err)) &
            call read_h5(nid,"rrrgam",rrrgam,h5in,h5err)
          if (obj_exists(nid,"zzzgam",h5err)) &
            call read_h5(nid,"zzzgam",zzzgam,h5in,h5err)
          if (obj_exists(nid,"aa1gam",h5err)) &
            call read_h5(nid,"aa1gam",aa1gam,h5in,h5err)
          if (obj_exists(nid,"aa2gam",h5err)) &
            call read_h5(nid,"aa2gam",aa2gam,h5in,h5err)
          if (obj_exists(nid,"aa3gam",h5err)) &
            call read_h5(nid,"aa3gam",aa3gam,h5in,h5err)
          if (obj_exists(nid,"aa4gam",h5err)) &
            call read_h5(nid,"aa4gam",aa4gam,h5in,h5err)
          if (obj_exists(nid,"aa5gam",h5err)) &
            call read_h5(nid,"aa5gam",aa5gam,h5in,h5err)
          if (obj_exists(nid,"aa6gam",h5err)) &
            call read_h5(nid,"aa6gam",aa6gam,h5in,h5err)
          if (obj_exists(nid,"aa7gam",h5err)) &
            call read_h5(nid,"aa7gam",aa7gam,h5in,h5err)
          if (obj_exists(nid,"iplots",h5err)) &
            call read_h5(nid,"iplots",iplots,h5in,h5err)
          if (obj_exists(nid,"kdomse",h5err)) &
            call read_h5(nid,"kdomse",kdomse,h5in,h5err)
          if (obj_exists(nid,"msebkp",h5err)) &
            call read_h5(nid,"msebkp",msebkp,h5in,h5err)
          if (obj_exists(nid,"msefitfun",h5err)) &
            call read_h5(nid,"msefitfun",msefitfun,h5in,h5err)
          if (obj_exists(nid,"mse_quiet",h5err)) &
            call read_h5(nid,"mse_quiet",mse_quiet,h5in,h5err)
          if (obj_exists(nid,"mse_spave_on",h5err)) &
            call read_h5(nid,"mse_spave_on",mse_spave_on,h5in,h5err)
          if (obj_exists(nid,"kwaitmse",h5err)) &
            call read_h5(nid,"kwaitmse",kwaitmse,h5in,h5err)
          if (obj_exists(nid,"dtmsefull",h5err)) &
            call read_h5(nid,"dtmsefull",dtmsefull,h5in,h5err)
          if (obj_exists(nid,"mse_strict",h5err)) &
            call read_h5(nid,"mse_strict",mse_strict,h5in,h5err)
          if (obj_exists(nid,"t_max_beam_off",h5err)) &
            call read_h5(nid,"t_max_beam_off",t_max_beam_off,h5in,h5err)
          if (obj_exists(nid,"ok_30rt",h5err)) &
            call read_h5(nid,"ok_30rt",ok_30rt,h5in,h5err)
          if (obj_exists(nid,"ok_210lt",h5err)) &
            call read_h5(nid,"ok_210lt",ok_210lt,h5in,h5err)
          if (obj_exists(nid,"mse_usecer",h5err)) &
            call read_h5(nid,"mse_usecer",mse_usecer,h5in,h5err)
          if (obj_exists(nid,"mse_certree",h5err)) &
            call read_h5(nid,"mse_certree",mse_certree,h5in,h5err)
          if (obj_exists(nid,"mse_use_cer330",h5err)) &
            call read_h5(nid,"mse_use_cer330",mse_use_cer330,h5in,h5err)
          if (obj_exists(nid,"mse_use_cer210",h5err)) &
            call read_h5(nid,"mse_use_cer210",mse_use_cer210,h5in,h5err)
          if (obj_exists(nid,"tgammauncor",h5err)) &
            call read_h5(nid,"tgammauncor",tgammauncor,h5in,h5err)
          if (obj_exists(nid,"v30lt",h5err)) &
            call read_h5(nid,"v30lt",v30lt,h5in,h5err)
          if (obj_exists(nid,"v30rt",h5err)) &
            call read_h5(nid,"v30rt",v30rt,h5in,h5err)
          if (obj_exists(nid,"v210lt",h5err)) &
            call read_h5(nid,"v210lt",v210lt,h5in,h5err)
          if (obj_exists(nid,"v210rt",h5err)) &
            call read_h5(nid,"v210rt",v210rt,h5in,h5err) 
          call close_group("ins",nid,h5err)
        endif
   
        call test_group(sid,"in_msels",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"in_msels",nid,h5err)
          if (obj_exists(nid,"bmsels",h5err)) &
            call read_h5(nid,"bmsels",bmsels,h5in,h5err)
          if (obj_exists(nid,"sbmsels",h5err)) &
            call read_h5(nid,"sbmsels",sbmsels,h5in,h5err)
          if (obj_exists(nid,"fwtbmsels",h5err)) &
            call read_h5(nid,"fwtbmsels",fwtbmsels,h5in,h5err)
          if (obj_exists(nid,"rrmsels",h5err)) &
            call read_h5(nid,"rrmsels",rrmsels,h5in,h5err)
          if (obj_exists(nid,"zzmsels",h5err)) &
            call read_h5(nid,"zzmsels",zzmsels,h5in,h5err)
          if (obj_exists(nid,"l1msels",h5err)) &
            call read_h5(nid,"l1msels",l1msels,h5in,h5err)
          if (obj_exists(nid,"l2msels",h5err)) &
            call read_h5(nid,"l2msels",l2msels,h5in,h5err)
          if (obj_exists(nid,"l4msels",h5err)) &
            call read_h5(nid,"l4msels",l4msels,h5in,h5err)
          if (obj_exists(nid,"emsels",h5err)) &
            call read_h5(nid,"emsels",emsels,h5in,h5err)
          if (obj_exists(nid,"semsels",h5err)) &
            call read_h5(nid,"semsels",semsels,h5in,h5err)
          if (obj_exists(nid,"fwtemsels",h5err)) &
            call read_h5(nid,"fwtemsels",fwtemsels,h5in,h5err)
          if (obj_exists(nid,"kdomsels",h5err)) &
            call read_h5(nid,"kdomsels",kdomsels,h5in,h5err)
          if (obj_exists(nid,"fmlscut",h5err)) &
            call read_h5(nid,"fmlscut",fmlscut,h5in,h5err)
          if (obj_exists(nid,"synmsels",h5err)) &
            call read_h5(nid,"synmsels",synmsels,h5in,h5err)
          if (obj_exists(nid,"avemsels",h5err)) &
            call read_h5(nid,"avemsels",avemsels,h5in,h5err)
          call close_group("in_msels",nid,h5err)
        endif

        call test_group(sid,"ina",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"ina",nid,h5err)
          if (obj_exists(nid,"spatial_avg_gam",h5err)) &
            call read_h5(nid,"spatial_avg_gam",spatial_avg_gam,h5in,h5err)
          call close_group("ina",nid,h5err)
        endif

        call test_group(sid,"inece",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"inece",nid,h5err)
          if (obj_exists(nid,"necein",h5err)) &
            call read_h5(nid,"necein",necein,h5in,h5err)
          if (obj_exists(nid,"teecein0",h5err)) &
            call read_h5(nid,"teecein0",teecein0,h5in,h5err)
          if (obj_exists(nid,"feece0",h5err)) &
            call read_h5(nid,"feece0",feece0,h5in,h5err)
          if (obj_exists(nid,"errorece0",h5err)) &
            call read_h5(nid,"errorece0",errorece0,h5in,h5err)
          if (obj_exists(nid,"fwtece0",h5err)) &
            call read_h5(nid,"fwtece0",fwtece0,h5in,h5err)
          if (obj_exists(nid,"fwtecebz0",h5err)) &
            call read_h5(nid,"fwtecebz0",fwtecebz0,h5in,h5err)
          if (obj_exists(nid,"ecefit",h5err)) &
            call read_h5(nid,"ecefit",ecefit,h5in,h5err)
          if (obj_exists(nid,"ecebzfit",h5err)) &
            call read_h5(nid,"ecebzfit",ecebzfit,h5in,h5err)
          if (obj_exists(nid,"kfitece",h5err)) &
            call read_h5(nid,"kfitece",kfitece,h5in,h5err)
          if (obj_exists(nid,"kinputece",h5err)) &
            call read_h5(nid,"kinputece",kinputece,h5in,h5err)
          if (obj_exists(nid,"kcallece",h5err)) &
            call read_h5(nid,"kcallece",kcallece,h5in,h5err)
          if (obj_exists(nid,"nharm",h5err)) &
            call read_h5(nid,"nharm",nharm,h5in,h5err)
          if (obj_exists(nid,"kfixro",h5err)) &
            call read_h5(nid,"kfixro",kfixro,h5in,h5err)
          if (obj_exists(nid,"rteo",h5err)) &
            call read_h5(nid,"rteo",rteo,h5in,h5err)
          if (obj_exists(nid,"zteo",h5err)) &
            call read_h5(nid,"zteo",zteo,h5in,h5err)
          if (obj_exists(nid,"kfixrece",h5err)) &
            call read_h5(nid,"kfixrece",kfixrece,h5in,h5err)
          if (obj_exists(nid,"rtep",h5err)) &
            call read_h5(nid,"rtep",rtep,h5in,h5err)
          if (obj_exists(nid,"rtem",h5err)) &
            call read_h5(nid,"rtem",rtem,h5in,h5err)
          if (obj_exists(nid,"rpbit",h5err)) &
            call read_h5(nid,"rpbit",rpbit,h5in,h5err)
          if (obj_exists(nid,"rmbit",h5err)) &
            call read_h5(nid,"rmbit",rmbit,h5in,h5err)
          if (obj_exists(nid,"robit",h5err)) &
            call read_h5(nid,"robit",robit,h5in,h5err)
          if (obj_exists(nid,"nfit",h5err)) &
            call read_h5(nid,"nfit",nfit,h5in,h5err)
          if (obj_exists(nid,"kcmin",h5err)) &
            call read_h5(nid,"kcmin",kcmin,h5in,h5err)
          if (obj_exists(nid,"fwtnow",h5err)) &
            call read_h5(nid,"fwtnow",fwtnow,h5in,h5err)
          if (obj_exists(nid,"kdoece",h5err)) &
            call read_h5(nid,"kdoece",kdoece,h5in,h5err)
          if (obj_exists(nid,"mtxece",h5err)) &
            call read_h5(nid,"mtxece",mtxece,h5in,h5err)
          if (obj_exists(nid,"nconstr",h5err)) &
            call read_h5(nid,"nconstr",nconstr,h5in,h5err)
          if (obj_exists(nid,"eceiter",h5err)) &
            call read_h5(nid,"eceiter",eceiter,h5in,h5err)
          if (obj_exists(nid,"eceerror",h5err)) &
            call read_h5(nid,"eceerror",eceerror,h5in,h5err)
          call close_group("inece",nid,h5err)
        endif

        call test_group(sid,"edgep",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"edgep",nid,h5err)
          if (obj_exists(nid,"symmetrize",h5err)) &
            call read_h5(nid,"symmetrize",symmetrize_int,h5in,h5err)
          if (symmetrize_int.ne.0) symmetrize=.true.  ! default off
          if (obj_exists(nid,"rpress",h5err)) &
            call read_h5(nid,"rpress",rpress,h5in,h5err)
          if (obj_exists(nid,"pressr",h5err)) &
            call read_h5(nid,"pressr",pressr,h5in,h5err)
          if (obj_exists(nid,"sigpre",h5err)) &
            call read_h5(nid,"sigpre",sigpre,h5in,h5err)
          if (obj_exists(nid,"npress",h5err)) &
            call read_h5(nid,"npress",npress,h5in,h5err)
          if (obj_exists(nid,"kprfit",h5err)) &
            call read_h5(nid,"kprfit",kprfit,h5in,h5err)
          if (obj_exists(nid,"kpressb",h5err)) &
            call read_h5(nid,"kpressb",kpressb,h5in,h5err)
          if (obj_exists(nid,"ndokin",h5err)) &
            call read_h5(nid,"ndokin",ndokin,h5in,h5err)
          if (obj_exists(nid,"kppfnc",h5err)) &
            call read_h5(nid,"kppfnc",kppfnc,h5in,h5err)
          if (obj_exists(nid,"kfffnc",h5err)) &
            call read_h5(nid,"kfffnc",kfffnc,h5in,h5err)
          if (obj_exists(nid,"kffcur",h5err)) &
            call read_h5(nid,"kffcur",kffcur,h5in,h5err)
          if (obj_exists(nid,"kppcur",h5err)) &
            call read_h5(nid,"kppcur",kppcur,h5in,h5err)
          if (obj_exists(nid,"mxiter",h5err)) &
            call read_h5(nid,"mxiter",mxiter,h5in,h5err)
          if (obj_exists(nid,"error",h5err)) &
            call read_h5(nid,"error",error,h5in,h5err)
          if (obj_exists(nid,"errmin",h5err)) &
            call read_h5(nid,"errmin",errmin,h5in,h5err)
          if (obj_exists(nid,"keecur",h5err)) &
            call read_h5(nid,"keecur",keecur,h5in,h5err)
          call close_group("edgep",nid,h5err)
        endif

        call test_group(sid,"iner",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"iner",nid,h5err)
          if (obj_exists(nid,"keecur",h5err)) &
            call read_h5(nid,"keecur",keecur,h5in,h5err)
          if (obj_exists(nid,"ecurbd",h5err)) &
            call read_h5(nid,"ecurbd",ecurbd,h5in,h5err)
          if (obj_exists(nid,"keefnc",h5err)) &
            call read_h5(nid,"keefnc",keefnc,h5in,h5err)
          if (obj_exists(nid,"eetens",h5err)) &
            call read_h5(nid,"eetens",eetens,h5in,h5err)
          if (obj_exists(nid,"keebdry",h5err)) &
            call read_h5(nid,"keebdry",keebdry,h5in,h5err)
          if (obj_exists(nid,"kee2bdry",h5err)) &
            call read_h5(nid,"kee2bdry",kee2bdry,h5in,h5err)
          if (obj_exists(nid,"eebdry",h5err)) &
            call read_h5(nid,"eebdry",eebdry,h5in,h5err)
          if (obj_exists(nid,"ee2bdry",h5err)) &
            call read_h5(nid,"ee2bdry",ee2bdry,h5in,h5err)
          if (obj_exists(nid,"eeknt",h5err)) &
            call read_h5(nid,"eeknt",eeknt,h5in,h5err)
          if (obj_exists(nid,"keeknt",h5err)) &
            call read_h5(nid,"keeknt",keeknt,h5in,h5err)
          if (obj_exists(nid,"keehord",h5err)) &
            call read_h5(nid,"keehord",keehord,h5in,h5err)
          call close_group("iner",nid,h5err)
        endif

        call test_group(sid,"insxr",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"insxr",nid,h5err)
          if (obj_exists(nid,"ksxr0",h5err)) &
            call read_h5(nid,"ksxr0",ksxr0,h5in,h5err)
          if (obj_exists(nid,"ksxr2",h5err)) &
            call read_h5(nid,"ksxr2",ksxr2,h5in,h5err)
          if (obj_exists(nid,"idosxr",h5err)) &
            call read_h5(nid,"idosxr",idosxr,h5in,h5err)
          call close_group("insxr",nid,h5err)
        endif

        call test_group(sid,"inms",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"inms",nid,h5err)
          if (obj_exists(nid,"xmprcg",h5err)) &
            call read_h5(nid,"xmprcg",xmprcg,h5in,h5err)
          if (obj_exists(nid,"xmp_k",h5err)) &
            call read_h5(nid,"xmp_k",xmp_k,h5in,h5err)
          if (obj_exists(nid,"vresxmp",h5err)) &
            call read_h5(nid,"vresxmp",vresxmp,h5in,h5err)
          if (obj_exists(nid,"t0xmp",h5err)) &
            call read_h5(nid,"t0xmp",t0xmp,h5in,h5err)
          if (obj_exists(nid,"psircg",h5err)) &
            call read_h5(nid,"psircg",psircg,h5in,h5err)
          if (obj_exists(nid,"psi_k",h5err)) &
            call read_h5(nid,"psi_k",psi_k,h5in,h5err)
          if (obj_exists(nid,"vrespsi",h5err)) &
            call read_h5(nid,"vrespsi",vrespsi,h5in,h5err)
          if (obj_exists(nid,"t0psi",h5err)) &
            call read_h5(nid,"t0psi",t0psi,h5in,h5err)
          if (obj_exists(nid,"fcrcg",h5err)) &
            call read_h5(nid,"fcrcg",fcrcg,h5in,h5err)
          if (obj_exists(nid,"fc_k",h5err)) &
            call read_h5(nid,"fc_k",fc_k,h5in,h5err)
          if (obj_exists(nid,"vresfc",h5err)) &
            call read_h5(nid,"vresfc",vresfc,h5in,h5err)
          if (obj_exists(nid,"t0fc",h5err)) &
            call read_h5(nid,"t0fc",t0fc,h5in,h5err)
          if (obj_exists(nid,"ercg",h5err)) &
            call read_h5(nid,"ercg",ercg,h5in,h5err)
          if (obj_exists(nid,"e_k",h5err)) &
            call read_h5(nid,"e_k",e_k,h5in,h5err)
          if (obj_exists(nid,"vrese",h5err)) &
            call read_h5(nid,"vrese",vrese,h5in,h5err)
          if (obj_exists(nid,"t0e",h5err)) &
            call read_h5(nid,"t0e",t0e,h5in,h5err)
          if (obj_exists(nid,"bcrcg",h5err)) &
            call read_h5(nid,"bcrcg",bcrcg,h5in,h5err)
          if (obj_exists(nid,"bc_k",h5err)) &
            call read_h5(nid,"bc_k",bc_k,h5in,h5err)
          if (obj_exists(nid,"vresbc",h5err)) &
            call read_h5(nid,"vresbc",vresbc,h5in,h5err)
          if (obj_exists(nid,"t0bc",h5err)) &
            call read_h5(nid,"t0bc",t0bc,h5in,h5err)
          if (obj_exists(nid,"prcg",h5err)) &
            call read_h5(nid,"prcg",prcg,h5in,h5err)
          if (obj_exists(nid,"p_k",h5err)) &
            call read_h5(nid,"p_k",p_k,h5in,h5err)
          if (obj_exists(nid,"vresp",h5err)) &
            call read_h5(nid,"vresp",vresp,h5in,h5err)
          if (obj_exists(nid,"t0p",h5err)) &
            call read_h5(nid,"t0p",t0p,h5in,h5err)
          if (obj_exists(nid,"bti322in",h5err)) &
            call read_h5(nid,"bti322in",bti322in,h5in,h5err)
          if (obj_exists(nid,"curc79in",h5err)) &
            call read_h5(nid,"curc79in",curc79in,h5in,h5err)
          if (obj_exists(nid,"curc139in",h5err)) &
            call read_h5(nid,"curc139in",curc139in,h5in,h5err)
          if (obj_exists(nid,"curc199in",h5err)) &
            call read_h5(nid,"curc199in",curc199in,h5in,h5err)
          if (obj_exists(nid,"devxmpin",h5err)) &
            call read_h5(nid,"devxmpin",devxmpin,h5in,h5err)
          if (obj_exists(nid,"rnavxmpin",h5err)) &
            call read_h5(nid,"rnavxmpin",rnavxmpin,h5in,h5err)
          if (obj_exists(nid,"devpsii",h5err)) &
            call read_h5(nid,"devpsii",devpsii,h5in,h5err)
          if (obj_exists(nid,"rnavpsiin",h5err)) &
            call read_h5(nid,"rnavpsiin",rnavpsiin,h5in,h5err)
          if (obj_exists(nid,"devfcin",h5err)) &
            call read_h5(nid,"devfcin",devfcin,h5in,h5err)
          if (obj_exists(nid,"rnavfcin",h5err)) &
            call read_h5(nid,"rnavfcin",rnavfcin,h5in,h5err)
          if (obj_exists(nid,"devein",h5err)) &
            call read_h5(nid,"devein",devein,h5in,h5err)
          if (obj_exists(nid,"rnavecin",h5err)) &
            call read_h5(nid,"rnavecin",rnavecin,h5in,h5err)
          call close_group("inms",nid,h5err)
        endif
   
        call test_group(sid,"inwant",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"inwant",nid,h5err)
          if (obj_exists(nid,"psiwant",h5err)) &
            call read_h5(nid,"psiwant",psiwant,h5in,h5err)
          if (obj_exists(nid,"vzeroj",h5err)) &
            call read_h5(nid,"vzeroj",vzeroj,h5in,h5err)
          if (obj_exists(nid,"fwtxxj",h5err)) &
            call read_h5(nid,"fwtxxj",fwtxxj,h5in,h5err)
          if (obj_exists(nid,"fbetap",h5err)) &
            call read_h5(nid,"fbetap",fbetap,h5in,h5err)
          if (obj_exists(nid,"fbetan",h5err)) &
            call read_h5(nid,"fbetan",fbetan,h5in,h5err)
          if (obj_exists(nid,"fli",h5err)) &
            call read_h5(nid,"fli",fli,h5in,h5err)
          if (obj_exists(nid,"fq95",h5err)) &
            call read_h5(nid,"fq95",fq95,h5in,h5err)
          if (obj_exists(nid,"fqsiw",h5err)) &
            call read_h5(nid,"fqsiw",fqsiw,h5in,h5err)
          if (obj_exists(nid,"jbeta",h5err)) &
            call read_h5(nid,"jbeta",jbeta,h5in,h5err)
          if (obj_exists(nid,"jli",h5err)) &
            call read_h5(nid,"jli",jli,h5in,h5err)
          if (obj_exists(nid,"alpax",h5err)) &
            call read_h5(nid,"alpax",alpax,h5in,h5err)
          if (obj_exists(nid,"gamax",h5err)) &
            call read_h5(nid,"gamax",gamax,h5in,h5err)
          if (obj_exists(nid,"jwantm",h5err)) &
            call read_h5(nid,"jwantm",jwantm,h5in,h5err)
          if (obj_exists(nid,"fwtxxq",h5err)) &
            call read_h5(nid,"fwtxxq",fwtxxq,h5in,h5err)
          if (obj_exists(nid,"fwtxxb",h5err)) &
            call read_h5(nid,"fwtxxb",fwtxxb,h5in,h5err)
          if (obj_exists(nid,"fwtxli",h5err)) &
            call read_h5(nid,"fwtxli",fwtxli,h5in,h5err)
          if (obj_exists(nid,"znose",h5err)) &
            call read_h5(nid,"znose",znose,h5in,h5err)
          if (obj_exists(nid,"fwtbdry",h5err)) &
            call read_h5(nid,"fwtbdry",fwtbdry,h5in,h5err)
          if (obj_exists(nid,"nqwant",h5err)) &
            call read_h5(nid,"nqwant",nqwant,h5in,h5err)
          if (obj_exists(nid,"siwantq",h5err)) &
            call read_h5(nid,"siwantq",siwantq,h5in,h5err)
          if (obj_exists(nid,"n_write",h5err)) &
            call read_h5(nid,"n_write",n_write,h5in,h5err)
          if (obj_exists(nid,"kccoils",h5err)) &
            call read_h5(nid,"kccoils",kccoils,h5in,h5err)
          if (obj_exists(nid,"ccoils",h5err)) &
            call read_h5(nid,"ccoils",ccoils,h5in,h5err)
          if (obj_exists(nid,"rexpan",h5err)) &
            call read_h5(nid,"rexpan",rexpan,h5in,h5err)
          if (obj_exists(nid,"xcoils",h5err)) &
            call read_h5(nid,"xcoils",xcoils,h5in,h5err)
          if (obj_exists(nid,"kcloops",h5err)) &
            call read_h5(nid,"kcloops",kcloops,h5in,h5err)
          if (obj_exists(nid,"cloops",h5err)) &
            call read_h5(nid,"cloops",cloops,h5in,h5err)
          if (obj_exists(nid,"xloops",h5err)) &
            call read_h5(nid,"xloops",xloops,h5in,h5err)
          if (obj_exists(nid,"currc79",h5err)) &
            call read_h5(nid,"currc79",currc79,h5in,h5err)
          if (obj_exists(nid,"currc139",h5err)) &
            call read_h5(nid,"currc139",currc139,h5in,h5err)
          if (obj_exists(nid,"nccoil",h5err)) &
            call read_h5(nid,"nccoil",nccoil,h5in,h5err)
          if (obj_exists(nid,"sizeroj",h5err)) &
            call read_h5(nid,"sizeroj",sizeroj,h5in,h5err)
          if (obj_exists(nid,"fitdelz",h5err)) &
            call read_h5(nid,"fitdelz",fitdelz_int,h5in,h5err)
          if (fitdelz_int.ne.0) fitdelz=.true.  ! default off
          if (obj_exists(nid,"ndelzon",h5err)) &
            call read_h5(nid,"ndelzon",ndelzon,h5in,h5err)
          if (obj_exists(nid,"relaxdz",h5err)) &
            call read_h5(nid,"relaxdz",relaxdz,h5in,h5err)
          if (obj_exists(nid,"stabdz",h5err)) &
            call read_h5(nid,"stabdz",stabdz,h5in,h5err)
          if (obj_exists(nid,"writepc",h5err)) &
            call read_h5(nid,"writepc",writepc_int,h5in,h5err)
          if (writepc_int.ne.0) writepc=.true.  ! default off
          if (obj_exists(nid,"table_dir",h5err)) &
            call read_h5(nid,"table_dir",table_dir,h5in,h5err)
          if (obj_exists(nid,"errdelz",h5err)) &
            call read_h5(nid,"errdelz",errdelz,h5in,h5err)
          if (obj_exists(nid,"oldccomp",h5err)) &
            call read_h5(nid,"oldccomp",oldccomp_int,h5in,h5err)
          if (oldccomp_int.ne.0) oldccomp=.true.  ! default off
          if (obj_exists(nid,"nicoil",h5err)) &
            call read_h5(nid,"nicoil",nicoil,h5in,h5err)
          if (obj_exists(nid,"oldcomp",h5err)) &
            call read_h5(nid,"oldcomp",oldcomp_int,h5in,h5err)
          if (oldcomp_int.ne.0) oldcomp=.true.  ! default off
          if (obj_exists(nid,"currc199",h5err)) &
            call read_h5(nid,"currc199",currc199,h5in,h5err)
          if (obj_exists(nid,"curriu30",h5err)) &
            call read_h5(nid,"curriu30",curriu30,h5in,h5err)
          if (obj_exists(nid,"curriu90",h5err)) &
            call read_h5(nid,"curriu90",curriu90,h5in,h5err)
          if (obj_exists(nid,"curriu150",h5err)) &
            call read_h5(nid,"curriu150",curriu150,h5in,h5err)
          if (obj_exists(nid,"curril30",h5err)) &
            call read_h5(nid,"curril30",curril30,h5in,h5err)
          if (obj_exists(nid,"curril90",h5err)) &
            call read_h5(nid,"curril90",curril90,h5in,h5err)
          if (obj_exists(nid,"curril150",h5err)) &
            call read_h5(nid,"curril150",curril150,h5in,h5err)
          if (obj_exists(nid,"ifitdelz",h5err)) &
            call read_h5(nid,"ifitdelz",ifitdelz,h5in,h5err)
          if (obj_exists(nid,"scaledz",h5err)) &
            call read_h5(nid,"scaledz",scaledz,h5in,h5err)
          call close_group("inwant",nid,h5err)
        endif
   
        call test_group(sid,"invt",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"invt",nid,h5err)
          if (obj_exists(nid,"omegat",h5err)) &
            call read_h5(nid,"omegat",omegat,h5in,h5err)
          if (obj_exists(nid,"nomegat",h5err)) &
            call read_h5(nid,"nomegat",nomegat,h5in,h5err)
          if (obj_exists(nid,"enw",h5err)) &
            call read_h5(nid,"enw",enw,h5in,h5err)
          if (obj_exists(nid,"emw",h5err)) &
            call read_h5(nid,"emw",emw,h5in,h5err)
          if (obj_exists(nid,"betapw0",h5err)) &
            call read_h5(nid,"betapw0",betapw0,h5in,h5err)
          if (obj_exists(nid,"kvtor",h5err)) &
            call read_h5(nid,"kvtor",kvtor,h5in,h5err)
          if (obj_exists(nid,"kwwcur",h5err)) &
            call read_h5(nid,"kwwcur",kwwcur,h5in,h5err)
          if (obj_exists(nid,"rvtor",h5err)) &
            call read_h5(nid,"rvtor",rvtor,h5in,h5err)
          if (obj_exists(nid,"wcurbd",h5err)) &
            call read_h5(nid,"wcurbd",wcurbd,h5in,h5err)
          if (obj_exists(nid,"preswb",h5err)) &
            call read_h5(nid,"preswb",preswb,h5in,h5err)
          if (obj_exists(nid,"fwtprw",h5err)) &
            call read_h5(nid,"fwtprw",fwtprw,h5in,h5err)
          if (obj_exists(nid,"npresw",h5err)) &
            call read_h5(nid,"npresw",npresw,h5in,h5err)
          if (obj_exists(nid,"presw",h5err)) &
            call read_h5(nid,"presw",presw,h5in,h5err)
          if (obj_exists(nid,"sigprw",h5err)) &
            call read_h5(nid,"sigprw",sigprw,h5in,h5err)
          if (obj_exists(nid,"rpresw",h5err)) &
            call read_h5(nid,"rpresw",rpresw,h5in,h5err)
          if (obj_exists(nid,"zpresw",h5err)) &
            call read_h5(nid,"zpresw",zpresw,h5in,h5err)
          if (obj_exists(nid,"kplotp",h5err)) &
            call read_h5(nid,"kplotp",kplotp,h5in,h5err)
          if (obj_exists(nid,"sbetaw",h5err)) &
            call read_h5(nid,"sbetaw",sbetaw,h5in,h5err)
          if (obj_exists(nid,"nsplot",h5err)) &
            call read_h5(nid,"nsplot",nsplot,h5in,h5err)
          if (obj_exists(nid,"comega",h5err)) &
            call read_h5(nid,"comega",comega,h5in,h5err)
          if (obj_exists(nid,"kcomega",h5err)) &
            call read_h5(nid,"kcomega",kcomega,h5in,h5err)
          if (obj_exists(nid,"xomega",h5err)) &
            call read_h5(nid,"xomega",xomega,h5in,h5err)
          if (obj_exists(nid,"kdovt",h5err)) &
            call read_h5(nid,"kdovt",kdovt,h5in,h5err)
          if (obj_exists(nid,"romegat",h5err)) &
            call read_h5(nid,"romegat",romegat,h5in,h5err)
          if (obj_exists(nid,"zomegat",h5err)) &
            call read_h5(nid,"zomegat",zomegat,h5in,h5err)
          if (obj_exists(nid,"sigome",h5err)) &
            call read_h5(nid,"sigome",sigome,h5in,h5err)
          if (obj_exists(nid,"scalep",h5err)) &
            call read_h5(nid,"scalep",scalep,h5in,h5err)
          if (obj_exists(nid,"kwwfnc",h5err)) &
            call read_h5(nid,"kwwfnc",kwwfnc,h5in,h5err)
          if (obj_exists(nid,"kwwknt",h5err)) &
            call read_h5(nid,"kwwknt",kwwknt,h5in,h5err)
          if (obj_exists(nid,"wwknt",h5err)) &
            call read_h5(nid,"wwknt",wwknt,h5in,h5err)
          if (obj_exists(nid,"wwtens",h5err)) &
            call read_h5(nid,"wwtens",wwtens,h5in,h5err)
          call close_group("invt",nid,h5err)
        endif

        call test_group(sid,"profile_ext",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"profile_ext",nid,h5err)
          if (obj_exists(nid,"npsi_ext",h5err)) &
            call read_h5(nid,"npsi_ext",npsi_ext,h5in,h5err)
          if (obj_exists(nid,"pprime_ext",h5err)) &
            call read_h5(nid,"pprime_ext",pprime_ext,h5in,h5err)
          if (obj_exists(nid,"ffprim_ext",h5err)) &
            call read_h5(nid,"ffprim_ext",ffprim_ext,h5in,h5err)
          if (obj_exists(nid,"psin_ext",h5err)) &
            call read_h5(nid,"psin_ext",psin_ext,h5in,h5err)
          if (obj_exists(nid,"geqdsk_ext",h5err)) &
            call read_h5(nid,"geqdsk_ext",geqdsk_ext,h5in,h5err)
          if (obj_exists(nid,"sign_ext",h5err)) &
            call read_h5(nid,"sign_ext",sign_ext,h5in,h5err)
          if (obj_exists(nid,"scalepp_ext",h5err)) &
            call read_h5(nid,"scalepp_ext",scalepp_ext,h5in,h5err)
          if (obj_exists(nid,"scaleffp_ext",h5err)) &
            call read_h5(nid,"scaleffp_ext",scaleffp_ext,h5in,h5err)
          if (obj_exists(nid,"shape_ext",h5err)) &
            call read_h5(nid,"shape_ext",shape_ext_int,h5in,h5err)
          if (shape_ext_int.ne.0) shape_ext=.true.  ! default off
          if (obj_exists(nid,"dr_ext",h5err)) &
            call read_h5(nid,"dr_ext",dr_ext,h5in,h5err)
          if (obj_exists(nid,"dz_ext",h5err)) &
            call read_h5(nid,"dz_ext",dz_ext,h5in,h5err)
          if (obj_exists(nid,"rc_ext",h5err)) &
            call read_h5(nid,"rc_ext",rc_ext,h5in,h5err)
          if (obj_exists(nid,"zc_ext",h5err)) &
            call read_h5(nid,"zc_ext",zc_ext,h5in,h5err)
          if (obj_exists(nid,"a_ext",h5err)) &
            call read_h5(nid,"a_ext",a_ext,h5in,h5err)
          if (obj_exists(nid,"eup_ext",h5err)) &
            call read_h5(nid,"eup_ext",eup_ext,h5in,h5err)
          if (obj_exists(nid,"elow_ext",h5err)) &
            call read_h5(nid,"elow_ext",elow_ext,h5in,h5err)
          if (obj_exists(nid,"dup_ext",h5err)) &
            call read_h5(nid,"dup_ext",dup_ext,h5in,h5err)
          if (obj_exists(nid,"dlow_ext",h5err)) &
            call read_h5(nid,"dlow_ext",dlow_ext,h5in,h5err)
          if (obj_exists(nid,"setlim_ext",h5err)) &
            call read_h5(nid,"setlim_ext",setlim_ext,h5in,h5err)
          if (obj_exists(nid,"reflect_ext",h5err)) &
            call read_h5(nid,"reflect_ext",reflect_ext,h5in,h5err)
          if (obj_exists(nid,"fixpp",h5err)) &
            call read_h5(nid,"fixpp",fixpp_int,h5in,h5err)
          if (fixpp_int.ne.0) fixpp=.true.  ! default off
          call close_group("profile_ext",nid,h5err)
        endif

        !H5: TODO: read neqdsk (from geqdsk_ext)?

        call test_group(sid,"inlibim",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"inlibim",nid,h5err)
          if (obj_exists(nid,"tlibim",h5err)) &
            call read_h5(nid,"tlibim",tlibim,h5in,h5err)
          if (obj_exists(nid,"slibim",h5err)) &
            call read_h5(nid,"slibim",slibim,h5in,h5err)
          if (obj_exists(nid,"fwtlib",h5err)) &
            call read_h5(nid,"fwtlib",fwtlib,h5in,h5err)
          if (obj_exists(nid,"rrrlib",h5err)) &
            call read_h5(nid,"rrrlib",rrrlib,h5in,h5err)
          if (obj_exists(nid,"zzzlib",h5err)) &
            call read_h5(nid,"zzzlib",zzzlib,h5in,h5err)
          if (obj_exists(nid,"aa1lib",h5err)) &
            call read_h5(nid,"aa1lib",aa1lib,h5in,h5err)
          if (obj_exists(nid,"aa8lib",h5err)) &
            call read_h5(nid,"aa8lib",aa8lib,h5in,h5err)
          call close_group("inlibim",nid,h5err)
        endif

        !H5: not used efitin and auxquant?

        call close_group(trim(line),sid,h5err)
        call close_group("time_slice",tid,h5err)
        call close_group("parameters",pid,h5err)
        call close_group("code",cid,h5err)
        call close_group("equilibrium",eqid,h5err)
        call close_h5file(fileid,rootgid,h5err)

#else
        ! this code should not be reachable
        call errctrl_msg('data_input','HDF5 needs to be linked')
        stop
#endif
      endif ! kdata.eq.9

!----------------------------------------------------------------------- 
!--   post-process inputs                                             --
!-----------------------------------------------------------------------
!--warn that idebug and jdebug inputs are depreciated
      if (idebug.ne.0) write(*,*) &
      "idebug input variable is depreciated, set cmake variable instead"
      if (jdebug.ne."NONE") write(*,*) &
      "jdebug input variable is depreciated, set cmake variable instead"
!--   roundoff differences can throw off zlim if limiter corners
!--   are too close to grid points (maybe zlim needs fixing...)
      do i=1,limitr
        ylim(i)=ylim(i)-1.e-10_dp
      enddo
!--   protect against underflow in fitting weights 
      if (abs(fwtdlc).le.1.e-30_dp)  fwtdlc=0.0
      if (abs(fwtcur).le.1.e-30_dp)  fwtcur=0.0
      do i=1,nfcoil
        if (abs(fwtfc(i)).le.1.e-30_dp)  fwtfc(i)=0.0
      enddo
      do i=1,nesum
        if (abs(fwtec(i)).le.1.e-30_dp)  fwtec(i)=0.0
      enddo
      do i=1,magpri
        if (abs(fwtmp2(i)).le.1.e-30_dp)  fwtmp2(i)=0.0
      enddo
      do i=1,nsilop
        if (abs(fwtsi(i)).le.1.e-30_dp)  fwtsi(i)=0.0
      enddo
      do i=1,nmtark
        if (abs(fwtgam(i)).le.1.e-30_dp)  fwtgam(i)=0.0
      enddo
      do i=1,nmsels
        if (abs(fwtbmsels(i)).le.1.e-30_dp)  fwtbmsels(i)=0.0
        if (abs(fwtemsels(i)).le.1.e-30_dp)  fwtemsels(i)=0.0
      enddo
      do i=1,nnece
        if (abs(fwtece0(i)).le.1.e-30_dp)  fwtece0(i)=0.0
      enddo
      if (abs(fwtecebz0).le.1.e-30_dp)  fwtecebz0=0.0
!
      if (nbdryp==-1) nbdryp=nbdry 

!--   Read msels_all.dat if needed
      if (kdomsels.gt.0) then 
        if (rrmsels(1).lt.-5.0) then 
          open(unit=nin,status='old',file='msels_all.dat') 
          do i=1,nmsels 
            read(nin,91008,iostat=istat) bmsels(i),sbmsels(i),rrmsels(i), & 
                 zzmsels(i), l1msels(i),l2msels(i),l4msels(i),emsels(i), & 
                 semsels(i),iemsels(i) 
            iermselt(jtime,i)=iemsels(i) 
          enddo 
          close(unit=nin) 
          if (istat>0) then 
            im1=i-1 
            write (6,*) 'msels_all i=',i,' bmsels(i-1)= ',bmsels(im1) 
            stop 
          endif
        endif 
      endif 
91008 format(9e12.5,i2) 

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
! TODO: devbcin, rnavbcin, devpin, and rnavpin are never defined or set...
!      devbc(jtime)=devbcin
!      rnavbc(jtime)=rnavbcin
!      devp(jtime)=devpin
!      rnavp(jtime)=rnavpin 
! 
!---------------------------------------------------------------------- 
!--   Input FF', P' arrays
!---------------------------------------------------------------------- 
! H5: should this come from OMAS file? (need to replace go to 11777...)
      if (geqdsk_ext.ne.'none') then 
        open(unit=neqdsk,status='old',file=geqdsk_ext) 
        read (neqdsk,11775) (case_ext(i),i=1,6),nh_ext,nw_ext,nh_ext 
        npsi_ext=nw_ext 
#ifdef DEBUG_LEVEL1
        write (nttyo,*) 'npsi_ext,nw_ext=',npsi_ext,nw_ext
#endif
        do i = 1,2 
          read (neqdsk,11773) 
        enddo 
        read (neqdsk,11776) plasma_ext,c_ext,c_ext,c_ext,c_ext 
        if (plasma_ext > 0.0) then 
          sign_ext = -1.0 
        endif 
        read (neqdsk,11773) 
        read (neqdsk,11776) (ffprim_ext(i),i=1,nw_ext) 
        read (neqdsk,11776) (pprime_ext(i),i=1,nw_ext) 
        prbdry=pprime_ext(nw_ext) 
        read (neqdsk,11776) (ffprim_ext(i),i=1,nw_ext) 
        read (neqdsk,11776) (pprime_ext(i),i=1,nw_ext) 
        read (neqdsk,11776) ((psirz_ext,i=1,nw_ext),j=1,nh_ext) 
        read (neqdsk,11776,err=11778) (qpsi_ext(i),i=1,nw_ext) 
        read (neqdsk,11774,err=11778) nbdry_ext,limitr_ext 
        read (neqdsk,11776,err=11778) (rbdry_ext(i),zbdry_ext(i),i=1,nbdry_ext) 
        read (neqdsk,11776,err=11778) (xlim_ext(i),ylim_ext(i),i=1,limitr_ext) 
11773   format (a) 
11774   format (2i5) 
11775   format (6a8,3i4) 
11776   format (5e16.9)
!      if (geqdsk_ext.ne.'none') then 
!        npsi_ext=nw_ext 
!        if (idebug /= 0) write (nttyo,*) 'npsi_ext,nw_ext=',npsi_ext, & 
!          nw_ext 
!        if (plasma_ext > 0.0) then 
!          sign_ext = -1.0 
!        endif 
!        prbdry=pprime_ext(nw_ext) 
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
            ylim(i)=ylim_ext(i)-1.e-10_dp 
          enddo 
        endif 
      endif 
11778 continue
#ifdef DEBUG_LEVEL1
      write (nttyo,*) 'npsi_ext=',npsi_ext
#endif
      if (npsi_ext > 0) then 
#ifdef DEBUG_LEVEL1
        write (nttyo,*) 'scalepp_ext,pprime_ext= ',scalepp_ext, & 
           pprime_ext(1) 
        write (nttyo,*) 'scaleffp_ext,ffpprim_ext= ',scaleffp_ext, & 
           ffprim_ext(1) 
#endif 
        pprime_ext = pprime_ext*darea*sign_ext*scalepp_ext 
        ffprim_ext = ffprim_ext*darea/twopi/tmu*sign_ext*scaleffp_ext 
        prbdry=prbdry*scalepp_ext*scalepp_ext 
#ifdef DEBUG_LEVEL1
        write (nttyo,*) 'scalepp_ext,pprime_ext= ',scalepp_ext, & 
           pprime_ext(1) 
        write (nttyo,*) 'scaleffp_ext,ffpprim_ext= ',scaleffp_ext, & 
           ffprim_ext(1) 
#endif 
 
        if (psin_ext(1) < 0) then 
          do i = 1, npsi_ext 
            psin_ext(i) = real(i-1,dp)/real(npsi_ext-1,dp) 
          enddo 
        endif 
        call zpline(npsi_ext,psin_ext,pprime_ext,bpp_ext,cpp_ext,dpp_ext) 
        call zpline(npsi_ext,psin_ext,ffprim_ext,bfp_ext,cfp_ext,dfp_ext) 
      endif 
!---------------------------------------------------------------------- 
!--   Scale boundary points                                          -- 
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
!--   Reflection, Lower = -Upper                                     -- 
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
        else
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
!--   Reflection, Upper = - Lower                                    -- 
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
        else
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
!--   Scale limiter                                                  -- 
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
!--   recalculate length of default directories in case any change   -- 
!---------------------------------------------------------------------- 
!      call set_table_dir 
!      call efit_read_tables 
 
!--------------------------------------------------------------------- 
!--   specific choice of current profile                            -- 
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
      if (mse_usecer .eq. 1) keecur = 0
      if (mse_usecer .eq. 2 .and. keecur .eq. 0) then
        keecur=2
        keefnc=0
        itek=5
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
!--   save fitting weights for FILE mode                            -- 
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
!--   adjust fit parameters based on basis function selected          -- 
!----------------------------------------------------------------------- 
      select case (kppfnc)
      case (3)
        kppcur = 4 * (kppknt - 1) 
      case (4)
        kppcur = 4 * (kppknt - 1) 
      case (5)
        kppcur = kppcur * (kppknt - 1) 
      case (6)
        kppcur = kppknt * 2 
      end select
      select case (kfffnc)
      case (3)
        kffcur = 4 * (kffknt - 1) 
      case (4)
        kffcur = 4 * (kffknt - 1) 
      case (5)
        kffcur = kffcur * (kffknt - 1) 
      case (6)
        kffcur = kffknt * 2 
      end select 
      select case (kwwfnc)
      case (3)
        kwwcur = 4 * (kwwknt - 1) 
      case (4)
        kwwcur = 4 * (kwwknt - 1) 
      case (5)
        kwwcur = kwwcur * (kwwknt - 1) 
      case (6)
        kwwcur = kwwknt * 2 
      end select 
      if (keecur.gt.0) then 
        select case (keefnc)
        case (3)
          keecur = 4 * (keeknt - 1) 
        case (4)
          keecur = 4 * (keeknt - 1) 
        case (5)
          keecur = keecur * (keeknt - 1) 
        case (6)
          keecur = keeknt * 2 
        end select 
      endif 
! 
      if (fbetan.gt.0.) brsp(nfcoil+jbeta)=alpax(jbeta)*darea 
      if (fli.gt.0.) brsp(nfcoil+kppcur+jli)=gamax(jli)*darea 
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
!--   itek > 100, write out PLTOUT.OUT individually                      --
!--
!--   TODO: nin is closed, what are read statements supposed to do here? 
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
        if (fli.gt.0.0) itell=2
        if (nqwant.gt.0.0) itell=2
        if ((symmetrize).and.(nbdry.gt.1)) itell=3 
      endif 
      if ((iconvr.ne.3).and.(qvfit.gt.0.0)) qenp=qvfit 
      if (ishot.eq.-1) then
        kerror = 1
        call errctrl_msg('data_input','shot number not set')
        return
      endif
      if ((limitr.gt.0).and.(xlim(1).le.-1.0)) then
!        read (nin,5000) (xlim(i),ylim(i),i=1,limitr) 
        ylim(i)=ylim(i)-1.e-10_dp
      endif
!      if ((nbdry.gt.0).and.(rbdry(1).le.-1.0)) read (nin,5020) & 
!          (rbdry(i),zbdry(i),i=1,nbdry) 
      if (kprfit.gt.0) then 
        do i=1,npress 
          premea(i)=pressr(i) 
        enddo 
      endif 
      if (kprfit.gt.0.and.sigpre(1).lt.0.0) then 
        scalep=abs(sigpre(1)) 
        scalemin=abs(sigpre(2)) 
        do i=1,npress 
          sigpre(i)=scalep*pressr(i) 
          sigpre(i)=max(sigpre(i),scalemin) 
        enddo
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
!-------------------------------------------------------------- 
!--   option to symmetrize added 8/14/91 eal                 -- 
!-------------------------------------------------------------- 
      if (symmetrize) then  ! symmetrize the fixed boundery 
        print *, 'option to symmetrize added 8/14/91 eal  ' 
        isetfb=0  ! be sure vertical feedback is off 
        zelip=0 ! must be symmetric about midplane 
        if (nbdry.gt.1) then  ! remove duplicate point 
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
            if (zbdry(i).gt.0.) then
              aup=atan2(zbdry(i),rbdry(i)-rbar) 
              iup=i 
              close=1e30 
              do j=1,nbdry 
                if (zbdry(j).lt.0.) then
                  adn=atan2(zbdry(j),rbdry(j)-rbar) 
                  val=abs(aup+adn) 
                  if (val.lt.close) then
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
!--   Symmetrize the limiter positions for fixed boundary if request-- 
!--   set LIMITR=1000+points for this option                        -- 
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
          open(unit=nin,status='old',file=edatname)
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
          open(unit=nin,status='old',file=edatname)
          bfract=-1. 
          if (tethom(1).lt.0.0) bfract=-tethom(1) 
          read (nin,edat) 
          close(unit=nin) 
        endif 
        if (nbeam.gt.0) then 
          do i=1,nbeam 
            dnbeam(i)=dnbeam(i)*1.e-19_dp
          enddo
        endif

        print *, 'reorder TS data points '   
!--------------------------------------------------------------------- 
!--     reorder TS data points                                      -- 
!--------------------------------------------------------------------- 
        call tsorder(npteth,zteth,dnethom,tethom,sgneth,sgteth) 
        if (sgnemin.lt.0.0) sgnemin=abs(sgnemin)*dnethom(1)*1.e-19_dp & 
                                                *co2cor 
        do i=1,npneth 
          dnethom(i)=dnethom(i)*1.e-19_dp*co2cor 
          sgneth(i)=sgneth(i)*1.e-19_dp*sgnethi*co2cor 
          sgneth(i)=max(sgneth(i),sgnemin) 
        enddo
        if (sgtemin.lt.0.0) sgtemin=abs(sgtemin)*tethom(1) 
        temax=tethom(1) 
        demax=dnethom(1) 
        do i=1,npteth 
          sgteth(i)=sgteth(i)*sgtethi 
          if (bfract.gt.0.0) then 
            tethom(i)=tethom(i)*bfract 
            sgteth(i)=sgteth(i)*bfract 
          endif 
          sgteth(i)=max(sgteth(i),sgtemin) 
          temax=max(temax,tethom(i)) 
          demax=max(demax,dnethom(i)) 
        enddo
        if (cstabte.lt.0.0) cstabte=abs(cstabte)*100./temax 
        if (cstabne.lt.0.0) cstabne=abs(cstabne)*100./demax 
        if (nption.lt.0) then 
          nptionf=-nption 
          if (nptionf.lt.100) then 
            call getfnmu(itimeu,'k',ishot,itime,edatname) 
            edatname='edat_'//edatname(2:7)//'_'//edatname(9:13)//'.cer' 
            open(unit=nin,status='old',file=edatname) 
            bfract=-1. 
            if (tionex(1).lt.0.0) bfract=-tionex(1) 
            read (nin,edat) 
            close(unit=nin) 
            do i=1,nption 
              sigti(i)=sigti(i)*sigtii 
              if (bfract.gt.0.0) then 
                sigti(i)=sigti(i)*bfract 
                tionex(i)=tionex(i)*bfract 
              endif 
            enddo
          endif 
          if (nptionf.gt.100) then 
            nptionf=nptionf-100 
            nption=npteth 
            nptionf=nptef 
            bfract=tionex(1) 
            do i=1,nption 
              sigti(i)=sgteth(i) 
              tionex(i)=tethom(i)*bfract 
              rion(i)=rteth(i) 
              zion(i)=zteth(i) 
            enddo
          endif 
        endif 
      endif

      print *, 'read in limiter data' 
!----------------------------------------------------------------------- 
!---- read in limiter data                                            -- 
!----------------------------------------------------------------------- 
      call getlim(1,xltype,xltype_180) 
! 
      if (iconvr.lt.0) then
        iecurr=1
        ivesel=1
      endif
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
      do i=1,nsilop 
        silopt(jtime,i)=coils(i) 
      enddo
      do i=1,nfcoil 
        fccurt(jtime,i)=brsp(i) 
      enddo
      do i=1,nmtark 
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
      enddo
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
      do i=1,nmsels 
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
      enddo
!---------------------------------------------------------------------- 
!     give the constraint value for matrix routine 
!---------------------------------------------------------------------- 
      do i=1,nnece 
        brspece(jtime,i)=ecefit(i) 
      enddo 
      brspecebz(jtime)=ecebzfit 
      do i=1,magpri 
        expmpi(jtime,i)=expmp2(i) 
      enddo 
!------------------------------------------------------------------------ 
!--   New E-coil connections                   LLao, 95/07/11          --
!------------------------------------------------------------------------ 
      if (size(ecurrt).ge.3) then
        if (ecurrt(3).le.-1.e10_dp) ecurrt(3)=ecurrt(1)
      endif
      if (size(ecurrt).ge.5) then
        if (ecurrt(5).le.-1.e10_dp) ecurrt(5)=ecurrt(1) 
      endif
      if (size(ecurrt).ge.4) then
        if (ecurrt(4).le.-1.e10_dp) ecurrt(4)=ecurrt(2) 
      endif
      if (size(ecurrt).ge.6) then
        if (ecurrt(6).le.-1.e10_dp) ecurrt(6)=ecurrt(2)
      endif 
      do i=1,nesum 
        eccurt(jtime,i)=ecurrt(i) 
      enddo
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
!--   + 0.01 to take care of truncation problem      03/16/91         --
!----------------------------------------------------------------------- 
      timeus=itimeu 
      timems=itime 
      time(jtime)=timems+timeus/1000. 
      bcentr(jtime)=btor 
      do i=1,nco2v 
        denvt(jtime,i)=denv(i) 
      enddo
      do i=1,nco2r 
        denrt(jtime,i)=denr(i) 
      enddo
      do i=1,nacoil 
        accurt(jtime,i)=acoilc(i) 
        caccurt(jtime,i)=acoilc(i) 
      enddo
      vloopt(jtime)=vloop 
      psiref(jtime)=siref 
! 
      gammap=1./gammap 
      gammaf=gammap 
      psibry0=psibry 
#ifdef DEBUG_LEVEL2
      write (6,*) 'DATA_INPUT PSIBRY0= ', psibry0
#endif
! 
      xlmint=xlmin 
      xlmaxt=xlmax 
!
      call zlim(zero,nw,nh,limitr,xlim,ylim,rgrid,zgrid,limfag)

      endif ! (kdata.eq.2).or.(kdata.eq.9)
!----------------------------------------------------------------------- 
!--   set up parameters for all modes                                 --
!----------------------------------------------------------------------- 
      if (kfffnc.eq.8) then 
        rkec=pi/(2.0*dpsiecn) 
      endif 
      chigam=0.0 
      tchimls=0.0 
!----------------------------------------------------------------------- 
!--   DIII-D shot > 112000 use a new set of table for magnetic probes --
!--          shot > 156000 new 2014 set                               --
!--          shot >= 168191 new 2017 set                              --
!----------------------------------------------------------------------- 
      if (pasmat(jtime).le.-1.e3_dp) then 
        negcur=1 
      else 
        negcur=0 
      endif 
      if (abs(pasmat(jtime)).le.cutip.and.iconvr.ge.0) then
        if (iconsi.eq.-1) iconsi=55 
        if (ivesel.gt.10) iconsi=0 
        iexcal=1 
        ivacum=1 
        ibunmn=0 
        ierchk=0 
        nxiter=1 
        iconvr=3 
      endif
!--------------------------------------------------------------------- 
!--   correction to 322 degree probes due to N1 coil                -- 
!--------------------------------------------------------------------- 
      if (oldccomp) then 
        if (n1coil.eq.2.and.ishot.le.108281) then 
          open(unit=60,file=input_dir(1:lindir)//'n1coil.ddd', & 
               status='old')
          j=jtime
          do i=30,magpri67+magpri322 
            read(60,*) namedum,xxxdum,signn1(i) 
            expmpi(j,i)=expmpi(j,i)-signn1(i)*curtn1(j) 
          enddo
          close(unit=60) 
        endif 
      endif 
!--------------------------------------------------------------------- 
!--   correction to 322 and 67 degree probes due to C coil          -- 
!--------------------------------------------------------------------- 
      if (nccoil.eq.1.and.oldccomp) then 
        open(unit=60,file=input_dir(1:lindir)//'ccoil.ddd', & 
             status='old') 
        j=jtime
        do i=30,magpri67+magpri322 
          read(60,*) namedum,signc139 
          expmpi(j,i)=expmpi(j,i)-signc139*curc139(j) 
        enddo 
        do i=1,29 
          read(60,*) namedum,signc79 
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
          do m=1,nsilop 
            silopt(jtime,m)=silopt(jtime,m)+psiref(jtime) 
          enddo
          psirefs(jtime)=psiref(jtime) 
          psiref(jtime)=0. 
        endif 
      endif 
      kconvr=iconvr 
      do kk=1,nwnh 
        www(kk)=zero(kk) 
      enddo
! 
      fwtref=fwtsi(iabs(nslref)) 
      if ((kersil.ne.2).and.(kersil.ne.3)) then
      do m=1,nsilop
        tdata1=serror*abs(silopt(jtime,m)) 
        tdata2=abs(psibit(m))*vbit 
        tdata=max(tdata1,tdata2) 
        sigmafl0(m)=tdata 
        if (tdata.gt.1.0e-10_dp) fwtsi(m)=fwtsi(m)/tdata**nsq 
        if (tdata.le.1.0e-10_dp) fwtsi(m)=0.0 
      enddo
!---------------------------------------------------------------------- 
!--   signal at psi loop # NSLREF is used as reference               -- 
!---------------------------------------------------------------------- 
      if (abs(psibit(iabs(nslref))).le.1.0e-10_dp) then
        coilmx=abs(silopt(jtime,1))
        do i=2,nsilop 
          abcoil=abs(silopt(jtime,i)) 
          coilmx=max(abcoil,coilmx) 
        enddo
        fwtsi(iabs(nslref))=1.0/coilmx**nsq/serror**nsq*fwtref
      endif 
!----------------------------------------------------------------------- 
!--   Fourier expansion of vessel sgments                             -- 
!----------------------------------------------------------------------- 
      else !kersil.eq.2.or.kersil.eq.3
      if (ifitvs.gt.0. .and. nfourier.gt.1) then
        do i=1,nvesel 
          if (rvs(i).ge.1.75_dp.and.zvs(i).ge.0.) & 
            thetav(i)=dasin(zvs(i)/sqrt ((rvs(i)-1.75_dp)**2+(zvs(i))**2)) 
          if (rvs(i).lt.1.75_dp.and.zvs(i).ge.0.) & 
            thetav(i)=pi-dasin(zvs(i)/sqrt ((rvs(i)-1.75_dp)**2+(zvs(i))**2)) 
          if (rvs(i).lt.1.75_dp.and.zvs(i).lt.0.) & 
            thetav(i)=pi-dasin(zvs(i)/sqrt ((rvs(i)-1.75_dp)**2+(zvs(i))**2)) 
          if (rvs(i).ge.1.75_dp.and.zvs(i).lt.0.) & 
            thetav(i)=2*pi+dasin(zvs(i)/sqrt ((rvs(i)-1.75_dp)**2+(zvs(i))**2)) 
          do j=1,nfourier 
            sinta(j,i)=sin(thetav(i)*j) 
            costa(j,i)=cos(thetav(i)*j) 
          enddo 
        enddo 
        do i=1,(2*nfourier+1) 
          do j=1,nvesel 
            if (i.eq.1) vecta(i,j)=1.0 
            if (i.gt.1.and.i.le.(nfourier+1)) vecta(i,j)=costa(i,j) 
            if (i.gt.(nfourier+1)) vecta(i,j)=sinta(i,j) 
          enddo 
        enddo 
      endif
! 
!     if (brsptu(1).le.-1.e-20_dp) & 
!        brsp(1:nfcoil)=brsptu(1:nfcoil)*turnfc(1:nfcoil) 
      if (brsptu(1).gt.-1.e-20_dp) & 
         brsp(1:nfcoil)=brsptu(1:nfcoil)*turnfc(1:nfcoil) 
      reflux=silopt(jtime,iabs(nslref)) 
      do m=1,nsilop 
        tdata1=errsil*abs(silopt(jtime,m)-reflux) 
        tdata2=sicont*rsi(m)*abs(pasmat(jtime)) 
        tdata=max(tdata1,tdata2) 
        tdata2=abs(psibit(m))*vbit 
        tdata=max(tdata,tdata2) 
        sigmafl0(m)=tdata 
        if (tdata.gt.1.0e-10_dp) fwtsi(m)=fwtsi(m)/tdata**nsq 
        if (tdata.le.1.0e-10_dp) fwtsi(m)=0.0 
      enddo
!---------------------------------------------------------------------- 
!--   signal at psi loop #8 is set to zero and used as reference     -- 
!---------------------------------------------------------------------- 
      if (kersil.ne.3) then
        fwtsi(iabs(nslref))=1.0/ersil8**nsq*fwtref
      else
!---------------------------------------------------------------------- 
!--     New option for reference flux loop uncertainty               -- 
!---------------------------------------------------------------------- 
        m=iabs(nslref) 
        tdata1=errsil*abs(silopt(jtime,m)) 
        tdata2=sicont*rsi(m)*abs(pasmat(jtime)) 
        tdata=max(tdata1,tdata2) 
        tdata2=abs(psibit(m))*vbit 
        tdata=max(tdata,tdata2) 
        !SEK: ??? this was commented out but I getting out of bounds error below
!       sigmafl0(m)=tdata 
        !SEK: ???
        if (tdata.gt.1.0e-10_dp) fwtsi(m)=fwtref/tdata**nsq 
        if (tdata.le.1.0e-10_dp) fwtsi(m)=0.0 
      endif 
!SEK: ???      sigmafl0(m)=tdata 
! 
      endif !kersil.eq.2.or.kersil.eq.3
      fwtref=fwtsi(iabs(nslref)) 
      do m=1,magpri 
        tdata1=serror*abs(expmpi(jtime,m)) 
        tdata2=abs(bitmpi(m))*vbit 
        tdata=max(tdata1,tdata2) 
        sigmamp0(m)=tdata 
        if (tdata.gt.1.0e-10_dp) fwtmp2(m)=fwtmp2(m)/tdata**nsq 
        if (tdata.le.1.0e-10_dp) fwtmp2(m)=0.0 
      enddo
      do m=1,nstark 
        tdata=abs(siggam(jtime,m)) 
        if (tdata.gt.1.0e-10_dp) fwtgam(m)=fwtgam(m)/tdata**nsq 
        if (tdata.le.1.0e-10_dp) fwtgam(m)=0.0 
      enddo
! 
#ifdef DEBUG_LEVEL2
      write (6,*) 'DATA fwtbmselt = ',(fwtbmselt(jtime,i),i=1,nmsels) 
      write (6,*) 'DATA sbmselt = ',(sbmselt(jtime,i),i=1,nmsels) 
#endif 
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
#ifdef DEBUG_LEVEL2
      write (6,*) 'DATA fwtbmselt = ', (fwtbmselt(jtime,i),i=1,nmsels)
#endif
! 
      do m=1,nfcoil 
        tdata1=serror*abs(fccurt(jtime,m)) 
        tdata2=abs(bitfc(m))*vbit 
        tdata=max(tdata1,tdata2) 
        sigmaf0(m)=tdata 
        if (tdata.gt.1.0e-10_dp) fwtfc(m)=fwtfc(m)/tdata**nsq 
        if (tdata.le.1.0e-10_dp) fwtfc(m)=0.0 
      enddo
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
!--   diamagetic flux                                                -- 
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
      do i=1,nsilop 
        if (fwtsi(i).gt.0.0) ipsi(jtime)=ipsi(jtime)+1 
      enddo
      ifc(jtime)=0 
      do i=1,nfcoil 
        if (fwtfc(i).gt.0.0) ifc(jtime)=ifc(jtime)+1 
      enddo
      iec(jtime)=0 
      do i=1,nesum 
        if (fwtec(i).gt.0.0) iec(jtime)=iec(jtime)+1 
      enddo 
      imag2(jtime)=0 
      do i=1,magpri 
        if (fwtmp2(i).gt.0.0) imag2(jtime)=imag2(jtime)+1 
      enddo
      kmtark=0 
      klibim=0 
      do i=1,nmtark 
        if (fwtgam(i).gt.0.0) kmtark=kmtark+1 
      enddo
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
#ifdef DEBUG_MSELS
      write (6,*) 'DATA mmbmsels = ', mmbmsels 
      write (6,*) 'bmselt ',(bmselt(1,i),i=1,nmsels) 
      write (6,*) 'iermselt ',(iermselt(1,i),i=1,nmsels) 
#endif 
! 
      iplasm(jtime)=0 
      if (fwtcur.gt.0.0) iplasm(jtime)=1 
      idlopc(jtime)=0 
      if (fwtdlc.gt.0.0) idlopc(jtime)=1 
      cpasma(jtime)=pasmat(jtime) 
      if (iconvr.eq.3.and.ivesel.eq.1) then 
        do i=1,nvesel 
          cpasma(jtime)=cpasma(jtime)-vcurrt(i) 
        enddo
      endif 
      do n=1,nsilop 
        csilop(n,jtime)=silopt(jtime,n) 
      enddo
      itime=time(jtime) 
      timems=itime 
      timeus=(time(jtime)-timems)*1000. 
      timeus=timeus+0.4_dp 
      itimeu=timeus 
!----------------------------------------------------------------------- 
!--   correction for truncation                                       -- 
!----------------------------------------------------------------------- 
      if (itimeu.ge.990) then 
        itime=itime+1 
        itimeu=0 
        time(jtime)=itime 
      endif 
      csiref=psiref(jtime) 
      do i=1,nesum 
        ecurrt(i)=eccurt(jtime,i) 
        cecurr(i)=ecurrt(i) 
      enddo
!      if (kdata.ne.2) then 
!      if ((iecurr.le.0).or.(idodo.gt.0)) go to 520 
!      open(unit=nrsppc,status='old',form='unformatted', & 
!           file=table_dir(1:ltbdir)//'re'//trim(ch1)//trim(ch2)//'.ddd') 
!      read (nrsppc) rsilec 
!      read (nrsppc) rmp2ec 
!      read (nrsppc) gridec 
!      close(unit=nrsppc) 
!      idodo=1 
! 520  continue 
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
      do i=1,nfcoil 
        if (rsisfc(i).le.-1.0) & 
        rsisfc(i)=turnfc(i)**2*twopi*rf(i)/wf(i)/hf(i)*zetafc 
      enddo
!  525 continue 
!----------------------------------------------------------------------- 
!--   read in the advance divertor coil response functions if needed  -- 
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
!--   optional vertical feedback                                  -- 
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
!--     DO DOUBLE PRECISION SUM OF GRDFDB IN M=1 CONFIGURATION      -- 
!--------------------------------------------------------------------- 
        do i=1,nw
          do j=1,nh
            kk=(i-1)*nh+j 
            gsum=0 
            do mmf=1,nfbcoil/2 
              gsum=gsum+grdfdb(kk,mmf) 
              gsum=gsum-grdfdb(kk,mmf+nfbcoil/2) 
            enddo
            grdfdb(kk,1)=gsum 
          enddo
        enddo
        idofb=1 
      endif 
!--------------------------------------------------------------------- 
!--   polarimetry ?                                                 -- 
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
!--   ECE --set kece, kecebz=0 before call setece                   -- 
!--------------------------------------------------------------------- 
      kece=0 
      kecebz=0 
      sigrid(1)=0.0 
      sigrid(nw)=1.0 
      do i=2,nw-1 
        sigrid(i)=1./real(nw-1,dp)*(i-1) 
      enddo
!-------------------------------------------------------------------- 
!--   kinputece=1, get Te, fe, error array from ECE data routine 
!--     hecedata.for with get_hece.for,fftabl.11,fts.pst( copy from 
!--     /u/austin/efit/hecefit/     (MAX AUSTIN)) 
!--     necein,teecein0(necein),feece0(necein),errorece0(necein) 
!--     fe(GHz), Te(Kev) 
!--   ( when kinputece>1,data from K-file ) 
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
!--   toroidal rotation ? Then set up geometric parameters          -- 
!--------------------------------------------------------------------- 
      if (kvtor.gt.0) then 
        do i=1,nw 
          rgrvt(i)=(rgrid(i)/rvtor)**2 
          rgrvt(i)=rgrvt(i)-1. 
          rgsvt(i)=rgrid(i)*rgrvt(i) 
        enddo 
      endif 
!---------------------------------------------------------------------- 
!--   make filement Green's tables only                              -- 
!---------------------------------------------------------------------- 
      if ((iconvr.lt.0).and.(iconvr.gt.-20)) then
      mx=iabs(iconvr) 
      do k=1,mx 
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
        do m=1,nsilop 
          rsilpf(m,k)=gsilpc(m,kk) 
        enddo
        do m=1,magpri 
          rmp2pf(m,k)=gmp2pc(m,kk) 
        enddo
        do ii=1,nw 
          do jj=1,nh 
            kkkk=(ii-1)*nh+jj 
            mj=iabs(j-jj)+1 
            mk=(i-1)*nh+mj 
            gridpf(kkkk,k)=gridpc(mk,ii) 
          enddo
        enddo
      enddo
      mw=nw 
      mh=nh 
      open(unit=nffile,status='old',form='unformatted', & 
           file='rpfxx.dat',iostat=ioerr)
      if (ioerr.eq.0) close(unit=nffile,status='delete')
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
! 
      elseif (iconvr.le.-20) then
      if (aelip.gt.0.0) then 
        do i=1,nw 
          do j=1,nh 
            kk=(i-1)*nh+j 
            erho=sqrt((rgrid(i)-relip)**2+((zgrid(j)-zelip)/eelip)**2) 
            xpsi(kk)=(erho/aelip)**2 
          enddo
        enddo
      endif 
      do m=1,nsilop 
        wsilpc(m)=0.0 
      enddo
      do m=1,magpri 
        wmp2pc(m)=0.0 
      enddo
      do m=1,nfcoil 
        wfcpc(m)=0.0 
      enddo
      do m=1,nesum 
        wecpc(m)=0.0 
      enddo
      do m=1,nvesel 
        wvspc(m)=0.0 
      enddo
      do m=1,nwnh 
        wgridpc(m)=0.0 
      enddo
      wpcpc=0.0 
      npc=0 
! 
      open(unit=nffile,status='old',form='unformatted', & 
           file=table_di2(1:ltbdi2)//'fc'//trim(ch1)//trim(ch2)//'.ddd')
      read (nffile) rfcfc 
      read (nffile) rfcpc 
      close(unit=nffile) 
! 
      do i=1,nw 
        do j=1,nh 
          kk=(i-1)*nh+j 
          if (xpsi(kk).lt.0.0.or.xpsi(kk).gt.1.0) cycle
          npc=npc+1 
          do m=1,nsilop 
            wsilpc(m)=wsilpc(m)+gsilpc(m,kk) 
          enddo
          do m=1,magpri 
            wmp2pc(m)=wmp2pc(m)+gmp2pc(m,kk) 
          enddo
          do m=1,nfcoil 
            wfcpc(m)=wfcpc(m)+rfcpc(m,kk) 
          enddo
          do m=1,nesum 
            wecpc(m)=wecpc(m)+gridec(kk,m) 
          enddo
          do m=1,nvesel 
            wvspc(m)=wvspc(m)+gridvs(kk,m) 
          enddo
          do ii=1,nw 
            do jj=1,nh 
              kkkk=(ii-1)*nh+jj 
              mj=iabs(j-jj)+1 
              mk=(i-1)*nh+mj 
              wgridpc(kkkk)=wgridpc(kkkk)+gridpc(mk,ii) 
              if (xpsi(kkkk).lt.0.0.or.xpsi(kkkk).gt.1.0) cycle
              wpcpc=wpcpc+gridpc(mk,ii) 
            enddo
          enddo
        enddo
      enddo
      xnpc=real(npc,dp) 
      do m=1,nsilop 
        wsilpc(m)=wsilpc(m)/xnpc 
      enddo
      do m=1,magpri 
        wmp2pc(m)=wmp2pc(m)/xnpc 
      enddo
      do m=1,nfcoil 
        wfcpc(m)=wfcpc(m)/xnpc 
      enddo
      do m=1,nesum 
        wecpc(m)=wecpc(m)/xnpc 
      enddo
      do m=1,nvesel 
        wvspc(m)=wvspc(m)/xnpc 
      enddo
      do m=1,nwnh 
        wgridpc(m)=wgridpc(m)/xnpc 
      enddo
      wpcpc=wpcpc/xnpc**2 
! 
      open(unit=nffile,status='old',form='unformatted', & 
           file='rpcxx.dat',iostat=ioerr)
      if (ioerr.eq.0) close(unit=nffile,status='delete')
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
! 
      else ! iconvr.ge.0
      if (ivesel.gt.0) call vescur(jtime) 
      if (nbdry.gt.0) then
      if (islve.gt.0) then
!------------------------------------------------------------------------------ 
!--     Solve equilibrium                                                    -- 
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
        if (kvtor.le.0) then 
!--------------------------------------------------------------------- 
!--       no rotation                                               -- 
!--------------------------------------------------------------------- 
          scc1=sqrt(2./sbeta)/srm**2 
          if (ssrm.lt.0.0) saaa=xlmint/srm/sqrt(1.-2.*scc1) 
          if (ssrm.gt.0.0) saaa=xlmax/srm/sqrt(1.+2.*scc1) 
          srma=srm*saaa 
          dth=twopi/real(nbdry,dp) 
          do i=1,nbdry 
            th=(i-1)*dth 
            rbdry(i)=srma*sqrt(1.-2.*scc1*cos(th)) 
            zbdry(i)=sin(th) 
            zbdry(i)=saaa*zbdry(i)*seee 
          enddo
        else 
!---------------------------------------------------------------------- 
!--       toroidal rotation                                          -- 
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
      endif
!----------------------------------------------------------------------------- 
!--   set up plasma response                                                -- 
!----------------------------------------------------------------------------- 
      do m=1,nbdry 
        do i=1,nw 
          rdif=rbdry(m)-rgrid(i) 
          do j=1,nh 
            k=(i-1)*nh+j 
            zdif=zbdry(m)-zgrid(j) 
            rsum=rdif**2+zdif**2 
            if (rsum.gt.dselsum) then
              rbdrpc(m,k)=psical(rbdry(m),rgrid(i),zdif)*tmu 
            else
!              mk=(i-1)*nh+1 
!              rbdrpc(m,k)=gridpc(mk,i) 
              zdif=dselsum 
              rselsum=rgrid(i)-dselsum 
              rbdrpc(m,k)=psical(rbdry(m),rselsum,zdif)*tmu 
            endif
          enddo
        enddo
      enddo
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
!                 mk=(i-1)*nh+1 
!                 rsolpc(m,k)=gridpc(mk,i) 
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
!--   set up parameters for fixed boundary calculations               -- 
!----------------------------------------------------------------------- 
      if (ifref.eq.-1) ifref=1 
      if (nbdry.gt.1) then 
        delx2=(rbdry(1)-rbdry(nbdry))**2 
        dely2=(zbdry(1)-zbdry(nbdry))**2 
        if ((delx2+dely2).le.1.0e-08_dp) nbdry=nbdry-1 
      endif 
      if (nbdry.ge.10) then
        xmin=rbdry(1) 
        xmax=xmin 
        ymin=zbdry(1) 
        ymax=ymin 
        do i=2,nbdry 
          xmin=min(xmin,rbdry(i)) 
          xmax=max(xmax,rbdry(i)) 
          ymin=min(ymin,zbdry(i)) 
          ymax=max(ymax,zbdry(i)) 
        enddo
        relip=(xmin+xmax)/2. 
        zelip=(ymin+ymax)/2. 
        aelip=(xmax-xmin)/2. 
        eelip=(ymax-ymin)/(xmax-xmin) 
      endif
      if (cfcoil.lt.0.) cfcoil=100./pasmat(jtime)*abs(cfcoil) 
      if (cupdown.lt.0.) cupdown=100./pasmat(jtime)*abs(cupdown) 
!----------------------------------------------------------------------- 
!--   symmetrize  F coil responses if needed                          -- 
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
!--   interpolate to get boundary response functions, first F coils   -- 
!----------------------------------------------------------------------- 
      do n=1,nfcoil
        call sets2d(gridfc(1,n),c,rgrid,nw,bkx,lkx,zgrid,nh,bky, & 
                    lky,wk,ier) 
        do i=1,nbdry 
          call seva2d(bkx,lkx,bky,lky,c,rbdry(i),zbdry(i),pds,ier,n111) 
          rbdrfc(i,n)=pds(1) 
        enddo
        if (nsol.gt.0) then 
          do i=1,nsol 
            call seva2d(bkx,lkx,bky,lky,c,rsol(i),zsol(i),pds,ier,n111) 
            rsolfc(i,n)=pds(1) 
          enddo 
        endif 
      enddo
!---------------------------------------------------------------------- 
!--   make sure interpolations are symmetrized                       -- 
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
!--   advance divertor coil                                           -- 
!----------------------------------------------------------------------- 
      if (iacoil.gt. 0) then 
        do n=1,nacoil 
        call sets2d(gridac(1,n),c,rgrid,nw,bkx,lkx,zgrid,nh,bky, & 
                    lky,wk,ier) 
          do i=1,nbdry 
            call seva2d(bkx,lkx,bky,lky,c,rbdry(i),zbdry(i),pds,ier,n111) 
            rbdrac(i,n)=pds(1) 
          enddo
        enddo
      endif 
!----------------------------------------------------------------------- 
!--   Ohmic coils                                                     -- 
!----------------------------------------------------------------------- 
      if (iecurr.gt.0) then
        do n=1,nesum 
          call sets2d(gridec(1,n),c,rgrid,nw,bkx,lkx,zgrid,nh,bky, & 
                    lky,wk,ier) 
          do i=1,nbdry 
            call seva2d(bkx,lkx,bky,lky,c,rbdry(i),zbdry(i),pds,ier,n111) 
            rbdrec(i,n)=pds(1) 
          enddo
          if (nsol.gt.0) then 
            do i=1,nsol 
              call seva2d(bkx,lkx,bky,lky,c,rsol(i),zsol(i),pds,ier,n111) 
              rsolec(i,n)=pds(1) 
            enddo 
          endif 
        enddo
      endif
!----------------------------------------------------------------------- 
!--   now vessel                                                      -- 
!----------------------------------------------------------------------- 
      if (ivesel.gt.0) then
        do n=1,nvesel 
          call sets2d(gridvs(1,n),c,rgrid,nw,bkx,lkx,zgrid,nh,bky, & 
                      lky,wk,ier) 
          do i=1,nbdry 
            call seva2d(bkx,lkx,bky,lky,c,rbdry(i),zbdry(i),pds,ier,n111) 
            rbdrvs(i,n)=pds(1) 
          enddo
        enddo
      endif
! 
      endif ! nbgry.gt.0 
      endif ! iconvr.ge.0
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
