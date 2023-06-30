#include "config.f"
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          data sets up the magnetic data and weighting arrays.    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**********************************************************************
      subroutine data_input(jtime,kconvr,ktime,kerror)
      use commonblocks,only: c,wk,bkx,bky,wgridpc,rfcpc
      use set_kinds, only: i4
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none

      integer*4, intent(in) :: jtime,ktime
      integer*4, intent(out) :: kconvr,kerror
      integer*4 i,j,k,ii,jj,kk,kkkk,m,n
      integer*4 idoac,mcontr,ktear
      !integer*4 idodo,idovs
      integer*4 istat,lshot,im1,loc,nbdry0,nbdry1,nbdry2,nbdry3,nbryup, &
                iup,nbranch,nqpsi,idn,limupper,lwant,mmemsels,idofb, &
                mw,mh,mmf,mx,ix,mj,mk,npc,npack,nskip,kkl,kku,ier,jupper
      integer*4 nrmin_e,nrmax_e,nzmin_e,nzmax_e,nsw,nsh
      integer*4 n_write !unused
      integer*4 nbabs_ext,jb_ext
      integer*4 ishot_save,itime_save,nbdry_save,iconvr_save, &
                mxiter_save,nxiter_save,ibunmn_save,npress_save
      integer*4 iemsels(nmsels)
      integer*4 ilower(mbdry)
      integer*8 n1d(1),n2d(2)
      !real*4 spatial_avg_ham(nmselp,ngam_vars,ngam_u,ngam_w)
      real*8 currn1,co2cor,fq95
      real*8 plasma,btor,dflux,xltype,xltype_180,vloop,siref,sgnemin, &
             idfila,sigtii,pnbeam,sgtemin,sgnethi,sgtethi, &
             currc79,currc139,currc199, &
             curriu30,curriu90,curriu150,curril30,curril90,curril150, &
             bremsigi,near,adn,val,rnow,znow,bfract, &
             temax,demax,timeus,timems,xlmint,xlmaxt,xxxdum, &
             signc79,signc139,tdata,tdata1,tdata2,coilmx,reflux, &
             zdif,erho,xnpc,ssrm,dth,xrm2,xrvt,th,siwant,rdif,rsum, &
             drgrids,dzgrids,psical,rselsum
      real*8 sicont,s1edge,scalep,scalemin,delx2,dely2,rbar,aup
      real*8 rbmin_e,rbmax_e,zbmin_e,zbmax_e
      real*8 plasma_ext,c_ext,dr_ext,dz_ext,rc_ext,zc_ext,a_ext
      real*8 eup_ext,elow_ext,dup_ext,dlow_ext,setlim_ext
      real*8 plasma_save,btor_save,rcentr_save,fwtcur_save,error_save
      real*8 r0min,r0max,z0min,z0max,zr0min,zr0max,rz0min,rz0max
      real*8 r0ave,z0ave,a0ave,e0top,e0bot,d0top,d0bot,dnmin
      real*8 gridpf(nwnh,mfila),gwork(nsilop,nwnh),rgrids(nw),zgrids(nh)
      real*8 coils(nsilop),expmp2(magpri),acoilc(nacoil), &
             tgamma(nmselp),sgamma(nmselp),rrrgam(nmselp),zzzgam(nmselp), &
             aa1gam(nmselp),aa2gam(nmselp),aa3gam(nmselp),aa4gam(nmselp), &
             aa5gam(nmselp),aa6gam(nmselp),aa7gam(nmselp), &
             tgammauncor(nmselp)
      real*8 bmsels(nmsels),sbmsels(nmsels),fwtbmsels(nmsels), &
             rrmsels(nmsels),zzmsels(nmsels), &
             l1msels(nmsels),l2msels(nmsels),l4msels(nmsels), &
             emsels(nmsels),semsels(nmsels),fwtemsels(nmsels)
      real*8 tlibim(libim),slibim(libim),rrrlib(libim),zzzlib(libim), &
             aa1lib(libim),aa8lib(libim),fwtlib(libim)
      real*8 pds(6),denr(nco2r),denv(nco2v)
      real*8 brsptu(nfsum),brsp_save(nrsmat)
      real*8, dimension(:), allocatable :: expmp2_save,coils_save, &
                                             rbdry_save,zbdry_save, &
                                             fwtsi_save,pressr_save, &
                                             rpress_save,sigpre_save
      real*8, dimension(:,:), allocatable :: psirz_temp
      character(1000) line
      character*10 case_ext(6)
      character*50 edatname
      character(256) table_save
      character*10 namedum,tindex,probeind
      character*2 reflect_ext
      logical shape_ext
      logical file_stat
      logical writepc ! unused variable
      integer*4, parameter :: m_ext=101,nsq=1
      real*8, parameter :: aaslop=0.6_dp,drslop=0.003_dp
      real*8, parameter :: ersil8=1.0e-03_dp,ten2m3=1.0e-03_dp,errorq=1.0e-03_dp
      !data idodo/0/,idovs/0/
      data idoac/0/,mcontr/35/

      namelist/in1/ishot,itime,plasma,itek,itrace,nxiter,fwtcur,kffcur, &
           coils,fwtsi,expmp2,fwtmp2,kppcur,mxiter,ierchk,fwtqa,qemp,error, &
           limitr,xlim,ylim,serror,nbdry,rbdry,zbdry,psibry,nslref,ibunmn, &
           btor,psibit,bitmpi,bitip,icurrt,icinit,brsp,iweigh,qenp,fwtbp, &
           relip,zelip,aelip,eelip,qvfit,fwtdlc,betap0,emp,enp,iconvr,icprof, &
           nextra,ixstrt,scrape,errmin,rbound,npnef,nptef,fwacoil,itimeu, &
           rcentr,rzero,gammap,cfcoil,fczero,fcsum,islve,icntour,iprobe, &
           salpha,srm,sbeta,ifref,isumip,n1coil,ifcurr,iecurr,ecurrt,iecoil, &
           co2cor,vcurrt,dflux,sigdlc,iplim,kinput,limfag,sigprebi,fwtxx, &
           kprfit,pressr,rpress,zpress,sigpre,npress,tethom,rteth,keqdsk, &
           zteth,sgteth,npteth,tionex,rion,zion,sigti,nption,dnethom,zeffvs, &
           rneth,zneth,sgneth,npneth,pbeam,sibeam,nbeam,rzeroj,xalpa,cgama, &
           ivesel,iexcal,iconsi,fwtfc,xltype,kcalpa,kcgama,calpa,iacoil, &
           limid,irfila,jzfila,vloop,iqplot,siref,denr,denv,xgama,sgnemin, &
           nptionf,currn1,ifitvs,bitfc,idfila,relax,saimin,icutfp,acoilc, &
           sigtii,cutip,iavem,pnbeam,xltype_180,sgtemin,sgprmin,elomin,dnmin, &
           sgnethi,fcurbd,pcurbd,prbdry,sgtethi,ndokin,zlowimp,kskipvs,limvs, &
           vcurfb,kpressb,pressbi,prespb,sigppb,kzeroj,rminvs,rmaxvs,errbry, &
           fwtpre,ibtcomp,klabel,zmaxvs,dnbeam,dmass,nmass,condin,iaveus, &
           sgtimin,kwripre,kbound,alphafp,kframe,zbound,vsdamp,zminvs,saicon, &
           kppfnc,kppknt,ppknt,pptens,kfffnc,kffknt,ffknt,fftens,fwtbdry, &
           kwwfnc,kwwknt,wwknt,wwtens,fwtec,fitsiref,bitec,scalepr,scalesir, &
           ppbdry,kppbdry,pp2bdry,kpp2bdry,scalea,sigrbd,sigzbd,nbskip, &
           ffbdry,kffbdry,ff2bdry,kff2bdry,errsil,vbit, &
           wwbdry,kwwbdry,ww2bdry,kww2bdry,f2edge,fe_width,fe_psin,kedgef, &
           ktear,kersil,iout,ixray,pedge,kedgep,pe_width,pe_psin, &
           table_dir,input_dir,store_dir,kautoknt,akchiwt,akerrwt, &
           kakloop,aktol,kakiter,akgamwt,akprewt, &
           kpphord,kffhord,keehord,psiecn,dpsiecn,fitzts,isolve,iplcout, &
           imagsigma,errmag,ksigma,errmagb,brsptu,fitfcsum,fwtfcsum,appendsnap, &
           nbdrymx,nsol,rsol,zsol,fwtsol,efitversion,kbetapr,nbdryp, &
           idebug,jdebug,ifindopt,tolbndpsi,siloplim,use_previous, &
           req_valid,chordv,chordr,nw_sub,nh_sub
      namelist/inwant/psiwant,vzeroj,fwtxxj,fbetap,fbetan,fli,fq95,fqsiw, &
           jbeta,jli,alpax,gamax,jwantm,fwtxxq,fwtxxb,fwtxli,znose, &
           fwtbdry,nqwant,siwantq,n_write,kccoils,ccoils,rexpan, &
           xcoils,kcloops,cloops,xloops,currc79,currc139,nccoil,sizeroj, &
           fitdelz,ndelzon,relaxdz,stabdz,writepc,table_dir,errdelz, &
           oldccomp,nicoil,oldcomp,currc199,curriu30,curriu90, &
           curriu150,curril30,curril90,curril150,ifitdelz,scaledz
      namelist/ink/isetfb,ioffr,ioffz,ishiftz,gain,gainp,idplace, &
           symmetrize,backaverage,lring,cupdown
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
      namelist/inece/necein,teecein0,feece0,errorece0,fwtece0,fwtecebz0, &
           ecefit,ecebzfit,kfitece,kinputece,kcallece,nharm, &
           kfixro,rteo,zteo,kfixrece,rtep,rtem,rpbit,rmbit,robit, &
           nfit,kcmin,fwtnow,kdoece,mtxece,nconstr,eceiter,eceerror
      namelist/iner/keecur,ecurbd,keefnc,eetens,keebdry,kee2bdry, &
           eebdry,ee2bdry,eeknt,keeknt,keehord
      namelist/insxr/ksxr0,ksxr2,idosxr
      namelist/edgep/symmetrize, &
           rpress,pressr,sigpre,npress,kprfit,kpressb,ndokin, &
           kppfnc,kfffnc,kffcur,kppcur,mxiter,error,errmin,keecur
      namelist/edat/nption,tionex,sigti,rion,zion, &
           npneth,dnethom,sgneth,rneth,zneth, &
           npteth,tethom,sgteth,rteth,zteth, &
           nbrmcrd,bremin,bremsig,brmrtan,brmzelev,ivbcuse, &
           sigtii,sgnethi,sgtethi,bremsigi, &
           npress,rpress,zpress,pressr,sigpre
      namelist/invt/omegat,nomegat,enw,emw,betapw0,kvtor,kwwcur,rvtor, &
           wcurbd,preswb,fwtprw,npresw,presw,sigprw,rpresw, &
           zpresw,kplotp,sbetaw,nsplot,comega,kcomega,xomega, &
           kdovt,romegat,zomegat,sigome,scalepw, &
           kwwfnc,kwwknt,wwknt,wwtens
      namelist/profile_ext/npsi_ext,pprime_ext,ffprim_ext,psin_ext, &
           geqdsk_ext,sign_ext,scalepp_ext,scaleffp_ext,shape_ext,dr_ext, &
           dz_ext,rc_ext,zc_ext,a_ext,eup_ext,elow_ext,dup_ext,dlow_ext, &
           setlim_ext,reflect_ext,fixpp
      namelist/out1/ishot,itime,betap0,rzero,qenp,enp,emp,plasma, &
           expmp2,coils,btor,rcentr,brsp,icurrt,rbdry,zbdry, &
           nbdry,fwtsi,fwtcur,mxiter,nxiter,limitr,xlim,ylim,error, &
           iconvr,ibunmn,pressr,rpress,nqpsi,npress,sigpre

      ! set defaults for all modes within loop
      kerror=0
      errorm=1.
      brsp_save=brsp
      brsp=0.0
      brsp(1)=-1.e-20_dp
      brsptu=0.0
      brsptu(1)=-1.e-20_dp
      sicont=tmu*drslop/aaslop

      ! initialize variables with zeros
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
      xgama=0.0
      sigtii=0.0
      sgnethi=0.0
      sgtethi=0.0
      zlowimp=0.0
      pressbi=0.0
      prespb=0.0
      dnbeam=0.0
      dmass=0.0
      sgtimin=0.0
      scalepr=0.0
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
      feece0=0.0
      errorece0=0.0
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
      rbdry=0.0
      zbdry=0.0
      fpol=0.0
      pres=0.0
      ffprim=0.0
      pprime=0.0
      saisref=0.0
      cdelz=0.0
      if (jtime.eq.0 .or. isicinit.ge.0 .or. isicinit.eq.-3) then
        cdeljsum=0.0
        psi=0.0
      endif
      iconsi=-1
      iconvr=2
      nxiter=1
!---------------------------------------------------------------------- 
!--   SNAP Mode
!--   Restore fitting weights
!--------------------------------------------------------------------- 
      snap: if ((kdata.ne.1).and.(kdata.ne.2)) then
      fwtdlc=swtdlc
      fwtcur=swtcur
      if(fwtqa.ne.0.0) fwtqa=1.
      if(fwtbp.ne.0.0) fwtbp=1.
      fwtfc=swtfc
      fwtec=swtec
      fwtmp2=swtmp2
      fwtsi=swtsi
      fwtgam=swtgam
!----------------------------------------------------------------------- 
!--   Set edge pedestal tanh paramters                                -- 
!----------------------------------------------------------------------- 
      if (fitzts.eq.'te'.and.ztserr(jtime)) then 
        nbdry=1 
        rbdry(1)=1.94_dp 
        zbdry(1)=ztssym(jtime)+0.5_dp*ztswid(jtime) 
      endif 
!---------------------------------------------------------------------- 
!--   File Mode
!--   initialize inputs
!---------------------------------------------------------------------- 
      else snap
      denr=0.
      denv=0.
      co2cor=1.0
      currn1=0.0
      fwtbmsels=0.0
      fwtemsels=0.0
      ktear=0
      pnbeam=0.0 
      rrmsels(1)=-10. 
      sgnemin=0.0 
      sgtemin=0.0 
      siref=0. 
      xltype=0.0 
      xltype_180=0.0 
      vloop=0. 

      fwtlib=0.0
      rrrlib=0.0
      zzzlib=0.0

      a_ext=-10. 
      dr_ext=0.0 
      dz_ext=0.0 
      dup_ext=-10. 
      dlow_ext=-10. 
      eup_ext=-10. 
      elow_ext=-10. 
      rc_ext=-10. 
      zc_ext=-10. 
      reflect_ext='no' 
      setlim_ext=-10. 
      shape_ext=.false. 

      table_save=table_dir

      call set_defaults

      file_type: if (kdata.eq.2) then 
!---------------------------------------------------------------------- 
!--   Read ascii input files                                         -- 
!---------------------------------------------------------------------- 
      open(unit=nin,status='old',file=ifname(jtime),iostat=istat)
      if (istat.ne.0) then
        call errctrl_msg('data_input',trim(ifname(jtime))//' not found')
        stop
      endif
      read (nin,in1,iostat=istat)
      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist in1: '//trim(line)
        kerror=1
        return
      endif

      rewind(nin)
      read (nin,ink,iostat=istat)
      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist ink: '//trim(line)
        stop
      endif

      rewind(nin)
      read (nin,ins,iostat=istat)
      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist ink: '//trim(line)
        stop
      endif

      rewind(nin)
      read (nin,in_msels,iostat=istat) 
      if (istat>0) then 
        backspace(nin) 
        read(nin,fmt='(A)') line 
        write(*,'(A)') 'Invalid line in namelist in_msels: '//trim(line)
        stop
      endif

      rewind(nin)
      read (nin,ina,iostat=istat)
      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist ina: '//trim(line)
        stop
      endif

      rewind(nin)
      read (nin,inece,iostat=istat) 
      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist inece: '//trim(line)
        stop
      endif

      rewind(nin)
      read (nin,edgep,iostat=istat) 
      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist edgep: '//trim(line)
        stop
      endif

      rewind(nin)
      read (nin,iner,iostat=istat) 
      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist iner: '//trim(line)
        stop
      endif

      rewind(nin)
      read (nin,insxr,iostat=istat) 
      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist insxr: '//trim(line)
        stop
      endif

      rewind(nin)
      read (nin,inwant,iostat=istat) 
      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist inwant: '//trim(line)
        stop
      endif

      rewind(nin)
      read (nin,invt,iostat=istat)
      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist invt: '//trim(line)
        stop
      endif
      close(unit=nin)

!--   Input FF', P' arrays
#ifdef DEBUG_LEVEL1
      write (nttyo,*) 'DATA_INPUT geqdsk_ext=',geqdsk_ext
#endif
      open(unit=nin,status='old',file=ifname(jtime))
      read(nin,profile_ext,iostat=istat)
      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist profile_ext: '//trim(line)
        stop
      endif
      close(unit=nin)
#ifdef DEBUG_LEVEL1
      write (nttyo,*) 'DATA_INPUT geqdsk_ext=',geqdsk_ext
#endif
      if (geqdsk_ext.ne.'none') then
        open(unit=neqdsk,status='old',file=geqdsk_ext)
        read (neqdsk,11775) (case_ext(i),i=1,6),nh_ext,nw_ext,nh_ext
        allocate(psirz_temp(nw_ext,nh_ext),psirz_ext(nw_ext*nh_ext), &
                 pprime_ext(nw_ext),ffprim_ext(nw_ext),qpsi_ext(nw_ext))
        read (neqdsk,11773)
        read (neqdsk,11776) c_ext,c_ext,simag_ext,psibry_ext,c_ext
        read (neqdsk,11776) plasma_ext,c_ext,c_ext,c_ext,c_ext
        read (neqdsk,11773)
        read (neqdsk,11776) (ffprim_ext(i),i=1,nw_ext)
        read (neqdsk,11776) (pprime_ext(i),i=1,nw_ext)
        prbdry=pprime_ext(nw_ext) 
        read (neqdsk,11776) (ffprim_ext(i),i=1,nw_ext)
        read (neqdsk,11776) (pprime_ext(i),i=1,nw_ext)
        read (neqdsk,11776) ((psirz_temp(i,j),i=1,nw_ext),j=1,nh_ext)
        psirz_ext=0.0
        read (neqdsk,11776,err=11777) (qpsi_ext(i),i=1,nw_ext)
        read (neqdsk,11774,err=11777) nbdry_ext,limitr_ext
        allocate(rbdry_ext(nbdry_ext),zbdry_ext(nbdry_ext), &
                 xlim_ext(limitr_ext),ylim_ext(limitr_ext))
        read (neqdsk,11776,err=11777) (rbdry_ext(i),zbdry_ext(i),i=1,nbdry_ext)
        read (neqdsk,11776,err=11777) (xlim_ext(i),ylim_ext(i),i=1,limitr_ext)
11773   format (a)
11774   format (2i5)
11775   format (6a8,3i4)
11776   format (5e16.9)
        if ((icinit.eq.-3).or.(icinit.eq.-4 .and. jtime.eq.1)) then
          if ((nw.ne.nw_ext) .or. (nh.ne.nh_ext)) then
            call errctrl_msg('data_input', &
                             'restart file must have same dimensions')
            stop
          endif
          do i=1,nw
            do j=1,nh
              kk=(i-1)*nh+j
              psirz_ext(kk)=psirz_temp(i,j)
            enddo
          enddo
!--       Read fcoil currents
!--       Avoid overwriting other variables that will impact restart
          allocate(expmp2_save(magpri),coils_save(nsilop), &
                   rbdry_save(mbdry),zbdry_save(mbdry), &
                   fwtsi_save(nsilop),pressr_save(mpress), &
                   rpress_save(mpress),sigpre_save(mpress))
          ishot_save=ishot
          itime_save=itime
          plasma_save=plasma
          expmp2_save=expmp2
          coils_save=coils
          btor_save=btor
          rcentr_save=rcentr
          brsp_save=brsp
          nbdry_save=nbdry
          rbdry_save=rbdry
          zbdry_save=zbdry
          fwtsi_save=fwtsi
          fwtcur_save=fwtcur
          nxiter_save=nxiter
          mxiter_save=mxiter
          error_save=error
          iconvr_save=iconvr
          ibunmn_save=ibunmn
          pressr_save=pressr
          rpress_save=rpress
          npress_save=npress
          sigpre_save=sigpre
          read (neqdsk,out1,iostat=istat)
          if (istat>0) then 
            backspace(nin)
            read(nin,fmt='(A)') line
            write(*,'(A)') 'Invalid line in namelist out1: '//trim(line)
            stop
          endif
          allocate(fcoil_ext(nfsum))
          fcoil_ext=brsp(1:nfsum)
          ishot=ishot_save
          itime=itime_save
          plasma=plasma_save
          expmp2=expmp2_save
          coils=coils_save
          btor=btor_save
          rcentr=rcentr_save
          brsp=brsp_save
          nbdry=nbdry_save
          rbdry=rbdry_save
          zbdry=zbdry_save
          fwtsi=fwtsi_save
          fwtcur=fwtcur_save
          nxiter=nxiter_save
          mxiter=mxiter_save
          error=error_save
          iconvr=iconvr_save
          ibunmn=ibunmn_save
          pressr=pressr_save
          rpress=rpress_save
          npress=npress_save
          sigpre=sigpre_save
          deallocate(expmp2_save,coils_save,rbdry_save,zbdry_save, &
                     fwtsi_save,pressr_save,rpress_save,sigpre_save)
        endif
      endif
11777 continue

!--   Read Li beam data
      open(unit=nin,status='old',file=ifname(jtime))
      read (nin,inlibim,iostat=istat)
      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist inlibim: '//trim(line)
        stop
      endif
      close(unit=nin)
      else file_type
!---------------------------------------------------------------------- 
!--     HDF5 file mode                                               -- 
!---------------------------------------------------------------------- 
#if defined(USE_HDF5)
        inquire(file=trim(ifname(1)),exist=file_stat)
        if (.not. file_stat) then
          call errctrl_msg('data_input',trim(line)//' not found')
          kerror=1
          return
        endif
        call fch5init
        call open_oldh5file(trim(ifname(1)),fileid,rootgid,h5in,h5err)
        call test_group(rootgid,"equilibrium",file_stat,h5err)
        if (.not. file_stat) then
          call errctrl_msg('data_input','equilibrium group not found')
          kerror=1
          return
        endif
        call open_group(rootgid,"equilibrium",eqid,h5err)
        call test_group(eqid,"code",file_stat,h5err)
        if (.not. file_stat) then
          call errctrl_msg('data_input','code group not found')
          kerror=1
          return
        endif
        call open_group(eqid,"code",cid,h5err)
        call test_group(cid,"parameters",file_stat,h5err)
        if (.not. file_stat) then
          call errctrl_msg('data_input','parameters group not found')
          kerror=1
          return
        endif
        call open_group(cid,"parameters",pid,h5err)
        call test_group(pid,"time_slice",file_stat,h5err)
        if (.not. file_stat) then
          call errctrl_msg('data_input','time_slice group not found')
          kerror=1
          return
        endif
        call open_group(pid,"time_slice",tid,h5err)
        write(tindex,"(I0)") jtime-1+rank*ktime
        call test_group(tid,trim(tindex),file_stat,h5err)
        if (.not. file_stat) then
          call errctrl_msg('data_input', &
                           trim(tindex)//' group not found')
          stop
        endif
        call open_group(tid,trim(tindex),sid,h5err)

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
          call read_h5_ex(nid,"coils",coils,h5in,h5err)
          call read_h5_ex(nid,"fwtsi",fwtsi,h5in,h5err)
          call read_h5_ex(nid,"expmp2",expmp2,h5in,h5err)
          call read_h5_ex(nid,"fwtmp2",fwtmp2,h5in,h5err)
          call read_h5_ex(nid,"kppcur",kppcur,h5in,h5err)
          call read_h5_ex(nid,"mxiter",mxiter,h5in,h5err)
          call read_h5_ex(nid,"ierchk",ierchk,h5in,h5err)
          call read_h5_ex(nid,"fwtqa",fwtqa,h5in,h5err)
          call read_h5_ex(nid,"qemp",qemp,h5in,h5err)
          call read_h5_ex(nid,"error",error,h5in,h5err)
          call read_h5_ex(nid,"limitr",limitr,h5in,h5err)
          call read_h5_ex(nid,"xlim",xlim,h5in,h5err)
          call read_h5_ex(nid,"ylim",ylim,h5in,h5err)
          call read_h5_ex(nid,"serror",serror,h5in,h5err)
          call read_h5_ex(nid,"nbdry",nbdry,h5in,h5err)
          call read_h5_ex(nid,"rbdry",rbdry,h5in,h5err)
          call read_h5_ex(nid,"zbdry",zbdry,h5in,h5err)
          call read_h5_ex(nid,"psibry",psibry,h5in,h5err)
          call read_h5_ex(nid,"nslref",nslref,h5in,h5err)
          call read_h5_ex(nid,"ibunmn",ibunmn,h5in,h5err)
          call read_h5_ex(nid,"btor",btor,h5in,h5err)
          call read_h5_ex(nid,"psibit",psibit,h5in,h5err)
          call read_h5_ex(nid,"bitmpi",bitmpi,h5in,h5err)
          call read_h5_ex(nid,"bitip",bitip,h5in,h5err)
          call read_h5_ex(nid,"icurrt",icurrt,h5in,h5err)
          call read_h5_ex(nid,"icinit",icinit,h5in,h5err)
          call read_h5_ex(nid,"brsp",brsp,h5in,h5err)
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
          call read_h5_ex(nid,"rzero",rzero,h5in,h5err)
          call read_h5_ex(nid,"gammap",gammap,h5in,h5err)
          call read_h5_ex(nid,"cfcoil",cfcoil,h5in,h5err)
          call read_h5_ex(nid,"fczero",fczero,h5in,h5err)
          call read_h5_ex(nid,"fcsum",fcsum,h5in,h5err)
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
          call read_h5_ex(nid,"ecurrt",ecurrt,h5in,h5err)
          call read_h5_ex(nid,"iecoil",iecoil,h5in,h5err)
          call read_h5_ex(nid,"co2cor",co2cor,h5in,h5err)
          call read_h5_ex(nid,"vcurrt",vcurrt,h5in,h5err)
          call read_h5_ex(nid,"dflux",dflux,h5in,h5err)
          call read_h5_ex(nid,"sigdlc",sigdlc,h5in,h5err)
          call read_h5_ex(nid,"iplim",iplim,h5in,h5err)
          call read_h5_ex(nid,"kinput",kinput,h5in,h5err)
          call read_h5_ex(nid,"limfag",limfag,h5in,h5err)
          call read_h5_ex(nid,"sigprebi",sigprebi,h5in,h5err)
          call read_h5_ex(nid,"fwtxx",fwtxx,h5in,h5err)
          call read_h5_ex(nid,"kprfit",kprfit,h5in,h5err)
          call read_h5_ex(nid,"pressr",pressr,h5in,h5err)
          call read_h5_ex(nid,"rpress",rpress,h5in,h5err)
          call read_h5_ex(nid,"zpress",zpress,h5in,h5err)
          call read_h5_ex(nid,"sigpre",sigpre,h5in,h5err)
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
          call read_h5_ex(nid,"pbeam",pbeam,h5in,h5err)
          call read_h5_ex(nid,"sibeam",sibeam,h5in,h5err)
          call read_h5_ex(nid,"nbeam",nbeam,h5in,h5err)
          call read_h5_ex(nid,"rzeroj",rzeroj,h5in,h5err)
          call read_h5_ex(nid,"xalpa",xalpa,h5in,h5err)
          ! need to ensure array is square
          call read_h5_sq(nid,"cgama",cgama,h5in,h5err)
          call read_h5_ex(nid,"ivesel",ivesel,h5in,h5err)
          call read_h5_ex(nid,"iexcal",iexcal,h5in,h5err)
          call read_h5_ex(nid,"iconsi",iconsi,h5in,h5err)
          call read_h5_ex(nid,"fwtfc",fwtfc,h5in,h5err)
          call read_h5_ex(nid,"xltype",xltype,h5in,h5err)
          call read_h5_ex(nid,"kcalpa",kcalpa,h5in,h5err)
          call read_h5_ex(nid,"kcgama",kcgama,h5in,h5err)
          ! need to ensure array is square
          call read_h5_sq(nid,"calpa",calpa,h5in,h5err)
          call read_h5_ex(nid,"iacoil",iacoil,h5in,h5err)
          call read_h5_ex(nid,"limid",limid,h5in,h5err)
          call read_h5_ex(nid,"irfila",irfila,h5in,h5err)
          call read_h5_ex(nid,"jzfila",jzfila,h5in,h5err)
          call read_h5_ex(nid,"vloop",vloop,h5in,h5err)
          call read_h5_ex(nid,"iqplot",iqplot,h5in,h5err)
          call read_h5_ex(nid,"siref",siref,h5in,h5err)
          call read_h5_ex(nid,"denr",denr,h5in,h5err)
          call read_h5_ex(nid,"denv",denv,h5in,h5err)
          call read_h5_ex(nid,"xgama",xgama,h5in,h5err)
          call read_h5_ex(nid,"sgnemin",sgnemin,h5in,h5err)
          call read_h5_ex(nid,"nptionf",nptionf,h5in,h5err)
          call read_h5_ex(nid,"currn1",currn1,h5in,h5err)
          call read_h5_ex(nid,"ifitvs",ifitvs,h5in,h5err)
          call read_h5_ex(nid,"bitfc",bitfc,h5in,h5err)
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
          call read_h5_ex(nid,"sgnethi",sgnethi,h5in,h5err)
          call read_h5_ex(nid,"fcurbd",fcurbd,h5in,h5err)
          call read_h5_ex(nid,"pcurbd",pcurbd,h5in,h5err)
          call read_h5_ex(nid,"prbdry",prbdry,h5in,h5err)
          call read_h5_ex(nid,"sgtethi",sgtethi,h5in,h5err)
          call read_h5_ex(nid,"ndokin",ndokin,h5in,h5err)
          call read_h5_ex(nid,"zlowimp",zlowimp,h5in,h5err)
          call read_h5_ex(nid,"kskipvs",kskipvs,h5in,h5err)
          call read_h5_ex(nid,"limvs",limvs,h5in,h5err)
          call read_h5_ex(nid,"vcurfb",vcurfb,h5in,h5err)
          call read_h5_ex(nid,"kpressb",kpressb,h5in,h5err)
          call read_h5_ex(nid,"pressbi",pressbi,h5in,h5err)
          call read_h5_ex(nid,"prespb",prespb,h5in,h5err)
          call read_h5_ex(nid,"sigppb",sigppb,h5in,h5err)
          call read_h5_ex(nid,"kzeroj",kzeroj,h5in,h5err)
          call read_h5_ex(nid,"rminvs",rminvs,h5in,h5err)
          call read_h5_ex(nid,"rmaxvs",rmaxvs,h5in,h5err)
          call read_h5_ex(nid,"errbry",errbry,h5in,h5err)
          call read_h5_ex(nid,"fwtpre",fwtpre,h5in,h5err)
          call read_h5_ex(nid,"ibtcomp",ibtcomp,h5in,h5err)
          call read_h5_ex(nid,"klabel",klabel,h5in,h5err)
          call read_h5_ex(nid,"zmaxvs",zmaxvs,h5in,h5err)
          call read_h5_ex(nid,"dnbeam",dnbeam,h5in,h5err)
          call read_h5_ex(nid,"dmass",dmass,h5in,h5err)
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
          call read_h5_ex(nid,"ppknt",ppknt,h5in,h5err)
          call read_h5_ex(nid,"pptens",pptens,h5in,h5err)
          call read_h5_ex(nid,"kfffnc",kfffnc,h5in,h5err)
          call read_h5_ex(nid,"kffknt",kffknt,h5in,h5err)
          call read_h5_ex(nid,"ffknt",ffknt,h5in,h5err)
          call read_h5_ex(nid,"fftens",fftens,h5in,h5err)
          call read_h5_ex(nid,"fwtbdry",fwtbdry,h5in,h5err)
          call read_h5_ex(nid,"kwwfnc",kwwfnc,h5in,h5err)
          call read_h5_ex(nid,"kwwknt",kwwknt,h5in,h5err)
          call read_h5_ex(nid,"wwknt",wwknt,h5in,h5err)
          call read_h5_ex(nid,"wwtens",wwtens,h5in,h5err)
          call read_h5_ex(nid,"fwtec",fwtec,h5in,h5err)
          call read_h5_ex(nid,"fitsiref",fitsiref,h5in,h5err)
          call read_h5_ex(nid,"bitec",bitec,h5in,h5err)
          call read_h5_ex(nid,"scalepr",scalepr,h5in,h5err)
          call read_h5_ex(nid,"scalesir",scalesir,h5in,h5err)
          call read_h5_ex(nid,"ppbdry",ppbdry,h5in,h5err)
          call read_h5_ex(nid,"kppbdry",kppbdry,h5in,h5err)
          call read_h5_ex(nid,"pp2bdry",pp2bdry,h5in,h5err)
          call read_h5_ex(nid,"kpp2bdry",kpp2bdry,h5in,h5err)
          call read_h5_ex(nid,"scalea",scalea,h5in,h5err)
          call read_h5_ex(nid,"sigrbd",sigrbd,h5in,h5err)
          call read_h5_ex(nid,"sigzbd",sigzbd,h5in,h5err)
          call read_h5_ex(nid,"nbskip",nbskip,h5in,h5err)
          call read_h5_ex(nid,"ffbdry",ffbdry,h5in,h5err)
          call read_h5_ex(nid,"kffbdry",kffbdry,h5in,h5err)
          call read_h5_ex(nid,"ff2bdry",ff2bdry,h5in,h5err)
          call read_h5_ex(nid,"kff2bdry",kff2bdry,h5in,h5err)
          call read_h5_ex(nid,"errsil",errsil,h5in,h5err)
          call read_h5_ex(nid,"vbit",vbit,h5in,h5err)
          call read_h5_ex(nid,"wwbdry",wwbdry,h5in,h5err)
          call read_h5_ex(nid,"kwwbdry",kwwbdry,h5in,h5err)
          call read_h5_ex(nid,"ww2bdry",ww2bdry,h5in,h5err)
          call read_h5_ex(nid,"kww2bdry",kww2bdry,h5in,h5err)
          call read_h5_ex(nid,"f2edge",f2edge,h5in,h5err)
          call read_h5_ex(nid,"fe_width",fe_width,h5in,h5err)
          call read_h5_ex(nid,"fe_psin",fe_psin,h5in,h5err)
          call read_h5_ex(nid,"kedgef",kedgef,h5in,h5err)
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
          call read_h5_ex(nid,"fitzts",fitzts,h5in,h5err)
          call read_h5_ex(nid,"isolve",isolve,h5in,h5err)
          call read_h5_ex(nid,"iplcout",iplcout,h5in,h5err)
          call read_h5_ex(nid,"brsptu",brsptu,h5in,h5err)
          call read_h5_ex(nid,"fitfcsum",fitfcsum,h5in,h5err)
          call read_h5_ex(nid,"fwtfcsum",fwtfcsum,h5in,h5err)
          call read_h5_ex(nid,"appendsnap",appendsnap,h5in,h5err)
          call read_h5_ex(nid,"idebug",idebug,h5in,h5err)
          call read_h5_ex(nid,"nbdrymx",nbdrymx,h5in,h5err)
          call read_h5_ex(nid,"nsol",nsol,h5in,h5err)
          call read_h5_ex(nid,"rsol",rsol,h5in,h5err)
          call read_h5_ex(nid,"zsol",zsol,h5in,h5err)
          call read_h5_ex(nid,"fwtsol",fwtsol,h5in,h5err)
          call read_h5_ex(nid,"efitversion",efitversion,h5in,h5err)
          call read_h5_ex(nid,"kbetapr",kbetapr,h5in,h5err)
          call read_h5_ex(nid,"nbdryp",nbdryp,h5in,h5err)
          call read_h5_ex(nid,"jdebug",jdebug,h5in,h5err)
          call read_h5_ex(nid,"ifindopt",ifindopt,h5in,h5err)
          call read_h5_ex(nid,"tolbndpsi",tolbndpsi,h5in,h5err)
          call read_h5_ex(nid,"siloplim",siloplim,h5in,h5err)
          call read_h5_ex(nid,"use_previous",use_previous,h5in,h5err)
          call read_h5_ex(nid,"nw_sub",nw_sub,h5in,h5err)
          call read_h5_ex(nid,"nh_sub",nh_sub,h5in,h5err)
          call close_group("in1",nid,h5err)
        endif
   
        call test_group(sid,"ink",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"ink",nid,h5err)
          call read_h5_ex(nid,"isetfb",isetfb,h5in,h5err)
          call read_h5_ex(nid,"ioffr",ioffr,h5in,h5err)
          call read_h5_ex(nid,"ioffz",ioffz,h5in,h5err)
          call read_h5_ex(nid,"ishiftz",ishiftz,h5in,h5err)
          call read_h5_ex(nid,"gain",gain,h5in,h5err)
          call read_h5_ex(nid,"gainp",gainp,h5in,h5err)
          call read_h5_ex(nid,"idplace",idplace,h5in,h5err)
          call read_h5_ex(nid,"symmetrize",symmetrize,h5in,h5err)
          call read_h5_ex(nid,"backaverage",backaverage,h5in,h5err)
          call read_h5_ex(nid,"lring",lring,h5in,h5err)
          call read_h5_ex(nid,"cupdown",cupdown,h5in,h5err)
          call close_group("ink",nid,h5err)
        endif
   
        call test_group(sid,"ins",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"ins",nid,h5err)
          call read_h5_ex(nid,"tgamma",tgamma,h5in,h5err)
          call read_h5_ex(nid,"sgamma",sgamma,h5in,h5err)
          call read_h5_ex(nid,"fwtgam",fwtgam,h5in,h5err)
          call read_h5_ex(nid,"rrrgam",rrrgam,h5in,h5err)
          call read_h5_ex(nid,"zzzgam",zzzgam,h5in,h5err)
          call read_h5_ex(nid,"aa1gam",aa1gam,h5in,h5err)
          call read_h5_ex(nid,"aa2gam",aa2gam,h5in,h5err)
          call read_h5_ex(nid,"aa3gam",aa3gam,h5in,h5err)
          call read_h5_ex(nid,"aa4gam",aa4gam,h5in,h5err)
          call read_h5_ex(nid,"aa5gam",aa5gam,h5in,h5err)
          call read_h5_ex(nid,"aa6gam",aa6gam,h5in,h5err)
          call read_h5_ex(nid,"aa7gam",aa7gam,h5in,h5err)
          call read_h5_ex(nid,"iplots",iplots,h5in,h5err)
          call read_h5_ex(nid,"kdomse",kdomse,h5in,h5err)
          call read_h5_ex(nid,"msebkp",msebkp,h5in,h5err)
          call read_h5_ex(nid,"msefitfun",msefitfun,h5in,h5err)
          call read_h5_ex(nid,"mse_quiet",mse_quiet,h5in,h5err)
          call read_h5_ex(nid,"mse_spave_on",mse_spave_on,h5in,h5err)
          call read_h5_ex(nid,"kwaitmse",kwaitmse,h5in,h5err)
          call read_h5_ex(nid,"dtmsefull",dtmsefull,h5in,h5err)
          call read_h5_ex(nid,"mse_strict",mse_strict,h5in,h5err)
          call read_h5_ex(nid,"t_max_beam_off",t_max_beam_off,h5in,h5err)
          call read_h5_ex(nid,"ok_30rt",ok_30rt,h5in,h5err)
          call read_h5_ex(nid,"ok_210lt",ok_210lt,h5in,h5err)
          call read_h5_ex(nid,"mse_usecer",mse_usecer,h5in,h5err)
          call read_h5_ex(nid,"mse_certree",mse_certree,h5in,h5err)
          call read_h5_ex(nid,"mse_use_cer330",mse_use_cer330,h5in,h5err)
          call read_h5_ex(nid,"mse_use_cer210",mse_use_cer210,h5in,h5err)
          call read_h5_ex(nid,"tgammauncor",tgammauncor,h5in,h5err)
          call read_h5_ex(nid,"v30lt",v30lt,h5in,h5err)
          call read_h5_ex(nid,"v30rt",v30rt,h5in,h5err)
          call read_h5_ex(nid,"v210lt",v210lt,h5in,h5err)
          call read_h5_ex(nid,"v210rt",v210rt,h5in,h5err)
          call close_group("ins",nid,h5err)
        endif
   
        call test_group(sid,"in_msels",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"in_msels",nid,h5err)
          call read_h5_ex(nid,"bmsels",bmsels,h5in,h5err)
          call read_h5_ex(nid,"sbmsels",sbmsels,h5in,h5err)
          call read_h5_ex(nid,"fwtbmsels",fwtbmsels,h5in,h5err)
          call read_h5_ex(nid,"rrmsels",rrmsels,h5in,h5err)
          call read_h5_ex(nid,"zzmsels",zzmsels,h5in,h5err)
          call read_h5_ex(nid,"l1msels",l1msels,h5in,h5err)
          call read_h5_ex(nid,"l2msels",l2msels,h5in,h5err)
          call read_h5_ex(nid,"l4msels",l4msels,h5in,h5err)
          call read_h5_ex(nid,"emsels",emsels,h5in,h5err)
          call read_h5_ex(nid,"semsels",semsels,h5in,h5err)
          call read_h5_ex(nid,"fwtemsels",fwtemsels,h5in,h5err)
          call read_h5_ex(nid,"kdomsels",kdomsels,h5in,h5err)
          call read_h5_ex(nid,"fmlscut",fmlscut,h5in,h5err)
          call read_h5_ex(nid,"synmsels",synmsels,h5in,h5err)
          call read_h5_ex(nid,"avemsels",avemsels,h5in,h5err)
          call close_group("in_msels",nid,h5err)
        endif

        call test_group(sid,"ina",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"ina",nid,h5err)
          call read_h5_ex(nid,"spatial_avg_gam",spatial_avg_gam,h5in,h5err)
          call close_group("ina",nid,h5err)
        endif

        call test_group(sid,"inece",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"inece",nid,h5err)
          call read_h5_ex(nid,"necein",necein,h5in,h5err)
          call read_h5_ex(nid,"teecein0",teecein0,h5in,h5err)
          call read_h5_ex(nid,"feece0",feece0,h5in,h5err)
          call read_h5_ex(nid,"errorece0",errorece0,h5in,h5err)
          call read_h5_ex(nid,"fwtece0",fwtece0,h5in,h5err)
          call read_h5_ex(nid,"fwtecebz0",fwtecebz0,h5in,h5err)
          call read_h5_ex(nid,"ecefit",ecefit,h5in,h5err)
          call read_h5_ex(nid,"ecebzfit",ecebzfit,h5in,h5err)
          call read_h5_ex(nid,"kfitece",kfitece,h5in,h5err)
          call read_h5_ex(nid,"kinputece",kinputece,h5in,h5err)
          call read_h5_ex(nid,"kcallece",kcallece,h5in,h5err)
          call read_h5_ex(nid,"nharm",nharm,h5in,h5err)
          call read_h5_ex(nid,"kfixro",kfixro,h5in,h5err)
          call read_h5_ex(nid,"rteo",rteo,h5in,h5err)
          call read_h5_ex(nid,"zteo",zteo,h5in,h5err)
          call read_h5_ex(nid,"kfixrece",kfixrece,h5in,h5err)
          call read_h5_ex(nid,"rtep",rtep,h5in,h5err)
          call read_h5_ex(nid,"rtem",rtem,h5in,h5err)
          call read_h5_ex(nid,"rpbit",rpbit,h5in,h5err)
          call read_h5_ex(nid,"rmbit",rmbit,h5in,h5err)
          call read_h5_ex(nid,"robit",robit,h5in,h5err)
          call read_h5_ex(nid,"nfit",nfit,h5in,h5err)
          call read_h5_ex(nid,"kcmin",kcmin,h5in,h5err)
!          call read_h5_ex(nid,"fwtnow",fwtnow,h5in,h5err) ! unused
!          call read_h5_ex(nid,"kdoece",kdoece,h5in,h5err) ! unused
          call read_h5_ex(nid,"mtxece",mtxece,h5in,h5err)
          call read_h5_ex(nid,"nconstr",nconstr,h5in,h5err)
          call read_h5_ex(nid,"eceiter",eceiter,h5in,h5err)
          call read_h5_ex(nid,"eceerror",eceerror,h5in,h5err)
          call close_group("inece",nid,h5err)
        endif

        call test_group(sid,"edgep",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"edgep",nid,h5err)
          call read_h5_ex(nid,"symmetrize",symmetrize,h5in,h5err)
          call read_h5_ex(nid,"rpress",rpress,h5in,h5err)
          call read_h5_ex(nid,"pressr",pressr,h5in,h5err)
          call read_h5_ex(nid,"sigpre",sigpre,h5in,h5err)
          call read_h5_ex(nid,"npress",npress,h5in,h5err)
          call read_h5_ex(nid,"kprfit",kprfit,h5in,h5err)
          call read_h5_ex(nid,"kpressb",kpressb,h5in,h5err)
          call read_h5_ex(nid,"ndokin",ndokin,h5in,h5err)
          call read_h5_ex(nid,"kppfnc",kppfnc,h5in,h5err)
          call read_h5_ex(nid,"kfffnc",kfffnc,h5in,h5err)
          call read_h5_ex(nid,"kffcur",kffcur,h5in,h5err)
          call read_h5_ex(nid,"kppcur",kppcur,h5in,h5err)
          call read_h5_ex(nid,"mxiter",mxiter,h5in,h5err)
          call read_h5_ex(nid,"error",error,h5in,h5err)
          call read_h5_ex(nid,"errmin",errmin,h5in,h5err)
          call read_h5_ex(nid,"keecur",keecur,h5in,h5err)
          call close_group("edgep",nid,h5err)
        endif

        call test_group(sid,"iner",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"iner",nid,h5err)
          call read_h5_ex(nid,"keecur",keecur,h5in,h5err)
          call read_h5_ex(nid,"ecurbd",ecurbd,h5in,h5err)
          call read_h5_ex(nid,"keefnc",keefnc,h5in,h5err)
          call read_h5_ex(nid,"eetens",eetens,h5in,h5err)
          call read_h5_ex(nid,"keebdry",keebdry,h5in,h5err)
          call read_h5_ex(nid,"kee2bdry",kee2bdry,h5in,h5err)
          call read_h5_ex(nid,"eebdry",eebdry,h5in,h5err)
          call read_h5_ex(nid,"ee2bdry",ee2bdry,h5in,h5err)
          call read_h5_ex(nid,"eeknt",eeknt,h5in,h5err)
          call read_h5_ex(nid,"keeknt",keeknt,h5in,h5err)
          call read_h5_ex(nid,"keehord",keehord,h5in,h5err)
          call close_group("iner",nid,h5err)
        endif

        call test_group(sid,"insxr",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"insxr",nid,h5err)
          call read_h5_ex(nid,"ksxr0",ksxr0,h5in,h5err)
          call read_h5_ex(nid,"ksxr2",ksxr2,h5in,h5err)
          call read_h5_ex(nid,"idosxr",idosxr,h5in,h5err)
          call close_group("insxr",nid,h5err)
        endif

        call test_group(sid,"inwant",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"inwant",nid,h5err)
          call read_h5_ex(nid,"psiwant",psiwant,h5in,h5err)
          call read_h5_ex(nid,"vzeroj",vzeroj,h5in,h5err)
          call read_h5_ex(nid,"fwtxxj",fwtxxj,h5in,h5err)
          call read_h5_ex(nid,"fbetap",fbetap,h5in,h5err)
          call read_h5_ex(nid,"fbetan",fbetan,h5in,h5err)
          call read_h5_ex(nid,"fli",fli,h5in,h5err)
          call read_h5_ex(nid,"fqsiw",fqsiw,h5in,h5err)
          call read_h5_ex(nid,"jbeta",jbeta,h5in,h5err)
          call read_h5_ex(nid,"jli",jli,h5in,h5err)
          call read_h5_ex(nid,"alpax",alpax,h5in,h5err)
          call read_h5_ex(nid,"gamax",gamax,h5in,h5err)
          call read_h5_ex(nid,"jwantm",jwantm,h5in,h5err)
          call read_h5_ex(nid,"fwtxxq",fwtxxq,h5in,h5err)
          call read_h5_ex(nid,"fwtxxb",fwtxxb,h5in,h5err)
          call read_h5_ex(nid,"fwtxli",fwtxli,h5in,h5err)
          call read_h5_ex(nid,"znose",znose,h5in,h5err)
          call read_h5_ex(nid,"fwtbdry",fwtbdry,h5in,h5err)
          call read_h5_ex(nid,"nqwant",nqwant,h5in,h5err)
          call read_h5_ex(nid,"siwantq",siwantq,h5in,h5err)
          call read_h5_ex(nid,"kccoils",kccoils,h5in,h5err)
          ! need to ensure array is square
          call read_h5_sq(nid,"ccoils",ccoils,h5in,h5err)
          call read_h5_ex(nid,"rexpan",rexpan,h5in,h5err)
          call read_h5_ex(nid,"xcoils",xcoils,h5in,h5err)
          call read_h5_ex(nid,"kcloops",kcloops,h5in,h5err)
          call read_h5_ex(nid,"cloops",cloops,h5in,h5err)
          call read_h5_ex(nid,"xloops",xloops,h5in,h5err)
          call read_h5_ex(nid,"currc79",currc79,h5in,h5err)
          call read_h5_ex(nid,"currc139",currc139,h5in,h5err)
          call read_h5_ex(nid,"nccoil",nccoil,h5in,h5err)
          call read_h5_ex(nid,"sizeroj",sizeroj,h5in,h5err)
          call read_h5_ex(nid,"fitdelz",fitdelz,h5in,h5err)
          call read_h5_ex(nid,"relaxdz",relaxdz,h5in,h5err)
          call read_h5_ex(nid,"stabdz",stabdz,h5in,h5err)
          call read_h5_ex(nid,"writepc",writepc,h5in,h5err)
          call read_h5_ex(nid,"table_dir",table_dir,h5in,h5err)
          call read_h5_ex(nid,"errdelz",errdelz,h5in,h5err)
          call read_h5_ex(nid,"oldccomp",oldccomp,h5in,h5err)
          call read_h5_ex(nid,"nicoil",nicoil,h5in,h5err)
          call read_h5_ex(nid,"oldcomp",oldcomp,h5in,h5err)
          call read_h5_ex(nid,"currc199",currc199,h5in,h5err)
          call read_h5_ex(nid,"curriu30",curriu30,h5in,h5err)
          call read_h5_ex(nid,"curriu90",curriu90,h5in,h5err)
          call read_h5_ex(nid,"curriu150",curriu150,h5in,h5err)
          call read_h5_ex(nid,"curril30",curril30,h5in,h5err)
          call read_h5_ex(nid,"curril90",curril90,h5in,h5err)
          call read_h5_ex(nid,"curril150",curril150,h5in,h5err)
          call read_h5_ex(nid,"ifitdelz",ifitdelz,h5in,h5err)
          call read_h5_ex(nid,"scaledz",scaledz,h5in,h5err)
          call close_group("inwant",nid,h5err)
        endif
   
        call test_group(sid,"invt",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"invt",nid,h5err)
          call read_h5_ex(nid,"omegat",omegat,h5in,h5err)
          call read_h5_ex(nid,"nomegat",nomegat,h5in,h5err)
          call read_h5_ex(nid,"enw",enw,h5in,h5err)
          call read_h5_ex(nid,"emw",emw,h5in,h5err)
          call read_h5_ex(nid,"betapw0",betapw0,h5in,h5err)
          call read_h5_ex(nid,"kdovt",kdovt,h5in,h5err)
          call read_h5_ex(nid,"kvtor",kvtor,h5in,h5err)
          call read_h5_ex(nid,"kwwcur",kwwcur,h5in,h5err)
          call read_h5_ex(nid,"rvtor",rvtor,h5in,h5err)
          call read_h5_ex(nid,"wcurbd",wcurbd,h5in,h5err)
          call read_h5_ex(nid,"preswb",preswb,h5in,h5err)
          call read_h5_ex(nid,"fwtprw",fwtprw,h5in,h5err)
          call read_h5_ex(nid,"npresw",npresw,h5in,h5err)
          call read_h5_ex(nid,"presw",presw,h5in,h5err)
          call read_h5_ex(nid,"sigprw",sigprw,h5in,h5err)
          call read_h5_ex(nid,"rpresw",rpresw,h5in,h5err)
          call read_h5_ex(nid,"zpresw",zpresw,h5in,h5err)
          call read_h5_ex(nid,"kplotp",kplotp,h5in,h5err)
          call read_h5_ex(nid,"nsplot",nsplot,h5in,h5err)
          call read_h5_ex(nid,"sbetaw",sbetaw,h5in,h5err)
          call read_h5_ex(nid,"comega",comega,h5in,h5err)
          call read_h5_ex(nid,"kcomega",kcomega,h5in,h5err)
          call read_h5_ex(nid,"xomega",xomega,h5in,h5err)
          call read_h5_ex(nid,"romegat",romegat,h5in,h5err)
          call read_h5_ex(nid,"zomegat",zomegat,h5in,h5err)
          call read_h5_ex(nid,"sigome",sigome,h5in,h5err)
          call read_h5_ex(nid,"scalepw",scalepw,h5in,h5err)
          call read_h5_ex(nid,"kwwfnc",kwwfnc,h5in,h5err)
          call read_h5_ex(nid,"kwwknt",kwwknt,h5in,h5err)
          call read_h5_ex(nid,"wwknt",wwknt,h5in,h5err)
          call read_h5_ex(nid,"wwtens",wwtens,h5in,h5err)
          call close_group("invt",nid,h5err)
        endif

        call test_group(sid,"profile_ext",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"profile_ext",nid,h5err)
          call read_h5_ex(nid,"npsi_ext",npsi_ext,h5in,h5err)
          call read_h5_ex(nid,"pprime_ext",pprime_ext,h5in,h5err)
          call read_h5_ex(nid,"ffprim_ext",ffprim_ext,h5in,h5err)
          call read_h5_ex(nid,"psin_ext",psin_ext,h5in,h5err)
          call read_h5_ex(nid,"geqdsk_ext",geqdsk_ext,h5in,h5err)
          call read_h5_ex(nid,"sign_ext",sign_ext,h5in,h5err)
          call read_h5_ex(nid,"scalepp_ext",scalepp_ext,h5in,h5err)
          call read_h5_ex(nid,"scaleffp_ext",scaleffp_ext,h5in,h5err)
          call read_h5_ex(nid,"shape_ext",shape_ext,h5in,h5err)
          call read_h5_ex(nid,"dr_ext",dr_ext,h5in,h5err)
          call read_h5_ex(nid,"dz_ext",dz_ext,h5in,h5err)
          call read_h5_ex(nid,"rc_ext",rc_ext,h5in,h5err)
          call read_h5_ex(nid,"zc_ext",zc_ext,h5in,h5err)
          call read_h5_ex(nid,"a_ext",a_ext,h5in,h5err)
          call read_h5_ex(nid,"eup_ext",eup_ext,h5in,h5err)
          call read_h5_ex(nid,"elow_ext",elow_ext,h5in,h5err)
          call read_h5_ex(nid,"dup_ext",dup_ext,h5in,h5err)
          call read_h5_ex(nid,"dlow_ext",dlow_ext,h5in,h5err)
          call read_h5_ex(nid,"setlim_ext",setlim_ext,h5in,h5err)
          call read_h5_ex(nid,"reflect_ext",reflect_ext,h5in,h5err)
          call read_h5_ex(nid,"fixpp",fixpp,h5in,h5err)
          call close_group("profile_ext",nid,h5err)
        endif

        call test_group(sid,"inlibim",file_stat,h5err)
        if (file_stat) then
          call open_group(sid,"inlibim",nid,h5err)
          call read_h5_ex(nid,"tlibim",tlibim,h5in,h5err)
          call read_h5_ex(nid,"slibim",slibim,h5in,h5err)
          call read_h5_ex(nid,"fwtlib",fwtlib,h5in,h5err)
          call read_h5_ex(nid,"rrrlib",rrrlib,h5in,h5err)
          call read_h5_ex(nid,"zzzlib",zzzlib,h5in,h5err)
          call read_h5_ex(nid,"aa1lib",aa1lib,h5in,h5err)
          call read_h5_ex(nid,"aa8lib",aa8lib,h5in,h5err)
          call close_group("inlibim",nid,h5err)
        endif

        call close_group(trim(tindex),sid,h5err)
        call close_group("time_slice",tid,h5err)
        call close_group("parameters",pid,h5err)
        call close_group("code",cid,h5err)

        ! read previous solution if requested
        if ((geqdsk_ext.ne.'none').or.(icinit.eq.-3).or.(icinit.eq.-4)) then
          call test_group(eqid,"time_slice",file_stat,h5err)
          if (.not. file_stat) then
            call errctrl_msg('data_input','time_slice group not found')
            stop
          endif
          call open_group(eqid,"time_slice",tid,h5err)
          call test_group(tid,trim(tindex),file_stat,h5err)
          if (.not. file_stat) then
            call errctrl_msg('data_input', &
                             trim(tindex)//' group not found')
            stop
          endif
          call open_group(tid,trim(tindex),sid,h5err)
          call test_group(sid,"profiles_2d",file_stat,h5err)
          if (.not. file_stat) then
            call errctrl_msg('data_input','profiles_2d group not found')
            stop
          endif
          call open_group(sid,"profiles_2d",cid,h5err)
          call test_group(cid,"0",file_stat,h5err)
          if (.not. file_stat) then
            call errctrl_msg('data_input','0 group not found')
            stop
          endif
          call open_group(cid,"0",nid,h5err)
          call read_dims(nid,"psi",n2d,h5in,h5err)
          nw_ext=int(n2d(1))
          nh_ext=int(n2d(2))
          allocate(psirz_ext(nw_ext*nh_ext),pprime_ext(nw_ext), &
                   ffprim_ext(nw_ext),qpsi_ext(nw_ext))
          call read_h5_sq(nid,"psi",psirz_ext,h5in,h5err)
          psirz_ext=psirz_ext/twopi
          call close_group("0",nid,h5err)
          call close_group("profiles_2d",cid,h5err)

          ! read in important scalars
          call test_group(sid,"global_quantities",file_stat,h5err)
          if (.not. file_stat) then
            call errctrl_msg('data_input', &
                             'global_quantities group not found')
            stop
          endif
          call open_group(sid,"global_quantities",nid,h5err)
          call read_h5_ex(nid,"psi_axis",simag_ext,h5in,h5err)
          simag_ext=simag_ext/twopi
          call read_h5_ex(nid,"psi_boundary",psibry_ext,h5in,h5err)
          psibry_ext=psibry_ext/twopi
          call read_h5_ex(nid,"ip",plasma_ext,h5in,h5err) ! TODO: could be missing twopi...
          call close_group("global_quantities",nid,h5err)

          ! read in boundary points
          call test_group(sid,"boundary",file_stat,h5err)
          if (.not. file_stat) then
            call errctrl_msg('data_input','boundary group not found')
            stop
          endif
          call open_group(sid,"boundary",cid,h5err)
          call test_group(cid,"outline",file_stat,h5err)
          if (.not. file_stat) then
            call errctrl_msg('data_input','outline group not found')
            stop
          endif
          call open_group(cid,"outline",nid,h5err)
          call read_dims(nid,"r",n1d,h5in,h5err)
          nbdry_ext=int(n1d(1))
          allocate(rbdry_ext(nbdry_ext),zbdry_ext(nbdry_ext))
          call read_h5_ex(nid,"r",rbdry_ext,h5in,h5err)
          call read_h5_ex(nid,"z",zbdry_ext,h5in,h5err)
          call close_group("outline",nid,h5err)
          call close_group("boundary",cid,h5err)

          ! read in p', FF', and q
          call test_group(sid,"profiles_1d",file_stat,h5err)
          if (.not. file_stat) then
            call errctrl_msg('data_input','profiles_1d group not found')
            stop
          endif
          call open_group(sid,"profiles_1d",nid,h5err)
          call read_h5_ex(nid,"dpressure_dpsi",pprime_ext,h5in,h5err)
          pprime_ext=pprime_ext*twopi
          call read_h5_ex(nid,"f_df_dpsi",ffprim_ext,h5in,h5err)
          ffprim_ext=ffprim_ext*twopi
          call read_h5_ex(nid,"q",qpsi_ext,h5in,h5err)
          call close_group("profiles_1d",nid,h5err)

          ! read in fcoil currents (skip ecoils)
          if ((icinit.eq.-3).or.(icinit.eq.-4)) then
            allocate(fcoil_ext(nfsum))
            call test_group(sid,"constraints",file_stat,h5err)
            if (.not. file_stat) then
              call errctrl_msg('data_input', &
                               'constraints group not found')
              stop
            endif
            call open_group(sid,"constraints",cid,h5err)
            call test_group(cid,"pf_current",file_stat,h5err)
            if (.not. file_stat) then
              call errctrl_msg('data_input', &
                               'pf_current group not found')
              stop
            endif
            call open_group(cid,"pf_current",nid,h5err)
            do i=nesum,nesum+nfsum-1
              write(probeind,"(I0)") i
              call test_group(nid,trim(probeind),file_stat,h5err)
              if (.not. file_stat) then
                call errctrl_msg('data_input', &
                                 trim(probeind)//' group not found')
                stop
              endif
              call open_group(nid,trim(probeind),fid,h5err)
              call read_h5_ex(fid,"reconstructed",fcoil_ext(i-nesum+1), &
                              h5in,h5err)
              call close_group(trim(probeind),fid,h5err)
            enddo
            call close_group("pf_current",nid,h5err)
            call close_group("constraints",cid,h5err)
          endif 
          call close_group(trim(tindex),sid,h5err)
          call close_group("time_slice",tid,h5err)
          call close_group("equilibrium",eqid,h5err)

          ! read in limiter position
          call test_group(rootgid,"wall",file_stat,h5err)
          if (.not. file_stat) then
            call errctrl_msg('data_input','wall group not found')
            stop
          endif
          call open_group(rootgid,"wall",eqid,h5err)
          call test_group(eqid,"description_2d",file_stat,h5err)
          if (.not. file_stat) then
            call errctrl_msg('data_input', &
                             'description_2d group not found')
            stop
          endif
          call open_group(eqid,"description_2d",cid,h5err)
          call test_group(cid,"0",file_stat,h5err)
          if (.not. file_stat) then
            call errctrl_msg('data_input','0 group not found')
            stop
          endif
          call open_group(cid,"0",pid,h5err)
          call test_group(pid,"limiter",file_stat,h5err)
          if (.not. file_stat) then
            call errctrl_msg('data_input','limiter group not found')
            stop
          endif
          call open_group(pid,"limiter",tid,h5err)
          call test_group(tid,"unit",file_stat,h5err)
          if (.not. file_stat) then
            call errctrl_msg('data_input','unit group not found')
            stop
          endif
          call open_group(tid,"unit",sid,h5err)
          call test_group(sid,"0",file_stat,h5err)
          if (.not. file_stat) then
            call errctrl_msg('data_input','0 group not found')
            stop
          endif
          call open_group(sid,"0",fid,h5err)
          call test_group(fid,"outline",file_stat,h5err)
          if (.not. file_stat) then
            call errctrl_msg('data_input','outline group not found')
            stop
          endif
          call open_group(fid,"outline",nid,h5err)
          call read_dims(nid,"r",n1d,h5in,h5err)
          limitr_ext = int(n1d(1))
          allocate(xlim_ext(limitr_ext),ylim_ext(limitr_ext))
          call read_h5_ex(nid,"r",xlim_ext,h5in,h5err)
          call read_h5_ex(nid,"z",ylim_ext,h5in,h5err)
          call close_group("outline",nid,h5err)
          call close_group("0",fid,h5err)
          call close_group("unit",sid,h5err)
          call close_group("limiter",tid,h5err)
          call close_group("0",pid,h5err)
          call close_group("description_2d",cid,h5err)
          call close_group("wall",eqid,h5err)

          call open_group(rootgid,"equilibrium",eqid,h5err)
        endif

        call close_group("equilibrium",eqid,h5err)
        call close_h5file(fileid,rootgid,h5err)
#else
        ! this code should not be reachable
        call errctrl_msg('data_input','HDF5 needs to be linked')
        stop
#endif
      endif file_type
!----------------------------------------------------------------------- 
!--   post-process inputs                                             --
!-----------------------------------------------------------------------
      if (link_efit(1:1).ne.'') then
        ! use directories specified in efit.setup if available
        table_dir=trim(link_efit)//'green/'
        input_dir=trim(link_efit)
      elseif (table_dir.ne.table_save .and. link_efit.eq.'') then
        ! error if table_dir changes (global_allocs already set)
        call errctrl_msg('data_input', &
          'changing machine during run is not supported (table_dir)')
        kerror=1
        return
      endif
      if(link_store(1:1).ne.'') store_dir=trim(link_store)
      ! replace parameters with in0 namelist versions
      if (ierchk_prior .ne. 99) then
        ! these parameters must be set together because logicals cannot
        ! be overloaded...
        ierchk = ierchk_prior
        req_valid = req_valid_prior
      endif
      if (iout_prior .ne. -1) then
        iout = iout_prior
      endif
      if (iplcout_prior .ne. -1) then
        iplcout = iplcout_prior
      endif
      
!--   warn that idebug, jdebug, and ktear inputs are deprecated
      if (idebug.ne.0) write(*,*) &
      "idebug input variable is deprecated, set cmake variable instead"
      if (jdebug.ne."NONE") write(*,*) &
      "jdebug input variable is deprecated, set cmake variable instead"
      if(ktear.ne.0) write(*,*) &
      "tearing calculations don't exist, ktear is deprecated"
!--   roundoff differences can throw off zlim if limiter corners
!--   are too close to grid points (maybe zlim needs fixing...)
      do i=1,limitr
        ylim(i)=ylim(i)-1.e-10_dp
      enddo
!--   protect against underflow in fitting weights 
      if(abs(fwtdlc).le.1.e-30_dp) fwtdlc=0.0
      if(abs(fwtcur).le.1.e-30_dp) fwtcur=0.0
      do i=1,nfsum
        if(abs(fwtfc(i)).le.1.e-30_dp) fwtfc(i)=0.0
      enddo
      do i=1,nesum
        if(abs(fwtec(i)).le.1.e-30_dp) fwtec(i)=0.0
      enddo
      do i=1,magpri
        if(abs(fwtmp2(i)).le.1.e-30_dp) fwtmp2(i)=0.0
      enddo
      do i=1,nsilop
        if(abs(fwtsi(i)).le.1.e-30_dp) fwtsi(i)=0.0
      enddo
      do i=1,nmselp
        if(abs(fwtgam(i)).le.1.e-30_dp) fwtgam(i)=0.0
      enddo
      do i=1,nmsels
        if(abs(fwtbmsels(i)).le.1.e-30_dp) fwtbmsels(i)=0.0
        if(abs(fwtemsels(i)).le.1.e-30_dp) fwtemsels(i)=0.0
      enddo
      do i=1,nnece
        if(abs(fwtece0(i)).le.1.e-30_dp) fwtece0(i)=0.0
      enddo
      if(abs(fwtecebz0).le.1.e-30_dp) fwtecebz0=0.0
!
      if(nbdryp==-1) nbdryp=nbdry

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
            kerror=1
            return
          endif
        endif 
      endif 
91008 format(9e12.5,i2) 

!----------------------------------------------------------------------
!--   Setup FF', P' arrays
!----------------------------------------------------------------------
      read_geqdsk: if (geqdsk_ext.ne.'none'.and.icinit.ne.-3 &
                                           .and.icinit.ne.-4) then
        npsi_ext=nw_ext 
#ifdef DEBUG_LEVEL1
        write (nttyo,*) 'npsi_ext,nw_ext=',npsi_ext,nw_ext 
#endif
        if (plasma_ext > 0.0) then 
          sign_ext = -1.0 
        endif 
        write_boundary: if (nbdry.le.0) then 
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
          fwtbdry(1:nbdry)=1. 
          do i=1,nbdry
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
        endif write_boundary
! 
        if (limitr.le.0) then
          xlim(1:limitr_ext)=xlim_ext(1:limitr_ext)
          ylim(1:limitr_ext)=ylim_ext(1:limitr_ext)-1.e-10_dp
        endif
      endif read_geqdsk
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
        loc=minloc(rbdry0(1:nbdry),1)
        r0min=rbdry0(loc)
        zr0min=zbdry0(loc)
        loc=maxloc(rbdry0(1:nbdry),1)
        r0max=rbdry0(loc)
        zr0max=zbdry0(loc)
        loc=minloc(zbdry0(1:nbdry),1)
        z0min=zbdry0(loc)
        rz0min=rbdry0(loc)
        loc=maxloc(zbdry0(1:nbdry),1)
        z0max=zbdry0(loc)
        rz0max=rbdry0(loc)
        r0ave=0.5*(r0min+r0max) 
        z0ave=0.5*(zr0min+zr0max) 
        a0ave=0.5*(r0max-r0min) 
        e0top=(z0max-z0ave)/a0ave 
        e0bot=(z0ave-z0min)/a0ave 
        d0top=(r0ave-rz0max)/a0ave 
        d0bot=(r0ave-rz0min)/a0ave 
        if(rc_ext.le.-10.0) rc_ext=r0ave
        if(zc_ext.le.-10.0) zc_ext=z0ave
        if(a_ext.le.-10.0) a_ext=a0ave
        if(eup_ext.le.-10.0) eup_ext=e0top
        if(elow_ext.le.-10.0) elow_ext=e0bot
        if(dup_ext.le.-10.0) dup_ext=d0top
        if(dlow_ext.le.-10.0) dlow_ext=d0bot
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
        loc=minloc(rbdry0(1:nbdry),1)
        r0min=rbdry0(loc)
        zr0min=zbdry0(loc)
        loc=maxloc(rbdry0(1:nbdry),1)
        r0max=rbdry0(loc)
        zr0max=zbdry0(loc)
        loc=minloc(zbdry0(1:nbdry),1)
        z0min=zbdry0(loc)
        rz0min=rbdry0(loc)
        loc=maxloc(zbdry0(1:nbdry),1)
        z0max=zbdry0(loc)
        rz0max=rbdry0(loc)
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
      call set_basis_params
      if (kedgep.gt.0) then
        s1edge=(1.0-pe_psin)/pe_width
        tpedge=tanh(s1edge)
        s1edge=(1.0-fe_psin)/fe_width
        tfedge=tanh(s1edge)
      endif
#ifdef DEBUG_LEVEL1
      write(*,*)  'before save fitting weights'
#endif
!--------------------------------------------------------------------- 
!--   save fitting weights
!--------------------------------------------------------------------- 
      swtdlc=fwtdlc
      swtcur=fwtcur
      swtfc=fwtfc
      swtec=fwtec
      swtmp2=fwtmp2
      swtsi=fwtsi
      swtgam(1:nmselp)=fwtgam(1:nmselp)
      swtbmsels=fwtbmsels
      swtemsels=fwtemsels
      swtece=fwtece0
      swtecebz=fwtecebz0
      iexcals=iexcal
      ibunmns=ibunmn
      ierchks=ierchk
#ifdef DEBUG_LEVEL1
      write(*,*)'adjust fit parameters based on basis function selected'
#endif
      if(fbetan.gt.0.) brsp(nfsum+jbeta)=alpax(jbeta)*darea
      if(fli.gt.0.) brsp(nfsum+kppcur+jli)=gamax(jli)*darea
      if(kedgep.gt.0) pedge=pedge*darea
      if(kedgef.gt.0) f2edge=f2edge*darea
      if (npress.lt.0) then
        kdopre=-npress
        npress=0
      endif
      if (itek.lt.0) then
        ixray=1
        itek=iabs(itek)
      endif
      if(psiwant.le.0.0) psiwant=1.e-5_dp

#ifdef DEBUG_LEVEL1
      write(*,*) 'itek > 100, write out PLTOUT.OUT individually '
#endif
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
        if(fitdelz) itell=4 
      endif 
      if (nxiter.lt.0) then 
        nxiter=-nxiter 
        itell=1 
        if(fbetan.gt.0.0) itell=2
        if(fli.gt.0.0) itell=2
        if(nqwant.gt.0.0) itell=2
        if((symmetrize).and.(nbdry.gt.1)) itell=3 
      endif 
      if((iconvr.ne.3).and.(qvfit.gt.0.0)) qenp=qvfit 
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
        ! save input values before being overwritten
        premea(1:npress)=pressr(1:npress)
        premew(1:npresw)=presw(1:npresw)
      endif
      if (kprfit.gt.0.and.sigpre(1).lt.0.0) then 
        scalep=abs(sigpre(1)) 
        scalemin=abs(sigpre(2)) 
        sigpre(1:npress)=scalep*pressr(1:npress) 
        do i=1,npress 
          sigpre(i)=max(sigpre(i),scalemin) 
        enddo
      endif 
      if (kprfit.gt.0.and.scalepr(1).gt.0.0) then 
        pressr(1:npress)=pressr(1:npress)*scalepr(1:npress) 
        sigpre(1:npress)=sigpre(1:npress)*scalepr(1:npress) 
      endif 
      if (kprfit.gt.0.and.scalepw(1).gt.0.0) then 
        if (npresw.gt.0) then 
          presw(1:npresw)=presw(1:npresw)*scalepw(1:npresw) 
          sigprw(1:npresw)=sigprw(1:npresw)*scalepw(1:npresw) 
        elseif (nomegat.gt.0) then 
          omegat(1:nomegat)=omegat(1:nomegat)*scalepw(1:nomegat) 
          sigome(1:nomegat)=sigome(1:nomegat)*scalepw(1:nomegat) 
        endif 
      endif
!-------------------------------------------------------------- 
!--   option to symmetrize added 8/14/91 eal                 -- 
!-------------------------------------------------------------- 
      sym: if (symmetrize) then
#ifdef DEBUG_LEVEL1
        write(*,*) 'option to symmetrize added 8/14/91 eal  '
#endif 
        isetfb=0  ! be sure vertical feedback is off 
        zelip=0 ! must be symmetric about midplane 
        fixed_bdry_1: if (nbdry.gt.1) then  ! remove duplicate point 
          nbryup=0 
          delx2=(rbdry(1)-rbdry(nbdry))**2 
          dely2=(zbdry(1)-zbdry(nbdry))**2 
          if((delx2+dely2).le.1.0e-08_dp) nbdry=nbdry-1 
          rbar=sum(rbdry(1:nbdry))/nbdry
          do i=1,nbdry 
            ilower(i)=-1 
            if (zbdry(i).gt.0.) then
              aup=atan2(zbdry(i),rbdry(i)-rbar) 
              iup=i 
              near=1e30 
              do j=1,nbdry 
                if (zbdry(j).lt.0.) then
                  adn=atan2(zbdry(j),rbdry(j)-rbar) 
                  val=abs(aup+adn) 
                  if (val.lt.near) then
                    near=val 
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
        endif fixed_bdry_1 
      endif sym
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
      if(kpressb.eq.2) pcurbd=0.0 
      if (kzeroj.gt.0) then 
        pcurbd=0.0 
        fcurbd=0.0 
      endif 
      close(unit=nin) 
      kinetic: if (kprfit.eq.1) then 
        if (npress.lt.0) then 
          call setfnmeq(itimeu,'k',ishot,itime,edatname) 
          edatname='edat_'//edatname(2:7)// & 
                       '_'//edatname(9:13)//'.pressure' 
          open(unit=nin,status='old',file=edatname)
          read (nin,edat)
          close(unit=nin) 
        endif
      elseif (kprfit.eq.2) then kinetic
        if (npteth.lt.0) then 
          nptef=-npteth 
          npnef=-npneth 
          call setfnmeq(itimeu,'k',ishot,itime,edatname) 
          edatname='edat_'//edatname(2:7)// & 
                       '_'//edatname(9:13)//'.thomson' 
          open(unit=nin,status='old',file=edatname)
          bfract=-1. 
          if (tethom(1).lt.0.0) bfract=-tethom(1) 
          read (nin,edat) 
          close(unit=nin) 
        endif 
        if(nbeam.gt.0) dnbeam(1:nbeam)=dnbeam(1:nbeam)*1.e-19_dp

#ifdef DEBUG_LEVEL1
        write(*,*) 'reorder TS data points '
#endif
!--------------------------------------------------------------------- 
!--     reorder TS data points                                      -- 
!--------------------------------------------------------------------- 
        call tsorder(npteth,zteth,dnethom,tethom,sgneth,sgteth)
        if(sgnemin.lt.0.0) sgnemin=abs(sgnemin)*dnethom(1)*1.e-19_dp & 
                                                *co2cor
        dnethom(1:npneth)=dnethom(1:npneth)*1.e-19_dp*co2cor 
        sgneth(1:npneth)=sgneth(1:npneth)*1.e-19_dp*sgnethi*co2cor 
        do i=1,npneth
          sgneth(i)=max(sgneth(i),sgnemin) 
        enddo
        if (sgtemin.lt.0.0) sgtemin=abs(sgtemin)*tethom(1) 
        temax=tethom(1) 
        demax=dnethom(1) 
        sgteth(1:npneth)=sgteth(1:npneth)*sgtethi 
        if (bfract.gt.0.0) then 
          tethom(1:npneth)=tethom(1:npneth)*bfract 
          sgteth(1:npneth)=sgteth(1:npneth)*bfract 
        endif 
        do i=1,npteth 
          sgteth(i)=max(sgteth(i),sgtemin) 
        enddo
        temax=max(temax,maxval(tethom(1:npneth))) 
        demax=max(demax,maxval(dnethom(1:npneth))) 
        if(cstabte.lt.0.0) cstabte=abs(cstabte)*100./temax
        if(cstabne.lt.0.0) cstabne=abs(cstabne)*100./demax
        if (nption.lt.0) then
          nptionf=-nption
          if (nptionf.lt.100) then
            call setfnmeq(itimeu,'k',ishot,itime,edatname)
            edatname='edat_'//edatname(2:7)//'_'//edatname(9:13)//'.cer'
            open(unit=nin,status='old',file=edatname)
            bfract=-1.
            if(tionex(1).lt.0.0) bfract=-tionex(1)
            read (nin,edat)
            close(unit=nin)
            sigti(1:nption)=sigti(1:nption)*sigtii
            if (bfract.gt.0.0) then
              sigti(1:nption)=sigti(1:nption)*bfract
              tionex(1:nption)=tionex(1:nption)*bfract
            endif
          endif
          if (nptionf.gt.100) then
            nptionf=nptionf-100
            nption=npteth
            nptionf=nptef
            bfract=tionex(1)
            sigti(1:nption)=sgteth(1:nption)
            tionex(1:nption)=tethom(1:nption)*bfract
            rion(1:nption)=rteth(1:nption)
            zion(1:nption)=zteth(1:nption)
          endif
        endif
      endif kinetic
!----------------------------------------------------------------------- 
!---- read in limiter data                                            -- 
!----------------------------------------------------------------------- 
#ifdef DEBUG_LEVEL1
      write(*,*) 'read in limiter data'
#endif
      ! this is only needed for jtime=1 unless profile_ext variables are
      ! used... (can a check for this be added?)
      call getlim(1,xltype,xltype_180,shape_ext)
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
        if(kwaitmse.ne.0) write (neqdsk,ina)
        write (nout,inece)
        write (nout,insxr)
        write (nout,invt)
        write (nout,inlibim)
      endif
      if(islve.gt.0) nbdry=40
! 
      diamag(jtime)=1.0e-03_dp*dflux
      sigdia(jtime)=1.0e-03_dp*abs(sigdlc)
      pbinj(jtime)=pnbeam
      silopt(jtime,:)=coils
      fccurt(jtime,:)=brsp(1:nfsum)
      tangam(jtime,1:nmselp)=tgamma 
      tangam_uncor(jtime,1:nmselp)=tgammauncor 
      siggam(jtime,1:nmselp)=sgamma 
      rrgam(jtime,1:nmselp)=rrrgam 
      zzgam(jtime,1:nmselp)=zzzgam 
      a1gam(jtime,1:nmselp)=aa1gam 
      a2gam(jtime,1:nmselp)=aa2gam 
      a3gam(jtime,1:nmselp)=aa3gam 
      a4gam(jtime,1:nmselp)=aa4gam 
      a5gam(jtime,1:nmselp)=aa5gam 
      a6gam(jtime,1:nmselp)=aa6gam 
      a7gam(jtime,1:nmselp)=aa7gam 
      a8gam(jtime,1:nmselp)=0.0 
      tangam(jtime,nmselp+1:nstark)=tlibim
      siggam(jtime,nmselp+1:nstark)=slibim
      rrgam(jtime,nmselp+1:nstark)=rrrlib
      zzgam(jtime,nmselp+1:nstark)=zzzlib
      a1gam(jtime,nmselp+1:nstark)=aa1lib
      a2gam(jtime,nmselp+1:nstark)=0.0
      a3gam(jtime,nmselp+1:nstark)=0.0
      a4gam(jtime,nmselp+1:nstark)=0.0
      a5gam(jtime,nmselp+1:nstark)=0.0
      a6gam(jtime,nmselp+1:nstark)=0.0
      a7gam(jtime,nmselp+1:nstark)=0.0
      a8gam(jtime,nmselp+1:nstark)=aa8lib
      fwtgam(nmselp+1:nstark)=fwtlib
      swtgam(nmselp+1:nstark)=fwtlib
! 
      bmselt(jtime,:)=bmsels
      sbmselt(jtime,:)=sbmsels
      fwtbmselt(jtime,:)=fwtbmsels
      rrmselt(jtime,:)=rrmsels
      zzmselt(jtime,:)=zzmsels
      l1mselt(jtime,:)=l1msels
      l2mselt(jtime,:)=l2msels
      l3mselt(jtime,:)=1.-l1msels
      l4mselt(jtime,:)=l4msels
      emselt(jtime,:)=emsels
      semselt(jtime,:)=semsels
      fwtemselt(jtime,:)=fwtemsels
!---------------------------------------------------------------------- 
!     give the constraint value for matrix routine 
!---------------------------------------------------------------------- 
      brspece(jtime,:)=ecefit
      brspecebz(jtime)=ecebzfit
      expmpi(jtime,:)=expmp2
!------------------------------------------------------------------------ 
!--   New E-coil connections                   LLao, 95/07/11          --
!------------------------------------------------------------------------ 
      if (nesum.ge.3) then
        if(ecurrt(3).le.-1.e10_dp) ecurrt(3)=ecurrt(1)
      endif
      if (nesum.ge.5) then
        if(ecurrt(5).le.-1.e10_dp) ecurrt(5)=ecurrt(1) 
      endif
      if (nesum.ge.4) then
        if(ecurrt(4).le.-1.e10_dp) ecurrt(4)=ecurrt(2) 
      endif
      if (nesum.ge.6) then
        if(ecurrt(6).le.-1.e10_dp) ecurrt(6)=ecurrt(2)
      endif
      eccurt(jtime,:)=ecurrt
      ipmeas(jtime)=plasma
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
      denvt(jtime,:)=denv 
      denrt(jtime,:)=denr
      accurt(jtime,:)=acoilc
      caccurt(jtime,:)=acoilc
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

      ! need to match snap mode saved variables
      zelipss=zelip

      endif snap
!----------------------------------------------------------------------- 
!--   set up parameters for all modes                                 --
!----------------------------------------------------------------------- 
      ! setup error calls
      call errctrl_setstate(rank,time(jtime))

      ! restore previous solution
      if(jtime.gt.1 .and. & 
         (use_previous .or. (isicinit.lt.0 .and. isicinit.ne.-3))) &
        brsp = brsp_save

      ! legacy options...
      ibunmn=ibunmns
      if(ibunmn.eq.3) ibunmn=1
      if(ibunmn.eq.4) ibunmn=2

      ! check grid subsampling
      if (abs(nw_sub).gt.nw) then
        write(*,*)  &
          "Higher resolution output request than run, ignoring nw_sub"
        nw_sub=nw
      elseif (nw_sub.le.0) then
        write(*,*)  &
          "Negative resolution output not possible, ignoring nw_sub"
        nw_sub=nw
      else
        nsw=(nw-1)/(nw_sub-1)
        nw_sub=(nw-1)/nsw+1
      endif
      if (abs(nh_sub).gt.nh) then
        write(*,*)  &
          "Higher resolution output request than run, ignoring nh_sub"
        nh_sub=nh
      elseif (nh_sub.le.0) then
        write(*,*)  &
          "Negative resolution output not possible, ignoring nh_sub"
        nh_sub=nh
      else
        nsh=(nh-1)/(nh_sub-1)
        nh_sub=(nh-1)/nsh+1
      endif

      if(kfffnc.eq.8) rkec=pi/(2.0*dpsiecn) 
      chigam=0.0 
      tchimls=0.0 
!
      if (ipmeas(jtime).le.-1.e3_dp) then 
        negcur=1 
      else 
        negcur=0 
      endif 
      iexcal=iexcals
      ivacum=0
      ierchk=ierchks
      if (abs(ipmeas(jtime)).le.cutip.and.iconvr.ge.0) then
        if(iconsi.eq.-1) iconsi=55 
        iexcal=1 
        ivacum=1 
        ibunmn=0 
        ierchk=0 
        nxiter=1 
        iconvr=3 
      endif

! Check that the grid sizes are compatible with cyclic reduction
! (single or double)
      if (ibunmn.ne.0) then
        if (nh.ne.0) then
          select case (nh)
          case (3,5,9,17,33,65,129,257,513,1025,2049)
            ! all good
          case default
            call errctrl_msg('data_input', &
                 'Chosen grid dimensions cannot be run')
            stop
          end select
        endif
        if (nh.eq.0 .or. isolve.eq.0) then
          select case (nw)
          case (3,5,9,17,33,65,129,257,513,1025,2049)
            ! all good
          case default
            call errctrl_msg('data_input', &
                 'Chosen grid dimensions cannot be run')
            stop
          end select
        endif
      endif
!--------------------------------------------------------------------- 
!--   DIIID correction to 322 degree probes due to N1 coil
!--------------------------------------------------------------------- 
      if (oldccomp) then 
        if (n1coil.eq.2.and.ishot.le.108281) then 
          open(unit=60,file=input_dir(1:lindir)//'n1coil.ddd', & 
               status='old')
          j=jtime
          do i=30,60 
            read(60,*) namedum,xxxdum,signn1(i) 
            expmpi(j,i)=expmpi(j,i)-signn1(i)*curtn1(j) 
          enddo
          close(unit=60) 
        endif 
      endif 
!--------------------------------------------------------------------- 
!--   DIIID correction to 322 and 67 degree probes due to C coil
!--------------------------------------------------------------------- 
      if (nccoil.eq.1.and.oldccomp) then 
        open(unit=60,file=input_dir(1:lindir)//'ccoil.ddd', & 
             status='old') 
        j=jtime
        do i=30,60
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
      if(ifitvs.eq.1) ivesel=3 
      if (.not.fitsiref) then 
        if (iecurr.gt.0.or.nslref.lt.0) then 
          silopt(jtime,1:nsilop)=silopt(jtime,1:nsilop)+psiref(jtime) 
          psiref(jtime)=0. 
        endif 
      endif 
      kconvr=iconvr
      www=zero 
!---------------------------------------------------------------------- 
!--   signal at psi loop # NSLREF is used as reference               -- 
!---------------------------------------------------------------------- 
      fwtref=fwtsi(iabs(nslref)) 
      kersil_23: if ((kersil.ne.2).and.(kersil.ne.3)) then
      do m=1,nsilop
        tdata1=serror*abs(silopt(jtime,m)) 
        tdata2=abs(psibit(m))*vbit 
        tdata=max(tdata1,tdata2) 
        sigsil(m)=tdata 
        if (tdata.gt.1.0e-10_dp) then
          fwtsi(m)=fwtsi(m)/tdata**nsq
        else
          fwtsi(m)=0.0
        endif
      enddo
      if (abs(psibit(iabs(nslref))).le.1.0e-10_dp) then
        coilmx=maxval(abs(silopt(jtime,1:nsilop))) 
        sigsil(iabs(nslref))=coilmx*serror
        fwtsi(iabs(nslref))=1.0/coilmx**nsq/serror**nsq*fwtref
      endif 
      else kersil_23
!----------------------------------------------------------------------- 
!--   Fourier expansion of vessel sgments                             -- 
!----------------------------------------------------------------------- 
      if (ivesel.eq.3 .and. nfourier.gt.1) then
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
            if(i.eq.1) vecta(i,j)=1.0 
            if(i.gt.1.and.i.le.(nfourier+1)) vecta(i,j)=costa(i-1,j) 
            if(i.gt.(nfourier+1)) vecta(i,j)=sinta(i-nfourier-1,j) 
          enddo 
        enddo 
      endif
! 
!     if(brsptu(1).le.-1.e-20_dp) & 
!       brsp(1:nfsum)=brsptu(1:nfsum)*turnfc(1:nfsum) 
      if((brsptu(1).gt.-1.e-20_dp).and.(icinit.ne.-3).and.(icinit.ne.-4)) &
        brsp(1:nfsum)=brsptu(1:nfsum)*turnfc(1:nfsum) 
      reflux=silopt(jtime,iabs(nslref)) 
      do m=1,nsilop 
        tdata1=errsil*abs(silopt(jtime,m)-reflux) 
        tdata2=sicont*rsi(m)*abs(ipmeas(jtime)) 
        tdata=max(tdata1,tdata2) 
        tdata2=abs(psibit(m))*vbit 
        tdata=max(tdata,tdata2) 
        sigsil(m)=tdata 
        if (tdata.gt.1.0e-10_dp) then
          fwtsi(m)=fwtsi(m)/tdata**nsq
        else
          fwtsi(m)=0.0
        endif
      enddo
      if (kersil.ne.3) then
!---------------------------------------------------------------------- 
!--     signal at psi loop #8 in D-III is set to zero and used as ref
!---------------------------------------------------------------------- 
        sigsil(iabs(nslref))=ersil8
        fwtsi(iabs(nslref))=1.0/ersil8**nsq*fwtref
      else
!---------------------------------------------------------------------- 
!--     Default option for reference flux loop uncertainty
!---------------------------------------------------------------------- 
        m=iabs(nslref)
        tdata1=errsil*abs(silopt(jtime,m))
        tdata2=sicont*rsi(m)*abs(ipmeas(jtime))
        tdata=max(tdata1,tdata2)
        tdata2=abs(psibit(m))*vbit
        tdata=max(tdata,tdata2)
        sigsil(m)=tdata
        if (tdata.gt.1.0e-10_dp) then
          fwtsi(m)=fwtref/tdata**nsq
        else
          fwtsi(m)=0.0
        endif
      endif
      endif kersil_23
      sigref=sigsil(iabs(nslref)) 
      fwtref=fwtsi(iabs(nslref)) 
! 
      do m=1,magpri 
        tdata1=serror*abs(expmpi(jtime,m)) 
        tdata2=abs(bitmpi(m))*vbit 
        tdata=max(tdata1,tdata2)
        sigmpi(m)=tdata
        if (tdata.gt.1.0e-10_dp) then
          fwtmp2(m)=fwtmp2(m)/tdata**nsq
        else
          fwtmp2(m)=0.0
        endif
      enddo
      do m=1,nstark 
        tdata=abs(siggam(jtime,m))
        if (tdata.gt.1.0e-10_dp) then
          fwtgam(m)=fwtgam(m)/tdata**nsq
        else
          fwtgam(m)=0.0
        endif
      enddo
! 
#ifdef DEBUG_LEVEL2
      write (6,*) 'DATA fwtbmselt = ',(fwtbmselt(jtime,i),i=1,nmsels) 
      write (6,*) 'DATA sbmselt = ',(sbmselt(jtime,i),i=1,nmsels) 
#endif 
      do i=1,nmsels 
        tdata=abs(sbmselt(jtime,i)) 
        if (tdata.gt.1.0e-10_dp) then
          fwtbmselt(jtime,i)=fwtbmselt(jtime,i)/tdata**nsq
        else
          fwtbmselt(jtime,i)=0.0
        endif
      enddo
      do i=1,nmsels 
        tdata=abs(semselt(jtime,i)) 
        if (tdata.gt.1.0e-10_dp) then
          fwtemselt(jtime,i)=fwtemselt(jtime,i)/tdata**nsq
        else
          fwtemselt(jtime,i)=0.0
        endif
      enddo 
#ifdef DEBUG_LEVEL2
      write (6,*) 'DATA fwtbmselt = ', (fwtbmselt(jtime,i),i=1,nmsels)
#endif
! 
      do m=1,nfsum 
        tdata1=serror*abs(fccurt(jtime,m)) 
        tdata2=abs(bitfc(m))*vbit 
        tdata=max(tdata1,tdata2) 
        sigfcc(m)=tdata
        if (tdata.gt.1.0e-10_dp) then
          fwtfc(m)=fwtfc(m)/tdata**nsq
        else
          fwtfc(m)=0.0
        endif
      enddo
      ecurrt=eccurt(jtime,:) 
      cecurr=ecurrt 
      if (iecurr.eq.2) then 
        do m=1,nesum 
          tdata1=serror*abs(ecurrt(m)) 
          tdata2=abs(bitec(m))*vbit 
          tdata=max(tdata1,tdata2) 
          sigecc(m)=tdata 
          if (tdata.gt.1.0e-10_dp) then
            fwtec(m)=fwtec(m)/tdata**nsq
          else
            fwtec(m)=0.0
          endif
        enddo 
      endif 
      tdata1=serror*abs(ipmeas(jtime)) 
      tdata2=abs(bitip)*vbit 
      tdata=max(tdata1,tdata2) 
      sigpasma=tdata
      if (tdata.gt.1.0e-10_dp) then 
        fwtcur=fwtcur/tdata**nsq
      else
        fwtcur=0.0
      endif
!---------------------------------------------------------------------- 
!--   diamagetic flux                                                -- 
!---------------------------------------------------------------------- 
      tdata=abs(sigdia(jtime)) 
      if(tdata.gt.1.0e-10_dp) fwtdlc=fwtdlc/tdata**nsq 
! 
      if(sidif.le.-1.0e+10_dp) sidif=tmu*ipmeas(jtime)*rcentr/2.0 
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
      nfnpcr=nfsum+kpcurn 
      nbase=nfsum+kppcur 
      nfnwcr=nfnpcr 
      if (kvtor.gt.0) then 
        nfnwcr=nfnwcr+kwwcur 
        kwcurn=kpcurn+kwwcur 
      else 
        kwcurn=kpcurn 
      endif 
      nqaxis=0 
      if(fwtqa.gt.1.0e-03_dp) nqaxis=1 
      nparam=nfnwcr 
      if(kprfit.gt.0) nparam=nparam+1 
      if(fitdelz) nparam=nparam+1 
      if(fitsiref) nparam=nparam+1 
      if(kedgep.gt.0) nparam=nparam+1 
      if(kedgef.gt.0) nparam=nparam+1 
      if(fwtqa.gt.0.0) fwtqa=fwtqa/errorq 
      if(fwtbp.gt.0.0) fwtbp=fwtbp/errorq 
      if(fbetap.gt.0.0) betap0=fbetap 
! 
      ipsi(jtime)=0 
      do i=1,nsilop 
        if(fwtsi(i).gt.0.0) ipsi(jtime)=ipsi(jtime)+1 
      enddo
      ifc(jtime)=0 
      do i=1,nfsum 
        if(fwtfc(i).gt.0.0) ifc(jtime)=ifc(jtime)+1 
      enddo
      iec(jtime)=0 
      do i=1,nesum 
        if(fwtec(i).gt.0.0) iec(jtime)=iec(jtime)+1 
      enddo 
      imag2(jtime)=0 
      do i=1,magpri 
        if(fwtmp2(i).gt.0.0) imag2(jtime)=imag2(jtime)+1 
      enddo
      kmtark=0 
      klibim=0 
      do i=1,nmselp 
        if(fwtgam(i).gt.0.0) kmtark=kmtark+1 
      enddo
      do i=nmselp+1,nstark 
        if(fwtgam(i).gt.0.0) klibim=klibim+1 
      enddo 
      kstark=kmtark+klibim 
! 
      mmbmsels=0 
      mmemsels=0 
      do i=1,nmsels 
        if(fwtbmselt(jtime,i).gt.1.e-06_dp) mmbmsels=mmbmsels+1 
        if(fwtemselt(jtime,i).gt.1.e-06_dp) mmemsels=mmemsels+1 
      enddo 
#ifdef DEBUG_MSELS
      write (6,*) 'DATA mmbmsels = ', mmbmsels 
      write (6,*) 'bmselt ',(bmselt(1,i),i=1,nmsels) 
      write (6,*) 'iermselt ',(iermselt(1,i),i=1,nmsels) 
#endif 
! 
      iplasm(jtime)=0 
      if(fwtcur.gt.0.0) iplasm(jtime)=1 
      idlopc(jtime)=0 
      if(fwtdlc.gt.0.0) idlopc(jtime)=1 
      ipmhd(jtime)=ipmeas(jtime) 
      if (iconvr.eq.3.and.ivesel.eq.1) then 
        do i=1,nvesel 
          ipmhd(jtime)=ipmhd(jtime)-vcurrt(i) 
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
!      if ((ivesel.eq.0).or.(idovs.gt.0)) go to 525 
!      open(unit=nrsppc,status='old',form='unformatted', & 
!           file=table_dir(1:ltbdir)//'rv'//trim(ch1)//trim(ch2)//'.ddd') 
!      read (nrsppc) rsilvs 
!      read (nrsppc) rmp2vs 
!      read (nrsppc) gridvs 
!      close(unit=nrsppc) 
!      idovs=1 
!      go to  525 
!      open(unit=nffile,status='old',form='unformatted', & 
!           file=table_dir(1:ltbdir)//'fc'//trim(ch1)//trim(ch2)//'.ddd') 
!      read (nffile) rfcfc 
!      close(unit=nffile) 
!      endif 
! 
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
      if(kstark.gt.0.or.kdomse.gt.0) call setstark(jtime) 
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
      do k=1,necein 
        kk=necein-k+1 
        if(kfitece.eq.3) kk = k 
        feece(kk)=feece0(k) 
        teecein(kk)=teecein0(k) 
        errorece(kk)=errorece0(k) 
      enddo 
!--------------------------------------------------------------------- 
!--   toroidal rotation - set up geometric parameters               -- 
!--------------------------------------------------------------------- 
      if (kvtor.gt.0) then 
        do i=1,nw 
          rgrvt(i)=(rgrid(i)/rvtor)**2 
          rgrvt(i)=rgrvt(i)-1. 
          rgsvt(i)=rgrid(i)*rgrvt(i) 
        enddo 
      endif 
!---------------------------------------------------------------------- 
!--   make filament Green's tables only                              -- 
!---------------------------------------------------------------------- 
      solution_mode: if ((iconvr.lt.0).and.(iconvr.gt.-20)) then
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
          if(k.gt.(mx/2+1)) ix=2 
          i=(rmx(k)-rgrid(1))/drgrid+1 
          j=(zmx(k)-zgrid(1))/dzgrid+ix 
          zdif=zmx(k)-zgrid(j) 
          if (abs(zdif).gt.0.6_dp*dzgrid) then 
            if(zdif.gt.0.0) j=j+1 
            if(zdif.lt.0.0) j=j-1 
          endif 
        endif 
        irfila(k)=i 
        jzfila(k)=j 
        rmx(k)=rgrid(i) 
        zmx(k)=zgrid(j) 
        kk=(i-1)*nh+j 
        rsilpf(:,k)=gsilpc(:,kk) 
        rmp2pf(:,k)=gmp2pc(:,kk) 
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
           file='rpfxx.dat',iostat=istat)
      if(istat.eq.0) close(unit=nffile,status='delete')
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
      elseif (iconvr.le.-20) then solution_mode
      if (aelip.gt.0.0) then 
        do i=1,nw 
          do j=1,nh 
            kk=(i-1)*nh+j 
            erho=sqrt((rgrid(i)-relip)**2+((zgrid(j)-zelip)/eelip)**2) 
            xpsi(kk)=(erho/aelip)**2 
          enddo
        enddo
      endif 
      wsilpc=0.0 
      wmp2pc=0.0 
      wfcpc=0.0 
      wecpc=0.0 
      wvspc=0.0 
      wgridpc=0.0 
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
          if(xpsi(kk).lt.0.0.or.xpsi(kk).gt.1.0) cycle
          npc=npc+1 
          wsilpc=wsilpc+gsilpc(:,kk) 
          wmp2pc=wmp2pc+gmp2pc(:,kk) 
          wfcpc=wfcpc+rfcpc(:,kk) 
          wecpc=wecpc+gridec(kk,:) 
          wvspc=wvspc+gridvs(kk,:) 
          do ii=1,nw 
            do jj=1,nh 
              kkkk=(ii-1)*nh+jj 
              mj=iabs(j-jj)+1 
              mk=(i-1)*nh+mj 
              wgridpc(kkkk)=wgridpc(kkkk)+gridpc(mk,ii) 
              if(xpsi(kkkk).lt.0.0.or.xpsi(kkkk).gt.1.0) cycle
              wpcpc=wpcpc+gridpc(mk,ii) 
            enddo
          enddo
        enddo
      enddo
      xnpc=real(npc,dp) 
      wsilpc=wsilpc/xnpc 
      wmp2pc=wmp2pc/xnpc 
      wfcpc=wfcpc/xnpc 
      wecpc=wecpc/xnpc 
      wvspc=wvspc/xnpc 
      wgridpc=wgridpc/xnpc 
      wpcpc=wpcpc/xnpc**2 
! 
      open(unit=nffile,status='old',form='unformatted', & 
           file='rpcxx.dat',iostat=istat)
      if(istat.eq.0) close(unit=nffile,status='delete')
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
      else solution_mode
      if(ivesel.eq.2) vcurrt=vloopt(jtime)/rsisvs
      fixed_bdry_2: if (nbdry.gt.0) then
      Solovev: if (islve.gt.0) then
!------------------------------------------------------------------------------ 
!--     Solovev equilibrium                                                    -- 
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
        rotation: if (kvtor.le.0) then 
!--------------------------------------------------------------------- 
!--       no rotation                                               -- 
!--------------------------------------------------------------------- 
          scc1=sqrt(2./sbeta)/srm**2 
          if(ssrm.lt.0.0) saaa=xlmint/srm/sqrt(1.-2.*scc1) 
          if(ssrm.gt.0.0) saaa=xlmax/srm/sqrt(1.+2.*scc1) 
          srma=srm*saaa 
          dth=twopi/real(nbdry,dp) 
          do i=1,nbdry 
            th=(i-1)*dth 
            rbdry(i)=srma*sqrt(1.-2.*scc1*cos(th)) 
            zbdry(i)=sin(th) 
            zbdry(i)=saaa*zbdry(i)*seee 
          enddo
        else rotation
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
                      rnow,znow,negcur,kerror,1) 
          if(kerror.gt.0) return 
          xmin=xout(1) 
          xmax=xmin 
          do i=2,nfound 
            if(xout(i).lt.xmin) xmin=xout(i) 
            if(xout(i).gt.xmax) xmax=xout(i) 
          enddo 
          if(ssrm.lt.0.0) saaa=xlmint/xmin 
          if(ssrm.gt.0.0) saaa=xlmax/xmax 
          nskip=nfound/mbdry+1 
          j=0 
          do i=1,nfound,nskip 
            j=j+1 
            rbdry(j)=xout(i)*saaa 
            zbdry(j)=yout(i)*saaa 
          enddo 
          nbdry=j ! nbdry change inside of .gt.0 block
          srma=srm*saaa 
          rvtor=srma 
          rbetaw=sbetaw/sbeta 
        endif rotation
      endif Solovev
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
      SOL: if (nsol.gt.0) then 
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
      endif SOL
!----------------------------------------------------------------------- 
!--   set up parameters for fixed boundary calculations               -- 
!----------------------------------------------------------------------- 
      if(ifref.eq.-1) ifref=1 
      if (nbdry.gt.1) then ! nbdry changed above
        delx2=(rbdry(1)-rbdry(nbdry))**2 
        dely2=(zbdry(1)-zbdry(nbdry))**2 
        if((delx2+dely2).le.1.0e-08_dp) nbdry=nbdry-1 
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
      if(cfcoil.lt.0.) cfcoil=100./ipmeas(jtime)*abs(cfcoil) 
      if(cupdown.lt.0.) cupdown=100./ipmeas(jtime)*abs(cupdown) 
!----------------------------------------------------------------------- 
!--   symmetrize  F coil responses if needed                          -- 
!----------------------------------------------------------------------- 
      if ((symmetrize).and.(nbdry.gt.1)) then 
        do i=1,nw 
          do j=1,nh 
            kkl=(i-1)*nh+j 
            kku=i*nh-j+1 
            do m=nfsum/2+1,nfsum 
              gridfc(kkl,m)=gridfc(kku,m-nfsum/2) 
            enddo 
          enddo 
        enddo 
      endif 
!----------------------------------------------------------------------- 
!--   interpolate to get boundary response functions, first F coils   -- 
!----------------------------------------------------------------------- 
      do n=1,nfsum
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
      if ((symmetrize).and.(nbdry.gt.1)) then ! nbdry changed above
        do i=1,nbryup 
          if (ilower(i).ne.-1) then 
            do j=nfsum/2 +1, nfsum 
              jupper=j-nfsum/2 
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
      endif fixed_bdry_2
      endif solution_mode
! 
      return 
 5000 format (2e12.6) 
 5020 format (1x,4e16.9) 
 6550 format (/,' ir = ',10i4) 
 6555 format (' iz = ',10i4) 
 6557 format (/,' npc = ',i4) 
10200 format (6e12.6) 
10220 format (5e10.4) 
      end subroutine data_input
