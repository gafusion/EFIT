!**********************************************************************
!>
!!    sets default values for the efitin and in1 namelists
!!      this also initializes variables before a new read
!!
!**********************************************************************
      subroutine set_defaults()
      use efit_bdata, only: iunit,m_write
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
!--------------------------------------------------------------------- 
!--   All Modes (efitin and in1 namelists)
!--------------------------------------------------------------------- 
      alphafp=0.
      avemsels=10.
      bitmpi=0.0
      condin=1.0e-06_dp
      dpsiecn=0.0
      dtmsefull=0.0
      ecurbd=0.0
      elomin=0.90_dp
      errdelz=0.06_dp
      errmin=1.e-02_dp
      error=1.e-02_dp
      fcurbd=1.
      fe_psin=0.98_dp
      fe_width=0.02_dp
      fitdelz=.false.
      fitsiref=.false.
      fitzts='no'
      fmlscut=1.e-6_dp
      fwacoil=0.
      fwtbdry=1.0
      fwtec=0.
      fwtfc=0.
      fwtgam=0.0
      fwtmp2=0.
      fwtnow=0.001_dp ! input value is never used (hard coded)
      fwtqa=0.
      fwtsi=0.
      fwtsol=1.0
      f2edge=0.0
      gain=.22_dp
      gainp=0.75_dp
      iacoil=0
      iavem=5
      iaveus=0
      ibatch=0 ! never used in code (just gets output)
      ibtcomp=1
      icinit=2
      icprof=0
      idebug=0 ! deprecated
      ierchk=1
      ifcurr=0
      ifindopt=2
      ifitdelz=1
      ioffr=-7
      ioffz=7
!jal 04/23/2004
      iplcout=0
      irfila=0.0
      ishiftz=0
      itek=0
      itrace=1
      jdebug='NONE' ! deprecated
      jzfila=0.0
      kbound=0
      kcallece=2
      kcalpa=0
      kcgama=0
      kcmin=0
      kdomsels=0
      kedgef=0
      kedgep=0
      keecur=0
      keefnc=0
      kersil=3
      keqdsk=1
      kfffnc=0
      kfitece=0
      kfixrece=0
      kfixro=0
      kframe=0
      kinputece=0 ! never used in code
      klabel=0
      kppfnc=0
      kskipvs=0
      kwaitmse=0
      kwripre=0
      kwwfnc=0
      kzeroj=0
      limid=33
      lring=0
      msebkp=0
      msefitfun=1
      mse_certree=0
      mse_spave_on=0
      mse_strict=0
      mse_usecer=0
      mse_use_cer210=0
      mse_use_cer330=0
      mtxece=0
      mxiter=-25
      nbdrymx=110
      nbeam=0
      nbskip=2
      nccoil=1
      necein=0
      nextra=1
      nfit=0
      nharm=1
      nh_sub=nh
      nicoil=1
      nmass=0
      nslref=1
      nw_sub=nw
      n1coil=0
      ok_30rt=0
      ok_210lt=0
      oldccomp=.false.
      oldcomp=.false.
      pcurbd=1.
      pedge=0.0
      pe_psin=0.98_dp
      pe_width=0.02_dp
      psibit=0.0
      psiecn=0.0
      psiwant=1.0
      qvfit=0.95_dp
      rbdry=0.0
      rcentr=1.6955_dp
      relax=1.00
      relaxdz=1.0
      req_valid=.false.
      rexpan=0.010_dp
      rmbit=0.0
      robit=0.0
      rpbit=0.0
      rtep=0.0
      rteo=0.0
      zteo=0.0
      rtem=0.0
      rzeroj=0.0 
      saicon=80.0
      saimin=80.0
      scaledz=1.e-03_dp
      scrape=0.030_dp
      serror=0.03_dp
      sigrbd=1.e10_dp
      sigzbd=1.e10_dp
      siloplim=0.0
      sizeroj=0.0
      stabdz=-1.e-4_dp
      synmsels='SYN'
      teecein0=0.0
      tolbndpsi=1.0e-12_dp
      t_max_beam_off=0.0
      use_previous=.false.
      vbit=10.
      vsdamp=0
      vzeroj=0.0
      v30lt=0.0
      v30rt=0.0
      v210lt=0.0
      v210rt=0.0
      zbdry=0.0
      znose=-1.276_dp

      if (nesum.gt.0) ecurrt(1)=0.0
      if (nesum.gt.1) ecurrt(2)=0.0
      if (nesum.gt.2) ecurrt(3)=-1.e10_dp
      if (nesum.gt.3) ecurrt(4)=-1.e10_dp
      if (nesum.gt.4) ecurrt(5)=-1.e10_dp
      if (nesum.gt.5) ecurrt(6)=-1.e10_dp
!---------------------------------------------------------------------- 
!--   All SNAP Modes (efitin namelist)
!--------------------------------------------------------------------- 
      if ((kdata.ne.1).and.(kdata.ne.2)) then
        alternate_pointname_file='/link/efit/pcsnames.dat'
        do_spline_fit=.true.
        iaved=5
        iavev=10
        ircfact=0
        kcaldia=0
        lookfw=1
        nvtime=-1
        sizeroj(1)=-1.0
        use_alternate_pointnames=0
      endif
!---------------------------------------------------------------------- 
!--   K-File Modes (efitin namelist)
!--------------------------------------------------------------------- 
      if ((kdata.eq.5).and.(kdata.eq.6)) then
        appendsnap='KG'
        fwtcur=1.
        fwtbp=1.
        iout=1                 ! default - write fitout.dat
        kffcur=3
        kppcur=3
        limid=33
        tangam_uncor=0.0
!---------------------------------------------------------------------- 
!--   Full Solution Modes (efitin and in1 namelists)
!---------------------------------------------------------------------- 
      else
        aelip=0.60_dp
        backaverage=.false.
        betap0=0.50_dp
        betapw0=0.0
        bitip=0.0
        bzmse=0.0
        cfcoil=-1.
        chigamt=0.0
        chilibt=0.0
        chipre=0.0
        condno=0.0
        cstabne=0.0 ! hardcoded option
        cstabte=0.0 ! hardcoded option
        cupdown=-100000.
        curmid=0.0
        cutip=80000.
        dfsqe=0.0
        eceerror=0.03_dp
        eceiter='pair'
        eebdry=0.0
        ee2bdry=0.0
        eeknt=0.0
        eelip=1.2_dp
        eetens=5.0_dp
        emaxis=1.3_dp
        emf=1.00
        emp=1.00
        enf=1.00
        enp=1.00
        errbry=1.0e-04_dp
        errsil=0.03_dp
        fbetap=0.0
        fbetat=0.0
        fco2ne=1.0
        fcsum=1.0
        fczero=1.0
        ffbdry=0.0
        ff2bdry=0.0
        ffknt=0.0
        fftens=0.0
        fixpp=.false.
        fli=0.0
        fwtbp=0.0
        fwtcur=0.0
        fwtdlc=0.0
        fwtece0=0.0
        fwtecebz0=0.0
        fwtxli=1.
        fwtxx=0.2_dp
        fwtxxb=1.
        fwtxxj=1.
        fwtxxq=1.
        gammap=1.0e+10_dp
        ibound=0 ! hardcoded variable, not part of any namelist...
        ibtcomp=1
        ibunmn=1
        icalbet=1
        icntour=0
        icondn=-1
        icurrt=2
        idplace=0
        iecoil=0
        iexcal=0
        ifref=1
        isolve=1
        iplim=0
        iplots=1
        iprobe=0
        iqplot=1
        isetfb=0
        islve=0
        isumip=0
        itimeu=0
        iunit=35
        iweigh=0
        ixray=0
        ixstrt=1
        jbeta=1
        jli=2
        jwantm=3
        kccoils=0
        kcomega=0
        kdofit=0
        kdomse=0
        kdopre=0
        kdovt=0
        keebdry=0
        kee2bdry=0
        kffbdry=0
        kff2bdry=0
        kffcur=1
        kinput=0
        kplotp=1
        kppcur=3
        kppbdry=0
        kpp2bdry=0
        kprfit=0
        ksxr0=0
        ksxr2=0
        kvtor=0
        kwwbdry=0.0
        kww2bdry=0.0
        kwwcur=2
        limfag=2
        limloc='N/A'
        limitr=-33
        limvs=1
        m_write=1 ! TODO: output file type option should not be hard coded...
        nbdry=0
        nconstr=1
        ndokin=1
        npsi_ext=-1
        nqwant=0
        nsol=0
        nsplot=4 
        n1coil=0
        ppbdry=0.0
        pp2bdry=0.0
        ppknt=0.0
        pptens=0.0
        pressr=0.0
        pressw=0.0
        preswb=0.0
        psibry=0.0 
        qemp=0.0 
        qenp=0.95_dp 
        qpsi=0.0
        rbound=0.0
        rdjdz=0.0
        relip=1.68_dp
        rmaxvs=100.
        rminvs=0
        rmp2ac=0.0
        rpress=0.0
        rsisec=-1.
        rvtor=1.70_dp
        rzero=1.6955_dp
        chiecc=0.0
        scalea=.false.
        scalesir=1.0e-3_dp
        sidif=-1.0e+10_dp
        sigppb=1000.
        sigpre=0.0
        symmetrize=.false. 
        rmaxis=rzero
        vcurrt=0.0
        wwbdry=0.0
        ww2bdry=0.0
        wcurbd=0.0
        wwknt=0.0
        wwtens=0.0
        xlim=0.0
        ylim=0.0
        zbound=0.0
        zelip=0.0
        zmaxvs=100.
        zminvs=-100.

        ! DIII-D positions
        if(nco2v.gt.0) chordv(1)=1.486_dp
        if(nco2v.gt.1) chordv(2)=1.945_dp
        if(nco2v.gt.2) chordv(3)=2.098_dp
        if(nco2r.gt.0) chordr(1)=0._dp
        if(nco2r.gt.1) chordr(2)=0.1524_dp
        idosxr=1
      endif
!---------------------------------------------------------------------- 
!--   File Mode (in1 namelist)
!---------------------------------------------------------------------- 
      mode: if ((kdata.eq.1).or.(kdata.eq.2)) then
        akchiwt=1.0
        akgamwt=0.0
        akerrwt=0.0
        akprewt=0.0
        aktol=0.1_dp
        alpax(1)=1.e4_dp
        bitfc=0.0
        fwtpre=1.
        gamax(2)=-1.e6_dp 
        iout=1                 ! default - write fitout.dat
        kautoknt=0
        kakloop=1
        kakiter=25
        kbetapr=0
        kfcurb=0
        kplotpr=1
        kpcurb=0
        kpressb=0
        ncstne=1
        ncstte=1
        prbdry=0.
        salpha=1._dp/40._dp
        sbeta=1._dp/8._dp
        sbetaw=0.0
        scalepr(1)=-1.
        scalepw(1)=-1.
        sgprmin=0.0
        srm=-3.5_dp
        vcurfb(1)=0.0
        vcurfb(2)=500.
        vcurfb(3)=500.

        xlim(1)=-1.0
        rbdry(1)=-1.0
        nbdryp=-1
        ecefit=0.0
        ecebzfit=0.0

        geqdsk_ext='none'
        nw_ext=0
        nh_ext=0
        psin_ext(1)=-1000.0 
        sign_ext=1.0
        cratio_ext=1.0
        cratiop_ext=1.0
        cratiof_ext=1.0
        scalepp_ext=1.0
        scaleffp_ext=1.0

        cgama=0.0
        calpa=0.0
        bitec=0.0
!---------------------------------------------------------------------- 
!--   SNAP Modes (efitin namelist)
!--------------------------------------------------------------------- 
      elseif ((kdata.ne.5).and.(kdata.ne.6)) then mode
        ecurrt=0.0
        idite=0
        gammap=1./gammap
        gammaf=gammap
        snapfile='none'
        nsnapf=66
! -- Qilong Ren
        write_kfile=.false.
        fitfcsum=.false.
!----------------------------------------------------------------------
!--     istore=0 : Leave EFIT results in the run directory. Otherwise
!--       collect them in store_dir. 
!----------------------------------------------------------------------
        istore=0
      endif mode
      end subroutine set_defaults
