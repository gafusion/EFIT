      subroutine getsets(ktime,kwake,mtear,kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          getsets performs inputing and initialization.           **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**          24/07/85..........revised                               **
!**          04/03/93..........revised name for NCAR                 **
!**          23/04/04...JAL iplcout added to namelist used in weqdsk **
!**          01/08/07...DPB namelist for mag uncertainty added       **
!**                                                                  **
!**********************************************************************
      use set_kinds
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      include 'basiscomdu.inc'
! MPI >>>
#if defined(USEMPI)
      include 'mpif.h'
#endif
! MPI <<<
      logical lopened
      character filenm*15,ishotime*12,news*72, &
                eqdsk*20,comfile*15,prefix1*1,header*42,fit_type*3
      dimension coils(nsilop),expmp2(magpri), &
                denr(nco2r),denv(nco2v), &
                tgamma(nmtark),sgamma(nmtark),rrrgam(nmtark), &
                zzzgam(nmtark),aa1gam(nmtark),aa2gam(nmtark), &
                aa3gam(nmtark),aa4gam(nmtark),aa5gam(nmtark), &
                aa6gam(nmtark),aa7gam(nmtark), &
                tgammauncor(nmtark)
      real*8 bmsels(nmsels),sbmsels(nmsels),fwtbmsels(nmsels), &
         rrmsels(nmsels),zzmsels(nmsels),l1msels(nmsels),l2msels(nmsels), &
         l4msels(nmsels),emsels(nmsels),semsels(nmsels),fwtemsels(nmsels)
      dimension tlibim(libim),slibim(libim),rrrlib(libim)
      dimension devxmpin(magpri),rnavxmpin(magpri) &
               ,devpsiin(nsilop),rnavpsiin(nsilop) &
               ,devfcin(nfcoil),rnavfcin(nfcoil) &
               ,devein(nesum),rnavecin(nesum)
      character*82 snap_ext
!vasorg      character*82 snap_file
      namelist/in1/ishot,itime,itimeu,qvfit,plasma,expmp2,coils,btor, &
           fwtsi,fwtcur,limitr,fwtmp2,kffcur,kppcur,fwtqa,ierchk, &
           fwtbp,serror,nextra,scrape,itrace,itek,xltype,rcentr,bitip, &
           psibit,bitmpi,denr,denv,siref,fwtfc,brsp,bitfc,iecurr,iplim, &
           ecurrt,ifitvs,vloop,dflux,ifcurr,iavem,icprof,currn1,n1coil, &
           pnbeam,error,errmin,mxiter,xltype_180,icutfp,keqdsk,ibtcomp, &
           fcurbd,pcurbd,kbound,alphafp,kskipvs,vsdamp,kframe,zelip, &
           fwtdlc,sigdlc,elomin,kcalpa,kcgama,calpa,cgama,xalpa,xgama, &
           kzeroj,rzeroj,iaveus,relax,fwtec,bitec,fitsiref, &
           kppfnc,kppknt,ppknt,pptens,kfffnc,kffknt,ffknt,fftens, &
           kwwfnc,kwwknt,wwknt,wwtens,nbdry,rbdry,zbdry, &
           ppbdry,kppbdry,pp2bdry,kpp2bdry, &
           ffbdry,kffbdry,ff2bdry,kff2bdry, &
           wwbdry,kwwbdry,ww2bdry,kww2bdry &
           ,ktear,kersil,iout,ixray,table_dir,input_dir,store_dir &
           ,kpphord,kffhord,keehord,psiecn,dpsiecn,fitzts,isolve &
           ,iplcout,imagsigma,errmag,saimin,errmagb,fitfcsum,fwtfcsum &
           ,appendsnap,vbit,nbdrymx,efitversion,ifindopt,tolbndpsi
      namelist/inwant/psiwant,vzeroj,nccoil,currc79,currc139,rexpan, &
           znose,sizeroj,fitdelz,relaxdz,errdelz,oldccomp,nicoil, &
           oldcomp,currc199,curriu30,curriu90, &
           curriu150,curril30,curril90,curril150,ifitdelz,scaledz
      namelist/inms/xmprcg,xmp_k,vresxmp,t0xmp,psircg,psi_k,vrespsi &
           ,t0psi,fcrcg,fc_k,vresfc,t0fc,ercg,e_k,vrese,t0e,bcrcg &
           ,bc_k,vresbc,t0bc,prcg,p_k,vresp,t0p,bti322in,curc79in &
           ,curc139in,curc199in,devxmpin,rnavxmpin,devpsiin,rnavpsiin &
           ,devfcin,rnavfcin,devein,rnavecin
      namelist/ink/isetfb,ioffr,ioffz,ishiftz,gain,gainp,idplace &
           ,symmetrize,backaverage
      namelist/ins/tgamma,sgamma,fwtgam,rrrgam,zzzgam,aa1gam,aa2gam, &
                   aa3gam,aa4gam,aa5gam,aa6gam,aa7gam,msebkp, &
            msefitfun,mse_quiet,mse_spave_on,kwaitmse, &
            dtmsefull,mse_strict,t_max_beam_off,ok_30rt,ok_210lt,&
            mse_usecer,mse_certree,mse_use_cer330,mse_use_cer210,&
            tgammauncor,v30lt,v30rt,v210lt,v210rt
      namelist/in_msels/bmsels,sbmsels,fwtbmsels,rrmsels,zzmsels, &
            l1msels,l2msels,l4msels,emsels,semsels,fwtemsels,kdomsels, &
            fmlscut,synmsels,avemsels
      namelist/ina/spatial_avg_gam
      namelist/inece/necein,teecein0,feece0,errorece0,fwtece0,fwtecebz0 &
           ,ecefit,ecebzfit,kfitece,kinputece,kcallece,nharm &
           ,kfixro,rteo,zteo,kfixrece,rtep,rtem,rpbit,rmbit,robit &
           ,nfit,kcmin,fwtnow,mtxece
      namelist/iner/keecur,ecurbd,keefnc,eetens,keebdry,kee2bdry, &
                    eebdry,ee2bdry,eeknt,keeknt,keehord
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
           psiecn,dpsiecn,relaxdz,fitzts,isolve,stabdz &
           ,iplcout,errdelz,imagsigma,errmag,ksigma,saimin,errmagb &
           ,write_Kfile,fitfcsum,fwtfcsum,appendsnap &
           ,mse_quiet,mse_spave_on,kwaitmse,dtmsefull &
           ,mse_strict,t_max_beam_off,ifitdelz,scaledz &
           ,mse_usecer,mse_certree,mse_use_cer330,mse_use_cer210 &
           ,ok_30rt,ok_210lt,vbit,nbdrymx,fwtbmsels,fwtemsels,idebug,jdebug &
           ,synmsels,avemsels,kwritime,v30lt,v30rt,v210lt,v210rt,ifindopt,tolbndpsi
      namelist/efitink/isetfb,ioffr,ioffz,ishiftz,gain,gainp,idplace &
           ,symmetrize,backaverage,lring
      data mcontr/35/,lfile/36/,ifpsi/0/
      data currn1/0.0/,currc79/0.0/,currc139/0.0/,currc199/0.0/ &
                      ,curriu30/0.0/,curriu90/0.0/,curriu150/0.0/ &
                      ,curril30/0.0/,curril90/0.0/,curril150/0.0/
      logical exists
      integer, intent(inout) :: kerror

      kerror = 0
      table_di2 = table_dir
! --- find length of default directories
      ltbdir=0
      lindir=0
      lstdir=0
      do i=1,len(table_dir)
         if (table_dir(i:i).ne.' ') ltbdir=ltbdir+1
         if (input_dir(i:i).ne.' ') lindir=lindir+1
         if (store_dir(i:i).ne.' ') lstdir=lstdir+1
      enddo
      ltbdi2=ltbdir
!
      mdoskip=0
      iout=1                 ! default - write fitout.dat
      appendsnap='KG'
      snapextin='none'
      if (kwake.eq.1) go to 10
      patmp2(1)=-1.
      tmu0=twopi*tmu
      tmu02=tmu0*2.0
      errorm=1.
      ibatch=0
      ilaser=0
!----------------------------------------------------------------------
!-- news and help information                                        --
!----------------------------------------------------------------------
! MPI >>>
      ! ONLY root process displays EFIT news and help information
      if (rank == 0) then
        open(unit=80,status='old', &
             file=input_dir(1:lindir)//'efithelp.txt',err=83220)
        do i=1,100
          read (80,83210,end=83220,err=83220) news
          write (nttyo,83210) news
        enddo
        close(unit=80)
      endif
83210 format (a)
! MPI <<<
!
83220 continue
! OPT_INPUT >>>
! MPI >>>
      ! ONLY root process allowed to interface with terminal
      if (rank == 0) then
        if (use_opt_input .eqv. .false.) then
          write (nttyo,5500) (mfvers(i),i=1,2)
          write (nttyo,6000)
          read (ntty,*) kdata
        else
          kdata = mode_in
        endif
      endif
#if defined(USEMPI)
      if (nproc > 1) then
        call MPI_BCAST(kdata,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      endif
#endif
! MPI <<<
! OPT_INPUT <<<
      if (kdata.lt.0) then
        kdata=-kdata
        ilaser=1
      endif
!---------------------------------------------------------------------
!--  KDATA=16, wake-up mode driven by file WAKEMFIT.DAT consisting  --
!--  of shot # and time, -shot # for quit                           --
!---------------------------------------------------------------------
      if (kdata.eq.16) then
        kwake=1
        kdata=3
        jwake=kwake
        mdoskip=1
      endif
!----------------------------------------------------------------------
!--   Changed kdata <= 7 to kdata < 7            --
!--   Snap-Extension mode = 7              --
!----------------------------------------------------------------------
   10 if (kdata.ge.5.and.kdata.lt.7) go to 3000
      if (kdata.eq.8) go to 3000
      if (kwake.eq.1.and.mdoskip.eq.0.and.(iand(iout,1).ne.0)) close(unit=nout)
! MPI >>>
! ONLY root process can check for existence of fitout.dat file
      if (rank == 0) then
        ! Delete fitout.dat if already exists
        open(unit=nout,status='old',file='fitout.dat',err=12913)
        close(unit=nout,status='delete')
12913   continue
      endif
      if (iand(iout,1).ne.0) then
#if defined(USEMPI)
        if (nproc > 1) then
          if (rank == 0) then
            open(unit=nout,status='new',file='fitout.dat')
          endif
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          if (rank > 0) then
            open(unit=nout,status='old',file='fitout.dat')
          endif
        else
          open(unit=nout,status='new',file='fitout.dat')
        endif
#else
        open(unit=nout,status='new',file='fitout.dat')
#endif
      endif
! MPI <<<
      if (kwake.eq.1.and.mdoskip.eq.0) go to 10999
      if (kdata.eq.2) go to 200
!----------------------------------------------------------------------
!--  look up magnetic data directly                                  --
!----------------------------------------------------------------------
      do 50 i=1,nsilop
        psibit(i)=0.0
        fwtsi(i)=0.0
   50 continue
      do 60 i=1,magpri
        bitmpi(i)=0.0
        fwtmp2(i)=0.0
   60 continue
      do 70 i=1,nstark
        fwtgam(i)=0.0
   70 continue
      do 71 i=1,nnece
        fwtece0(i)=0.0
   71 continue
      fwtecebz0=0.0
      backaverage=.false.
      bitip=0.0
      betap0=0.50_dp
      brsp(1)=-1.e+20_dp
      cfcoil=-1.
      cutip=80000.
      do 185 i=1,nesum
        ecurrt(i)=0.0
        rsisec(i)=-1.
  185 continue
      emf=1.00
      emp=1.00
      enf=1.00
      enp=1.00
      error=1.0e-03_dp
      fbetap=0.0
      fbetat=0.0
      fcurbd=1.
      do 190 i=1,nfcoil
        fcsum(i)=1.0
        fczero(i)=1.0
        fwtfc(i)=0.
        rsisfc(i)=-1.
  190 continue
      do i=1,nesum
        fwtec(i)=0.0
      enddo
      do i=1,mbdry
       fwtbdry(i)=1.0
       fwtsol(i)=1.0
       sigrbd(i)=1.e10_dp
       sigzbd(i)=1.e10_dp
      enddo
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
      ifcurr=0
!jal 04/23/2004
      iplcout=0
      ifitvs=0
      ifref=-1
      itimeu=0
      iplim=0
      iprobe=0
      iqplot=1
      isetfb=0
      idplace=0
      ishot=-1
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
      kffcur=1
      kinput=0
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
      write_Kfile = .false.
      fitfcsum = .false.
      ifindopt = 2
      tolbndpsi = 1.0e-12_dp
!----------------------------------------------------------------------
!--   Snap-Extension mode              --
!--   Initialize istore = 0                                          --
!--   Central directory to collect EFIT results is the default       --
!--   directory. Otherwise in store_dir (default to /link/store/)    --
!----------------------------------------------------------------------
      istore = 0
!
   75 continue
      if (kdata.eq.3) then
      open(unit=neqdsk,status='old', &
           file='efit_snap.dat',err=80)
      snapfile='efit_snap.dat'
      go to 95
   80 continue
      open(unit=neqdsk,status='old',       &
           file= input_dir(1:lindir)//'efit_snap.dat'         )
      snapfile=input_dir(1:lindir)//'efit_snap.dat'
      endif
      if (kdata.eq.4) then
      open(unit=neqdsk,status='old', &
           file='efit_time.dat',err=85)
      go to 95
   85 continue
      open(unit=neqdsk,status='old', &
           file= input_dir(1:lindir)//'efit_time.dat'         )
      endif
!----------------------------------------------------------------------
!--    Snap-Extension mode                                           --
!----------------------------------------------------------------------
      if (kdata.eq.7) then
! OPT_INPUT >>>
! MPI >>>
! ONLY root process interaces with terminal
         if (rank == 0) then
           if (use_opt_input .eqv. .false.) then
             write (nttyo,6617)
             read (ntty,6620) snap_ext
           else
             snap_ext = snapext_in
           endif
         endif
         snapextin=snap_ext
#if defined(USEMPI)
         if (nproc > 1) then
           call MPI_BCAST(snap_ext,82,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
         endif
#endif
! MPI <<<
! OPT_INPUT <<<
   j = 0
   do i = 1, 82
      if (snap_ext(i:i) .ne. ' ') then
         j = j + 1
         snap_ext(j:j) = snap_ext(i:i)
      endif
   enddo
   snap_file = 'efit_snap.dat_'//snap_ext
         open(unit=neqdsk,status='old', &
           file= snap_file,err=81)
         snapfile=snap_file
         go to 95
 81    continue
         open(unit=neqdsk,status='old', &
           file= input_dir(1:lindir)//snap_file,err=83  )
          snapfile=input_dir(1:lindir)//snap_file
         go to 95
  83   continue
       snap_file = snap_ext
         open(unit=neqdsk,status='old', &
           file= snap_file       )
         snapfile=snap_file
      endif


95    continue
      read (neqdsk,efitin,end=108)
 108   continue
      read (neqdsk,efitink,err=96,end=109)
 109   continue
   96 close(unit=neqdsk)
!----------------------------------------------------------------------
!--   writes out the efitin namelist. Flag iout = 32.                --
!----------------------------------------------------------------------
      if (iand(iout,32).ne.0) then
         open(unit=nin,status='unknown',file='efit_snap.dat_out', &
              err=11231)
         write(nin,efitin)
11231    close(unit=nin)
      endif
!---------------------------------------------------------------------
!---- recalculate length of default directories in case any change  --
!---------------------------------------------------------------------
      ltbdir=0
      lindir=0
      lstdir=0
      do i=1,len(table_dir)
         if (table_dir(i:i).ne.' ') ltbdir=ltbdir+1
         if (input_dir(i:i).ne.' ') lindir=lindir+1
         if (store_dir(i:i).ne.' ') lstdir=lstdir+1
      enddo
      table_di2 = table_dir
      ltbdi2 = ltbdir
      iteks=itek
      mxiters=mxiter
      zelipss=zelip
      n1coils=n1coil
!
      ktime=mtime
      mtear=ktear
      qenp=qvfit
      if (itek.lt.0) then
        ixray=1
        itek=iabs(itek)
      endif
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
      if (imagsigma.gt.0) then
         do_spline_fit=.false.
         saimin=300.
      endif
!---------------------------------------------------------------------
!-- adjust fit parameters based on basis function selected          --
!---------------------------------------------------------------------
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
      if (kzeroj.eq.1.and.sizeroj(1).lt.0.0) sizeroj(1)=psiwant
!---------------------------------------------------------------------
!-- wakeup mode KDATA=16                                            --
!---------------------------------------------------------------------
10999 continue
      if (kwake.eq.1) then
28000   inquire(file='wakeefit.dat',opened=lopened)
        if (lopened) close(unit=neqdsk)
        open(unit=neqdsk, file='wakeefit.dat', status='old', err=28002)
        go to 28005
28002   call lib$wait(10.0)
        go to 28000
28005   iread=0
28010   read (neqdsk,*,end=28020,err=28000) ishot,timeb,dtime,ktime
        iread=iread+1
        if (iread.ge.ireadold+1) go to 28020
        go to 28010
28020   close(unit=neqdsk)
        if (iread.le.ireadold) then
          call lib$wait(10.0)
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
      if (kwake.eq.0) then
!
! OPT_INPUT >>>
! MPI >>>
! ONLY root process interfaces with terminal
        if (rank == 0) then
          if (use_opt_input .eqv. .false.) then
            write (nttyo,6040)
            read (ntty,*) ishot,timeb,dtime,ktime
          else
            ishot = shot_in
            timeb = starttime_in
            dtime = deltatime_in
            ktime = steps_in
          endif
        endif
#if defined(USEMPI)
        if (nproc > 1) then
          call MPI_BCAST(ishot,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(ktime,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(timeb,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(dtime,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        endif
#endif
! MPI <<<
! OPT_INPUT <<<
      endif
! -- Qilong Ren
      iishot = ishot
      ttimeb = timeb
      ddtime = dtime
      kktime = ktime
!----------------------------------------------------------------------
!--   Set proper Green's directory table_dir based on shot number    --
!----------------------------------------------------------------------
      if (ishot.ge.112000) then
        if (ishot.lt.156000) then
          table_di2 = table_di2(1:ltbdi2)//'112000/'
        elseif (ishot.lt.168191) then
          if (kdata.ne.2) then
            table_di2 = table_di2(1:ltbdi2)//'156014/'
          else
            if (efitversion <= 20140331) then
               table_di2 = table_di2(1:ltbdi2)//'112000/'
            else
               table_di2 = table_di2(1:ltbdi2)//'156014/'
            endif
          endif
        elseif (ishot.lt.181292) then
          table_di2 = table_di2(1:ltbdi2)//'168191/'
        elseif (ishot.ge.181292) then
          table_di2 = table_di2(1:ltbdi2)//'181292/'
        endif
        ltbdi2=ltbdi2+7
      endif
!-------------------------------------------------------------------------------
!--  Set bit noise for ishot > 152000                                         --
!-------------------------------------------------------------------------------
      if (ishot.gt.152000) vbit = 80
!-------------------------------------------------------------------------------
!-- read in limiter data                                                      --
!-------------------------------------------------------------------------------
      call getlim(1,xltype,xltype_180)
!
  100 continue
      if (lookfw.ge.0) then
        do 102 i=1,magpri
        rwtmp2(i)=0.0
  102   continue
        do i=1,nsilop
         rwtsi(i)=0.0
        enddo
        open(unit=neqdsk,status='old', &
             file=table_di2(1:ltbdi2)//'fitweight.dat'         )
  105   read (neqdsk,*,end=107) irshot
        if (irshot.gt.ishot) go to 107
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
  107   continue
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
! MPI >>>
#if defined(USEMPI)
! MK 2020.10.08 TODO All control paths *should* be identical here
! MK i.e. only call getpts_mpi, regardless of nproc
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
      do 142 i=1,nstark
        swtgam(i)=fwtgam(i)
        if (fwtgam(i).gt.1.e-06_dp) mmstark=mmstark+1
  142 continue
      if (mmstark.gt.0) then
! MPI >>>
#if defined(USEMPI)
        if (nproc == 1) then
          call getstark(ktime)
        else
          call getstark_mpi(ktime)
        endif
#else
        call getstark(ktime)
#endif
! MPI <<<
      endif
!
      do jtime=1,ktime
        do 99142 i=1,nmsels
          swtbmselt(jtime,i)=fwtbmsels(i)
          swtemselt(jtime,i)=fwtemsels(i)
99142   continue
      enddo
      mmbmsels=0
      mmemsels=0
      do 99144 i=1,nmsels
        swtbmsels(i)=fwtbmsels(i)
        if (fwtbmsels(i).gt.1.e-06_dp) mmbmsels=mmbmsels+1
        if (fwtemsels(i).gt.1.e-06_dp) mmemsels=mmemsels+1
99144 continue
      if (mmbmsels.gt.0) then
! MPI >>>
#if defined(USEMPI)
        if (nproc == 1) then
          call getmsels(ktime)
        else
          call getmsels(ktime)
        endif
#else
        call getmsels(ktime)
#endif
! MPI <<<
      endif
!
      do 145 i=1,ktime
        time(i)=time(i)*1000.
  145 continue
!-----------------------------------------------------------------------
!-- Get edge pedestal tanh paramters                                  --
!-----------------------------------------------------------------------
      if (fitzts.eq.'te') then
        call gettanh(ishot,fitzts,ktime,time,ztssym,ztswid, &
                          ptssym,ztserr)
      endif
!----------------------------------------------------------------------
!-- save fitting weights for SNAP modes                              --
!----------------------------------------------------------------------
      swtdlc=fwtdlc
      swtcur=fwtcur
      do i=1,nfcoil
        swtfc(i)=fwtfc(i)
      enddo
      do i=1,nesum
        swtec(i)=fwtec(i)
      enddo
      do i=1,magpri
        if (lookfw.gt.0) then
           if (fwtmp2(i).gt.0.0) fwtmp2(i)=rwtmp2(i)
        endif
        swtmp2(i)=fwtmp2(i)
      enddo
      do i=1,nsilop
        if (lookfw.gt.0) then
           if (fwtsi(i).gt.0.0) fwtsi(i)=rwtsi(i)
        endif
        swtsi(i)=fwtsi(i)
      enddo
      go to 1000
!
  200 continue
! OPT_INPUT >>>
      ifname(:) = ''
! ONLY root process interfaces with terminal
      if (rank == 0) then
        if (use_opt_input .eqv. .false.) then
          write (nttyo,6200)
          read (ntty,*) ktime
          write (nttyo,6220)
          do 210 i=1,ktime
            write (nttyo,6230)
            read (ntty,6240) ifname(i)
  210     continue
        else
          ktime = steps_in
          do i=1,ktime
             ifname(i) = inpfile_in(i)
          enddo
        endif
      endif
#if defined(USEMPI)
! Distribute steps among ALL processes if necessary
      if (nproc > 1) then
        dist_data(:) = 0
        dist_data_displs(:) = 0
        if (rank == 0) then
! Compute number of steps per process
          i = 1
          do while (i <= ktime)
            do j=1,nproc
              if (i <= ktime) then
                dist_data(j) = dist_data(j)+1
                i = i+1
              endif
            enddo
          enddo
! Compute array displacements
          do i=2,nproc
            do j=1,i-1
! Input filenames are up to 80 characters and displacements given as number of bytes
              dist_data_displs(i) = dist_data_displs(i)+dist_data(j)*80
            enddo
          enddo
        endif
! Explicitly synchronize processes
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
! Distribute time step and filename information to ALL processes
        call MPI_SCATTER(dist_data,1,MPI_INTEGER,ktime,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! Recall each filename 80 characters
        if (rank == 0) then
          dist_data(:) = dist_data(:)*80
          call MPI_SCATTERV(ifname,dist_data,dist_data_displs,MPI_CHARACTER,MPI_IN_PLACE,dist_data(rank+1),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
        else
          call MPI_SCATTERV(ifname,dist_data,dist_data_displs,MPI_CHARACTER,ifname,ktime*80,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
        endif
      endif
#endif
! OPT_INPUT <<<
 1000 continue
      open(unit=mcontr,status='old',form='unformatted', &
           file=table_dir(1:ltbdir)//'ec'//trim(ch1)//trim(ch2)//'.ddd')
      read (mcontr) mw,mh
      if(.not.allocated(rgrid)) then
       allocate(rgrid(mw),stat=iallocate_stat)
      if(iallocate_stat/=0) stop "*** Not enough space for rgrid ***"
      endif
      if(.not.allocated(zgrid)) then
       allocate(zgrid(mh),stat=iallocate_stat)
      if(iallocate_stat/=0) stop "*** Not enough space for zgrid ***"
      endif
      if(.not.allocated(gridfc)) then
       allocate(gridfc(mw*mh,nfcoil),stat=iallocate_stat)
      if(iallocate_stat/=0) stop "*** Not enough space for gridfc ***"
      endif
!vas-------
      read (mcontr) rgrid,zgrid
      read (mcontr) gridfc
      if(.not.allocated(gridpc)) then
       allocate(gridpc(mw*mh,mw),stat=iallocate_stat)
      if(iallocate_stat/=0) stop "*** Not enough space for gridpc ***"
      endif
      read (mcontr) gridpc
      close(unit=mcontr)
!----------------------------------------------------------------------
!-- read in the f coil response functions                            --
!----------------------------------------------------------------------
      if(.not.allocated(rsilfc)) then
       allocate(rsilfc(nsilop,nfcoil),stat=iallocate_stat)
      if(iallocate_stat/=0) stop "*** Not enough space for rsilfc ***"
      endif
      if(.not.allocated(rmp2fc)) then
       allocate(rmp2fc(magpri,nfcoil),stat=iallocate_stat)
      if(iallocate_stat/=0) stop "*** Not enough space for rmp2fc ***"
      endif
      if(.not.allocated(gsilpc)) then
       allocate(gsilpc(nsilop,mw*mh),stat=iallocate_stat)
      if(iallocate_stat/=0) stop "*** Not enough space for gsilpc ***"
      endif
      if(.not.allocated(gmp2pc)) then
       allocate(gmp2pc(magpri,mw*mh),stat=iallocate_stat)
      if(iallocate_stat/=0) stop "*** Not enough space for gmp2pc ***"
      endif
!vas-------
      open(unit=nrspfc,form='unformatted', &
           status='old',file=table_dir(1:ltbdir)//'rfcoil.ddd')
      read (nrspfc) rsilfc
      read (nrspfc) rmp2fc
      close(unit=nrspfc)
!----------------------------------------------------------------------
!-- read in the plasma response function                             --
!----------------------------------------------------------------------
      open(unit=nrsppc,status='old',form='unformatted', &
           file=table_dir(1:ltbdir)//'ep'//trim(ch1)//trim(ch2)//'.ddd')
      read (nrsppc) gsilpc
      read (nrsppc) gmp2pc
      close(unit=nrsppc)
!
      if (kdata.ne.2) &
      call zlim(zero,nw,nh,limitr,xlim,ylim,rgrid,zgrid,limfag)
      drgrid=rgrid(2)-rgrid(1)
      dzgrid=zgrid(2)-zgrid(1)
      darea=drgrid*dzgrid
      tmu2=-pi*tmu*dzgrid/drgrid
!
      return
 3000 continue
      call write_K(ksstime,kerror)
      ktime = ksstime
      if (kerror.gt.0) return
      return

 4042 format (1x,a42,1x,a3)
 4958 format ('#!/bin/csh -f')
 4960 format ('      runefit.sc k',a12)
 4962 format ('#',/,'exit')
 4970 format (2x,a,1x,a,1x,a,1x,a)
 4980 format (i5)
 5000 format (2e12.6)
 5500 format (/,10x,'EFITD Version  ',2a5,/)
 6000 format (/,1x,'type mode (2=file, 3=snap, 4=time', &
               ', 5=input, 6=com file, 7=snap_ext,',    &
               ' 8=pefit):')
 6040 format (/,1x,'type shot #, start time(ms), time step(ms), steps' &
        ,'(<1001):')
 6080 format (/,1x,'type limiter position (cm, 0=ptdata):')
 6090 format(' enter number of extra field lines to trace:')
 6091 format(' enter scrape off depth(m),'/ &
       '       sense of tracing (+1 for down, -1 for up),'/ &
       '       ixstrt (+1 for start on outside, -1' &
       ' for start inside):')
 6100 format(/,1x,'type plot mode (0=none, 1=tektronix, 2=versatec, 3=qms, -=x ray):')
 6200 format (/,1x,'number of time slices?')
 6220 format (/,1x,'type input file names:')
 6230 format (1x,'#')
 6240 format (a)
 6600 format (/,1x,'good shot list file name ( 0=tty) ?')
 6610 format (/,1x,'command file name ( 0=none) ?')
 6617 format (/,1x,'type snap file extension (def for default):')
 6620 format (a)
 6700 format (a1,a12)
30000 format (i9)
30200 format (10f3.0)
      end
