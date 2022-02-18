#include "config.f"
!**********************************************************************
!>
!!    Subroutine for writing K file in K file mode Qilong Ren         --
!!    
!!    
!!
!!    @param ksstime :
!!
!!    @param kerror : Error flag
!!
!**********************************************************************
      subroutine write_K(ksstime,kerror)
      use set_kinds
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
#if defined(USEMPI)
      include 'mpif.h'
#endif

      character filenm*15,ishotime*12,news*72,table_s*72, &
                eqdsk*20,comfile*15,prefix1*1,header*42,fit_type*3
      integer*4 ltbdis,versiondate
      real*8,dimension(:),allocatable :: coils,expmp2, &
                denr,denv, tgamma,sgamma,rrrgam, &
                zzzgam,aa1gam,aa2gam, aa3gam,aa4gam,aa5gam, &
                aa6gam,aa7gam,tgammauncor
      real*8,dimension(:),allocatable :: bmsels,sbmsels,fwtbmsels, &
                rrmsels,zzmsels,l1msels,l2msels, &
                l4msels,emsels,semsels,fwtemsels
      real*8,dimension(:),allocatable :: tlibim,slibim,rrrlib
      character*82 snap_ext
      namelist/machinein/nsilds,nsilol,nfcoil,nrogow,nacoil,mfcoil,necoil,nvesel, &
      mpress,nesum,magpri67,magpri322,magprirdp,magudom,maglds,mse315,mse45, &
      mse15,mse1h,mse315_2,mse210,libim,nmsels,nnece,nnecein,neceo,nnnte, &
      ngam_vars,ngam_u,ngam_w,nlimit,nlimbd,nangle,ntangle,nfbcoil,mccoil, &
      micoil,ndata,nwwcur,nffcur,nppcur,nercur,ntime,ndim,kxiter,mqwant, &
      mbdry,mbdry1,nxtram,nxtlim,nco2v,nco2r,modef,modep,modew,kubics, &
      icycred_loopmax,nfourier
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
           kwwfnc,kwwknt,wwknt,wwtens,nbdry,rbdry,zbdry,vbit, nbdrymx, &
           ppbdry,kppbdry,pp2bdry,kpp2bdry, &
           ffbdry,kffbdry,ff2bdry,kff2bdry, &
           wwbdry,kwwbdry,ww2bdry,kww2bdry &
           ,ktear,kersil,iout,ixray,table_dir,input_dir,store_dir &
           ,kpphord,kffhord,keehord,psiecn,dpsiecn,fitzts,isolve &
           ,iplcout,imagsigma,errmag,saimin,errmagb,fitfcsum,fwtfcsum,efitversion &
           ,kwripre,ifindopt,tolbndpsi
      namelist/inwant/psiwant,vzeroj,nccoil,currc79,currc139,rexpan, &
           znose,sizeroj,fitdelz,relaxdz,errdelz,oldccomp,nicoil, &
           oldcomp,currc199,curriu30,curriu90, &
           curriu150,curril30,curril90,curril150,ifitdelz,scaledz
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
      integer*4, intent(out) :: ksstime
      integer*4, intent(inout) :: kerror

      ALLOCATE(coils(nsilop),expmp2(magpri), &
                denr(nco2r),denv(nco2v), &
                tgamma(nmtark),sgamma(nmtark),rrrgam(nmtark), &
                zzzgam(nmtark),aa1gam(nmtark),aa2gam(nmtark), &
                aa3gam(nmtark),aa4gam(nmtark),aa5gam(nmtark), &
                aa6gam(nmtark),aa7gam(nmtark),tgammauncor(nmtark))
      ALLOCATE(bmsels(nmsels),sbmsels(nmsels),fwtbmsels(nmsels), &
         rrmsels(nmsels),zzmsels(nmsels),l1msels(nmsels),l2msels(nmsels), &
         l4msels(nmsels),emsels(nmsels),semsels(nmsels),fwtemsels(nmsels))
      ALLOCATE(tlibim(libim),slibim(libim),rrrlib(libim))

      kerror = 0
!---------------------------------------------------------------------
!--  generate input files and command file for running EFIT-AI      --
!---------------------------------------------------------------------
      if (kdata .eq. 5 .or. kdata .eq. 6) then
        if (rank == 0) then
          if (use_opt_input .eqv. .false.) then
            write (nttyo,6610)
            read (ntty,6620) comfile
            write (nttyo,6600)
            read (ntty,6620) filenm
            if (filenm.eq.'0' .or. filenm.eq.'') then
              write (nttyo,6040)
              read (ntty,*) lshot,timeb,dtime,ktime
            endif
          else
            comfile = cmdfile_in
            filenm = shotfile_in
            if (filenm.eq.'0' .or. filenm.eq.'') then
              lshot = shot_in
              timeb = starttime_in
              dtime = deltatime_in
              ktime = steps_in
            endif
          endif
        endif
#if defined(USEMPI)
        if (nproc > 1) then
          call MPI_BCAST(comfile,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(filenm,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
          if (filenm.eq.'0' .or. filenm.eq.'') then
            call MPI_BCAST(lshot,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(ktime,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(timeb,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(dtime,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
          endif
        endif
#endif
        if (comfile.ne.'0' .and. comfile.ne.'') then
          open(unit=nffile,status='old',file=comfile,iostat=ioerr)
          if (ioerr.ne.0) close(unit=nffile,status='delete')
          open(unit=nffile,status='new',file=comfile)
          write (nffile,4958)
        endif
        if (filenm.ne.'0' .and. filenm.ne.'') then
          open(unit=nout,status='old',file=filenm,iostat=ioerr)
          if (ioerr.ne.0) open(unit=nout,status='old', &
            file='phys_data:[d3phys.diiid.gsl]'//filenm,iostat=ioerr)
          if (ioerr.ne.0) then
            call errctrl_msg('write_K', &
                             'could not open good shot list file')
            stop
          endif
        endif
      endif
      open(unit=neqdsk,status='old',file='efit_snap.dat',iostat=ioerr)
      if (ioerr.eq.0) then
        snapfile='efit_snap.dat'
      else
        open(unit=neqdsk,status='old', &
           file= input_dir(1:lindir)//'efit_snap.dat')
        snapfile=input_dir(1:lindir)//'efit_snap.dat'
      endif
      do i=1,nfcoil
        fwtfc(i)=0.
      enddo
      do i=1,nesum
        fwtec(i)=0.
      enddo
      do i=1,mbdry
       fwtbdry(i)=1.0
       fwtsol(i)=1.0
       sigrbd(i)=1.e10_dp
       sigzbd(i)=1.e10_dp
      enddo
      do i=1,magpri
        fwtmp2(i)=0.
        bitmpi(i)=0.0
      enddo
      do i=1,nstark
        fwtgam(i)=0.0
      enddo
      do i=1,nsilop
        fwtsi(i)=0.
        psibit(i)=0.0
      enddo
      error=1.e-03_dp
      fcurbd=1.
      fwtcur=1.
      fwtbp=1.
      fwtqa=0.
      iaved=5
      iavem=5
      iavev=10
      icprof=0
      icutfp=0
      iecurr=1
      ierchk=1
      ifcurr=0
      ifitvs=0
      itek=0
      iout=1                 ! default - write fitout.dat
      itrace=1
!jal 04/23/2004
      iplcout=0
      keqdsk=1
      kffcur=3
      kppcur=3
      limid=33
      lookfw=1
      mxiter=25
      nextra=1
      pcurbd=1.
      scrape=0.030_dp
      qvfit=0.95_dp
      serror=0.03_dp
      xltype=0.
      xltype_180=0.0
      ifindopt = 2
      tolbndpsi = 1.0e-12_dp
!
      read (neqdsk,efitin)
      read (neqdsk,efitink,iostat=ioerr)
      if (ioerr.ne.0) close(unit=neqdsk)
      if (imagsigma.eq.1) then
        do_spline_fit=.false.
        saimin=150.
      endif
#ifdef DEBUG_LEVEL2
      write (6,*) 'WRITE_K fwtbmsels= ',(fwtbmsels(i),i=1,nmsels)
#endif
!----------------------------------------------------------------------
!--   recalculate length of default directories in case any change   --
!----------------------------------------------------------------------
      ltbdir=0
      lindir=0
      lstdir=0
      do i=1,len(table_dir)
        if (table_dir(i:i).ne.' ') ltbdir=ltbdir+1
        if (input_dir(i:i).ne.' ') lindir=lindir+1
        if (store_dir(i:i).ne.' ') lstdir=lstdir+1
      enddo
      ! TODO: jtime is undefined here...
!      if (lshot.ge.112000.and.jtime.le.1) then
      if (lshot.ge.112000) then
        if (lshot.lt.156000) then
          table_di2 = table_dir(1:ltbdir)//'112000/'
        elseif  (ishot.lt.168191) then
          if ((kdata.ne.1).and.(kdata.ne.2)) then
            table_di2 = table_dir(1:ltbdir)//'156014/'
          else
            read(efitversion,*,iostat=ioerr) versiondate
            if (ioerr.ne.0) versiondate = 99999999 ! EFIT-AI hash string
            if (versiondate <= 20140331) then
              table_di2 = table_dir(1:ltbdir)//'112000/'
            else
              table_di2 = table_dir(1:ltbdir)//'156014/'
            endif
          endif
        elseif (ishot.lt.181292) then
          table_di2 = table_dir(1:ltbdir)//'168191/'
        elseif (ishot.ge.181292) then
          table_di2 = table_dir(1:ltbdir)//'181292/'
        endif
        ltbdi2=ltbdir+7
      endif
      efitversion = efitvers//" "
      ksstime = ktime
!---------------------------------------------------------------------
!--    specific choice of current profile                           --
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
      if (mse_usecer .eq. 1) keecur = 0
      if (mse_usecer .eq. 2 .and. keecur .eq. 0) then
        keecur = 2
        keefnc = 0
        itek = 5
      endif
!
      limitr=-limid
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
      swtcur=fwtcur
      do i=1,nstark
        swtgam(i)=fwtgam(i)
      enddo
      do i=1,nmsels
        swtbmsels(i)=fwtbmsels(i)
        swtemsels(i)=fwtemsels(i)
      enddo
 3046 continue
#ifdef DEBUG_LEVEL2
      write (6,*) 'WRITE_K fwtbmsels= ',(fwtbmsels(i),i=1,nmsels)
      write (6,*) 'WRITE_K swtbmsels= ',(swtbmsels(i),i=1,nmsels)
#endif
!
      if (filenm.ne.'0' .and. filenm.ne.'') then
        read (nout,4970,iostat=ioerr) ishotime
        if (ioerr.eq.0) then
          read (ishotime,fmt='(i6,1x,i5)',err=3046) ishot,itime
          times=itime/1000.
          delt=0.002_dp
          ktime=1
          timeb=itime
          dtime=0.
        endif
      else
        times=timeb/1000.
        delt=dtime/1000.
        ishot=lshot
        ioerr=0
      endif
      if (ioerr.eq.0) then
      if (ifitvs.gt.0) then
        istop=-1
      else
        istop=0
      endif
      if (ishot.gt.108281) n1coil = 0
!----------------------------------------------------------------------
!--   Fetch data                                                     --
!----------------------------------------------------------------------
#if defined(USEMPI)
      if (nproc == 1) then
        call getpts(ishot,times,delt,ktime,istop)
      else
        call getpts_mpi(ishot,times,delt,ktime,istop)
      endif
#else
      call getpts(ishot,times,delt,ktime,istop)
#endif
      if (istop.gt.0) then
        write (6,20000)
        if (filenm.ne.'0' .and. filenm.ne.'') then
          go to 3046
        else
          kerror = 1
          call errctrl_msg('write_K', 'shot data not on disk')
          return
        endif
      endif
!
      mmstark=0
      do i=1,nstark
        if (swtgam(i).gt.1.e-06_dp) mmstark=mmstark+1
      enddo
      if (mmstark.gt.0) then
#if defined(USEMPI)
        if (nproc == 1) then
          call getstark(ktime)
        else
          call getstark_mpi(ktime)
        endif
#else
        call getstark(ktime)
#endif
      endif
!
      mmbmsels=0
      mmemsels=0
      do i=1,nmsels
        if (swtbmsels(i).gt.1.e-06_dp) mmbmsels=mmbmsels+1
        if (swtemsels(i).gt.1.e-06_dp) mmemsels=mmemsels+1
      enddo
#ifdef DEBUG_LEVEL2
      write (6,*) 'WRITE_K mmbmsels= ', mmbmsels
#endif
      if (mmbmsels.gt.0) then
#if defined(USEMPI)
        if (nproc == 1) then
          call getmsels(ktime)
        else
          call getmsels(ktime)
        endif
#else
        call getmsels(ktime)
#endif
      endif
!----------------------------------------------------------------------
!--   get xltype, xltype_180 without reading limiters                --
!----------------------------------------------------------------------
      call getlim(0,xltype,xltype_180)
!----------------------------------------------------------------------
!--   adjust fitting weights based on FITWEIGHT.DAT                  --
!----------------------------------------------------------------------
      if (lookfw.ge.0) then
        do i=1,magpri
          rwtmp2(i)=0.0
        enddo
        do i=1,nsilop
         rwtsi(i)=0.0
        enddo
        open(unit=neqdsk,status='old', &
             file=table_di2(1:ltbdi2)//'fitweight.dat')
 3050   read (neqdsk,*,iostat=ioerr) irshot
        if ((ioerr.ne.0).and.(irshot.le.ishot)) then
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
          go to 3050
        endif
        close(unit=neqdsk)
      endif
!
      do i=1,ktime
        time(i)=time(i)*1000.
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
!-----------------------------------------------------------------------
!--   Get edge pedestal tanh paramters                                --
!-----------------------------------------------------------------------
      if (fitzts.eq.'te') then
        call gettanh(ishot,fitzts,ktime,time,ztssym,ztswid, &
                     ptssym,ztserr)
      endif
!
      do jtime=1,ktime
        itime=time(jtime)
        timems=itime
        timeus=(time(jtime)-timems)*1000.
        itimeu=timeus
!-----------------------------------------------------------------------
!--     correction for truncation                                     --
!-----------------------------------------------------------------------
        if (itimeu.ge.990) then
          itime=itime+1
          itimeu=0
          time(jtime)=itime
        endif
        do i=1,nsilop
          coils(i)=silopt(jtime,i)
          fwtsi(i)=swtsi(i)
          if (lookfw.gt.0.and.fwtsi(i).gt.0.0) fwtsi(i)=rwtsi(i)
          if (ierpsi(i).ne.0) fwtsi(i)=0.0
        enddo
        do i=1,magpri
          expmp2(i)=expmpi(jtime,i)
          fwtmp2(i)=swtmp2(i)
          if (lookfw.gt.0.and.fwtmp2(i).gt.0.0) fwtmp2(i)=rwtmp2(i)
          if (iermpi(i).ne.0) fwtmp2(i)=0.0
        enddo
        do i=1,nfcoil
          brsp(i)=fccurt(jtime,i)
          fwtfc(i)=swtfc(i)
          if (ierfc(i).ne.0) fwtfc(i)=0.0
        enddo
        do i=1,nesum
          ecurrt(i)=eccurt(jtime,i)
          fwtec(i)=swtec(i)
          if (ierec(i).ne.0) fwtec(i)=0.0
        enddo
        do i=1,nco2r
          denr(i)=denrt(jtime,i)
        enddo
        do i=1,nco2v
          denv(i)=denvt(jtime,i)
        enddo
        if (mmstark.gt.0) then
          do i=1,nmtark
            tgamma(i)=tangam(jtime,i)
            tgammauncor(i)=tangam_uncor(jtime,i)
            sgamma(i)=siggam(jtime,i)
            rrrgam(i)=rrgam(jtime,i)
            zzzgam(i)=zzgam(jtime,i)
            aa1gam(i)=a1gam(jtime,i)
            aa2gam(i)=a2gam(jtime,i)
            aa3gam(i)=a3gam(jtime,i)
            aa4gam(i)=a4gam(jtime,i)
            aa5gam(i)=a5gam(jtime,i)
            aa6gam(i)=a6gam(jtime,i)
            aa7gam(i)=a7gam(jtime,i)
            fwtgam(i)=swtgam(i)
            if (iergam(i).ne.0) fwtgam(i)=0.0
          enddo
        endif
!
        if (mmbmsels.gt.0) then
          do i=1,nmsels
            bmsels(i)=bmselt(jtime,i)
            sbmsels(i)=sbmselt(jtime,i)
            fwtbmsels(i)=swtbmsels(i)
            rrmsels(i)=rrmselt(jtime,i)
            zzmsels(i)=zzmselt(jtime,i)
            l1msels(i)=l1mselt(jtime,i)
            l2msels(i)=l2mselt(jtime,i)
            l4msels(i)=l4mselt(jtime,i)
            emsels(i)=emselt(jtime,i)
            semsels(i)=semselt(jtime,i)
            fwtemsels(i)=swtemsels(i)
            if (iermselt(jtime,i).ne.0) then
              fwtbmsels(i)= 0.0
              fwtemsels(i)= 0.0
            endif
          enddo
        endif
!
        fwtcur=swtcur
        if (ierpla.ne.0) fwtcur=0.0
        btor=bcentr(jtime)
        plasma=pasmat(jtime)
        siref=psiref(jtime)
        vloop=vloopt(jtime)
        dflux=1.0e+03_dp*diamag(jtime)
        sigdlc=1.0e+03_dp*sigdia(jtime)
        pnbeam=pbinj(jtime)
        if (n1coil.gt.0) currn1=curtn1(jtime)
        if (nccoil.gt.0) then
          currc79=curc79(jtime)
          currc139=curc139(jtime)
          currc199=curc199(jtime)
        endif
        curriu30=curiu30(jtime)
        curriu90=curiu90(jtime)
        curriu150=curiu150(jtime)
        curril30=curil30(jtime)
        curril90=curil90(jtime)
        curril150=curil150(jtime)
!-----------------------------------------------------------------------
!--     Set edge pedestal tanh paramters                              --
!-----------------------------------------------------------------------
        if (fitzts.eq.'te'.and.ztserr(jtime)) then
          nbdry=1
          rbdry(1)=1.94_dp
          zbdry(1)=ztssym(jtime)+0.5_dp*ztswid(jtime)
        endif
        call getfnmu(itimeu,'k',ishot,itime,eqdsk)
        open(unit=neqdsk,status='old',file=eqdsk,iostat=ioerr, &
             recl=72,delim='APOSTROPHE')
        if (ioerr.eq.0) close(unit=neqdsk,status='delete')
!-----------------------------------------------------------------------
!--     Set bit noise for ishot > 152000                              --
!-----------------------------------------------------------------------
!
!05-03-2016 Bob Johnson - problem seen on iris but not sun
!The next line gets rid of valgrind complaints on the
!following if statement checking ishot.  Also, ktear was printing
!as garbage when ktear is not defined in the input file.
!This statement fixes the problem - so there could still be
!an underlying bug or it could be a fortran issue with namelists.
!Note that ktear is part of 2 namelists.
!If ktear IS defined in the input file, everything is fine
!without this statement.
        if (ktear.gt.100) write(6,*) 'weird value of ktear=',ktear
        if (ishot.gt.152000) vbit = 80
        open(unit=neqdsk,file=eqdsk,status='new', &
             recl=72,delim='quote')
!             recl=72,delim='APOSTROPHE')
        write (neqdsk,machinein)
        write (neqdsk,in1)
        write (neqdsk,inwant)
        if (isetfb.ne.0) write (neqdsk,ink)
        if (mmstark.gt.0) write (neqdsk,ins)
        if (mmbmsels.gt.0) write (neqdsk,in_msels)
        if (kwaitmse.ne.0) write (neqdsk,ina)
        if (kfitece.gt.0) write (neqdsk,inece)
        if (keecur.gt.0) write (neqdsk,iner)
!-----------------------------------------------------------------------
!--     fitting type flag                                             --
!-----------------------------------------------------------------------
        if ((kprfit.gt.0).and.(mmstark.gt.0)) then
          fit_type = 'KIM'
        elseif (kprfit.gt.0) then
          fit_type = 'KIN'
        elseif (mmstark.gt.0) then
          fit_type = 'MSE'
          if (mmbmsels.gt.0) fit_type = 'MSL'
        elseif (mmbmsels.gt.0) then
          fit_type = 'MLS'
        else
          fit_type = 'MAG'
        endif
!
        call db_header(ishot,itime,header)
        write (neqdsk,4042) header,fit_type
!---------------------------------------------------------------------
!--     Append SNAP file                                            --
!---------------------------------------------------------------------
        if (appendsnap.eq.'K'.or.appendsnap.eq.'KG') then
          if (snapfile/='none') then
            open(unit=nsnapf,status='old', &
                 file=snapfile,iostat=ioerr)
            if (ioerr.eq.0) then
              do i=1,1000000
                read (nsnapf,9991,iostat=ioerr) tmpdata
                if (ioerr.ne.0) exit
                if (INDEX(tmpdata,'&efitin')/=0) exit
              enddo
            endif
            if (ioerr.eq.0) then
              do i=1,1000000
                write (neqdsk,9991) tmpdata
                read (nsnapf,9991,iostat=ioerr) tmpdata
                if (ioerr.ne.0) exit
                if (INDEX(tmpdata,'/')/=0) then
                  write (neqdsk,9991) tmpdata
                  exit
                endif
              enddo
            endif
            if (ioerr.eq.0) close (unit=nsnapf)
 9991       format (a)
          endif
        endif
!
        close(unit=neqdsk)
        if (filenm.eq.'0' .or. filenm.eq.'') then
          read (eqdsk,6700) prefix1,ishotime
        endif
        if (comfile.ne.'0' .and. comfile.ne.'') &
          write (nffile,4960) ishotime
      enddo
      if (filenm.ne.'0' .and. filenm.ne.'') go to 3046
      endif ! ioerr.eq.0
      if (comfile.ne.'0' .and. comfile.ne.'') then
        write (nffile,4962)
        close(unit=nffile)
      endif
      if (filenm.ne.'0' .and. filenm.ne.'') close(unit=nout)
 4042 format (1x,a42,1x,a3)
 4958 format ('#!/bin/csh -f')
 4960 format ('      runefit.sc k',a12)
 4962 format ('#',/,'exit')
 4970 format (2x,a,1x,a,1x,a,1x,a)
 6040 format (/,1x,'type shot #, start time(ms), time step(ms), steps' &
        ,'(<1001):')
 6600 format (/,1x,'good shot list file name ( 0=tty) ?')
 6610 format (/,1x,'command file name ( 0=none) ?')
 6620 format (a)
 6700 format (a1,a12)
20000 format (/,1x,'shot data not on disk')
      end


!**********************************************************************
!>
!!    Subroutine for writing K file in SNAP mode 2014/04/16
!!    
!!
!!    @param jtime : time index
!!
!!    @param kerror : error flag
!!
!**********************************************************************
      subroutine write_K2(jtime,kerror)
      use set_kinds
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
#if defined(USEMPI)
      include 'mpif.h'
#endif

      character table_s*72,eqdsk*20,header*42,fit_type*3
      integer*4 ltbdis,limitrss
      real*8,dimension(:),allocatable :: coils,expmp2, &
                denr,denv, tgamma,sgamma,rrrgam, &
                zzzgam,aa1gam,aa2gam, aa3gam,aa4gam,aa5gam, &
                aa6gam,aa7gam,tgammauncor
      real*8,dimension(:),allocatable :: tlibim,slibim,rrrlib
      character*82 snap_ext
      namelist/machinein/nsilds,nsilol,nfcoil,nrogow,nacoil,mfcoil,necoil,nvesel, &
           mpress,nesum,magpri67,magpri322,magprirdp,magudom,maglds,mse315,mse45, &
           mse15,mse1h,mse315_2,mse210,libim,nmsels,nnece,nnecein,neceo,nnnte, &
           ngam_vars,ngam_u,ngam_w,nlimit,nlimbd,nangle,ntangle,nfbcoil,mccoil, &
           micoil,ndata,nwwcur,nffcur,nppcur,nercur,ntime,ndim,kxiter,mqwant, &
           mbdry,mbdry1,nxtram,nxtlim,nco2v,nco2r,modef,modep,modew,kubics, &
           icycred_loopmax,nfourier
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
           kwwfnc,kwwknt,wwknt,wwtens,nbdry,rbdry,zbdry,vbit, nbdrymx, &
           ppbdry,kppbdry,pp2bdry,kpp2bdry, &
           ffbdry,kffbdry,ff2bdry,kff2bdry, &
           wwbdry,kwwbdry,ww2bdry,kww2bdry &
           ,ktear,kersil,iout,ixray,table_dir,input_dir,store_dir &
           ,kpphord,kffhord,keehord,psiecn,dpsiecn,fitzts,isolve &
           ,iplcout,imagsigma,errmag,saimin,errmagb,fitfcsum,fwtfcsum,efitversion &
           ,kwripre,ifindopt,tolbndpsi
      namelist/inwant/psiwant,vzeroj,nccoil,currc79,currc139,rexpan, &
           znose,sizeroj,fitdelz,relaxdz,errdelz,oldccomp,nicoil, &
           oldcomp,currc199,curriu30,curriu90, &
           curriu150,curril30,curril90,curril150,ifitdelz,scaledz
      namelist/ink/isetfb,ioffr,ioffz,ishiftz,gain,gainp,idplace &
           ,symmetrize,backaverage
      namelist/ins/tgamma,sgamma,fwtgam,rrrgam,zzzgam,aa1gam,aa2gam, &
                   aa3gam,aa4gam,aa5gam,aa6gam,aa7gam,msebkp, &
            msefitfun,mse_quiet,mse_spave_on,kwaitmse, &
            dtmsefull,mse_strict,t_max_beam_off,ok_30rt,ok_210lt, &
            mse_usecer,mse_certree,mse_use_cer330,mse_use_cer210, &
            tgammauncor,v30lt,v30rt,v210lt,v210rt
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
      data mcontr/35/,lfile/36/,ifpsi/0/
      data currn1/0.0/,currc79/0.0/,currc139/0.0/,currc199/0.0/ &
                      ,curriu30/0.0/,curriu90/0.0/,curriu150/0.0/ &
                      ,curril30/0.0/,curril90/0.0/,curril150/0.0/
      integer*4, intent(in) :: jtime
      integer*4, intent(inout) :: kerror


      ALLOCATE(coils(nsilop),expmp2(magpri), &
                denr(nco2r),denv(nco2v), &
                tgamma(nmtark),sgamma(nmtark),rrrgam(nmtark), &
                zzzgam(nmtark),aa1gam(nmtark),aa2gam(nmtark), &
                aa3gam(nmtark),aa4gam(nmtark),aa5gam(nmtark), &
                aa6gam(nmtark),aa7gam(nmtark),tgammauncor(nmtark))
      ALLOCATE(tlibim(libim),slibim(libim),rrrlib(libim))

      kerror = 0
      efitversion = efitvers//" "
!----------------------------------------------------------------
!--   recover the value of table_dir for mode 3 or 7           --
!----------------------------------------------------------------
      if (ishot.ge.112000) then
        ltbdis =ltbdir
        ltbdir = ltbdir-7
        table_s = table_dir
        table_dir = table_dir(1:ltbdir)  ! necessary??
      endif
!
      itime=time(jtime)
      timems=itime
      timeus=(time(jtime)-timems)*1000.
      itimeu=timeus
      siref=psirefs(jtime)
      do i=1,nsilop
        coils(i)=silopt(jtime,i)-siref
        fwtsi(i)=swtsi(i)
      enddo
      do i=1,magpri
        expmp2(i)=expmpi(jtime,i)
        fwtmp2(i)=swtmp2(i)
      enddo
      brspss=brsp
      brsp=0.0
      do i=1,nfcoil
        brsp(i)=fccurt(jtime,i)
        fwtfc(i)=swtfc(i)
      enddo
      do i=1,nesum
        ecurrt(i)=eccurt(jtime,i)
        fwtec(i)=swtec(i)
      enddo
      do i=1,nco2r
        denr(i)=denrt(jtime,i)
      enddo
      do i=1,nco2v
        denv(i)=denvt(jtime,i)
      enddo
      mmstark=0
      do i=1,nstark
        if (swtgam(i).gt.1.e-06_dp) mmstark=mmstark+1
      enddo
      if (mmstark.gt.0) then
        do i=1,nmtark
          tgamma(i)=tangam(jtime,i)
          tgammauncor(i)=tangam_uncor(jtime,i)
          sgamma(i)=siggam(jtime,i)
          rrrgam(i)=rrgam(jtime,i)
          zzzgam(i)=zzgam(jtime,i)
          aa1gam(i)=a1gam(jtime,i)
          aa2gam(i)=a2gam(jtime,i)
          aa3gam(i)=a3gam(jtime,i)
          aa4gam(i)=a4gam(jtime,i)
          aa5gam(i)=a5gam(jtime,i)
          aa6gam(i)=a6gam(jtime,i)
          aa7gam(i)=a7gam(jtime,i)
          fwtgam(i)=swtgam(i)
        enddo
      endif
      fwtcur=swtcur
      btor=bcentr(jtime)
      plasma=pasmat(jtime)
      rbdryss=rbdry
      zbdryss=zbdry
      nbdryss=nbdry
      rbdry=0.0
      zbdry=0.0
      if (fitzts.eq.'te'.and.ztserr(jtime)) then
        nbdry=1
        rbdry(1)=1.94_dp
        zbdry(1)=ztssym(jtime)+0.5_dp*ztswid(jtime)
      endif
      ppbdryss=ppbdry
      pp2bdryss=pp2bdry
      ffbdryss=ffbdry
      ff2bdryss=ff2bdry
      wwbdryss=wwbdry
      ww2bdryss=ww2bdry
      kppbdryss=kppbdry
      kpp2bdryss=kpp2bdry
      kffbdryss=kffbdry
      kff2bdryss=kff2bdry
      kwwbdryss=kwwbdry
      kww2bdryss=kww2bdry
!
      ppbdry=0.0
      pp2bdry=0.0
      ffbdry=0.0
      ff2bdry=0.0
      wwbdry=0.0
      ww2bdry=0.0
      kppbdry=0.0
      kpp2bdry=0.0
      kffbdry=0.0
      kff2bdry=0.0
      kwwbdry=0.0
      kww2bdry=0.0
      limitrss=limitr
      limitr=-limid
      vloop=vloopt(jtime)
      dflux=1.0e+03_dp*diamag(jtime)
      sigdlc=1.0e+03_dp*sigdia(jtime)
      pnbeam=pbinj(jtime)
      currn1=curtn1(jtime)
      currc79=curc79(jtime)
      currc139=curc139(jtime)
      currc199=curc199(jtime)
      curriu30=curiu30(jtime)
      curriu90=curiu90(jtime)
      curriu150=curiu150(jtime)
      curril30=curil30(jtime)
      curril90=curil90(jtime)
      curril150=curil150(jtime)
      itekt=itek
      mxitert=mxiter
      n1coilt=n1coil
      zeliptt=zelip
      itek=iteks
      mxiter=mxiters
      n1coil=n1coils
      zelip=zelipss
      if (ishot.gt.108281) n1coil = 0
      call getfnmu(itimeu,'k',ishot,itime,eqdsk)
      open(unit=neqdsk,status='old',file=eqdsk,iostat=ioerr, &
           recl=72,delim='APOSTROPHE')
      if (ioerr.eq.0) close(unit=neqdsk,status='delete')
!-----------------------------------------------------------------------
!--   Write K file                                                    --
!-----------------------------------------------------------------------
      open(unit=neqdsk,file=eqdsk,status='new', &
           recl=72,delim='quote')
!           recl=72,delim='APOSTROPHE')
      write (neqdsk,machinein)
      write (neqdsk,in1)
      write (neqdsk,inwant)
      if (isetfb.ne.0) write (neqdsk,ink)
      if (mmstark.gt.0) write (neqdsk,ins)
      if (kwaitmse.ne.0) write (neqdsk,ina)
      if (kfitece.gt.0) write (neqdsk,inece)
      if (keecur.gt.0) write (neqdsk,iner)
!-----------------------------------------------------------------------
!--   fitting type flag                                               --
!-----------------------------------------------------------------------
      if ((kprfit.gt.0).and.(mmstark.gt.0)) then
        fit_type = 'KIM'
      elseif (kprfit.gt.0) then
        fit_type = 'KIN'
      elseif (mmstark.gt.0) then
        fit_type = 'MSE'
      else
        fit_type = 'MAG'
      endif
!
      call db_header(ishot,itime,header)
      write (neqdsk,4042) header,fit_type
!----------------------------------------------------------------------
!--   Restore variables                                              --
!----------------------------------------------------------------------
      limitr=limitrss
      if (ishot.ge.112000) then
        ltbdir =ltbdis
        table_dir = table_s
      endif
      rbdry=rbdryss
      zbdry=zbdryss
      nbdry=nbdryss
      brsp=brspss
      ppbdry=ppbdryss
      pp2bdry=pp2bdryss
      ffbdry=ffbdryss
      ff2bdry=ff2bdryss
      wwbdry=wwbdryss
      ww2bdry=ww2bdryss
      kppbdry=kppbdryss
      kpp2bdry=kpp2bdryss
      kffbdry=kffbdryss
      kff2bdry=kff2bdryss
      kwwbdry=kwwbdryss
      kww2bdry=kww2bdryss
      itek=itekt
      mxiter=mxitert
      n1coil=n1coilt
      zelip=zeliptt
!---------------------------------------------------------------------
!--   Append SNAP file                                              --
!---------------------------------------------------------------------
      if (kdata.eq.7) then
        open(unit=nsnapf,status='old', &
             file=snap_file,iostat=ioerr)
        if (ioerr.eq.0) then
          snapfile=snap_file
        else
          open(unit=nsnapf,status='old', &
               file=input_dir(1:lindir)//snap_file,iostat=ioerr)
          if (ioerr.eq.0) then
            snapfile=input_dir(1:lindir)//snap_file
          else
            open(unit=nsnapf,status='old', &
                 file=snapextin)
            snapfile=snapextin
          endif
        endif
      else
        open(unit=nsnapf,status='old', &
             file='efit_snap.dat',iostat=ioerr)
        if (ioerr.eq.0) then
          snapfile='efit_snap.dat'
        else
          open(unit=nsnapf,status='old', &
               file= input_dir(1:lindir)//'efit_snap.dat')
          snapfile=input_dir(1:lindir)//'efit_snap.dat'
        endif
      endif
      if (appendsnap.eq.'K'.or.appendsnap.eq.'KG') then
        do i=1,1000000
          read (nsnapf,9991,iostat=ioerr) tmpdata
          if (ioerr.ne.0) exit
          if (INDEX(tmpdata,'&efitin')/=0) exit
        enddo
        if (ioerr.eq.0) then
          do i=1,1000000
            write (neqdsk,9991) tmpdata
            read (nsnapf,9991,iostat=ioerr) tmpdata
            if (ioerr.ne.0) exit
            if (INDEX(tmpdata,'/')/=0) then
              write (neqdsk,9991) tmpdata
              exit
            endif
          enddo
        endif
        if (ioerr.eq.0) close (unit=nsnapf)
 9991   format (a)
      endif
!
      close(unit=neqdsk)
      ltbdir=ltbdis
      ! TODO: which is best?
      table_dir = table_s
!      table_dir(1:ltbdis) = table_s(1:ltbdis)
!      table_dir = table_s(1:ltbdis) !compilers do not like...
 4042 format (1x,a42,1x,a3)
      end

!**********************************************************************
!>
!!    wtime writes out time history of plasma
!!    parameters in o-file format.         --
!!    
!!    
!!
!!    @param kwtime : number of time elements 
!!
!**********************************************************************      
      subroutine wtime(kwtime)
      use set_kinds
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      character*30 ofname, chname
      Character*28 xxtitle,yytitle,zztitle
      data czero/0.0/
!
      if (kwtime.le.1) return
      if (kwritime.ne.2) return
!
      xdum=0.0
      ydum=0.0
      itime00=time(1)
      call getfnm2('ot',ishot,itime00,ofname)
!
      ofname=ofname(1:14)//'_chi2mag'
      open(unit=74,status='old',file=ofname,iostat=ioerr)
      if (ioerr.eq.0) close(unit=74,status='delete')
      open(unit=74,status='new',file=ofname)
      xxtitle='Time(ms)'
      yytitle='Magnetic Chi Square'
      zztitle='EFIT'
      write (74,93024) xxtitle
      write (74,93024) yytitle
      write (74,93024) zztitle
      do i=1,kwtime
        if (kerrot(i).eq.0) then
          write (74,92924) time(i),tsaisq(i),xdum,xdum
        endif
      enddo
      close(unit=74)
!
      if (kstark.gt.0) then
        ofname=ofname(1:14)//'_chi2mse'
        open(unit=74,status='old',file=ofname,iostat=ioerr)
        if (ioerr.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=ofname)
        xxtitle='Time(ms)'
        yytitle='MSE Chi Square'
        zztitle='EFIT'
        write (74,93024) xxtitle
        write (74,93024) yytitle
        write (74,93024) zztitle
        do i=1,kwtime
          if (kerrot(i).eq.0) then
            write (74,92924) time(i),chi2gamt(i),xdum,xdum
          endif
        enddo
        close(unit=74)
      endif
!
      if (mmbmsels.gt.0) then
        ofname=ofname(1:14)//'_chi2mls'
        open(unit=74,status='old',file=ofname,iostat=ioerr)
        if (ioerr.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=ofname)
        xxtitle='Time(ms)'
        yytitle='MSE-LS Chi Square'
        zztitle='EFIT'
        write (74,93024) xxtitle
        write (74,93024) yytitle
        write (74,93024) zztitle
        do i=1,kwtime
          if (kerrot(i).eq.0) then
            write (74,92924) time(i),chi2mls(i),xdum,xdum
          endif
        enddo
        close(unit=74)
!
        do j=1,nmsels
          if (sinmls(j).gt.0.0) then
            ichan00=j
            call getfnm2('ot',ishot,ichan00,chname)
            ofname=ofname(1:14)//'_echan'//chname(13:14)
            open(unit=74,status='old',file=ofname,iostat=ioerr)
            if (ioerr.eq.0) close(unit=74,status='delete')
            open(unit=74,status='new',file=ofname)
            xxtitle='Time(ms)'
            yytitle='Measured B MSE-LS (T)'
            zztitle='MSE-LS Channel '//chname(13:14)
            write (74,93024) xxtitle
            write (74,93024) yytitle
            write (74,93024) zztitle
            do i=1,kwtime
              if (kerrot(i).eq.0) then
                iges=i
                ydum=sbmselt(iges,j)
                if (swtbmselt(iges,j).gt.1.e-06_dp) ydum=ydum/ &
                                                       swtbmselt(iges,j)
                write (74,92924) time(i),bmselt(iges,j),xdum,ydum
              endif
            enddo
            close(unit=74)
          endif
        enddo
!
        do j=1,nmsels
          if (sinmls(j).gt.0.0) then
            ichan00=j
            call getfnm2('ot',ishot,ichan00,chname)
            ofname=ofname(1:14)//'_cchan'//chname(13:14)
            open(unit=74,status='old',file=ofname,iostat=ioerr)
            if (ioerr.eq.0) close(unit=74,status='delete')
            open(unit=74,status='new',file=ofname)
            xxtitle='Time(ms)'
            yytitle='Computed B MSE-LS (T)'
            zztitle='MSE-LS Channel '//chname(13:14)
            write (74,93024) xxtitle
            write (74,93024) yytitle
            write (74,93024) zztitle
            do i=1,kwtime
              if (kerrot(i).eq.0) then
                iges=i
                write (74,92924) time(i),cmmls(iges,j),xdum,xdum
              endif
            enddo
            close(unit=74)
          endif
        enddo
!
      endif
      return
92924 format (4(1pe12.5,1x))
93024 format (a28)
      end
