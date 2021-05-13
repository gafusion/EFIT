      program efitd
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD EQUILIBRIUM ANALYSIS                      **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          efit is the main driver for equilibrium analysis.       **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1) L.L. Lao, H. St. John, R.D. Stambaugh,              **
!**              A.G. Kellman, and W. Pfeiffer, Nuclear Fusion       **
!**              25 (1985) 1611.                                     **
!**          (2) L.L. Lao, H. St. John, R.D. Stambaugh, and          **
!**              W. Pfeiffer, Nuclear Fusion 25 (1985) 1421.         **
!**          (3) L.L. Lao, J.R. Ferron, R.J. Groebner, W. Howl,      **
!**              H. St. John, E.J. Strait, and T.S. Taylor           **
!**              Nuclear Fusion 30 (1990) 1035.                      **
!**          (4) L.L. Lao and T.H. Jensen Nuclear Fusion             **
!**              31 (1991) 1909.                                     **
!**          (5) L.L. Lao, H. St. John, et al, Fusion Sci. Technol.  **
!**              48 (2005) 968.                                      **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**          13/02/86..........revised                               **
!**          06/12/91..........revised                               **
!**          15/02/93..........revised                               **
!**          11/11/93..........revised                               **
!**          03/01/94..........revised to include rotation           **
!**          97/04/28..........ER from MSE                           **
!**          97/10/01..........took off the DIIID specific data part **
!**          98/05/01..........Add hyperbolic tangent representation **
!**          98/05/12..........Q.P. added Toroidal X-ray overlay     **
!**          99/09/15..........Q.P. use IOUT to determine whether to **
!**                            write rstarkxx.dat(8) or esave.dat(16)**
!**                            default to NO.                        **
!**        2001/01/18..........Q.P. added saipre2 to preserve saipre **
!**        2005/03/25..........New magnetic uncertainty              **
!**        2006/01/12..........New Magnetic Uncertainty              **
!**        2007/08/01..........Mag. Uncnty. namelist added to k file **
!**        2012/05/17..........MPI, SNAP revision                    **
!**        2020/09/18..........R.S. Bug fix, changed mpi_abort to    **
!**                            mpi_finalize at end to allow all      **
!**                            processes to complete.                **
!**        2020/09/18..........R.S. changed "shape" to "shapesurf",  **
!**                            shape is an intrinsic procedure name  **
!**        2020/09/18..........R.S. changed some Hollerith to quotes **
!**                                                                  **
!**********************************************************************
     use commonblocks
     use set_kinds
     include 'eparmdud129.inc'
     include 'modules2.inc'
     include 'modules1.inc'
     implicit integer*4 (i-n), real*8 (a-h,o-z)
! MPI >>>
#ifdef USEMPI
     include 'mpif.h'
#endif
! MPI <<<
     data kwake/0/
     parameter (krord=4,kzord=4)
     character inp1*4,inp2*4
     integer :: nargs, iargc, finfo, kerror, terr

     integer :: iend1, iend2
     character*80 :: cmdline

     kerror = 0
!------------------------------------------------------------------------------
!--   Set paths
!------------------------------------------------------------------------------ 
      call set_expath
      ! If environment variable exists, then override values from set_expath
      call getenv("link_efit",link_efitx)
      call getenv("link_store",link_storex)
      if (link_efitx(1:1).ne.' ') then
          table_dir=trim(link_efitx)//'green/'
          input_dir=trim(link_efitx)
      endif
      if (link_storex(1:1).ne.' ')  store_dir=trim(link_storex)           
!------------------------------------------------------------------------------
!--   MPI if we have it
!------------------------------------------------------------------------------ 
! MPI >>>
#if defined(USEMPI)
! Initialize MPI environment
      call MPI_INIT_THREAD(MPI_THREAD_SINGLE,terr,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
#else
      rank  = 0
#endif
! MPI <<<

      ! Set global constants for each rank
      call set_constants()

!----------------------------------------------------------------------
!-- Read in grid size from command line and set global variables     --
!-- ONLY root process reads command-line arguments                   --
!----------------------------------------------------------------------
      nw = 0
      nh = 0
      if (rank == 0) then
       nargs = iargc()
! Using mpirun command so will have different number of arguments than serial case
#if defined(LF95)
       call getcl(cmdline)
       cmdline = adjustl(cmdline)
       iend1 = scan(cmdline,' ')
       inp1 = cmdline(1:iend-1)//'  '
       cmdline = adjustl(cmdline(iend1:len(cmdline)))
       iend2 = scan(cmdline,' ')
       inp2 = cmdline(1:iend2-1)//'  '
#else
       call getarg(1,inp1)
       call getarg(2,inp2)
#endif
       read (inp1,'(i4)',err=9876) nw
       read (inp2,'(i4)',err=9876) nh
9876   continue
       if (nh == 0) nh = nw

      endif

! MPI >>>
#if defined(USEMPI)
! Distribute command-line information (meaning grid dimensions) to all processes if necessary
      if (nproc > 1) then
       call MPI_BCAST(nw,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(nh,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      endif
#endif
! MPI <<<
      if (nw == 0 .or. nh == 0) then
        if (rank == 0) then
          call errctrl_msg('efitd','Must specify grid dimensions as arguments')
        endif
! MPI >>>
#if defined(USEMPI)
        deallocate(dist_data,dist_data_displs,fwtgam_mpi)
        call mpi_finalize(ierr)
#endif
        STOP
! MPI <<<
      endif

      IF (nw .le. 129) THEN
        npoint=800
      ELSE
        npoint=3200
      ENDIF
       nwnh=nw*nh
       nh2=2*nh
       nwrk=2*(nw+1)*nh
!      nwwf=2*nw
       nwwf=3*nw
       nwf=nwwf
       kubicx = 4
       kubicy = 4
       lubicx = nw - kubicx + 1
       lubicy = nh - kubicy + 1
       kujunk = kubicx*kubicy*lubicx*lubicy
       boundary_count=2*nh+2*(nw-2)
       lr0=nw-krord+1
       lz0=nh-kzord+1
       nxtrap=npoint
       mfila = 10

      
      call read_efitin
      call inp_file_ch(nw,nh,ch1,ch2)

      call get_opt_input(ktime)
      ntime = ktime
      call get_eparmdud_defaults()

      if (kdata==2) then
        call read_shot(ifname(1))     !this assume machine is always the same
        call read_eparmdud(ifname(1)) !this assume machine is always the same
      elseif(kdata==7) then
        call read_eparmdud('efit_snap.dat_'//adjustl(snapext_in))
      else
        call read_eparmdud('efit_snap.dat')
      endif
      call get_eparmdud_dependents()

!----------------------------------------------------------------------
!-- Global Allocations                                               --
!----------------------------------------------------------------------
  
      include 'global_allocs.f90'
      call set_ecom_mod1_arrays
      call set_ecom_mod2_arrays

!----------------------------------------------------------------------
!-- get data                                                         --
!----------------------------------------------------------------------
#if defined(USEMPI)
! Arrays can only be allocated after MPI has been initialized because dimension is # of processes
      allocate(dist_data(nproc),dist_data_displs(nproc),fwtgam_mpi(nstark,nproc))
#endif
    !call set_table_dir
      !call efit_read_tables
      print *, 'Entering getsets'
  20  call getsets(ktime,kwake,mtear,kerror)
      print * ,'exiting getsets'
! MPI >>>
#if defined(USEMPI)
      if (nproc > 1) then
        call MPI_ALLREDUCE(kerror,MPI_IN_PLACE,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
      endif
      if (kerror.gt.0) then
        call errctrl_msg('efitd','Aborting due to fatal error in getsets')
        call mpi_abort(MPI_COMM_WORLD,ierr) ! kill all processes, something is wrong with the setup.
      endif
#else
      if (kerror.gt.0) then
        stop
      end if
#endif

! MPI <<<

! Looping (below) on the number of time slices depends on the number of ranks.
! Time slices are assigned to ranks in groups as follows:
! slice: 100, 120, 140, 160, 180, 200, 220, 240, ...
! rank:    0,   0,   0,   1,   1,   1,   2,   2, ...
! However, not all ranks necessarily have the same number of slices.
! The number of slices per rank = floor(nslices/nranks), BUT
!   ranks 0 to N also have ONE additional slice, where
!   N = nslices - nranks*floor(nslices/nranks) - 1
!   nslices = total number of slices for all ranks

!----------------------------------------------------------------------
!-- start simulation for KTIME time slices per rank                  --
!----------------------------------------------------------------------
      k=0
  100 k=k+1
        ks=k ! ks=1,2,3... in serial, but ks=1,1,1,... in parallel
!----------------------------------------------------------------------
!--  set up data                                                     --
!----------------------------------------------------------------------
        
        if(idebug>=2) write(6,*) ' Entering prtoutheader subroutine'
        call prtoutheader()
        if(idebug>=2) write(6,*) ' Entering data_input subroutine'

        call data_input(ks,iconvr,ktime,mtear,kerror)

        if(idebug>=2) write(6,*) ' Entering errctrl_setstate'
        call errctrl_setstate(rank,time(ks))
        if (kerror.gt.0) go to 500
        if (iconvr.lt.0) go to 500
        if (kautoknt .eq. 1) then
           call autoknot(ks,iconvr,ktime,mtear,kerror)
        else
!----------------------------------------------------------------------
!--  initialize current profile                                      --
!----------------------------------------------------------------------

           if(idebug>=2) write(6,*) 'Entering inicur subroutine'
           call inicur(ks)
!----------------------------------------------------------------------
!--  get equilibrium                                                 --
!----------------------------------------------------------------------

           if(idebug>=2) write(6,*) 'Entering fit subroutine'
           call fit(ks,kerror)
           if (kerror.gt.0) go to 500
        endif
!----------------------------------------------------------------------
!--  post processing for graphic and text outputs                    --
!----------------------------------------------------------------------
        if(idebug>=2) write(6,*) 'Entering shapesurf'
        call shapesurf(ks,ktime,kerror)
        if (kerror.gt.0) go to 500
!DEPRECATED        if (mtear.ne.0) call tearing(ks,mtear,kerror)
        if (kerror.gt.0) go to 500
        if (idebug /= 0) write (6,*) 'Main/PRTOUT ks/kerror = ', ks, kerror
        call prtout(ks)
        if ((kwaitmse.ne.0).and.(kmtark.gt.0)) call fixstark(-ks,kerror)
!----------------------------------------------------------------------
!--  write A and G EQDSKs                                            --
!----------------------------------------------------------------------
        if (idebug /= 0) write (6,*) 'Main/WQEDSK ks/kerror = ', ks, kerror
        call weqdsk(ks)
        if (iconvr.ge.0) then
           call shipit(ktime,ks,ks)
!DEPRECATED           call wtear(mtear,ks)
        endif
#if HAVE_NETCDF
        if (idebug /= 0) write (6,*) 'Main/wmeasure ks/kerror = ', ks, kerror
        call wmeasure(ktime,ks,ks,1)
#endif
!----------------------------------------------------------------------
! -- write Kfile if needed                                           --
!----------------------------------------------------------------------
      if (kdata.eq.3 .or. kdata.eq.7) then
         if (write_Kfile) then
           call write_K2(ks,kerror)
         endif
      endif
  500 if (k.lt.ktime) then
        kerrot(ks)=kerror
        go to 100
      endif
      if (kwake.ne.0) go to 20
#if HAVE_NETCDF
      call wmeasure(ktime,1,ktime,2)
#endif
      call wtime(ktime)

! MPI >>>
#if defined(USEMPI)
      ! Finalize MPI
      if (allocated(dist_data)) deallocate(dist_data)
      if (allocated(dist_data_displs)) deallocate(dist_data_displs)
      if (allocated(fwtgam_mpi)) deallocate(fwtgam_mpi)
      call errctrl_msg('efitd','Done processing',3)
      call mpi_finalize(ierr)
#endif
      stop
! MPI <<<
      end program efitd

!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          betali computes betas and li.                           **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**          22/07/96 revised  by Q.Peng to add the surface area of  **
!**                            the last closed flux surface (psurfa) **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      subroutine betali(jtime,rgrid,zgrid,idovol,kerror)
      use commonblocks,only: c,wk,copy,bkx,bky,cw,wkw,copyw,bwx, &
                  bwy,sifprw,bwprw,cwprw,dwprw,sfprw,sprwp
      include 'eparmdud129.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      common/cww/lwx,lwy
      dimension pds(6),rgrid(*),zgrid(*)
      real*8,dimension(:),allocatable :: worksi,workrm,bwork, &
             cwork,dwork,x,y,dpleng
      dimension xsier(nercur)
      integer, intent(inout) :: kerror
      data inorm/3/,ibtcal/2/

      kerror = 0

      ALLOCATE(worksi(nw),workrm(nw),bwork(nw), &
         cwork(nw),dwork(nw),x(nw),y(nh),dpleng(npoint))

      if (ivacum.gt.0) return
      sumbp2=0.0
      select case (licalc)
      case (1)
        go to 20
      case (2)
        go to 60
      end select

   20 continue
      call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
      do i=1,nw
        do j=1,nh
          kk=(i-1)*nh+j
          call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n333)
          bpolsq=(pds(2)**2+pds(3)**2)/rgrid(i)
          sumbp2=sumbp2+bpolsq*www(kk)
        enddo
      enddo
      sumbp2=sumbp2*twopi*darea
      go to 90
   60 continue
      do 70 kk=1,nwnh
      sumbp2=sumbp2+(psi(kk)-psibry)*pcurrt(kk)*www(kk)
   70 continue
      sumbp2=sumbp2*twopi*tmu*twopi
   90 continue
!
      psurfa(jtime)=0.0
      plengt(1)=0.0
      do 100 i=1,nfound-1
        ip1=i+1
        dli=sqrt((xout(ip1)-xout(i))**2+(yout(ip1)-yout(i))**2)
        plengt(ip1)=plengt(i)+dli
        dpleng(i)=dli
        psurfa(jtime)=psurfa(jtime)+dli*0.5_dp*(xout(i)+xout(ip1))
  100 continue
      psurfa(jtime)=psurfa(jtime)*twopi
!
      s1(jtime)=0.0
      s2(jtime)=0.0
      s3(jtime)=0.0
      qout(jtime)=0.0
      sbp=0.0
      sbpi=0.0
      sumbzz=0.0
      sumr2=0.0
      sumz=0.0
      do 500 i=1,nfound-1
        ip1=i+1
        axout=(xout(i)+xout(ip1))/2.
        ayout=(yout(i)+yout(ip1))/2.
        erhor=axout-rcentr
        erhoz=ayout
        rho=sqrt(erhor**2+erhoz**2)
        enz=(xout(ip1)-xout(i))/2.
        enr=-(yout(ip1)-yout(i))/2.
        erhor=erhor/rho
        erhoz=erhoz/rho
        rlnr=sqrt(enr**2+enz**2)
        enr=enr/rlnr
        enz=enz/rlnr
        sdnrho=enr*erhor+enz*erhoz
        abpol=(bpol(ip1)+bpol(i))/2.
        sumbp=dpleng(i)*abpol
        sumbzz=sumbzz+(bpolz(ip1)+bpolz(i))/2.*(yout(ip1)-yout(i))
        if (kbetapr.eq.0) then
          sum=dpleng(i)*abpol**2*axout
        else
          if (kbetapr.eq.1) pxxxxx=prbdry
          if (kbetapr.eq.2) pxxxxx=pressb
          sum=dpleng(i)*(abpol**2+4.0*pi*tmu*pxxxxx)*axout
        endif
        s1(jtime)=s1(jtime)+sum*sdnrho*rho
        s2(jtime)=s2(jtime)+sum*enr
        s3(jtime)=s3(jtime)+sum*ayout*enz
        sbp=sbp+sumbp
        sumr2=sumr2+sumbp*axout**2
        sumz=sumz+sumbp*ayout
        dlbpi=dpleng(i)/abpol
        sbpi=sbpi+dlbpi
        qout(jtime)=qout(jtime)+dlbpi/axout**2
  500 continue
      bp2flx=sbp/sbpi
      if (inorm.eq.2) bp2flx=2.*rcentr/vout(jtime)*(pi*tmu* &
                                cpasma(jtime))**2*1.0E+06_dp
      if (inorm.eq.3) bp2flx=(tmu*2.0*pi*cpasma(jtime)/ &
                             plengt(nfound))**2
      if (inorm.eq.4) bp2flx=2.*(tmu*cpasma(jtime)/aout(jtime))**2 &
                        /(eout(jtime)**2+1.)*1.0e+04_dp
      bpolav(jtime)=sqrt(bp2flx)
      rcurrt(jtime)=sqrt(sumr2/twopi/tmu/abs(cpasma(jtime)))*100.
      zcurrt(jtime)=sumz/twopi/tmu/abs(cpasma(jtime))*100.
      const=twopi/vout(jtime)*1.0e+06_dp/bp2flx
      s1(jtime)=const*s1(jtime)
      s2(jtime)=const*rcentr*s2(jtime)
      s3(jtime)=const*s3(jtime)
      sumbzz=abs(sumbzz)/tmu/abs(cpasma(jtime))/twopi
      rcurrm=rcentr
      ali(jtime)=sumbp2/vout(jtime)/bp2flx*1.0e+06_dp
      ali3(jtime)=(tmu*2.0*pi*cpasma(jtime))**2*rout(jtime)/200.0
      ali3(jtime)=sumbp2/ali3(jtime)
      betap(jtime)=s1(jtime)/4.+s2(jtime)/4.*(1.+rcurrm/rcentr) &
           -ali(jtime)/2.
      betat(jtime)=betap(jtime)*bp2flx/bcentr(jtime)**2*100.
      betat(jtime)=betat(jtime)*(rout(jtime)/100./rcentr)**2
      qout(jtime)=qout(jtime)*abs(bcentr(jtime))*rcentr/twopi
      wplasm(jtime)=1.5_dp*betap(jtime)*bp2flx/2./tmu/2./pi*vout(jtime) &
                     /1.0e+06_dp
!---------------------------------------------------------------------
!-- calculations of current moment y2                               --
!---------------------------------------------------------------------
      sumy2=0.0
      rcccc=rcurrt(jtime)/100.
      rmajz0(nw)=0.0
      do 520 i=1,nfound-1
        ip1=i+1
        axout=(xout(i)+xout(ip1))/2./rcccc
        ayout=(yout(i)+yout(ip1))/2./rcccc
        abpol=(bpol(ip1)+bpol(i))/2.
        sumbp=dpleng(i)*abpol
        axout=axout**2
        ayout=ayout**2
        f222=axout*(axout-4.*ayout)
        ff222=f222-2.*axout+1.
        sumy2=sumy2+sumbp*ff222
        if (yout(i)*yout(ip1).le.0.0.and.xout(i).gt.rcentr) then
          rmajz0(nw)=xout(i)-yout(i)*(xout(ip1)-xout(i))/ &
             (yout(ip1)-yout(i))
        endif
  520 continue
      yyy2(jtime)=sumy2/(tmu*2.0*pi*cpasma(jtime))/4. &
                  *(rcurrt(jtime)/aout(jtime))**2
!
      dsi=(psibry-simag)/(nw-1)
      mx=(rmaxis-rgrid(1))/drgrid+1
      my=(zmaxis-zgrid(1))/dzgrid+1
      mkk=(mx-1)*nh+my+1
      sicut=psi(mkk)
      volp(nw)=vout(jtime)/1.0e+06_dp
      volp(1)=0.0
      if (abs(zmaxis).le.0.001_dp) then
         rmajz0(1)=rmaxis
      else
         rmajz0(1)=0.0
      endif
      r1surf(1)=1./rmaxis
      r2surg(1)=1./rmaxis**2
      r1surf(nw)=r1bdry
      r2surg(nw)=r2bdry
      bpolss(1)=0.0
      bpolss(nw)=bpolav(jtime)
!----------------------------------------------------------------------
!--  rotational terms                                                --
!----------------------------------------------------------------------
      if (kvtor.gt.0) then
        sumprt=0.0
        sumprw=0.0
        n1set=1
        ypsi=0.5_dp
        prew0=pwcur4(n1set,ypsi,kwwcur)
        pres0=prcur4(n1set,ypsi,kppcur)
        n1set=0
        do i=1,nw
          pwprim(i)=sprwp(i)
          pressw(i)=sfprw(i)
        enddo
        do i=1,nw
        do j=1,nh
           kk=(i-1)*nh+j
           ypsi=xpsi(kk)
           pres0=prcur4(n1set,ypsi,kppcur)-prbdry
           prew0=pwcur4(n1set,ypsi,kwwcur)-preswb
           if (kvtor.eq.1) then
             presst(kk)=pres0+prew0*rgrvt(i)
           elseif (kvtor.eq.2) then
             if (abs(pres0).gt.1.e-10_dp) then
               pwop0=prew0/pres0
               pwp0r2=pwop0*rgrvt(i)
             else
               pwop0=0.0
               pwp0r2=0.0
             endif
             presst(kk)=pres0*(1.+0.5_dp*pwp0r2**2)+prew0*rgrvt(i)
           elseif (kvtor.eq.11.or.kvtor.eq.3) then
             if (abs(pres0).gt.1.e-10_dp) then
               pwop0=prew0/pres0
               ptop0=exp(pwop0*rgrvt(i))
             else
               ptop0=1.0
             endif
             presst(kk)=pres0*ptop0
           endif
           sumprt=sumprt+rgrid(i)*presst(kk)*www(kk)
           prewx=prew0*(rgrid(i)/rvtor)**2
           sumprw=sumprw+rgrid(i)*prewx*www(kk)
        enddo
        enddo
        sumprt=sumprt*darea*twopi
        sumprw=sumprw*darea*twopi
        call sets2d(presst,cw,rgrid,nw,bwx,lwx,zgrid,nh,bwy,lwy,wkw,ier)
      endif
!
      sumpre=0.0
      rzzmax(1)=rmaxis
      zzmax(1)=zmaxis
      rhovn(1)=0.0
      nnn=1
      d11=30.
      d22=0.03_dp
      d33=0.01_dp
      nzz=0
      zzz=0.0
      do 800 i=2,nw-1
        ii=nw-i+1
        siii=(i-1)*dsi
        siwant=psibry-siii
        siii=1.0_dp-1.0_dp/(nw-1)*(i-1)
        if (idovol.gt.1) go to 790
        rzzmax(ii)=-99.0
        call surfac(siwant,psi,nw,nh,rgrid,zgrid,bpol,bpolz,nfind, &
                    npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nnn, &
                    rmaxis,zmaxis,negcur,kerror)
        if (kerror.gt.0) return
        if (nfind.le.40.and.icntour.eq.0) then
        if (idebug >= 2) write (6,*) ' SHAPE/BETALI kerror,i,nfind = ', & 
                            kerror,i,nfind
        call cntour(rmaxis,zmaxis,siwant,xcmin,xcmax,ycmin,ycmax, &
                    yxcmin,yxcmax,xycmin,xycmax,d11,drgrid,d22, &
                    d33 ,d33 ,xmin,xmax,ymin,ymax,nzz,iautoc, &
                    bpol,bpolz,nfind,rgrid,nw,zgrid,nh, &
                    c,n222,nh2,nttyo,npoint, &
                    negcur,bkx,lkx,bky,lky,kerror)
        if (idebug >= 2) write (6,*) ' BETALI/CNTOUR kerror,nfind = ', &
                            kerror,nfind
        if (kerror /= 0) return
        endif
        if (nfind.lt.10) go to 750
        r2surf(ii)=0.0
        r1surf(ii)=0.0
        volp(ii)=0.0
        rhovn(ii)=0.0
        xym=bpol(1)*bpolz(1)
!------------------------------------------------------------------
!-- integration over z from 0 to bpolz                           --
!------------------------------------------------------------------
        xww=bpol(1)
        dyww=bpolz(1)/(nh-1)
        bpolzs=0.5*fpol(ii)
        do kz=2,nh-1
           yww=dyww*(kz-1)
          call seva2d(bkx,lkx,bky,lky,c,xww,yww,pds,ier,n111)
           ypsz=(simag-pds(1))/sidif
           fnow=seval(nw,ypsz,sigrid,fpol,bbfpol,ccfpol,ddfpol)
           bpolzs=bpolzs+fnow
        enddo
        yww=0.0
          call seva2d(bkx,lkx,bky,lky,c,xww,yww,pds,ier,n111)
        ypsz=(simag-pds(1))/sidif
        fnow=seval(nw,ypsz,sigrid,fpol,bbfpol,ccfpol,ddfpol)
        bpolzs=bpolzs+fnow*0.5
        yoxm=bpolzs/bpol(1)*dyww
!
        izzmax=1
        zzmax(ii)=bpolz(1)
        rmajz0(ii)=0.0
        ssdlbp=0.0
        sdlbpol=0.0
        circum=0.0
        do 700 k=2,nfind
          km1=k-1
          xyp=bpol(k)*bpolz(k)
!------------------------------------------------------------------
!-- integration over z from 0 to bpolz                           --
!------------------------------------------------------------------
          xww=bpol(k)
          dyww=bpolz(k)/(nh-1)
          bpolzs=0.5*fpol(ii)
          do kz=2,nh-1
           yww=dyww*(kz-1)
          call seva2d(bkx,lkx,bky,lky,c,xww,yww,pds,ier,n111)
           ypsz=(simag-pds(1))/sidif
           fnow=seval(nw,ypsz,sigrid,fpol,bbfpol,ccfpol,ddfpol)
           bpolzs=bpolzs+fnow
          enddo
          yww=0.0
          call seva2d(bkx,lkx,bky,lky,c,xww,yww,pds,ier,n111)
          ypsz=(simag-pds(1))/sidif
          fnow=seval(nw,ypsz,sigrid,fpol,bbfpol,ccfpol,ddfpol)
          bpolzs=bpolzs+fnow*0.5
          yoxp=bpolzs/bpol(k)*dyww
!
          dx=bpol(k)-bpol(km1)
          volp(ii)=volp(ii)+(xyp+xym)/2.0*dx
          rhovn(ii)=rhovn(ii)+(yoxp+yoxm)/2.0*dx
          xym=xyp
          yoxm=yoxp
          dli=sqrt((bpol(k)-bpol(km1))**2+(bpolz(k)-bpolz(km1))**2)
          xww=(bpol(k)+bpol(km1))/2.
          yww=(bpolz(k)+bpolz(km1))/2.
          call seva2d(bkx,lkx,bky,lky,c,xww,yww,pds,ier,n333)
          bpolsq=sqrt(pds(2)**2+pds(3)**2)/xww
          sdlbp=dli/bpolsq
          sdlbpol=sdlbpol+dli*bpolsq
          circum=circum+dli
          r2surf(ii)=r2surf(ii)+sdlbp/xww**2
          r1surf(ii)=r1surf(ii)+sdlbp/xww
          ssdlbp=ssdlbp+sdlbp
          if (bpolz(k).gt.zzmax(ii)) then
            izzmax=k
            zzmax(ii)=bpolz(k)
          endif
          if (bpolz(km1)*bpolz(k).le.0.0.and.bpol(km1).gt.rcentr) then
            rmajz0(ii)=bpol(km1)-bpolz(km1)*(bpol(k)-bpol(km1))/ &
             (bpolz(k)-bpolz(km1))
          endif
  700   continue
        r1surf(ii)=r1surf(ii)/ssdlbp
        r2surg(ii)=r2surf(ii)/ssdlbp
        bpolss(ii)=sdlbpol/circum
        call qfit(n333,bpol(izzmax-1),bpol(izzmax),bpol(izzmax+1), &
                  bpolz(izzmax-1),bpolz(izzmax),bpolz(izzmax+1), &
                  zaaa,zbbb,zccc,ierr)
        if (ierr.ne.0) then
          kerror = 1
          return
        end if
        rzzmax(ii)=-zbbb/2./zaaa
        zzmax(ii)=zaaa*rzzmax(ii)**2+zbbb*rzzmax(ii)+zccc
        volp(ii)=abs(volp(ii))*twopi
        rhovn(ii)=rhovn(ii)/rhovn(nw)
        go to 790
  750   continue
        bpolss(ii)=0.0
        rmajz0(ii)=0.0
        sisi=simag-siwant
        bbsi=sqrt(sisi/siaz)
        aasi=sqrt(sisi/siar)
        volp(ii)=twopi*rmaxis*pi*aasi*bbsi
        rhovn(ii)=pi*aasi*bbsi/rmaxis/rhovn(nw)*fpol(1)
!
        eesi=bbsi/aasi
        select case (icurrt)
        case (1)
          go to 760
        case (2)
          go to 775
        case (3)
          go to 780
        case (4)
          go to 770
        case (5)
          go to 775
        end select
  760   continue
        rdiml=rmaxis/srma
        cjmaxi=cratio*(sbeta*rdiml+2.*salpha/rdiml)/darea
        if (kvtor.gt.0) then
          cjmaxi=cjmaxi+cratio/darea*sbetaw*rdiml*(rdiml**2-1.)
        endif
        go to 785
  770   continue
        rdiml=rmaxis/rzero
        cjmaxi=cratio/darea*(rdiml+rbetap/rdiml)
        cjmaxi=cratio/darea*(rdiml+rbetap/rdiml)
        if (kvtor.eq.1) then
          rgmvt=(rmaxis/rvtor)**2-1.
          cjmaxi=cjmaxi+cratio/darea*rdiml*rbetaw*rgmvt
        elseif (kvtor.eq.11) then
           ypsm=0.0
           n1set=1
           pres0=prcur4(n1set,ypsm,kppcur)
           prew0=pwcur4(n1set,ypsm,kwwcur)
           rgmvt=(rmaxis/rvtor)**2-1.
           pwop0=prew0/pres0
           ptop0=exp(pwop0*rgmvt)
           pp0= 1.-pwop0*rgmvt
           ppw=rbetaw*rgmvt
           cjmaxi=cjmaxi+(pp0+ppw)*rdiml*ptop0
        endif
        go to 785
  775   continue
        cjmaxi=(rmaxis*ppcurr(zzz,kppcur) &
             +fpcurr(zzz,kffcur)/rmaxis)*cratio/darea
        if (kvtor.gt.0) then
          cjmaxi=0.0
          do j=1,kppcur
            cjmaxi=rjjjx(j)*brsp(nfcoil+j)+cjmaxi
          enddo
          do j=1,kffcur
            cjmaxi=rjjfx(j)*brsp(nfcoil+kppcur+j)+cjmaxi
          enddo
          do j=1,kwwcur
            cjmaxi=rjjwx(j)*brsp(nfcoil+kpcurn+j)+cjmaxi
          enddo
          cjmaxi=cjmaxi/darea
        endif
        go to 785
  780   continue
  785   continue
        r2surf(ii)=(eesi**2+1.)/rmaxis**2/tmu/cjmaxi
        r1surf(ii)=1./rmaxis
  790   continue
        sumpre=sumpre+volp(ii)*pprime(ii)
  800 continue
      sumpre=sumpre+volp(nw)*pprime(nw)/2.
      if (ibtcal.le.1) return
      if (kbetapr.eq.0) then
        sumpre=-sumpre*dsi/volp(nw)
      else
        if (kbetapr.eq.1) pxxxxx=prbdry
        if (kbetapr.eq.2) pxxxxx=pressb
        sumpre=-sumpre*dsi/volp(nw)+pxxxxx
      endif
      betat(jtime)=sumpre*2.0*twopi*tmu/bcentr(jtime)**2
      betap(jtime)=sumpre*2.0*twopi*tmu/bp2flx
      betat(jtime)=100.*betat(jtime)*(rout(jtime)/100./rcentr)**2
      wplasm(jtime)=1.5_dp*betap(jtime)*bp2flx/2./tmu/2./pi*vout(jtime) &
                     /1.0e+06_dp
      pasman=cpasma(jtime)/1.e4_dp/aout(jtime)/abs(bcentr(jtime))
      pasman=pasman*rout(jtime)/100./rcentr
      betatn=betat(jtime)/pasman
      dmui=1.0e+06_dp*diamag(jtime)*4.*pi*bcentr(jtime)*rcentr &
            /bp2flx/vout(jtime)
      betapd(jtime)=s1(jtime)/2.+s2(jtime)/2.*(1.-rcurrm/rcentr)-dmui
      betatd(jtime)=betapd(jtime)*bp2flx/bcentr(jtime)**2*100.
      betatd(jtime)=betatd(jtime)*(rout(jtime)/100./rcentr)**2
      wplasmd(jtime)=1.5_dp*betapd(jtime)*bp2flx/2./tmu/2./pi*vout(jtime) &
                    /1.0e+06_dp
!-----------------------------------------------------------------------
!-- rotational terms                                                  --
!-----------------------------------------------------------------------
      if (kvtor.gt.0) then
        sumprt=sumprt/volp(nw)
        sumprw=sumprw/volp(nw)
        betat(jtime)=sumprt/sumpre*betat(jtime)
        betap(jtime)=sumprt/sumpre*betap(jtime)
        wplasm(jtime)=sumprt/sumpre*wplasm(jtime)
        betapw(jtime)=sumprw/sumprt*betap(jtime)
        betatw(jtime)=sumprw/sumprt*betat(jtime)
        wplasw(jtime)=betapw(jtime)*bp2flx/2./tmu/2./pi &
                      *vout(jtime)/1.0e+06_dp
      endif
!----------------------------------------------------------------------
!-- get normalized radial coordinate square root of toroidal flux    --
!-- at uniform poloidal flux grid sigrid                             --
!----------------------------------------------------------------------
      rhovn(nw)=1.0
      do i=2,nw-1
       rhovn(i)=sqrt(abs(rhovn(i)))
      enddo
      call zpline(nw,sigrid,rhovn,brhovn,crhovn,drhovn)
      if (kstark.gt.0.or.kdomse.gt.0) then
      do i=1,nstark
        sixxx=sigam(i)
        if (sixxx.le.1.0) then
          rhogam(i)=seval(nw,sixxx,sigrid,rhovn,brhovn,crhovn,drhovn)
        else
          rhogam(i)=1.1_dp
        endif
      enddo
      endif
!----------------------------------------------------------------------
!-- get electrostatic potential                                      --
!----------------------------------------------------------------------
      if (keecur.gt.0) then
        do i=1,nw
          epoten(i)=0.0
          call seter(sigrid(i),xsier)
          do j=1,keecur
            epoten(i)=epoten(i)+cerer(j)*xsier(j)
          enddo
        enddo
      endif
      if (idovol.gt.1) return
!-----------------------------------------------------------------------
!--  compute beta*, taking P(1)=0.0    10/25/90                       --
!-----------------------------------------------------------------------
      sumpr2=0.0
      do 1003 i=2,nw-1
        sumpr2=sumpr2+volp(i)*pprime(i)*pres(i)
 1003 continue
      sumpr2=-2.*sumpr2*dsi
      sumpr2=sumpr2/volp(nw)
      if (sumpr2.ge.0.0) then
        sumpr2=sqrt(sumpr2)
      else
        sumpr2=0.0
      endif
      betat2=sumpr2*2.0*twopi*tmu/bcentr(jtime)**2
      betat2=100.*betat2*(rout(jtime)/100./rcentr)**2
!
      do 1020 i=1,nw-1
        qpsi(i)=abs(fpol(i))/twopi*r2surf(i)
 1020 continue
      qpsi(nw)=qout(jtime)
      qpsi(1)=qmaxis
      qqmin=qpsi(1)
      iiqmin=1
      do i=2,nw
        if (qpsi(i).lt.qqmin) then
         qqmin=qpsi(i)
         iiqmin=i
        endif
      enddo
      rqqmin=sqrt(volp(iiqmin)/volp(nw))
!
      btaxp(jtime)=fpol(1)/rmaxis
      btaxv(jtime)=fpol(nw)/rmaxis
      vbeta0=pres(1)/sumpre*betat(jtime)
!
      do 1200 i=1,nw
        workrm(i)=rmaxis+(xmax-rmaxis)/(nw-1)*(i-1)
        call seva2d(bkx,lkx,bky,lky,c,workrm(i),zmaxis,pds,ier,n111)
        worksi(i)=(pds(1)-simag)/(psibry-simag)
 1200 continue
      call zpline(nw,worksi,workrm,bwork,cwork,dwork)
      do 1220 i=1,nw
        sixxx=1.0_dp/(nw-1)*(i-1)
        rpres(i)=seval(nw,sixxx,worksi,workrm,bwork,cwork,dwork)
 1220 continue
!
      do 1500 i=1,nw
        ppnow=pprime(i)
        fpnow=ffprim(i)/twopi/tmu
        cjor(i)=(ppnow+fpnow*r2surg(i))/r1surf(i)
        fpnow=ffprec(i)/twopi/tmu
        cjorec(i)=fpnow*r2surg(i)/r1surf(i)
 1500 continue
!
      DEALLOCATE(worksi,workrm,bwork,cwork,dwork,x,y,dpleng)
!
      return
 1980 format (1x,i6)
 2000 format (1x,6e12.5)
      end subroutine betali

!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          betsli computes betas and li.                           **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**                                                                  **
!**********************************************************************
      subroutine betsli(jtime,rgrid,zgrid,kerror)
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparmdud129.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      dimension pds(6),rgrid(*),zgrid(*)
      real*8,dimension(:),allocatable :: x,y,dpleng,xxs,yys
      data inorm/3/,ibtcal/2/
! MPI >>>
      integer, intent(inout) :: kerror
      kerror = 0
! MPI <<<
!
      ALLOCATE(x(nw),y(nh),dpleng(npoint),xxs(npoint),yys(npoint))
!
      if (ivacum.gt.0) return
      sumbp2=0.0
      call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
      do i=1,nw
        do j=1,nh
          kk=(i-1)*nh+j
          call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n333)
          bpolsq=(pds(2)**2+pds(3)**2)/rgrid(i)
          sumbp2=sumbp2+bpolsq*www(kk)
        enddo
      enddo
      sumbp2=sumbp2*twopi*darea
!
      plengt(1)=0.0
      do 100 i=1,nfound-1
        ip1=i+1
        dli=sqrt((xout(ip1)-xout(i))**2+(yout(ip1)-yout(i))**2)
        plengt(ip1)=plengt(i)+dli
        dpleng(i)=dli
  100 continue
!
      bp2flx=(tmu*2.0*pi*cpasma(jtime)/ &
                             plengt(nfound))**2
      const=twopi/vout(jtime)*1.0e+06_dp/bp2flx
      ali(jtime)=sumbp2/vout(jtime)/bp2flx*1.0e+06_dp
      ali3(jtime)=(tmu*2.0*pi*cpasma(jtime))**2*rout(jtime)/200.0
      ali3(jtime)=sumbp2/ali3(jtime)
!
      dsi=(psibry-simag)/(nw-1)
      volp(nw)=vout(jtime)/1.0e+06_dp
      volp(1)=0.0
      if (icurrt.ne.1) go to 540
      pprime(1)=cratio*sbeta/darea/srma
      pprime(nw)=pprime(1)
  540 continue
      if (icurrt.ne.2.and.icurrt.ne.5) go to 550
      pprime(nw)=ppcurr(x111,kppcur)/darea
      pprime(1)=ppcurr(x000,kppcur)/darea
  550 continue
      if (icurrt.ne.4) go to 600
      call currnt(n222,jtime,n222,n222,kerror)
      if (kerror.gt.0) return
      pprime(1)=cratio/darea/rzero
      pprime(nw)=pprime(1)*gammap
  600 continue
      sumpre=0.0
      nnn=1
      d11=30.
      d22=0.03_dp
      d33=0.01_dp
      nzz=0
      zzz=0.0
      do 800 i=2,nw-1
        ii=nw-i+1
        siii=(i-1)*dsi
        siwant=psibry-siii
        siii=1.0_dp-1.0_dp/(nw-1)*(i-1)
        call surfac(siwant,psi,nw,nh,rgrid,zgrid,xxs,yys,nfind, &
                    npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nnn, &
                    rmaxis,zmaxis,negcur,kerror)
        if (kerror.gt.0) return
        if (nfind.le.40.and.icntour.eq.0) then
        call cntour(rmaxis,zmaxis,siwant,xcmin,xcmax,ycmin,ycmax, &
                    yxcmin,yxcmax,xycmin,xycmax,d11,drgrid,d22, &
                    d33 ,d33 ,xmin,xmax,ymin,ymax,nzz,iautoc, &
                    xxs,yys,nfind,rgrid,nw,zgrid,nh, &
                    c,n222,nh2,nttyo,npoint, &
                    negcur,bkx,lkx,bky,lky,kerror)
          if (kerror.gt.0) return
        endif
        if (nfind.lt.10) go to 750
        volp(ii)=0.0
        xym=xxs(1)*yys(1)
        do 700 k=2,nfind
          km1=k-1
          xyp=xxs(k)*yys(k)
          dx=xxs(k)-xxs(km1)
          volp(ii)=volp(ii)+(xyp+xym)/2.0*dx
          xym=xyp
  700   continue
        volp(ii)=abs(volp(ii))*twopi
        go to 790
  750   continue
        sisi=simag-siwant
        bbsi=sqrt(sisi/siaz)
        aasi=sqrt(sisi/siar)
        volp(ii)=twopi*rmaxis*pi*aasi*bbsi
!
  790   continue
        if (icurrt.ne.2.and.icurrt.ne.5) go to 792
        pprime(ii)=ppcurr(siii,kppcur)/darea
  792   continue
        if (icurrt.ne.4) go to 794
        pprime(ii)=(1.-siii**enp)**emp*(1.-gammap)+gammap
        pprime(ii)=pprime(1)*pprime(ii)
  794   continue
        if (icurrt.ne.1) go to 796
        pprime(ii)=pprime(1)
  796   continue
        sumpre=sumpre+volp(ii)*pprime(ii)
  800 continue
      sumpre=sumpre+volp(nw)*pprime(nw)/2.
      if (kbetapr.eq.0) then
        sumpre=-sumpre*dsi/volp(nw)
      else
        if (kbetapr.eq.1) pxxxxx=prbdry
        if (kbetapr.eq.2) pxxxxx=pressb
        sumpre=-sumpre*dsi/volp(nw)+pxxxxx
      endif
      betat(jtime)=sumpre*2.0*twopi*tmu/bcentr(jtime)**2
      betap(jtime)=sumpre*2.0*twopi*tmu/bp2flx
      betat(jtime)=100.*betat(jtime)*(rout(jtime)/100./rcentr)**2
      pasman=cpasma(jtime)/1.e4_dp/aout(jtime)/abs(bcentr(jtime))
      pasman=pasman*rout(jtime)/100./rcentr
      betatn=betat(jtime)/pasman
!
      DEALLOCATE(x,y,dpleng,xxs,yys)
!
      return
      end subroutine betsli
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          chisqr computes the figure of merit for fitting         **
!**          chisq.                                                  **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      subroutine chisqr(jtime)
      include 'eparmdud129.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
!
      if (idebug>=2) then
        write (6,*) 'CHISQR, jtime= ',jtime
      endif
!
      if (ivesel.gt.10) return
      if (nbdry.le.0) go to 200
!----------------------------------------------------------------------
!-- read in the plasma response function                             --
!----------------------------------------------------------------------

      open(unit=nrsppc,status='old',form='unformatted', &
           file=table_di2(1:ltbdi2)//'ep'//trim(ch1)//trim(ch2)//'.ddd')
      read (nrsppc) gsilpc
      read (nrsppc) gmp2pc
      close(unit=nrsppc)
!
  200 continue
      nsq=2
      saisq=0.0
      chi2rm(jtime)=0.0
      do 4600 m=1,nsilop
        cm=0.0
        do 4520 n=1,nfcoil
          cm=cm+rsilfc(m,n)*brsp(n)
 4520   continue
        do 4540 n=1,nwnh
          cm=cm+gsilpc(m,n)*pcurrt(n)
 4540   continue
        if (ivesel.le.0) go to 4548
        do 4544 n=1,nvesel
          cm=cm+rsilvs(m,n)*vcurrt(n)
 4544   continue
 4548   continue
        if (iecurr.ne.1) go to 4560
        do 4550 n=1,nesum
          cm=cm+rsilec(m,n)*ecurrt(n)
 4550   continue
 4560   continue
        if (iacoil.le.0) go to 4571
        do 4569 n=1,nacoil
          cm=cm+rsilac(m,n)*caccurt(jtime,n)
 4569   continue
 4571   continue
        if (iecurr.eq.2) then
          do n=1,nesum
            cm=cm+rsilec(m,n)*cecurr(n)
          enddo
        endif
        if (fitsiref) then
            cm=cm-csiref
        endif
        if (swtsi(m).ne.0.0) then
        saisil(m)=(fwtsi(m)/swtsi(m))**nsq*(silopt(jtime,m)-cm)**2
        else
        saisil(m)=0.0
        endif
        saisq=saisq+saisil(m)
        csilop(m,jtime)=cm
        chi2rm(jtime)=max(chi2rm(jtime),saisil(m))
 4600 continue
!
      do 4700 m=1,magpri
        cm=0.0
        do 4620 n=1,nfcoil
          cm=cm+rmp2fc(m,n)*brsp(n)
 4620   continue
        do 4640 n=1,nwnh
          cm=cm+gmp2pc(m,n)*pcurrt(n)
 4640   continue
        if (ivesel.le.0) go to 4648
        do 4644 n=1,nvesel
          cm=cm+rmp2vs(m,n)*vcurrt(n)
 4644   continue
 4648   continue
        if (iecurr.ne.1) go to 4660
        do 4650 n=1,nesum
          cm=cm+rmp2ec(m,n)*ecurrt(n)
 4650   continue
 4660   continue
        if (iacoil.le.0) go to 4671
        do 4669 n=1,nacoil
          cm=cm+rmp2ac(m,n)*caccurt(jtime,n)
 4669   continue
 4671   continue
        if (iecurr.eq.2) then
          do n=1,nesum
            cm=cm+rmp2ec(m,n)*cecurr(n)
          enddo
        endif
        if (swtmp2(m).ne.0.0) then
        saimpi(m)=(fwtmp2(m)/swtmp2(m))**nsq*(expmpi(jtime,m)-cm)**2
        else
        saimpi(m)=0.0
        endif
        saisq=saisq+saimpi(m)
        cmpr2(m,jtime)=cm
        chi2rm(jtime)=max(chi2rm(jtime),saimpi(m))
 4700 continue
!-------------------------------------------------------------------
!--    calculate ECE chisqr (chiece, tchiece for psi(R+)=psi(R-))
!-------------------------------------------------------------------
       tchiece=0.0
       do 4880 m=1,nece
        cm=0.0
        do 4820 n=1,nfcoil
          cm=cm+recefc(m,n)*brsp(n)
 4820   continue
        do 4840 n=1,kwcurn
          cm=cm+recepc(m,n)*brsp(nfcoil+n)
 4840   continue
        if (ivesel.le.0) go to 4848
        do 4844 n=1,nvesel
          cm=cm+recevs(m,n)*vcurrt(n)
 4844   continue
 4848   continue
        if (iecurr.ne.1) go to 4855
        do 4850 n=1,nesum
          cm=cm+receec(m,n)*ecurrt(n)
 4850   continue
 4855   continue
       if (iacoil.le.0) go to 4871
        do 4869 n=1,nacoil
          cm=cm+receac(m,n)*caccurt(jtime,n)
 4869   continue
 4871   continue
        if (iecurr.eq.2) then
          do n=1,nesum
            cm=cm+receec(m,n)*cecurr(n)
          enddo
        endif
        if (swtece(m).ne.0.0) then
        chiece(m)=fwtece(m)**nsq*(brspece(jtime,m)-cm)**2
        chiece(m)=chiece(m)/swtece(m)**nsq
        else
        chiece(m)=0.0
        endif
        cmece(m,jtime)=cm
        tchiece=tchiece+chiece(m)
 4880 continue
!-------------------------------------------------------------------
!--    calculate ECE chisqr (chiecebz for Bz(receo)=0)
!-------------------------------------------------------------------
        cm=0.0
        do 84820 n=1,nfcoil
          cm=cm+recebzfc(n)*brsp(n)
84820   continue
        do 84840 n=1,kwcurn
          cm=cm+recebzpc(n)*brsp(nfcoil+n)
84840   continue
        if (ivesel.le.0) go to 84848
        do 84844 n=1,nvesel
          cm=cm+recevs(m,n)*vcurrt(n)
84844   continue
84848   continue
        if (iecurr.ne.1) go to 84855
        do 84850 n=1,nesum
          cm=cm+recebzec(n)*ecurrt(n)
84850   continue
84855   continue
        if (iacoil.le.0) go to 84871
        do 84869 n=1,nacoil
          cm=cm+receac(m,n)*caccurt(jtime,n)
84869   continue
84871   continue
        if (iecurr.eq.2) then
          do n=1,nesum
            cm=cm+recebzec(n)*cecurr(n)
          enddo
        endif
        if (swtecebz.ne.0.0) then
        chiecebz=fwtecebz**nsq*(brspecebz(jtime)-cm)**2
        chiecebz=chiecebz/swtecebz**nsq
        else
        chiecebz=0.0
        endif
        cmecebz(jtime)=cm
!
      cm=0.0
      do 4900 n=1,nwnh
        cm=cm+pcurrt(n)
 4900 continue
      cpasma(jtime)=cm
      if (ifitvs.gt.0.or.ivesel.gt.0) then
        do 4910 i=1,nvesel
          cm=cm+vcurrt(i)
 4910   continue
      endif
      if (swtcur.ne.0.0) then
      saiip=(fwtcur/swtcur)**nsq*(pasmat(jtime)-cm)**2
      else
      saiip=0.0
      endif
      saisq=saisq+saiip
      chi2rm(jtime)=max(chi2rm(jtime),saiip)
!
      tsaifc=0.0
      do 4920 i=1,nfcoil
        saifc(i)=0.0
        if (fwtfc(i).gt.0.0) then
          saifc(i)=fwtfc(i)**nsq*(brsp(i)-fccurt(jtime,i))**2
          saifc(i)=saifc(i)/swtfc(i)**nsq
        endif
        saisq=saisq+saifc(i)
        tsaifc=tsaifc+saifc(i)
        chi2rm(jtime)=max(chi2rm(jtime),saifc(i))
 4920 continue
      if (iecurr.eq.2) then
      do i=1,nesum
        saiec(i)=0.0
        if (fwtec(i).gt.0.0) then
          saiec(i)=fwtec(i)**nsq*(cecurr(i)-ecurrt(i))**2
          saiec(i)=saiec(i)/swtec(i)**nsq
        endif
        saisq=saisq+saiec(i)
        chi2rm(jtime)=max(chi2rm(jtime),saiec(i))
      enddo
      endif
      if (fitsiref) then
        saisref=0.0
        if (fwtref.gt.0.0) then
          saisref=fwtref**nsq*(psiref(jtime)-csiref)**2
          saisref=saisref/swtsi(nslref)**nsq
        endif
        saisq=saisq+saisref
        chi2rm(jtime)=max(chi2rm(jtime),saisref)
      endif
!
      do 4950 n=1,nfcoil
        ccbrsp(n,jtime)=brsp(n)
 4950 continue
!
      tsaisq(jtime)=saisq
      if (iand(iout,1).ne.0) then
      write (nout,7400) time(jtime),tsaisq(jtime),cpasma(jtime)
      write (nout,7420)
      write (nout,7450) (saisil(m),m=1,nsilop)
      write (nout,7430)
      write (nout,7450) (saimpi(m),m=1,magpri)
      write (nout,7460) saiip
      write (nout,7480) tsaifc
      write (nout,7450) (saifc(m),m=1,nfcoil)
      write (nout,7482) saisref
      write (nout,7485)
      write (nout,7450) (saiec(m),m=1,nesum)
      if (kprfit.gt.0) write (nout,7470) chipre
      if (kprfit.gt.0) write (nout,7450) (saipre(m),m=1,npress)
      if (kecebz.gt.0) write (nout,7486) chiecebz
      if (kece.gt.0) write (nout,7487) tchiece
      if (kece.gt.0) then
       write(nout,7488)
       write (nout,7450) (chiece(m),m=1,nece)
      endif
      endif
!-------------------------------------------------------------------------
!-- compute signals at MSE locations if requested                       --
!-------------------------------------------------------------------------
      if (kdomse.gt.0) then
      nnn=0
      call green(nnn,jtime,n222)
      do 4800 m=1,nstark
        cmbr=0.0
        cmbz=0.0
        do 4720 n=1,nfcoil
          cmbr=cmbr+rbrfc(m,n)*brsp(n)
          cmbz=cmbz+rbzfc(m,n)*brsp(n)
 4720   continue
        do 4740 n=1,kwcurn
          cmbr=cmbr+rbrpc(m,n)*brsp(nfcoil+n)
          cmbz=cmbz+rbzpc(m,n)*brsp(nfcoil+n)
 4740   continue
        if (kedgep.gt.0) then
          cmbr=cmbr+rbrpe(m)*pedge
          cmbz=cmbz+rbzpe(m)*pedge
        endif
        if (kedgef.gt.0) then
          cmbr=cmbr+rbrfe(m)*f2edge
          cmbz=cmbz+rbzfe(m)*f2edge
        endif
        if (ivesel.le.0) go to 4748
        do 4744 n=1,nvesel
          cmbr=cmbr+rbrvs(m,n)*vcurrt(n)
          cmbz=cmbz+rbzvs(m,n)*vcurrt(n)
 4744   continue
 4748   continue
        if (iecurr.ne.1) go to 4755
        do 4750 n=1,nesum
          cmbr=cmbr+rbrec(m,n)*ecurrt(n)
          cmbz=cmbz+rbzec(m,n)*ecurrt(n)
 4750   continue
 4755   continue
        if (iecurr.eq.2) then
          do n=1,nesum
            cmbr=cmbr+rbrec(m,n)*cecurr(n)
            cmbz=cmbz+rbzec(m,n)*cecurr(n)
          enddo
        endif
        if (iacoil.le.0) go to 4761
        do 4759 n=1,nacoil
          cmbr=cmbr+rbrac(m,n)*caccurt(jtime,n)
          cmbz=cmbz+rbzac(m,n)*caccurt(jtime,n)
 4759   continue
 4761   continue
        cm=a2gam(jtime,m)*btgam(m)+a3gam(jtime,m)*cmbr+a4gam(jtime,m) &
             *cmbz
        bzmsec(m)=cmbz
        if (keecur.le.0) then
        cm=a1gam(jtime,m)*cmbz/cm
        else
        ce1rbz=0.0
        ce2rbz=0.0
        ce3rbr=0.0
        do n=1,keecur
         ce1rbz=ce1rbz+e1rbz(m,n)*cerer(n)
         ce2rbz=ce2rbz+e2rbz(m,n)*cerer(n)
         ce3rbr=ce3rbr+e3rbr(m,n)*cerer(n)
        enddo
        cm=cm-ce2rbz-ce3rbr
        cm=(a1gam(jtime,m)*cmbz-ce1rbz)/cm
        endif
        cmgam(m,jtime)=cm
 4800 continue
      endif
!-------------------------------------------------------------------------
!-- compute signals at MSE-LS locations if requested                    --
!-------------------------------------------------------------------------
      if (kdomsels.gt.0) then
      call green(nnn,jtime,n222)
      do 94800 m=1,nmsels
        cmbr=0.0
        do 94720 n=1,kpcurn
          cmbr=cmbr+rmlspc(m,n)*brsp(nfcoil+n)
94720   continue
        cmbr=cmbr-rhsmls(jtime,m)
        cmmls(jtime,m)=sqrt(cmbr)
        cmmls2(jtime,m)=l1mselt(jtime,m)*btmls(m)**2+l2mselt(jtime,m)* &
           brmls(m)*btmls(m)+l3mselt(jtime,m)*brmls(m)**2+             &
          (l1mselt(jtime,m)+l3mselt(jtime,m))*bzmls(m)**2
        cmmls2(jtime,m)=sqrt(cmmls2(jtime,m))
94800 continue
      write (nout,7585)
      write (nout,7450) (cmmls(jtime,m),m=1,nmsels)
      write (nout,7588)
      write (nout,7450) (cmmls2(jtime,m),m=1,nmsels)
      endif
!
      if (mmbmsels.gt.0.or.kdomsels.gt.0) then
        do m=1,nmsels
          cmmlsv(jtime,m)=0.0
          if (rrmselt(jtime,m).gt.0.0) then
            cmmlsv(jtime,m)=abs(bcentr(jtime)*rcentr/rrmselt(jtime,m))
            cmmlsv(jtime,m)=cmmlsv(jtime,m)*sqrt(l1mselt(jtime,m))
          endif
        enddo 
        write (nout,7590)
        write (nout,7450) (cmmlsv(jtime,m),m=1,nmsels)
      endif
!
      return
 7400 format (/,2x,'time = ',e12.5,2x,'chisq = ',e12.5,2x,'current = ',e12.5)
 7420 format (10x,'chi psi loops:')
 7430 format (10x,'chi inner magnetic probes:')
 7450 format (8(1x,e12.5:,1x))
 7460 format (10x,'chi ip:',/,15x,e12.5)
 7470 format (10x,'chi pressure:         ',/,1x,e12.5)
 7480 format (10x,'chi F-coils:          ',/,10x,e12.5)
 7482 format (10x,'chi psiref:',/,15x,e12.5)
 7485 format (10x,'chi E-coils:          ',/,1x,e12.5)
 7486 format (10x,'chi ecebz:            ',/,1x,e12.5)
 7487 format (10x,'chi total eceR+R-:    ',/,1x,e12.5)
 7488 format (10x,'chi eceR+R-:          ')
 7585 format (10x,'Simulated MSE-LS (T): ',/)
 7588 format (10x,'Simulated MSE-LS2 (T): ',/)
 7590 format (10x,'Simulated MSE-LSV (T): ',/)
      end subroutine chisqr

!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          erpote computes the stream function for the             **
!**          radial electric field. eradial computes the             **
!**          radial electric field.                                  **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          97/04/24..........first created                         **
!**                                                                  **
!**********************************************************************
      function erpote(ypsi,nnn)
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparmdud129.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      dimension xpsii(nercur)
      data init/0/
!
      if (abs(ypsi).gt.1.0) then
        erpote=0.0
        return
      endif
      erpote=0.0
      call seter(ypsi,xpsii)
      do 1400 iiij=1,nnn
        erpote=erpote+cerer(iiij)*xpsii(iiij)
 1400 continue
      return
!
      entry erppote(ypsi,nnn)
      if (abs(ypsi).gt.1.0) then
        erppote=0.0
        return
      endif
      erppote=0.0
      call seterp(ypsi,xpsii)
      do 2400 iiij=1,nnn
        erppote=erppote+cerer(iiij)*xpsii(iiij)
 2400 continue
      erppote=-erppote/sidif
      return
      end function erpote

!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          eradial computes the radial electric field.             **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          97/04/24..........first created                         **
!**                                                                  **
!**********************************************************************
      function eradial(ypsi,nnn,reee,zeee)
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparmdud129.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      dimension pds(6)
      data init/0/
!
      if (abs(ypsi).ge.1.0) then
        eradial=0.0
        return
      endif
        call seva2d(bkx,lkx,bky,lky,c,reee,zeee,pds,ier,n333)
      eradial=erpote(ypsi,nnn)
      eradial=-pds(2)*eradial
      return
!
      entry esradial(ypsi,nnn,reee,zeee)
      if (abs(ypsi).ge.1.0) then
        esradial=0.0
        return
      endif
      call seva2d(bkx,lkx,bky,lky,c,reee,zeee,pds,ier,n333)
      esradial=erppote(ypsi,nnn)
      esradial=-(reee*pds(2))**2*esradial
      fnow=seval(nw,ypsi,sigrid,fpol,bbfpol,ccfpol,ddfpol)
      bbnow=sqrt(fnow**2+pds(2)**2+pds(3)**2)/reee
      esradial=esradial/bbnow
      return
      end function

!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          autoknot minimizes chi-squared as a function of knot    **
!**          location                                                **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          20/01/2000........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine autoknot(ks,lconvr,ktime,mtear,kerror)
      include 'eparmdud129.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      external ppakfunc,ffakfunc,wwakfunc,eeakfunc
      integer, intent(inout) :: kerror
      kerror = 0

!
!    Store the values away so the functions can restore them 
!    after reading the k file. 
!
      ks_a = ks
      lconvr_a = lconvr
      ktime_a = ktime
      mtear_a = mtear
      kerror_a = kerror_a
      tol = aktol
      saveiter = mxiter
      call store_autoknotvals
      mxiter_a = kakiter
   
!
!     Minimize chi^2 and error for pp knot location
!
      if(kppfnc .eq. 6)then
        kakpploop = kakloop
        if(kppknt .le. 3) kakpploop = 1
        do j= 1, kakpploop
          do i= 2,kppknt-1
            kappknt = kppknt
            kadknt = i
            dzero = appknt(i-1)
            done = appknt(i+1)
            prevknt = appknt(kadknt)
            appknt(kadknt) = fmin(dzero,done,ppakfunc,tol)
            write(6,*) 'pp knot ',kadknt,' set to ',appknt(kadknt)
          enddo
          if(abs(prevknt - appknt(kadknt)) .le. aktol)then
            write(6,*)'Last PPknot changed by less that use tolerance'
            go to 10
          endif
        enddo
      endif

!
!     Minimize chi^2 and error for ff knot location
!
10    continue
      if(kfffnc .eq. 6)then
      kakffloop = kakloop
      if(kffknt .le. 3) kakffloop = 1
      do j= 1, kakffloop
      do i= 2,kffknt-1
         kaffknt = kffknt
         kadknt = i
         dzero = affknt(i-1)
         done = affknt(i+1)
         prevknt = affknt(kadknt)
         affknt(kadknt) = fmin(dzero,done,ffakfunc,tol)
         write(6,*) 'ff knot ',kadknt,' set to ',affknt(kadknt)
      enddo
         if(abs(prevknt - affknt(kadknt)) .le. aktol)then
           write(6,*)'Last FFknot changed by less that use tolerance'
           go to 20
         endif
      enddo
      endif

!
!     Minimize chi^2 and error for ww knot location
!

20    continue
      if(kwwfnc .eq. 6)then
         kakwwloop = kakloop
         if(kwwknt .le. 3) kakwwloop = 1
         do j= 1, kakwwloop
         do i= 2,kwwknt-1
            kawwknt = kwwknt
            kadknt = i
            dzero = awwknt(i-1)
            done = awwknt(i+1)
            prevknt = awwknt(kadknt)
            awwknt(kadknt) = fmin(dzero,done,wwakfunc,tol)
            write(6,*) 'ww knot ',kadknt,' set to ',awwknt(kadknt)
         enddo
            if(abs(prevknt - awwknt(kadknt)) .le. aktol)then
              write(6,*)'Last WWknot changed by less that use tolerance'
              go to 30
            endif
         enddo
      endif


!
!     Minimize chi^2 and error for ee knot location
!

30    continue
      if(keefnc .eq. 6)then
         kakeeloop = kakloop
         if(keeknt .le. 3) kakeeloop = 1
         do j= 1, kakeeloop
         do i= 2,keeknt-1
            kaeeknt = keeknt
            kadknt = i
            dzero = aeeknt(i-1)
            done = aeeknt(i+1)
            prevknt = aeeknt(kadknt)
            aeeknt(kadknt) = fmin(dzero,done,eeakfunc,tol)
            write(6,*) 'ee knot ',kadknt,' set to ',aeeknt(kadknt)
         enddo
            if(abs(prevknt - aeeknt(kadknt)) .le. aktol)then
              write(6,*)'Last EEknot changed by less that use tolerance'
              go to 40
            endif
         enddo
      endif
!
!     Now do the final fit with adjusted knots and full iterations
!

40    continue
      call data_input(ks_a,lconvr_a,ktime_a,mtear_a,kerror)
      if (kerror.gt.0) return
      if (lconvr_a.lt.0) return
      mxiter_a = saveiter
      call restore_autoknotvals
      call inicur(ks_a)
      call fit(ks_a,kerror_a)
      if (kerror_a /= 0) then
        kerror = 1
        return
      endif

      return
      end subroutine autoknot
!
!    store values read from k file into autoknot variables
!
      subroutine restore_autoknotvals
      include 'eparmdud129.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      do i = 1,npcurn
          ppknt(i) = appknt(i)
          ffknt(i) = affknt(i)
          wwknt(i) = awwknt(i)
          eeknt(i) = aeeknt(i)
      enddo
      kppknt = kappknt
      kffknt = kaffknt
      kwwknt = kawwknt
      keeknt = kaeeknt
      mxiter = mxiter_a
      return
      end subroutine restore_autoknotvals
!
!     store autoknot variables into standard efit names
!     for example knot locations
!
      subroutine store_autoknotvals
      include 'eparmdud129.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      do i = 1,npcurn
          appknt(i) = ppknt(i)
          affknt(i) = ffknt(i)
          awwknt(i) = wwknt(i)
          aeeknt(i) = eeknt(i)
      enddo
      kappknt = kppknt
      kaffknt = kffknt
      kawwknt = kwwknt
      kaeeknt = keeknt
      mxiter_a = mxiter
      return
      end subroutine store_autoknotvals
!
! used by autoknot which passes the routine to a minimization subroutine
! which calls it to evaulate the function being minimized
!
      function ppakfunc(xknot) ! TODO: kerror is not returned
      include 'eparmdud129.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      kerror = 0
      ppakfunc = 1000.0
      write(6,*)
      write(6,*)' trying pp knot ',kadknt,' at location ',xknot, &
        ' out of ',kappknt,' knots'
      call data_input(ks_a,lconvr_a,ktime_a,mtear_a,kerror)
      if (kerror.gt.0) return
      if (lconvr_a.lt.0) return
      call restore_autoknotvals
      ppknt(kadknt) = xknot
      call inicur(ks_a)
      call fit(ks_a,kerror_a)
      if(kerror_a .gt. 0) return
      ppakfunc = akchiwt * tsaisq(ks_a)  + akerrwt * errorm  &
        + akgamwt * chigamt + akprewt * chipre
      return
      end function ppakfunc
!
! used by autoknot which passes the routine to a minimization subroutine
! which calls it to evaulate the function being minimized
!
      function ffakfunc(xknot) ! TODO: kerror is not returned
      include 'eparmdud129.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      kerror = 0
      ffakfunc = 1000.0
      write(6,*)
      write(6,*)' trying ff knot ',kadknt,' at location ',xknot, &
        ' out of ',kaffknt,' knots'
      call data_input(ks_a,lconvr_a,ktime_a,mtear_a,kerror)
      if (kerror.gt.0) return
      if (lconvr_a.lt.0) return
      call restore_autoknotvals
      ffknt(kadknt) = xknot
      call inicur(ks_a)
      call fit(ks_a,kerror_a)
      if(kerror_a .gt. 0) return
      ffakfunc = akchiwt * tsaisq(ks_a)  + akerrwt * errorm  &
        + akgamwt * chigamt + akprewt * chipre
      return
      end function ffakfunc
!
! used by autoknot which passes the routine to a minimization subroutine
! which calls it to evaulate the function being minimized
!
      function wwakfunc(xknot) ! TODO: kerror is not returned
      include 'eparmdud129.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      kerror = 0
      wwakfunc = 1000.0
      write(6,*)
      write(6,*)' trying ww knot ',kadknt,' at location ',xknot, &
                 ' out of ',kawwknt,' knots'
      call data_input(ks_a,lconvr_a,ktime_a,mtear_a,kerror)
      if (kerror.gt.0) return
      if (lconvr_a.lt.0) return
      call restore_autoknotvals
      wwknt(kadknt) = xknot
      call inicur(ks_a)
      call fit(ks_a,kerror_a)
      if(kerror_a .gt. 0) return
      wwakfunc = akchiwt * tsaisq(ks_a)  + akerrwt * errorm &
        + akgamwt * chigamt + akprewt * chipre
      return
      end function wwakfunc
!
! used by autoknot which passes the routine to a minimization subroutine
! which calls it to evaulate the function being minimized
!
      function eeakfunc(xknot) ! TODO: kerror is not returned
      include 'eparmdud129.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      kerror = 0
      eeakfunc = 1000.0
      write(6,*)
      write(6,*)' trying ee knot ',kadknt,' at location ',xknot, &
                 ' out of ',kaeeknt,' knots'
      call data_input(ks_a,lconvr_a,ktime_a,mtear_a,kerror)
      if (kerror.gt.0) return
      if (lconvr_a.lt.0) return
      call restore_autoknotvals
      eeknt(kadknt) = xknot
      call inicur(ks_a)
      call fit(ks_a,kerror_a)
      if(kerror_a .gt. 0) return
      eeakfunc = akchiwt * tsaisq(ks_a)  + akerrwt * errorm &
        + akgamwt * chigamt + akprewt * chipre
      return
      end function

!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          fit carries out the fitting and equilibrium iterations. **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**          94/03/11..........revised                               **
!**                                                                  **
!**********************************************************************
      subroutine fit(jtime,kerror)
      include 'eparmdud129.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      data nzero/0/
      save nzero
!-----------------------------------------------------------------------
!--  inner equilibrium do loop and outer current profile parameters   --
!--  do loop                                                          --
!-----------------------------------------------------------------------
      if (idebug /= 0) write (6,*) 'Enter FIT'
      kerror=0
      jerror(jtime)=0
      nitera=0
      ndelzon=999
      iwantk=0
      iend=mxiter+1
      if (iconvr.eq.3) iend=1
      idosigma=1
      do i=1,iend
        ix=i
        if (i.gt.1) then
          iwantk=iwantk+1
          !------------------------------------------------------------------
          !--  nitera.ge.kcallece, then call setece                        --
          !--  mtxece=1 call setece every iteration                        --
          !--  mtxece>1 call setece per mtece time iteration               --
          !------------------------------------------------------------------
          if ((nitera.ge.kcallece).and.(kfitece.gt.0)) then
            nleft=abs(mxiter)-nitera
            if(nleft .ge. mtxece*nconstr) then
              if (idebug.ge.2) then
                write (6,*) 'Call FIT/SETECE',ix
                write (6,*) '  nitera/kcallece/kfitece/nleft/mtxece/nconstr = ', &
                  nitera,kcallece,kfitece,nleft,mtxece,nconstr
              endif
              call setece(jtime,kerror)
              if (kerror /= 0) then
                jerror(jtime) = 1
                return
              endif
            endif
          endif
          call green(nzero,jtime,nitera)
          if (kprfit.gt.0.and.iwantk.eq.ndokin) then
            call presur(jtime,nitera,kerror)
            if (kerror /= 0) then
              jerror(jtime) = 1
              return
            endif
            iwantk=0
          endif
          if (kprfit.ge.3) call presurw(jtime,nitera)
          if (errorm.lt.errmagb) then
            if ((imagsigma.eq.1) .AND. (errorm > errmag) ) &
              call getsigma(jtime,nitera)
            if ((imagsigma.eq.2).and.(idosigma.eq.1)) then
              call getsigma(jtime,nitera)
              idosigma=2
            endif
          endif
          if (idebug.ge.2) write (6,*) 'Call FIT/MATRIX',ix
          call matrix(jtime,ix,ichisq,nitera,kerror)
           
          if (kerror /= 0) then
            jerror(jtime) = 1
            return
          endif
          if ((iconvr.eq.2).and.(ichisq.gt.0)) then
            call errctrl_msg('fit','not converged properly',2)
            go to 2020
          end if
        end if

        do in=1,nxiter
          ixnn=in
          nitera=nitera+1

          if (idebug.ge.2) write(6,*) 'Entering currnt' 
          call currnt(ix,jtime,ixnn,nitera,kerror)
          if (kerror /= 0) then
            jerror(jtime) = 1
            return
          endif
           
          if (ivesel.ge.2) then
             if (idebug.ge.2) write(6,*) 'Entering vescur'
             call vescur(jtime)
          endif
          if ((i.le.1).or.(in.gt.1)) then
              if (idebug.ge.2) write(6,*) 'Entering fcurrt'
             call fcurrt(jtime,ix,nitera,kerror)
          endif 
          if (kerror /= 0) then
            jerror(jtime) = 1
            return
          endif
          if (idebug.ge.2) write(6,*) 'Entering pflux' 
          call pflux(ix,ixnn,nitera,jtime,kerror)
          if (kerror.gt.0) then
            jerror(jtime) = 1
            return
          endif

          if (idebug.ge.2) write(6,*) 'Entering steps'
          call steps(ixnn,nitera,ix,jtime,kerror)
          if (kerror.gt.0) then
            jerror(jtime) = 1
            return
          endif
          if (kmtark.gt.0) then
            if (kwaitmse.ne.0 .and. i.ge.kwaitmse) then
              if (idebug.ge.2) write(6,*) 'Entering fixstark'
              call fixstark(jtime,kerror)
              if (kerror.gt.0) then
                jerror(jtime) = 1
                return
              endif
            end if
          endif

          if (idebug.ge.2) write(6,*) 'Entering residu'
          call residu(nitera,jtime)
          if ((nitera.lt.kcallece).and.(kfitece.gt.0.0)) exit
          if ((in.eq.1).and.(idone.gt.0).and.(tsaisq(jtime).le.saimin)) then
            go to 2020
          end if
          if (idone.gt.0) exit
          if (i.eq.mxiter+1) exit
        end do ! in
      end do ! i
      call errctrl_msg('fit','not converged, reached max iterations',2)
2020  continue
      !---------------------------------------------------------------------
      !--  update pressure if needed                                      --
      !---------------------------------------------------------------------
      if (kprfit.gt.1) then

        if (idebug.ge.2) write(6,*) 'Entering presur'
        call presur(jtime,nitera,kerror)
        if (kerror /= 0) then
          jerror(jtime) = 1
          return
        endif
      end if
      return
      end subroutine fit

!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          fpcurr computes the radial derivative                   **
!**          of the poloidal current ff. ffcurr computes             **
!**          the poloidal current F=twopi RBt/mu0                    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**          94/03/08..........revised                               **
!**                                                                  **
!**********************************************************************
      function fpcurr(upsi,nnn)
      include 'eparmdud129.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      dimension xpsii(nffcur)

      if (icutfp.gt.0) then
        ypsi=upsi*xpsimin
      else
        ypsi=upsi
      endif
      if (abs(ypsi).gt.1.0) then
        fpcurr=0.0
        return
      endif
      if (npsi_ext > 0) then
        fpcurr = seval(npsi_ext,ypsi,psin_ext,ffprim_ext,bfp_ext,cfp_ext,dfp_ext)
        fpcurr = fpcurr * cratiof_ext
        return
      endif
      fpcurr=0.0
      call setfp(ypsi,xpsii)
      do 1400 iiij=nbase+1,nbase+nnn
        iijj=iiij-nbase
        fpcurr=fpcurr+brsp(iiij)*xpsii(iijj)
 1400 continue
!----------------------------------------------------------------------
!-- edge hyperbolic tangent component                                --
!----------------------------------------------------------------------
      if (kedgef.eq.0) return
      siedge=(ypsi-fe_psin)/fe_width
      f0edge=f2edge
      f0back=f0edge/fe_width/sidif
      fpcurr=fpcurr+f0back/cosh(siedge)**2
      return
!
      entry fpecrr(upsi,nnn)
      if (icutfp.gt.0) then
        ypsi=upsi*xpsimin
      else
        ypsi=upsi
      endif
      if (abs(ypsi).gt.1.0) then
        fpecrr=0.0
        return
      endif
      fpecrr=0.0
      call setfp(ypsi,xpsii)
      do 1500 iiij=nbase+nnn,nbase+nnn
        iijj=iiij-nbase
        fpecrr=fpecrr+brsp(iiij)*xpsii(iijj)
 1500 continue
      return
!
      entry ffcurr(upsi,nnn)
      if (icutfp.gt.0) then
        ypsi=upsi*xpsimin
        xsidif=-sidif/xpsimin
      else
        ypsi=upsi
        xsidif=-sidif
      endif
      if (abs(ypsi).ge.1.0) then
        ffcurr=fbrdy
        return
      endif
      ffcurr=0.0
      call setff(ypsi,xpsii)
      do 1600 i=nbase+1,nbase+nnn
        nn=i-nbase
        ffcurr=ffcurr+brsp(i)*xpsii(nn)
 1600 continue
      fb22=fbrdy**2
      ff22=fb22+xsidif*ffcurr/constf2
      if (ff22.lt.0.0) ff22=fb22
      ffcurr=sqrt(ff22)*fbrdy/abs(fbrdy)
!----------------------------------------------------------------------
!-- edge hyperbolic tangent component                                --
!----------------------------------------------------------------------
      if (kedgef.eq.0) return
      siedge=(ypsi-fe_psin)/fe_width
      f0edge=f2edge/constf2
      ff22=ff22+f0edge*(tfedge-tanh(siedge))
      if (ffcurr.lt.0.0) ffcurr=fb22
      ffcurr=sqrt(ff22)*fbrdy/abs(fbrdy)
      return
      end function fpcurr
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          FLUXAV does the flux surface average.                   **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**          f     : function array                                  **
!**          si    : flux array                                      **
!**          x,y,n : n (R,Z) coordinates of contour                  **
!**                                                                  **
!**     RETURN ARGUMENTS:                                            **
!**          fave  : int(dl f/Bp) / sdlobp                           **
!**          sdlbp : int(dl Bp)                                      **
!**          sdlobp: int(dl/Bp)                                      **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          14/10/87..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine fluxav(f,x,y,n,si,rx,msx,ry,msy,fave,ns,sdlobp,sdlbp)
      use set_kinds
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparmdud129.inc'
      use var_cwork3, only:lkx,lky
!vas
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      dimension f(1),x(1),y(1),si(1),pds(6),rx(1),ry(1)
!
      if (ns.eq.0) go to 2000
      call sets2d(si,c,rx,msx,bkx,lkx,ry,msy,bky,lky,wk,ier)
!
 2000 continue
!------------------------------------------------------------------
!--  flux surface average of f                                   --
!------------------------------------------------------------------
      fave=0.0
      fnorm=0.0
      sdlbp=0.0
      do 2050 i=2,n
        xnow=0.5_dp*(x(i-1)+x(i))
        ynow=0.5_dp*(y(i-1)+y(i))
        fnow=0.5_dp*(f(i-1)+f(i))
        dxnow=x(i)-x(i-1)
        dynow=y(i)-y(i-1)
        dl=sqrt(dxnow**2+dynow**2)
        call seva2d(bkx,lkx,bky,lky,c,xnow,ynow,pds,ier,n333)
        bpol = sqrt(pds(2)**2+pds(3)**2)/xnow
        dlbpol = dl/bpol
        fnorm = fnorm + dlbpol
        fave = fave + dlbpol*fnow
        sdlbp = sdlbp + dl*bpol
 2050 continue
      fave = fave/fnorm
      sdlobp = fnorm
      return
      end subroutine fluxav
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          green set up the appropriate response functions for use **
!**          with the routine matrix.                                **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**          94/03/08..........revised                               **
!**        2015/03/27..........revised MSE-LS                        **  
!**                                                                  **
!**********************************************************************
      subroutine green(ifag,jtime,niter)
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparmdud129.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      dimension xpsii(nwcurn),xpsis(nwcurn),xpsisb(nwcurn)
      dimension xsier(nercur)
      dimension pds(6),xnsi(nppcur)
!
 1100 continue
      kffp1=kffcur+1
      do 1300 jj=1,kwcurn
        if (ifag.eq.1) go to 1280
        do 1250 m=1,nsilop
          rsilpc(m,jj)=0.0
 1250   continue
        do 1260 m=1,magpri
          rmp2pc(m,jj)=0.0
 1260   continue
        if (kstark.gt.0.or.kdomse.gt.0) then
        do 1270 m=1,nstark
          rbrpc(m,jj)=0.0
          rbzpc(m,jj)=0.0
 1270   continue
        endif
        if (mmbmsels.gt.0.or.kdomsels.gt.0) then
          do m=1,nmsels
            rmlspc(m,jj)=0.0
          enddo
        endif
        if (kece.gt.0) then
        do 1290 m=1,nnece
          recepc(m,jj)=0.0
 1290   continue
        endif
        if (kecebz.gt.0) then
        recebzpc(jj)=0.0
        endif
        if (nbdry.gt.0) then
           do m=1,nbdry
             gbdrpc(m,jj)=0.0
           enddo
        endif
 1280   continue
        fgowpc(jj)=0.0
 1300 continue
      if (fitdelz.and.niter.ge.ndelzon) then
          do m=1,nsilop
            gsildz(m)=0.0
          enddo
          do m=1,magpri
            gmp2dz(m)=0.0
          enddo
          fgowdz=0.0
          if (kstark.gt.0) then
            do m=1,nstark
              gbrdz(m)=0.0
              gbzdz(m)=0.0
            enddo
          endif
          if (kece.gt.0) then
            do m=1,nnece
              gecedz(m)=0.0
            enddo
          endif
          if (kecebz.gt.0) then
              gecebzdz=0.0
          endif
          if (nbdry.gt.0) then
             do m=1,nbdry
               gbdrdz(m)=0.0
             enddo
          endif
      endif
      if (ifag.ne.1) then
       do jj=1,kffcur
         rspdlc(jj)=0.0
       enddo
      endif
!
      if (fwtdlc.gt.0.0) then
        upsi1=1.0
        call setff(upsi1,xpsisb)
      endif
!------------------------------------------------------------------
!-- Hyperbolic tangent term                                      --
!------------------------------------------------------------------
      if (kedgep.gt.0) then
       if (ifag.ne.1) then
         do m=1,nsilop
           rsilpe(m)=0.0
         enddo
         do m=1,magpri
           rmp2pe(m)=0.0
         enddo
         if (kstark.gt.0.or.kdomse.gt.0) then
           do m=1,nstark
             rbrpe(m)=0.0
             rbzpe(m)=0.0
           enddo
         endif
         if (nbdry.gt.0) then
           do m=1,nbdry
             gbdrpe(m)=0.0
           enddo
         endif
       endif
         fgowpe=0.0
      endif
!
      if (kedgef.gt.0) then
       if (ifag.ne.1) then
         do m=1,nsilop
           rsilfe(m)=0.0
         enddo
         do m=1,magpri
           rmp2fe(m)=0.0
         enddo
         if (kstark.gt.0.or.kdomse.gt.0) then
           do m=1,nstark
             rbrfe(m)=0.0
             rbzfe(m)=0.0
           enddo
         endif
         if (nbdry.gt.0) then
           do m=1,nbdry
             gbdrfe(m)=0.0
           enddo
         endif
         rdlcfe=0.0
       endif
         fgowfe=0.0
      endif
!
      do 2000 i=1,nw
      do 2000 j=1,nh
        kk=(i-1)*nh+j
        if ((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.1.0)) then
          go to 1582
        end if
        call setpp(xpsi(kk),xpsii)
        do 1580 jj=1,kppcur
          factor=xpsii(jj)*rgrid(i)*www(kk)
!------------------------------------------------------------------------
!-- correction for rotation, p0' terms                                 --
!------------------------------------------------------------------------
          if (niter.gt.1) then
          if (kvtor.eq.2) then
            prew0=pwcurr(xpsi(kk),kwwcur)
            pres0=prcurr(xpsi(kk),kppcur)
            if (abs(pres0).gt.1.e-10_dp) then
               pwop0=prew0/pres0
               pwp0r2=pwop0*rgrvt(i)
            else
               pwop0=0.0
               pwp0r2=0.0
            endif
            factor=factor*(1.-0.5_dp*pwp0r2**2)
          elseif (kvtor.eq.3) then
            prew0=pwcurr(xpsi(kk),kwwcur)
            pres0=prcurr(xpsi(kk),kppcur)
            if (abs(pres0).gt.1.e-10_dp) then
               pwop0=prew0/pres0
               pwp0r2=pwop0*rgrvt(i)
               ptop0=exp(pwp0r2)
            else
               pwop0=0.0
               pwp0r2=0.0
               ptop0=1.0
            endif
            factor=factor*ptop0*(1.-pwp0r2)
          endif
          endif
          if (ifag.eq.1) go to 1520
          do 1350 m=1,nsilop
            rsilpc(m,jj)=rsilpc(m,jj) + gsilpc(m,kk)*factor
 1350     continue
          do 1400 m=1,magpri
            rmp2pc(m,jj)=rmp2pc(m,jj) + gmp2pc(m,kk)*factor
 1400     continue
          if (kstark.gt.0.or.kdomse.gt.0) then
          do 1500 m=1,nstark
            rbrpc(m,jj)=rbrpc(m,jj) + gbrpc(m,kk)*factor
            rbzpc(m,jj)=rbzpc(m,jj) + gbzpc(m,kk)*factor
 1500     continue
          endif
          if (kece.gt.0) then
          do 1510 m=1,nnece
            recepc(m,jj)=recepc(m,jj) + gecepc(m,kk)*factor
 1510     continue
          endif
          if (kecebz.gt.0) then
            recebzpc(jj)=recebzpc(jj) + gecebzpc(kk)*factor
          endif
!---------------------------------------------------------------------------
!-- Boundary constraints                                                  --
!---------------------------------------------------------------------------
          if (nbdry.gt.0) then
           do m=1,nbdry
             gbdrpc(m,jj)=gbdrpc(m,jj) + rbdrpc(m,kk)*factor
           enddo
          endif
 1520     continue
          fgowpc(jj)=fgowpc(jj) + factor
 1580   continue
!-------------------------------------------------------------------------
!-- Hyperbolic tangent term for P'                                      --
!-------------------------------------------------------------------------
        if (kedgep.gt.0) then
          siedge=(xpsi(kk)-pe_psin)/pe_width
          xpsnow=1./pe_width/sidif/cosh(siedge)**2
          factor=xpsnow*rgrid(i)*www(kk)
          if (ifag.eq.1) go to 19580
          do m=1,nsilop
            rsilpe(m)=rsilpe(m) + gsilpc(m,kk)*factor
          enddo
          do m=1,magpri
            rmp2pe(m)=rmp2pe(m) + gmp2pc(m,kk)*factor
          enddo
          if (kstark.gt.0.or.kdomse.gt.0) then
          do m=1,nstark
            rbrpe(m)=rbrpe(m) + gbrpc(m,kk)*factor
            rbzpe(m)=rbzpe(m) + gbzpc(m,kk)*factor
          enddo
          endif
          if (nbdry.gt.0) then
           do m=1,nbdry
             gbdrpe(m)=gbdrpe(m) + rbdrpc(m,kk)*factor
           enddo
          endif
19580     fgowpe=fgowpe + factor
        endif
 1582   continue
!-------------------------------------------------------------------------
!--  attached current for ff' ?                                         --
!-------------------------------------------------------------------------
        if (icutfp.eq.0) then
          upsi=xpsi(kk)
          wwwww=www(kk)
        else
          upsi=xpsi(kk)*xpsimin
          wwwww=zero(kk)
        endif
        if ((upsi.lt.0.0).or.(upsi.gt.1.0)) go to 1890
        call setfp(upsi,xpsii)
        if (fwtdlc.gt.0.0) then
          call setff(upsi,xpsis)
        endif
        do 1880 jj=kppcur+1,kpcurn
          jjk=jj-kppcur
          factor=xpsii(jjk)/rgrid(i)*wwwww
          if (ifag.eq.1) go to 1820
          do 1600 m=1,nsilop
            rsilpc(m,jj)=rsilpc(m,jj) + gsilpc(m,kk)*factor
 1600     continue
          do 1700 m=1,magpri
            rmp2pc(m,jj)=rmp2pc(m,jj) + gmp2pc(m,kk)*factor
 1700     continue
          if (kstark.gt.0.or.kdomse.gt.0) then
          do 1800 m=1,nstark
            rbrpc(m,jj)=rbrpc(m,jj) + gbrpc(m,kk)*factor
            rbzpc(m,jj)=rbzpc(m,jj) + gbzpc(m,kk)*factor
 1800     continue
          endif
          if (kece.gt.0) then
          do 1810 m=1,nece
            recepc(m,jj)=recepc(m,jj) + gecepc(m,kk)*factor
 1810     continue
          endif
          if (kecebz.gt.0) then
            recebzpc(jj)=recebzpc(jj) + gecebzpc(kk)*factor
          endif
!---------------------------------------------------------------------------
!-- Boundary constraints                                                  --
!---------------------------------------------------------------------------
          if (nbdry.gt.0) then
           do m=1,nbdry
             gbdrpc(m,jj)=gbdrpc(m,jj) + rbdrpc(m,kk)*factor
           enddo
          endif
 1820     continue
          fgowpc(jj)=fgowpc(jj)+factor
          if (ifag.eq.1) go to 1880
          if (fwtdlc.gt.0.0) then
            xpsdd=xpsis(jjk)
            xpsdb=xpsisb(jjk)
            rspdlc(jjk)=(xpsdd-xpsdb) &
                      /rgrid(i)*www(kk)+rspdlc(jjk)
          endif
 1880   continue
!----------------------------------------------------------------------
!-- Hyperbolic tangent term                                          --
!----------------------------------------------------------------------
        if (kedgef.gt.0) then
          siedge=(upsi-fe_psin)/fe_width
          xpsnow=1./fe_width/sidif/cosh(siedge)**2
          if (icutfp.gt.0) xpsnow=xpsnow*xpsimin
          factor=xpsnow/rgrid(i)*wwwww
          if (ifag.eq.1) go to 11820
          do m=1,nsilop
            rsilfe(m)=rsilfe(m) + gsilpc(m,kk)*factor
          enddo
          do m=1,magpri
            rmp2fe(m)=rmp2fe(m) + gmp2pc(m,kk)*factor
          enddo
          if (kstark.gt.0.or.kdomse.gt.0) then
          do m=1,nstark
            rbrfe(m)=rbrfe(m) + gbrpc(m,kk)*factor
            rbzfe(m)=rbzfe(m) + gbzpc(m,kk)*factor
          enddo
          endif
          if (nbdry.gt.0) then
           do m=1,nbdry
             gbdrfe(m)=gbdrfe(m) + rbdrpc(m,kk)*factor
           enddo
          endif
          if (fwtdlc.gt.0.0) then
            xpsdd=tanh(siedge)
            xpsdb=tfedge
            rdlcfe=(xpsdd-xpsdb) &
                      /rgrid(i)*www(kk)+rdlcfe
          endif
11820     fgowfe=fgowfe + factor
        endif
!----------------------------------------------------------------------
!--  toroidal rotation terms                                         --
!----------------------------------------------------------------------
 1890   continue
        if (kvtor.le.0) go to 1982
        if ((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.1.0)) go to 1982
        call setpwp(xpsi(kk),xpsii)
        do 1980 jj=kpcurn+1,kwcurn
          jjii=jj-kpcurn
          factor=xpsii(jjii)*rgsvt(i)*www(kk)
          if (niter.gt.1) then
          if (kvtor.eq.2) then
            factor=factor*(1.+pwp0r2)
          elseif (kvtor.eq.3) then
            factor=factor*ptop0
          endif
          endif
          if (ifag.eq.1) go to 1930
          do 1900 m=1,nsilop
            rsilpc(m,jj)=rsilpc(m,jj) + gsilpc(m,kk)*factor
 1900     continue
          do 1910 m=1,magpri
            rmp2pc(m,jj)=rmp2pc(m,jj) + gmp2pc(m,kk)*factor
 1910     continue
          if (kstark.gt.0.or.kdomse.gt.0) then
           do 1920 m=1,nstark
            rbrpc(m,jj)=rbrpc(m,jj) + gbrpc(m,kk)*factor
            rbzpc(m,jj)=rbzpc(m,jj) + gbzpc(m,kk)*factor
 1920      continue
          endif
          if (kece.gt.0) then
           do 1925 m=1,nece
            recepc(m,jj)=recepc(m,jj) + gecepc(m,kk)*factor
 1925      continue
          endif
          if (kecebz.gt.0) then
            recebzpc(jj)=recebzpc(jj) + gecebzpc(kk)*factor
          endif
!---------------------------------------------------------------------------
!-- Boundary constraints                                                  --
!---------------------------------------------------------------------------
          if (nbdry.gt.0) then
           do m=1,nbdry
             gbdrpc(m,jj)=gbdrpc(m,jj) + rbdrpc(m,kk)*factor
           enddo
          endif
 1930     continue
          fgowpc(jj)=fgowpc(jj) + factor
 1980   continue
 1982   continue
 2000 continue
!----------------------------------------------------------------------
!--  diamagnetic flux                                                --
!----------------------------------------------------------------------
      if (fwtdlc.le.0.0) go to 2020
      do 2010 jj=1,kffcur
        rspdlc(jj)=rspdlc(jj)*sidif*twopi/fbrdy
        if (icutfp.gt.0) rspdlc(jj)=rspdlc(jj)/xpsimin
 2010 continue
      rdlcfe=rdlcfe*twopi/fbrdy
 2020 continue
!-----------------------------------------------------------------------
!--  axial q constraint                                               --
!-----------------------------------------------------------------------
      qconst=twopi*emaxis*rmaxis**2/(emaxis**2+1.)/abs(fcentr)/darea
      xpzero=0.0
      call setpp(xpzero,xpsii)
      rmvtor=qconst*rmaxis
      rmvjj=rmaxis
      if (kvtor.gt.0) rmggvt=(rmaxis/rvtor)**2-1.
      if (niter.gt.1) then
      if (kvtor.eq.2) then
        prew0=pwcurr(xpzero,kwwcur)
        pres0=prcurr(xpzero,kppcur)
        if (abs(pres0).gt.1.e-10_dp) then
           pwop0=prew0/pres0
           pwp0r2=pwop0*rmggvt
        else
           pwop0=0.0
           pwp0r2=0.0
        endif
        rmvtor=rmvtor*(1.-0.5_dp*pwp0r2**2)
        rmvjj=rmvjj*(1.-0.5_dp*pwp0r2**2)
      endif
      if (kvtor.eq.3) then
        prew0=pwcurr(xpzero,kwwcur)
        pres0=prcurr(xpzero,kppcur)
        if (abs(pres0).gt.1.e-10_dp) then
           pwop0=prew0/pres0
           pwp0r2=pwop0*rmggvt
           ptop0=exp(pwp0r2)
        else
           pwop0=0.0
           pwp0r2=0.0
           ptop0=1.0
        endif
        rmvtor=rmvtor*ptop0*(1.-pwp0r2)
        rmvjj=rmvjj*ptop0*(1.-pwp0r2)
      endif
      endif
      rqajtor = rmvtor
      do jj=1,kppcur
          factor= xpsii(jj)
          rqajx(jj)=rmvtor*factor
          rjjjx(jj)=rmvjj*factor
      enddo
      if (kedgep.gt.0) then
        siedge=(xpzero-pe_psin)/pe_width
        rqapetor=rmvtor/pe_width/sidif
        rqape=rqapetor/cosh(siedge)**2
      endif
!
      call setfp(xpzero,xpsii)
      rmvtor=qconst/rmaxis
      rqaftor = rmvtor
      do jj=1,kffcur
          factor= xpsii(jj)
          rqafx(jj)=rmvtor*factor
          rjjfx(jj)=factor/rmaxis
      enddo
      if (kedgef.gt.0) then
        siedge=(xpzero-fe_psin)/fe_width
        rqafetor=rmvtor/fe_width/sidif
        rqafe=rqafetor/cosh(siedge)**2
      endif
!
      if (kvtor.gt.0) then
        call setpwp(xpzero,xpsii)
        rmvnow=rmaxis*rmggvt
        if (niter.gt.1) then
        if (kvtor.eq.2) then
          rmvnow=rmvnow*(1.+pwp0r2)
        endif
        if (kvtor.eq.3) then
          rmvnow=rmvnow*ptop0
        endif
        endif
        rmvtor=qconst*rmvnow
        do jj=1,kwwcur
          factor= xpsii(jj)
          rqawx(jj)=rmvtor*factor
          rjjwx(jj)=rmvnow*factor
        enddo
      endif
!----------------------------------------------------------------------
!-- response functions for MSE                                       --
!----------------------------------------------------------------------
      if (kstark.gt.0.or.kdomse.gt.0) then
        do 2040 m=1,nstark
        if (fwtgam(m).le.0.0) go to 32036
        do 2030 jj=1,kwcurn
          rgampc(m,jj)=(a3gam(jtime,m)*tangam(jtime,m) &
                       - a8gam(jtime,m))*rbrpc(m,jj) &
            +(a4gam(jtime,m)*tangam(jtime,m)-a1gam(jtime,m))*rbzpc(m,jj)
 2030 continue
        do 2032 jj=1,nfcoil
          rgamfc(m,jj)=(a3gam(jtime,m)*tangam(jtime,m) &
                       - a8gam(jtime,m))*rbrfc(m,jj) &
            +(a4gam(jtime,m)*tangam(jtime,m)-a1gam(jtime,m))*rbzfc(m,jj)
 2032 continue
        do 2034 jj=1,nesum
          rgamec(m,jj)=(a3gam(jtime,m)*tangam(jtime,m) &
                       - a8gam(jtime,m))*rbrec(m,jj) &
            +(a4gam(jtime,m)*tangam(jtime,m)-a1gam(jtime,m))*rbzec(m,jj)
 2034   continue
        if (ifitvs.gt.0) then
        do 2036 jj=1,nvesel
          rgamvs(m,jj)=(a3gam(jtime,m)*tangam(jtime,m) &
                       - a8gam(jtime,m))*rbrvs(m,jj) &
            +(a4gam(jtime,m)*tangam(jtime,m)-a1gam(jtime,m))*rbzvs(m,jj)
 2036   continue
        endif
        rhsgam(jtime,m)=-tangam(jtime,m)*btgam(m)*a2gam(jtime,m)
        if (kedgep.gt.0) then
          rgampe(m)=(a3gam(jtime,m)*tangam(jtime,m) &
                       - a8gam(jtime,m))*rbrpe(m) &
            +(a4gam(jtime,m)*tangam(jtime,m)-a1gam(jtime,m))*rbzpe(m)
        endif
        if (kedgef.gt.0) then
          rgamfe(m)=(a3gam(jtime,m)*tangam(jtime,m) &
                       - a8gam(jtime,m))*rbrfe(m) &
            +(a4gam(jtime,m)*tangam(jtime,m)-a1gam(jtime,m))*rbzfe(m)
        endif
32036   continue
        if (keecur.gt.0) then
          rmsnow=rrgam(jtime,m)
          zmsnow=zzgam(jtime,m)
        call seva2d(bkx,lkx,bky,lky,c,rmsnow,zmsnow,pds,ier,n333)
          xmsinow=(simag-pds(1))/sidif
          call seter(xmsinow,xsier)
          do jj=1,keecur
           e1rbz(m,jj)=a5gam(jtime,m)*pds(2)*xsier(jj)
           e2rbz(m,jj)=a7gam(jtime,m)*pds(2)*xsier(jj)
           e3rbr(m,jj)=a6gam(jtime,m)*pds(3)*xsier(jj)
           rgamer(m,jj)=-(e2rbz(m,jj) + e3rbr(m,jj))*tangam(jtime,m) &
                        + e1rbz(m,jj)
          enddo
        endif
 2040 continue
      endif
!----------------------------------------------------------------------
!-- response functions for MSE-LS                                    --
!----------------------------------------------------------------------
      if (jdebug.eq.'MSEL') write (6,*) 'GREEN mmbmsels = ',mmbmsels
      if (mmbmsels.gt.0.or.kdomsels.gt.0) then
        do 92040 m=1,nmsels
        if ((fwtbmselt(jtime,m).le.0.0).and.(kdomsels.eq.0)) go to 92038
        if (rrmselt(jtime,m).le.0.0) go to 92040
        call seva2d(bkx,lkx,bky,lky,c,rrmselt(jtime,m),       &
                          zzmselt(jtime,m),pds,ier,n333)
        xn=(simag-pds(1))/sidif
        if (xn.ge.1.0) then
          do jj=kppcur+1,kpcurn
            rmlspc(m,jj)=fmlscut
          enddo
        else
          call setff(xn,xnsi)
          if (idebug>=2) then
            write (6,*) 'GREEN MSELS m,xn,xnsi,tmu02,dsi,dA= ',m,xn, &
                          xnsi(1),xnsi(2),xnsi(3),tmu02,sidif,darea
          endif
          do 92030 jj=kppcur+1,kpcurn
            mjj=jj-kppcur
            rmlspc(m,jj)=-l1mselt(jtime,m)*sidif/darea*xnsi(mjj)*tmu02       
            rmlspc(m,jj)=rmlspc(m,jj)/rrmselt(jtime,m)**2
92030     continue
        endif
        rhsmls(jtime,m)=-l1mselt(jtime,m)*                   &
          (rcentr*bcentr(jtime)/rrmselt(jtime,m))**2
        ecorrect=1.0
        if (mmemsels.gt.0) then
          epotp=erpote(ypsi,keecur)
          ecorrect=ecorrect-l4mselt(jtime,m)*rrmselt(jtime,m)*epotp
        endif
        rhsmls(jtime,m)=rhsmls(jtime,m)-l2mselt(jtime,m)*ecorrect*        &
          brmls(m)*btmls(m)
        rhsmls(jtime,m)=rhsmls(jtime,m)-l3mselt(jtime,m)*ecorrect**2*     &
          brmls(m)**2
        rhsmls(jtime,m)=rhsmls(jtime,m)-(l3mselt(jtime,m)*ecorrect**2+    &
          l1mselt(jtime,m))*bzmls(m)**2
92038   continue
        if (mmemsels.gt.0) then
          call seter(xn,xsier)
          do jj=1,keecur                                            
            relser(m,jj)=xsier(jj)                                   
          enddo
        endif
92040   continue
        if (idebug>=2.or.jdebug.eq.'MSEL') then
          m=3
          write (6,*) 'GREEN MSELS m,rmlspc,rhsmls= ',m,rmlspc(m,3), &
                          rhsmls(1,m)
        endif
      endif
!-----------------------------------------------------------------------
!--  response functions for q constraints at PSIWANT                  --
!-----------------------------------------------------------------------
      if (nqwant.gt.0) then
       do 87999 k=1,nqwant
       do 3550 jj=1,kwcurn
        fgowsw(jj,k)=0.0
 3550  continue
       do 4000 i=1,nw
       do 4000 j=1,nh
        kk=(i-1)*nh+j
        if ((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.siwantq(k))) &
                  go to 3582
        call setpp(xpsi(kk),xpsii)
        do 3580 jj=1,kppcur
          factor=xpsii(jj)*rgrid(i)*www(kk)
          fgowsw(jj,k)=fgowsw(jj,k) + factor
 3580   continue
 3582   continue
        if (icutfp.eq.0) then
          upsi=xpsi(kk)
          wwwww=www(kk)
        else
          upsi=xpsi(kk)*xpsimin
          wwwww=zero(kk)
        endif
        if ((upsi.lt.0.0).or.(upsi.gt.psiwant)) go to 4000
        call setfp(upsi,xpsii)
        do 3880 jj=kppcur+1,kpcurn
          jji=jj-kppcur
          factor=xpsii(jji)/rgrid(i)*wwwww
          fgowsw(jj,k)=fgowsw(jj,k)+factor
 3880   continue
 4000  continue
87999  continue
      endif
      if (ifag.eq.1) return
!--------------------------------------------------------------------------
!-- deltaz in fitting                                                    --
!--------------------------------------------------------------------------
      if (fitdelz.and.niter.ge.ndelzon) then
      do 12000 i=1,nw
      do 12000 j=1,nh
        kk=(i-1)*nh+j
        rdjdz(kk)=0.0
        if ((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.1.0)) go to 11582
        call setppp(xpsi(kk),xpsii)
        call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n333)
        do 11580 jj=1,kppcur
          factor=xpsii(jj)*rgrid(i)*www(kk)*pds(3)*brsp(nfcoil+jj)
          factor=-factor/sidif
          rdjdz(kk)=rdjdz(kk)+factor
          do 11350 m=1,nsilop
            gsildz(m)=gsildz(m) + gsilpc(m,kk)*factor
11350     continue
          do 11400 m=1,magpri
            gmp2dz(m)=gmp2dz(m) + gmp2pc(m,kk)*factor
11400     continue
          if (kstark.gt.0) then
          do 11500 m=1,nstark
            gbrdz(m)=gbrdz(m) + gbrpc(m,kk)*factor
            gbzdz(m)=gbzdz(m) + gbzpc(m,kk)*factor
11500     continue
          endif
          if (kece.gt.0) then
          do 11550 m=1,nece
            recedz(m)=recedz(m) + gecepc(m,kk)*factor
11550     continue
          endif
          if (kecebz.gt.0) then
            recebzdz=recebzdz + gecebzpc(kk)*factor
          endif
!---------------------------------------------------------------------------
!-- Boundary constraints                                                  --
!---------------------------------------------------------------------------
          if (nbdry.gt.0) then
           do m=1,nbdry
             gbdrdz(m)=gbdrdz(m) + rbdrpc(m,kk)*factor
           enddo
          endif
          fgowdz=fgowdz + factor
11580   continue
11582   continue
!-------------------------------------------------------------------------
!--  attached current for ff' ?                                         --
!-------------------------------------------------------------------------
        if (icutfp.eq.0) then
          upsi=xpsi(kk)
          wwwww=www(kk)
        else
          upsi=xpsi(kk)*xpsimin
          wwwww=zero(kk)
        endif
        if ((upsi.lt.0.0).or.(upsi.gt.1.0)) go to 11890
        call setfpp(upsi,xpsii)
        do 11880 jj=kppcur+1,kpcurn
          jjk=jj-kppcur
          factor=xpsii(jjk)/rgrid(i)*wwwww*pds(3)*brsp(nfcoil+jj)
          factor=-factor/sidif
          rdjdz(kk)=rdjdz(kk)+factor
          do 11600 m=1,nsilop
            gsildz(m)=gsildz(m) + gsilpc(m,kk)*factor
11600     continue
          do 11700 m=1,magpri
            gmp2dz(m)=gmp2dz(m) + gmp2pc(m,kk)*factor
11700     continue
          if (kstark.gt.0) then
          do 11800 m=1,nstark
            gbrdz(m)=gbrdz(m) + gbrpc(m,kk)*factor
            gbzdz(m)=gbzdz(m) + gbzpc(m,kk)*factor
11800     continue
          endif
          if (kece.gt.0) then
          do 11810 m=1,nece
            recedz(m)=recedz(m) + gecepc(m,kk)*factor
11810     continue
          endif
          if (kecebz.gt.0) then
            recebzdz=recebzdz + gecebzpc(kk)*factor
          endif
!---------------------------------------------------------------------------
!-- Boundary constraints                                                  --
!---------------------------------------------------------------------------
          if (nbdry.gt.0) then
           do m=1,nbdry
             gbdrdz(m)=gbdrdz(m) + rbdrpc(m,kk)*factor
           enddo
          endif
          fgowdz=fgowdz + factor
11880   continue
11890   continue
!----------------------------------------------------------------------
!--  toroidal rotation contributions                                 --
!----------------------------------------------------------------------
        if (kvtor.le.0) go to 11982
        if ((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.1.0)) go to 11982
        call setpwpp(xpsi(kk),xpsii)
        do 11980 jj=kpcurn+1,kwcurn
          jjii=jj-kpcurn
          factor=xpsii(jjii)*rgsvt(i)*www(kk)*pds(3)*brsp(nfcoil+jj)
          factor=-factor/sidif
          rdjdz(kk)=rdjdz(kk)+factor
          do 11900 m=1,nsilop
            gsildz(m)=gsildz(m) + gsilpc(m,kk)*factor
11900     continue
          do 11910 m=1,magpri
            gmp2dz(m)=gmp2dz(m) + gmp2pc(m,kk)*factor
11910     continue
          if (kstark.gt.0) then
           do 11920 m=1,nstark
            gbrdz(m)=gbrdz(m) + gbrpc(m,kk)*factor
            gbzdz(m)=gbzdz(m) + gbzpc(m,kk)*factor
11920      continue
          endif
          if (kece.gt.0) then
          do 11930 m=1,nece
            recedz(m)=recedz(m) + gecepc(m,kk)*factor
11930     continue
          endif
          if (kecebz.gt.0) then
            recebzdz=recebzdz + gecebzpc(kk)*factor
          endif
!---------------------------------------------------------------------------
!-- Boundary constraints                                                  --
!---------------------------------------------------------------------------
          if (nbdry.gt.0) then
           do m=1,nbdry
             gbdrdz(m)=gbdrdz(m) + rbdrpc(m,kk)*factor
           enddo
          endif
          fgowdz=fgowdz + factor
11980   continue
11982   continue
12000 continue
!----------------------------------------------------------------------
!-- response functions for MSE                                       --
!----------------------------------------------------------------------
      if (kstark.gt.0) then
        do 12040 m=1,nstark
        if (fwtgam(m).le.0.0) go to 12040
         rgamdz(m)=(a3gam(jtime,m)*tangam(jtime,m) &
                   -a8gam(jtime,m))*gbrdz(m) &
            +(a4gam(jtime,m)*tangam(jtime,m)-a1gam(jtime,m))*gbzdz(m)
12040   continue
      endif
      endif
!
      return
      end subroutine green
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          inicur initializes the current density distribution.    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**           ks..............time slice number                      **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine inicur(ks)
      include 'eparmdud129.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      save isicinit,zelips
!
      if (ivacum.gt.0) return
!----------------------------------------------------------------------
!--  icinit=1 uniform and elliptical flux surfaces                   --
!--         2 parabolic and elliptical                               --
!--        -2 use ESAVE.DAT at first slice and previous slice later  --
!--       -12 parabolic and elliptical first slice and previous      --
!--           slice subsequently                                     --
!----------------------------------------------------------------------
      if (ks.eq.1) isicinit=icinit
      if (isicinit.lt.0)   then
        if (ks.gt.1) then
          icinit=isicinit
          return
        else
          if (isicinit.eq.-12) then
            icinit=2
            go to 1200
          endif
          if (icinit.eq.-2)  go to 1100
        endif
      endif
        select case (abs(icinit))
        case (1)
          go to 100
        case (2)
          go to 1100
        end select
      return
  100 continue
      if (aelip.gt.0.0) go to 150
      aelip=0.50_dp
      do 140 i=1,limitr
        erho=sqrt((xlim(i)-relip)**2+((ylim(i)-zelip)/eelip)**2)
        aelip=min(aelip,erho)
  140 continue
  150 continue
      delcur=1.
      sumi=0.0
      do 300 i=1,nw
      do 300 j=1,nh
        kk=(i-1)*nh+j
        erho=sqrt((rgrid(i)-relip)**2+((zgrid(j)-zelip)/eelip)**2)
        if (erho.gt.aelip) go to 300
        pcurrt(kk)=delcur*zero(kk)/rgrid(i)
        sumi=sumi+pcurrt(kk)
  300 continue
      cratio=pasmat(ks)/sumi
      do 320 i=1,nwnh
        pcurrt(i)=pcurrt(i)*cratio*zero(i)
  320 continue
      return
!
 1100 continue
      if (icinit.gt.0) go to 1200
      open(unit=nsave,form='unformatted',file='esave.dat', &
           status='old',err=1200)
      read (nsave) mw,mh
      read (nsave) xpsi
      read (nsave) brsp
      read (nsave) www
      read (nsave) emaxis, rmaxis, fcentr
      close(unit=nsave)
      return
 1200 continue
      if (ks.eq.1)  zelips=zelip
      if (zelip.gt.1.e5_dp .or. zelips.gt.1.e5_dp) then
        zelip=1.447310_dp*(expmpi(ks,37)-expmpi(ks,43)) &
             +0.692055_dp*(expmpi(ks,57)-expmpi(ks,53)) &
             +0.728045_dp*(silopt(ks,27)-silopt(ks,37)) &
             +2.047150_dp*(silopt(ks,2) -silopt(ks,11))
        zelip=zelip*1.e6_dp/pasmat(ks)
        zbound=zelip
        eelip=1.5_dp
      endif
!----------------------------------------------------------------
!-- set zelip=0.0 if bad signals              96/06/24         --
!----------------------------------------------------------------
      if (abs(fwtmp2(37)).le.1.0e-30_dp)  zelip=0.0
      if (abs(fwtmp2(43)).le.1.0e-30_dp)  zelip=0.0
      if (abs(fwtmp2(57)).le.1.0e-30_dp)  zelip=0.0
      if (abs(fwtmp2(53)).le.1.0e-30_dp)  zelip=0.0
      if (abs(fwtsi(27)).le.1.0e-30_dp)  zelip=0.0
      if (abs(fwtsi(37)).le.1.0e-30_dp)  zelip=0.0
      if (abs(fwtsi(2)).le.1.0e-30_dp)  zelip=0.0
      if (abs(fwtsi(11)).le.1.0e-30_dp)  zelip=0.0
!
      do 1300 i=1,nw
      do 1300 j=1,nh
        kk=(i-1)*nh+j
        erho=sqrt((rgrid(i)-relip)**2+((zgrid(j)-zelip)/eelip)**2)
        xpsi(kk)=(erho/aelip)**2
 1300 continue
      return
      end subroutine inicur

!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          pflux computes the poloidal fluxes on the r-z grid.     **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**          09/03/90..........Lazarus vertical feedback             **
!**          28/03/93..........fix relax                             **
!**                                                                  **
!**********************************************************************
      subroutine pflux(niter,nnin,ntotal,jtime,kerror)
!vas  f90 modifi
      use var_bunemn
      use commonblocks,only: c,wk,copy,bkx,bky,psiold,psipold, &
                             psipp,work
      include 'eparmdud129.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      integer initresult
      dimension pds(6)
      real*8,dimension(:),allocatable :: psikkk,gfbsum
      data initfb/0/,init/0/

      kerror = 0
      ALLOCATE(psikkk(nwnh),gfbsum(nwnh))

      vfeed=(isetfb.ne.0).and.(init.ne.0).and.(niter.gt.2.or.nnin.gt.2)
      if (ivesel.gt.10) return
!----------------------------------------------------------------------------
!--  save flux from current iterations before update                       --
!----------------------------------------------------------------------------
      do 2100 kk=1,nwnh
        psiold(kk)=psi(kk)
        psipold(kk)=psipla(kk)
 2100 continue
!
      if (ibunmn.eq.1) go to 2000
      if ((ibunmn.eq.2).and.(errorm.gt.errcut)) go to 2000
      if (ibunmn.eq.3) go to 2000
      if ((ibunmn.eq.4).and.(errorm.gt.errcut)) go to 2000
!-----------------------------------------------------------------------------
!--  ibunmn=0, and 2,4 when errorm less than errcut                         --
!--  Green's integral method of obtaining flux, can be computationally      --
!--  intensive                                                              --
!-----------------------------------------------------------------------------
      do 1000 i=1,nw
      do 1000 j=1,nh
        kk=(i-1)*nh+j
        psi(kk)=0.0
        do 300 m=1,nfcoil
          psi(kk)=psi(kk)+gridfc(kk,m)*brsp(m)
  300   continue
        if (ivesel.le.0) go to 340
        do 335 m=1,nvesel
          psi(kk)=psi(kk)+gridvs(kk,m)*vcurrt(m)
  335   continue
  340   continue
        if (iecurr.ne.1) go to 400
        do 350 m=1,nesum
          psi(kk)=psi(kk)+gridec(kk,m)*ecurrt(m)
  350   continue
  400   continue
        if (iecurr.ne.2) go to 405
        do  m=1,nesum
          psi(kk)=psi(kk)+gridec(kk,m)*cecurr(m)
        enddo
  405   continue
        if (iacoil.gt.0) then
         do 421 m=1,nacoil
          psi(kk)=psi(kk)+gridac(kk,m)*caccurt(jtime,m)
  421    continue
        endif
        psipla(kk)=psi(kk)
        if (ivacum.gt.0) go to 1000
        do 500 ii=1,nw
        do 500 jj=1,nh
          kkkk=(ii-1)*nh+jj
          mj=abs(j-jj)+1
          mk=(i-1)*nh+mj
          psi(kk)=psi(kk)+gridpc(mk,ii)*pcurrt(kkkk)
  500   continue
        psipla(kk)=psi(kk)-psipla(kk)
 1000 continue
!
      ibunmn=0
!------------------------------------------------------------------------------
!--  relaxation before return                                                --
!------------------------------------------------------------------------------
      go to 6000
 2000 continue

!-----------------------------------------------------------------------------
!-- Buneman's method of obtaining psi at the inner grid points              --
!-- only plasma contribution                                                --
!-----------------------------------------------------------------------------
      if (init.gt.0) go to 2020
      nww=nw-1
      nhh=nh-1
      nbmdim=max(nw,nh)+1
!-----These must be brought into the integrals 
      rgrid1=rgrid(1)
      delrgrid=rgrid(2)-rgrid(1)
      delz=zgrid(2)-zgrid(1)
      drdz2=(delrgrid/delz)**2
!---------------------------------------------------------------------
!--  optional vertical feedback control                             --
!---------------------------------------------------------------------
      if (isetfb.ne.0) then
        if(nw.gt.30) ioffr=ioffr*((nw+1)/33)
        if(nh.gt.30) ioffz=ioffz*((nh+1)/33)
        imd=(nw*nh)/2+1+nh*ioffr+ishiftz
        kct1=imd-ioffz
        kct2=imd+ioffz
        deltaz=zgrid(ioffz+nh/2+1)-zgrid(-ioffz+nh/2+1)
      endif
      init=1
 2020 continue
!----------------------------------------------------------------------------
!--  obtain fluxes by inverting del*                                       --
!----------------------------------------------------------------------------
      if (vfeed) then
        dpsip_last=dpsip
        psic1=psi(kct1)
        psic2=psi(kct2)
        dpsip=psic2-psic1
        psibar=(psic2+psic1)/2
      endif
      if (abs(vcurfb(1)).gt.1.e-6_dp) then
        iinow=vcurfb(3)
        if (ntotal.lt.iinow) then
          if (ntotal.eq.1) then
            do 2027 kk=1,nwnh
              psikkk(kk)=0.0
 2027       continue
          endif
        else
          icurfb=vcurfb(2)
          if (mod(ntotal-iinow,icurfb).eq.0) then
            do 2033 kk=1,nwnh
              psikkk(kk)=psiold(kk)
 2033       continue
          endif
        endif
      endif

!-----------------------------------------------------------------------
!-- boundary terms                                                    --
!-----------------------------------------------------------------------
 2050 continue 
      do 2200 j=1,nh
        kk=(nw-1)*nh+j
        psipp(j)=0.
        psipp(kk)=0.
        do 2170 ii=1,nw
        do 2170 jj=1,nh
          kkkk=(ii-1)*nh+jj
          mj=abs(j-jj)+1
          mk=(nw-1)*nh+mj
          psipp(j)=psipp(j)-gridpc(mj,ii)*pcurrt(kkkk)
          psipp(kk)=psipp(kk)-gridpc(mk,ii)*pcurrt(kkkk)
          psi(j)=psipp(j)
          psi(kk)=psipp(kk)
 2170   continue
 2200 continue
      do 2400 i=2,nw-1
        kk1=(i-1)*nh
        kknh=kk1+nh
        kk1=kk1+1
        psipp(kk1)=0.
        psipp(kknh)=0.
        do 2370 ii=1,nw
        do 2370 jj=1,nh
          kkkk=(ii-1)*nh+jj
          mj1=abs(jj-1)+1
          mjnh=abs(nh-jj)+1
          mk1=(i-1)*nh+mj1
          mknh=(i-1)*nh+mjnh
          psipp(kk1)=psipp(kk1)-gridpc(mk1,ii)*pcurrt(kkkk)
          psipp(kknh)=psipp(kknh)-gridpc(mknh,ii)*pcurrt(kkkk)
          psi(kk1)=psipp(kk1)
          psi(kknh)=psipp(kknh)
 2370   continue
 2400 continue
!-------------------------------------------------------------------------
!-- get flux at inner points by inverting del*, only plasma flux        --
!-- first set up plasma currents, single cyclic method gets factor of 2 --
!-------------------------------------------------------------------------
      if (isolve.eq.0) then
!       original buneman solver method
        do 2600 i=2,nw-1
        do 2600 j=2,nh-1
          kk=(i-1)*nh+j
          psi(kk)=tmu2*pcurrt(kk)*rgrid(i)
 2600   continue
        call buneto(psi,nw,nh,work)
      else 
!       new faster single cyclic reduction method
        do 2700 i=2,nw-1
        do 2700 j=2,nh-1
          kk=(i-1)*nh+j
          psi(kk)=tmu2*pcurrt(kk)*rgrid(i)*2.0
 2700   continue
        call pflux_cycred(psi,work,kerror)
        if (kerror.gt.0) return
      endif
      do 3000 i=1,nwnh
        psi(i)=-psi(i)
 3000 continue

!----------------------------------------------------------------------------
!--  optional symmetrized solution                                         --
!----------------------------------------------------------------------------
      if (symmetrize) then
        difpsi=0.0
        do 3100 i=1,nw
        do 3100 j=1,nh/2
          k1=(i-1)*nh+j
          k2=i*nh-j+1
          val=(psi(k1)+psi(k2))/2
          if (sidif.ne.0.0) then
            difnow=(psi(k1)-psi(k2))/sidif
            difnow=abs(difnow)
            difpsi=max(difnow,difpsi)
          endif
          psi(k1)=val
          psi(k2)=val
3100    continue
      endif ! end symmetrize loop
!------------------------------------------------------------------------
!--  optional damping out the m=1 vertical eigen mode                  --
!------------------------------------------------------------------------
      if (abs(vcurfb(1)).lt.1.e-6_dp) go to 43000
      if (ntotal.lt.5) go to 43000
!-----------------------------------------------------------------------
!-- sum Green's functions for m=1 eigen mode                          --
!-----------------------------------------------------------------------
      if (initfb.eq.0) then
      idelr=nw/8
      idelz=nh/12
      jtop=(nh+1)/2+nh/4-1
      do 42200 i=1,nw
      do 42200 j=1,nh
        kk=(i-1)*nh+j
        gfbsum(kk)=0.0
        do 42170 ii=1+idelr,nw-idelr
        do 42170 jj=jtop-idelz,jtop+idelz
          mj=abs(j-jj)+1
          mk=(nw-1)*nh+mj
          gfbsum(kk)=gfbsum(kk)+gridpc(mj,ii)
42170   continue
42200 continue
      jlow=(nh+1)/2-nh/4+1
      do 42500 i=1,nw
      do 42500 j=1,nh
        kk=(i-1)*nh+j
        do 42400 ii=1+idelr,nw-idelr
        do 42400 jj=jlow-idelz,jlow+idelz
          mj=abs(j-jj)+1
          mk=(nw-1)*nh+mj
          gfbsum(kk)=gfbsum(kk)-gridpc(mj,ii)
42400   continue
42500 continue
      initfb=-1
      endif
!---------------------------------------------------------------------
!--  get damping currents                                           --
!---------------------------------------------------------------------
      tvfbrt(ntotal)=0.0
      do 42600  i=1+idelr,nw-idelr
      do 42600  j=jtop-idelz,jtop+idelz
          kk=(i-1)*nh+j
          tvfbrt(ntotal)=tvfbrt(ntotal)+(psikkk(kk)-psiold(kk))
42600 continue
      do 42620  i=1+idelr,nw-idelr
      do 42620  j=jlow-idelz,jlow+idelz
          kk=(i-1)*nh+j
          tvfbrt(ntotal)=tvfbrt(ntotal)-(psikkk(kk)-psiold(kk))
42620 continue
      if (initfb.eq.-1) then
        vcurfi=vcurfb(1)*cpasma(jtime)/abs(tvfbrt(ntotal))
        initfb=1
      endif
      tvfbrt(ntotal)=tvfbrt(ntotal)*vcurfi
      do 42800 kk=1,nwnh
        psi(kk)=psi(kk)+gfbsum(kk)*tvfbrt(ntotal)
42800 continue
43000 continue
!--------------------------------------------------------------------
!--  optional vertical feedback control                            --
!--    psi(kct1)  !psi at lower point                              --
!--    psi(kct2) ! at upper point                                  --
!--    psilu     !psi (per amp) at lower point due to upper coils  --
!--    brfb(1) is lower current                                    --
!--------------------------------------------------------------------
      if (.not.vfeed) go to 3401
      psiul_psiuu=-grdfdb(kct2,1)
      psilu_psill=grdfdb(kct1,1)
      zdwn_now=0
      zup_now=0
      zcontr=0
      zcurnow=0
      sumc=0
      rcurnow=0
      do 12300 icontr=1,nfound
        zcontr=zcontr+yout(icontr)/nfound
12300 continue
      do 12320 i=1,nw
      do 12320 j=1,nh
        kk=(i-1)*nh+j
        zcurnow=zcurnow+pcurrt(kk)*zgrid(j)
        rcurnow=rcurnow+pcurrt(kk)*rgrid(j)
        sumc=sumc+pcurrt(kk)
12320 continue
      zcurnow=zcurnow/sumc
      rcurnow=rcurnow/sumc
      if (zelip.eq.0.) then
        brfb(1)=-gainp*dpsip/(psiul_psiuu+psilu_psill) &
        -gain*(dpsip-dpsip_last)/(psiul_psiuu+psilu_psill)
      else
        brfb(1)=-gainp* &
        (dpsip-psibar*((zelip-ishiftz*delz)/deltaz)) &
        /(psiul_psiuu+psilu_psill) &
        -gain*(dpsip-dpsip_last)/(psiul_psiuu+psilu_psill)
      endif
      brfb(2)=-brfb(1)
      brfbc(ntotal)=brfb(1)
      if (.not.vfeed) then
          brfb(1)=0.
          brfb(2)=0.
      endif
 3401 continue
!----------------------------------------------------------------------------
!--  add flux from external coils                                          --
!----------------------------------------------------------------------------
      do 3600 kk=1,nwnh
        psipla(kk)=psi(kk)
        if (ivesel.le.0) go to 3140
        do 3135 m=1,nvesel
          psi(kk)=psi(kk)+gridvs(kk,m)*vcurrt(m)
 3135   continue
 3140   continue
        if (iecurr.ne.1) go to 3160
        do 3150 m=1,nesum
          psi(kk)=psi(kk)+gridec(kk,m)*ecurrt(m)
 3150   continue
 3160   continue
        if (iecurr.ne.2) go to 3165
        do  m=1,nesum
          psi(kk)=psi(kk)+gridec(kk,m)*cecurr(m)
        enddo
 3165   continue
        if (vfeed) then
          psi(kk)=psi(kk)+grdfdb(kk,1)*brfb(2)
        endif
        do 3200 m=1,nfcoil
          psi(kk)=psi(kk)+gridfc(kk,m)*brsp(m)
 3200   continue
        if (iacoil.gt.0) then
         do 3535 m=1,nacoil
          psi(kk)=psi(kk)+gridac(kk,m)*caccurt(jtime,m)
 3535    continue
        endif
 3600   continue
 6000 continue
!----------------------------------------------------------------------------
!-- rigid vertical shift correction ?                                      --
!----------------------------------------------------------------------------
      if (ifitdelz.eq.1) then
      if (fitdelz.and.ntotal.ge.ndelzon) then
      cdelznow=cdelz(ntotal-1)/100.
      do i=1,nw
      do j=1,nh
        kk=(i-1)*nh+j
        call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n333)
        psi(kk)=psi(kk)+cdelznow*pds(3)
      enddo
      enddo
      endif
      endif
!----------------------------------------------------------------------------
!--  relax flux if needed                                                  --
!----------------------------------------------------------------------------
      if (ntotal.le.1) go to 7000
      if (abs(relax-1.0).lt.1.0e-03_dp) go to 7000
      do 6500 kk=1,nwnh
        psi(kk)=relax*psi(kk)+(1.-relax)*psiold(kk)
        psipla(kk)=relax*psipla(kk)+(1.-relax)*psipold(kk)
 6500 continue
 7000 continue
!
      DEALLOCATE(psikkk,gfbsum)
!
      return
      end subroutine pflux

!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**         residu computes the flux variations on the r-z grid.     **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**          28/03/93..........fixe relax                            **
!**                                                                  **
!**********************************************************************
      subroutine residu(nx,jtime)
      use commonblocks,only: psiold,psipold,psipp
      include 'eparmdud129.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!      include 'ecomdu1.f90'

      if (ivacum.gt.0) return
      errold=errorm
      errave=0.0
      errorm=0.0
      do i=1,nw
        do j=1,nh
          kk=(i-1)*nh+j
          change=abs(psi(kk)-psiold(kk))
          errorm=max(errorm,change)
          errave=errave+change
          if (errorm.le.change) then
            iermax(nx)=i
            jermax(nx)=j
          end if
        end do
      end do

      errorm=errorm/abs(sidif)

      aveerr(nx)=errave/abs(sidif)/nwnh
      cerror(nx)=errorm
      idone=0
      if (errorm.le.error) idone=1
      !----------------------------------------------------------------------
      !--  Turn on vertical stabilization if error small                   --
      !----------------------------------------------------------------------
      if ((errorm.le.errdelz).and.fitdelz) then
        ndelzon = 3
      else
        ndelzon = 999
      endif
      !----------------------------------------------------------------------
      !--  vertical stabilization and iteration information                --
      !----------------------------------------------------------------------
      if (itell.gt.0.and.isetfb.ge.0) then
        if (itell.eq.1) then
          !if (nx.eq.1) write (nttyo,10017) itime, rank
          if (nx.eq.1) write (nttyo,'(x)')
          if (nsol.eq.0) then
            if (mmbmsels.eq.0) then
              write (nttyo,10019) rank,itime,nx,tsaisq(jtime),zmaxis,errorm,delzmm, &
                sum(chigam)
            else
              write (nttyo,90019) rank,itime,nx,tsaisq(jtime),zmaxis,errorm,delzmm, &
                sum(chigam),tchimls
            endif
          else
            write (nttyo,10020) rank,itime,nx,tsaisq(jtime),zmaxis,errorm,delzmm, &
              sum(chigam),erbmax,erbsmax
          endif
        elseif (itell.eq.2) then
          write (nttyo,10021) rank,itime,nx,ali(jtime),abs(betatn),errorm,qsiw(1)
        elseif (itell.eq.3) then
          write (nttyo,10023) rank,itime,nx,difpsi,zmaxis,errorm,delzmm
        elseif (itell.eq.4) then
          if (nx.eq.1) then
            !write (nttyo,10017) itime, rank
            write (nttyo,'(x)')
            if (kstark.gt.0) then
              write (nttyo,80019) rank,itime,nx,tsaisq(jtime),zmaxis,errorm,delzmm, &
                cdelz(1),cdeljsum,sum(chigam)
            else
              write (nttyo,10025) rank,itime,nx,tsaisq(jtime),zmaxis,errorm &
                ,delzmm,cdelz(1),cdeljsum
            endif
          else
            if (kstark.gt.0) then
              write (nttyo,80019) rank,itime,nx,tsaisq(jtime),zmaxis,errorm,delzmm, &
                cdelz(nx-1),cdeljsum,sum(chigam)
            else
              write (nttyo,10025)  rank,itime,nx,tsaisq(jtime),zmaxis,errorm &
                ,delzmm,cdelz(nx-1),cdeljsum
            endif
          endif
        endif
      endif
      if (isetfb.ne.0) then
        write (4,10009) rank,itime,nx,tsaisq(jtime),zmaxis,errorm,delzmm &
          ,brfb(1)
        if (isetfb.lt.0) &
          write (6,10009) rank,itime,nx,tsaisq(jtime),zmaxis,errorm,delzmm &
          ,brfb(1)
      elseif (eelip.gt.2.25_dp .and. itell.eq.0) then
        write (6,10009) rank,itime,nx,tsaisq(jtime),zmaxis,errorm,delzmm &
          ,brfb(1)
      endif
      if (idebug /= 0) then
        write (nttyo,*) 'cratio,cratio_ext,cratiop_ext,cratiof_ext= ', &
          cratio,cratio_ext,cratiop_ext,cratiof_ext
        write (nttyo,*) 'scalepp_ext,scaleffp_ext= ', &
          scalepp_ext,scaleffp_ext
      endif
      return
10009 format (x,'r=',i3,1x,'t=',i6,1x,'iter',i3.3, &
      ' chsq=',1pe8.2,' zmag=',1pe9.2,' err=',1pe8.2,' dz=',1pe10.3, &
      ' Ifb=',1pe9.2)
!10017 format (/,x,' ----- time =',i6,' ms ----- (',i2,')')
10019 format (x,'r=',i3,1x,'t=',i6,1x,'it=',i3, &
      ' chi2=',1pe8.2,' zm=',1pe9.2,' err=',1pe9.3, &
      ' dz=',1pe10.3,' chigam=',1pe9.2)
80019 format (x,'r=',i3,1x,'t=',i6,1x,'it=',i3, &
      ' chi2=',1pe8.2,' zm=',1pe9.2,' err=',1pe9.3, &
      ' dz=',1pe10.3,' delz=',1pe10.3,' dj=',1pe9.3,' chigam=',1pe9.2)
90019 format (x,'r=',i3,1x,'t=',i6,1x,'it=',i3, &
      ' chi2=',1pe8.2,' zm=',1pe9.2,' err=',1pe9.3, &
      ' dz=',1pe10.3,' chigam=',1pe9.2,' chimls=',1pe9.2)
10020 format (x,'r=',i3,1x,'t=',i6,1x,'it=',i3, &
      ' chi2=',1pe8.2,' zm=',1pe9.2,' err=',1pe9.3, &
      ' dz=',1pe10.3,' chigam=',1pe9.2,' errb=',1pe9.2,' errbs=',1pe9.2)
10021 format (x,'r=',i3,1x,'t=',i6,1x,'it=',i3, &
      ' li=',1pe9.3,' betan=',1pe9.3,' err=',1pe9.3, &
      ' qs=',1pe9.3)
10023 format (x,'r=',i3,1x,'t=',i6,1x,'it=',i3, &
      ' dpsi=',1pe10.3,' zm=',1pe9.2,' err=',1pe9.3, &
      ' dz=',1pe10.3)
10025 format (x,'r=',i3,1x,'t=',i6,1x,'it=',i3, &
      ' chi2=',1pe8.2,' zm=',1pe9.2,' err=',1pe9.3, &
      ' dz=',1pe10.3,' delz=',1pe10.3,' dj=',1pe9.3)
      end subroutine residu

!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**                                                                  **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine splitc(is,rs,zs,cs,rc,zc,wc,hc,ac,ac2,cc)
      use global_constants
      use set_kinds
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension rs(1),zs(1),cs(1)

      if(ac+ac2.eq.0.) go to 100
      if(ac.ne.0.) go to 200
      if(ac2.ne.0.) go to 300
!----------------------------------------------------------------------
!-- rectangle                                                        --
!----------------------------------------------------------------------
  100 continue
      wdelt=wc/is
      hdelt=hc/is
      rstrt=rc-wc/2.+wdelt/2.
      zstrt=zc-hc/2.+hdelt/2.
      zz=zstrt
      ic=0
      c=cc/(is*is)
      do 120 ii=1,is
      rr=rstrt
      do 110 jj=1,is
      ic=ic+1
      zs(ic)=zz
      rs(ic)=rr
      cs(ic)=c
      rr=rr+wdelt
  110 continue
      zz=zz+hdelt
  120 continue
      go to 900
!----------------------------------------------------------------------
!-- ac .ne. 0                                                        --
!----------------------------------------------------------------------
  200 continue
      side=tan(radeg*ac)*wc
      hdelt=hc/is
      wdelt=wc/is
      zdelt=tan(radeg*ac)*wdelt
      rstrt=rc-wc/2.+wdelt/2.
      tsid=hc+side
      zstrt =zc-tsid/2.+tsid/2.*1./is
      rr=rstrt
      ic=0
      c=cc/(is*is)
      do 220 ii=1,is
      zz=zstrt+(ii-1)*zdelt
      do 210 jj=1,is
      ic=ic+1
      zs(ic)=zz
      rs(ic)=rr
      cs(ic)=c
      zz=zz+hdelt
  210 continue
      rr=rr+wdelt
  220 continue
      go to 900
!----------------------------------------------------------------------
!-- ac2 .ne. 0                                                       --
!----------------------------------------------------------------------
  300 continue
!
  340 continue
      side=hc/tan(radeg*ac2)
      hdelt=hc/is
      wdelt=wc/is
      zstrt=zc-hc/2.+hdelt/2.
      rdelt=hdelt/tan(radeg*ac2)
      rstrt=rc-side/2.-wc/2.+rdelt/2.+wdelt/2.
      side=hc/tan(radeg*ac2)
      wtot=side+wc
      whaf=(side+wc)/2.
      rcorn=rc-whaf
      rcornr=rc+whaf
      rcorn2=rcorn+wtot/is
      rstrt=(rcorn+rcorn2)/2.
      zz=zstrt
      ic=0
      c=cc/(is*is)
      do 320 ii=1,is
      rr=rstrt+(ii-1)*rdelt
      do 310 jj=1,is
      ic=ic+1
      zs(ic)=zz
      rs(ic)=rr
      cs(ic)=c
      rr=rr+wdelt
  310 continue
      zz=zz+hdelt
  320 continue
      go to 900
!
  900 continue
      return
      end subroutine splitc
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          step computes the dimensionless poloidal fluxes for     **
!**          the r-z grid.                                           **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**                                                                  **
!**********************************************************************
      subroutine steps(ix,ixt,ixout,jtime,kerror)
      use commonblocks,only: c,wk,copy,bkx,bky,zeros,xouts,youts, &
           rsplt,zsplt,csplt
      include 'eparmdud129.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      common/wwork1/xlims(5),ylims(5),limtrs,xlmins
      dimension pds(6)
      integer iii
      real :: zmaxis_last = 0.0
      data isplit/8/,psitol/1.0e-04_dp/
      save xguess, yguess, xltrac, radbou
!
      if (ivacum.gt.0) return
      if (ixt.gt.1) go to 100
      xguess=(rgrid(1)+rgrid(nw))/2.
      yguess=(zgrid(1)+zgrid(nh))/2.
      if (zbound.ne.0.0) yguess=zbound
      if (rbound.ne.0.0) xguess=rbound
      xltrac=xlmin
      if (ibound.eq.-1) xltrac=xlmax
      radbou=(xguess+xltrac)/2.
!
  100 continue
!----------------------------------------------------------------------
!-- first set up bi-cubic spline interpolation in findax             --
!----------------------------------------------------------------------
      m10=10

      if (idebug >= 2) write (6,*) 'Entering findax'
      !print *, 'nw,nh,rgrid,zgrid',nw,nh,rgrid,zgrid
      !print *, 'rmaxis,zmaxis,simag', rmaxis,zmaxis,simag
      !print *, 'psibry,rseps(1,jtime),zseps(1,jtime),m10', psibry,rseps(1,jtime),zseps(1,jtime),m10
      !print *, 'xout,yout,nfound,psi', xout,yout,nfound,psi
      !print *, 'xmin,xmax,ymin,ymax',xmin,xmax,ymin,ymax
      !print *,  'zxmin,zxmax,rymin,rymax' , zxmin,zxmax,rymin,rymax
      !print *,  'dpsi,bpol,bpolz' , dpsi,bpol,bpolz
      !print *, 'limitr,xlim,ylim,limfag', limitr,xlim,ylim,limfag
      !print *, 'ixt,jtime,kerror', ixt,jtime,kerror
      call findax(nw,nh,rgrid,zgrid,rmaxis,zmaxis,simag, &
                  psibry,rseps(1,jtime),zseps(1,jtime),m10, &
                  xout,yout,nfound,psi,xmin,xmax,ymin,ymax, &
                  zxmin,zxmax,rymin,rymax,dpsi,bpol,bpolz, &
                  limitr,xlim,ylim,limfag,ixt,jtime,kerror)

      if (kerror.gt.0) return
      if (nsol.gt.0) then
   
          write (6,*) 'STEPS R,Z,Si,Err = ', rsol(1),zsol(1),wsisol,ier
          call seva2d(bkx,lkx,bky,lky,c,rbdry(1),zbdry(1),pds,ier,n111)
          write (6,*) 'STEPS R,Z,Si,Err = ', rbdry(1),zbdry(1),pds(1),ier
          call seva2d(bkx,lkx,bky,lky,c,rbdry(nbdry),zbdry(nbdry),pds,ier &
             ,n111)
          write (6,*) 'STEPS R,Z,Si,Err = ', rbdry(nbdry),zbdry(nbdry) &
            ,pds(1),ier

            write(6,*) 'rsplt(kk),zsplt(kk)',rbdry(nbdry),zbdry(nbdry)
            write(6,*) 'lkx, lky',lkx,lky
            write(6,*) 'pds,ier,n111', pds,ier,n111

        if (idebug >= 2) then

          call seva2d(bkx,lkx,bky,lky,c,xout(1),yout(1),pds,ier &
             ,n111)
          write (6,*) 'STEPS simag,psibry,n1,si = ',simag,psibry,n111,pds(1)
          write (6,*) 'STEPS Z,R = ', zgrid(33),(rgrid(i),i=45,45)
          write (6,*) 'STEPS si = ',(psi((i-1)*65+33),i=45,45)
          call seva2d(bkx,lkx,bky,lky,c,rgrid(45),zgrid(33),pds,ier &
             ,n111)
          write (6,*) 'STEPS R,Z,si = ', rgrid(45),zgrid(33),pds(1)
        endif
      endif
!-----------------------------------------------------------------------
!--  Trace boundary, first check for counter beam injection           --
!-----------------------------------------------------------------------
      if (pasmat(jtime).lt.-1.e3_dp) then
        nnerr=10000
      else
        nnerr=0
      endif
      call bound(psi,nw,nh,nwnh,psibry,xmin,xmax,ymin,ymax, &
                 zero,rgrid,zgrid,xguess,yguess,ixt,limitr,xlim,ylim, &
                 xout,yout,nfound,xltrac,npoint,rymin,rymax,dpsi, &
                 zxmin,zxmax,nnerr,ishot,itime, &
                 limfag,radbou,kbound,tolbndpsi)
      if (nnerr.gt.0) then
        kerror=1
        return
      endif
!----------------------------------------------------------------------
!--  find magnetic axis and poloidal flux at axis simag              --
!----------------------------------------------------------------------
      m20=20
      call findax(nw,nh,rgrid,zgrid,rmaxis,zmaxis,simag, &
                  psibry,rseps(1,jtime),zseps(1,jtime),m20, &
                  xout,yout,nfound,psi,xmin,xmax,ymin,ymax, &
                  zxmin,zxmax,rymin,rymax,dpsi,bpol,bpolz, &
                  limitr,xlim,ylim,limfag,ixt,jtime,kerror)
      if (kerror.gt.0) return
      sidif=simag-psibry
      eouter=(ymax-ymin)/(xmax-xmin)
      zplasm=(ymin+ymax)/2.
      aouter=(xmax-xmin)/2.

!-----------------------------------------------------------------------
!--   force free current in the scrape-off layer                      --
!-----------------------------------------------------------------------
      if (icutfp.eq.2) then
        xvsmaxo=xvsmax
        xvsmin=1.e10_dp
        xvsmax=-1.e10_dp
        if (limvs.eq.0) then
        itot=isplit*isplit
        do 51000 k=1,nvesel
          call splitc(isplit,rsplt,zsplt,csplt, &
                  rvs(k),zvs(k),wvs(k),hvs(k),avs(k),avs2(k),cdum)
          do 50900 kk=2,itot
            call seva2d(bkx,lkx,bky,lky,c,rsplt(kk),zsplt(kk), &
                        pds,ier,n111)
            write(6,*) 'rgrid(i),zgrid(j)',rsplt(kk),zsplt(kk)
            write(6,*) 'lkx, lky',lkx,lky
            write(6,*) 'pds,ier,n111', pds,ier,n111
            xvsmin=min(xvsmin,pds(1))
            xvsmax=max(xvsmax,pds(1))
50900     continue
51000   continue
        else
        do 51009 k=1,limitr-1
           delx=xlim(k+1)-xlim(k)
           dely=ylim(k+1)-ylim(k)
           dels=sqrt(delx**2+dely**2)
           nn=dels/0.002_dp
           nn=max(5,nn)
           delx=delx/(nn-1)
           dely=dely/(nn-1)
           do 51007 kk=2,nn
            xww=xlim(k)+delx *(kk-1)
            yww=ylim(k)+dely *(kk-1)
            call seva2d(bkx,lkx,bky,lky,c,xww,yww,pds,ier,n111)
            xvsmin=min(xvsmin,pds(1))
            xvsmax=max(xvsmax,pds(1))
51007      continue
51009   continue
        endif
!--------------------------------------------------------------
!--  exclude private region flux                             --
!--------------------------------------------------------------
        xvsmax=psibry
!--------------------------------------------------------------
!--  possible second separatrix                              --
!--------------------------------------------------------------
        rsepex=-999.
        yvs2=1000.
        if (kskipvs.eq.0) go to 10000
        avebp=cpasma(jtime)*tmu/aouter
        bpmin=avebp
        sibpold=sibpmin
        do 51100 i=1,nw
        do 51090 j=1,nh
          kk=(i-1)*nh+j
          if (zero(kk).gt.0.0005_dp.and.www(kk).lt.0.1_dp) then
           call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n333)
           write(6,*) 'rgrid(i),zgrid(j)',rgrid(i),zgrid(j)
           write(6,*) 'lkx, lky',lkx,lky
           write(6,*) 'pds,ier,n111', pds,ier,n333
           bpnow=sqrt(pds(2)**2+pds(3)**2)/rgrid(i)
           if (bpnow.le.bpmin) then
           if ((abs(dpsi).le.psitol).or.((abs(dpsi).gt.psitol).and. &
               (zgrid(j)*zseps(1,jtime).lt.0.0))) then
              bpmin=bpnow
              xs=rgrid(i)
              ys=zgrid(j)
              sibpmin=pds(1)
           endif
           endif
!          endif
          endif
51090   continue
51100   continue

        if (bpmin.eq.avebp) go to 9320
        relsi=abs((sibpmin-psibry)/sidif)
        if (bpmin.le.0.10_dp*avebp.and.relsi.gt.0.005_dp) then
!------------------------------------------------------------------
!-- find second separatrix                                       --
!------------------------------------------------------------------
          do 9300 j=1,40
          call seva2d(bkx,lkx,bky,lky,c,xs,ys,pds,ier,n666)
          write(6,*) 'xs,ys',xs, ys
          write(6,*) 'lkx, lky',lkx,lky
          write(6,*) 'pds,ier,n111', pds,ier,n111

          det=pds(5)*pds(6)-pds(4)*pds(4)
          if (abs(det).lt.1.0e-15_dp) go to 9305
          xerr=(-pds(2)*pds(6)+pds(4)*pds(3))/det
          yerr=(-pds(5)*pds(3)+pds(2)*pds(4))/det
          xs=xs+xerr
          ys=ys+yerr
          if (xerr*xerr+yerr*yerr.lt.1.0e-12_dp) go to 9310
 9300     continue
 9305     continue
          epssep=xerr*xerr+yerr*yerr
          write (nttyo,11001) epssep,ixt
          if (iand(iout,1).ne.0) write (nout,11001) epssep,ixt
          if (epssep.lt.1.0e-10_dp) go to 9310
          go to 9320
 9310     continue
         sibpmin=pds(1)
         yvs2=ys
         rsepex=xs
         relsi=abs((sibpmin-psibry)/sidif)
         if (relsi.gt.0.005_dp) then
         if (ixt.gt.1) sibpmin=sibpmin*(1.-vsdamp)+sibpold*vsdamp
         xvsmin=max(xvsmin,sibpmin)
         endif
        endif
10000   continue
 9320   continue
        if (alphafp.ge.0.0) then
           xpsimin=xvsmin+alphafp*(xvsmax-xvsmin)
        else
           xpsimino=xpsimins
           xpsimin=abs(alphafp)*xvsmax
           xpsimins=xpsimin
           if (ixt.gt.1) xpsimin=xpsimin*(1.-vsdamp)+xpsimino*vsdamp
           alphamu=(xpsimin-xvsmin)/(xvsmax-xvsmin)
        endif
        xpsialp=xpsimin
        xpsimin=sidif/(simag-xpsimin)
      endif 
!-----------------------------------------------------------------------
!-- get normalized flux function XPSI                                 --
!-----------------------------------------------------------------------
      do 1000 i=1,nw
      do 1000 j=1,nh
        kk=(i-1)*nh+j
        if (icutfp.eq.0) then
          xpsi(kk)=1.1_dp
          if ((rgrid(i).lt.xmin).or.(rgrid(i).gt.xmax)) go to 1000
          if ((zgrid(j).lt.ymin).or.(zgrid(j).gt.ymax)) go to 1000
          xpsi(kk)=(simag-psi(kk))/sidif
        else
          if (zero(kk).gt.0.0005_dp) then
            xpsi(kk)=(simag-psi(kk))/sidif
            if (xpsi(kk)*xpsimin.le.1.0.and.xpsi(kk)*xpsimin.ge.0.0) then
              if ((rgrid(i).lt.rminvs).or.(rgrid(i).gt.rmaxvs)) &
                   xpsi(kk)=1000.
              if ((zgrid(j).lt.zminvs).or.(zgrid(j).gt.zmaxvs)) then
                   xpsi(kk)=1000.
              endif
              if (xpsi(kk).lt.1.0.and.zgrid(j).lt.ymin) xpsi(kk)=1000.
              if (xpsi(kk).lt.1.0.and.zgrid(j).gt.ymax) xpsi(kk)=1000.
              if (abs(zgrid(j)).gt.abs(yvs2).and. &
                zgrid(j)*yvs2.gt.0.0) xpsi(kk)=1000.
            endif
          else
            xpsi(kk)=1000.
          endif
        endif
 1000 continue
!-----------------------------------------------------------------------
!-- get SOL flux if needed                                            --
!-----------------------------------------------------------------------
      if (nsol.gt.0) then
        call seva2d(bkx,lkx,bky,lky,c,rsol(1),zsol(1), &
                     pds,ier,n111)
        wsisol=pds(1)
        if (idebug >= 2) then
          write (6,*) 'STEPS R,Z,Si,Err = ', rsol(1),zsol(1),wsisol,ier 
          call seva2d(bkx,lkx,bky,lky,c,rbdry(1),zbdry(1),pds,ier,n111)
          write (6,*) 'STEPS R,Z,Si,Err = ', rbdry(1),zbdry(1),pds(1),ier
          call seva2d(bkx,lkx,bky,lky,c,rbdry(nbdry),zbdry(nbdry),pds,ier &
             ,n111)
          write (6,*) 'STEPS R,Z,Si,Err = ', rbdry(nbdry),zbdry(nbdry) &
            ,pds(1),ier
        endif

      endif
      if (idebug >= 2) then
        call seva2d(bkx,lkx,bky,lky,c,xout(1),yout(1),pds,ier &
             ,n111)
        write (6,*) 'STEPS simag,psibry,n1,si = ',simag,psibry,n111,pds(1)
        write (6,*) 'STEPS Z,R = ', zgrid(33),(rgrid(i),i=45,45)
        write (6,*) 'STEPS si = ',(psi((i-1)*nw+33),i=45,45)
        call seva2d(bkx,lkx,bky,lky,c,rgrid(45),zgrid(33),pds,ier &
             ,n111)
        write (6,*) 'STEPS R,Z,si = ', rgrid(45),zgrid(33),pds(1)
        write (6,*) 'STEPS lkx,lky = ',lkx,lky

      endif

!-----------------------------------------------------------------------
!-- get weighting function                                            --
!-----------------------------------------------------------------------

      call weight(rgrid,zgrid)
!-----------------------------------------------------------------------
!--  get response functions for MSE                                   --
!-----------------------------------------------------------------------
      if (kstark.gt.0.or.kdomse.gt.0) then
      do 50299 k=1,nstark
        if (rrgam(jtime,k).le.0.0) go to 50299
        call seva2d(bkx,lkx,bky,lky,c,rrgam(jtime,k) &
                    ,zzgam(jtime,k),pds,ier,n111)
        sistark(k)=pds(1)
        sisinow=(simag-pds(1))/sidif
        sigam(k)=sisinow
        fpnow=ffcurr(sisinow,kffcur)
        btgam(k)=fpnow*tmu/rrgam(jtime,k)
50299 continue
      endif

!-----------------------------------------------------------------------
!--  get response functions for MSE-LS                                --
!-----------------------------------------------------------------------
      if (mmbmsels.gt.0.or.kdomsels.gt.0) then
      do 60299 k=1,nmsels
        if (idebug >= 2) then
          write (6,*) 'STEPS MSE-LS k,rrmselt= ', k,rrmselt(jtime,k)
        endif
        if (rrmselt(jtime,k).le.0.0) go to 60299
        call seva2d(bkx,lkx,bky,lky,c,rrmselt(jtime,k) &
                    ,zzmselt(jtime,k),pds,ier,n333)
        simls(k)=pds(1)
        sisinow=(simag-pds(1))/sidif
        sinmls(k)=sisinow
        fpnow=ffcurr(sisinow,kffcur)
        btmls(k)=fpnow*tmu/rrmselt(jtime,k)
        brmls(k)=-pds(3)/rrmselt(jtime,k)
        bzmls(k)=pds(2)/rrmselt(jtime,k)
        if (idebug >= 2) then
          write (6,*) 'STEPS MSE-LS k,rrmselt,br,bz,bt= ', k,rrmselt(jtime,k) &
                       ,brmls(k),bzmls(k),btmls(k)
        endif
60299 continue

      endif
!
      do 51200 k=1,nfound
        xouts(k)=1./xout(k)**2
51200 continue
      nzz=0
      call fluxav(xouts,xout,yout,nfound,psi,rgrid,nw,zgrid,nh, &
                  r2bdry,nzz,sdlobp,sdlbp)
      do 51210 k=1,nfound
        xouts(k)=1./xout(k)
51210 continue
      call fluxav(xouts,xout,yout,nfound,psi,rgrid,nw,zgrid,nh, &
                  r1bdry,nzz,sdlobp,sdlbp)
!-----------------------------------------------------------------------
!--  get metric elements at PSIWANT for edge constraint if needed     --
!----------------------------------------------------------------------
      r1sdry(1)=r1bdry
      r2sdry(1)=r2bdry
      nnn=1
      if (abs(sizeroj(1)-1.0).le.1.e-05_dp.and.kzeroj.eq.1) go to 51977
      if (kzeroj.gt.0) then
       do i=1,kzeroj
       if (sizeroj(i).ge.1.0) sizeroj(i)=0.99999_dp
       siwant=simag+sizeroj(i)*(psibry-simag)
       call surfac(siwant,psi,nw,nh,rgrid,zgrid,rsplt,zsplt,nfounc, &
                    npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nnn, &
                    rmaxis,zmaxis,negcur,kerror)
       if (kerror.gt.0) return
       do 51900 k=1,nfounc
        csplt(k)=1./rsplt(k)**2
51900  continue
       call fluxav(csplt,rsplt,zsplt,nfounc,psi,rgrid,nw,zgrid,nh, &
                  r2sdry(i),nzz,sdlobp,sdlbp)
       do 51920 k=1,nfounc
        csplt(k)=1./rsplt(k)
51920  continue
       call fluxav(csplt,rsplt,zsplt,nfounc,psi,rgrid,nw,zgrid,nh, &
                  r1sdry(i),nzz,sdlobp,sdlbp)
       r2surs = r2sdry(i)*sdlobp
       fpnow = ffcurr(psiwant,kffcur)
       fpnow = fpnow*tmu
       enddo
      endif
51977 continue
!-----------------------------------------------------------------------
!--  get metric elements at PSIWANT for q constraint if needed        --
!-----------------------------------------------------------------------
      if (nqwant.gt.0) then
      do 53999 i=1,nqwant
       siwant=simag+siwantq(i)*(psibry-simag)
       call surfac(siwant,psi,nw,nh,rgrid,zgrid,rsplt,zsplt,nfounc, &
                    npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nnn, &
                    rmaxis,zmaxis,negcur,kerror)
       if (kerror.gt.0) return
       do 53810 k=1,nfounc
        csplt(k)=1./rsplt(k)**2
53810  continue
       call fluxav(csplt,rsplt,zsplt,nfounc,psi,rgrid,nw,zgrid,nh, &
                  r2qdry,nzz,sdlobp,sdlbp)
       do 53820 k=1,nfounc
        csplt(k)=1./rsplt(k)
53820  continue
       call fluxav(csplt,rsplt,zsplt,nfounc,psi,rgrid,nw,zgrid,nh, &
                  r1qdry,nzz,sdlobp,sdlbp)
       r2surq = r2qdry*sdlobp
       fpnow = ffcurr(siwantq(i),kffcur)
       fpnow = fpnow*tmu
       qsiw(i)= abs(fpnow)/twopi*r2surq
       pasmsw(i)=sdlbp/tmu/twopi
53999 continue
      endif
!
      cvolp(ixt)=0.0
      carea=0.0
      xym=xout(1)*yout(1)
      xyma=yout(1)
      do 1450 i=2,nfound
        xyp=xout(i)*yout(i)
        xypa=yout(i)
        dx=xout(i)-xout(i-1)
        cvolp(ixt)=cvolp(ixt)+twopi*(xyp+xym)/2.0*dx
        carea=carea+(xyma+xypa)/2.0*dx
        xym=xyp
        xyma=xypa
 1450 continue
      carea=abs(carea)
      if (iconvr.eq.3) go to 1470
      if ((ix.gt.1).or.(ixout.le.1)) call chisqr(jtime)
 1470 continue
      cvolp(ixt)=abs(cvolp(ixt))*1.0e+06_dp
      vout(jtime)=cvolp(ixt)
      rout(jtime)=(xmax+xmin)/2.*100.
      aout(jtime)=100.*(xmax-xmin)/2.0
      csimag(ixt)=simag
      csibry(ixt)=psibry
      crmaxi(ixt)=rmaxis*100.
      czmaxi(ixt)=zmaxis*100.
      cemaxi(ixt)=emaxis
      cchisq(ixt)=tsaisq(jtime)
      csumip(ixt)=sumip
      tratio(ixt)=cratio
!---------------------------------------------------------------------
!--  get beta and li for constraints                                --
!---------------------------------------------------------------------
      if (fli.gt.0.0.or.fbetan.gt.0.0) then
           call betsli(jtime,rgrid,zgrid,kerror)
           if (kerror.gt.0) return
      endif
!----------------------------------------------------------------------
!--  vertical stabilization information                              --
!----------------------------------------------------------------------
      if (itell.gt.0.and.isetfb.ge.0) then
        delzmm = zmaxis - zmaxis_last
        zmaxis_last=zmaxis
      endif
      if (isetfb.ne.0) then
        fb_plasma(jtime)=abs(brfb(1)*nfbcoil/cpasma(jtime))
        delzmm = zmaxis - zmaxis_last
        zmaxis_last=zmaxis
      elseif (eelip.gt.2.25_dp .and. itell.eq.0) then
        delzmm = zmaxis - zmaxis_last
        zmaxis_last=zmaxis
      endif
!----------------------------------------------------------------------
!--  magnetic axis parameters if needed                              --
!----------------------------------------------------------------------
      if (icinit.gt.0) then
        if ((iconvr.ne.3).and.(ixout.le.1)) go to 1580
      endif
        select case (icurrt)
        case (1)
          go to 1570
        case (2)
          go to 1590
        case (3)
          go to 1595
        case (4)
          go to 1580
        case (5)
          go to 1590
        end select
!
 1570 continue
      rdiml=rmaxis/srma
      cjmaxi=cratio*(sbeta*rdiml+2.*salpha/rdiml)/darea
      if (kvtor.gt.0) then
          cjmaxi=cjmaxi+cratio/darea*sbetaw*rdiml*(rdiml**2-1.)
      endif
      go to 1600
 1580 continue
      nqend=1
      n22=2
      if (errorm.lt.0.1_dp.and.icurrt.eq.4) nqend=nqiter
      do 1585 i=1,nqend
      if (i.gt.1) then
        call currnt(n22,jtime,n22,n22,kerror)
        if (kerror.gt.0) return
      end if

      fcentr=fbrdy**2+sidif*dfsqe
      if (fcentr.lt.0.0) fcentr=fbrdy**2
      fcentr=sqrt(fcentr)*fbrdy/abs(fbrdy)
      rdiml=rmaxis/rzero
      cjmaxi=cratio/darea*(rdiml+rbetap/rdiml)
      if (kvtor.eq.1) then
        rgmvt=(rmaxis/rvtor)**2-1.
        cjmaxi=cjmaxi+cratio/darea*rdiml*rbetaw*rgmvt
      elseif (kvtor.eq.11) then
         ypsm=0.0
         n1set=1
         pres0=prcur4(n1set,ypsm,kppcur)
         prew0=pwcur4(n1set,ypsm,kwwcur)
         rgmvt=(rmaxis/rvtor)**2-1.
         pwop0=prew0/pres0
         ptop0=exp(pwop0*rgmvt)
         pp0= 1.-pwop0*rgmvt
         ppw=rbetaw*rgmvt
         cjmaxi=cjmaxi+(pp0+ppw)*rdiml*ptop0
      endif
      cqmaxi(ixt)=(emaxis**2+1.)*abs(fcentr)/twopi/emaxis &
                  /rmaxis**2/abs(cjmaxi)
      qmaxis=cqmaxi(ixt)
      if (icurrt.eq.4) then
      if (qenp.gt.0.0) enp=enp*qenp/qmaxis
      if (qemp.gt.0.0) emp=emp*qmaxis/qemp
      enf=enp
      emf=emp
      endif
 1585 continue

      return
 1590 continue
      if ((ixt.le.1).and.(icinit.gt.0)) go to 1580
      fcentr=ffcurr(x000,kffcur)
      if (kvtor.eq.0) then
        cjmaxi=(rmaxis*ppcurr(x000,kppcur) &
             +fpcurr(x000,kffcur)/rmaxis)*cratio/darea
      else
        cjmaxi=0.0
        do j=1,kppcur
          cjmaxi=rjjjx(j)*brsp(nfcoil+j)+cjmaxi
        enddo
        do j=1,kffcur
          cjmaxi=rjjfx(j)*brsp(nfcoil+kppcur+j)+cjmaxi
        enddo
        do j=1,kwwcur
          cjmaxi=rjjwx(j)*brsp(nfnpcr+j)+cjmaxi
        enddo
        cjmaxi=cjmaxi/darea
      endif
      go to 1600
 1595 continue
 1600 continue
      cqmaxi(ixt)=(emaxis**2+1.)*abs(fcentr)/twopi/emaxis &
                   /rmaxis**2/abs(cjmaxi)
      qmaxis=cqmaxi(ixt)
      return
 2100 format(/,1x,'shot',i6,' at ',i6,' ms ','** Problem in BOUND **')
11001 format(/,1x,'** 2nd seperatrix **',2x,e10.3,2x,i4)
      end subroutine steps
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          this subroutine reorders the z profile data to be       **
!**          in ascending order and sets the ne and te data to       **
!**          correspond to the new order.                            **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          27/11/90..........first created, T. Carlstrom           **
!**          08/07/91..........revised for EFIT, L. Lao              **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      subroutine tsorder(mbox,zprof,nemprof,temprof,nerprof,terprof)
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension ztemp(40),temtemp(40),tertemp(40),zinc(40), &
                  zprof(1), temprof(1), terprof(1)
      real*8 nemtemp(40),nertemp(40),nemprof(1),nerprof(1)
      !---------------------------------------------------------------------
      !-- copy zprof to ztemp (temporary work space)                  --
      !---------------------------------------------------------------------
      do 100 j=1,mbox
        ztemp(j)=zprof(j)
100   continue
      !---------------------------------------------------------------------
      !-- find min z in ztemp                                         --
      !---------------------------------------------------------------------
      do 1000 j=1,mbox
        zmin=999.
        do 800 i=1,mbox-j+1
          if(ztemp(i).lt.zmin) then
            zmin=ztemp(i)
            kmin=i
          end if
800     continue
        !---------------------------------------------------------------------
        !-- put zmin into new vectors                                   --
        !---------------------------------------------------------------------
        zinc(j)=zmin
        nemtemp(j)=nemprof(kmin)
        temtemp(j)=temprof(kmin)
        nertemp(j)=nerprof(kmin)
        tertemp(j)=terprof(kmin)
        !---------------------------------------------------------------------
        !-- create new ztemp with remaining data                        --
        !---------------------------------------------------------------------
        k=0
        do 900 i=1,mbox-j+1
          if(zmin.ne.ztemp(i))then
            k=k+1
            ztemp(k)=ztemp(i)
          end if
900     continue
1000  continue
      !---------------------------------------------------------------------
      !-- rewrite new vectors                                         --
      !---------------------------------------------------------------------
      do 1200 j=1,mbox
        zprof(j)=zinc(j)
        nemprof(j)=nemtemp(j)
        temprof(j)=temtemp(j)
        nerprof(j)=nertemp(j)
        terprof(j)=tertemp(j)
1200  continue
      !
      return
      end subroutine tsorder

!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          vescur computes the currents induced in the vessel      **
!**          segments due to E coils and F coils.                    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/11/85..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine vescur(jtime)
      include 'eparmdud129.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
!
      if (ivesel.eq.5) return
!
 2000 continue
      if (ivesel.eq.1) return
      do 2100 i=1,nvesel
 2100 vcurrt(i)=vloopt(jtime)/rsisvs(i)
      return
      end subroutine vescur

!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          weight computes the weighting function w.               **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**                                                                  **
!**********************************************************************
      subroutine weight(x,y)
      include 'eparmdud129.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      dimension x(1),y(1)
!
      if (iweigh.le.0) go to 5000
!----------------------------------------------------------------------
!-- find a rectangle to search based on output from subroutine bound --
!----------------------------------------------------------------------
      do 100 j=1,nwnh
  100 www(j)=0.0
      do 110 i=1,nw-1
      ileft=i
  110 if((x(i).le.xmin).and.(x(i+1).gt.xmin))go to 120
  120 do 130 i=ileft,nw-1
      iright=i
  130 if((x(i).lt.xmax).and.(x(i+1).ge.xmax))go to 140
  140 do 150 j=1,nh-1
      jbotm=j
  150 if((y(j).le.ymin).and.(y(j+1).gt.ymin))go to 160
  160 do 170 j=jbotm,nh-1
      jtop=j
  170 if((y(j).lt.ymax).and.(y(j+1).ge.ymax))go to 180
  180 jtop=min0(jtop,nh)
      jbotm=max0(jbotm,1)
      ileft=max0(ileft,1)
      iright=min0(iright,nw)
      do 320 i = ileft,iright
      do 320 j = jbotm,jtop
      kk = j+nh*(i-1)
      a = psi(kk-nh)
      b = psi(kk+nh)
      c = psi(kk-1)
      d = psi(kk+1)
      if (i.eq.1) a = 0.
      if (i.eq.nw) b = 0.
      if (j.eq.1) c = 0.
      if (j.eq.nh) d = 0.
      psil = psibry
      in = 0
      p1 = 0.
      p2 = 0.
      p3 = 0.
      p4 = 0.
      a1 = 0.
      a2 = 0.
      a3 = 0.
      a4 = 0.
      if (a.ge.psil) go to 1
      in = in+1
      a1 = a-psil
      go to 2
    1 p1 = a-psil
    2 if (b.ge.psil) go to 3
      in = in+1
      a2 = b-psil
      go to 4
    3 p2 = b-psil
    4 if (c .ge. psil) go to 5
      in = in+1
      a3 = c-psil
      go to 6
    5 p3 = c-psil
    6 if (d.ge.psil) go to 7
      in = in+1
      a4 = d-psil
      go to 8
    7 p4 = d-psil
    8 in = in+1
      select case (in)
      case (1)
        go to 10
      case (2)
        go to 11
      case (3)
        go to 12
      case (4)
        go to 13
      case (5)
        go to 14
      end select
   10 www(kk) = 1.
      go to 20
   11 xxw = (p1+p2+p3+p4)/4.
      yyw = (a1+a2+a3+a4)
      yyw = (yyw/(xxw-yyw))**2
      www(kk) = 1.-0.5_dp*yyw
      go to 20
   12 xxw = (p1+p2+p3+p4)
      yyw = (a1+a2+a3+a4)
      www(kk) = xxw/(xxw-yyw)
      go to 20
   13 xxw = (p1+p2+p3+p4)
      yyw = (a1+a2+a3+a4)/4.
      xxw = (xxw/(xxw-yyw))**2
      www(kk) = 0.5_dp*xxw
      go to 20
   14 www(kk) = 0.
   20 continue
!     do 15 kk = 1,nwnh
      www(kk) = www(kk)*zero(kk)
   15 continue
  320 continue
      return
!
 5000 continue
      do 5500 kk=1,nwnh
        www(kk)=0.0
        if ((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.1.0)) go to 5500
        www(kk)=zero(kk)
 5500 continue
      return
      end subroutine weight

