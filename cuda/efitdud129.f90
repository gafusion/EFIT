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
!**                                                                  **
!**********************************************************************
     use commonblocks
     include 'eparmdud129.f90'
     include 'modules2.f90'
     include 'modules1.f90'
     include 'basiscomdu.f90'
! MPI >>>
#if defined(USEMPI)
     include 'mpif.h'
#endif
! MPI <<<
     data kwake/0/
     parameter (krord=4,kzord=4)
     character inp1*4,inp2*4
     integer :: nargs, iargc, finfo, kerror, terr
! OPT_INPUT >>>
     logical input_flag
     integer*4 mode, shot, steps
     character cmdfile*15, shotfile*15, snapext*82
     real*8 starttime, deltatime
     character(80),dimension(ntime) :: inpfile
     namelist/optin/mode,cmdfile,shotfile,shot,starttime,deltatime,steps,snapext,inpfile
! OPT_INPUT <<<
     kerror = 0
! MPI >>>
#if defined(USEMPI)
! Initialize MPI environment
     call MPI_INIT_THREAD(MPI_THREAD_SINGLE,terr,ierr)
     call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
     call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
! Arrays can only be allocated after MPI has been initialized because dimension is # of processes
     allocate(dist_data(nproc),dist_data_displs(nproc),fwtgam_mpi(nstark,nproc))
#else
     rank  = 0
#endif
! MPI <<<

!----------------------------------------------------------------------
!-- Read in grid size from command line and set global variables     --
!-- ONLY root process reads command-line arguments                   --
!----------------------------------------------------------------------
     if (rank == 0) then
       nargs = iargc()
! Using mpirun command so will have different number of arguments than serial case
#if defined(LF95)
       integer :: iend1, iend2
       character*80 :: cmdline
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
       nw = 0
       nh = 0
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
         print *, 'ERROR: Must specify grid dimensions as arguments'
       endif
! MPI >>>
#if defined(USEMPI)
       deallocate(dist_data,dist_data_displs,fwtgam_mpi)
       call MPI_FINALIZE(ierr)
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
!----------------------------------------------------------------------
!-- Global Allocations                                               --
!----------------------------------------------------------------------
      include 'global_allocs.f90'
!----------------------------------------------------------------------
!-- get data                                                         --
!----------------------------------------------------------------------
      call inp_file_ch(nw,nh,ch1,ch2)
      
! OPT_INPUT >>>
     use_opt_input = .false.
! MPI >>>
! ONLY root process check for existence of input file
     if (rank == 0) then
       inquire(file='efit.input',exist=input_flag)
     endif
#if defined(USEMPI)
! Distribute file existence flag to all processes
     if (nproc > 1) then
       call MPI_BCAST(input_flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
     endif
#endif
     if (input_flag .eqv. .true.) then
       if (rank == 0) then
         print *, ' Using efit.input file...'
       endif
! ALL processes open and read efit.input file
       open(unit=nin,status='old',file='efit.input')
       read(nin,optin)
       close(nin)
       use_opt_input = .true.
       mode_in = mode
       cmdfile_in = cmdfile
       shotfile_in = shotfile
       shot_in = shot
       starttime_in = starttime
       deltatime_in = deltatime
       steps_in = steps
       snapext_in = snapext
       inpfile_in = inpfile
     endif
! OPT_INPUT <<<

   20 call getsets(ktime,kwake,mtear,kerror)

! MPI >>>
#if defined(USEMPI)
    if (nproc > 1) then
! NOTE : Want all processes to exit if any encounter error condition
      call MPI_ALLREDUCE(kerror,MPI_IN_PLACE,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    endif
    if (kerror /= 0) then
      call mpi_stop
    endif
#endif
! MPI <<<
!----------------------------------------------------------------------
!-- start simulation for KTIME timeslices                            --
!----------------------------------------------------------------------
      k=0
  100 k=k+1
        ks=k
!----------------------------------------------------------------------
!--  set up data                                                     --
!----------------------------------------------------------------------
        call data_input(ks,iconvr,ktime,mtear)
        if (iconvr.lt.0) go to 500
        if (kautoknt .eq. 1) then
           call autoknot(ks,iconvr,ktime,mtear,kerror)
        else
!----------------------------------------------------------------------
!--  initialize current profile                                      --
!----------------------------------------------------------------------
           call inicur(ks)
!----------------------------------------------------------------------
!--  get equilibrium                                                 --
!----------------------------------------------------------------------
           call fit(ks,kerror)
           if (kerror.gt.0.and.k.lt.ktime) go to 500
        endif
!----------------------------------------------------------------------
!--  post processing for graphic and text outputs                    --
!----------------------------------------------------------------------
        call shape(ks,ktime,kerror)
        if (mtear.ne.0) call tearing(ks,mtear)
        if (kerror.gt.0) go to 500
        if (idebug /= 0) write (6,*) 'Main/PRTOUT ks/kerror = ', ks, kerror
        call prtout(ks)
        if (kwaitmse.ne.0) call fixstark(-ks,kerror)
!----------------------------------------------------------------------
!--  write A and G EQDSKs                                            --
!----------------------------------------------------------------------
        if (idebug /= 0) write (6,*) 'Main/WQEDSK ks/kerror = ', ks, kerror
        call weqdsk(ks)
        if (iconvr.ge.0) then
           call shipit(ktime,ks,ks)
        call wtear(mtear,jtime)
        endif
        call wmeasure(ktime,ks,ks,1)
!----------------------------------------------------------------------
! -- write Kfile if needed                                           --
!----------------------------------------------------------------------
      if (kdata.eq.3 .or. kdata.eq.7) then
         if (write_Kfile) then
           call write_K2(ks,kerror)
         endif
      endif
  500 if (k.lt.ktime) go to 100
      if (kwake.ne.0) go to 20
      call wmeasure(ktime,1,ktime,2)

! MPI >>>
#if defined(USEMPI)
      ! Finalize MPI
      if (allocated(dist_data)) deallocate(dist_data)
      if (allocated(dist_data_displs)) deallocate(dist_data_displs)
      if (allocated(fwtgam_mpi)) deallocate(fwtgam_mpi)
      if (rank == 0) then
        print *, 'FORTRAN STOP'
      endif
      call MPI_FINALIZE(ierr)
#else
      stop
#endif
! MPI <<<
      end
      subroutine betali(jtime,rgrid,zgrid,idovol,kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          betali computes betas and li.                           **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**          22/07/96 revised  by Q.Peng to add the surface area of  **
!**                            the last closed flux surface (psurfa) **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky,cw,wkw,copyw,bwx, &
                  bwy,sifprw,bwprw,cwprw,dwprw,sfprw,sprwp
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
      common/cwork3/lkx,lky
      common/cww/lwx,lwy
      dimension pds(6),rgrid(1),zgrid(1)
      real*8,dimension(:),allocatable :: worksi,workrm,bwork, &
             cwork,dwork,x,y,dpleng
      dimension xsier(nercur)
! MPI >>>
      integer, intent(inout) :: kerror
      kerror = 0
! MPI <<<
      data inorm/3/,ibtcal/2/
!
      ALLOCATE(worksi(nw),workrm(nw),bwork(nw), &
         cwork(nw),dwork(nw),x(nw),y(nh),dpleng(npoint))
!
      if (ivacum.gt.0) return
      sumbp2=0.0
      go to (20,60) licalc
   20 continue
      call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
      do 50 i=1,nw
      do 50 j=1,nh
        kk=(i-1)*nh+j
        call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n333)
        bpolsq=(pds(2)**2+pds(3)**2)/rgrid(i)
        sumbp2=sumbp2+bpolsq*www(kk)
   50 continue
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
        psurfa(jtime)=psurfa(jtime)+dli*0.5*(xout(i)+xout(ip1))
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
                                cpasma(jtime))**2*1.0E+06
      if (inorm.eq.3) bp2flx=(tmu*2.0*pi*cpasma(jtime)/ &
                             plengt(nfound))**2
      if (inorm.eq.4) bp2flx=2.*(tmu*cpasma(jtime)/aout(jtime))**2 &
                        /(eout(jtime)**2+1.)*1.0e+04
      bpolav(jtime)=sqrt(bp2flx)
      rcurrt(jtime)=sqrt(sumr2/twopi/tmu/abs(cpasma(jtime)))*100.
      zcurrt(jtime)=sumz/twopi/tmu/abs(cpasma(jtime))*100.
      const=twopi/vout(jtime)*1.0e+06/bp2flx
      s1(jtime)=const*s1(jtime)
      s2(jtime)=const*rcentr*s2(jtime)
      s3(jtime)=const*s3(jtime)
      sumbzz=abs(sumbzz)/tmu/abs(cpasma(jtime))/twopi
      rcurrm=rcentr
      ali(jtime)=sumbp2/vout(jtime)/bp2flx*1.0e+06
      ali3(jtime)=(tmu*2.0*pi*cpasma(jtime))**2*rout(jtime)/200.0
      ali3(jtime)=sumbp2/ali3(jtime)
      betap(jtime)=s1(jtime)/4.+s2(jtime)/4.*(1.+rcurrm/rcentr) &
           -ali(jtime)/2.
      betat(jtime)=betap(jtime)*bp2flx/bcentr(jtime)**2*100.
      betat(jtime)=betat(jtime)*(rout(jtime)/100./rcentr)**2
      qout(jtime)=qout(jtime)*abs(bcentr(jtime))*rcentr/twopi
      wplasm(jtime)=1.5*betap(jtime)*bp2flx/2./tmu/2./pi*vout(jtime) &
                     /1.0e+06
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
      dsi=(psibry-simag)/float(nw-1)
      mx=(rmaxis-rgrid(1))/drgrid+1
      my=(zmaxis-zgrid(1))/dzgrid+1
      mkk=(mx-1)*nh+my+1
      sicut=psi(mkk)
      volp(nw)=vout(jtime)/1.0e+06
      volp(1)=0.0
      if (abs(zmaxis).le.0.001) then
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
        ypsi=0.5
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
             if (abs(pres0).gt.1.e-10) then
               pwop0=prew0/pres0
               pwp0r2=pwop0*rgrvt(i)
             else
               pwop0=0.0
               pwp0r2=0.0
             endif
             presst(kk)=pres0*(1.+0.5*pwp0r2**2)+prew0*rgrvt(i)
           elseif (kvtor.eq.11.or.kvtor.eq.3) then
             if (abs(pres0).gt.1.e-10) then
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
      d22=0.03
      d33=0.01
      nzz=0
      zzz=0.0
      do 800 i=2,nw-1
        ii=nw-i+1
        siii=(i-1)*dsi
        siwant=psibry-siii
        siii=1.0-1./float(nw-1)*(i-1)
        if (idovol.gt.1) go to 790
        rzzmax(ii)=-99.0
        call surfac(siwant,psi,nw,nh,rgrid,zgrid,bpol,bpolz,nfind &
                    ,npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nnn, &
                    rmaxis,zmaxis,negcur)
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
        if (kerror /= 0) then
           return
        endif
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
        dyww=bpolz(1)/float(nh-1)
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
          dyww=bpolz(k)/float(nh-1)
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
        if (ierr.eq.0) then
        rzzmax(ii)=-zbbb/2./zaaa
        zzmax(ii)=zaaa*rzzmax(ii)**2+zbbb*rzzmax(ii)+zccc
        endif
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
        go to (760,775,780,770,775) icurrt
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
      wplasm(jtime)=1.5*betap(jtime)*bp2flx/2./tmu/2./pi*vout(jtime) &
                     /1.0e+06
      pasman=cpasma(jtime)/1.e4/aout(jtime)/abs(bcentr(jtime))
      pasman=pasman*rout(jtime)/100./rcentr
      betatn=betat(jtime)/pasman
      dmui=1.0e+06*diamag(jtime)*4.*pi*bcentr(jtime)*rcentr &
            /bp2flx/vout(jtime)
      betapd(jtime)=s1(jtime)/2.+s2(jtime)/2.*(1.-rcurrm/rcentr)-dmui
      betatd(jtime)=betapd(jtime)*bp2flx/bcentr(jtime)**2*100.
      betatd(jtime)=betatd(jtime)*(rout(jtime)/100./rcentr)**2
      wplasmd(jtime)=1.5*betapd(jtime)*bp2flx/2./tmu/2./pi*vout(jtime) &
                    /1.0e+06
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
                      *vout(jtime)/1.0e+06
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
          rhogam(i)=1.1
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
        workrm(i)=rmaxis+(xmax-rmaxis)/float(nw-1)*(i-1)
        call seva2d(bkx,lkx,bky,lky,c,workrm(i),zmaxis,pds,ier,n111)
        worksi(i)=(pds(1)-simag)/(psibry-simag)
 1200 continue
      call zpline(nw,worksi,workrm,bwork,cwork,dwork)
      do 1220 i=1,nw
        sixxx=1./float(nw-1)*(i-1)
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
      end
      subroutine betsli(jtime,rgrid,zgrid,kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          betsli computes betas and li.                           **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
      common/cwork3/lkx,lky
      dimension pds(6),rgrid(1),zgrid(1)
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
      do 50 i=1,nw
      do 50 j=1,nh
        kk=(i-1)*nh+j
        call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n333)
        bpolsq=(pds(2)**2+pds(3)**2)/rgrid(i)
        sumbp2=sumbp2+bpolsq*www(kk)
   50 continue
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
      const=twopi/vout(jtime)*1.0e+06/bp2flx
      ali(jtime)=sumbp2/vout(jtime)/bp2flx*1.0e+06
      ali3(jtime)=(tmu*2.0*pi*cpasma(jtime))**2*rout(jtime)/200.0
      ali3(jtime)=sumbp2/ali3(jtime)
!
      dsi=(psibry-simag)/float(nw-1)
      volp(nw)=vout(jtime)/1.0e+06
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
      pprime(1)=cratio/darea/rzero
      pprime(nw)=pprime(1)*gammap
  600 continue
      sumpre=0.0
      nnn=1
      d11=30.
      d22=0.03
      d33=0.01
      nzz=0
      zzz=0.0
      do 800 i=2,nw-1
        ii=nw-i+1
        siii=(i-1)*dsi
        siwant=psibry-siii
        siii=1.0-1./float(nw-1)*(i-1)
        call surfac(siwant,psi,nw,nh,rgrid,zgrid,xxs,yys,nfind &
                    ,npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nnn, &
                    rmaxis,zmaxis,negcur)
        if (nfind.le.40.and.icntour.eq.0) then
        call cntour(rmaxis,zmaxis,siwant,xcmin,xcmax,ycmin,ycmax, &
                    yxcmin,yxcmax,xycmin,xycmax,d11,drgrid,d22, &
                    d33 ,d33 ,xmin,xmax,ymin,ymax,nzz,iautoc, &
                    xxs,yys,nfind,rgrid,nw,zgrid,nh, &
                    c,n222,nh2,nttyo,npoint, &
                    negcur,bkx,lkx,bky,lky,kerror)
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
      pasman=cpasma(jtime)/1.e4/aout(jtime)/abs(bcentr(jtime))
      pasman=pasman*rout(jtime)/100./rcentr
      betatn=betat(jtime)/pasman
!
      DEALLOCATE(x,y,dpleng,xxs,yys)
!
      return
      end
      subroutine chisqr(jtime)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          chisqr computes the figure of merit for fitting         **
!**          chisq.                                                  **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
!
      if (ivesel.gt.10) return
      if (nbdry.le.0) go to 200
!----------------------------------------------------------------------
!-- read in the plasma response function                             --
!----------------------------------------------------------------------
!vas
!      print*,'file name : ', 'ep'//trim(ch1)// &  
!                         trim(ch2)//'.ddd' 
      open(unit=nrsppc,status='old',form='unformatted', &
           file=table_dir(1:ltbdir)//'ep'//trim(ch1)//trim(ch2)//'.ddd')
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
!
      return
 7400 format (/,2x,7htime = ,e12.5,2x,8hchisq = ,e12.5, &
              2x,10hcurrent = ,e12.5)
 7420 format (10x,14hchi psi loops:)
 7430 format (10x,26hchi inner magnetic probes:)
 7450 format (8(1x,e12.5,1x))
 7460 format (10x,7hchi ip:,/,15x,e12.5)
 7470 format (10x,22hchi pressure:         ,/,1x,e12.5)
 7480 format (10x,22hchi F-coils:          ,/,10x,e12.5)
 7482 format (10x,11hchi psiref:,/,15x,e12.5)
 7485 format (10x,22hchi E-coils:          ,/,1x,e12.5)
 7486 format (10x,22hchi ecebz:            ,/,1x,e12.5)
 7487 format (10x,22hchi total eceR+R-:    ,/,1x,e12.5)
 7488 format (10x,22hchi eceR+R-:          )
      end
      subroutine currnt(iter,jtime,ixn,nitett,kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          currnt computes the current density on the r-z mesh.    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      include 'basiscomdu.f90'
      common/cwork3/lkx,lky
!               lkx,lky
      dimension pds(6)
      dimension alipc(npcur3,nwcurn),xpspp(nppcur),xpsfp(nffcur) &
          ,wlipc(nwcurn),work(nwcur2),xrsp(npcur3),xpspwp(nwwcur)
      dimension crsp(4*(npcurn-2)+6+npcurn*npcurn,nrsmat)
      dimension b(nrsmat),z(4*(npcurn-2)+6+npcurn*npcurn)
      real*8 :: tcurrt, tcurrtpp, tcurrtffp
      data initc/0/
      data ten24/1.e4/
! MPI >>>
      kerror = 0
! MPI <<<
      initc=initc+1
      if (ivacum.gt.0) return
      if ((nitett.le.1).and.(icinit.eq.1)) return
      if (icinit.gt.0) then
        if ((iconvr.ne.3).and.(iter.le.1)) go to 3100
      endif
      go to (100,1100,2100,3100,5100) icurrt
  100 continue
!------------------------------------------------------------------------------
!--  uniform current for Solove equilibrium                                  --
!------------------------------------------------------------------------------
      tcurrt=0.0
      do 130 i=1,nw
      do 130 j=1,nh
        kk=(i-1)*nh+j
        pcurrt(kk)=0.0
        pcurrw(kk)=0.0
        if ((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.1.0)) go to 130
        rdiml=rgrid(i)/srma
        pcurrt(kk)=sbeta*rdiml+2.*salpha/rdiml
        if (kvtor.gt.0) then
          pcurrw(kk)=sbetaw*rdiml*(rdiml**2-1.)
          pcurrt(kk)=pcurrt(kk)+pcurrw(kk)
        endif
        pcurrt(kk)=pcurrt(kk)*www(kk)
        pcurrw(kk)=pcurrw(kk)*www(kk)
        tcurrt=tcurrt+pcurrt(kk)
  130 continue
      cratio=cpasma(jtime)/tcurrt
      do 140 kk=1,nwnh
        pcurrt(kk)=pcurrt(kk)*cratio
        pcurrw(kk)=pcurrw(kk)*cratio
  140 continue
      cj0=cratio/darea/2.
      fcentr=twopi*2.*saaa**2*cj0
      fbrdy=fcentr*sqrt(1.-4.*salpha)
      bcentr(jtime)=fbrdy*tmu/rcentr
      return
!
 1100 continue
!----------------------------------------------------------------------
!--  polynomial current profile                                      --
!----------------------------------------------------------------------
      if ((nitett.le.1).and.(icinit.gt.0)) go to 3100
      nnn=1
      call green(nnn,jtime,nitett)
      if ((nitett.le.1).and.(icinit.lt.0)) go to 1800
      if (iconvr.eq.3) then
        if (kcgama.gt.0.or.kcalpa.gt.0) go to 1200
      endif
      if ((qenp.le.0.0).or.(nitett.le.1)) go to 1800
      if ((qvfit.le.0.0).or.(fwtqa.le.0.0)) go to 1800
      if ((kedgep.gt.0).or.(kedgef.gt.0)) go to 1800
      tz = 0.0
      aqax=1.0/qvfit/(rqajtor*ppcurr(tz,kppcur) &
                   +rqaftor*fpcurr(tz,kffcur))
      aqax = abs(aqax)
      do i=1,kppcur
      	brsp(nfcoil+i)=aqax*brsp(nfcoil+i)
      enddo
      do i=1,kffcur
      	brsp(nbase+i)=aqax*brsp(nbase+i)
      enddo
      cwant0=cpasma(jtime)-fgowpc(1)*brsp(nfcoil+1)-fgowpc(kppcur+1) &
             *brsp(nbase+1)
      cwant1=0.0
      do 1150 i=2,kpcurn
        if (i.eq.kppcur+1) go to 1150
        cwant1=cwant1+fgowpc(i)*brsp(nfcoil+i)
 1150 continue
      if (abs(cwant1).gt.1.0e-10) then
        cwant1=cwant0/cwant1
      else
        kerror=1
        return
      endif
      do 1170 i=2,kpcurn
        if (i.eq.kppcur+1) go to 1170
        brsp(nfcoil+i)=cwant1*brsp(nfcoil+i)
 1170 continue
      go to 1800
!----------------------------------------------------------------------
!--  Adjust current profile to keep q(0), I, J(1), and others fixed  --
!----------------------------------------------------------------------
 1200 continue
      nj=0
      if (fwtqa.gt.0.0) then
        nj=nj+1
        do 1210 j=1,kppcur
          alipc(nj,j)=rqajx(j)*fwtqa
 1210   continue
        do j=1,kffcur
          alipc(nj,kppcur+j)=rqafx(j)*fwtqa
        enddo
        nnow=kpcurn
        if (kedgep.gt.0) then
          nnow=nnow+1
          alipc(nj,nnow)=fwtqa*rqape
        endif
        if (kedgef.gt.0) then
          nnow=nnow+1
          alipc(nj,nnow)=fwtqa*rqafe
        endif
        xrsp(nj)=fwtqa/qvfit*pasmat(jtime)/abs(pasmat(jtime))
      endif
!
      if (fwtcur.gt.0.0) then
      fwtcux=fwtcur
        nj=nj+1
        do 1220 j=1,kpcurn
          alipc(nj,j)=fwtcux*fgowpc(j)
 1220   continue
        nnow=kpcurn
        if (kedgep.gt.0) then
          nnow=nnow+1
          alipc(nj,nnow)=fwtcux*fgowpe
        endif
        if (kedgef.gt.0) then
          nnow=nnow+1
          alipc(nj,nnow)=fwtcux*fgowfe
        endif
        xrsp(nj)=fwtcux*pasmat(jtime)
      endif
!-----------------------------------------------------------------------
!--  constraints on q at psiwant by successive iterations             --
!-----------------------------------------------------------------------
      if (nqwant.gt.0)   then
        do 13226 i=1,nqwant
        nj=nj+1
        if (initc.ge.jwantm) then
          fwtqqq=fwtxxq/0.001/pasmsw(i)
          do 13220 j=1,kpcurn
            alipc(nj,j)=fwtqqq*fgowsw(j,i)
13220     continue
          xrsp(nj)=fwtqqq/fqsiw(i)*qsiw(i)*pasmsw(i)
        endif
13226   continue
      endif
!-----------------------------------------------------------------------
!--  J at PSIWANT constraint                                          --
!-----------------------------------------------------------------------
      brspmin=max(ten24,abs(brsp(nfcoil+1)))
      fwtxxx=fwtxxj*1000./brspmin
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
!----------------------------------------------------------------------
!--  local or flux surface averaged J constraint ?                   --
!----------------------------------------------------------------------
        if (rzeroj(i).gt.0.0) rxxx=rzeroj(i)
        if (rzeroj(i).lt.0.0) rxxx=rseps(1,jtime)/100.
        rxxxf=rxxx
        if (rzeroj(i).eq.0.0) then
             rxxx=1./r1sdry(i)
             rxxxf=r1sdry(i)/r2sdry(i)
        endif
        if (rxxx.le.0.0) then
            rxxx=rcentr
            rxxxf=rxxx
        endif
!
        ysiwant=sizeroj(i)
        call setfp(ysiwant,xpsfp)
        call setpp(ysiwant,xpspp)
        do 1230 j=1,kpcurn
          if (j.le.kppcur) then
            xjj=xpspp(j)
            alipc(nj,j)=rxxx*fwtxxx*xjj
          else
            xjj=xpsfp(j-kppcur)
            alipc(nj,j)=fwtxxx/rxxxf*xjj
          endif
 1230   continue
        nnow=kpcurn
        if (kedgep.gt.0) then
          nnow=nnow+1
          siedge=(sizeroj(i)-pe_psin)/pe_width
          alipc(nj,nnow)=fwtxxx*rxxx/cosh(siedge)**2/pe_width/sidif
        endif
        if (kedgef.gt.0) then
          nnow=nnow+1
          siedge=(sizeroj(i)-fe_psin)/fe_width
          alipc(nj,nnow)=fwtxxx/rxxxf/cosh(siedge)**2/fe_width/sidif
        endif
        xrsp(nj)=fwtxxx*vzeroj(i)*darea*pasmat(jtime)/carea
       enddo
      endif
!-----------------------------------------------------------------------
!--  constraint on betan by successive iterations                     --
!-----------------------------------------------------------------------
      if (fbetan.gt.0.0) then
        brspmin=max(ten24,abs(brsp(nfcoil+jbeta)))
        fwtxxn=fwtxxb*1000./brspmin
        nj=nj+1
        do 13230 j=1,kpcurn
          alipc(nj,j)=0.0
13230   continue
        calpao=brsp(nfcoil+jbeta)
        alipc(nj,jbeta)=fwtxxn
        if (initc.ge.jwantm) then
          xrsp(nj)=fwtxxn/abs(betatn)*fbetan*calpao
        else
          xrsp(nj)=fwtxxn*calpao
        endif
      endif
!-----------------------------------------------------------------------
!--  constraint on li successive iterations                           --
!-----------------------------------------------------------------------
      if (fli.gt.0.0) then
        brspmin=max(ten24,abs(brsp(nfcoil+jli)))
        fwtxxi=fwtxli*1000./brspmin
        nj=nj+1
        do 13240 j=1,kpcurn
          alipc(nj,j)=0.0
13240   continue
        cgamao=brsp(nbase+jli)
        alipc(nj,kppcur+jli)=fwtxxi
        if (initc.ge.jwantm) then
         if (cgamao.lt.0.0) then
          xrsp(nj)=fwtxxi*ali(jtime)/fli*cgamao
         else
          xrsp(nj)=fwtxxi/ali(jtime)*fli*cgamao
         endif
        else
          xrsp(nj)=fwtxxi*cgamao
        endif
      endif
!-----------------------------------------------------------------------
!-- constraints on P' and FF'                                         --
!-----------------------------------------------------------------------
      if (kcalpa.gt.0) then
        do 1240 j=1,kcalpa
        nj=nj+1
        do 1236 i=1,kppcur
          alipc(nj,i)=calpa(i,j)*fwtxxx
 1236   continue
        do 1238 i=kppcur+1,kpcurn
          alipc(nj,i)=0.
 1238   continue
        nnow=kpcurn
        if (kedgep.gt.0) then
          nnow=nnow+1
          alipc(nj,nnow)=calpa(kppcur+1,j)*fwtxxx
        endif
        if (kedgef.gt.0) then
          nnow=nnow+1
          alipc(nj,nnow)=0.0
        endif
        xrsp(nj)=xalpa(j)*darea*fwtxxx
 1240   continue
      endif
!
      if (kcgama.gt.0) then
        do 1250 j=1,kcgama
        nj=nj+1
        do 1244 i=1,kppcur
          alipc(nj,i)=0.0
 1244   continue
        do 1246 i=kppcur+1,kpcurn
          alipc(nj,i)=cgama(i-kppcur,j)*fwtxxx
 1246   continue
        nnow=kpcurn
        if (kedgep.gt.0) then
          nnow=nnow+1
          alipc(nj,nnow)=0.0
        endif
        if (kedgef.gt.0) then
          nnow=nnow+1
          alipc(nj,nnow)=cgama(kffcur+1,j)*fwtxxx
        endif
        xrsp(nj)=xgama(j)*darea*fwtxxx
        if (kfffnc.eq.6) then
           xrsp(nj)=xrsp(nj)/twopi/tmu
        endif
 1250   continue
      endif
!
      nnn=1
      kknow=kcalpa+kcgama
      if (fwtcur.gt.0.0) kknow=kknow+1
      if (fwtqa.gt.0.0) kknow=kknow+1
      if (kzeroj.gt.0) kknow=kknow+kzeroj
      nownow=kpcurn
      if (kedgep.gt.0) nownow=nownow+1
      if (kedgef.gt.0) nownow=nownow+1
      ncrsp = 0
      if (kknow.lt.nownow) then
      nzzzz = 0
      call ppcnst(ncrsp,crsp,z,nzzzz)
      call ffcnst(ncrsp,crsp,z,nzzzz)
      endif
      if (ncrsp .le. 0) then
      call sdecm(alipc,npcur3,nj,nownow,xrsp,npcur3,nnn,wlipc,work,ier)
      if (ier.ne.129) go to 1560
      write (nttyo,8000) ier
! MPI >>>
#if defined(USEMPI)
      ! ERROR_FIX >>>
      !kerror = 1
      !return
      ! <<< >>>
      call mpi_stop
      ! ERROR_FIX <<<
#else
      stop
#endif
! MPI <<<
 1560 continue
      cond=ier
      toler=1.0e-06*wlipc(1)
      do 1570 i=1,nownow
        t=0.0
        if (wlipc(i).gt.toler) t=xrsp(i)/wlipc(i)
        work(i)=t
 1570 continue
      do 1575 i=1,nownow
        brsp(nfcoil+i)=0.0
        do 1575 j=1,nownow
          brsp(nfcoil+i)=brsp(nfcoil+i)+alipc(i,j)*work(j)
 1575 continue
      else
        do j=1,nj
            b(j) = xrsp(j)
        enddo
        call dgglse(nj,nownow,ncrsp,alipc,npcur3,crsp,4*(npcurn-2)+6+ &
                   npcurn*npcurn,b,z,xrsp,work,nrsma2,info,condno)
        do i=1,nownow
          brsp(nfcoil+i)=xrsp(i)
        enddo
        if (info.ne.0) then
        write (nttyo,8000) info
! MPI >>>
#if defined(USEMPI)
        call mpi_stop
#else
        stop
#endif
! MPI <<<
        endif
      endif
      nownow=kpcurn
      if (kedgep.gt.0) then
         nownow=nownow+1
         pedge=brsp(nfcoil+nownow)
      endif
      if (kedgef.gt.0) then
         nownow=nownow+1
         f2edge=brsp(nfcoil+nownow)
      endif
!-----------------------------------------------------------------------
!-- update total plasma current if needed to                          --
!-----------------------------------------------------------------------
      if (abs(fwtcur).le.1.e-30.and.nqwant.gt.0) then
        cm=0.0
        do 4999 n=nfcoil+1,nfnpcr
         cm=cm+brsp(n)*fgowpc(n-nfcoil)
 4999   continue
        if (kedgep.gt.0) then
           cm=cm+pedge*fgowpe
        endif
        if (kedgef.gt.0) then
           cm=cm+f2edge*fgowfe
        endif
        cpasma(jtime)=cm
        pasmat(jtime)=cm
      endif
!
 1800 continue
      tcurrt=0.0
      tcurrp=0.0
      tcurrtpp=0.0
      do 2000 i=1,nw
      do 2000 j=1,nh
        kk=(i-1)*nh+j
        pcurrt(kk)=0.0
        pcurrtpp(kk)=pcurrt(kk)
        if (icutfp.eq.0) then
          if ((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.1.0)) go to 2000
          pcurrt(kk)=rgrid(i)*ppcurr(xpsi(kk),kppcur)
          pcurrtpp(kk)=pcurrt(kk)
          pcurrt(kk)=pcurrt(kk) &
                     +fpcurr(xpsi(kk),kffcur)/rgrid(i)
          pcurrt(kk)=pcurrt(kk)*www(kk)
          pcurrtpp(kk)=pcurrtpp(kk)*www(kk)
        else
          if ((xpsi(kk).ge.0.0).and.(xpsi(kk).le.1.0)) &
            pcurrt(kk)=rgrid(i)*ppcurr(xpsi(kk),kppcur)
          pcurrtpp(kk)=pcurrt(kk)
          upsi=xpsi(kk)*xpsimin
          if ((upsi.ge.0.0).and.(upsi.le.1.0)) &
           pcurrt(kk)=pcurrt(kk)+fpcurr(xpsi(kk),kffcur)/rgrid(i)
          pcurrt(kk)=pcurrt(kk)*zero(kk)
          pcurrtpp(kk)=pcurrtpp(kk)*zero(kk)
        endif
        tcurrt=tcurrt+pcurrt(kk)
        tcurrtpp=tcurrtpp+pcurrtpp(kk)
 2000 continue
      tcurrtffp=tcurrt-tcurrtpp
      if ((nitett.le.1).and.(icinit.lt.0)) then
        cratio=1.0
        return
      endif
      if (.not.fixpp) then
      cratio=cpasma(jtime)/tcurrt
      if (abs(cpasma(jtime)).le.1.e-3) cratio=1.0
      cratio_ext = cratio * cratio_ext
      cratiop_ext = cratio_ext
      cratiof_ext = cratio_ext
      do 2010 kk=1,nwnh
        pcurrt(kk)=pcurrt(kk)*cratio
        pcurrtpp(kk)=pcurrtpp(kk)*cratio
        if ((xpsi(kk).ge.0.0).and.(xpsi(kk).le.1.0)) then
           tcurrp=tcurrp+pcurrt(kk)
        endif
 2010 continue
      do 2020 i=nfcoil+1,nfnpcr
        brsp(i)=cratio*brsp(i)
 2020 continue
      if (npsi_ext > 0) then
        prbdry=prbdry*cratio*cratio
      endif
      if (kedgep.gt.0) then
         pedge=pedge*cratio
      endif
      if (kedgef.gt.0) then
         f2edge=f2edge*cratio
      endif
      else
      cratio=1.0
      cratiop_ext = 1.0
      cratiof = (cpasma(jtime)-tcurrtpp)/tcurrtffp
      cratiof_ext = cratiof * cratiof_ext
      pcurrt(1:nwnh)=pcurrtpp(1:nwnh)+(pcurrt(1:nwnh)-pcurrtpp(1:nwnh))*cratiof
      do i=nbase+kppcur+1,nfnpcr
        brsp(i)=cratiof*brsp(i)
      enddo
      if (kedgef.gt.0) then
         f2edge=f2edge*cratiof
      endif
      endif
!----------------------------------------------------------------------------
!-- rigid vertical shift correction ?                                      --
!----------------------------------------------------------------------------
      if (ifitdelz.eq.3) then
      if (fitdelz.and.nitett.ge.ndelzon) then
      call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
      cdelznow=cdelz(nitett-1)/100.
      cdeljsum=0.0
      do i=1,nw
      do j=1,nh
        kk=(i-1)*nh+j
        call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n333)
!sri-feb09
!        cdeljnow=cdelznow*pds(3)*rdjdz(kk)
        cdeljnow=cdelznow*rdjdz(kk)
        pcurrt(kk)=pcurrt(kk)+cdeljnow
        cdeljsum=cdeljsum+abs(cdeljnow)
      enddo
      enddo
      cdeljsum=abs(cdeljsum/tcurrp)
      endif
      endif
      return
!
 2100 continue
      return
!
 3100 continue
!------------------------------------------------------------------------
!--  GAQ type current profile                                          --
!------------------------------------------------------------------------
      if (kvtor.eq.11) then
        n1set=1
        ypsi=0.5
        pres0=prcur4(n1set,ypsi,kppcur)
        prew0=pwcur4(n1set,ypsi,kwwcur)
        n1set=0
      endif
      tcurrt=0.0
      do 3300 i=1,nw
      do 3300 j=1,nh
        kk=(i-1)*nh+j
        pcurrt(kk)=0.0
        if ((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.1.0)) go to 3300
        rdiml=rgrid(i)/rzero
        pp0       = (1.-xpsi(kk)**enp)**emp*(1.-gammap)+gammap
        pcurrt(kk)=  rbetap/rdiml*((1.-xpsi(kk)**enf)**emf &
                *(1.-gammaf)+gammaf)
!-----------------------------------------------------------------
!--  toroidal rotation ?                                        --
!-----------------------------------------------------------------
        if (kvtor.eq.0) then
         pcurrt(kk)=pcurrt(kk)+pp0*rdiml
        elseif (kvtor.ge.1.and.kvtor.le.3) then
         ppw= (1.-xpsi(kk)**enw)**emw*(1.-gammaw)+gammaw
         ppw=rbetaw*ppw*rgrvt(i)
         pcurrt(kk)=pcurrt(kk)+(pp0+ppw)*rdiml
        elseif (kvtor.eq.11) then
         ypsi=xpsi(kk)
         pres0=prcur4(n1set,ypsi,kppcur)
         prew0=pwcur4(n1set,ypsi,kwwcur)
         if (abs(pres0).gt.1.e-10) then
          pwop0=prew0/pres0
          ptop0=exp(pwop0*rgrvt(i))
         else
          ptop0=1.0
          pwop0=0.0
         endif
         pp0=pp0*(1.-pwop0*rgrvt(i))
         ppw= (1.-xpsi(kk)**enw)**emw*(1.-gammaw)+gammaw
         ppw=rbetaw*ppw*rgrvt(i)
         pcurrt(kk)=pcurrt(kk)+(pp0+ppw)*rdiml*ptop0
        endif
        pcurrt(kk)=pcurrt(kk)*www(kk)
        tcurrt=tcurrt+pcurrt(kk)
 3300 continue
      if (abs(tcurrt).gt.1.0e-10) then
        cratio=cpasma(jtime)/tcurrt
      else
        kerror=1
        return
      endif
      do 4000 kk=1,nwnh
        pcurrt(kk)=pcurrt(kk)*cratio
 4000 continue
      dfsqe=0.5
      ddpsi=1./float(nw-1)
      sinow=0.0
      do 4020 i=2,nw-1
        sinow=sinow+ddpsi
        dfsqe=dfsqe+((1.-sinow**enf)**emf*(1.-gammaf)+gammaf)
 4020 continue
      dfsqe=dfsqe*2.*twopi/tmu*rzero*cratio*rbetap/darea*ddpsi
      if (nitett.gt.1) go to 4100
      if (icurrt.ne.2.and.icurrt.ne.5) go to 4100
      if (fwtbp.le.0.0) go to 4100
      do 4050 i=1,kppcur
        brsp(nfcoil+i)=0.0
 4050 continue
      do 4060 i=1,kffcur
        brsp(nbase+i)=0.0
 4060 continue
      brsp(nfcoil+1)=cratio/rzero
      brsp(nbase+1)=cratio*rbetap*rzero
      brsp(nfcoil+2)=-brsp(nfcoil+1)
      brsp(nbase+2)=-brsp(nbase+1)
 4100 continue
      return
!---------------------------------------------------------------------
!--  toroidal rotation, node points, bases :  ICURRT=5              --
!---------------------------------------------------------------------
 5100 continue
      if ((nitett.le.1).and.(icinit.gt.0)) go to 3100
      nnn=1
      call green(nnn,jtime,nitett)
      if ((nitett.le.1).and.(icinit.lt.0)) go to 5800
      if (iconvr.ne.3) go to 5800
!----------------------------------------------------------------------
!--  Adjust current profile to keep q(0), I, J(1), and others fixed  --
!----------------------------------------------------------------------
 5200 continue
      nj=0
      if (fwtqa.gt.0.0) then
        nj=nj+1
        do 5210 j=1,kppcur
          alipc(nj,j)=rqajx(j)*fwtqa
 5210   continue
        do j=1,kffcur
          alipc(nj,kppcur+j)=rqafx(j)*fwtqa
        enddo
        if (kvtor.gt.0) then
          do j=1,kwwcur
            alipc(nj,kpcurn+j)=rqawx(j)*fwtqa
          enddo
        endif
        xrsp(nj)=fwtqa/qvfit*pasmat(jtime)/abs(pasmat(jtime))
      endif
!
      if (fwtcur.gt.0.0) then
      fwtcux=fwtcur
        nj=nj+1
        do 5220 j=1,kpcurn
          alipc(nj,j)=fwtcux*fgowpc(j)
 5220   continue
        if (kvtor.gt.0) then
          do j=kpcurn+1,kwcurn
            alipc(nj,j)=fwtcux*fgowpc(j)
          enddo
        endif
        xrsp(nj)=fwtcux*pasmat(jtime)
      endif
!-----------------------------------------------------------------------
!--  constraints on q at psiwant by successive iterations             --
!-----------------------------------------------------------------------
      if (nqwant.gt.0)   then
        do 15226 i=1,nqwant
        nj=nj+1
        if (initc.ge.jwantm) then
          fwtqqq=fwtxxq/0.001/pasmsw(i)
          do 15220 j=1,kpcurn
            alipc(nj,j)=fwtqqq*fgowsw(j,i)
15220     continue
          if (kvtor.gt.0) then
            do j=kpcurn+1,kwcurn
              alipc(nj,j)=fwtqqq*fgowsw(j,i)
            enddo
          endif
          xrsp(nj)=fwtqqq/fqsiw(i)*qsiw(i)*pasmsw(i)
        endif
15226   continue
      endif
!-----------------------------------------------------------------------
!--  J at PSIWANT constraint                                          --
!-----------------------------------------------------------------------
      brspmin=max(ten24,abs(brsp(nfcoil+1)))
      fwtxxx=fwtxxj*1000./brspmin
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
!----------------------------------------------------------------------
!--  local or flux surface averaged J constraint ?                   --
!----------------------------------------------------------------------
        if (rzeroj(i).gt.0.0) rxxx=rzeroj(i)
        if (rzeroj(i).lt.0.0) rxxx=rseps(1,jtime)/100.
        rxxxf=rxxx
        if (rzeroj(i).eq.0.0) then
             rxxx=1./r1sdry(i)
             rxxxf=r1sdry(i)/r2sdry(i)
        endif
        if (rxxx.le.0.0) then
            rxxx=rcentr
            rxxxf=rxxx
        endif
!
        ysiwant=sizeroj(i)
        call setfp(ysiwant,xpsfp)
        call setpp(ysiwant,xpspp)
!---------------------------------------------------------------------------
!--  rotational term                                                      --
!---------------------------------------------------------------------------
        if (kvtor.gt.0) then
          call setpwp(ysiwant,xpspwp)
          if (kvtor.eq.2.or.kvtor.eq.3) then
            prew0=pwcurr(ysiwant,kwwcur)
            pres0=prcurr(ysiwant,kppcur)
            if (abs(pres0).gt.1.e-10) then
               pwop0=prew0/pres0
            else
               pwop0=0.0
            endif
          endif
          if (rzeroj(i).gt.0.0) then
              rxx2=(rzeroj(i)/rvtor)**2-1.
              rxxw=rxx2*rzeroj(i)
              if (kvtor.eq.2) then
                rxxw=rxxw*(1.+pwop0*rxx2)
                rxxx=rxxx*(1.-0.5*(pwop0*rxx2)**2)
              endif
              if (kvtor.eq.3) then
                pwp0r2=pwop0*rxx2
                ptop0=exp(pwp0r2)
                rxxw=rxxw*ptop0
                rxxx=rxxx*ptop0*(1.-pwp0r2)
              endif
          endif
          if (rzeroj(i).lt.0.0) then
              rxxw=rseps(1,jtime)/100.
              rxx2=(rxxw/rvtor)**2-1.
              rxxw=rxx2*rxxw
              if (kvtor.eq.2) then
                rxxw=rxxw*(1.+pwop0*rxx2)
                rxxx=rxxx*(1.-0.5*(pwop0*rxx2)**2)
              endif
              if (kvtor.eq.3) then
                pwp0r2=pwop0*rxx2
                ptop0=exp(pwp0r2)
                rxxw=rxxw*ptop0
                rxxx=rxxx*ptop0*(1.-pwp0r2)
              endif
          endif
          if (rzeroj(i).eq.0.0) then
              rxx2=  r2wdry/rvtor**2-1.
              rxxw=  rxx2/r1sdry(i)
              if (kvtor.eq.2) then
                rxxw=rxxw+pwop0*r4wdry/r1sdry(i)
                rxxx=rxxx-0.5*pwop0**2*r4wdry/r1sdry(i)
              endif
              if (kvtor.eq.3) then
                rxxx=(rpwdry-pwop0*rp2wdry)/r1sdry(i)
                rxxw=rp2wdry/r1sdry(i)
              endif
          endif
        endif
!
        do 5230 j=1,kwcurn
          if (j.le.kppcur) then
            xjj=xpspp(j)
            alipc(nj,j)=rxxx*fwtxxx*xjj
          elseif (j.le.kpcurn) then
            xjj=xpsfp(j-kppcur)
            alipc(nj,j)=fwtxxx/rxxxf*xjj
          elseif (kvtor.gt.0) then
            xjj=xpspwp(j-kpcurn)
            alipc(nj,j)=rxxw*fwtxxx*xjj
          endif
 5230   continue
        xrsp(nj)=fwtxxx*vzeroj(i)*darea*pasmat(jtime)/carea
       enddo
      endif
!-----------------------------------------------------------------------
!--  constraint on betan by successive iterations                     --
!-----------------------------------------------------------------------
      if (fbetan.gt.0.0) then
        brspmin=max(ten24,abs(brsp(nfcoil+jbeta)))
        fwtxxn=fwtxxb*1000./brspmin
        nj=nj+1
        do 15230 j=1,kpcurn
          alipc(nj,j)=0.0
15230   continue
        calpao=brsp(nfcoil+jbeta)
        alipc(nj,jbeta)=fwtxxn
        if (initc.ge.jwantm) then
          xrsp(nj)=fwtxxn/abs(betatn)*fbetan*calpao
        else
          xrsp(nj)=fwtxxn*calpao
        endif
      endif
!-----------------------------------------------------------------------
!--  constraint on li successive iterations                           --
!-----------------------------------------------------------------------
      if (fli.gt.0.0) then
        brspmin=max(ten24,abs(brsp(nfcoil+jli)))
        fwtxxi=fwtxli*1000./brspmin
        nj=nj+1
        do 15240 j=1,kpcurn
          alipc(nj,j)=0.0
15240   continue
        cgamao=brsp(nbase+jli)
        alipc(nj,kppcur+jli)=fwtxxi
        if (initc.ge.jwantm) then
         if (cgamao.lt.0.0) then
          xrsp(nj)=fwtxxi*ali(jtime)/fli*cgamao
         else
          xrsp(nj)=fwtxxi/ali(jtime)*fli*cgamao
         endif
        else
          xrsp(nj)=fwtxxi*cgamao
        endif
      endif
!-----------------------------------------------------------------------
!-- constraints on P' and FF'                                         --
!-----------------------------------------------------------------------
      if (kcalpa.gt.0) then
        do 5240 j=1,kcalpa
        nj=nj+1
        do 5236 i=1,kppcur
          alipc(nj,i)=calpa(i,j)*fwtxxx
 5236   continue
        do 5238 i=kppcur+1,kpcurn
          alipc(nj,i)=0.
 5238   continue
        xrsp(nj)=xalpa(j)*darea*fwtxxx
 5240   continue
      endif
!
      if (kcgama.gt.0) then
        do 5250 j=1,kcgama
        nj=nj+1
        do 5244 i=1,kppcur
          alipc(nj,i)=0.0
 5244   continue
        do 5246 i=kppcur+1,kpcurn
          alipc(nj,i)=cgama(i-kppcur,j)*fwtxxx
 5246   continue
        xrsp(nj)=xgama(j)*darea*fwtxxx
 5250   continue
      endif
!
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          do i=1,kwcurn
            alipc(nj,i)=comega(i,j)*fwtxxx
          enddo
          xrsp(nj)=xomega(j)*darea*fwtxxx
        enddo
      endif
!
      if (nj.le.0) go to 5800
      nnn=1
      call sdecm(alipc,npcur3,nj,kwcurn,xrsp,npcur3,nnn,wlipc,work,ier)
      if (ier.ne.129) go to 5560
      write (nttyo,8000) ier
! MPI >>>
#if defined(USEMPI)
      ! ERROR_FIX >>>
      !kerror = 1
      !return
      ! <<< >>>
      call mpi_stop
      ! ERROR_FIX <<<
#else
      stop
#endif
! MPI <<<
 5560 continue
      cond=ier
      toler=1.0e-06*wlipc(1)
      do 5570 i=1,kwcurn
        t=0.0
        if (wlipc(i).gt.toler) t=xrsp(i)/wlipc(i)
        work(i)=t
 5570 continue
      do 5575 i=1,kwcurn
        brsp(nfcoil+i)=0.0
        do 5575 j=1,kwcurn
          brsp(nfcoil+i)=brsp(nfcoil+i)+alipc(i,j)*work(j)
 5575 continue
!-----------------------------------------------------------------------
!-- update total plasma current if needed to                          --
!-----------------------------------------------------------------------
      if (abs(fwtcur).le.1.e-30.and.nqwant.gt.0) then
        cm=0.0
        do 5999 n=nfcoil+1,nfnwcr
         cm=cm+brsp(n)*fgowpc(n-nfcoil)
 5999   continue
        cpasma(jtime)=cm
        pasmat(jtime)=cm
      endif
!
 5800 continue
      tcurrt=0.0
      tcurrp=0.0
      do 6000 i=1,nw
      do 6000 j=1,nh
        kk=(i-1)*nh+j
        pcurrt(kk)=0.0
        pcurrw(kk)=0.0
        if (icutfp.eq.0) then
!----------------------------------------------------------------------
!--  no attached current                                             --
!----------------------------------------------------------------------
          if ((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.1.0)) go to 6000
          pp0=ppcurr(xpsi(kk),kppcur)
          pcurrt(kk)=fpcurr(xpsi(kk),kffcur)/rgrid(i)
          if (kvtor.eq.0) then
             pcurrt(kk)=pcurrt(kk)+pp0*rgrid(i)
          else
!-----------------------------------------------------------------------
!--      rotation                                                     --
!-----------------------------------------------------------------------
             rdimw=rgrid(i)/rvtor
             ppw=pwpcur(xpsi(kk),kwwcur)
             ppw=ppw*(rdimw**2-1.)
             pcurrw(kk)=ppw*rgrid(i)
             if (kvtor.eq.2) then
               prew0=pwcurr(xpsi(kk),kwwcur)
               pres0=prcurr(xpsi(kk),kppcur)
               if (abs(pres0).gt.1.e-10) then
                 pwop0=prew0/pres0
                 pwp0r2=pwop0*rgrvt(i)
               else
                 pwop0=0.0
                 pwp0r2=0.0
               endif
               pcurrw(kk)=pcurrw(kk)*(1.+pwp0r2)
               pp0=pp0*(1.-0.5*pwp0r2**2)
             endif
             if (kvtor.eq.3) then
               prew0=pwcurr(xpsi(kk),kwwcur)
               pres0=prcurr(xpsi(kk),kppcur)
               if (abs(pres0).gt.1.e-10) then
                 pwop0=prew0/pres0
                 pwp0r2=pwop0*rgrvt(i)
                 ptop0=exp(pwp0r2)
               else
                 pwop0=0.0
                 pwp0r2=0.0
                 ptop0=1.0
               endif
               pcurrw(kk)=pcurrw(kk)*ptop0
               pp0=pp0*ptop0*(1.-pwp0r2)
             endif
             pcurrt(kk)=pcurrt(kk)+pp0*rgrid(i)+pcurrw(kk)
          endif
          pcurrt(kk)=pcurrt(kk)*www(kk)
          pcurrw(kk)=pcurrw(kk)*www(kk)
        else
!------------------------------------------------------------------------
!--  attached current                                                  --
!------------------------------------------------------------------------
          if ((xpsi(kk).ge.0.0).and.(xpsi(kk).le.1.0)) then
             pp0=ppcurr(xpsi(kk),kppcur)
             if (kvtor.eq.0) then
               pcurrt(kk)=rgrid(i)*pp0
             else
!------------------------------------------------------------------------
!--        rotation                                                    --
!------------------------------------------------------------------------
               rdimw=rgrid(i)/rvtor
               ppw=pwpcur(xpsi(kk),kwwcur)
               ppw=ppw*(rdimw**2-1.)
               pcurrw(kk)=ppw*rgrid(i)
               if (kvtor.eq.2) then
                 prew0=pwcurr(xpsi(kk),kwwcur)
                 pres0=prcurr(xpsi(kk),kppcur)
                 if (abs(pres0).gt.1.e-10) then
                   pwop0=prew0/pres0
                   pwp0r2=pwop0*rgrvt(i)
                 else
                   pwop0=0.0
                   pwp0r2=0.0
                 endif
                 pcurrw(kk)=pcurrw(kk)*(1.+pwp0r2)
                 pp0=pp0*(1.-0.5*pwp0r2**2)
               endif
               if (kvtor.eq.3) then
                 prew0=pwcurr(xpsi(kk),kwwcur)
                 pres0=prcurr(xpsi(kk),kppcur)
                 if (abs(pres0).gt.1.e-10) then
                   pwop0=prew0/pres0
                   pwp0r2=pwop0*rgrvt(i)
                   ptop0=exp(pwp0r2)
                 else
                   pwop0=0.0
                   pwp0r2=0.0
                   ptop0=1.0
                 endif
                 pcurrw(kk)=pcurrw(kk)*ptop0
                 pp0=pp0*ptop0*(1.-pwp0r2)
               endif
               pcurrt(kk)=pp0*rgrid(i)+pcurrw(kk)
             endif
          endif
          upsi=xpsi(kk)*xpsimin
          if ((upsi.ge.0.0).and.(upsi.le.1.0)) &
            pcurrt(kk)=pcurrt(kk)+fpcurr(xpsi(kk),kffcur)/rgrid(i)
          pcurrt(kk)=pcurrt(kk)*zero(kk)
        endif
        tcurrt=tcurrt+pcurrt(kk)
 6000 continue
      if ((nitett.le.1).and.(icinit.lt.0)) then
        cratio=1.0
        return
      endif
!-------------------------------------------------------------------
!--  adjust to match total current                                --
!-------------------------------------------------------------------
      cratio=cpasma(jtime)/tcurrt
      do 6010 kk=1,nwnh
        pcurrt(kk)=pcurrt(kk)*cratio
        pcurrw(kk)=pcurrw(kk)*cratio
        if ((xpsi(kk).ge.0.0).and.(xpsi(kk).le.1.0)) then
           tcurrp=tcurrp+pcurrt(kk)
        endif
 6010 continue
      do 6020 i=nfcoil+1,nfcoil+kwcurn
        brsp(i)=cratio*brsp(i)
 6020 continue
      return
!
 8000 format (/,'  ** Problem in Decomposition **',i10)
      end
      subroutine data_input(jtime,kconvr,ktime,mtear)
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
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
      include 'basiscomdu.f90'
      parameter(mfila=10)
      parameter (m_ext=101)
      common/cwork3/lkx,lky
      common/gwork1/rmx(mfila),zmx(mfila),rsilpf(nsilop,mfila), &
        rmp2pf(magpri,mfila), &
        irfila(mfila),jzfila(mfila) &
        ,wsilpc(nsilop),wmp2pc(magpri),wfcpc(nfcoil), &
         wecpc(nesum),wvspc(nvesel),wpcpc
      real*8,allocatable :: gridpf(:,:),gwork(:,:),rgrids(:),zgrids(:) 
      dimension coils(nsilop),expmp2(magpri),acoilc(nacoil) &
         ,tgamma(nmtark),sgamma(nmtark),rrrgam(nmtark),zzzgam(nmtark) &
         ,aa1gam(nmtark),aa2gam(nmtark),aa3gam(nmtark),aa4gam(nmtark) &
         ,aa5gam(nmtark),aa6gam(nmtark),aa7gam(nmtark)
      dimension tlibim(libim),slibim(libim),rrrlib(libim) &
         ,zzzlib(libim),aa1lib(libim),aa8lib(libim),fwtlib(libim)
      dimension pds(6),denr(nco2r),denv(nco2v)
      dimension ilower(mbdry)
      dimension devxmpin(magpri),rnavxmpin(magpri) &
               ,devpsiin(nsilop),rnavpsiin(nsilop) &
               ,devfcin(nfcoil),rnavfcin(nfcoil) &
               ,devein(nesum),rnavecin(nesum)
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
      ,idebug,nbdrymx,nsol,rsol,zsol,fwtsol,efitversion,kbetapr,nbdryp
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
            dtmsefull,mse_strict,t_max_beam_off,ok_30rt,ok_210lt
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
          ,ok_30rt,ok_210lt,vbit,nbdrymx
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
      dimension brsptu(nfcoil)
      character*50 edatname
      character*82 table_nam
      character*10 namedum
      character*2 :: reflect_ext
      logical :: shape_ext
      real*4 spatial_avg_ham(nmtark,ngam_vars,ngam_u,ngam_w)
      data nsq/1/
      data ersil8/1.0e-03/,currn1/0.0/
      data idodo/0/,idovs/0/,zetafc/2.5e-08/
      data co2cor/1.0/,idoac/0/,fq95/0.0/
      data mcontr/35/
      data ten2m3/1.0e-03/
      data idtime/0/,itimeb/0/,brsptu(1)/-1.e-20/
      save idodo, idovs, idoac
!
      ALLOCATE(gridpf(nwnh,mfila),gwork(nbwork,nwnh), &
         rgrids(nw),zgrids(nh))
!
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
  175 continue
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
        rbdry(1)=1.94
        zbdry(1)=ztssym(jtime)+0.5*ztswid(jtime)
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
      alpax(jbeta)=1.e4
      backaverage=.false.
      bitip=0.0
      betap0=0.50
      brsp(1)=-1.e+20
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
      error=1.0e-03
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
      do 11299 i=1,nnece
       fwtece0(i)=0.0
11299 continue
       fwtecebz0=0.0
      do 12399 i=1,mbdry
       fwtbdry(i)=1.0
       sigrbd(i)=1.e10
       sigzbd(i)=1.e10
       fwtsol(i)=1.0
12399 continue
      akchiwt=1.0
      akprewt=0.0
      akgamwt=0.0
      akerrwt=0.0
      aktol=0.1
      fwtqa=0.0
      gammap=1.0e+10
      gamax(jli)=-1.e6
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
      qenp=0.95
      qvfit=0.95
      rzeroj(1)=0.0
      salpha=1./40.
      sbeta=1./8.
      sbetaw=0.0
      scrape=0.030
      serror=0.03
      sgnemin=0.0
      sgprmin=0.0
      sgtemin=0.0
      sidif=-1.0e+10
      srm=-3.5
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
!----------------------------------------------------------------------
!--   Read input file for KDATA = 2                                  --
!----------------------------------------------------------------------
      open(unit=nin,status='old',file=ifname(jtime) &
                                 )
!
      xlim(1)=-1.0
      rbdry(1)=-1.0
      ishot=-1
      itimeu=0
      table_nam = table_dir
      nbdryp=-1
      read (nin,in1)
      if (nbdryp==-1) nbdryp=nbdry
      read (nin,ink,err=11111,end=101)
101    continue
11111 close(unit=nin)
      open(unit=nin,status='old',file=ifname(jtime) &
                                 )
      read (nin,ins,err=11113,end=103)
103    continue
11113 close(unit=nin)
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
      shape_ext=.F.
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
11773   format (a)
11775   format (6a8,3i4)
11776   format (5e16.9)
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
            psin_ext(i) = float(i-1)/float(npsi_ext-1)
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
        ylim(1)=0.5*(z0min+z0max)
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
      ltbdir=0
      lindir=0
      lstdir=0
      do i=1,len(table_dir)
         if (table_dir(i:i).ne.' ') ltbdir=ltbdir+1
         if (input_dir(i:i).ne.' ') lindir=lindir+1
         if (store_dir(i:i).ne.' ') lstdir=lstdir+1
      enddo
!---------------------------------------------------------------------
!-- Set up proper Green's tables area                               --
!--        shot > 156000 new 2014 set                               --
!---------------------------------------------------------------------
      if (ishot.ge.112000) then
        if (ishot.lt.156000) then
          table_dir = table_dir(1:ltbdir)//'112000/'
        else
          if (kdata.ne.2) then
            table_dir = table_dir(1:ltbdir)//'156014/'
          else
            if (efitversion <= 20140331) then
               table_dir = table_dir(1:ltbdir)//'112000/'
            else
               table_dir = table_dir(1:ltbdir)//'156014/'
            endif
          endif
        endif
        ltbdir=ltbdir+7
        ltbdi2=ltbdir
        table_di2 = table_dir
        if (rank == 0) then
          write(*,*) 'table_dir = <',table_dir(1:ltbdir),'>'
        endif
      endif
!---------------------------------------------------------------------
!-- Re-read Green's tables from table_dir if necessary              --
!---------------------------------------------------------------------
      if ((table_dir.eq.table_nam).and.(jtime.gt.1)) go to 11337 
      open(unit=mcontr,status='old',form='unformatted', &
           file=table_dir(1:ltbdir)//'ec'//trim(ch1)//trim(ch2)//'.ddd')
      read (mcontr) mw,mh
      read (mcontr) rgrid,zgrid
      read (mcontr) gridfc
      read (mcontr) gridpc
      close(unit=mcontr)
!----------------------------------------------------------------------
!-- read in the f coil response functions                            --
!----------------------------------------------------------------------
      open(unit=mcontr,form='unformatted', &
           status='old',file=table_dir(1:ltbdir)//'rfcoil.ddd')
      read (mcontr) rsilfc
      read (mcontr) rmp2fc
      close(unit=mcontr)
!----------------------------------------------------------------------
!-- read in the plasma response function                             --
!----------------------------------------------------------------------
      open(unit=mcontr,status='old',form='unformatted', &
           file=table_dir(1:ltbdir)//'ep'//trim(ch1)//trim(ch2)//'.ddd')
      read (mcontr) gsilpc
      read (mcontr) gmp2pc
      close(unit=mcontr)
!----------------------------------------------------------------------
!-- read in the E-coil response function                             --
!----------------------------------------------------------------------
      if (iecurr.gt.0) then
        open(unit=mcontr,status='old',form='unformatted', &
           file=table_dir(1:ltbdir)//'re'//trim(ch1)//trim(ch2)//'.ddd')
        read (mcontr) rsilec
        read (mcontr) rmp2ec
        read (mcontr) gridec
        close(unit=mcontr)
      endif
!----------------------------------------------------------------------
!-- read in the vessel response function                             --
!----------------------------------------------------------------------
      if (ivesel.gt.0) then
        open(unit=mcontr,status='old',form='unformatted', &
           file=table_dir(1:ltbdir)//'rv'//trim(ch1)//trim(ch2)//'.ddd')
        read (mcontr) rsilvs
        read (mcontr) rmp2vs
        read (mcontr) gridvs
        close(unit=mcontr)
      endif
!
      open(unit=mcontr,status='old',  &
          file=table_di2(1:ltbdi2)//'dprobe.dat')
      rsi(1)=-1.
      read (mcontr,in3)
      read (mcontr,10200) (rf(i),zf(i),wf(i),hf(i),af(i),af2(i), &
                i=1,mfcoil)
      if (rsi(1).lt.0.) &
        read (mcontr,10200) (rsi(i),zsi(i),wsi(i),hsi(i),as(i),as2(i), &
                i=1,nsilop)
      read (mcontr,10220) (re(i),ze(i),we(i),he(i),ecid(i), &
                                        i=1,necoil)
      if (ifitvs.gt.0.or.icutfp.eq.2) then
        read (mcontr,10200) (rvs(i),zvs(i),wvs(i),hvs(i), &
                                        avs(i),avs2(i),i=1,nvesel)
      endif
      close(unit=mcontr)
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
        calpa(1,1)=0.1
        calpa(2,1)=0.1
        calpa(3,1)=0.1
        xalpa(1)=0.0
        kcgama=1
        cgama(1,1)=0.1
        cgama(2,1)=0.1
        cgama(3,1)=0.1
        xgama(1)=0.0
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
      do i=1,nnece
        swtece(i)=fwtece0(i)
      enddo
      swtecebz=fwtecebz0
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
      if (psiwant.le.0.0) psiwant=1.e-5
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
      if (ishot.eq.-1) go to 10000
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
!--------------------------------------------------------------
!-- option to symmetrize added 8/14/91 eal                   --
!--------------------------------------------------------------
      if (symmetrize) then	! symmetrize the fixed boundery
        isetfb=0	!	be sure vertical feedback is off
        zelip=0	! must be symmetric about midplane
        if(nbdry.gt.1)then	! remove duplicate point
          nbryup=0
          delx2=(rbdry(1)-rbdry(nbdry))**2
          dely2=(zbdry(1)-zbdry(nbdry))**2
          if ((delx2+dely2).le.1.0e-08) nbdry=nbdry-1
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
      endif	!	end boundary symmertization
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
           dnbeam(i)=dnbeam(i)*1.e-19
43901     continue
        endif
!---------------------------------------------------------------------
!--  reorder TS data points                                         --
!---------------------------------------------------------------------
        call tsorder(npteth,zteth,dnethom,tethom,sgneth,sgteth)
        if (sgnemin.lt.0.0) sgnemin=abs(sgnemin)*dnethom(1)*1.e-19 &
                                    *co2cor
        do 40020 i=1,npneth
          dnethom(i)=dnethom(i)*1.e-19*co2cor
          sgneth(i)=sgneth(i)*1.e-19*sgnethi*co2cor
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
!-----------------------------------------------------------------------
!---- read in limiter data                                            --
!-----------------------------------------------------------------------
      call getlim(1,xltype,xltype_180)
!
  212 continue
      if (iconvr.ge.0) go to 214
      iecurr=1
      ivesel=1
  214 continue
      if (iand(iout,1).ne.0) then
      write (nout,in1)
      write (nout,inwant)
      write (nout,ink)
      write (nout,ins)
      if (kwaitmse.ne.0) write (neqdsk,ina)
      write (nout,inece)
      write (nout,insxr)
      write (nout,invt)
      write (nout,inlibim)
      endif
      if (islve.gt.0) nbdry=40
  216 continue
!
      diamag(jtime)=1.0e-03*dflux
      sigdia(jtime)=1.0e-03*abs(sigdlc)
      pbinj(jtime)=pnbeam
      do 220 i=1,nsilop
        silopt(jtime,i)=coils(i)
  220 continue
      do 230 i=1,nfcoil
        fccurt(jtime,i)=brsp(i)
  230 continue
      do 240 i=1,nmtark
        tangam(jtime,i)=tgamma(i)
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
      if (ecurrt(3).le.-1.e10) ecurrt(3)=ecurrt(1)
      if (ecurrt(5).le.-1.e10) ecurrt(5)=ecurrt(1)
      if (ecurrt(4).le.-1.e10) ecurrt(4)=ecurrt(2)
      if (ecurrt(6).le.-1.e10) ecurrt(6)=ecurrt(2)
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
!-----------------------------------------------------------------------
!-- DIII-D shot > 112000 use a new set of table for magnetic probes   --
!--        shot > 156000 new 2014 set                                 --
!-----------------------------------------------------------------------
      if (kdata.ne.2) then
      if (ishot.ge.112000.and.jtime.le.1) then
        if (ishot.lt.156000) then
          table_dir = table_dir(1:ltbdir)//'112000/'
        else
          table_dir = table_dir(1:ltbdir)//'156014/'
        endif
        ltbdir=ltbdir+7
        ltbdi2=ltbdir
        table_di2 = table_dir
        if (rank == 0) then
          write(*,*) 'table_dir = <',table_dir(1:ltbdir),'>'
        endif
!----------------------------------------------------------------------
!-- read in the f coil response functions                            --
!----------------------------------------------------------------------
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
      endif
      endif
!
      if (pasmat(jtime).le.-1.e3) then
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
  292 continue
!
      fwtref=fwtsi(iabs(nslref))
      if (kersil.eq.2) go to 325
      if (kersil.eq.3) go to 325
      do 300 m=1,nsilop
        tdata1=serror*abs(silopt(jtime,m))
        tdata2=abs(psibit(m))*vbit
        tdata=max(tdata1,tdata2)
        sigmafl0(m)=tdata
        if (tdata.gt.1.0e-10) fwtsi(m)=fwtsi(m)/tdata**nsq
        if (tdata.le.1.0e-10) fwtsi(m)=0.0
  300 continue
!----------------------------------------------------------------------
!--  signal at psi loop # NSLREF is used as reference                --
!----------------------------------------------------------------------
      if (abs(psibit(iabs(nslref))).gt.1.0e-10) go to 340
      coilmx=abs(silopt(jtime,1))
      do 320 i=2,nsilop
        abcoil=abs(silopt(jtime,i))
        coilmx=max(abcoil,coilmx)
  320 continue
      fwtsi(iabs(nslref))=1.0/coilmx**nsq/serror**nsq*fwtref
      go to 340
!
  325 continue
      if (kdata.ne.2) then
      if (jtime.le.1) then
        open(unit=80,status='old',file=table_di2(1:ltbdi2)//'dprobe.dat')
        rsi(1)=-1.
        read (80,in3)
        read (80,10200) (rf(i),zf(i),wf(i),hf(i),af(i),af2(i), &
                i=1,mfcoil)
        if (rsi(1).lt.0.) &
        read (80,10200) (rsi(i),zsi(i),wsi(i),hsi(i),as(i),as2(i), &
                i=1,nsilop)
        read (80,10220) (re(i),ze(i),we(i),he(i),ecid(i), &
                                        i=1,necoil)
        if (ifitvs.gt.0.or.icutfp.eq.2) then
          read (80,10200) (rvs(i),zvs(i),wvs(i),hvs(i), &
                                        avs(i),avs2(i),i=1,nvesel)
        endif
        close(unit=80)
      endif
      endif
!-----------------------------------------------------------------------
!--  Fourier expansion of vessel sgments                              --
!-----------------------------------------------------------------------
      if (ifitvs.gt.0. .and. nfourier.gt.1) then
      do i=1,nvesel
      if(rvs(i).ge.1.75.and.zvs(i).ge.0.) &
      thetav(i)=dasin(zvs(i)/sqrt ((rvs(i)-1.75)**2+(zvs(i))**2))
      if(rvs(i).lt.1.75.and.zvs(i).ge.0.) &
      thetav(i)=pi-dasin(zvs(i)/sqrt ((rvs(i)-1.75)**2+(zvs(i))**2))
      if(rvs(i).lt.1.75.and.zvs(i).lt.0.) &
      thetav(i)=pi-dasin(zvs(i)/sqrt ((rvs(i)-1.75)**2+(zvs(i))**2))
      if(rvs(i).ge.1.75.and.zvs(i).lt.0.) &
      thetav(i)=2*pi+dasin(zvs(i)/sqrt ((rvs(i)-1.75)**2+(zvs(i))**2))
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
!     if (brsptu(1).le.-1.e-20) &
!        brsp(1:nfcoil)=brsptu(1:nfcoil)*turnfc(1:nfcoil)
      if (brsptu(1).gt.-1.e-20) &
         brsp(1:nfcoil)=brsptu(1:nfcoil)*turnfc(1:nfcoil)
      reflux=silopt(jtime,iabs(nslref))
      do 330 m=1,nsilop
        tdata1=errsil*abs(silopt(jtime,m)-reflux)
        tdata2=sicont*rsi(m)*abs(pasmat(jtime))
        tdata=max(tdata1,tdata2)
        tdata2=abs(psibit(m))*vbit
        tdata=max(tdata,tdata2)
        sigmafl0(m)=tdata
        if (tdata.gt.1.0e-10) fwtsi(m)=fwtsi(m)/tdata**nsq
        if (tdata.le.1.0e-10) fwtsi(m)=0.0
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
        if (tdata.gt.1.0e-10) fwtsi(m)=fwtref/tdata**nsq
        if (tdata.le.1.0e-10) fwtsi(m)=0.0
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
        if (tdata.gt.1.0e-10) fwtmp2(m)=fwtmp2(m)/tdata**nsq
        if (tdata.le.1.0e-10) fwtmp2(m)=0.0
  350 continue
      do 400 m=1,nstark
        tdata=abs(siggam(jtime,m))
        if (tdata.gt.1.0e-10) fwtgam(m)=fwtgam(m)/tdata**nsq
        if (tdata.le.1.0e-10) fwtgam(m)=0.0
  400 continue
      do 402 m=1,nfcoil
        tdata1=serror*abs(fccurt(jtime,m))
        tdata2=abs(bitfc(m))*vbit
        tdata=max(tdata1,tdata2)
        sigmaf0(m)=tdata
        if (tdata.gt.1.0e-10) fwtfc(m)=fwtfc(m)/tdata**nsq
        if (tdata.le.1.0e-10) fwtfc(m)=0.0
  402 continue
      if (iecurr.eq.2) then
      do m=1,nesum
        tdata1=serror*abs(ecurrt(m))
        tdata2=abs(bitec(m))*vbit
        tdata=max(tdata1,tdata2)
        sigmae0(m)=tdata
        if (tdata.gt.1.0e-10) fwtec(m)=fwtec(m)/tdata**nsq
        if (tdata.le.1.0e-10) fwtec(m)=0.0
      enddo
      endif
      tdata1=serror*abs(pasmat(jtime))
      tdata2=abs(bitip)*vbit
      tdata=max(tdata1,tdata2)
      sigmaip0=tdata
      if (tdata.gt.1.0e-10) fwtcur=fwtcur/tdata**nsq
      if (tdata.le.1.0e-10) fwtcur=0.0
!----------------------------------------------------------------------
!-- diamagetic flux                                                  --
!----------------------------------------------------------------------
      tdata=abs(sigdia(jtime))
      if (tdata.gt.1.0e-10) fwtdlc=fwtdlc/tdata**nsq
!
      if (sidif.le.-1.0e+10) sidif=tmu*pasmat(jtime)*rcentr/2.0
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
      if (fwtqa.gt.1.0e-03) nqaxis=1
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
      timeus=timeus+0.4
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
      if (kdata.ne.2) then
      if ((iecurr.le.0).or.(idodo.gt.0)) go to 520
      open(unit=nrsppc,status='old',form='unformatted', &
           file=table_dir(1:ltbdir)//'re'//trim(ch1)//trim(ch2)//'.ddd')
      read (nrsppc) rsilec
      read (nrsppc) rmp2ec
      read (nrsppc) gridec
      close(unit=nrsppc)
      idodo=1
  520 continue
!
      if ((ivesel.le.0).or.(idovs.gt.0)) go to 525
      open(unit=nrsppc,status='old',form='unformatted', &
           file=table_dir(1:ltbdir)//'rv'//trim(ch1)//trim(ch2)//'.ddd')
      read (nrsppc) rsilvs
      read (nrsppc) rmp2vs
      read (nrsppc) gridvs
      close(unit=nrsppc)
      idovs=1
      if (ivesel.le.10) go to  525
      open(unit=nffile,status='old',form='unformatted', &
           file=table_dir(1:ltbdir)//'fc'//trim(ch1)//trim(ch2)//'.ddd')
      read (nffile) rfcfc
      close(unit=nffile)
      endif
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
             file=table_dir(1:ltbdir)//'ra'//trim(ch1)//trim(ch2)//'.ddd')
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
             file=table_dir(1:ltbdir)//'ef'//trim(ch1)//trim(ch2)//'.ddd')
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
        sigrid(i)=1./float(nw-1)*(i-1)
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
        th=twopi*(k-1)/float(mx)
        rmx(k)=relip-aelip*cos(th)
        zmx(k)=zelip+eelip*aelip*sin(th)
        ix=1
        if (k.gt.(mx/2+1)) ix=2
        i=(rmx(k)-rgrid(1))/drgrid+1
        j=(zmx(k)-zgrid(1))/dzgrid+ix
        zdif=zmx(k)-zgrid(j)
        if (abs(zdif).gt.0.6*dzgrid) then
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
               file=table_dir(1:ltbdir)//'fc'//trim(ch1)//trim(ch2)//'.ddd')
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
        xnpc=float(npc)
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
!--  Solove equilibrium                                                      --
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
        dth=twopi/float(nbdry)
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
        saaa=0.50
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
        rnow=0.5*(rgrids(1)+rgrids(nw))
        znow=0.0
        call surfac(siwant,psi,nw,nh,rgrids,zgrids,xout,yout,nfound, &
                    npoint,drgrids,dzgrids,xmin,xmax,ymin,ymax,npack, &
                    rnow,znow,negcur)
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
          if (rsum.gt.1.0e-08) go to 72054
          mk=(i-1)*nh+1
          rbdrpc(m,k)=gridpc(mk,i)
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
              if (rsum.le.1.0e-08) then
                mk=(i-1)*nh+1
                rsolpc(m,k)=gridpc(mk,i)
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
        if ((delx2+dely2).le.1.0e-08) nbdry=nbdry-1
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
 1250 continue
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
 4980 format (i5)
 5000 format (2e12.6)
 5020 format (1x,4e16.9)
 6550 format (/,' ir = ',10i4)
 6555 format (' iz = ',10i4)
 6557 format (/,' npc = ',i4)
10000 continue
      open(unit=40,file='errfil.out',status='unknown',access='append' &
                                 )
      write (40,10020) ishot,itime
      write (nttyo,10020) ishot,itime
      close(unit=40)
10020 format (///,1x,4hshot,i6,4h at ,i6,4h ms ,'%% Bad Input Data %%')
! MPI >>>
#if defined(USEMPI)
      call mpi_stop
#else
      stop
#endif
! MPI <<<
10140 format (i7,2x,i4,4(1x,f7.3),9(1x,f8.3))
10200 format (6e12.6)
10220 format (5e10.4)
      end
      subroutine dslant(x,y,np,xmin,xmax,ymin,ymax,x1,y1,x2,y2,dismin)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          dslant finds the minimum distance between a curve       **
!**          represented by (x,y) and the line given by (x1,y1)      **
!**          and (x2,y2).                                            **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          12/08/86..........first created                         **
!**                                                                  **
!**********************************************************************
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      data nn/30/
      dimension x(1),y(1)
!
      dismin=1.0e+20
      delx=x2-x1
      dely=y2-y1
      dels=sqrt(delx**2+dely**2)
      nn=dels/0.002
      nn=max(5,nn)
      delx=delx/float(nn-1)
      dely=dely/float(nn-1)
      do 1000 i=1,nn
        xw=x1+delx *(i-1)
        yw=y1+dely *(i-1)
        do 900 m=1,np
          if (x(m).lt.xmin) go to 900
          if (x(m).gt.xmax) go to 900
          if (y(m).lt.ymin) go to 900
          if (y(m).gt.ymax) go to 900
            disw=sqrt((xw-x(m))**2+(yw-y(m))**2)
            dismin=min(dismin,disw)
  900 continue
 1000 continue
      dismin=dismin*100.
      return
      end
      function erpote(ypsi,nnn)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          erpote computes the stream function for the             **
!**          radial electric field. eradial computes the             **
!**          radial electric field.                                  **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          97/04/24..........first created                         **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
       common/cwork3/lkx,lky
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
      end
      function eradial(ypsi,nnn,reee,zeee)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          eradial computes the radial electric field.             **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          97/04/24..........first created                         **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
      common/cwork3/lkx,lky
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
      end
      subroutine fcurrt(jtime,iter,itertt,kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          fcurrt computes the currents in the f coils.            **
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
!**          16/08/90..........revised                               **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: rfcpc
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'

      parameter(nfnves=nfcoil+nvesel)
      dimension afma(nfcoil,nfcoil),ifmatr(nfcoil),wfmatr(nfcoil) &
           ,wbry(msbdry),work(msbdr2),ut(msbdry,msbdry)
      dimension abry(msbdry,nfnves),bbry(msbdry),ainbry(nfnves,msbdry)
      dimension fcref(nfcoil)
      dimension pbry(msbdry),psbry(msbdry)
      data iskip/0/,idoit/0/
! MPI >>>
      integer, intent(inout) :: kerror
      kerror = 0
! MPI <<<
!
      if (ifcurr.gt.0) return
      if (itertt.le.1.and.icinit.lt.0) return
      if (islpfc.ne.1) go to 2000
!-----------------------------------------------------------------------
!--  flux loop on F coils, seldom used for DIII-D                     --
!-----------------------------------------------------------------------
      if (iskip.gt.0) go to 1200
      iskip=1
!vas
!vas      print*,'file name : ','fc'//trim(ch1)// &  
!vas                         trim(ch2)//'.ddd' 
      open(unit=nffile,status='old',form='unformatted', &
		   file=table_dir(1:ltbdir)//'fc'//trim(ch1)//trim(ch2)//'.ddd')
	      read (nffile) rfcfc
	      read (nffile) rfcpc
	      close(unit=nffile)
	!vas
	!vas      print*,'file name : ','fm'//trim(ch1)// &  
	!vas                         trim(ch2)//'.ddd' 
	      open(unit=nffile,status='old',form='unformatted', &
		   file=table_dir(1:ltbdir)//'fm'//trim(ch1)//trim(ch2)//'.ddd')
      read (nffile,err=1150) afma
      read (nffile) ifmatr
      close(unit=nffile)
      go to 1200
 1150 continue
      do 1170 i=1,nfcoil
      do 1170 j=1,nfcoil
        afma(i,j)=rfcfc(i,j)
 1170 continue
      m111=-1.0
      call decomp(nfcoil,nfcoil,afma,m111,ifmatr,wfmatr)
!vas
!vas      print*,'file name : ','fm'//trim(ch1)// &  
!vas                         trim(ch2)//'.ddd' 
      open(unit=nffile,status='old',form='unformatted', &
           file='fm'//trim(ch1)//trim(ch2)//'.ddd' & 
              ,err=12927)
      close(unit=nffile,status='delete')
12927 continue
!vas
!vas      print*,'file name : ','fm'//trim(ch1)// & 
!vas                          trim(ch2)//'.ddd'
      open(unit=nffile,status='new',form='unformatted', &
           file='fm'//trim(ch1)//trim(ch2)//'.ddd')
      write (nffile) afma
      write (nffile) ifmatr
      close(unit=nffile)
 1200 continue
      do 1320 i=1,nfcoil
        brsp(i)=0.0
        if (ivacum.gt.0) go to 1310
        do 1308 j=1,nwnh
        brsp(i)=brsp(i)+rfcpc(i,j)*pcurrt(j)
 1308   continue
 1310   continue
        if (iecurr.le.0) go to 1314
        do 1312 j=1,nesum
          brsp(i)=brsp(i)+rfcec(i,j)*ecurrt(j)
 1312   continue
 1314   continue
        if (ivesel.le.0.or.ifitvs.gt.0) go to 1318
        do 1316 j=1,nvesel
          brsp(i)=brsp(i)+rfcvs(i,j)*vcurrt(j)
 1316   continue
 1318   continue
        brsp(i)=csilop(i,jtime)-brsp(i)
        if (fitsiref) brsp(i)=brsp(i)+psiref(jtime)
 1320 continue
      call solve(nfcoil,nfcoil,afma,brsp,ifmatr)
      return
!-----------------------------------------------------------------------
!--  standard option, flux loops away from F coils                    --
!-----------------------------------------------------------------------
 2000 continue
      if (idoit.gt.0.and.itertt.gt.1) go to 3000
!-----------------------------------------------------------------------
!--  set up response matrix once for all                              --
!-----------------------------------------------------------------------
      idoit=1
      wsibry=psibry
!-----------------------------------------------------------------------
!--   get fixed boundary response from plasma                         --
!-----------------------------------------------------------------------
      if (nbdry.le.0.or.iconvr.ne.3) go to 2098
!-----------------------------------------------------------------------
!-- set up boundary fitting weights                                   --
!-----------------------------------------------------------------------
      z04=1.0e-04
      fwtbdr=abs(errbry)*max(abs(sidif),z04)
      fwtbdr=1.0/fwtbdr
      do 2080 i=1,nbdry
        fwtbry(i)=fwtbdr*fwtbdry(i)
 2080 continue
 2098 continue
!-----------------------------------------------------------------------
!--  set up response matrix, first fixed boundary, then flux loops    --
!--  and F coil currents                                              --
!-----------------------------------------------------------------------
      do 2500 nk=1,nfcoil
      nj=0
      if (nbdry.le.0.or.iconvr.ne.3) go to 2200
      do 2100 m=1,nbdry
        abry(m,nk)=fwtbry(m)*rbdrfc(m,nk)
 2100 continue
      nj=nbdry
 2200 continue
!-----------------------------------------------------------------------
!--  flux loop option                                                 --
!-----------------------------------------------------------------------
      do 2220 m=1,nsilop
        if (fwtsi(m).le.0.0) go to 2220
        nj=nj+1
        abry(nj,nk)=fwtsi(m)*rsilfc(m,nk)
 2220 continue
!-----------------------------------------------------------------------
!--  magnetic probes used only for vacuum error field analysis        --
!-----------------------------------------------------------------------
      if (ivacum.gt.0) then
      do 2222 m=1,magpri
        if (fwtmp2(m).le.0.0) go to 2222
        nj=nj+1
        abry(nj,nk)=fwtmp2(m)*rmp2fc(m,nk)
 2222 continue
      endif
!-----------------------------------------------------------------------
!--  F coil currents option                                           --
!-----------------------------------------------------------------------
      do 2225 m=1,nfcoil
        if (fwtfc(m).le.0.0) go to 2225
        nj=nj+1
        abry(nj,nk)=0.
        if (m.eq.nk) abry(nj,nk)=fwtfc(m)
 2225 continue
!-----------------------------------------------------------------------
!--  minimize coil current oscillations for fixed boundary option     --
!-----------------------------------------------------------------------
      if (nbdry.le.0.or.iconvr.ne.3) go to 2250
      do 2240 m=1,nfcoil
        nj=nj+1
        abry(nj,nk)=0.0
        if (nk.eq.m) abry(nj,nk)=cfcoil*fczero(m)
 2240 continue
!-----------------------------------------------------------------------
!--  constraint F coil currents to be symmetric if requested          --
!-----------------------------------------------------------------------
      if ((symmetrize).and.(cupdown.ne.0.0)) then
      do m=1,nfcoil/2
        nj=nj+1
        abry(nj,nk)=0.0
        if (nk.eq.m) abry(nj,nk)=cupdown
        if (nk.eq.(m+nfcoil/2)) abry(nj,nk)=-cupdown
      enddo
      endif
 2250 continue
!-----------------------------------------------------------------------
!--  F coil constraints Cf=x                                          --
!-----------------------------------------------------------------------
      if (kccoils.gt.0) then
      do m=1,kccoils
        nj=nj+1
        abry(nj,nk)=ccoils(nk,m)
      enddo
      endif
!-----------------------------------------------------------------------
!--  SOL constraints 2014/05/20 LL                                    --
!-----------------------------------------------------------------------
      if (nsol.gt.0) then
!-----------------------------------------------------------------------
!-- set up SOL fitting weights and response matrix                    --
!-----------------------------------------------------------------------
        z04s=1.0e-04
        fwtsols=abs(errbry)*max(abs(sidif),z04s)
        fwtsols=1.0/fwtsols
        do i=1,nsol
          fwtsolw(i)=fwtsols*fwtsol(i)
          nj=nj+1
          abry(nj,nk)=fwtsolw(i)*rsolfc(i,nk)
        enddo
        if (idebug >= 2) write (6,*) 'FCURRT nsol,fwtsolw = ', nsol,(fwtsolw(i),i=1,nsol)
      endif
 2500 continue
!-----------------------------------------------------------------------
!--  optional vessel current model                                    --
!-----------------------------------------------------------------------
      neqn=nfcoil
      if (ifitvs.le.0) go to 2610
! 
      if (nfourier.gt.1) then
        nuuu=nfourier*2+1
      else
        nuuu=nvesel
      endif 
      do 2600 nkk=1,nuuu
      nk=nkk+nfcoil
      nj=0
      if (nbdry.le.0.or.iconvr.ne.3) go to 2510
      do 2505 m=1,nbdry
          if(nfourier.gt.1)then
          temp=0.
          do i=1,nvesel
          temp=temp+rbdrvs(m,i)*vecta(nkk,i)
          enddo
            abry(m,nk)=fwtbry(m)*temp
          else
            abry(m,nk)=fwtbry(m)*rbdrvs(m,nkk)
          endif
 2505 continue
      nj=nbdry
 2510 continue
      do 2520 m=1,nsilop
        if (fwtsi(m).le.0.0) go to 2520
        nj=nj+1
        if(nfourier.gt.1)then
          temp=0.
          do i=1,nvesel
          temp=temp+rsilvs(m,i)*vecta(nkk,i)
          enddo
          abry(nj,nk)=fwtsi(m)*temp
        else
           abry(nj,nk)=fwtsi(m)*rsilvs(m,nkk)
        endif
 2520 continue
!-----------------------------------------------------------------------
!-- magnetic probes used only for vacuum error analysis               --
!-----------------------------------------------------------------------
      if (ivacum.gt.0) then
      do 2522 m=1,magpri
        if (fwtmp2(m).le.0.0) go to 2522
        nj=nj+1
        if(nfourier.gt.1)then
          temp=0.
          do i=1,nvesel
          temp=temp+rmp2vs(m,i)*vecta(nkk,i)
          enddo
           abry(nj,nk)=fwtmp2(m)*temp
        else
           abry(nj,nk)=fwtmp2(m)*rmp2vs(m,nkk)
        endif
 2522 continue
      endif
      do 2525 m=1,nfcoil
        if (fwtfc(m).le.0.0) go to 2525
        nj=nj+1
        abry(nj,nk)=0.
 2525 continue
      if (nbdry.le.0.or.iconvr.ne.3) go to 2550
      do 2540 m=1,nfcoil
        nj=nj+1
        abry(nj,nk)=0.0
 2540 continue
 2550 continue
 2600 continue
      neqn=neqn+nvesel
 2610 continue
!-----------------------------------------------------------------------
!--  set up working matrix                                            --
!-----------------------------------------------------------------------
      do 2800 i=1,msbdry
      do 2800 j=1,msbdry
        ut(i,j)=0.0
 2800 continue
      do 2870 i=1,msbdry
 2870 ut(i,i) = 1.0
!-----------------------------------------------------------------------
!--  compute inverse of abry in least squares sence with singular     --
!--  value decomposition                                              --
!-----------------------------------------------------------------------
      if (idebug >= 2) write (6,*) 'FCURRT nj,neqn = ', nj,neqn
      call sdecm(abry,msbdry,nj,neqn,ut,msbdry,nj,wbry,work,ier)
      if (ier.ne.129) go to 2880
      write (nttyo,8000) ier
! MPI >>>
#if defined(USEMPI)
      call mpi_stop
#else
      stop
#endif
! MPI <<<
 2880 continue
      cond=ier
      toler=1.0e-06*wbry(1)
      do 2890 i=1,neqn
      do 2890 j=1,nj
        ainbry(i,j)=0.0
        do 2890 k=1,neqn
          t=0.0
          if (wbry(k).gt.toler) t=1.0/wbry(k)
          ainbry(i,j) = ainbry(i,j) + abry(i,k)*ut(k,j)*t
 2890   continue
!-----------------------------------------------------------------------
!--  compute and store coil reference currents, not need for          --
!--  non-fixed boundary when IFREF=-1                                 --
!-----------------------------------------------------------------------
      if (ifref.le.0.or.iconvr.ne.3) goto 3000
      do 2900 i=1,nj
        if (i.le.nbdry) then
          bbry(i)=fwtbry(i)
        else
          bbry(i)=0.0
        endif
 2900 continue 
      do 2952 i=1,nfcoil
        fcref(i)=0.0
        do 2952 j=1,nj
          fcref(i)=fcref(i)+ainbry(i,j)*bbry(j)
 2952   continue
!-----------------------------------------------------------------------
!--  RHS of response matrix, start here if got inverse matrix already --
!-----------------------------------------------------------------------
 3000 continue
      nj=0
      if (nbdry.le.0.or.iconvr.ne.3) go to 3600
!-----------------------------------------------------------------------
!-- fixed boundary portion                                            --
!-----------------------------------------------------------------------
      do 3596 m=1,nbdry
        bbry(m)=0.0
        if (ivesel.le.0.or.ifitvs.gt.0) go to 3570
        do 3565 k=1,nvesel
          bbry(m)=bbry(m)+rbdrvs(m,k)*vcurrt(k)
 3565   continue
 3570   continue
        if (iecurr.le.0) go to 3590
        do 3580 k=1,nesum
          bbry(m)=bbry(m)+rbdrec(m,k)*ecurrt(k)
 3580   continue
 3590   continue
        if (iacoil.gt.0) then
          do 3592 k=1,nacoil
            bbry(m)=bbry(m)+rbdrac(m,k)*caccurt(jtime,k)
 3592     continue
        endif
        do 3594 k=1,nwnh
         bbry(m)=bbry(m)+rbdrpc(m,k)*pcurrt(k)
 3594   continue
        pbry(m)=bbry(m)
        bbry(m)=fwtbry(m)*(wsibry-bbry(m))
 3596 continue
      nj=nbdry
!-----------------------------------------------------------------------
!--  flux loops portion                                               --
!-----------------------------------------------------------------------
 3600 continue
      do 3630 m=1,nsilop
        if (fwtsi(m).le.0.0) go to 3630
        nj=nj+1
        bbry(nj)=0.0
        do 3610 j=1,nwnh
        bbry(nj)=bbry(nj)+gsilpc(m,j)*pcurrt(j)
 3610   continue
        if (iecurr.le.0) go to 3614
        do 3612 j=1,nesum
          bbry(nj)=bbry(nj)+rsilec(m,j)*ecurrt(j)
 3612   continue
 3614   continue
!-----------------------------------------------------------------------
!-- specify vessel currents ?                                         --
!-----------------------------------------------------------------------
        if (ivesel.le.0.or.ifitvs.gt.0) go to 3618
        do 3616 j=1,nvesel
          bbry(nj)=bbry(nj)+rsilvs(m,j)*vcurrt(j)
 3616   continue
 3618   continue
!-----------------------------------------------------------------------
!-- use advance divertor coil ?                                       --
!-----------------------------------------------------------------------
        if (iacoil.gt.0) then
          do 3621 j=1,nacoil
            bbry(nj)=bbry(nj)+rsilac(m,j)*caccurt(jtime,j)
 3621     continue
        endif
        if (fitsiref) then
        bbry(nj)=fwtsi(m)*(silopt(jtime,m)+psiref(jtime)-bbry(nj))
        else
        bbry(nj)=fwtsi(m)*(silopt(jtime,m)-bbry(nj))
        endif
 3630 continue
!-----------------------------------------------------------------------
!--  magnetic probes used only for vacuum error field analysis        --
!-----------------------------------------------------------------------
      if (ivacum.gt.0) then
      do 3650 m=1,magpri
        if (fwtmp2(m).le.0.0) go to 3650
        nj=nj+1
        bbry(nj)=0.0
        do 3640 j=1,nwnh
        bbry(nj)=bbry(nj)+gmp2pc(m,j)*pcurrt(j)
 3640   continue
        if (iecurr.le.0) go to 3644
        do 3642 j=1,nesum
          bbry(nj)=bbry(nj)+rmp2ec(m,j)*ecurrt(j)
 3642   continue
 3644   continue
!-----------------------------------------------------------------------
!-- specify vessel currents ?                                         --
!-----------------------------------------------------------------------
        if (ivesel.le.0.or.ifitvs.gt.0) go to 3648
        do 3646 j=1,nvesel
          bbry(nj)=bbry(nj)+rmp2vs(m,j)*vcurrt(j)
 3646   continue
 3648   continue
        if (iacoil.gt.0) then
        do 3649 j=1,nacoil
          bbry(nj)=bbry(nj)+rmp2ac(m,j)*caccurt(jtime,j)
 3649   continue
        endif
        bbry(nj)=fwtmp2(m)*(expmpi(jtime,m)-bbry(nj))
 3650 continue
      endif
!-----------------------------------------------------------------------
!-- F-coil currents specification                                     --
!-----------------------------------------------------------------------
      do 3660 m=1,nfcoil
        if (fwtfc(m).le.0.0) go to 3660
        nj=nj+1
        bbry(nj)=fccurt(jtime,m)*fwtfc(m)
 3660 continue
!-----------------------------------------------------------------------
!--  F coil current minimization                                      --
!-----------------------------------------------------------------------
      if (nbdry.le.0.or.iconvr.ne.3) go to 3665
      do 3662 m=1,nfcoil
        nj=nj+1
        bbry(nj)=0.0
 3662 continue
!-----------------------------------------------------------------------
!--  symmetrize F coil currents ?                                     --
!-----------------------------------------------------------------------
      if ((symmetrize).and.(cupdown.ne.0.0)) then
      do m=1,nfcoil/2
        nj=nj+1
        bbry(nj)=0.0
      enddo
      endif
 3665 continue
!-----------------------------------------------------------------------
!--  F coil constraints ?                                             --
!-----------------------------------------------------------------------
      if (kccoils.gt.0) then
        do m=1,kccoils
          nj=nj+1
          bbry(nj)=xcoils(m)
        enddo
      endif
!-----------------------------------------------------------------------
!-- SOL Constraints RHS                                               --
!-----------------------------------------------------------------------
      if (nsol.gt.0) then
        if (idebug >= 2) write (6,*) 'FCURRT wsisol = ', wsisol
        do m=1,nsol
          nj=nj+1
          bbry(nj)=0.0
          if (ivesel.gt.0.and.ifitvs.lt.0) then
            do k=1,nvesel
              bbry(nj)=bbry(nj)+rsolvs(m,k)*vcurrt(k)
            enddo
          endif
          if (iecurr.gt.0) then
             do k=1,nesum
                bbry(nj)=bbry(nj)+rsolec(m,k)*ecurrt(k)
             enddo
           endif
           do k=1,nwnh
             bbry(nj)=bbry(nj)+rsolpc(m,k)*pcurrt(k)
           enddo
           psbry(m)=bbry(nj)
           bbry(nj)=fwtsolw(m)*(wsisol-bbry(nj))
        enddo
      endif
!-----------------------------------------------------------------------
!--  now get F coil currents from precomputed inverse matrix          --
!-----------------------------------------------------------------------
      do 3670 i=1,nfcoil
        brsp(i)=0.0
        do 3670 j=1,nj
          brsp(i)=brsp(i)+ainbry(i,j)*bbry(j)
 3670 continue
      if (ifitvs.le.0) go to 3685
      do 3680 ii=1,nvesel
        i=ii+nfcoil
        vcurrt(ii)=0.0
        do 3680 j=1,nj
          vcurrt(ii)=vcurrt(ii)+ainbry(i,j)*bbry(j)
 3680 continue
!
       if(nfourier.gt.1) then
        do j=1,nvesel
        temp=0.
        do 2671 i=1,(nfourier*2+1)
          temp=temp+vcurrt(i)*vecta(i,j)
 2671   continue
        vcurrt(j)=temp
        enddo
       endif
!
 3685 continue
!-----------------------------------------------------------------------
!--  adjustments for various boundary flux options                    --
!-----------------------------------------------------------------------
      if (ifref.le.0.or.iconvr.ne.3) go to 4000
      sumif = 0.0
      sumifr = 0.0
!-----------------------------------------------------------------------
!-- sum of inner F-coils 1-5 A and B zero IFREF=1                     --
!-----------------------------------------------------------------------
      if (ifref.ne.1) go to 3700
      do 3690 i=1,5
 3690 sumif = sumif + brsp(i) + brsp(i+9)
      sumif = sumif + brsp(8) + brsp(17)
      do 3695 i=1,5
 3695 sumifr = sumifr + fcref(i) + fcref(i+9)
      sumifr = sumifr + fcref(8) + fcref(17)
      go to 3750
 3700 continue
!-----------------------------------------------------------------------
!--  sum of F coils selected through FCSUM vanish IFREF=2             --
!-----------------------------------------------------------------------
      if (ifref.ne.2) go to 3720
      do 3705 i=1,nfcoil
 3705 sumif=sumif+fcsum(i)*brsp(i)/turnfc(i)
      do 3710 i=1,nfcoil
 3710 sumifr = sumifr + fcref(i) * fcsum(i)/turnfc(i)
      go to 3750
 3720 continue
!----------------------------------------------------------------------
!--  choose boundary flux by minimize coil currents IFREF=3          --
!----------------------------------------------------------------------
      if (ifref.eq.3) then
      do 3725 i=1,nfcoil
 3725 sumif=sumif+fcref(i)*brsp(i)*fczero(i)
      do 3730 i=1,nfcoil
 3730 sumifr = sumifr + fcref(i)**2*fczero(i)
      go to 3750
      endif
 3750 continue
!-----------------------------------------------------------------------
!--  update boundary flux for IFREF=1-3                               --
!-----------------------------------------------------------------------
      if (ifref.le.3) then
      ssiref = sumif/sumifr
      do 3759 m=1,nfcoil
      silopt(jtime,m)=silopt(jtime,m)-ssiref
 3759 brsp(m) = brsp(m) - ssiref*fcref(m)
      wsibry=wsibry-ssiref
      wsisol=wsisol-ssiref
      endif
!------------------------------------------------------------------------
!--  fixed boundary flux specified through PSIBRY, IFREF=4             --
!------------------------------------------------------------------------
      if (ifref.eq.4) then
      ssiref=psibry0-psibry
      do 3770 m=1,nfcoil
        brsp(m)=brsp(m)+ssiref*fcref(m)
        silopt(jtime,m)=silopt(jtime,m)+ssiref
 3770 continue
      wsibry=wsibry+ssiref
      wsisol=wsisol+ssiref
      endif
!-----------------------------------------------------------------------
!--  done, estimated errors for fixed boundary calculations           --
!-----------------------------------------------------------------------
 4000 continue
      if (nbdry.gt.0.and.iconvr.eq.3) then
      erbmax=0.0
      erbave=0.0
      do 4100 i=1,nbdry
        xsibry=pbry(i)
        do 4050 m=1,nfcoil
          xsibry=xsibry+rbdrfc(i,m)*brsp(m)
 4050   continue
        erbloc(i)=abs((wsibry-xsibry)/sidif)
        erbmax=max(erbloc(i),erbmax)
        erbave=erbave+erbloc(i)
 4100 continue
      erbave=erbave/float(nbdry)
      if (idebug /= 0) write (6,*) 'FCURRT erbmax,erbave,si = ', erbmax,erbave,wsibry
      endif
!
      if (nsol.gt.0.and.iconvr.eq.3) then
      erbsmax=0.0
      erbsave=0.0
      do i=1,nsol
        xsisol=psbry(i)
        do m=1,nfcoil
          xsisol=xsisol+rsolfc(i,m)*brsp(m)
        enddo
        erbsloc(i)=abs((wsisol-xsisol)/sidif)
        erbsmax=max(erbsloc(i),erbsmax)
        erbsave=erbsave+erbsloc(i)
      enddo
      erbsave=erbsave/float(nsol)
      if (idebug /= 0) write (6,*) 'FCURRT erbsmax,erbsave,si = ', erbsmax,erbsave,wsisol
      endif
!
 4500 return
 8000 format (/,'  ** Problem in Decomposition **',i10)
      end
      subroutine autoknot(ks,lconvr,ktime,mtear,kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  Autoknot locator for spline basis function    **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          autoknow minimizes chi-squared as a function of knot    **
!**          location                                                **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          20/01/2000........first created                         **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
      common/autok/ks_a,lconvr_a,ktime_a,mtear_a,kerror_a,kadknt, &
                   appknt(npcurn),kappknt, &
                   affknt(npcurn),kaffknt, &
                   awwknt(npcurn),kawwknt, &
                   aeeknt(npcurn),kaeeknt,mxiter_a
      include 'basiscomdu.f90'
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
		goto 10
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
           goto 20
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
              goto 30
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
              goto 40
            endif
         enddo
      endif
!
!     Now do the final fit with adjusted knots and full iterations
!

40    continue
      call data_input(ks_a,lconvr_a,ktime_a,mtear_a)
      if (lconvr_a.lt.0) go to 500
      mxiter_a = saveiter
      call restore_autoknotvals
      call inicur(ks_a)
      call fit(ks_a,kerror_a)
      if (kerror_a /= 0) then
        kerror = 1
      endif
500   return
      end
!
!    store values read from k file into autoknot variables
!
      subroutine restore_autoknotvals
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
      common/autok/ks_a,lconvr_a,ktime_a,mtear_a,kerror_a,kadknt, &
                   appknt(npcurn),kappknt, &
                   affknt(npcurn),kaffknt, &
                   awwknt(npcurn),kawwknt, &
                   aeeknt(npcurn),kaeeknt,mxiter_a
      include 'basiscomdu.f90'

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
      end
!
!     store autoknot variables into standard efit names
!     for example knot locations
!
      subroutine store_autoknotvals
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
      common/autok/ks_a,lconvr_a,ktime_a,mtear_a,kerror_a,kadknt, &
                   appknt(npcurn),kappknt, &
                   affknt(npcurn),kaffknt, &
                   awwknt(npcurn),kawwknt, &
                   aeeknt(npcurn),kaeeknt,mxiter_a
      include 'basiscomdu.f90'

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
      end
!
! used by autoknot which passes the routine to a minimization subroutine
! which calls it to evaulate the function being minimized
!

      function ppakfunc(xknot)          
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
      common/autok/ks_a,lconvr_a,ktime_a,mtear_a,kerror_a,kadknt, &
                   appknt(npcurn),kappknt, &
                   affknt(npcurn),kaffknt, &
                   awwknt(npcurn),kawwknt, &
                   aeeknt(npcurn),kaeeknt,mxiter_a
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      include 'basiscomdu.f90'
      ppakfunc = 1000.0
      write(6,*)
      write(6,*)' trying pp knot ',kadknt,' at location ',xknot, &
                 ' out of ',kappknt,' knots'
        call data_input(ks_a,lconvr_a,ktime_a,mtear_a)
        if (lconvr_a.lt.0) go to 500
        call restore_autoknotvals
        ppknt(kadknt) = xknot
        call inicur(ks_a)
        call fit(ks_a,kerror_a)
        if(kerror_a .gt. 0)goto 500
        ppakfunc = akchiwt * tsaisq(ks_a)  + akerrwt * errorm  &
                  + akgamwt * chigamt + akprewt * chipre
500      return
        end
!
! used by autoknot which passes the routine to a minimization subroutine
! which calls it to evaulate the function being minimized
!
      function ffakfunc(xknot)          
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
      common/autok/ks_a,lconvr_a,ktime_a,mtear_a,kerror_a,kadknt, &
                   appknt(npcurn),kappknt, &
                   affknt(npcurn),kaffknt, &
                   awwknt(npcurn),kawwknt, &
                   aeeknt(npcurn),kaeeknt,mxiter_a
!
      include 'basiscomdu.f90'
      ffakfunc = 1000.0
      write(6,*)
      write(6,*)' trying ff knot ',kadknt,' at location ',xknot, &
                 ' out of ',kaffknt,' knots'
        call data_input(ks_a,lconvr_a,ktime_a,mtear_a)
        if (lconvr_a.lt.0) go to 500
        call restore_autoknotvals
        ffknt(kadknt) = xknot
        call inicur(ks_a)
        call fit(ks_a,kerror_a)
        if(kerror_a .gt. 0)goto 500
        ffakfunc = akchiwt * tsaisq(ks_a)  + akerrwt * errorm  &
                  + akgamwt * chigamt + akprewt * chipre
500      return
        end
!
! used by autoknot which passes the routine to a minimization subroutine
! which calls it to evaulate the function being minimized
!
      function wwakfunc(xknot)          
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
      common/autok/ks_a,lconvr_a,ktime_a,mtear_a,kerror_a,kadknt, &
                   appknt(npcurn),kappknt, &
                   affknt(npcurn),kaffknt, &
                   awwknt(npcurn),kawwknt, &
                   aeeknt(npcurn),kaeeknt,mxiter_a
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      include 'basiscomdu.f90'
      wwakfunc = 1000.0
      write(6,*)
      write(6,*)' trying ww knot ',kadknt,' at location ',xknot, &
                 ' out of ',kawwknt,' knots'
        call data_input(ks_a,lconvr_a,ktime_a,mtear_a)
        if (lconvr_a.lt.0) go to 500
        call restore_autoknotvals
        wwknt(kadknt) = xknot
        call inicur(ks_a)
        call fit(ks_a,kerror_a)
        if(kerror_a .gt. 0)goto 500
        wwakfunc = akchiwt * tsaisq(ks_a)  + akerrwt * errorm &
                  + akgamwt * chigamt + akprewt * chipre
500      return
        end
!
! used by autoknot which passes the routine to a minimization subroutine
! which calls it to evaulate the function being minimized
!
      function eeakfunc(xknot)          
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
      common/autok/ks_a,lconvr_a,ktime_a,mtear_a,kerror_a,kadknt, &
                   appknt(npcurn),kappknt, &
                   affknt(npcurn),kaffknt, &
                   awwknt(npcurn),kawwknt, &
                   aeeknt(npcurn),kaeeknt,mxiter_a
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      include 'basiscomdu.f90'
      eeakfunc = 1000.0
      write(6,*)
      write(6,*)' trying ee knot ',kadknt,' at location ',xknot, &
                 ' out of ',kaeeknt,' knots'
!
        call data_input(ks_a,lconvr_a,ktime_a,mtear_a)
        if (lconvr_a.lt.0) go to 500
        call restore_autoknotvals
        eeknt(kadknt) = xknot
        call inicur(ks_a)
        call fit(ks_a,kerror_a)
        if(kerror_a .gt. 0)goto 500
        eeakfunc = akchiwt * tsaisq(ks_a)  + akerrwt * errorm &
                  + akgamwt * chigamt + akprewt * chipre
500      return
        end
      subroutine fit(jtime,kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          fit carries out the fitting and equilibrium iterations. **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**          94/03/11..........revised                               **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules1.f90'
      include 'basiscomdu.f90'
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
      do 2010 i=1,iend
        if (i.gt.1) iwantk=iwantk+1
        ix=i
        if (i.le.1) go to 500
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
        if ((iconvr.eq.2).and.(ichisq.gt.0)) go to 2020
  500   continue
        do 2000 in=1,nxiter
          ixnn=in
          nitera=nitera+1
          call currnt(ix,jtime,ixnn,nitera,kerror)
          if (kerror /= 0) then
            jerror(jtime) = 1
            return
          endif
          if (ivesel.ge.2) call vescur(jtime)
          if ((i.le.1).or.(in.gt.1)) call fcurrt(jtime,ix,nitera,kerror)
          if (kerror /= 0) then
            jerror(jtime) = 1
            return
          endif
          call pflux(ix,ixnn,nitera,jtime)
          call steps(ixnn,nitera,ix,jtime,kerror)
          if (kerror /= 0) then
            jerror(jtime) = 1
            return
          endif
           if (kwaitmse.ne.0 .and. i.ge.kwaitmse)  &
      	              call fixstark(jtime,kerror)
          call residu(nitera,jtime)
          if ((nitera.lt.kcallece).and.(kfitece.gt.0.0)) go to 2010
          if ((in.eq.1).and.(idone.gt.0)) then
            if (tsaisq(jtime).le.saimin) go to 2020
          endif
          if (idone.gt.0) go to 2010
          if (i.eq.mxiter+1) go to 2010
 2000   continue
 2010 continue
 2020 continue
!---------------------------------------------------------------------
!--  update pressure if needed                                      --
!---------------------------------------------------------------------
      if (kprfit.gt.1) call presur(jtime,nitera,kerror)
      if (kerror /= 0) then
        jerror(jtime) = 1
      endif
      return
      end
      function fpcurr(upsi,nnn)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          fpcurr computes the radial derivative                   **
!**          of the poloidal current ff. ffcurr computes             **
!**          the poloidal current F=twopi RBt/mu0                    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**          94/03/08..........revised                               **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules1.f90'
      dimension xpsii(nffcur)
      real*8, external :: linear
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
      end
      subroutine fluxav(f,x,y,n,si,rx,msx,ry,msy,fave,ns,sdlobp,sdlbp)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
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
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          14/10/87..........first created                         **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparmdud129.f90'
!vas
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      common/cwork3/lkx,lky
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
        xnow=0.5*(x(i-1)+x(i))
        ynow=0.5*(y(i-1)+y(i))
        fnow=0.5*(f(i-1)+f(i))
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
      end
      subroutine getbeam
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          GETBEAM gets the beam pressure.                         **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          15/09/87..........first created                         **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      dimension bwork(ndata),cwork(ndata),dwork(ndata)
!
      if (nbeam.lt.0) go to 2000
!-----------------------------------------------------------------
!--  interpolate beam pressure into Thomson grid                --
!-----------------------------------------------------------------
      call zpline(nbeam,sibeam,pbeam,bwork,cwork,dwork)
      do 3000 i=1,npress
        xn=-rpress(i)
        pbimth(i)=seval(nbeam,xn,sibeam,pbeam,bwork,cwork,dwork)
 3000 continue
      pbeamb=seval(nbeam,x111,sibeam,pbeam,bwork,cwork,dwork)
      pbimpb=speval(nbeam,x111,sibeam,pbeam,bwork,cwork,dwork)
!-----------------------------------------------------------------
!--  interpolate beam ion density into Thomson grid             --
!-----------------------------------------------------------------
      call zpline(nbeam,sibeam,dnbeam,bwork,cwork,dwork)
      do 4000 i=1,npress
        xn=-rpress(i)
        dnbthom(i)=seval(nbeam,xn,sibeam,dnbeam,bwork,cwork,dwork)
 4000 continue
      return
 2000 continue
!------------------------------------------------------------------
!--  compute beam pressure analytically                          --
!------------------------------------------------------------------
      return
      end
      subroutine geteceb(jtime,kerror)
!**********************************************************************
!**                                                                  **
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          geteceb obtains the receo, R+ R-                        **
!**          from ECE measurement data, (fitting T(B))               **
!**          if kfixro kfixrece = -1, called in setece               **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          1/99..........first created, Cheng Zhang                **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!
      parameter (nn=30)
      parameter (kbre=5)
      common/cwork3/lkx,lky

      dimension pds(6),nnout(kbre),bmink(kbre),bmaxk(kbre) 
      real*8,allocatable :: rrgrid(:,:),bfield(:),rrout(:,:), &
          bout(:,:),babs(:,:),bbb(:),ccc(:),ddd(:),btttt(:), &
          dsidr(:),ddsiddr(:),bx(:),ry(:),bbk(:),dbdr(:)
      dimension arspfit(nnecein,nn),brspfit(nnecein) &
           ,s(nn),tte(nnecein),x(nn),an(nnecein) &
           ,tebit(nnecein)
      dimension telowf(nnnte) &
         ,blowf(nnnte),bb(nnnte),cc(nnnte),dd(nnnte) &
         ,teece(nnece),pteprm(nnece),pteprp(nnece) &
         ,idestp(nnece),idestm(nnece),becem(nnece),becep(nnece) &
         ,dbdrp(nnece),dbdrm(nnece)
      integer, intent(inout) :: kerror
      kerror = 0
!-------------------------------------------------------------------
      ALLOCATE(rrgrid(kbre,nw),bfield(nw),rrout(kbre,nw), &
         bout(kbre,nw),babs(kbre,nw),bbb(nw),ccc(nw), &
         ddd(nw),btttt(nw),dsidr(nw),ddsiddr(nw),bx(nw), &
         ry(nw),bbk(nw),dbdr(nw))
!
      telowf=0 &
         ;blowf=0;bb=0;cc=0;dd=0 &
         ;teece=0;pteprm=0;pteprp=0 &
         ;idestp=0;idestm=0;becem=0;becep=0 &
         ;dbdrp=0;dbdrm=0

!-------------------------------------------------------------------
      do k=1,nnece
         fwtece0(k)=swtece(k)
      enddo
      fwtecebz0=swtecebz
      do k=1,kbre
       do i=1,nw
         babs(k,i)=0.0
         bout(k,i)=0.0
         rrout(k,i)=0.0
         rrgrid(k,i)=0.0
       enddo
      enddo
!-------------------------------------------------------------------
!--     input kgeteceb=0 from input file
!-------------------------------------------------------------------
      if (kgeteceb.gt.0) go to 539
      kgeteceb=kgeteceb+1
!---------------------------------------------------------------------
!--  Calculation of |B| array from fe array ( harmonic nharm)       --
!--     becein(necein),   fe(GHz),|B|(T)                            --
!--     !!! becein  from Low field to high field !!!                --
!---------------------------------------------------------------------
      do 500 k=1,necein
      becein(k)=0.001*6.0*9.1095*3.14159/4.8032*feece0(k)/float(nharm)
500   continue
!EALW      write(*,*)'becein'
!EALW      write(*,*)becein
!--------------------------------------------------------------------
!--   fitting data from teecein0,errorece0 and becein (nnecein)    --
!--     bbx=(B-b00)/baa                                            --
!--     Te=x(1)+x(2)*bbx+x(3)*bbx**2+...+x(nfit)*bbx**(nfit-1)     --
!--------------------------------------------------------------------
!heng          mm--nnecein   m---necein  nn--parameter, n--nfit input
      binmin=becein(1)
      binmax=becein(necein)
      baa=0.5*(binmax-binmin)
      b00=0.5*(binmax+binmin)
      do i=1,necein
        an(i)=(becein(i)-b00)/baa
        tebit(i)=max(errorece0(i),dble(1.e-4))
      enddo
      do 1110 nj=1,necein
      do 1100 nk=1,nfit
      if (nk.eq.1) then
        arspfit(nj,nk)=1./tebit(nj)
      else
        arspfit(nj,nk)=an(nj)**(nk-1)/tebit(nj)
      endif
1100   continue
1110   continue
!---------------------------------------------------------------------
!--   teecein0,errorece0  from low field to high field
!---------------------------------------------------------------------
      do 1200 nj=1,necein
      brspfit(nj)=teecein0(nj)/tebit(nj)
1200   continue
!
      mnow=necein
      if (kcmin.gt.0) then
      fwtnow=0.001
      fwtcm =1.0
      do j=1,nfit
      mnow=mnow+1
      do k=1,nfit
        if (j.ne.k) then
          arspfit(necein+j,k)=0.0
        else
          arspfit(necein+j,k)=fwtcm/fwtnow
        endif
      enddo
      brspfit(necein+j)=0.0
      enddo
      endif
!
      nnn1=1
      iieerr=0
      call sdecm(arspfit,nnecein,mnow,nfit,brspfit,nnecein &
           ,nnn1,s,wk,iieerr)
!EALW      write(*,*) iieerr
!EALW      write(*,*) 's'
!EALW      write(*,*) s
      toler=1.0e-06*s(1)
      DO 2010 I = 1,nfit
            T = 0.0
            IF (S(I).gt.toler) T = Brspfit(I)/S(I)
            Brspfit(I) = T
2010    CONTINUE
2015    DO 2025 I = 1, Nfit
            X(I) = 0.0
            DO 2020 J = 1,nfit
2020                X(I) = X(I) + Arspfit(I,J)*Brspfit(J)
2025    CONTINUE
      do k=1,nfit
      xfit(k)=x(k)
      enddo
!EALW      write(*,*)'x'
!EALW      write(*,*)x
      chisqfit=0.0
      do 2400 k=1,necein
      tte(k)=0.
      do 2500 nk=1,nfit
      tte(k)=tte(k)+x(nk)*an(k)**(nk-1)
2500   continue
      chisqfit=chisqfit+(tte(k)-teecein0(k))**2/tebit(k)
2400   continue
!EALW      write(*,*) 'chisqfit='
!EALW      write(*,*) chisqfit
!EALW      write(*,*)'tte'
!EALW      write(*,*)tte
!--------------------------------------------------------------------
!--  get Teeceb(bbf) in ECE data region, (Te(B)), bbf-B            --
!--------------------------------------------------------------------
      dbbf=(becein(necein)-becein(1))/float(nnnte-1)
      do 3016 i=1,nnnte
          bbf(i)=becein(1)+dbbf*float(i-1)
          bbx=(bbf(i)-b00)/baa
          teeceb(i)=0.
      do 3012 nk=1,nfit
         teeceb(i)=teeceb(i)+x(nk)*bbx**(nk-1)
3012  continue
3016  continue
!---------------------------------------------------------------------
!--   find  beceo which is the B value of Te peak point             --
!---------------------------------------------------------------------
      if (kfixro.eq.1) go to 3019
        teeceo=teeceb(1)
        iio=1
        do 3018 i=2,nnnte
          if (teeceb(i).gt.teeceo) then
             iio=i
             teeceo=teeceb(i)
          endif
3018    continue
         beceo=bbf(iio)
!EALW           write(*,*) 'find beceo, iio,bbf(iio),teeceo'
!EALW           write(*,*) iio,bbf(iio),teeceo
!EALW        write(*,*)'beceo'
!EALW        write(*,*)beceo
!--------------------------------------------------------------------
!--    find becein(idesto), it close to beceo                      --
!--       dTe on beceo from tebit(idesto)                          --
!--------------------------------------------------------------------
        desto=abs(beceo-becein(1))
        idesto=1
        do i=2,necein
        if (abs(beceo-becein(i)).lt.desto) then
        desto=abs(beceo-becein(i))
        idesto=i
        endif
        enddo
!EALW        write(*,*)'idesto'
!EALW        write(*,*)idesto
!--------------------------------------------------------------------
!--    get bobit=dB=sqrt(dTe/Te'')                                 --
!--    Te''-- (d2Te/dB2) at beceo--ppteppbo, dTe=tebit(idesto)     --
!--------------------------------------------------------------------
         bbx1=(bbf(iio+1)-b00)/baa
         bbx2=(bbf(iio-1)-b00)/baa
         ptpr1=0.
         ptpr2=0.
        do nk=2,nfit
           ptpr1=ptpr1+x(nk)*bbx1**(nk-2)
           ptpr2=ptpr2+x(nk)*bbx2**(nk-2)
        enddo
          ptpr1=ptpr1/baa
          ptpr2=ptpr2/baa
        ppteppbo=abs(0.5*(ptpr1-ptpr2)/dbbf)
        dtero=abs(tebit(idesto))
        bobit=sqrt(dtero/ppteppbo)
!EALW        write(*,*)'bobit'
!EALW        write(*,*)bobit
!---------------------------------------------------------------------
!-- take B+ (becep) from becein>beceo and get B- (becem)            --
!-- find nece (the number of B+)                                    --
!--   B+, B- from centre to edge                                    --
!---------------------------------------------------------------------
3019    continue
      if ((kfitece.eq.1).or.(kfixrece.eq.1)) go to 3069
        ii=0
        do k=1,necein
          if ((beceo-becein(k)).lt.0.) then
              ii=ii+1
              becep(ii)=becein(k)
          endif
        enddo
        nece=ii
        do k=1,nece
           bbx=(becep(k)-b00)/baa
           teece(k)=0.
          do 3020 nk=1,nfit
             teece(k)=teece(k)+x(nk)*bbx**(nk-1)
3020      continue
        enddo
!
        ii=0
        do 3025 i=1,nnnte
         if (bbf(i).lt.beceo) then
          ii=ii+1
          blowf(ii)=bbf(i)
          telowf(ii)=teeceb(i)
         endif
3025    continue
!
        nlowf=ii
        call zpline(nlowf,telowf,blowf,bb,cc,dd)
        do 3028 k=1,nece
          becem(k)=seval(nlowf,teece(k),telowf,blowf,bb,cc,dd)
         if ((becem(k).ge.becein(1)).and.(becem(k).lt.beceo)) &
           go to 3028
         fwtece0(k)=0.0
         becem(k)=1.E-6
3028    continue
!--------------------------------------------------------------------
!--   idestm(nece)- the point becein(idestm) close to B-(nece)
!---------------------------------------------------------------------
!
       do k=1,nece
        dest=abs(becem(k)-becein(1))
          idestm(k)=1
          do i=2,necein
           if (abs(becem(k)-becein(i)).lt.dest) then
            dest=abs(becem(k)-becein(i))
            idestm(k)=i
           endif
          enddo
       enddo
!---------------------------------------------------------------------
!--  iteration from 539, before 539 only once calculated            --
!--  Calculation of |B| on rgrid (z=zeceo)   bfield(nw)             --
!---------------------------------------------------------------------
 539  continue
      if (icurrt.ne.1) go to 540
      ffprim(1)=cratio*srma*2.*salpha/darea*twopi*tmu
      ffprim(nw)=ffprim(1)
  540 continue
      if (icurrt.ne.2.and.icurrt.ne.5) go to 550
      ffprim(nw)=fpcurr(x111,kffcur)/darea*twopi*tmu
      ffprim(1)=fpcurr(x000,kffcur)/darea*twopi*tmu
  550 continue
      if (icurrt.ne.4) go to 600
      call currnt(n222,jtime,n222,n222,kerror)
! MPI >>>
#if defined(USEMPI)
      ! NOTE : Serial code does NOT have this error check and does NOT return KERROR
      !if (kerror /= 0) then
      !  ! NOTE : Do NOT need to set KERROR return value because will be set to 1 by CURRNT if error occurred
      !  return
      !endif
#endif
! MPI <<<
      ffprim(1)=rbetap*cratio*rzero*twopi*tmu/darea
      ffprim(nw)=ffprim(1)*gammaf
  600 continue
      do 699 i=2,nw-1
      siii=sigrid(i)
        if (icurrt.ne.2.and.icurrt.ne.5) go to 692
        ffprim(i)=fpcurr(siii,kffcur)/darea*twopi*tmu
  692   continue
        if (icurrt.ne.4) go to 694
        ffprim(i)=ffprim(1)*(1.-siii**enp)**emp*(1.-gammap)+gammap
  694   continue
        if (icurrt.ne.1) go to 696
        ffprim(i)=ffprim(1)
  696   continue
  699 continue
      fpol(nw)=fbrdy*tmu
!EALW      write(*,*)'fpol(nw)'
!EALW      write(*,*)fpol(nw)
      sumf=fpol(nw)**2/2.
!EALW      write(*,*)'psibry'
!EALW      write(*,*)psibry
!EALW      write(*,*)'simag'
!EALW      write(*,*)simag
      delsi=-(psibry-simag)/float(nw-1)
      do 700 i=1,nw-1
        sumf=sumf+0.5*delsi*(ffprim(nw-i+1)+ffprim(nw-i))
        if(sumf .ge. 0.0) then
          fpol(nw-i)=sqrt(2.*sumf)*fpol(nw)/abs(fpol(nw))
        else
          fpol(nw-i)=fpol(nw)
        endif
 700  continue
!EALW      write(*,*)'fpol'
!EALW      write(*,*)fpol
!EALW      write(*,*)'sigrid'
!EALW      write(*,*)sigrid
      call zpline(nw,sigrid,fpol,bbb,ccc,ddd)
      do 702 iw=1,nw
            kk=(iw-1)*nh+jo
            if (xpsi(kk).gt.1.0.or.ivacum.gt.0) then
            fnow=fbrdy*tmu
!EALW            write(*,*)'iw, fnow'
!EALW            write(*,*)iw,fnow
            else
            fnow=seval(nw,xpsi(kk),sigrid,fpol,bbb,ccc,ddd)
!EALW            write(*,*)'iw, xpsi(kk),fnow'
!EALW            write(*,*)iw, xpsi(kk),fnow
            endif
            btttt(iw)=fnow/rgrid(iw)
702   continue
           ier = 0
      call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
         do 708 iw=1,nw
            rw=rgrid(iw)
            rh=zgrid(jo)
            kk=(iw-1)*nh+jo
            call seva2d(bkx,lkx,bky,lky,c,rw,rh,pds,ier,n555)
            if (ier.eq.0) go to 706
            write (nttyo,9910) ier,rw,rh
!EALW            write(*,*) 'ier,rw,rh'
!EALW            write(*,*) ier,rw,rh
            return
 706      continue
            bfield(iw)=sqrt(pds(2)**2+pds(3)**2)/rgrid(iw)
            dsidr(iw)=pds(2)
            ddsiddr(iw)=pds(5)
            bfield(iw)=sqrt(bfield(iw)**2+btttt(iw)**2)
 708   continue
!EALW       write(*,*)'bfield'
!EALW       write(*,*)bfield
 9910 format('  error in spline =',i4,' (r,z)= ( ',f5.2,',',f5.2,')')
!--------------------------------------------------------------------
!--     get d|B|/dr on rgrid
!--------------------------------------------------------------------
       do i=2,nw-1
          dbdr(i)=(bfield(i+1)-bfield(i-1))/(rgrid(i+1)-rgrid(i-1))
       enddo
       dbdr(1)=(bfield(2)-bfield(1))/(rgrid(2)-rgrid(1))
       dbdr(nw)=(bfield(nw)-bfield(nw-1))/(rgrid(nw)-rgrid(nw-1))
!--------------------------------------------------------------------
!--   find the number of bmink and bmaxk (local min and max)
!--------------------------------------------------------------------
      kmin=0
      kmax=0
      do k=2,nw-1
         if((bfield(k).lt.bfield(k+1)).and.(bfield(k).lt.bfield &
                (k-1))) then
            kmin=kmin+1
            bmink(kmin)=bfield(k)
         endif
         if((bfield(k).gt.bfield(k+1)).and.(bfield(k).gt.bfield &
                 (k-1))) then
            kmax=kmax+1
            bmaxk(kmax)=bfield(k)
         endif
      enddo
      if (kmin.eq.(kmax+1)) then
           kmax=kmax+1
           bmaxk(kmax)=bfield(nw)
      endif
      if (kmin.eq.(kmax-1)) then
           kmin=kmin+1
           bmink(1)=bfield(1)
           do k=2,kmin
             bbk(k)=bmink(k-1)
           enddo
           do k=2,kmin
             bmink(k)=bbk(k)
           enddo
      endif
!
      if (kmax.eq.0) go to 1012
!
!--------------------------------------------------------------------
!--   get babsk array (kmax+1) is strictly increasing order
!--------------------------------------------------------------------
      do k=1,kmax
         if(bmaxk(k).lt.bmink(k)) then
!EALW          write(*,*) 'stop by bmaxk1'
!EALW          write(*,*)k,bmaxk(k),bmink(k)
! MPI >>>
#if defined(USEMPI)
          call mpi_stop
#else
          stop
#endif
! MPI <<<
         endif
      enddo
!
      iout1=0
      do i=1,nw
        if (bfield(i).gt.bmaxk(1)) then
          iout1=iout1+1
          bout(1,iout1)=bfield(i)
          rrout(1,iout1)=rgrid(i)
        endif
      enddo
      nnout(1)=iout1
      do k=2,kmax
          ioutk=0
        do i=1,nw
         if((bfield(i).gt.bmaxk(k)).and.(bfield(i).lt.bmink &
             (k-1))) then
          ioutk=ioutk+1
          bout(k,ioutk)=bfield(i)
          rrout(k,ioutk)=rgrid(i)
         endif
        enddo
        nnout(k)=ioutk
      enddo
      ioutk1=0
      do i=1,nw
        if (bfield(i).lt.bmink(kmax)) then
          ioutk1=ioutk1+1
          bout(kmax+1,ioutk1)=bfield(i)
          rrout(kmax+1,ioutk1)=rgrid(i)
         endif    
       enddo
        nnout(kmax+1)=ioutk1
!
      do k=1,kmax+1
        if (nnout(k).gt.3) then
         n=nnout(k)
         do i=1,n
           babs(k,i)=bout(k,n-i+1)
           rrgrid(k,i)=rrout(k,n-i+1)
         enddo
        endif
      enddo
!EALW       write(*,*)'kmax,kmin'
!EALW       write(*,*)kmax,kmin
!EALW       write(*,*)'bmaxk,bmink'
!EALW       write(*,*)bmaxk,bmink
!EALW       write(*,*)'nnout'
!EALW       write(*,*)nnout
!EALW       write(*,*) 'babs, rrgrid'
      do k=1,kmax+1
       n=nnout(k)
       do i=1,n
!EALW         write(*,*) babs(k,i)
!EALW         write(*,*)rrgrid(k,i)
       enddo
      enddo
!-------------------------------------------------------------------
!--   get R-,R+,Ro  where |B| = B+,B-,Bo                          --
!-------------------------------------------------------------------
      do m=1,nece
          recem(m)=1.E-6
          recep(m)=1.E-6
      enddo
      receo=1.e-6
!
      if (nnout(1).gt.3) then
        n=nnout(1)
        do i=1,n
             bx(i)=babs(1,i)
             ry(i)=rrgrid(1,i)
        enddo
        call zpline(n,bx,ry,bbb,ccc,ddd)
!
        if (beceo.ge.bmaxk(1)) then
           receo=seval(n,beceo,bx,ry,bbb,ccc,ddd)
        endif
        do m=1,nece
         if (becep(m).ge.bmaxk(1)) then
          recem(m)=seval(n,becep(m),bx,ry,bbb,ccc,ddd)
         endif
         if (becem(m).ge.bmaxk(1)) then
          recep(m)=seval(n,becem(m),bx,ry,bbb,ccc,ddd)
         endif
        enddo
      endif
!
      do k=2,kmax
       if (nnout(k).gt.3) then
         n=nnout(k)
         do i=1,n
             bx(i)=babs(k,i)
             ry(i)=rrgrid(k,i)
         enddo
         call zpline(n,bx,ry,bbb,ccc,ddd)
        if ((beceo.ge.bmaxk(k)).and.(beceo.le.bmink(k-1))) then
           receo=seval(n,beceo,bx,ry,bbb,ccc,ddd)
        endif
        do m=1,nece
         if((becep(m).ge.bmaxk(k)).and.(becep(m).le.bmink(k-1)))then
          recem(m)=seval(n,becep(m),bx,ry,bbb,ccc,ddd)
         endif
         if((becem(m).ge.bmaxk(k)).and.(becem(m).le.bmink(k-1)))then
          recep(m)=seval(n,becem(m),bx,ry,bbb,ccc,ddd)
         endif
        enddo
       endif
      enddo
!
        if (nnout(kmax+1).gt.3) then
         n=nnout(kmax+1)
         do i=1,n
             bx(i)=babs(kmax+1,i)
             ry(i)=rrgrid(kmax+1,i)
         enddo
         call zpline(n,bx,ry,bbb,ccc,ddd)
        if (beceo.le.bmink(kmax)) then
           receo=seval(n,beceo,bx,ry,bbb,ccc,ddd)
        endif
         do m=1,nece
           if(becep(m).le.bmink(kmax))then
             recem(m)=seval(n,becep(m),bx,ry,bbb,ccc,ddd)
           endif
           if(becem(m).le.bmink(kmax))then
             recep(m)=seval(n,becem(m),bx,ry,bbb,ccc,ddd)
           endif
         enddo
        endif
!
      do m=1,nece
          if((recep(m).lt.1.E-5).or.(recem(m).lt.1.E-5)) then
            fwtece0(m)=0.0
            recep(m)=1.0E-6
            recem(m)=1.0E-6
          endif
      enddo
      if (receo.lt.1.E-5) fwtecebz0=0.0
      go to 1015
!
1012  do i=1,nw
          bx(i)=bfield(nw-i+1)
          ry(i)=rgrid(nw-i+1)
      enddo
      call zpline(nw,bx,ry,bbb,ccc,ddd)
       do m=1,nece
         recep(m)=seval(nw,becem(m),bx,ry,bbb,ccc,ddd)
         recem(m)=seval(nw,becep(m),bx,ry,bbb,ccc,ddd)
      enddo
      receo=seval(nw,beceo,bx,ry,bbb,ccc,ddd)
1015    continue
!EALW      write(*,*)'recem'
!EALW      write(*,*)recem
!EALW      write(*,*)'recep'
!EALW      write(*,*)recep
!EALW      write(*,*)'receo'
!EALW      write(*,*)receo
!EALW      write(*,*)'nece'
!EALW      write(*,*)nece
!------------------------------------------------------------------
!--   get dB/dr at receo (dbdro) and recep,recem (dbdrp,dbdrm)
!------------------------------------------------------------------
      call zpline(nw,rgrid,dbdr,bbb,ccc,ddd)
         if (fwtecebz0.gt.1.e-6) then
                 dbdro=seval(nw,receo,rgrid,dbdr,bbb,ccc,ddd)
         endif
         do k=1,nece
           if (fwtece0(k).gt.1.e-6) then
                 dbdrp(k)=seval(nw,recep(k),rgrid,dbdr,bbb,ccc,ddd)
                 dbdrm(k)=seval(nw,recem(k),rgrid,dbdr,bbb,ccc,ddd)
           endif
         enddo  
!--------------------------------------------------------------------
!--   get robit-- from bobit and dB/dr, (robit=dB/(dB/dr),bobit--dB)
!--------------------------------------------------------------------
         robit=bobit/dbdro
!--------------------------------------------------------------------
!--       get ecebit from sqrt(bitm**2+bitp**2) ,
!--         bit(m,p)=dTe *(dpsi/dR)/(dTe/dB)/(dB/dR)
!--         dTe(m)=tebit(idestm),  dTe(p)=tebit(idestp)
!--         (dTe/dB)*(dB/dR)--pteprp,pteprm
!--         (dpsi/dR)--dsidrp,dsidrm
!---------------------------------------------------------------------
        do k=1,nece
         if (fwtece0(k).gt.1.e-6) then 
           bbxp=(becep(k)-b00)/baa
           bbxm=(becem(k)-b00)/baa
           pteprp(k)=0.
           pteprm(k)=0.
          do 3030 nk=2,nfit
             pteprp(k)=pteprp(k)+x(nk)*bbx**(nk-2)
             pteprm(k)=pteprm(k)+x(nk)*bbx**(nk-2)
3030      continue
          pteprp(k)=pteprp(k)/baa
          pteprm(k)=pteprm(k)/baa
          pteprp(k)=pteprp(k)*dbdrp(k)
          pteprm(k)=pteprm(k)*dbdrm(k)
         endif
        enddo
!
        call zpline(nw,rgrid,dsidr,bbb,ccc,ddd)
      do 3036 k=1,nece
       if (fwtece0(k).gt.1.e-6) then
         dsidrm=seval(nw,recem(k),rgrid,dsidr,bbb,ccc,ddd)
         dsidrp=seval(nw,recep(k),rgrid,dsidr,bbb,ccc,ddd)
        if((abs(pteprm(k)).gt.1.E-10).and.(abs(pteprp(k)).gt.1.E-10)) &
           then
         imk=idestm(k)
          rmbit(k)=tebit(imk)/pteprm(k)
          bitm=rmbit(k)*dsidrm
         ipk=idestp(k)
          rpbit(k)=tebit(ipk)/pteprp(k)
          bitp=rpbit(k)*dsidrp
          ecebit(k)=sqrt(bitm**2+bitp**2)
         else
          fwtece0(k)=0.
         endif
        endif
3036  continue
!EALW       write(*,*)'rmbit'
!EALW       write(*,*)rmbit
!EALW       write(*,*)'rpbit'
!EALW       write(*,*)rpbit
3069  continue
!
      DEALLOCATE(rrgrid,bfield,rrout,bout,babs,bbb,ccc,ddd,btttt, &
          dsidr,ddsiddr,bx,ry,bbk,dbdr)
!
      return
        end
!
       subroutine getecer(jtime,kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          getecer obtains the receo, R+ R-                        **
!**          from ECE measurement data                               **
!**          if kfixro  kfixrece = 0 called in setece                **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          11/98..........first created, Cheng Zhang               **
!**     2013/08/07..........Update for real-space Ti                 **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!
      parameter (nn=30)
      parameter (kbre=5)
      common/cwork3/lkx,lky
      dimension pds(6),nnout(kbre),bmink(kbre),bmaxk(kbre)
      real*8,allocatable :: rrgrid(:,:),bfield(:),rrout(:,:), &
          bout(:,:),babs(:,:),bbb(:),ccc(:),ddd(:),btttt(:), &
          dsidr(:),ddsiddr(:),bx(:),ry(:),bbk(:)
      dimension arspfit(nnecein,nn),brspfit(nnecein) &
           ,s(nn),tte(nnecein),x(nn),an(nnecein) &
           ,tebit(nnecein)
      dimension telowf(nnnte) &
         ,rlowf(nnnte),bb(nnnte),cc(nnnte),dd(nnnte) &
         ,teece(nnece),pteprm(nnece),pteprp(nnece) &
         ,idestp(nnece),idestm(nnece)
      integer, intent(inout) :: kerror
!
      kerror = 0
      ALLOCATE(rrgrid(kbre,nw),bfield(nw),rrout(kbre,nw), &
         bout(kbre,nw),babs(kbre,nw),bbb(nw),ccc(nw), &
         ddd(nw),btttt(nw),dsidr(nw),ddsiddr(nw),bx(nw), &
         ry(nw),bbk(nw))
!
      do k=1,nnece
         fwtece0(k)=swtece(k) 
      enddo
      fwtecebz0=swtecebz
      do k=1,kbre
       do i=1,nw
         babs(k,i)=0.0
         bout(k,i)=0.0
         rrout(k,i)=0.0
         rrgrid(k,i)=0.0
       enddo
      enddo
!---------------------------------------------------------------------
!--  Calculation of |B| array from fe array ( harmonic nharm)       --
!--     becein(necein),   fe(GHz),|B|(T) becein form H.f to L.f     --
!---------------------------------------------------------------------
      do 500 k=1,necein
        becein(k)=0.001*6.0*9.1095*3.14159/4.8032*feece(k)/float(nharm)
500   continue
!---------------------------------------------------------------------
!--  Calculation of |B| on rgrid (z=zeceo)   bfield(nw)             --
!---------------------------------------------------------------------
      if (icurrt.ne.1) go to 540
      ffprim(1)=cratio*srma*2.*salpha/darea*twopi*tmu
      ffprim(nw)=ffprim(1)
  540 continue
      if (icurrt.ne.2.and.icurrt.ne.5) go to 550
      ffprim(nw)=fpcurr(x111,kffcur)/darea*twopi*tmu
      ffprim(1)=fpcurr(x000,kffcur)/darea*twopi*tmu
  550 continue
      if (icurrt.ne.4) go to 600
      call currnt(n222,jtime,n222,n222,kerror)
      ffprim(1)=rbetap*cratio*rzero*twopi*tmu/darea
      ffprim(nw)=ffprim(1)*gammaf
  600 continue
      do 699 i=2,nw-1
        siii=sigrid(i)
        if (icurrt.ne.2.and.icurrt.ne.5) go to 692
        ffprim(i)=fpcurr(siii,kffcur)/darea*twopi*tmu
  692   continue
        if (icurrt.ne.4) go to 694
        ffprim(i)=ffprim(1)*(1.-siii**enp)**emp*(1.-gammap)+gammap
  694   continue
        if (icurrt.ne.1) go to 696
        ffprim(i)=ffprim(1)
  696   continue
  699 continue
      fpol(nw)=fbrdy*tmu
      sumf=fpol(nw)**2/2.
      delsi=-(psibry-simag)/float(nw-1)
      do 700 i=1,nw-1
        sumf=sumf+0.5*delsi*(ffprim(nw-i+1)+ffprim(nw-i))
        if(sumf .ge. 0.0) then
          fpol(nw-i)=sqrt(2.*sumf)*fpol(nw)/abs(fpol(nw))
        else
          fpol(nw-i)=fpol(nw)
        endif
 700  continue
      call zpline(nw,sigrid,fpol,bbb,ccc,ddd)
      do 702 iw=1,nw
            kk=(iw-1)*nh+jo
            if (xpsi(kk).gt.1.0.or.ivacum.gt.0) then
            fnow=fbrdy*tmu
            else
            fnow=seval(nw,xpsi(kk),sigrid,fpol,bbb,ccc,ddd)
            endif
            btttt(iw)=fnow/rgrid(iw)
702   continue
!
      ier = 0
      call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
      do 708 iw=1,nw
            rw=rgrid(iw)
            rh=zgrid(jo)
            kk=(iw-1)*nh+jo
            call seva2d(bkx,lkx,bky,lky,c,rw,rh,pds,ier,n555)
            if (ier.eq.0) go to 706
            write (nttyo,9910) ier,rw,rh
            return
 706        continue
            bfield(iw)=sqrt(pds(2)**2+pds(3)**2)/rgrid(iw)
            dsidr(iw)=pds(2)
            ddsiddr(iw)=pds(5)            
            bfield(iw)=sqrt(bfield(iw)**2+btttt(iw)**2)
 708  continue
 9910 format('  error in spline =',i4,' (r,z)= ( ',f5.2,',',f5.2,')')
!--------------------------------------------------------------------
!--   find the number of bmink and bmaxk (local min and max)
!--------------------------------------------------------------------
      kmin=0
      kmax=0
      do k=2,nw-1
         if((bfield(k).lt.bfield(k+1)).and.(bfield(k).lt.bfield &
                (k-1))) then
            kmin=kmin+1
            bmink(kmin)=bfield(k)
         endif
         if((bfield(k).gt.bfield(k+1)).and.(bfield(k).gt.bfield &
                 (k-1))) then
            kmax=kmax+1
            bmaxk(kmax)=bfield(k)
         endif
      enddo
      if (kmin.eq.(kmax+1)) then
           kmax=kmax+1  
           bmaxk(kmax)=bfield(nw)
      endif
      if (kmin.eq.(kmax-1)) then
           kmin=kmin+1
           bmink(1)=bfield(1)
           do k=2,kmin
             bbk(k)=bmink(k-1)
           enddo
           do k=2,kmin
             bmink(k)=bbk(k)
           enddo
      endif 
!
      if (idebug.ge.3) write (6,*) 'GETECER, kmin/kmax = ', kmin, kmax
      if (kmax.eq.0) go to 1012
!--------------------------------------------------------------------
!--   get babsk array (kmax+1) is strictly increasing order
!-------------------------------------------------------------------- 
      do k=1,kmax
         if(bmaxk(k).lt.bmink(k)) then
           stop
         endif
      enddo
!
      iout1=0
      do i=1,nw
        if (bfield(i).gt.bmaxk(1)) then
          iout1=iout1+1
          bout(1,iout1)=bfield(i)
          rrout(1,iout1)=rgrid(i)
        endif
      enddo
      nnout(1)=iout1
      do k=2,kmax
        ioutk=0
        do i=1,nw
         if((bfield(i).gt.bmaxk(k)).and.(bfield(i).lt.bmink &
             (k-1))) then
          ioutk=ioutk+1
          bout(k,ioutk)=bfield(i)
          rrout(k,ioutk)=rgrid(i)
         endif
        enddo
        nnout(k)=ioutk
      enddo
      ioutk1=0
      do i=1,nw
        if (bfield(i).lt.bmink(kmax)) then
          ioutk1=ioutk1+1
          bout(kmax+1,ioutk1)=bfield(i)
          rrout(kmax+1,ioutk1)=rgrid(i)
         endif       
       enddo
        nnout(kmax+1)=ioutk1
!
      do k=1,kmax+1
        if (nnout(k).gt.3) then
         n=nnout(k)
         do i=1,n
           babs(k,i)=bout(k,n-i+1)
           rrgrid(k,i)=rrout(k,n-i+1)
         enddo
        endif
      enddo 
!-------------------------------------------------------------------
!--   get recein, at which |B| = becein                           --
!--     recein(mecein)                                            --
!-------------------------------------------------------------------
      kout=0
!
      if (nnout(1).gt.3) then
      n=nnout(1)
      do i=1,n
        bx(i)=babs(1,i)
        ry(i)=rrgrid(1,i)
      enddo
      call zpline(n,bx,ry,bbb,ccc,ddd)
      do m=1,necein
        if (becein(m).ge.bmaxk(1)) then
         kout=kout+1
         recein(kout)=seval(n,becein(m),bx,ry,bbb,ccc,ddd) 
         teeceinr(kout)=teecein(m)
          tebit(kout)=errorece(m)         
        endif
      enddo
      endif
!
      do k=2,kmax
       if (nnout(k).gt.3) then
         n=nnout(k)
         do i=1,n
             bx(i)=babs(k,i)
             ry(i)=rrgrid(k,i)
         enddo
         call zpline(n,bx,ry,bbb,ccc,ddd)
         do m=1,necein
          if((becein(m).ge.bmaxk(k)).and.(becein(m).le.bmink(k-1)))then
            kout=kout+1
            recein(kout)=seval(n,becein(m),bx,ry,bbb,ccc,ddd)
            teeceinr(kout)=teecein(m)
            tebit(kout)=errorece(m)
          endif
         enddo
       endif
      enddo
!
        if (nnout(kmax+1).gt.3) then
         n=nnout(kmax+1)
         do i=1,n
             bx(i)=babs(kmax+1,i)
             ry(i)=rrgrid(kmax+1,i)
         enddo
         call zpline(n,bx,ry,bbb,ccc,ddd)
         do m=1,necein
           if(becein(m).le.bmink(kmax))then
             kout=kout+1
             recein(kout)=seval(n,becein(m),bx,ry,bbb,ccc,ddd)
             teeceinr(kout)=teecein(m)
             tebit(kout)=errorece(m)
           endif
         enddo
        endif
! 
      mecein=kout
!
      go to 1015
!
1012    do  i=1,nw
          bx(i)=bfield(nw-i+1)
          ry(i)=rgrid(nw-i+1)
        enddo
       call zpline(nw,bx,ry,bbb,ccc,ddd)
       do m=1,necein
         recein(m)=seval(nw,becein(m),bx,ry,bbb,ccc,ddd)
         teeceinr(m)=teecein(m)
         tebit(m)=errorece(m)
      enddo
      mecein=necein
!
1015    continue
!--------------------------------------------------------------------
!--   fitting data from teeceinr,tebit and recein (nnecein)        --
!--     rx=(R-r00)/raa                                             --
!--     Te=x(1)+x(2)*rx+x(3)*rx**2+...+x(nfit)*rx**(nfit-1)        --
!--------------------------------------------------------------------
!heng          mm--nnecein   m---mecein  nn--parameter, n--nfit input
      rmin=recein(1)
      rmax=recein(mecein)
      raa=0.5*(rmax-rmin)
      r00=0.5*(rmax+rmin)
      do i=1,mecein
        an(i)=(recein(i)-r00)/raa
        tebit(i)=max(tebit(i),dble(1.e-4))
      enddo
      do 1110 nj=1,mecein
      do 1100 nk=1,nfit
      if (nk.eq.1) then
        arspfit(nj,nk)=1./tebit(nj)
      else
        arspfit(nj,nk)=an(nj)**(nk-1)/tebit(nj)
      endif
1100   continue
1110   continue
      do 1200 nj=1,mecein
      brspfit(nj)=teeceinr(nj)/tebit(nj)
1200   continue
!
      mnow=mecein
      if (kcmin.gt.0) then
      fwtnow=0.001
      fwtcm =1.0
      do j=1,nfit
      mnow=mnow+1
      do k=1,nfit
        if (j.ne.k) then
          arspfit(mecein+j,k)=0.0
        else
          arspfit(mecein+j,k)=fwtcm/fwtnow
        endif
      enddo
      brspfit(mecein+j)=0.0
      enddo
      endif
!
      nnn1=1
      iieerr=0
      call sdecm(arspfit,nnecein,mnow,nfit,brspfit,nnecein &
           ,nnn1,s,wk,iieerr)
      toler=1.0e-06*s(1)
      DO 2010 I = 1,nfit
            T = 0.0
            IF (S(I).gt.toler) T = Brspfit(I)/S(I)
            Brspfit(I) = T
2010  CONTINUE
2015  DO 2025 I = 1, Nfit
            X(I) = 0.0
            DO 2020 J = 1,nfit
2020                X(I) = X(I) + Arspfit(I,J)*Brspfit(J)
2025  CONTINUE
      do k=1,nfit
      xfit(k)=x(k)
      enddo
      chisqfit=0.0
      do 2400 k=1,mecein
      tte(k)=0.
      do 2500 nk=1,nfit
      tte(k)=tte(k)+x(nk)*an(k)**(nk-1)
2500  continue
      chisqfit=chisqfit+(tte(k)-teeceinr(k))**2/tebit(k)
2400  continue
!--------------------------------------------------------------------
!--  get Teecer(rrr) in ECE data region                            --
!--------------------------------------------------------------------
      drrr=(recein(mecein)-recein(1))/float(nnnte-1)
      do 3016 i=1,nnnte
          rrr(i)=recein(1)+drrr*float(i-1)
          rx=(rrr(i)-r00)/raa
          teecer(i)=0.
      do 3012 nk=1,nfit
         teecer(i)=teecer(i)+x(nk)*rx**(nk-1)
3012  continue
3016  continue  
!---------------------------------------------------------------------
!--   find receo  which is Te peak point                            --
!---------------------------------------------------------------------
      if (kfixro.eq.1) go to 3019      
        teeceo=teecer(1)
        iio=1
        do 3018 i=2,nnnte
          if (teecer(i).gt.teeceo) then
             iio=i
             teeceo=teecer(i)
          endif
3018    continue
        receo=rrr(iio)
!--------------------------------------------------------------------
!--    find recein(idesto), it close to receo                      --
!--       dTe on receo from tebit(idesto)                          --
!--------------------------------------------------------------------
        desto=abs(receo-recein(1))
        idesto=1
        do i=2,mecein
        if (abs(receo-recein(i)).lt.desto) then
        desto=abs(receo-recein(i))
        idesto=i
        endif
        enddo
!--------------------------------------------------------------------
!--    get robit,  robit=sqrt(dTe/Te'')                            --
!--    Te''-- (d2Te/dR2) at receo--ppteppro, dTe=tebit(idesto)     --
!--------------------------------------------------------------------
        rx1=(rrr(iio+1)-r00)/raa
        rx2=(rrr(iio-1)-r00)/raa
        ptpr1=0.
        ptpr2=0.
        do nk=2,nfit
           ptpr1=ptpr1+x(nk)*rx1**(nk-2)
           ptpr2=ptpr2+x(nk)*rx2**(nk-2)
        enddo
        ptpr1=ptpr1/raa
        ptpr2=ptpr2/raa
        ppteppro=abs(0.5*(ptpr1-ptpr2)/drrr)
        dtero=abs(tebit(idesto))
        robit=sqrt(dtero/ppteppro)
!---------------------------------------------------------------------
!-- take R- and get R+                                              --
!--         nece=the number of R-,  recem(nece), recep(nece)        --
!---------------------------------------------------------------------
3019  continue
      if ((kfitece.eq.1).or.(kfixrece.eq.1)) go to 3069
        ii=0
        do k=mecein,1,-1
          if ((receo-recein(k)).gt.0.) then
              ii=ii+1
              recem(ii)=recein(k)
          endif
        enddo
        nece=ii
        do k=1,nece
           rx=(recem(k)-r00)/raa
           teece(k)=0.
           pteprm(k)=0.
          do 3020 nk=1,nfit
             teece(k)=teece(k)+x(nk)*rx**(nk-1)
3020      continue
          do 3021 nk=2,nfit
             pteprm(k)=pteprm(k)+x(nk)*rx**(nk-2)
3021      continue          
          pteprm(k)=pteprm(k)/raa
        enddo
!
        ii=0
        do 3025 i=nnnte,1,-1
         if (rrr(i).gt.receo) then
          ii=ii+1
          rlowf(ii)=rrr(i)
          telowf(ii)=teecer(i)
         endif
3025    continue
        nlowf=ii
        call zpline(nlowf,telowf,rlowf,bb,cc,dd)
        do 3028 k=1,nece
          recep(k)=seval(nlowf,teece(k),telowf,rlowf,bb,cc,dd)
         if ((recep(k).gt.receo).and.(recep(k).lt.recein(mecein))) &
           go to 3028
         fwtece0(k)=0.0
3028    continue
!--------------------------------------------------------------------
!--   idestp(nece)- the point recein(idestp) close to R+(nece)
!--   idestm(nece)- the point recein(idestm) close to R-(nece)
!---------------------------------------------------------------------
       do k=1,nece
        dest=abs(recep(k)-recein(1))
          idestp(k)=1
          do i=2,mecein
           if (abs(recep(k)-recein(i)).lt.dest) then
            dest=abs(recep(k)-recein(i))
            idestp(k)=i
           endif
          enddo
       enddo 
! 
       do k=1,nece
        dest=abs(recem(k)-recein(1))
          idestm(k)=1
          do i=2,mecein
           if (abs(recem(k)-recein(i)).lt.dest) then
            dest=abs(recem(k)-recein(i))
            idestm(k)=i
           endif
          enddo
       enddo 
!--------------------------------------------------------------------
!--       get ecebit from sqrt(bitm**2+bitp**2) ,
!--         bit(m,p)=dTe *(dpsi/dR)/(dTe/dR)
!--         dTe(m)=tebit(idestm),  dTe(p)=tebit(idestp)
!---------------------------------------------------------------------
        do k=1,nece
           rx=(recep(k)-r00)/raa
           pteprp(k)=0.
          do 3030 nk=2,nfit
             pteprp(k)=pteprp(k)+x(nk)*rx**(nk-2)
3030      continue
          pteprp(k)=pteprp(k)/raa
        enddo
!
        call zpline(nw,rgrid,dsidr,bbb,ccc,ddd)
        do 3036 k=1,nece
         dsidrm=seval(nw,recem(k),rgrid,dsidr,bbb,ccc,ddd)
         dsidrp=seval(nw,recep(k),rgrid,dsidr,bbb,ccc,ddd)
         if((abs(pteprm(k)).gt.1.E-10).and.(abs(pteprp(k)).gt.1.E-10)) &
           then
            imk=idestm(k)
            rmbit(k)=tebit(imk)/pteprm(k)
            bitm=rmbit(k)*dsidrm
            ipk=idestp(k)
            rpbit(k)=tebit(ipk)/pteprp(k)  
            bitp=rpbit(k)*dsidrp
            ecebit(k)=sqrt(bitm**2+bitp**2) 
           else
            fwtece0(k)=0.
         endif
3036    continue
3069  continue
!
      DEALLOCATE(rrgrid,bfield,rrout,bout,babs,bbb,ccc,ddd,btttt, &
          dsidr,ddsiddr,bx,ry,bbk)
!
      return
      end
      subroutine gettir(jtime,kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          gettir obtains the receo, R+ R-                         **
!**          from Ti data                                            **
!**          kfixro = 0, kfixrece = 3, called from setece            **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**     2013/08/07..........Update for real-space Ti based on ECE/Te **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!
      parameter (nn=30)
      parameter (kbre=5)
      common/cwork3/lkx,lky
      dimension pds(6),nnout(kbre),bmink(kbre),bmaxk(kbre)
      real*8,allocatable :: rrgrid(:,:),bfield(:),rrout(:,:), &
          bout(:,:),babs(:,:),bbb(:),ccc(:),ddd(:),btttt(:), &
          dsidr(:),ddsiddr(:),bx(:),ry(:),bbk(:)
      dimension arspfit(nnecein,nn),brspfit(nnecein) &
           ,s(nn),tte(nnecein),x(nn),an(nnecein) &
           ,tebit(nnecein)
      dimension telowf(nnnte) &
         ,rlowf(nnnte),bb(nnnte),cc(nnnte),dd(nnnte) &
         ,teece(nnece),pteprm(nnece),pteprp(nnece) &
         ,idestp(nnece),idestm(nnece)
      integer, intent(inout) :: kerror
!
      if (idebug.ge.3) write (6,*) 'Enter GETTIR, kfitece/kfixrece = ',&
         kfitece, kfixrece
      kerror = 0
      ALLOCATE(rrgrid(kbre,nw),bfield(nw),rrout(kbre,nw), &
         bout(kbre,nw),babs(kbre,nw),bbb(nw),ccc(nw), &
         ddd(nw),btttt(nw),dsidr(nw),ddsiddr(nw),bx(nw), &
         ry(nw),bbk(nw))
!
      do k=1,nnece
         fwtece0(k)=swtece(k) 
      enddo
      fwtecebz0=swtecebz
      do k=1,kbre
       do i=1,nw
         babs(k,i)=0.0
         bout(k,i)=0.0
         rrout(k,i)=0.0
         rrgrid(k,i)=0.0
       enddo
      enddo
!---------------------------------------------------------------------
!--  Copy Ti  array                                                 --
!---------------------------------------------------------------------
      do 500 k=1,necein
        becein(k)=feece(k)
        recein(k)=becein(k)
        teeceinr(k)=teecein(k)
        tebit(k)=errorece(k)
500   continue
!---------------------------------------------------------------------
!--  Calculation of |B| on rgrid (z=zeceo)   bfield(nw)             --
!---------------------------------------------------------------------
      do 708 iw=1,nw
            rw=rgrid(iw)
            rh=zteo
            kk=(iw-1)*nh+jo
            call seva2d(bkx,lkx,bky,lky,c,rw,rh,pds,ier,n555)
            if (ier.eq.0) go to 706
              write (nttyo,9910) ier,rw,rh
            return
 706        continue
            dsidr(iw)=pds(2)
            ddsiddr(iw)=pds(5)            
 708  continue
 9910 format('  error in getecer/spline = ',i4,' (r,z)= ( ', &
                f5.2,',',f5.2,')')
!--------------------------------------------------------------------
!--   fitting data from teeceinr,tebit and recein (nnecein)        --
!--     rx=(R-r00)/raa                                             --
!--     Te=x(1)+x(2)*rx+x(3)*rx**2+...+x(nfit)*rx**(nfit-1)        --
!--     mm--nnecein   m---mecein  nn--parameter, n--nfit input     --
!--------------------------------------------------------------------
      mecein=necein
      rmin=recein(1)
      rmax=recein(mecein)
      raa=0.5*(rmax-rmin)
      r00=0.5*(rmax+rmin)
      do i=1,mecein
        an(i)=(recein(i)-r00)/raa
        tebit(i)=max(tebit(i),dble(1.e-4))
      enddo
      do 1110 nj=1,mecein
      do 1100 nk=1,nfit
        if (nk.eq.1) then
          arspfit(nj,nk)=1./tebit(nj)
        else
          arspfit(nj,nk)=an(nj)**(nk-1)/tebit(nj)
        endif
1100  continue
1110  continue
      do 1200 nj=1,mecein
      brspfit(nj)=teeceinr(nj)/tebit(nj)
1200   continue
!
      if (idebug.ge.3) write (6,*) 'GETTIR/SDECM, mecein/raa/r00 = ',&
         mecein, raa, r00, nfit
      if (idebug.ge.3) write (6,*) 'GETTIR teeceinr = ', &
           (teeceinr(i),i=1,mecein)
      if (idebug.ge.3) write (6,*) 'GETTIR tebit = ', &
           (tebit(i),i=1,mecein)
      nnn1=1
      iieerr=0
      mnow=mecein
      call sdecm(arspfit,nnecein,mnow,nfit,brspfit,nnecein &
           ,nnn1,s,wk,iieerr)
      toler=1.0e-06*s(1)
      DO 2010 I = 1,nfit
        T = 0.0
        IF (S(I).gt.toler) T = Brspfit(I)/S(I)
        Brspfit(I) = T
2010  CONTINUE
2015  DO 2025 I = 1, Nfit
        X(I) = 0.0
        DO 2020 J = 1,nfit
2020      X(I) = X(I) + Arspfit(I,J)*Brspfit(J)
2025  CONTINUE
!
      do k=1,nfit
        xfit(k)=x(k)
      enddo
      chisqfit=0.0
      do 2400 k=1,mecein
      tte(k)=0.
      do 2500 nk=1,nfit
        tte(k)=tte(k)+x(nk)*an(k)**(nk-1)
2500  continue
      chisqfit=chisqfit+(tte(k)-teeceinr(k))**2/tebit(k)
2400  continue
      mmmte = nnnte
      if (idebug.ge.3) write (6,*) 'GETTIR chisqfit/kfixro/mnow = ', &
         chisqfit, kfixro, mnow
      if (idebug.ge.3) write (6,*) 'GETTIR tte = ', &
           (tte(i),i=1,mecein)
!--------------------------------------------------------------------
!--  get Teecer(rrr) in ECE data region                            --
!--------------------------------------------------------------------
      drrr=(recein(mecein)-recein(1))/float(nnnte-1)
      do 3016 i=1,nnnte
          rrr(i)=recein(1)+drrr*float(i-1)
          rx=(rrr(i)-r00)/raa
          teecer(i)=0.
      do 3012 nk=1,nfit
         teecer(i)=teecer(i)+x(nk)*rx**(nk-1)
3012  continue
3016  continue  
!---------------------------------------------------------------------
!--   find receo  which is Te peak point                            --
!---------------------------------------------------------------------
      if (kfixro.eq.1) go to 3019
        teeceo=teecer(1)
        iio=1
        do 3018 i=2,nnnte
          if (teecer(i).gt.teeceo) then
             iio=i
             teeceo=teecer(i)
          endif
3018    continue
        receo=rrr(iio)
        if (idebug.ge.3) write (6,*) 'GETTIR teece, receo, iio = ', &
           teeceo, receo, iio
!--------------------------------------------------------------------
!--    find recein(idesto), it close to receo                      --
!--       dTe on receo from tebit(idesto)                          --
!--------------------------------------------------------------------
        desto=abs(receo-recein(1))
        idesto=1
        do i=2,mecein
        if (abs(receo-recein(i)).lt.desto) then
        desto=abs(receo-recein(i))
        idesto=i
        endif
        enddo
!--------------------------------------------------------------------
!--    get robit,  robit=sqrt(2dTe/Te'')                            --
!--    Te''-- (d2Te/dR2) at receo--ppteppro, dTe=tebit(idesto)     --
!--------------------------------------------------------------------
        rx1=(rrr(iio+1)-r00)/raa
        rx2=(rrr(iio-1)-r00)/raa
        ptpr1=0.
        ptpr2=0.
        do nk=2,nfit
           ptpr1=ptpr1+x(nk)*rx1**(nk-2)*(nk-1)
           ptpr2=ptpr2+x(nk)*rx2**(nk-2)*(nk-1)
        enddo
        ptpr1=ptpr1/raa
        ptpr2=ptpr2/raa
        ppteppro=abs(0.5*(ptpr1-ptpr2)/drrr)
        dtero=abs(tebit(idesto))
        robit=sqrt(2*dtero/ppteppro)
!---------------------------------------------------------------------
!-- take R- and get R+                                              --
!--         nece=the number of R-,  recem(nece), recep(nece)        --
!---------------------------------------------------------------------
3019  continue
      if (idebug.ge.3) write (6,*) 'GETTIR R-, kfitece/kfixrece  = ', &
        kfitece, kfixrece
      if ((kfitece.eq.1).or.(kfixrece.eq.1)) go to 3069
        ii=0
!       do k=mecein,1,-1
        do k=1,mecein
          if ((receo-recein(k)).gt.0.) then
              ii=ii+1
              recem(ii)=recein(k)
          endif
        enddo
        nece=ii
        do k=1,nece
          rx=(recem(k)-r00)/raa
          teece(k)=0.
          pteprm(k)=0.
          do 3020 nk=1,nfit
            teece(k)=teece(k)+x(nk)*rx**(nk-1)
3020      continue
          do 3021 nk=2,nfit
            pteprm(k)=pteprm(k)+x(nk)*rx**(nk-2)*(nk-1)
3021      continue          
          pteprm(k)=pteprm(k)/raa
        enddo
        if (idebug.ge.3) write (6,*) 'GETTIR R-, nece = ', &
          nece
        if (idebug.ge.3) write (6,*) 'GETTIR R-, recem = ', &
          (recem(i),i=1,nece)
        if (idebug.ge.3) write (6,*) 'GETTIR R-, teece = ', &
          (teece(i),i=1,nece)
        if (idebug.ge.3) write (6,*) 'GETTIR R-, pteprm = ', &
          (pteprm(i),i=1,nece)
!
        ii=0
        do 3025 i=nnnte,1,-1
!       do 3025 i=1,nnnte
         if (rrr(i).gt.receo) then
          ii=ii+1
          rlowf(ii)=rrr(i)
          telowf(ii)=teecer(i)
         endif
3025    continue
        nlowf=ii

!
        call zpline(nlowf,telowf,rlowf,bb,cc,dd)
        do 3028 k=1,nece
         recep(k)=seval(nlowf,teece(k),telowf,rlowf,bb,cc,dd)
         if ((recep(k).gt.receo).and.(recep(k).lt.recein(mecein))) &
           go to 3028
         fwtece0(k)=0.0
3028    continue
        if (idebug.ge.3) write (6,*) 'GETTIR R+, recep = ', &
          (recep(i),i=1,nece)
!--------------------------------------------------------------------
!--   idestp(nece)- the point recein(idestp) close to R+(nece)
!--   idestm(nece)- the point recein(idestm) close to R-(nece)
!---------------------------------------------------------------------
       do k=1,nece
        dest=abs(recep(k)-recein(1))
          idestp(k)=1
          do i=2,mecein
           if (abs(recep(k)-recein(i)).lt.dest) then
            dest=abs(recep(k)-recein(i))
            idestp(k)=i
           endif
          enddo
       enddo 
! 
       do k=1,nece
        dest=abs(recem(k)-recein(1))
          idestm(k)=1
          do i=2,mecein
           if (abs(recem(k)-recein(i)).lt.dest) then
            dest=abs(recem(k)-recein(i))
            idestm(k)=i
           endif
          enddo
       enddo 
!--------------------------------------------------------------------
!--       get ecebit from sqrt(bitm**2+bitp**2) ,
!--         bit(m,p)=dTe *(dpsi/dR)/(dTe/dR)
!--         dTe(m)=tebit(idestm),  dTe(p)=tebit(idestp)
!---------------------------------------------------------------------
        do k=1,nece
           rx=(recep(k)-r00)/raa
           pteprp(k)=0.
          do 3030 nk=2,nfit
             pteprp(k)=pteprp(k)+x(nk)*rx**(nk-2)*(nk-1)
3030      continue
          pteprp(k)=pteprp(k)/raa
        enddo
!
        call zpline(nw,rgrid,dsidr,bbb,ccc,ddd)
        do 3036 k=1,nece
         dsidrm=seval(nw,recem(k),rgrid,dsidr,bbb,ccc,ddd)
         dsidrp=seval(nw,recep(k),rgrid,dsidr,bbb,ccc,ddd)
         if((abs(pteprm(k)).gt.1.E-10).and.(abs(pteprp(k)).gt.1.E-10)) &
           then
            imk=idestm(k)
            rmbit(k)=tebit(imk)/pteprm(k)
            bitm=rmbit(k)*dsidrm
            ipk=idestp(k)
            rpbit(k)=tebit(ipk)/pteprp(k)  
            bitp=rpbit(k)*dsidrp
            ecebit(k)=sqrt(bitm**2+bitp**2) 
           else
            fwtece0(k)=0.
         endif
3036    continue
3069  continue
      if (idebug.ge.3) write (6,*) 'GETTIR, ecebit = ', &
          (ecebit(i),i=1,nece)

!
      DEALLOCATE(rrgrid,bfield,rrout,bout,babs,bbb,ccc,ddd,btttt, &
          dsidr,ddsiddr,bx,ry,bbk)
!
      return
      end
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
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      include 'basiscomdu.f90'
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
                aa6gam(nmtark),aa7gam(nmtark)
      dimension tlibim(libim),slibim(libim),rrrlib(libim)
      dimension devxmpin(magpri),rnavxmpin(magpri) &
               ,devpsiin(nsilop),rnavpsiin(nsilop) &
               ,devfcin(nfcoil),rnavfcin(nfcoil) &
               ,devein(nesum),rnavecin(nesum)
      character*82 snap_ext
      real pefitktime 
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
           ,appendsnap,vbit,nbdrymx,efitversion
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
            dtmsefull,mse_strict,t_max_beam_off,ok_30rt,ok_210lt
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
           ,ok_30rt,ok_210lt,vbit,nbdrymx
      namelist/efitink/isetfb,ioffr,ioffz,ishiftz,gain,gainp,idplace &
           ,symmetrize,backaverage,lring
      data mcontr/35/,lfile/36/,ifpsi/0/
      data currn1/0.0/,currc79/0.0/,currc139/0.0/,currc199/0.0/ &
                      ,curriu30/0.0/,curriu90/0.0/,curriu150/0.0/ &
                      ,curril30/0.0/,curril90/0.0/,curril150/0.0/
      logical exists
! MPI >>>
      integer, intent(inout) :: kerror
      kerror = 0
! MPI <<<

      radeg=pi/180.
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
      twopi=2.0*pi
      tmu0=twopi*tmu
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
!--   Changed kdata <= 7 to kdata < 7				     --
!--   Snap-Extension mode = 7					     --
!----------------------------------------------------------------------
   10 if (kdata.ge.5.and.kdata.lt.7) go to 3000
      if (kdata.eq.8) then 
      CALL system("/u/huangyao/EFIT-PEFIT/pefit_read_kfile/pre_read")   
      go to 3000
      endif
      if (kdata.eq. 9) go to 3001
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
      betap0=0.50
      brsp(1)=-1.e+20
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
      error=1.0e-03
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
       sigrbd(i)=1.e10
       sigzbd(i)=1.e10
      enddo
      fli=0.0
      fwtbp=0.0
      fwtdlc=0.0
      fwtqa=0.0
      gammap=1.0e+10
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
      qenp=0.95
      qvfit=0.95
      scrape=0.030
      serror=0.03
      sidif=-1.0e+10
      symmetrize=.false.
      xltype=0.0
      xltype_180=0.
      gammap=1./gammap
      gammaf=gammap
      rmaxis=rzero
      mtear=0
      snapfile='none'
      nsnapf=66
! -- Qilong Ren
      write_Kfile = .false.
      fitfcsum = .false.
!----------------------------------------------------------------------
!--   Snap-Extension mode					     --
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
        calpa(1,1)=0.1
        calpa(2,1)=0.1
        calpa(3,1)=0.1
        xalpa(1)=0.0
        kcgama=1
        cgama(1,1)=0.1
        cgama(2,1)=0.1
        cgama(3,1)=0.1
        xgama(1)=0.0
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
        open(unit=neqdsk,status='old',         err=28002, &
          file='wakeefit.dat'                       )
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
! MPI >>>
#if defined(USEMPI)
! NOTE : ALL processes reading same input file so no need to broadcast value of ISHOT
          call mpi_stop
#else
          stop
#endif
! MPI <<<
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
        else
          if (kdata.ne.2) then
            table_di2 = table_di2(1:ltbdi2)//'156014/'
          else
            if (efitversion <= 20140331) then
               table_di2 = table_di2(1:ltbdi2)//'112000/'
            else
               table_di2 = table_di2(1:ltbdi2)//'156014/'
            endif
          endif
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
      if (nproc == 1) then
        call getpts(ishot,times,delt,ktime,istop)
      else
        call getpts_mpi(ishot,times,delt,ktime,istop)
      endif
#else
      call getpts(ishot,times,delt,ktime,istop)
#endif
      if (istop.gt.0) then
        if (rank == 0) then
          write (6,20000)
        endif
! MPI >>>
#if defined(USEMPI)
        call mpi_stop
#else
        stop
#endif
! MPI <<<
      endif
!
      mmstark=0
      do 142 i=1,nstark
        swtgam(i)=fwtgam(i)
        if (fwtgam(i).gt.1.e-06) mmstark=mmstark+1
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
! MPI >>>
#if defined(USEMPI)
      if (kerror /= 0) then
        kerror = 1
        return
      endif
 3400 continue
 ! NOTE : Finished EFIT run so STOP execution
 call mpi_stop
#else
 3400 if (kdata.eq.5) stop
!----------------------------------------------------------------------
!--   PEFIT mode = 8 and 9                                                --
!--   Make system call to drive PEFIT on Linux GPU                   --
!----------------------------------------------------------------------
      if (kdata.eq.8) then
      pefitktime = ktime    
      open(unit=99,file='sliceno.dat',status='unknown')
      write (99,"(e15.7)") pefitktime 
      close (unit=99)
      CALL system("/u/huangyao/EFIT-PEFIT/pefit_read_kfile/read")
      CALL system("/u/huangyao/EFIT-PEFIT/pefit_257_spline/bin/test")
      stop ('P-EFIT DONE!')
      endif

 3001 continue
      if (kdata.eq.9) then
      CALL system("/u/huangyao/EFIT-PEFIT/pefit_read_kfile/pre_read")
      open(unit=99,file='sliceno.dat',status='unknown')
      write (99,"(e15.7)") 1.0 
      close (unit=99)
      CALL system("/u/huangyao/EFIT-PEFIT/pefit_read_kfile/read_9")
      CALL system("/u/huangyao/EFIT-PEFIT/pefit_257_spline/bin/test")
      stop ('P-EFIT DONE!')
      endif
#endif
! MPI <<<
!
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
               ' 8=pefit_snap, 9=pefit_kfile):')
 6040 format (/,1x,'type shot #, start time(ms), time step(ms), steps' &
        ,'(<401):')
 6080 format (/,1x,'type limiter position (cm, 0=ptdata):')
 6090 format(' enter number of extra field lines to trace:')
 6091 format(' enter scrape off depth(m),'/ &
       '       sense of tracing (+1 for down, -1 for up),'/ &
       '       ixstrt (+1 for start on outside, -1' &
       ' for start inside):')
 6100 format(/,1x,48htype plot mode (0=none, 1=tektronix, 2=versatec, &
           ,17h 3=qms, -=x ray):)
 6200 format (/,1x,22hnumber of time slices?)
 6220 format (/,1x,22htype input file names:)
 6230 format (1x,1h#)
 6240 format (a)
 6600 format (/,1x,'good shot list file name ( 0=tty) ?')
 6610 format (/,1x,'command file name ( 0=none) ?')
 6617 format (/,1x,'type snap file extension (def for default):')
 6620 format (a)
 6700 format (a1,a12)
20000 format (/,1x,'shot data not on disk')
30000 format (i9)
30200 format (10f3.0)
      end
      subroutine fixstark(jtime,kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          fixstark adjusts the internal pitch angles              **
!**          based on spatial averaging data                         **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          01/01/07..........first created                         **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: ct,wkt,bkrt,bkzt
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
      include 'basiscomdu.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      dimension pds(6)
      real*8,dimension(:),allocatable :: bt,br,bzt,bz,bwork, &
             cwork,dwork
      real*4 save_gam(ntime,nstark)
      real*4 save_tangam(ntime,nstark)
      common/cworkbt/lkrt,lkzt
      data ifirst /0/
      save ifirst,save_gam,save_tangam
! MPI >>>
      integer, intent(inout) :: kerror
      kerror = 0
! MPI <<<
!
      ALLOCATE(bt(nwnh),br(nwnh),bzt(nwnh),bz(nwnh), &
         bwork(nw),cwork(nw),dwork(nw))
!
      if (keecur .gt. 0 .and. ifirst .eq. 0) then
	write(6,*) "Spatial averaging correction of MSE data"
	write(6,*) "not supported with Er fit at this time."
	write(6,*) "Going ahead with correction anyway,at this"
	write(6,*) "time, on channels requested."
!	write(6,*) "Calculating but not applying correction."
!	where(mse_spave_on .ne. 0) mse_spave_on = -1
      endif

      ltime = jtime
      if (jtime .lt. 0) then
	 ltime = -jtime
!	 do ichan=1,nstark
!		tangam(ltime,ichan) = save_tangam(ltime,ichan) 
!	 enddo
      endif

      if ( ifirst .eq. 0) then
	 ifirst = 1
         write(6,*)"Calculate pitch angle corrections", &
      	 " using spatial averaging"

	 do ichan=1,nstark
	    save_gam(ltime,ichan) = atan(tangam(ltime,ichan))
	    save_tangam(ltime,ichan) = tangam(ltime,ichan)
	 enddo
 	 return
      endif

      ip_sign = - cpasma(ltime) / abs(cpasma(ltime))

       call sets2d(psi,ct,rgrid,nw,bkrt,lkrt,zgrid,nh,bkzt, &
      		 lkzt,wkt,ier)
  
      if (pasmat(ltime).gt.0.0) then
        ssimag=-simag
        ssibry=-psibry
      else
        ssimag=simag
        ssibry=psibry
      endif


      if (jtime .gt. 0.0) then
!---------------------------------------------------------------------
!-- set up P' and FF', then integration                             --
!-- ffprim = (RBt) * d/dpsi(RBt)                                    --
!---------------------------------------------------------------------
      if (icurrt.ne.1) go to 7540
      pprime(1)=cratio*sbeta/darea/srma
      ffprim(1)=cratio*srma*2.*salpha/darea*twopi*tmu
      pprime(nw)=pprime(1)
      ffprim(nw)=ffprim(1)
 7540 continue
      if (icurrt.ne.2.and.icurrt.ne.5) go to 7550
      pprime(nw)=ppcurr(x111,kppcur)/darea
      ffprim(nw)=fpcurr(x111,kffcur)/darea*twopi*tmu
      pprime(1)=ppcurr(x000,kppcur)/darea
      ffprim(1)=fpcurr(x000,kffcur)/darea*twopi*tmu
      if (kfffnc.eq.8) then
         ffprec(nw)=fpecrr(x111,kffcur)/darea*twopi*tmu
         ffprec(1)=fpecrr(x000,kffcur)/darea*twopi*tmu
      else
         ffprec(nw)=0.0
         ffprec(1)=0.0
      endif
 7550 continue
      if (icurrt.ne.4) go to 7600
      call currnt(n222,iges,n222,n222,kerror)
! MPI >>>
#if defined(USEMPI)
      ! NOTE : Serial codes does NOT have this error check and does NOT return KERROR
      !if (kerror /= 0) then
      !  return
      !endif
#endif
! MPI <<<
      pprime(1)=cratio/darea/rzero
      ffprim(1)=rbetap*cratio*rzero*twopi*tmu/darea
      ffprim(nw)=ffprim(1)*gammaf
      pprime(nw)=pprime(1)*gammap
 7600 continue
!
      do i=2,nw-1
        ii=nw-i+1
        siii=1.0-1./float(nw-1)*(i-1)
        sigrid(ii)=siii
        if (icurrt.ne.2.and.icurrt.ne.5) go to 7792
        pprime(ii)=ppcurr(siii,kppcur)/darea
        ffprim(ii)=fpcurr(siii,kffcur)/darea*twopi*tmu
         if (kfffnc.eq.8) then
           ffprec(ii)=fpecrr(siii,kffcur)/darea*twopi*tmu
         else
           ffprec(ii)=0.0
         endif
 7792   continue
        if (icurrt.ne.4) go to 7794
        pprime(ii)=(1.-siii**enp)**emp*(1.-gammap)+gammap
        ffprim(ii)=ffprim(1)*pprime(ii)
        pprime(ii)=pprime(1)*pprime(ii)
 7794   continue
        if (icurrt.ne.1) go to 7796
        pprime(ii)=pprime(1)
        ffprim(ii)=ffprim(1)
 7796   continue
      enddo

      endif
      fpol(nw)=fbrdy*tmu
      sumf=fpol(nw)**2/2.
      delsi=-(psibry+psimag)/float(nw-1)
      do i=1,nw-1
        sumf=sumf+0.5*delsi*(ffprim(nw-i+1)+ffprim(nw-i))
        if(sumf .ge. 0.0) then
          fpol(nw-i)=sqrt(2.*sumf)*fpol(nw)/abs(fpol(nw))
        else
          fpol(nw-i)=fpol(nw)
        endif
      enddo

      call zpline(nw,sigrid,fpol,bwork,cwork,dwork)


      do ichan = 1,nmtark
      if (mse_spave_on(ichan) .ne. 0) then
!	    print *,i,ichan,cerer(1),e1rbz(ichan,1),e2rbz(ichan,1),
!     .	    e3rbr(ichan,1)
!	   print *,spatial_avg_gam(ichan,7,2,3)
!	   print *,spatial_avg_gam(ichan,8,2,3)
!	   print *,spatial_avg_gam(ichan,9,2,3)

!
!   spatial_avg_gam(chan, var [1 = r, 2 = z, 3-9 = A1-7], 5,3)
!
         ttl = 0.0
 	 rl = rrgam(ltime,ichan)
 	 zl = zzgam(ltime,ichan)
         call seva2d(bkrt,lkrt,bkzt,lkzt,ct,rl,zl,pds,ier,n333)
	 brl = -pds(3) / rl
  	 bzl = pds(2) / rl
	 psi_norm = (ssimag -pds(1)/ip_sign)/(ssimag-ssibry)
         btl = seval(nw,abs(psi_norm),sigrid,fpol,bwork, &
       	                              cwork,dwork) / rl
 	 tglocal = (bzl * a1gam(ltime,ichan)) /  &
      	 (btl * a2gam(ltime,ichan) + brl * a3gam(ltime,ichan) &
            + bzl * a4gam(ltime,ichan))


	 do i = 1,ngam_u
	 do j = 1,ngam_w
	    rl = spatial_avg_gam(ichan,1,i,j)
	    zl = spatial_avg_gam(ichan,2,i,j)
            call seva2d(bkrt,lkrt,bkzt,lkzt,ct,rl,zl,pds,ier,n333)
  	    brl = -pds(3) / rl
  	    bzl = pds(2) / rl
	    psi_norm = (ssimag -pds(1)/ip_sign)/(ssimag-ssibry)
            btl = seval(nw,abs(psi_norm),sigrid,fpol,bwork, &
       	                                 cwork,dwork) / rl
            tl = 0.0
	    tl = tl + spatial_avg_gam(ichan,4,i,j) * btl
	    tl = tl + spatial_avg_gam(ichan,5,i,j) * brl
	    tl = tl + spatial_avg_gam(ichan,6,i,j) * bzl
	    tl = spatial_avg_gam(ichan,3,i,j) * bzl / tl
	    ttl = ttl + tl
! 	    if(jtime .lt. 0)
!     .	    write(7,'(I2,6F13.8)') ichan,rl,zl,btl,brl,bzl,tl
	 enddo
	 enddo
!            spatial_fix(ichan,ltime) = 
!    .  	   atan(cmgam(ichan,ltime)) - atan(ttl)  
             spatial_fix(ichan,ltime) =  &
        	   atan(tglocal) - atan(ttl)  

          if (jtime.gt.0.and.mse_spave_on(ichan) .eq. 1) then
  	     tangam(ltime,ichan) = tan(save_gam(ltime,ichan) &
      	      - spatial_fix(ichan,ltime)) 
          endif
	 

      endif

      
      enddo
          if (jtime .lt. 0) then
	      ifirst = 0
          endif
!
      DEALLOCATE(bt,br,bzt,bz,bwork,cwork,dwork)
!
      return
      end


      subroutine getsigma(jtimex,niterax)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          GETSIGMA is the control for getting the uncertainty     **
!**          in Magnetic Data                                        **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**          jtime - time slice number                               **
!**          nitera - iteration number                               **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**        2002/11/03..........First Created   EKS                   **
!**        2006/01/19..........Updated                               **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      common/cwork3/lkx,lky
      integer jtimex,niterax
      real*8             :: gradsdr,gradsdz,brdr,brdz,bzdr,bzdz,cost,sint &
                          ,oldfit
      dimension pds(6)
!----------------------------------------------------------------------
!--    BR=-1/R dpsi/dZ           BZ=1/R dpsi/dR                      --
!--            Bprobe= BR cost + BZ sint                             --
!--            gradsdr= dBRdR cost + dBZdR sint         dBprobe/dR   --
!--            gradsdz= dBRdZ cost + dBZdZ sint         dBrpobe/dZ   --
!--            gradsmp= sqrt (gradsdr**2 + gradsdz**2)               --
!--            gradsfl= |flux loop gradient|                         --
!--            bpermp= -BR sint + BZ cost                            --
!----------------------------------------------------------------------
      if (jtimex*niterax.eq.1) then
      open(unit=33,file='getsigma.log',form='formatted',status='unknown' &
            )
      open(unit=34,file='mprobe.err',form='formatted',status='unknown' &
            )
      open(unit=35,file='siloop.err',form='formatted',status='unknown' &
            )
      open(unit=36,file='fcoil.err',form='formatted',status='unknown' &
            )
      open(unit=37,file='ecoil.err',form='formatted',status='unknown' &
            )
      open(unit=38,file='bandcur.err',form='formatted',status='unknown' &
            )
      endif
      do i=1,magpri
         call seva2d(bkx,lkx,bky,lky,c,xmp2(i),ymp2(i),pds,ier,n666)
!----------------------------------------------------------------------
!-- Calculate dBR/dr, dBR/dz, dBZ/dr, dBZ/dz                         --
!----------------------------------------------------------------------
         brdr=(-pds(4)+pds(3)/xmp2(i))/xmp2(i)
         brdz=(-pds(6))/xmp2(i)
         bzdr=(pds(5)-(pds(2)/xmp2(i)))/xmp2(i)
         bzdz=pds(4)/xmp2(i)
!----------------------------------------------------------------------
!--  Form dBprobe/dR and dBprobe/dZ, then gradient                   --
!----------------------------------------------------------------------
         sinm=sin(radeg*amp2(i))
         cosm=cos(radeg*amp2(i))
         gradsdr=brdr*cosm + bzdr*sinm
         gradsdz=brdz*cosm + bzdz*sinm
         gradsmp(jtimex,i)=sqrt(gradsdr**2+gradsdz**2)
!----------------------------------------------------------------------
!-- Calculate B perpendicular to magnetic probe                      --
!----------------------------------------------------------------------
         bpermp(jtimex,i)=(pds(2)*cosm + pds(3)*sinm) &
                             /xmp2(i)
      enddo
!----------------------------------------------------------------------
!-- calc gradient of flux loops                                      --
!----------------------------------------------------------------------
      do i=1,nsilop
         call seva2d(bkx,lkx,bky,lky,c,rsi(i),zsi(i),pds,ier,n333)
         gradsfl(jtimex,i)=sqrt(pds(2)**2+pds(3)**2)
      enddo
      write(33,*) '#', ishot,time(jtimex),jtimex
      write(33,*)'#magprobe no.,gradsmp,bpermp'
      do i=1,magpri
         write(33,'(i5,2e13.5)')i,gradsmp(jtimex,i),bpermp(jtimex,i)
      enddo
      write(33,*)'#fluxloop no., gradsfl'
      do i=1,nsilop
         write(33,'(i5,2e13.5)') i,gradsfl(jtimex,i)
      enddo
      call magsigma(ishot,time(jtimex),jtimex,gradsmp,gradsfl, &
           bpermp,sigmaf,sigmab,sigmae,sigmaip,sigmafl,sigmamp)
!----------------------------------------------------------------------
!-- Set fitting weights                                              --
!----------------------------------------------------------------------
      do i=33,38
         write(i,*) '#', ishot,time(jtimex),jtimex
         write(i,*) '#errorm=', errorm
         write(i,*) '#errmag=', errmag
      enddo
      write(38,*) '#sigmab'
      write(38,'(e13.5)') sigmab(jtimex)
      write(38,*) '#sigmaip0,sigmaip'
      write(38, '(2e13.5)') sigmaip0,sigmaip(jtimex)
      write(34,*) '#sigmamp0,sigmamp'
      do i=1,magpri
         write(34,'(i5,2e13.5)') i,sigmamp0(i),sigmamp(jtimex,i)
      end do
      write(34,*) ' '
      write(34,*) '#t0mp,mp_k, mprcg, vresmp, devmp'
      do i=1,magpri
         write(34,'(6e13.5)') t0xmp(i),xmp_k(i),xmprcg(i),vresxmp(i), &
                 devxmp(jtimex,i),rnavxmp(jtimex,i)
      enddo
      write(34,*) ' '
      write(36,*) '#sigmaf0,sigmaf'
      do i=1,nfcoil
         write(36,'(i5,2e13.5)') i,sigmaf0(i),sigmaf(jtimex,i)
      end do
      write(36,*) ' '
      write(37,*) '#sigmae0,sigmae'
      do i=1,nesum
         write(37,'(i5, 2e13.5)') i, sigmae0(i),sigmae(jtimex,i)
      end do
      write(37,*) ' '
      write(35,*) '#sigmafl0,sigmafl'
      do i=1,nsilop
         write(35,'(i5,2e13.5)') i, sigmafl0(i),sigmafl(jtimex,i)
      enddo
      write(35,*) ' '
      write(35,*) '#t0psi,psi_k, psircg, vrespsi, devpsi'
      do i=1,nsilop
         write(35,'(6e13.5)') t0psi(i),psi_k(i),psircg(i),vrespsi(i), &
               devpsi(jtimex,i),rnavpsi(jtimex,i)
      enddo
      write(35,*) ' '
      do i=1,nsilop
         oldfit = fwtsi(i)
         if (sigmafl(jtimex,i)/=0.0) then
            fwtsi(i)=swtsi(i)/sigmafl(jtimex,i)
         else
            fwtsi(i)=0.0
         endif
         write (99,*) i, swtsi(i), oldfit, fwtsi(i)
      enddo
      do i=1,magpri
         oldfit = fwtmp2(i)
         if (sigmamp(jtimex,i)/=0.0) then
            fwtmp2(i)=swtmp2(i)/sigmamp(jtimex,i)
         else
            fwtmp2(i)=0.0
         endif
         write (99,*) i, swtmp2(i), oldfit, fwtmp2(i)
      enddo
      do i=1,nesum
         oldfit = fwtec(i)
         if (sigmae(jtimex,i)/=0.0) then
           fwtec(i)=swtec(i)/sigmae(jtimex,i)
         else
           fwtec(i)=0.0
         endif
         write (99,*) i, swtec(i), oldfit, fwtec(i)
      enddo
      do i=1,nfcoil
         oldfit = fwtfc(i)
         if (sigmaf(jtimex,i)/=0.0) then
           fwtfc(i)=swtfc(i)/sigmaf(jtimex,i)
         else
           fwtfc(i)=0.0
         endif
         write (99,*) i, swtfc(i), oldfit, fwtfc(i)
      enddo
      oldfit = fwtcur
      if (sigmaip(jtimex)/=0.0) then
           fwtcur=swtcur/sigmaip(jtimex)
      else
           fwtcur=0.0
      endif
      write (99,*) swtcur, oldfit, fwtcur
!
      return
      end
      subroutine getstark(ktime)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          getstark obtains the internal pitch angles              **
!**          from polarimetry measurement using Wroblewski's routine **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          23/03/90..........first created                         **
!**          93/04/23..........revised for double precision version  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      real*4 avem,tanham(ktime,nmtark),sigham(ktime,nmtark), &
         rrham(nmtark),zzham(nmtark), &
         sarkar,sarkaz,a1ham(nmtark), &
         a2ham(nmtark),a3ham(nmtark),a4ham(nmtark), &
         a5ham(nmtark),a6ham(nmtark),a7ham(nmtark),atime(ktime), &
         spatial_avg_ham(nmtark,ngam_vars,ngam_u,ngam_w), &
         hgain(nmtark),hslope(nmtark),hscale(nmtark), &
         hoffset(nmtark),max_beamOff
!
      do 10 i=1,ktime
        atime(i)=time(i)
   10 continue
      if(dtmsefull .gt. 0.0) then
         avem=dtmsefull / 1000.0
      else
         avem=2.0*iavem / 1000.0
      endif
       max_beamOff = t_max_beam_off / 1000.0
       call  set_mse_beam_logic(mse_strict,max_beamOff, &
                           ok_210lt,ok_30rt)
       tanham = 0.0
       sigham = 0.0
       call  stark2(ishot,atime,ktime,avem,msefitfun,tanham,sigham, &
             rrham,zzham,a1ham,a2ham,a3ham,a4ham,a5ham,a6ham,a7ham, &
             iergam,msebkp,mse_quiet)
       call get_mse_spatial_data(spatial_avg_ham)
       call get_mse_calibration(msefitfun,hgain, &
                    hslope,hscale,hoffset)
      kfixstark = 0
      do 200 n=1,nmtark
	rmse_gain(n) = hgain(n) 
	rmse_slope(n) = hslope(n) 
	rmse_scale(n) = hscale(n) 
	rmse_offset(n) = hoffset(n) 
	if(mse_spave_on(n) .ne. 0) kfixstark = 1
	do i=1,ngam_vars
	do j=1,ngam_u
	do k=1,ngam_w
	    spatial_avg_gam(n,i,j,k)= spatial_avg_ham(n,i,j,k)
	enddo
	enddo
	enddo
!
      do 100 i=1,ktime
        tangam(i,n)=tanham(i,n)
        siggam(i,n)=sigham(i,n)
        rrgam(i,n)=rrham(n)
        zzgam(i,n)=zzham(n)
        starkar(i,n)=sarkar
        starkaz(i,n)=sarkaz
        a1gam(i,n)=a1ham(n)
        a2gam(i,n)=a2ham(n)
        a3gam(i,n)=a3ham(n)
        a4gam(i,n)=a4ham(n)
        a5gam(i,n)=a5ham(n)
        a6gam(i,n)=a6ham(n)
        a7gam(i,n)=a7ham(n)
        a8gam(i,n)=0.0
        if (abs(tangam(i,n)).le.1.e-10.and. &
            abs(siggam(i,n)).le.1.e-10)then
            fwtgam(n)=0.0
            siggam(i,n)=0.0
        else if (abs(tangam(i,n)).le.1.e-10.and. &
                abs(siggam(i,n)).le.100.0) then
            fwtgam(n)=0.0
            siggam(i,n)=0.0
        else if (iergam(n).gt.0) then
            fwtgam(n)=0.0
            siggam(i,n)=0.0
       endif
  100 continue
  200 continue
!
! kfixstark is zero when none of the individual channels is turned on
! in this case set kwaitmse to zero which turns the mse spacial average
! off globally
!
      if (kfixstark .eq. 0) then
		kwaitmse = 0
      elseif (kwaitmse .eq. 0) then
		kwaitmse = 5
      endif
      return
      end
      subroutine gette(kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          GETTE gets the electron temperature                     **
!**          profiles.                                               **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          15/09/87..........first created                         **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      common/cwork2/arsp(ndata,nppcur),wrsp(nppcur),work(ndata), &
                    bdata(ndata),ematrix(nppcur,nppcur), &
                    einv(nppcur,nppcur)
! MPI >>>
      integer, intent(inout) :: kerror
      kerror = 0
! MPI <<<
!----------------------------------------------------------------------
!--  singular decomposition                                          --
!----------------------------------------------------------------------
      do 1100 nj=1,npress
        do 1000 nk=1,nptef
          xn=-rpress(nj)
          arsp(nj,nk)=xn**(nk-1)/sgteth(nj)
 1000   continue
        bdata(nj)=tethom(nj)/sgteth(nj)
 1100 continue
      ntedat=npress
      if (cstabte.gt.0.0) then
        nj=npress
        do 1120 jj=ncstte,nptef
         nj=nj+1
         do 1110 nk=1,nptef
          arsp(nj,nk)=0.0
          if (jj.eq.nk) arsp(nj,nk)=cstabte
 1110    continue
         bdata(nj)=0.0
 1120   continue
        ntedat=ntedat+nptef-ncstte+1
      endif
!---------------------------------------------------------------------
!-- form error matrix                                               --
!---------------------------------------------------------------------
      do 1900 i=1,nptef
      do 1900 j=1,nptef
        ematrix(i,j)=0.0
        do 1900 k=1,npress
          ematrix(i,j)=ematrix(i,j)+arsp(k,i)*arsp(k,j)
 1900   continue
!
      call sdecm(arsp,ndata,ntedat,nptef,bdata,ntedat,n111,wrsp, &
                 work,ier)
      if (ier.ne.129) go to 1200
      write (nttyo,8000) ier
! MPI >>>
#if defined(USEMPI)
      ! ERROR_FIX >>>
      !kerror = 1
      !return
      ! <<< >>>
      call mpi_stop
      ! ERROR_FIX <<<
#else
      stop
#endif
! MPI <<<
!
 1200 continue
      cond=ier
      toler=1.0e-06*wrsp(1)
      do 1600 i=1,nptef
        t=0.0
        if (wrsp(i).gt.toler) t=bdata(i)/wrsp(i)
        work(i)=t
 1600 continue
      do 1650 i=1,nptef
        tefit(i)=0.0
        do 1650 j=1,nptef
          tefit(i)=tefit(i)+arsp(i,j)*work(j)
 1650   continue
!------------------------------------------------------------------
!-- compute chi square                                           --
!------------------------------------------------------------------
      chisqte=0.0
      do 1700 i=1,npress
        tenow=0.0
        xn=-rpress(i)
        do 1670 j=1,nptef
 1670     tenow=tenow+tefit(j)*xn**(j-1)
        chisqte=chisqte+((tenow-tethom(i))/sgteth(i))**2
 1700 continue
!---------------------------------------------------------------------
!-- get inverse of error matrix                                     --
!---------------------------------------------------------------------
      call linv1f(ematrix,nptef,nppcur,einv,n444,work,ier)
!----------------------------------------------------------------------
!--  boundary values                                                 --
!----------------------------------------------------------------------
        tebdry=0.0
        stebdry=0.0
        tepbry=0.0
        sigtepb=0.0
        do 1850 j=1,nptef
          do 1840 i=1,nptef
            stebdry=stebdry+einv(i,j)
            sigtepb=sigtepb+(i-1)*(j-1)*einv(i,j)
 1840     continue
          tepbry=tepbry+(j-1)*tefit(j)
 1850     tebdry=tebdry+tefit(j)
        if (stebdry.gt.0.0) then
          stebdry=sqrt(stebdry)
        else
          stebdry=tebdry
        endif
        if (sigtepb.gt.0.0) then
          sigtepb=sqrt(sigtepb)
        else
          sigtepb=tepbry
        endif
!
      return
 8000 format (/,'  ** Problem in Decomposition **',i10)
      end
      subroutine gettion(kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          GETTION gets the ion temperature profile.               **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          15/09/87..........first created                         **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      common/cwork2/arsp(ndata,nppcur),wrsp(nppcur),work(ndata), &
                    bdata(ndata),ematrix(nppcur,nppcur), &
                    einv(nppcur,nppcur)
      common/cwork3/lkx,lky
      dimension pds(6),bwork(ndata),cwork(ndata),dwork(ndata)
! MPI >>>
      integer, intent(inout) :: kerror
! MPI <<<
      if (rion(2).lt.0.0) then
      do 2900 i=1,nption
       xsiion(i)=-rion(i)
       sigti(i)=sgtimin*tionex(i)
 2900 continue
      call zpline(nption,xsiion,tionex,bwork,cwork,dwork)
      do 3000 i=1,npress
        xn=-rpress(i)
        tithom(i)=seval(nption,xn,xsiion,tionex,bwork,cwork,dwork)
        stitho(i)=sgtimin*tithom(i)
 3000 continue
      tibdry=seval(nption,x111,xsiion,tionex,bwork,cwork,dwork)
      tipbry=speval(nption,x111,xsiion,tionex,bwork,cwork,dwork)
      stibdry=sgtimin*tibdry
      sigtipb=sgtimin*tipbry
      return
      endif
!----------------------------------------------------------------------
!--  singular decomposition                                          --
!----------------------------------------------------------------------
      do 200 i=1,nption
        call seva2d(bkx,lkx,bky,lky,c,rion(i),zion(i),pds,ier,n111)
        xsiion(i)=(simag-pds(1))/sidif
  200 continue
      need=nption+1
      xsiion(need)=1.0
      tionex(need)=tebdry
      sigti(need)=stebdry
!
      do 1100 nj=1,need
        do 1000 nk=1,nptionf
          xn=xsiion(nj)
          arsp(nj,nk)=xn**(nk-1)/sigti(nj)
 1000   continue
        bdata(nj)=tionex(nj)/sigti(nj)
 1100 continue
!---------------------------------------------------------------------
!-- form error matrix                                               --
!---------------------------------------------------------------------
      do 1900 i=1,nptionf
      do 1900 j=1,nptionf
        ematrix(i,j)=0.0
        do 1900 k=1,need
          ematrix(i,j)=ematrix(i,j)+arsp(k,i)*arsp(k,j)
 1900   continue
!
      nnn=1
      call sdecm(arsp,ndata,need,nptionf,bdata,need,nnn,wrsp,work,ier)
      if (ier.ne.129) go to 1200
      write (nttyo,8000) ier
! MPI >>>
#if defined(USEMPI)
      ! ERROR_FIX >>>
      !kerror = 1
      !return
      ! <<< >>>
      call mpi_stop
      ! ERROR_FIX <<<
#else
      stop
#endif
! MPI <<<
!
 1200 continue
      cond=ier
      toler=1.0e-06*wrsp(1)
      do 1600 i=1,nptionf
        t=0.0
        if (wrsp(i).gt.toler) t=bdata(i)/wrsp(i)
        work(i)=t
 1600 continue
      do 1650 i=1,nptionf
        tifit(i)=0.0
        do 1650 j=1,nptionf
          tifit(i)=tifit(i)+arsp(i,j)*work(j)
 1650   continue
!------------------------------------------------------------------
!-- compute chi square                                           --
!------------------------------------------------------------------
      chisqti=0.0
      do 1700 i=1,need
        tinow=0.0
        do 1670 j=1,nptionf
 1670     tinow=tinow+tifit(j)*xsiion(i)**(j-1)
        chisqti=chisqti+((tinow-tionex(i))/sigti(i))**2
 1700 continue
!---------------------------------------------------------------------
!-- get inverse of error matrix                                     --
!---------------------------------------------------------------------
      call linv1f(ematrix,nptionf,nppcur,einv,n444,work,ier)
!---------------------------------------------------------------------
!--  project ion temperature into Thompson flux grid                --
!---------------------------------------------------------------------
      do 1800 i=1,npress
        tithom(i)=0.0
        stitho(i)=0.0
        xn=-rpress(i)
        do 1750 j=1,nptionf
          do 1740 k=1,nptionf
          stitho(i)=stitho(i)+einv(k,j)*xn**(j-1)*xn**(k-1)
 1740     continue
 1750     tithom(i)=tithom(i)+tifit(j)*xn**(j-1)
        if (stitho(i).gt.0.0) then
          stitho(i)=sqrt(stitho(i))
        else
          stitho(i)=0.5*tithom(i)
        endif
 1800 continue
!----------------------------------------------------------------------
!--  boundary values                                                 --
!----------------------------------------------------------------------
        tibdry=0.0
        stibdry=0.0
        tipbry=0.0
        sigtipb=0.0
        do 1850 j=1,nptionf
          do 1840 i=1,nptionf
            stibdry=stibdry+einv(i,j)
            sigtipb=sigtipb+(i-1)*(j-1)*einv(i,j)
 1840     continue
          tipbry=tipbry+(j-1)*tifit(j)
 1850     tibdry=tibdry+tifit(j)
        if (stibdry.gt.0.0) then
          stibdry=sqrt(stibdry)
        else
          stibdry=tibdry
        endif
        if (sigtipb.gt.0.0) then
          sigtipb=sqrt(sigtipb)
        else
          sigtipb=tipbry
        endif
!
      return
 8000 format (/,'  ** Problem in Decomposition **',i10)
      end
      subroutine green(ifag,jtime,niter)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          green set up the appropriate response functions for use **
!**          with the routine matrix.                                **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**          94/03/08..........revised                               **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      common/cwork3/lkx,lky
      dimension xpsii(nwcurn),xpsis(nwcurn),xpsisb(nwcurn)
      dimension xsier(nercur)
      dimension pds(6)
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
        if ((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.1.0)) go to 1582
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
            if (abs(pres0).gt.1.e-10) then
               pwop0=prew0/pres0
               pwp0r2=pwop0*rgrvt(i)
            else
               pwop0=0.0
               pwp0r2=0.0
            endif
            factor=factor*(1.-0.5*pwp0r2**2)
          elseif (kvtor.eq.3) then
            prew0=pwcurr(xpsi(kk),kwwcur)
            pres0=prcurr(xpsi(kk),kppcur)
            if (abs(pres0).gt.1.e-10) then
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
 1560     continue
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
 1860     continue
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
        if (kvtor.le.0) goto 1982
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
        if (abs(pres0).gt.1.e-10) then
           pwop0=prew0/pres0
           pwp0r2=pwop0*rmggvt
        else
           pwop0=0.0
           pwp0r2=0.0
        endif
        rmvtor=rmvtor*(1.-0.5*pwp0r2**2)
        rmvjj=rmvjj*(1.-0.5*pwp0r2**2)
      endif
      if (kvtor.eq.3) then
        prew0=pwcurr(xpzero,kwwcur)
        pres0=prcurr(xpzero,kppcur)
        if (abs(pres0).gt.1.e-10) then
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
        if (kvtor.le.0) goto 11982
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
      end
      subroutine inicur(ks)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          inicur initializes the current density distribution.    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**           ks..............time slice number                      **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
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
      go to (100,1100) iabs(icinit)
      return
  100 continue
      if (aelip.gt.0.0) go to 150
      aelip=0.50
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
      if (zelip.gt.1.e5.or.zelips.gt.1.e5) then
        zelip=1.447310*(expmpi(ks,37)-expmpi(ks,43)) &
             +0.692055*(expmpi(ks,57)-expmpi(ks,53)) &
             +0.728045*(silopt(ks,27)-silopt(ks,37)) &
             +2.047150*(silopt(ks,2) -silopt(ks,11))
        zelip=zelip*1.e6/pasmat(ks)
        zbound=zelip
        eelip=1.5
      endif
!----------------------------------------------------------------
!-- set zelip=0.0 if bad signals              96/06/24         --
!----------------------------------------------------------------
      if (abs(fwtmp2(37)).le.1.0e-30)  zelip=0.0
      if (abs(fwtmp2(43)).le.1.0e-30)  zelip=0.0
      if (abs(fwtmp2(57)).le.1.0e-30)  zelip=0.0
      if (abs(fwtmp2(53)).le.1.0e-30)  zelip=0.0
      if (abs(fwtsi(27)).le.1.0e-30)  zelip=0.0
      if (abs(fwtsi(37)).le.1.0e-30)  zelip=0.0
      if (abs(fwtsi(2)).le.1.0e-30)  zelip=0.0
      if (abs(fwtsi(11)).le.1.0e-30)  zelip=0.0
!
      do 1300 i=1,nw
      do 1300 j=1,nh
        kk=(i-1)*nh+j
        erho=sqrt((rgrid(i)-relip)**2+((zgrid(j)-zelip)/eelip)**2)
        xpsi(kk)=(erho/aelip)**2
 1300 continue
      return
      end
      subroutine matrix(jtime,iter,ichisq,nniter,kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          matrix calculates the appropriate response matrix and   **
!**          invert it to get the plasma current strengths.  note    **
!**          that npcurn=nffcur+nppcur, nrsmat=nfcoil+npcurn+        **
!**          number of constraints.                                  **
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
!**          94/03/08..........revised                               **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      include 'basiscomdu.f90'
      dimension arsp(nrsmat,mfnpcr),wrsp(mfnpcr)
      dimension brsold(nrsmat),work(nrsma2),vcurrto(nvesel)
      dimension xpsfp(nffcur),xpspp(nppcur),xpspwp(nwwcur)
      dimension crsp(4*(npcurn-2)+6+npcurn*npcurn,nrsmat)
      dimension b(nrsmat),z(4*(npcurn-2)+6+npcurn*npcurn)
      common/cwork3/lkx,lky
      dimension pds(6)
      dimension rxxx(ndata),rxxxf(ndata),rxx2(ndata),rxxw(ndata)
!---------------------------------------------------------------------
!--   relax saimin=50 from 30               04/27/90                --
!--                60 from 50               03/31/93                --
!---------------------------------------------------------------------
      data iupdat/0/,minite/8/,ten24/1.e4/,z04/1.0e-04/
      save z04
!
      integer, intent(inout) :: kerror
      kerror = 0
      if (iconvr.eq.3) go to 6000
!----------------------------------------------------------------------
!-- Variable fitdelz                                                 --
!----------------------------------------------------------------------
      if (fitdelz) scadelz=scaledz
!----------------------------------------------------------------------
!--  set up fitting weight for boundary constraints                  --
!----------------------------------------------------------------------
      if (nbdry.gt.0) then
        fwtbdr=abs(errbry)*max(abs(sidif),z04)
        fwtbdr=1.0/fwtbdr
        do i=1,nbdry
          if (sigrbd(i).lt.1.e10.and.sigzbd(i).lt.1.e10) then
          call seva2d(bkx,lkx,bky,lky,c,rbdry(i),zbdry(i),pds,ier,n666)
          fwtbdr=sqrt((sigrbd(i)*pds(2))**2+(sigzbd(i)*pds(3))**2)
          fwtbdr=1.0/fwtbdr
          endif
          fwtbry(i)=fwtbdr*fwtbdry(i)
        enddo
      endif
!----------------------------------------------------------------------
!-- set up the response matrix arsp                                  --
!----------------------------------------------------------------------
      ichisq=0
      nsq=2
      saiold=saisq
      do 1000 i=1,nrsmat
        brsold(i)=brsp(i)
 1000 continue
      if (ifitvs.gt.0) then
        do 1010 i=1,nvesel
          vcurrto(i)=vcurrt(i)
 1010   continue
      endif
!----------------------------------------------------------------------
!--  singular decomposition, first F-coil currents, set up arsp      --
!----------------------------------------------------------------------
      do 2100 nk=1,nfcoil
        nj=0
        do 2020 m=1,nsilop
          if (fwtsi(m).le.0.0) go to 2020
          nj=nj+1
          arsp(nj,nk)=fwtsi(m)*rsilfc(m,nk)
 2020   continue
        do 2040 m=1,magpri
          if (fwtmp2(m).le.0.0) go to 2040
          nj=nj+1
          arsp(nj,nk)=fwtmp2(m)*rmp2fc(m,nk)
 2040   continue
        do 2060 m=1,nstark
          if (fwtgam(m).le.0.0) go to 2060
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgamfc(m,nk)
 2060   continue
        do 2070 m=1,nece
          if (fwtece(m).le.0.0) go to 2070
          nj=nj+1
          arsp(nj,nk)=fwtece(m)*recefc(m,nk)
 2070   continue
          if (fwtecebz.le.0.0) go to 2075
          nj=nj+1
          arsp(nj,nk)=fwtecebz*recebzfc(nk)
 2075   continue
        if (fwtcur.le.0.0) go to 2080
        nj=nj+1
        arsp(nj,nk)=0.0
 2080   continue
        if (fwtqa.le.0.0) go to 2085
        nj=nj+1
        arsp(nj,nk)=0.0
 2085 continue
      if (fwtbp.le.0.0) go to 2090
      do 2087 m=2,kffcur
        nj=nj+1
        arsp(nj,nk)=0.0
 2087 continue
 2090 continue
      if (fwtdlc.le.0.0) go to 2092
      nj=nj+1
      arsp(nj,nk)=0.0
 2092 continue
 2096 continue
      do 2097 m=1,nfcoil
        if (fwtfc(m).le.0.0) go to 2097
        nj=nj+1
        arsp(nj,nk)=0.
        if (m.eq.nk) arsp(nj,nk)=fwtfc(m)
 2097 continue
!--------------------------------------------------------------------------
!--  pressure data                                                       --
!--------------------------------------------------------------------------
      if (kprfit.le.0.or.kdofit.eq.0) go to 2099
      if (npress.le.0) goto 2099
      do 2098 m=1,npress
        if (fwtpre(m).le.0.0) goto 2098
        nj=nj+1
        arsp(nj,nk)=0.0
 2098 continue
 2099 continue
      if (kpressb.eq.2) then
        nj=nj+1
        arsp(nj,nk)=0.0
      endif
!---------------------------------------------------------------------------
!-- rotational pressure data                                              --
!---------------------------------------------------------------------------
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          if (fwtprw(m).gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
          endif
        enddo
      endif
!---------------------------------------------------------------------------
!--  J(PSIWANT) constraint                                                --
!---------------------------------------------------------------------------
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
        arsp(nj,nk)=0.0
       enddo
      endif
!----------------------------------------------------------------------------
!-- P', FF', and rotational constraints                                    --
!----------------------------------------------------------------------------
      if (kcalpa.gt.0) then
        do 30050 j=1,kcalpa
          nj=nj+1
          arsp(nj,nk)=0.0
30050   continue
      endif
      if (kcgama.gt.0) then
        do 30060 j=1,kcgama
          nj=nj+1
          arsp(nj,nk)=0.0
30060   continue
      endif
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
!------------------------------------------------------------------------------
!-- Boundary constraints, data                                               --
!------------------------------------------------------------------------------
      if (nbdry.gt.0) then
        do j=1,nbdry
          if (fwtbdry(j).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=fwtbry(j)*rbdrfc(j,nk)
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!--  E coil currents, data                                                   --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
       do m=1,nesum
        if (fwtec(m).gt.0.0) then
        nj=nj+1
        arsp(nj,nk)=0.
        endif
       enddo
      endif
!------------------------------------------------------------------------------
!--  Reference flux, data                                                    --
!------------------------------------------------------------------------------
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
      endif
!------------------------------------------------------------------------
!-- Summation of F-coils currents
!------------------------------------------------------------------------
      if (fitfcsum) then
         nj = nj + 1
         arsp(nj,nk) = fwtfcsum(nk)
      endif
 2100 continue
!-----------------------------------------------------------------------
!--  plasma current P', FF', and Pw', set up response matrix arsp     --
!-----------------------------------------------------------------------
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        ysiwant=sizeroj(i)
        call setfp(ysiwant,xpsfp)
        call setpp(ysiwant,xpspp)
!----------------------------------------------------------------------
!--  local or flux surface averaged J constraint ?                   --
!----------------------------------------------------------------------
        if (rzeroj(i).gt.0.0) rxxx(i)=rzeroj(i)
        if (rzeroj(i).lt.0.0) rxxx(i)=rseps(1,jtime)/100.
        rxxxf(i)=rxxx(i)
        if (rzeroj(i).eq.0.0) then
             rxxx(i)=1./r1sdry(i)
             rxxxf(i)=r1sdry(i)/r2sdry(i)
        endif
        if (rxxx(i).le.0.0) then
            rxxx(i)=rcentr
            rxxxf(i)=rxxx(i)
        endif
        brspmin=max(ten24,abs(brsp(nfcoil+1)))
        fwtxxzj=fwtxxj*1000./brspmin
        if (kvtor.gt.0) then
          call setpwp(ysiwant,xpspwp)
          if (nniter.gt.0) then
          if (kvtor.eq.2.or.kvtor.eq.3) then
            prew0=pwcurr(ysiwant,kwwcur)
            pres0=prcurr(ysiwant,kppcur)
            if (abs(pres0).gt.1.e-10) then
               pwop0=prew0/pres0
            else
               pwop0=0.0
            endif
          endif
          endif
          if (rzeroj(i).gt.0.0) then
              rxx2(i)=(rzeroj(i)/rvtor)**2-1.
              rxxw(i)=rxx2(i)*rzeroj(i)
              if (nniter.gt.0) then
              if (kvtor.eq.2) then
                rxxw(i)=rxxw(i)*(1.+pwop0*rxx2(i))
                rxxx(i)=rxxx(i)*(1.-0.5*(pwop0*rxx2(i))**2)
              endif
              if (kvtor.eq.3) then
                pwp0r2=pwop0*rxx2(i)
                ptop0=exp(pwp0r2)
                rxxw(i)=rxxw(i)*ptop0
                rxxx(i)=rxxx(i)*ptop0*(1.-pwp0r2)
              endif
              endif
          endif
          if (rzeroj(i).lt.0.0) then
              rxxw(i)=rseps(1,jtime)/100.
              rxx2(i)=(rxxw(i)/rvtor)**2-1.
              rxxw(i)=rxx2(i)*rxxw(i)
              if (nniter.gt.0) then
              if (kvtor.eq.2) then
                rxxw(i)=rxxw(i)*(1.+pwop0*rxx2(i))
                rxxx(i)=rxxx(i)*(1.-0.5*(pwop0*rxx2(i))**2)
              endif
              if (kvtor.eq.3) then
                pwp0r2=pwop0*rxx2(i)
                ptop0=exp(pwp0r2)
                rxxw(i)=rxxw(i)*ptop0
                rxxx(i)=rxxx(i)*ptop0*(1.-pwp0r2)
              endif
              endif
          endif
          if (rzeroj(i).eq.0.0) then
              rxx2(i)=  r2wdry/rvtor**2-1.
              rxxw(i)=  rxx2(i)/r1sdry(i)
              if (nniter.gt.0) then
              if (kvtor.eq.2) then
                rxxw(i)=rxxw(i)+pwop0*r4wdry/r1sdry(i)
                rxxx(i)=rxxx(i)-0.5*pwop0**2*r4wdry/r1sdry(i)
              endif
              if (kvtor.eq.3) then
                rxxx(i)=(rpwdry-pwop0*rp2wdry)/r1sdry(i)
                rxxw(i)=rp2wdry/r1sdry(i)
              endif
              endif
          endif
        endif
       enddo
      endif
!----------------------------------------------------------------------
!--  start loop for plasma parameters: P', FF', Pw'                  --
!----------------------------------------------------------------------
      do 2210 nk=nfcoil+1,nfnwcr
        n=nk-nfcoil
        nj=0
        do 2120 m=1,nsilop
          if (fwtsi(m).le.0.0) go to 2120
          nj=nj+1
          arsp(nj,nk)=fwtsi(m)*rsilpc(m,n)
 2120   continue
        do 2140 m=1,magpri
          if (fwtmp2(m).le.0.0) go to 2140
          nj=nj+1
          arsp(nj,nk)=fwtmp2(m)*rmp2pc(m,n)
 2140   continue
        do 2160 m=1,nstark
          if (fwtgam(m).le.0.0) go to 2160
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgampc(m,n)
 2160   continue
        do 2170 m=1,nece
          if (fwtece(m).le.0.0) go to 2170
          nj=nj+1
          arsp(nj,nk)=fwtece(m)*recepc(m,n)
 2170   continue
          if (fwtecebz.le.0.0) go to 2175
          nj=nj+1
          arsp(nj,nk)=fwtecebz*recebzpc(n)
 2175   continue
        if (fwtcur.le.0.0) go to 2180
        nj=nj+1
        arsp(nj,nk)=fwtcur*fgowpc(n)
 2180   continue
        if (fwtqa.le.0.0) go to 2190
        nj=nj+1
        arsp(nj,nk)=0.0
 2190 continue
      if (fwtbp.le.0.0) go to 2194
      do 2192 m=2,kffcur
        nj=nj+1
        arsp(nj,nk)=0.0
 2192 continue
 2194 continue
      if (fwtdlc.le.0.0) go to 2196
      nj=nj+1
      arsp(nj,nk)=0.0
 2196 continue
 2202 continue
      do 2203 m=1,nfcoil
        if (fwtfc(m).le.0.0) go to 2203
        nj=nj+1
        arsp(nj,nk)=0.
 2203 continue
!--------------------------------------------------------------------
!--  pressure                                                      --
!--------------------------------------------------------------------
      if (kprfit.le.0.or.kdofit.eq.0) go to 2206
      if (npress.le.0) goto 2206
      do 2204 m=1,npress
        if (fwtpre(m).le.0.0) goto 2204
        nj=nj+1
        arsp(nj,nk)=0.0
        if (n.le.kppcur) arsp(nj,nk)=rprepc(m,n)/sigpre(m)*fwtpre(m)
 2204 continue
 2206 continue
!--------------------------------------------------------------------
!-- P'(1)                                                          --
!--------------------------------------------------------------------
      if (kpressb.eq.2) then
        nj=nj+1
        if (n.le.kppcur) arsp(nj,nk)=1./sigppb/darea
      endif
!---------------------------------------------------------------------
!-- rotational pressure                                             --
!---------------------------------------------------------------------
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          if (fwtprw(m).gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
          if (n.gt.kpcurn) then
            nmw=n-kpcurn
            arsp(nj,nk)=rprwpc(m,nmw)/sigprw(m)*fwtprw(m)
          endif
          endif
        enddo
      endif
!----------------------------------------------------------------------
!--  J(PSIWANT) constraint                                           --
!----------------------------------------------------------------------
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
        ysiwant=sizeroj(i)
        if (n.le.kppcur) then
            call setpp(ysiwant,xpspp)
            xjj=xpspp(n)
            arsp(nj,nk)=rxxx(i)*fwtxxzj*xjj
        elseif (n.le.kpcurn) then
            call setfp(ysiwant,xpsfp)
            xjj=xpsfp(n-kppcur)
            arsp(nj,nk)=fwtxxzj/rxxxf(i)*xjj
        elseif (kvtor.gt.0) then
            xjj=xpspwp(n-kpcurn)
            arsp(nj,nk)=rxxw(i)*fwtxxzj*xjj
        endif
       enddo
      endif
!-------------------------------------------------------------------------
!--  p' and ff' constraints                                             --
!-------------------------------------------------------------------------
      if (kcalpa.gt.0) then
        brspmin=max(ten24,abs(brsp(nfcoil+1)))
        fwtxxa=fwtxx*1000./brspmin
        do 30140 j=1,kcalpa
          nj=nj+1
          if (n.le.kppcur) then
            arsp(nj,nk)=calpa(n,j)*fwtxxa
          else
            arsp(nj,nk)=0.0
          endif
30140   continue
      endif
      if (kcgama.gt.0) then
        brspmin=max(ten24,abs(brsp(nfcoil+1)))
        fwtxxf=fwtxx*1000./brspmin
        do 30150 j=1,kcgama
          nj=nj+1
          if (n.le.kppcur) then
            arsp(nj,nk)=0.0
          elseif (n.le.kpcurn) then
            arsp(nj,nk)=cgama(n-kppcur,j)*fwtxxf
          else
            arsp(nj,nk)=0.0
          endif
30150   continue
      endif
!---------------------------------------------------------------------
!-- rotational constraints                                          --
!---------------------------------------------------------------------
      if (kcomega.gt.0) then
        brspmin=max(ten24,abs(brsp(nfcoil+1)))
        fwtxxo=fwtxx*1000./brspmin
        do j=1,kcomega
          nj=nj+1
          if (n.le.kpcurn) then
            arsp(nj,nk)=0.0
          else
            arsp(nj,nk)=comega(n,j)*fwtxxo
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!-- Boundary constraints due to plasma contributions                         --
!------------------------------------------------------------------------------
      if (nbdry.gt.0) then
        do j=1,nbdry
          if (fwtbdry(j).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=fwtbry(j)*gbdrpc(j,n)
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!--  E coil currents                                                         --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
       do m=1,nesum
        if (fwtec(m).gt.0.0) then
        nj=nj+1
        arsp(nj,nk)=0.
        endif
       enddo
      endif
!------------------------------------------------------------------------------
!--  Reference flux                                                          --
!------------------------------------------------------------------------------
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
      endif
!-------------------------------------------------------------------------
!-- Summation of F-coils currents
!-------------------------------------------------------------------------
      if (fitfcsum) then
         nj = nj + 1
         arsp(nj,nk) = 0.0
      endif
 2210 continue
      need=nfnwcr
!----------------------------------------------------------------------
!-- fit vessel currents                                              --
!----------------------------------------------------------------------
      if (ifitvs.le.0) go to 2310
      if (nfourier.gt.1) then
       need=need+nfourier*2+1
      else
       need=need+nvesel
      endif
      do 2300 nk=nfnwcr+1,need
        mk=nk-nfnwcr
        nj=0
        do 2220 m=1,nsilop
          if (fwtsi(m).le.0.0) go to 2220
          nj=nj+1
          if(nfourier.gt.1)then
          temp=0.
          do i=1,nvesel
          temp=temp+rsilvs(m,i)*vecta(mk,i)
          enddo
           arsp(nj,nk)=fwtsi(m)*temp 
          else
          arsp(nj,nk)=fwtsi(m)*rsilvs(m,mk)
          endif
 2220   continue
        do 2240 m=1,magpri
          if (fwtmp2(m).le.0.0) go to 2240
          nj=nj+1
          if(nfourier.gt.1)then
          temp=0.
          do i=1,nvesel
          temp=temp+rmp2vs(m,i)*vecta(mk,i)
          enddo
           arsp(nj,nk)=fwtmp2(m)*temp        
          else
          arsp(nj,nk)=fwtmp2(m)*rmp2vs(m,mk)
          endif
 2240   continue
        do 2260 m=1,nstark
          if (fwtgam(m).le.0.0) go to 2260
          nj=nj+1
          if(nfourier.gt.1)then
          temp=0.
          do i=1,nvesel
          temp=temp+rgamvs(m,i)*vecta(mk,i)
          enddo
           arsp(nj,nk)=fwtgam(m)*temp
          else
          arsp(nj,nk)=fwtgam(m)*rgamvs(m,mk)
          endif
 2260   continue
        do 2270 m=1,nece
          if (fwtece(m).le.0.0) go to 2270  ! add by qian for ece 
          nj=nj+1
          if(nfourier.gt.1)then
          temp=0.
          do i=1,nvesel
          temp=temp+recevs(m,i)*vecta(mk,i)
          enddo
           arsp(nj,nk)=fwtece(m)*temp
          else
          arsp(nj,nk)=fwtece(m)*recevs(m,mk)
          endif
 2270   continue
          if (fwtecebz.le.0.0) go to 2275
          nj=nj+1
          arsp(nj,nk)=fwtecebz*recebzvs(mk)
 2275   continue
        if (fwtcur.le.0.0) go to 2280
        nj=nj+1
          if(nfourier.gt.1)then
          temp=0.
          do i=1,nvesel
          temp=temp+vecta(mk,i)
          enddo
           arsp(nj,nk)=fwtcur*temp
          else
          arsp(nj,nk)=fwtcur
          endif
 2280   continue
        if (fwtqa.le.0.0) go to 2285
        nj=nj+1
        arsp(nj,nk)=0.0
 2285 continue
      if (fwtbp.le.0.0) go to 2290
      do 2287 m=2,kffcur
        nj=nj+1
        arsp(nj,nk)=0.0
 2287 continue
 2290 continue
      if (fwtdlc.le.0.0) go to 2292
      nj=nj+1
      arsp(nj,nk)=0.0
 2292 continue
 2296 continue
      do 2297 m=1,nfcoil
        if (fwtfc(m).le.0.0) go to 2297
        nj=nj+1
        arsp(nj,nk)=0.
 2297 continue
      if (kprfit.le.0.or.kdofit.eq.0) go to 2299
      if (npress.le.0) goto 2299
      do 2298 m=1,npress
        if (fwtpre(m).le.0.0) goto 2298
        nj=nj+1
        arsp(nj,nk)=0.0
 2298 continue
 2299 continue
      if (kpressb.eq.2) then
        nj=nj+1
        arsp(nj,nk)=0.0
      endif
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          if (fwtprw(m).gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
          endif
        enddo
      endif
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
        arsp(nj,nk)=0.0
       enddo
      endif
      if (kcalpa.gt.0) then
        do 30190 j=1,kcalpa
          nj=nj+1
          arsp(nj,nk)=0.0
30190   continue
      endif
      if (kcgama.gt.0) then
        do 30200 j=1,kcgama
          nj=nj+1
          arsp(nj,nk)=0.0
30200   continue
      endif
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
!------------------------------------------------------------------------------
!-- Boundary constraints                                                     --
!------------------------------------------------------------------------------
      if (nbdry.gt.0) then
        do j=1,nbdry
          if (fwtbdry(j).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=fwtbry(j)*rbdrvs(j,mk)
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!--  E coil currents                                                         --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
       do m=1,nesum
        if (fwtec(m).gt.0.0) then
        nj=nj+1
        arsp(nj,nk)=0.
        endif
       enddo
      endif
!------------------------------------------------------------------------------
!--  Reference flux                                                          --
!------------------------------------------------------------------------------
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
      endif
!------------------------------------------------------------------------------
!--  add by qian  sum of vessel sement currents=IPV1-IP1                    --
!------------------------------------------------------------------------------
      if (isumvesel.gt.0) then
          nj=nj+1
          if(nfourier.gt.0)then
           temp=0.
            do i=1,nvesel
             temp=temp+vecta(mk,i)
            enddo
             arsp(nj,nk)=fwtcur*temp/100.
           else
             arsp(nj,nk)=1.0
          endif
      endif
 2300 continue
 2310 continue
!-----------------------------------------------------------------------
!-- boundary pressure term for kinetic fitting P(1)                   --
!-- P(1) is an additional fitting parameter in kinetic fitting        --
!-----------------------------------------------------------------------
        if (kprfit.le.0.or.kdofit.eq.0) go to 9210
        if (npress.le.0) goto 9210
        need=need+1
        nk=need
        nj=0
        do 9120 m=1,nsilop
          if (fwtsi(m).le.0.0) go to 9120
          nj=nj+1
          arsp(nj,nk)=0.0
 9120   continue
        do 9140 m=1,magpri
          if (fwtmp2(m).le.0.0) go to 9140
          nj=nj+1
          arsp(nj,nk)=0.0
 9140   continue
        do 9160 m=1,nstark
          if (fwtgam(m).le.0.0) go to 9160
          nj=nj+1
          arsp(nj,nk)=0.0
 9160   continue
        do 9170 m=1,nece
          if (fwtece(m).le.0.0) go to 9170
          nj=nj+1
          arsp(nj,nk)=0.0
 9170   continue
          if (fwtecebz.le.0.0) go to 9175
          nj=nj+1
          arsp(nj,nk)=0.0
 9175   continue
        if (fwtcur.le.0.0) go to 9180
        nj=nj+1
        arsp(nj,nk)=0.0
 9180   continue
        if (fwtqa.le.0.0) go to 9190
        nj=nj+1
        arsp(nj,nk)=0.0
 9190 continue
      if (fwtbp.le.0.0) go to 9194
      do 9192 m=2,kffcur
        nj=nj+1
        arsp(nj,nk)=0.0
 9192 continue
 9194 continue
      if (fwtdlc.le.0.0) go to 9196
      nj=nj+1
      arsp(nj,nk)=0.0
 9196 continue
 9202 continue
      do 9203 m=1,nfcoil
        if (fwtfc(m).le.0.0) go to 9203
        nj=nj+1
        arsp(nj,nk)=0.
 9203 continue
      do 9204 m=1,npress
        if (fwtpre(m).le.0.0) goto 9204
        nj=nj+1
        arsp(nj,nk)=1./sigpre(m)*fwtpre(m)
 9204 continue
!-----------------------------------------------------------------------
!--  boundary pressure constraint on P'(1), no coupling to P(1)       --
!-----------------------------------------------------------------------
      if (kpressb.eq.2) then
        nj=nj+1
        arsp(nj,nk)=0.0
      endif
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          if (fwtprw(m).gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
          endif
        enddo
      endif
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
        arsp(nj,nk)=0.0
       enddo
      endif
      if (kcalpa.gt.0) then
        do 30250 j=1,kcalpa
          nj=nj+1
          arsp(nj,nk)=0.0
30250   continue
      endif
      if (kcgama.gt.0) then
        do 30260 j=1,kcgama
          nj=nj+1
          arsp(nj,nk)=0.0
30260   continue
      endif
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
      if (nbdry.gt.0) then
        do j=1,nbdry
          if (fwtbdry(j).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.0
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!--  E coil currents                                                         --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
       do m=1,nesum
        if (fwtec(m).gt.0.0) then
        nj=nj+1
        arsp(nj,nk)=0.
        endif
       enddo
      endif
!------------------------------------------------------------------------------
!--  Reference flux                                                          --
!------------------------------------------------------------------------------
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
      endif
!-----------------------------------------------------------------------------
!-- Summation of Fcoils current
!-----------------------------------------------------------------------------
      if (fitfcsum) then
         nj = nj + 1
         arsp(nj,nk) = 0.0
      endif
 9210 continue
!-----------------------------------------------------------------------
!-- DELZ rigid vertical shift   96/01                                 --
!-----------------------------------------------------------------------
      if (fitdelz.and.nniter.ge.ndelzon) then
        need=need+1
        nsavdz=need
        nk=need
        nj=0
        do 39120 m=1,nsilop
          if (fwtsi(m).le.0.0) go to 39120
          nj=nj+1
          arsp(nj,nk)=fwtsi(m)*gsildz(m)*scadelz
39120   continue
        do 39140 m=1,magpri
          if (fwtmp2(m).le.0.0) go to 39140
          nj=nj+1
          arsp(nj,nk)=fwtmp2(m)*gmp2dz(m)*scadelz
39140   continue
        do 39160 m=1,nstark
          if (fwtgam(m).le.0.0) go to 39160
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgamdz(m)*scadelz
39160   continue
        do 39170 m=1,nece
          if (fwtece(m).le.0.0) go to 39170
          nj=nj+1
          arsp(nj,nk)=fwtece(m)*recedz(m)*scaledz
39170   continue
          if (fwtecebz.le.0.0) go to 39175
          nj=nj+1
          arsp(nj,nk)=fwtecebz*recebzdz*scaledz
39175   continue
        if (fwtcur.le.0.0) go to 39180
          nj=nj+1
          arsp(nj,nk)=fwtcur*fgowdz*scadelz
39180   continue
        if (fwtqa.le.0.0) go to 39190
          nj=nj+1
          arsp(nj,nk)=0.0
39190   continue
        if (fwtbp.le.0.0) go to 39194
          do 39192 m=2,kffcur
            nj=nj+1
            arsp(nj,nk)=0.0
39192     continue
39194   continue
        if (fwtdlc.le.0.0) go to 39196
          nj=nj+1
          arsp(nj,nk)=0.0
39196   continue
        do 39203 m=1,nfcoil
          if (fwtfc(m).le.0.0) go to 39203
          nj=nj+1
          arsp(nj,nk)=0.
39203   continue
        if (kprfit.le.0.or.kdofit.eq.0) go to 39206
        if (npress.le.0) goto 39206
        do 39204 m=1,npress
          if (fwtpre(m).le.0.0) goto 39204
          nj=nj+1
          arsp(nj,nk)=rpredz(m)/sigpre(m)*fwtpre(m)*scadelz
39204   continue
39206   continue
        if (kpressb.eq.2) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        if (kprfit.ge.3.and.npresw.gt.0) then
          do m=1,npresw
            if (fwtprw(m).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=rprwdz(m)/sigprw(m)*fwtprw(m)*scadelz
            endif
          enddo
        endif
        if (kzeroj.gt.0) then
          do i=1,kzeroj
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (kcalpa.gt.0) then
          do 39250 j=1,kcalpa
            nj=nj+1
            arsp(nj,nk)=0.0
39250     continue
        endif
        if (kcgama.gt.0) then
          do 39260 j=1,kcgama
            nj=nj+1
            arsp(nj,nk)=0.0
39260     continue
        endif
        if (kcomega.gt.0) then
          do j=1,kcomega
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (nbdry.gt.0) then
          do j=1,nbdry
            if (fwtbdry(j).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=fwtbry(j)*gbdrdz(j)*scadelz
            endif
          enddo
        endif
        if (iecurr.eq.2) then
         do m=1,nesum
          if (fwtec(m).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.
          endif
         enddo
        endif
        if (fitsiref) then
          if (fwtref.gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.0
          endif
        endif
!-----------------------------------------------------------------------------
!-- Summation of Fcoils current
!-----------------------------------------------------------------------------
        if (fitfcsum) then
           nj = nj + 1
           arsp(nj,nk) = 0.0
        endif
      endif
!-----------------------------------------------------------------------
!--  set up response matrix for advanced divertor coil                --
!--  lacking boundary constraints                                     --
!-----------------------------------------------------------------------
      if (iacoil.le.0) go to 9410
      nkb=need+1
      need=need+nacoil
      do 9400 nk=nkb,need
        mk=nk-nkb+1
        nj=0
        do 9320 m=1,nsilop
          if (fwtsi(m).le.0.0) go to 9320
          nj=nj+1
          arsp(nj,nk)=fwtsi(m)*rsilac(m,mk)
 9320   continue
        do 9340 m=1,magpri
          if (fwtmp2(m).le.0.0) go to 9340
          nj=nj+1
          arsp(nj,nk)=fwtmp2(m)*rmp2ac(m,mk)
 9340   continue
        do 9360 m=1,nstark
          if (fwtgam(m).le.0.0) go to 9360
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgamac(m,mk)
 9360   continue
        do 9370 m=1,nece
          if (fwtece(m).le.0.0) go to 9370
          nj=nj+1
          arsp(nj,nk)=fwtece(m)*receac(m,mk)
 9370   continue
          if (fwtecebz.le.0.0) go to 9375
          nj=nj+1
          arsp(nj,nk)=fwtecebz*recebzac(mk)
 9375   continue
        if (fwtcur.le.0.0) go to 9380
        nj=nj+1
        arsp(nj,nk)=0.0
 9380   continue
        if (fwtqa.le.0.0) go to 9385
        nj=nj+1
        arsp(nj,nk)=0.0
 9385 continue
      if (fwtbp.le.0.0) go to 9390
      do 9387 m=2,kffcur
        nj=nj+1
        arsp(nj,nk)=0.0
 9387 continue
 9390 continue
      if (fwtdlc.le.0.0) go to 9392
      nj=nj+1
      arsp(nj,nk)=0.0
 9392 continue
      do 9397 m=1,nfcoil
        if (fwtfc(m).le.0.0) go to 9397
        nj=nj+1
        arsp(nj,nk)=0.
 9397 continue
      if (kprfit.le.0.or.kdofit.eq.0) go to 9399
      if (npress.le.0) goto 9399
      do 9398 m=1,npress
        if (fwtpre(m).le.0.0) goto 9398
        nj=nj+1
        arsp(nj,nk)=0.0
 9398 continue
 9399 continue
      if (kpressb.eq.2) then
        nj=nj+1
        arsp(nj,nk)=0.0
      endif
!---------------------------------------------------------------------------
!-- rotational pressure                                                   --
!---------------------------------------------------------------------------
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          if (fwtprw(m).gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
          endif
        enddo
      endif
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
        arsp(nj,nk)=0.0
       enddo
      endif
      if (kcalpa.gt.0) then
        do 39390 j=1,kcalpa
          nj=nj+1
          arsp(nj,nk)=0.0
39390   continue
      endif
      if (kcgama.gt.0) then
        do 39400 j=1,kcgama
          nj=nj+1
          arsp(nj,nk)=0.0
39400   continue
      endif
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
!------------------------------------------------------------------------------
!--  E coil currents                                                         --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
       do m=1,nesum
        if (fwtec(m).gt.0.0) then
        nj=nj+1
        arsp(nj,nk)=0.
        endif
       enddo
      endif
!------------------------------------------------------------------------------
!--  Reference flux                                                          --
!------------------------------------------------------------------------------
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
      endif
!-----------------------------------------------------------------------------
!-- Summation of Fcoils current
!-----------------------------------------------------------------------------
      if (fitfcsum) then
         nj = nj + 1
         arsp(nj,nk) = 0.0
      endif
 9400 continue
 9410 continue
!-----------------------------------------------------------------------
!--  set up response matrix for E coils                               --
!-----------------------------------------------------------------------
      if (iecurr.ne.2) go to 48410
      nkb=need+1
      need=need+nesum
      do 48400 nk=nkb,need
        mk=nk-nkb+1
        nj=0
        do 48320 m=1,nsilop
          if (fwtsi(m).le.0.0) go to 48320
          nj=nj+1
          arsp(nj,nk)=fwtsi(m)*rsilec(m,mk)
48320   continue
        do 48340 m=1,magpri
          if (fwtmp2(m).le.0.0) go to 48340
          nj=nj+1
          arsp(nj,nk)=fwtmp2(m)*rmp2ec(m,mk)
48340   continue
        do 48360 m=1,nstark
          if (fwtgam(m).le.0.0) go to 48360
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgamec(m,mk)
48360   continue
        do 48370 m=1,nece
          if (fwtece(m).le.0.0) go to 48370
          nj=nj+1
          arsp(nj,nk)=fwtece(m)*receec(m,mk)
48370   continue
          if (fwtecebz.le.0.0) go to 48375
          nj=nj+1
          arsp(nj,nk)=fwtecebz*recebzec(mk)
48375   continue
        if (fwtcur.le.0.0) go to 48380
        nj=nj+1
        arsp(nj,nk)=0.0
48380   continue
        if (fwtqa.le.0.0) go to 48385
        nj=nj+1
        arsp(nj,nk)=0.0
48385 continue
      if (fwtbp.le.0.0) go to 48390
      do 48387 m=2,kffcur
        nj=nj+1
        arsp(nj,nk)=0.0
48387 continue
48390 continue
      if (fwtdlc.le.0.0) go to 48392
      nj=nj+1
      arsp(nj,nk)=0.0
48392 continue
      do 48397 m=1,nfcoil
        if (fwtfc(m).le.0.0) go to 48397
        nj=nj+1
        arsp(nj,nk)=0.
48397 continue
      if (kprfit.le.0.or.kdofit.eq.0) go to 48399
      if (npress.le.0) goto 48399
      do 48398 m=1,npress
        if (fwtpre(m).le.0.0) goto 48398
        nj=nj+1
        arsp(nj,nk)=0.0
48398 continue
48399 continue
      if (kpressb.eq.2) then
        nj=nj+1
        arsp(nj,nk)=0.0
      endif
!---------------------------------------------------------------------------
!-- rotational pressure                                                   --
!---------------------------------------------------------------------------
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          if (fwtprw(m).gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
          endif
        enddo
      endif
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
        arsp(nj,nk)=0.0
       enddo
      endif
      if (kcalpa.gt.0) then
        do 48490 j=1,kcalpa
          nj=nj+1
          arsp(nj,nk)=0.0
48490   continue
      endif
      if (kcgama.gt.0) then
        do 48495 j=1,kcgama
          nj=nj+1
          arsp(nj,nk)=0.0
48495   continue
      endif
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
!------------------------------------------------------------------------------
!-- Boundary constraints                                                     --
!------------------------------------------------------------------------------
      if (nbdry.gt.0) then
        do j=1,nbdry
          if (fwtbdry(j).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=fwtbry(j)*rbdrec(j,mk)
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!--  E coil currents                                                         --
!------------------------------------------------------------------------------
      do 48497 m=1,nesum
        if (fwtec(m).le.0.0) go to 48497
        nj=nj+1
        arsp(nj,nk)=0.
        if (m.eq.mk) arsp(nj,nk)=fwtec(m)
48497 continue
!------------------------------------------------------------------------------
!--  Reference flux                                                          --
!------------------------------------------------------------------------------
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
      endif
!-----------------------------------------------------------------------------
!-- Summation of Fcoils current
!-----------------------------------------------------------------------------
      if (fitfcsum) then
         nj = nj + 1
         arsp(nj,nk) = 0.0
      endif
48400 continue
48410 continue
!------------------------------------------------------------------------------
!--  fitting relative flux, set up response for fitted reference flux        --
!------------------------------------------------------------------------------
        if (.not.fitsiref) goto 52501
        need=need+1
        nj=0
        nk=need
        do 52020 m=1,nsilop
          if (fwtsi(m).le.0.0) go to 52020
          nj=nj+1
          arsp(nj,nk)=-fwtsi(m)*scalesir
52020   continue
        do 52040 m=1,magpri
          if (fwtmp2(m).le.0.0) go to 52040
          nj=nj+1
          arsp(nj,nk)=0.0
52040   continue
        do 52060 m=1,nstark
          if (fwtgam(m).le.0.0) go to 52060
          nj=nj+1
          arsp(nj,nk)=0.0
52060   continue
        do 52070 m=1,nece
          if (fwtece(m).le.0.0) go to 52070
          nj=nj+1
          arsp(nj,nk)=0.0
52070   continue
          if (fwtecebz.le.0.0) go to 52075
          nj=nj+1
          arsp(nj,nk)=0.0
52075   continue
        if (fwtcur.le.0.0) go to 52080
        nj=nj+1
        arsp(nj,nk)=0.0
52080   continue
        if (fwtqa.le.0.0) go to 52085
        nj=nj+1
        arsp(nj,nk)=0.0
52085 continue
      if (fwtbp.le.0.0) go to 52090
      do 52087 m=2,kffcur
        nj=nj+1
        arsp(nj,nk)=0.0
52087 continue
52090 continue
      if (fwtdlc.le.0.0) go to 52092
      nj=nj+1
      arsp(nj,nk)=0.0
52092 continue
52096 continue
      do 52097 m=1,nfcoil
        if (fwtfc(m).le.0.0) go to 52097
        nj=nj+1
        arsp(nj,nk)=0.
52097 continue
!--------------------------------------------------------------------------
!--  pressure                                                            --
!--------------------------------------------------------------------------
      if (kprfit.le.0.or.kdofit.eq.0) go to 52099
      if (npress.le.0) goto 52099
      do 52098 m=1,npress
        if (fwtpre(m).le.0.0) goto 52098
        nj=nj+1
        arsp(nj,nk)=0.0
52098 continue
52099 continue
      if (kpressb.eq.2) then
        nj=nj+1
        arsp(nj,nk)=0.0
      endif
!---------------------------------------------------------------------------
!-- rotational pressure                                                   --
!---------------------------------------------------------------------------
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          if (fwtprw(m).gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
          endif
        enddo
      endif
!---------------------------------------------------------------------------
!--  J(PSIWANT)                                                           --
!---------------------------------------------------------------------------
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
        arsp(nj,nk)=0.0
       enddo
      endif
!----------------------------------------------------------------------------
!-- P', FF', and rotational                                                --
!----------------------------------------------------------------------------
      if (kcalpa.gt.0) then
        do 53050 j=1,kcalpa
          nj=nj+1
          arsp(nj,nk)=0.0
53050   continue
      endif
      if (kcgama.gt.0) then
        do 53060 j=1,kcgama
          nj=nj+1
          arsp(nj,nk)=0.0
53060   continue
      endif
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
!------------------------------------------------------------------------------
!-- Boundary constraints                                                     --
!------------------------------------------------------------------------------
      if (nbdry.gt.0) then
        do j=1,nbdry
          if (fwtbdry(j).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.0
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!--  E coil currents                                                         --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
       do m=1,nesum
        if (fwtec(m).gt.0.0) then
        nj=nj+1
        arsp(nj,nk)=0.
        endif
       enddo
      endif
!------------------------------------------------------------------------------
!--  Reference flux                                                          --
!------------------------------------------------------------------------------
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=fwtref*scalesir
        endif
      endif
!-----------------------------------------------------------------------------
!-- Summation of Fcoils current
!-----------------------------------------------------------------------------
      if (fitfcsum) then
         nj = nj + 1
         arsp(nj,nk) = 0.0
      endif
52501 continue
!-----------------------------------------------------------------------
!-- set up response matrix for ER                                     --
!-----------------------------------------------------------------------
      if (keecur.le.0.or.kdomse.gt.0) goto 72111
        needs=need
        need=need+keecur
      do 72100 nk=needs+1,need
        nkk=nk-needs
        nj=0
        do 72020 m=1,nsilop
          if (fwtsi(m).le.0.0) go to 72020
          nj=nj+1
          arsp(nj,nk)=0.0
72020   continue
        do 72040 m=1,magpri
          if (fwtmp2(m).le.0.0) go to 72040
          nj=nj+1
          arsp(nj,nk)=0.0
72040   continue
        do 72060 m=1,nstark
          if (fwtgam(m).le.0.0) go to 72060
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgamer(m,nkk)
72060   continue
        do 72070 m=1,nece
          if (fwtece(m).le.0.0) go to 72070
          nj=nj+1
          arsp(nj,nk)=0.0
72070   continue
          if (fwtecebz.le.0.0) go to 72075
          nj=nj+1
          arsp(nj,nk)=0.0
72075   continue
        if (fwtcur.le.0.0) go to 72080
        nj=nj+1
        arsp(nj,nk)=0.0
72080   continue
        if (fwtqa.le.0.0) go to 72085
        nj=nj+1
        arsp(nj,nk)=0.0
72085 continue
      if (fwtbp.le.0.0) go to 72090
      do 72087 m=2,kffcur
        nj=nj+1
        arsp(nj,nk)=0.0
72087 continue
72090 continue
      if (fwtdlc.le.0.0) go to 72092
      nj=nj+1
      arsp(nj,nk)=0.0
72092 continue
      do 72097 m=1,nfcoil
        if (fwtfc(m).le.0.0) go to 72097
        nj=nj+1
        arsp(nj,nk)=0.
72097 continue
      if (kprfit.le.0.or.kdofit.eq.0) go to 72099
      if (npress.le.0) goto 72099
      do 72098 m=1,npress
        if (fwtpre(m).le.0.0) goto 72098
        nj=nj+1
        arsp(nj,nk)=0.0
72098 continue
72099 continue
      if (kpressb.eq.2) then
        nj=nj+1
        arsp(nj,nk)=0.0
      endif
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          if (fwtprw(m).gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
          endif
        enddo
      endif
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
        arsp(nj,nk)=0.0
       enddo
      endif
      if (kcalpa.gt.0) then
        do j=1,kcalpa
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
      if (kcgama.gt.0) then
        do j=1,kcgama
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
      if (nbdry.gt.0) then
        do j=1,nbdry
          if (fwtbdry(j).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.0
          endif
        enddo
      endif
      if (iecurr.eq.2) then
       do m=1,nesum
        if (fwtec(m).gt.0.0) then
        nj=nj+1
        arsp(nj,nk)=0.
        endif
       enddo
      endif
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
      endif
!-----------------------------------------------------------------------------
!-- Summation of Fcoils current
!-----------------------------------------------------------------------------
      if (fitfcsum) then
         nj = nj + 1
         arsp(nj,nk) = 0.0
      endif
72100 continue
72111 continue
!-----------------------------------------------------------------------
!-- Response for Pedge hyperbolic tangent                             --
!-----------------------------------------------------------------------
      if (kedgep.eq.0) go to 73311
        need=need+1
        nk=need
        npedge=need
        nj=0
        do 73120 m=1,nsilop
          if (fwtsi(m).le.0.0) go to 73120
          nj=nj+1
          arsp(nj,nk)=fwtsi(m)*rsilpe(m)
73120   continue
        do 73140 m=1,magpri
          if (fwtmp2(m).le.0.0) go to 73140
          nj=nj+1
          arsp(nj,nk)=fwtmp2(m)*rmp2pe(m)
73140   continue
        do 73160 m=1,nstark
          if (fwtgam(m).le.0.0) go to 73160
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgampe(m)
73160   continue
        if (fwtcur.le.0.0) go to 73180
        nj=nj+1
        arsp(nj,nk)=fwtcur*fgowpe
73180   continue
        if (fwtqa.le.0.0) go to 73190
        nj=nj+1
        arsp(nj,nk)=0.0
73190   continue
        if (fwtbp.le.0.0) go to 73194
        do 73192 m=2,kffcur
          nj=nj+1
          arsp(nj,nk)=0.0
73192   continue
73194   continue
        if (fwtdlc.le.0.0) go to 73196
        nj=nj+1
        arsp(nj,nk)=0.0
73196   continue
        do 73203 m=1,nfcoil
          if (fwtfc(m).le.0.0) go to 73203
          nj=nj+1
          arsp(nj,nk)=0.
73203   continue
!--------------------------------------------------------------------
!--  pressure                                                      --
!--------------------------------------------------------------------
        if (kprfit.le.0.or.kdofit.eq.0) go to 73206
        if (npress.le.0) goto 73206
        do 73204 m=1,npress
          if (fwtpre(m).le.0.0) goto 73204
          nj=nj+1
          arsp(nj,nk)=rprepe(m)/sigpre(m)*fwtpre(m)
73204 continue
73206 continue
!--------------------------------------------------------------------
!-- P'(1)                                                          --
!--------------------------------------------------------------------
      if (kpressb.eq.2) then
        nj=nj+1
        arsp(nj,nk)=1./sigppb/darea/cosh(s1edge)**2/pe_width/sidif
      endif
!---------------------------------------------------------------------
!-- rotational pressure                                             --
!---------------------------------------------------------------------
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          if (fwtprw(m).gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
          endif
        enddo
      endif
!----------------------------------------------------------------------
!--  J(PSIWANT) constraint                                           --
!----------------------------------------------------------------------
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
        siedge=(sizeroj(i)-pe_psin)/pe_width
        arsp(nj,nk)=rxxx(i)*fwtxxzj/cosh(siedge)**2/pe_width/sidif
       enddo
      endif
!-------------------------------------------------------------------------
!--  p' and ff' constraints                                             --
!-------------------------------------------------------------------------
      if (kcalpa.gt.0) then
        do 73240 j=1,kcalpa
          nj=nj+1
          arsp(nj,nk)=calpa(kppcur+1,j)*fwtxxa
73240   continue
      endif
      if (kcgama.gt.0) then
        do 73250 j=1,kcgama
          nj=nj+1
          arsp(nj,nk)=0.0
73250   continue
      endif
!---------------------------------------------------------------------
!-- rotational constraints                                          --
!---------------------------------------------------------------------
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
!------------------------------------------------------------------------------
!-- Boundary constraints due to plasma contributions                         --
!------------------------------------------------------------------------------
      if (nbdry.gt.0) then
        do j=1,nbdry
          if (fwtbdry(j).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=fwtbry(j)*gbdrpe(j)
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!--  E coil currents                                                         --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
       do m=1,nesum
        if (fwtec(m).gt.0.0) then
        nj=nj+1
        arsp(nj,nk)=0.
        endif
       enddo
      endif
!------------------------------------------------------------------------------
!--  Reference flux                                                          --
!------------------------------------------------------------------------------
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
      endif
!-----------------------------------------------------------------------------
!-- Summation of Fcoils current
!-----------------------------------------------------------------------------
      if (fitfcsum) then
         nj = nj + 1
         arsp(nj,nk) = 0.0
      endif
73311 continue
!-----------------------------------------------------------------------
!-- Response for f2edge hyperbolic tangent                            --
!-----------------------------------------------------------------------
      if (kedgef.eq.0) go to 74311
        need=need+1
        nk=need
        nfedge=need
        nj=0
        do 74120 m=1,nsilop
          if (fwtsi(m).le.0.0) go to 74120
          nj=nj+1
          arsp(nj,nk)=fwtsi(m)*rsilfe(m)
74120   continue
        do 74140 m=1,magpri
          if (fwtmp2(m).le.0.0) go to 74140
          nj=nj+1
          arsp(nj,nk)=fwtmp2(m)*rmp2fe(m)
74140   continue
        do 74160 m=1,nstark
          if (fwtgam(m).le.0.0) go to 74160
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgamfe(m)
74160   continue
        if (fwtcur.le.0.0) go to 74180
        nj=nj+1
        arsp(nj,nk)=fwtcur*fgowfe
74180   continue
        if (fwtqa.le.0.0) go to 74190
        nj=nj+1
        arsp(nj,nk)=0.0
74190   continue
        if (fwtbp.le.0.0) go to 74194
        do 74192 m=2,kffcur
          nj=nj+1
          arsp(nj,nk)=0.0
74192   continue
74194   continue
        if (fwtdlc.le.0.0) go to 74196
        nj=nj+1
        arsp(nj,nk)=0.0
74196   continue
        do 74203 m=1,nfcoil
          if (fwtfc(m).le.0.0) go to 74203
          nj=nj+1
          arsp(nj,nk)=0.
74203   continue
!--------------------------------------------------------------------
!--  pressure                                                      --
!--------------------------------------------------------------------
        if (kprfit.le.0.or.kdofit.eq.0) go to 74206
        if (npress.le.0) goto 74206
        do 74204 m=1,npress
          if (fwtpre(m).le.0.0) goto 74204
          nj=nj+1
          arsp(nj,nk)=0.0
74204 continue
74206 continue
!--------------------------------------------------------------------
!-- P'(1)                                                          --
!--------------------------------------------------------------------
      if (kpressb.eq.2) then
        nj=nj+1
        arsp(nj,nk)=0.0
      endif
!---------------------------------------------------------------------
!-- rotational pressure                                             --
!---------------------------------------------------------------------
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          if (fwtprw(m).gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
          endif
        enddo
      endif
!----------------------------------------------------------------------
!--  J(PSIWANT) constraint                                           --
!----------------------------------------------------------------------
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
        siedge=(sizeroj(i)-fe_psin)/fe_width
        arsp(nj,nk)=fwtxxzj/rxxxf(i)/cosh(siedge)**2/fe_width/sidif
       enddo
      endif
!-------------------------------------------------------------------------
!--  p' and ff' constraints                                             --
!-------------------------------------------------------------------------
      if (kcalpa.gt.0) then
        do 74240 j=1,kcalpa
          nj=nj+1
          arsp(nj,nk)=0.0
74240   continue
      endif
      if (kcgama.gt.0) then
        do 74250 j=1,kcgama
          nj=nj+1
          arsp(nj,nk)=cgama(kffcur+1,j)*fwtxxf
74250   continue
      endif
!---------------------------------------------------------------------
!-- rotational constraints                                          --
!---------------------------------------------------------------------
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
!------------------------------------------------------------------------------
!-- Boundary constraints due to plasma contributions                         --
!------------------------------------------------------------------------------
      if (nbdry.gt.0) then
        do j=1,nbdry
          if (fwtbdry(j).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=fwtbry(j)*gbdrfe(j)
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!--  E coil currents                                                         --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
       do m=1,nesum
        if (fwtec(m).gt.0.0) then
        nj=nj+1
        arsp(nj,nk)=0.
        endif
       enddo
      endif
!------------------------------------------------------------------------------
!--  Reference flux                                                          --
!------------------------------------------------------------------------------
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
      endif
!-----------------------------------------------------------------------------
!-- Summation of Fcoils current
!-----------------------------------------------------------------------------
      if (fitfcsum) then
         nj = nj + 1
         arsp(nj,nk) = 0.0
      endif
74311 continue
!-----------------------------------------------------------------------
!-- stabilization constraint for dz, setup here for all fitting paras --
!-----------------------------------------------------------------------
      if (fitdelz.and.nniter.ge.ndelzon) then
        if (stabdz.gt.0.0) then
        nj=nj+1
        do i=1,need
          arsp(nj,i )=0.0
        enddo
        arsp(nj,nsavdz)=stabdz
        endif
      endif
!-----------------------------------------------------------------------
!-- user defined constraint equations
!-----------------------------------------------------------------------
      if (kccoils.gt.0) then
         do m=1,kccoils
           nj=nj+1
           do i=1,nfcoil
             arsp(nj,i)=ccoils(i,m)
           enddo
         enddo
      endif
!-----------------------------------------------------------------------
!-- set up the corresponding data                                     --
!-----------------------------------------------------------------------
      nj=0
      do 2332 m=1,nsilop
        if (fwtsi(m).le.0.0) go to 2332
        nj=nj+1
        brsp(nj)=fwtsi(m)*silopt(jtime,m)
        if (iecurr.ne.1) go to 2320
        ework=0.0
        do 2315 ne=1,nesum
          ework=ework+rsilec(m,ne)*ecurrt(ne)
 2315   continue
        brsp(nj)=brsp(nj)-fwtsi(m)*ework
 2320   if (ivesel.le.0.or.ifitvs.gt.0) go to 2332
        ework=0.0
        do 2330 ne=1,nvesel
          ework=ework+rsilvs(m,ne)*vcurrt(ne)
 2330   continue
        brsp(nj)=brsp(nj)-fwtsi(m)*ework
 2332 continue
      do 2352 m=1,magpri
        if (fwtmp2(m).le.0.0) go to 2352
        nj=nj+1
        brsp(nj)=fwtmp2(m)*expmpi(jtime,m)
        if (iecurr.ne.1) go to 2340
        ework=0.0
        do 2335 ne=1,nesum
          ework=ework+rmp2ec(m,ne)*ecurrt(ne)
 2335   continue
        brsp(nj)=brsp(nj)-fwtmp2(m)*ework
 2340 continue
        if (ivesel.le.0.or.ifitvs.gt.0) go to 2352
        ework=0.0
        do 2345 ne=1,nvesel
          ework=ework+rmp2vs(m,ne)*vcurrt(ne)
 2345   continue
        brsp(nj)=brsp(nj)-fwtmp2(m)*ework
 2352 continue
      do 2360 m=1,nstark
        if (fwtgam(m).le.0.0) go to 2360
        nj=nj+1
        brsp(nj)=fwtgam(m)*rhsgam(jtime,m)
        if (iecurr.ne.1) go to 2356
        ework=0.0
        do 2354 ne=1,nesum
          ework=ework+rgamec(m,ne)*ecurrt(ne)
 2354   continue
        brsp(nj)=brsp(nj)-fwtgam(m)*ework
 2356 continue
        if (ivesel.le.0.or.ifitvs.gt.0) go to 2360
        ework=0.0
        do 2358 ne=1,nvesel
          ework=ework+rgamvs(m,ne)*vcurrt(ne)
 2358   continue
        brsp(nj)=brsp(nj)-fwtgam(m)*ework
 2360 continue
      do 2370 m=1,nece
        if (fwtece(m).le.0.0) go to 2370
        nj=nj+1
        brsp(nj)=fwtece(m)*brspece(jtime,m)
        if (iecurr.ne.1) go to 2366
        ework=0.0
        do 2364 ne=1,nesum
          ework=ework+receec(m,ne)*ecurrt(ne)
 2364   continue
        brsp(nj)=brsp(nj)-fwtece(m)*ework
 2366 continue
        if (ivesel.le.0.or.ifitvs.gt.0) go to 2370
        ework=0.0
        do 2368 ne=1,nvesel
          ework=ework+recevs(m,ne)*vcurrt(ne)
 2368   continue
        brsp(nj)=brsp(nj)-fwtece(m)*ework
 2370 continue
        if (fwtecebz.le.0.0) go to 2379
        nj=nj+1
        brsp(nj)=fwtecebz*brspecebz(jtime)
        if (iecurr.ne.1) go to 2376
        ework=0.0
        do 2374 ne=1,nesum
          ework=ework+recebzec(ne)*ecurrt(ne)
 2374   continue
        brsp(nj)=brsp(nj)-fwtecebz*ework
 2376 continue
        if (ivesel.le.0.or.ifitvs.gt.0) go to 2379
        ework=0.0
        do 2378 ne=1,nvesel
          ework=ework+recebzvs(ne)*vcurrt(ne)
 2378   continue
        brsp(nj)=brsp(nj)-fwtecebz*ework
 2379 continue
      if (fwtcur.le.0.0) go to 2380
      nj=nj+1
      brsp(nj)=fwtcur*pasmat(jtime)
 2380 continue
      if (fwtqa.le.0.0) go to 2400
      nj=nj+1
      do i=1,kppcur
        arsp(nj,nfcoil+i)=fwtqa*rqajx(i)
      enddo
      do i=1,kffcur
        arsp(nj,nbase+i)=fwtqa*rqafx(i)
      enddo
      if (kedgep.gt.0) then
        arsp(nj,npedge)=fwtqa*rqape
      endif
      if (kedgef.gt.0) then
        arsp(nj,nfedge)=fwtqa*rqafe
      endif
      if (kvtor.gt.0) then
        do i=1,kwwcur
          arsp(nj,nfcoil+kpcurn+i)=fwtqa*rqawx(i)
        enddo
      endif
      brsp(nj)=fwtqa/qvfit*pasmat(jtime)/abs(pasmat(jtime))
 2400 continue
      if (fwtbp.le.0.0) go to 2450
      fwtbpp=fwtbp/abs(brsold(nfcoil+1)*brsold(nbase+2))
      do 2420 i=2,kffcur
        nj=nj+1
        arsp(nj,nfcoil+1)=fwtbpp*brsold(nbase+i)
        arsp(nj,nfcoil+i)=-fwtbpp*brsold(nbase+1)
        arsp(nj,nbase+1)=-fwtbpp*brsold(nfcoil+i)
        arsp(nj,nbase+i)=fwtbpp*brsold(nfcoil+1)
        brsp(nj)=fwtbpp*(brsold(nfcoil+1)*brsold(nbase+i)- &
                        brsold(nfcoil+i)*brsold(nbase+1))
 2420 continue
 2450 continue
      if (fwtdlc.le.0.0) go to 2470
      nj=nj+1
      do 2460 i=1,kffcur
        arsp(nj,nbase+i)=fwtdlc*rspdlc(i)
 2460 continue
      if (kedgef.gt.0) then
        arsp(nj,nfedge)=fwtdlc*rdlcfe
      endif
      brsp(nj)=-fwtdlc*diamag(jtime)
 2470 continue
 2475 continue
      do 2476 i=1,nfcoil
        if (fwtfc(i).le.0.0) go to 2476
        nj=nj+1
        brsp(nj)=fccurt(jtime,i)*fwtfc(i)
 2476 continue
!------------------------------------------------------------------
!--  pressure                                                    --
!------------------------------------------------------------------
      if (kprfit.le.0.or.kdofit.eq.0) go to 2480
      if (npress.le.0) goto 2480
      do 2477 i=1,npress
        if (fwtpre(i).le.0.0) goto 2477
        nj=nj+1
        brsp(nj)=pressr(i)/sigpre(i)*fwtpre(i)
 2477 continue
 2480 continue
      if (kpressb.eq.2) then
        nj=nj+1
        brsp(nj)=prespb/sigppb
      endif
!-------------------------------------------------------------------
!--  rotational pressure                                          --
!-------------------------------------------------------------------
      if (kprfit.ge.3.and.npresw.gt.0) then
        do i=1,npresw
          if (fwtprw(i).gt.0.0) then
          nj=nj+1
          brsp(nj)=(presw(i)-preswb)/sigprw(i)*fwtprw(i)
          endif
        enddo
      endif
!--------------------------------------------------------------------
!-- J(PSIWANT)                                                     --
!--------------------------------------------------------------------
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
        brsp(nj)=fwtxxzj*vzeroj(i)*darea*pasmat(jtime)/carea
       enddo
      endif
!--------------------------------------------------------------------
!--  P' and FF'                                                    --
!--------------------------------------------------------------------
      if (kcalpa.gt.0) then
        do 30430 j=1,kcalpa
          nj=nj+1
          brsp(nj)=fwtxxa*xalpa(j)*darea
30430   continue
      endif
      if (kcgama.gt.0) then
        do 30440 j=1,kcgama
          nj=nj+1
          brsp(nj)=fwtxxf*xgama(j)*darea
30440   continue
      endif
!---------------------------------------------------------------------
!--  rotational                                                     --
!---------------------------------------------------------------------
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          brsp(nj)=xomega(j)*darea*fwtxxo
        enddo
      endif
!------------------------------------------------------------------------------
!-- Boundary constraints                                                     --
!------------------------------------------------------------------------------
      if (nbdry.gt.0) then
        do j=1,nbdry
          if (fwtbdry(j).gt.0.0) then
            nj=nj+1
            brsp(nj)=0.0
            if (iecurr.eq.1) then
              do k=1,nesum
                brsp(nj)=brsp(nj)+rbdrec(j,k)*ecurrt(k)
              enddo
            endif
            brsp(nj)=fwtbry(j)*(psibry-brsp(nj))
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!--  E coil currents                                                         --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
      do i=1,nesum
        if (fwtec(i).gt.0.0) then
         nj=nj+1
         brsp(nj)=ecurrt(i)*fwtec(i)
        endif
      enddo
      endif
!------------------------------------------------------------------------------
!--  fitting relative flux                                                   --
!------------------------------------------------------------------------------
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          brsp(nj)=fwtref*psiref(jtime)
        endif
      endif
!-----------------------------------------------------------------------
!-- Summation of F-coils current
!----------------------------------------------------------------------
      if (fitfcsum) then
         nj = nj + 1
         brsp(nj) = 0.0
      endif
!-----------------------------------------------------------------------
!-- stabilization constraint for dz                                   --
!-----------------------------------------------------------------------
      if (fitdelz.and.nniter.ge.ndelzon) then
        if (stabdz.gt.0.0) then
        nj=nj+1
        brsp(nj)=0.0
        endif
      endif
!-----------------------------------------------------------------------
!-- user defined constraint equations
!-----------------------------------------------------------------------
      if (kccoils.gt.0) then
         do j=1,kccoils
           nj=nj+1
           brsp(nj)=xcoils(j)
         enddo
      endif
!
      nnn=1
      ncrsp = 0
      nfffff = nfcoil
      call ppcnst(ncrsp,crsp,z,nfffff)
      call ffcnst(ncrsp,crsp,z,nfffff)
      call wwcnst(ncrsp,crsp,z,nfffff)
      if (keecur.gt.0.and.kdomse.eq.0) then
      nfffff = needs-kppcur-kffcur
      needer = needs
      call eecnst(ncrsp,crsp,z,nfffff)
      endif
!---------------------------------------------------------------------
!-- preconditioning A matrix if need                                --
!---------------------------------------------------------------------
      if (scalea) then
         call dgeequ(nj,need,arsp,nrsmat,rowscale,colscale, &
                     rowcnd,colcnd,arspmax,infosc)
         do j = 1,nj
         do i = 1,need
            arsp(j,i) = arsp(j,i) * colscale(i)
         enddo
         enddo
      endif
      if (ncrsp .le. 0) then
         call sdecm(arsp,nrsmat,nj,need,brsp,nrsmat, &
                      nnn,wrsp,work,ier)
         if (ier.ne.129) go to 2500
         write (nttyo,8000) ier
! MPI >>>
#if defined(USEMPI)
         ! ERROR_FIX >>>
         !kerror = 1
         !return
         ! <<< >>>
         call mpi_stop
         ! ERROR_FIX <<<
#else
         stop
#endif
! MPI <<<
 2500    continue
!-----------------------------------------------------------------------
!--  unfold fitting parameters                                        --
!-----------------------------------------------------------------------
         condno=wrsp(1)/wrsp(need)
         toler=condin*wrsp(1)
         do 2600 i=1,need
           t=0.0
          if (wrsp(i).gt.toler) t=brsp(i)/wrsp(i)
          work(i)=t
 2600   continue
        do 2650 i=1,need
          brsp(i)=0.0
          do 2650 j=1,need
            brsp(i)=brsp(i)+arsp(i,j)*work(j)
 2650   continue
      else
	do 2655 j=1,nrsmat
	    b(j) = brsp(j)
 2655	continue
	call dgglse(nj,need,ncrsp,arsp,nrsmat,crsp,4*(npcurn-2)+6+ &
                   npcurn*npcurn,b,z,brsp,work,nrsma2,info,condno)
	if (info.eq.0) goto 2656
	write (nttyo,8000) info
! MPI >>>
#if defined(USEMPI)
        ! ERROR_FIX >>>
        !kerror = 1
        !return
        ! <<< >>>
        call mpi_stop
        ! ERROR_FIX <<<
#else
	stop
#endif
! MPI <<<
 2656   continue
      endif
!----------------------------------------------------------------------
!--  rescale results if A is preconditioned                          --
!----------------------------------------------------------------------
      if (scalea) then
        do i=1, need
          brsp(i)=brsp(i)*colscale(i)
        enddo
      endif
      nload=nfnwcr
      if (ifitvs.gt.0) then
       if(nfourier.gt.1) then
        do j=1,nvesel
          temp=0.
         do 2671 i=1,(nfourier*2+1)
          temp=temp+brsp(i+nfnwcr)*vecta(i,j)
 2671    continue
         vcurrt(j)=temp
        enddo
       nload=nload+(nfourier*2+1)
       else
         do 2670 i=1,nvesel
          vcurrt(i)=brsp(i+nfnwcr)
 2670    continue
       nload=nload+nvesel
       endif
      endif
!
      if (kprfit.gt.0.and.kdofit.gt.0) then
        nload=nload+1
        prbdry=brsp(nload)
      endif
      if (fitdelz.and.nniter.ge.ndelzon) then
        nload=nload+1
        cdelz(nniter)=brsp(nload)*scaledz*100.*relaxdz
      else
        cdelz(nniter)=0.0
      endif
      if (iacoil.eq.0) then
          do 2711 i=1,nacoil
            caccurt(jtime,i)=0.0
 2711     continue
      else
          do 2717 i=1,nacoil
            caccurt(jtime,i)=brsp(nload+i)
 2717     continue
          nload=nload+nacoil
      endif
!------------------------------------------------------------------------------
!--  E coil currents                                                         --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
       do m=1,nesum
          nload=nload+1
          cecurr(m)=brsp(nload)
       enddo
      endif
!------------------------------------------------------------------------------
!--  Reference flux                                                          --
!------------------------------------------------------------------------------
      if (fitsiref) then
          nload=nload+1
          csiref=brsp(nload)*scalesir
      endif
!------------------------------------------------------------------------------
!--  Hyperbolic tangents                                                     --
!------------------------------------------------------------------------------
      if (kedgep.gt.0) then
        pedge=brsp(npedge)
      endif
      if (kedgef.gt.0) then
        f2edge=brsp(nfedge)
      endif
!
 3000 continue
!----------------------------------------------------------------------------
!-- Update poloidal flux due to vertical shift if necessary                --
!----------------------------------------------------------------------------
      if (ifitdelz.eq.2) then
      if (fitdelz.and.nniter.ge.ndelzon) then
      cdelznow=cdelz(nniter)/100.
      call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
      do i=1,nw
      do j=1,nh
        kk=(i-1)*nh+j
        call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n333)
        psi(kk)=psi(kk)+cdelznow*pds(3)
      enddo
      enddo
      endif
      endif
!----------------------------------------------------------------------
!-- calculate the fitting figure of merit saisq                      --
!----------------------------------------------------------------------
      saisq=0.0
      do 4600 m=1,nsilop
        cm=0.0
        do 4520 n=1,nfcoil
          cm=cm+rsilfc(m,n)*brsp(n)
 4520   continue
        if (ivesel.le.0) go to 4550
        do 4545 n=1,nvesel
          cm=cm+rsilvs(m,n)*vcurrt(n)
 4545   continue
 4550   continue
        if (iecurr.ne.1) go to 4565
        do 4560 n=1,nesum
          cm=cm+rsilec(m,n)*ecurrt(n)
 4560   continue
 4565   continue
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
        cmv=cm
        do 4540 n=1,kwcurn
          cm=cm+rsilpc(m,n)*brsp(nfcoil+n)
 4540   continue
        if (fitsiref) then
            cm=cm-csiref
        endif
        if (kedgep.gt.0) then
          cm=cm+rsilpe(m)*pedge
        endif
        if (kedgef.gt.0) then
          cm=cm+rsilfe(m)*f2edge
        endif
        if (swtsi(m).ne.0.0) then
        saisil(m)=fwtsi(m)**nsq*(silopt(jtime,m)-cm)**2
        saisil(m)=saisil(m)/swtsi(m)**nsq
        else
        saisil(m)=0.0
        endif
        saisq=saisq+saisil(m)
        csilop(m,jtime)=cm
        csilopv(m,jtime)=cmv
 4600 continue
!
      do 4700 m=1,magpri
        cm=0.0
        do 4620 n=1,nfcoil
          cm=cm+rmp2fc(m,n)*brsp(n)
 4620   continue
        if (ivesel.le.0) go to 4648
        do 4644 n=1,nvesel
          cm=cm+rmp2vs(m,n)*vcurrt(n)
 4644   continue
 4648   continue
        if (iecurr.ne.1) go to 4655
        do 4650 n=1,nesum
          cm=cm+rmp2ec(m,n)*ecurrt(n)
 4650   continue
 4655   continue
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
        cmv=cm
        do 4640 n=1,kwcurn
          cm=cm+rmp2pc(m,n)*brsp(nfcoil+n)
 4640   continue
        if (kedgep.gt.0) then
          cm=cm+rmp2pe(m)*pedge
        endif
        if (kedgef.gt.0) then
          cm=cm+rmp2fe(m)*f2edge
        endif
        if (swtmp2(m).ne.0.0) then
        saimpi(m)=fwtmp2(m)**nsq*(expmpi(jtime,m)-cm)**2
        saimpi(m)=saimpi(m)/swtmp2(m)**nsq
        else
        saimpi(m)=0.0
        endif
        saisq=saisq+saimpi(m)
        cmpr2(m,jtime)=cm
        cmpr2v(m,jtime)=cmv
 4700 continue
!
      if (kstark.gt.0) then
      chigamt=0.0
      do 4800 m=1,nstark
        chigam(m)=0.0
        cmgam(m,jtime)=0.0
        if (rrgam(jtime,m).le.0.0) goto 4800
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
        if (iacoil.le.0) go to 4761
        do 4759 n=1,nacoil
          cmbr=cmbr+rbrac(m,n)*caccurt(jtime,n)
          cmbz=cmbz+rbzac(m,n)*caccurt(jtime,n)
 4759   continue
 4761   continue
        if (iecurr.eq.2) then
          do n=1,nesum
            cmbr=cmbr+rbrec(m,n)*cecurr(n)
            cmbz=cmbz+rbzec(m,n)*cecurr(n)
          enddo
        endif
        if (kedgep.gt.0) then
          cmbr=cmbr+rbrpe(m)*pedge
          cmbz=cmbz+rbzpe(m)*pedge
        endif
        if (kedgef.gt.0) then
          cmbr=cmbr+rbrfe(m)*f2edge
          cmbz=cmbz+rbzfe(m)*f2edge
        endif
        cm=a2gam(jtime,m)*btgam(m)+a3gam(jtime,m)*cmbr+a4gam(jtime,m) &
             *cmbz
        bzmsec(m)=cmbz
        if (keecur.le.0) then
          bzmse(m)=(tangam(jtime,m)*cm-a8gam(jtime,m)*cmbr) &
                       /a1gam(jtime,m)
          cm=(a1gam(jtime,m)*cmbz+a8gam(jtime,m)*cmbr)/cm
        else
          erup=0.0
          erbot=0.0
          do i=1,keecur
            cerer(i)=brsp(needer+i)
            erup=erup+e1rbz(m,i)*cerer(i)
            erbot=erbot+(e2rbz(m,i)+e3rbr(m,i))*cerer(i)
          enddo
          cm= cm-erbot
          bzmse(m)=(tangam(jtime,m)*cm+erup-a8gam(jtime,m)*cmbr) &
                      /a1gam(jtime,m)
          cm=(a1gam(jtime,m)*cmbz+a8gam(jtime,m)*cmbr-erup)/cm
          ermse(m)=-erup/a5gam(jtime,m)
        endif
        if (swtgam(m).ne.0.0) then
        chigam(m)=fwtgam(m)**nsq*(tangam(jtime,m)-cm)**2
        chigam(m)=chigam(m)/swtgam(m)**nsq
        else
        chigam(m)=0.0
        endif
        cmgam(m,jtime)=cm
 4800 continue
      do m=1,nmtark
        chigamt=chigamt+chigam(m)
      enddo
      chilibt=0.0
      do m=nmtark+1,nstark
        chilibt=chilibt+chigam(m)
      enddo
      if (ishot.le.97400) then
        mcentral=15
      else
        mcentral=10
      endif
      do m=1,mcentral
        mp1=m+1
        drgam=rrgam(jtime,mp1)-rrgam(jtime,m)
        cjmse(m)=-(bzmse(mp1)-bzmse(m))/drgam/twopi/tmu
        cjmsec(m)=-(bzmsec(mp1)-bzmsec(m))/drgam/twopi/tmu
      enddo
      endif
!
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
!
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
      do 4900 n=nfcoil+1,nfnwcr
        cm=cm+brsp(n)*fgowpc(n-nfcoil)
 4900 continue
        if (kedgep.gt.0) then
          cm=cm+fgowpe*pedge
        endif
        if (kedgef.gt.0) then
          cm=cm+fgowfe*f2edge
        endif
      cpasma(jtime)=cm
      if (kfffnc.eq.8) then
        cjeccd=brsp(nfnwcr)*fgowpc(nfnwcr-nfcoil)/1000.
      else
        cjeccd=0.0
      endif
      if (ifitvs.gt.0) then
        do 30500 i=1,nvesel
          cm=cm+vcurrt(i)
30500   continue
      endif
      if (swtcur.ne.0.0) then
      saiip=(fwtcur/swtcur)**nsq*(pasmat(jtime)-cm)**2
      else
      saiip=0.0
      endif
      saisq=saisq+saiip
!
      tsaifc=0.0
      do 30510 i=1,nfcoil
        saifc(i)=0.0
        if (fwtfc(i).gt.0.0) then
          saifc(i)=fwtfc(i)**nsq*(brsp(i)-fccurt(jtime,i))**2
          saifc(i)=saifc(i)/swtfc(i)**nsq
        endif
        saisq=saisq+saifc(i)
        tsaifc=tsaifc+saifc(i)
30510 continue
!
      if (iecurr.eq.2) then
      do i=1,nesum
        saiec(i)=0.0
        if (fwtec(i).gt.0.0) then
          saiec(i)=fwtec(i)**nsq*(cecurr(i)-ecurrt(i))**2
          saiec(i)=saiec(i)/swtec(i)**nsq
        endif
        saisq=saisq+saiec(i)
      enddo
      endif
!
      if (fitsiref) then
        saisref=0.0
        if (fwtref.gt.0.0) then
          saisref=fwtref**nsq*(psiref(jtime)-csiref)**2
          saisref=saisref/swtsi(nslref)**nsq
        endif
        saisq=saisq+saisref
      endif
      tsaisq(jtime)=saisq
!
      chipre=0.0
      if (kprfit.le.0.or.kdofit.eq.0) go to 4910
      if (npress.le.0) goto 4910
      do 4908 m=1,npress
        cm=0.0
        do 4906 n=1,kppcur
          cm=cm+rprepc(m,n)*brsp(nfcoil+n)
 4906   continue
        if (kedgep.gt.0) then
          cm=cm+rprepe(m)*pedge
        endif
        cm=cm+prbdry
        if (fwtpre(m).gt.0.0) then
        saipre(m)=((cm-pressr(m))/sigpre(m))**2
        else
        saipre(m)=0.0
        endif
        saipre2(m)=saipre(m)  ! preserve saipre - changed later by pltout
        chipre=chipre+saipre(m)
        precal(m)=cm
 4908 continue
 4910 continue
!
      chiprw=0.0
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          cm=0.0
          do n=1,kwwcur
            cm=cm+rprwpc(m,n)*brsp(nfnpcr+n)
          enddo
        cm=cm+preswb
        if (fwtprw(m).gt.0.0) then
        saiprw(m)=((cm-presw(m))/sigprw(m))**2
        else
        saiprw(m)=0.0
        endif
        chiprw=chiprw+saiprw(m)
        prwcal(m)=cm
        enddo
      endif
!
      if ((nniter.lt.minite).and.((eouter.gt.elomin).or.(fwtdlc.gt.0.0) &
      ))  go to 4950
      if ((nniter.lt.kcallece).and.(kfitece.gt.0.0))  go to 4950
      if ((errorm.gt.errmin).and.((eouter.gt.elomin).or.(fwtdlc.gt.0.0) &
      ))  go to 4950
      if (saisq.gt.saicon) go to 4950
      if (iconvr.ne.2) go to 4950
      if (abs(saisq-saiold).le.0.10) go to 4918
      if (saisq.lt.saiold) go to 4950
 4918 continue
      ichisq=1
      do 4920 i=1,nrsmat
        brsp(i)=brsold(i)
 4920 continue
      if (ifitvs.gt.0) then
        do 30520 i=1,nvesel
          vcurrt(i)=vcurrto(i)
30520   continue
      endif
      saisq=saiold
 4950 continue
!
      if (iand(iout,1).ne.0) then
      write (nout,7400) time(jtime),chipre,cpasma(jtime), &
                        nniter,condno,saisq,chigamt
      write (nout,7445) need
      write (nout,7450) (brsp(i),i=1,need)
      write (nout,7450) (wrsp(i),i=1,need)
      endif
!     
      if (iupdat.gt.0) return
      if (saisq.gt.saimin) return
      tcrrnt=cpasma(jtime)
      iupdat=1
      return
!
 6000 continue
      return
 7400 format (/,2x,7htime = ,e12.5,2x,8hchipr = ,e12.5, &
              2x,10hcurrent = ,e12.5,/,2x,5hit = ,i5, &
              1x,8hcondn = ,1pe11.4, &
              1x,8hchisq = ,1pe11.4,1x,9hchigam = ,1pe11.4, &
              /,2x,11hchiecebz = ,1pe11.4,1x,14htchieceR+R- = , &
              1pe11.4)
 7445 format (10x,22hfitting parameters:   ,i5)
 7450 format (8(1x,e12.5,1x))
 7460 format (10x,7hchi ip:,/,15x,e12.5)
 8000 format (/,'  ** Problem in Decomposition **',i10)
      end

      subroutine pflux(niter,nnin,ntotal,jtime)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          pflux computes the poloidal fluxes on the r-z grid.     **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**          09/03/90..........Lazarus vertical feedback             **
!**          28/03/93..........fix relax                             **
!**                                                                  **
!**********************************************************************
!vas  f90 modifi
      use var_bunemn
      use commonblocks,only: c,wk,copy,bkx,bky,psiold,psipold, &
                             psipp,work
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      common/cwork3/lkx,lky
      integer initresult
      dimension pds(6)
      real*8,dimension(:),allocatable :: psikkk,gfbsum
      data initfb/0/,init/0/
!
      ALLOCATE(psikkk(nwnh),gfbsum(nwnh))
!
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
          mj=iabs(j-jj)+1
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
		if(nw.gt.30)ioffr=ioffr*((nw+1)/33)
		if(nh.gt.30)ioffz=ioffz*((nh+1)/33)
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
      if (abs(vcurfb(1)).gt.1.e-6) then
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
          mj=iabs(j-jj)+1
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
          mj1=iabs(jj-1)+1
          mjnh=iabs(nh-jj)+1
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
        call pflux_cycred(psi,work)
      endif
      do 3000 i=1,nwnh
        psi(i)=-psi(i)
 3000 continue

 3010 continue
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
      if (abs(vcurfb(1)).lt.1.e-6) go to 43000
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
          mj=iabs(j-jj)+1
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
          mj=iabs(j-jj)+1
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
      if (abs(relax-1.0).lt.1.0e-03) go to 7000
      do 6500 kk=1,nwnh
        psi(kk)=relax*psi(kk)+(1.-relax)*psiold(kk)
        psipla(kk)=relax*psipla(kk)+(1.-relax)*psipold(kk)
 6500 continue
 7000 continue
!
      DEALLOCATE(psikkk,gfbsum)
!
      return
      end
! The following set of routines are specific to the fast cyclic solver
! originated by Holger St. John, optimised to realtime format by John Ferron
! and then modified for use as an option here by Dylan Brennan.  These
! routines continue up to the mark END_CYCLIC_ROUTINES

! ======================================================================
! FUNCTION ef_init_cycred_data
! Initialize a set of precomputed data for the cyclic reduction algorithm.
! These data depend on the grid size and spacings.
! ----------------------------------------------------------------------

      integer function ef_init_cycred_data()

      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'

      u0 = 4.0 * 3.1415927 * 1.0e-7

! ----------------------------------------------------------------------
! Power of 2 that specifies the grid height.
      i=1
      nhpwr = -1
      do 10 j=0,14
        i = i*2
        if (i.eq.(nh - 1)) then
          nhpwr = j+1
          goto 15
        endif
   10 continue

   15 if (nhpwr.eq.-1) then
         write(6,*) ' (init_cycred_data): grid height must be a', &
                    ' power of 2 + 1.\n PROGRAMMING ERROR!\n'
!        default to nh=33 
         nhpwr = 5
         ef_init_cycred_data=1 
         return 
      endif  

! ----------------------------------------------------------------------
! Elements of the tridiagonal matrix.
! todo: replace these with the common block couterparts

      dr = rgrid(2) - rgrid(1)
      dz = zgrid(2) - zgrid(1)

      dzsq = dz*dz
      dzdrsq = (dz/dr) * (dz/dr)
      dumy = dzsq/(2.0*dr)
      dumy1 = 2.0 * (1.0 + dzdrsq)

!     All elements of the matrix diagonal have the same value
      diag = dumy1 
      do 20 i=0,nw-3
         denom = dumy/rgrid(i+2)
         diagl(i+1) = -1.0 * dzdrsq - denom
         diagu(i+1) = -1.0 * dzdrsq + denom
   20 continue 

! ----------------------------------------------------------------------
! Misc. constants used in computing the rhs vector.
!
! 
      rhs_a_dumy = dzdrsq + dzsq/(2.0 * rgrid(2) * dr)
      rhs_b_dumy = dzdrsq - dzsq/(2.0 * rgrid(nw-1) * dr)

! ----------------------------------------------------------------------
! Constants used during the forward reduction procedure.

      index = 0
!     nred steps are required
      nred = nhpwr - 1    
      jpow = 1   
!     k is the reduction step index
      do 90 k=1,nred   
!                2**(k-1)
        jpowm1 = jpow   
!                2**k
        jpow   = jpow * 2 
!                2**k + 1
        jstart = jpow + 1   
!                nh - 2**k
        jend   = nh - jpow 
!                2**k 
        jstep  = jpow    

        do 80 j=jstart,jend,jstep
           m1 = -1.0
           do 70 l=1,jpowm1
              m1 = -1.0 * m1
 
              index = index + 1
              if (index.gt.icycred_loopmax) then
                 write(6,*) 'rtefit (init_cycred_data):', &
                            'constant data index is ',index, &
                            'too large! PROGRAMMING ERROR! id:1', &
                            ' ',jstart,' ',jend,' ',jstep,' ',j, &
                            ' ',nred,' ',k,' ',jpowm1,' ',l
                 ef_init_cycred_data=1
                 return 
              endif 

              cosdii = -2.0 * cos( (2.0 * l - 1.0) * pi/jpow)
              sindi = sin( (2.0 * l - 1.0) * pi/jpow)
              alphab(index) = m1 * sindi/jpowm1
              diag1(index) = diag + cosdii
   70     continue 

          if (index.gt.icycred_loopmax) goto 90
              
   80   continue
   90 continue

! ----------------------------------------------------------------------
! Constants used during the back solution procedure.

!            nhpwr
      nred = nred + 1 
!            2**nhpwr
      jpow = 2**nred 

      do 120 k=nred,1,-1
!               2**k
        jstep = jpow
!              2**(k-1)
        jpow = jpow/2   
        jstart = jpow + 1
        jend = nh - jpow
        do 110 j=jstart,jend,jstep
          m1 = -1.0
          do 100 l=1,jpow
            m1 = -1.0 * m1

            index = index + 1
            if (index.gt.icycred_loopmax) then
              write(6,*) 'rtefit (init_cycred_data): ', &
                         'constant data index is ',index, &
                         'too large! PROGRAMMING ERROR! id:2'
              ef_init_cycred_data=1
              return
            endif 

            cosdii = -2.0 * cos( (2.0 * l - 1.0) * pi/(2**k)) 
            sindi = sin( (2.0 * l - 1.0) * pi/(2**k))
            alphab(index) = m1 * sindi/jpow
            diag1(index) = diag + cosdii
  100     continue 

          if(index.gt.icycred_loopmax) goto 120
             
  110   continue
  120 continue   

! ----------------------------------------------------------------------
! Vectors of precalculated values used by the tridiag routine.
! These vectors are various combinations of the diagonals of the matrices
! that are solved.

! At this point "index" holds the number of reduction loops that are used.
      do 130 i=1,index
        beti(i,1) = 1.0/diag1(i)
!       not actually used
        abeti(i,1) = 0.0 
!       not actually used
        wk1(i,1) = 0.0   

        do 130 j=2,nw-2
          wk1(i,j) = diagu(j-1) * beti(i,j-1)
          beti(i,j) =  &
            1.0/( diag1(i) - diagl(j) * wk1(i,j))
          abeti(i,j) =  &
            diagl(j) * beti(i,j)
  130 continue    

! ----------------------------------------------------------------------
! All done.

      ef_init_cycred_data=0
      end



! ======================================================================
! FUNCTION cyclic_reduction
! The core routine to compute the flux on the grid using Holger StJohn's
! single cyclic reduction algorithm.
! ----------------------------------------------------------------------

      subroutine cyclic_reduction(f)

      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'

      dimension f(1)

      constant=0.5
! ----------------------------------------------------------------------
! Forward reduction.

      nv = nw - 2

      index = 0 
!     nred steps are required.
      nred = nhpwr - 1
      jpow = 1   
      do 5 i=0,nv
        wk2(i+1) = 0.0
    5 continue

!     k is the reduction step index. 
      do 20 k=1,nred   
!                2**(k-1)
        jpowm1 = jpow 
!                2**k
        jpow   = jpow * 2  
!                2**k + 1 
        jstart = jpow + 1
!                nh - 2**k 
        jend   = nh - jpow 
!                2**k     
        jstep  = jpow  

        do 20 j=jstart,jend,jstep
       
!         Index of the first element of jth row 
          jd = (j-1) * nw + 1
!         Next row up
          jdp = jd + jpowm1*nw  
!         Next row down
          jdm = jd - jpowm1*nw 

          call ef_vadd_shrt(f(jdm+1),f(jdp+1),phi,nv)
          
          do 10 l=1,jpowm1
            index = index + 1
            call ef_tridiag2(f(jd+1),nv,index)
   10     continue

!         use the declared var constant, instead of 0.5 directly because
!         of change in alignment for different compilations

          call ef_vmul_const_shrt(f(jd+1),constant,f(jd+1),nv)

   20 continue

! ----------------------------------------------------------------------
! Back solution.

!            nhpwr
      nred = nred + 1
!            2**nhpwr 
      jpow = 2**nred 

      do 30 k=nred,1,-1
!               2**k
        jstep = jpow    
!              2**(k-1) 
        jpow = jpow/2  
        jstart = jpow + 1
        jend = nh - jpow

        do 30 j=jstart,jend,jstep
!         Index of the first element of jth row 
          jd = (j-1) * nw   
!         Next row up 
          jdp = jd + jpow*nw 
!         Next row down 
          jdm = jd - jpow*nw  
 
          call ef_vadd_shrt(f(jdm+2),f(jdp+2),phi,nv)
          do 30 l=1,jpow
            index = index + 1

            if (l.eq.1) then
              call ef_tridiag1(f(jd+2),nv,index)
            else
              call ef_tridiag2(f(jd+2),nv,index)
            endif
   30 continue

! All done.

      return
      end

! ======================================================================
! FUNCTION pflux_cycred
! Use the single cyclic reduction method of Holger StJohn and John Ferron 
! to get the flux on the grid resulting from the plasma current.
! ----------------------------------------------------------------------

      subroutine pflux_cycred(psigrid,sia)

      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'

      dimension psigrid(1),sia(1)

! ----------------------------------------------------------------------
! Call the initialisation and check the result
      initresult=ef_init_cycred_data()
      if (initresult.eq.1) then
! MPI >>>
#if defined(USEMPI)
        ! ERROR_FIX >>>
        !return
        ! <<< >>>
        call mpi_stop
        ! ERROR_FIX <<<
#else
        stop 'init failure ef_init_cycred_data'
#endif
! MPI <<<
      endif

! ----------------------------------------------------------------------
! transpose the matrix, solver needs transposition (same as buneto)

      do 2 i = 1,nw
        ii = (i-1)*nh 
        do 2 j = 1,nh
          sia(i+(j-1)*nw) = psigrid(ii+j)
    2 continue

! ----------------------------------------------------------------------
! Add the finite difference expressions for the second and nw-1 grid columns
! to the rhs. These terms are known boundary terms and
! hence do not form part of the tridiagonal matrix.

      ius = nw
      ivs = nw
      iys = nw
      nn = nh-2
      call vsma_(sia(nw+1), ius, &
          rhs_a_dumy, &
          sia(nw+2), ivs, &
          sia(nw+2), iys, &
          nn)
      call vsma_(sia(2*nw), ius, &
          rhs_b_dumy, &
          sia(2*nw-1), ivs, &
          sia(2*nw-1), iys, &
          nn)
! ----------------------------------------------------------------------
! Do the cyclic reduction to get psigrid.

      call cyclic_reduction(sia)

! ----------------------------------------------------------------------
! transpose the matrix back

      do 6 i = 1,nh
        ii = (i-1)*nw
        do 6 j = 1,nw
          psigrid(i+(j-1)*nh) = sia(ii+j)
    6 continue

! All done.

      return
      end

! ======================================================================
! ======================================================================
!  +++++++++++++++++++++++
!  SUBROUTINE: vsma_
!  +++++++++++++++++++++++
!
      subroutine vsma_(a, ia, b, c, ic, d, id, n)
      implicit integer*4 (i-n), real*8 (a-h, o-z)

      dimension a(1),c(1),d(1)

      if (n.le.0) then
        return
      endif

      iai = 1
      ici = 1
      idi = 1

      do 10 i=1,n 
        d(idi) = a(iai)*b + c(ici)
        iai = iai + ia
        ici = ici + ic
        idi = idi + id
   10 continue
      return
      end


! ef_vvmul: multiply two vectors together
! ef_tridiag1
! ef_tridiag2: triagonal matrix solution routines used by the cyclic
!              reduction algorithm for computing flux on the grid.
!
! ef_vadd_shrt: Add two vectors.  
!               Handles vectors one element at a time.  Good for
!      short, unaligned vectors but not optimized at all for long vectors.
! ef_vmul_const_shrt: Multiply a vector by a constant. 
!            Handles vectors one element
!      at a time.  Good for
!      short, unaligned vectors but not optimized at all for long vectors.
!----------------------------------------------------------------------
! Modifications:
! 12/30/97: created by J. Ferron from the C language equivalent comments
!           in rtefitutils.s
! 3/31/00:  converted to fortran by Dylan Brennan
!**********************************************************************


      subroutine ef_vvmul(vin1,vin2,out,nelements)
      implicit integer*4 (i-n), real*8 (a-h, o-z)

      dimension vin1(1),vin2(1),out(1)

      do 10 i=1,nelements
        out(i) = vin1(i) * vin2(i)
   10 continue
      return
      end


! ======================================================================
! ======================================================================
!  +++++++++++++++++++++++
!  SUBROUTINE: ef_tridiag2
!  +++++++++++++++++++++++
! 
!  Solve system of equations defined by a tridiagonal matrix.
!  This is a version of the routine in Numerical Recipies which
!  uses some values precalculated from the diagonals of the matrix.
! 
!  This routine also does a bunch of preparation to compute the rhs vector
!  and accumulate the result.
! 
!  con: array of structures giving precalculated values
!  wk1: vector of precalculated values. 
!  wk2,phi: vectors that when combined yield the right hand side vector
!  alphab: a constant that would be used in computing the right hand side
!         vector.
!  f: the vector in which the result is being accumulated.
!  v: temporary storage area for the result vector
!  n: number of elements in the vector. 
! 

      subroutine ef_tridiag2(f,n,index)

      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'

      dimension f(1)

      v(1)=(wk2(1)+alphab(index)*phi(1))* &
           beti(index,1)
      do 10 j=2,n
        v(j)=wk2(j)*beti(index,j)+  &
          alphab(index)*phi(j)*beti(index,j)-  &
          v(j-1) * abeti(index,j)
   10 continue

      usave = v(n)
      f(n) = usave + f(n)

      do 20 j=n-1,1,-1
        usave = v(j)-wk1(index,j+1)*usave
        f(j) = usave + f(j)
   20 continue 
      return
      end


! ======================================================================
! ======================================================================
!  +++++++++++++++++++++++
!  SUBROUTINE: ef_tridiag1
!  +++++++++++++++++++++++
! 
!  Solve system of equations defined by a tridiagonal matrix.
!  This is a version of the routine in Numerical Recipies which
!  uses some values precalculated from the diagonals of the matrix.
! 
!  This routine also does a bunch of preparation to compute the rhs vector
!  and accumulate the result.
! 
!  Basically the same as tridiag2 except:
!  f is the input instead of wk2.  At the end of the routine, f as it
!  was at the start of the routine has been copied into wk2.
! 
!  The result is written directly into f rather than adding the result
!  to the input value of f.
! 
!  con: array of structures giving precalculated values.
!  wk1: vector of precalculated values. 
!  wk2,phi: vectors that when combined yield the right hand side vector
!  alphab: a constant that would be used in computing the right hand side
!         vector.
!  f: the vector in which the result is being accumulated.
!  v: temporary storage area for the result vector
!  n: number of elements in the vector. 
! 

      subroutine ef_tridiag1(f,n,index)

      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'

      dimension f(1)
 
      v(1)=(f(1)+alphab(index)*phi(1))* &
                    beti(index,1)
 
      do 10 j=2,n
       v(j) = f(j)*beti(index,j)+  &
         alphab(index)*phi(j)*beti(index,j)-  &
         v(j-1)*abeti(index,j)
   10 continue
 
      wk2(n) = f(n)
      f(n) = v(n)
 
      do 20 j=n-1,1,-1
        wk2(j) = f(j)
        v(j)=v(j)-wk1(index,j+1)*v(j+1)
        f(j) = v(j)
   20 continue 
      return
      end


! ======================================================================
! ======================================================================
!  ++++++++++++++++++++++++
!  SUBROUTINE: ef_vadd_shrt
!  ++++++++++++++++++++++++
! 
!  Add two vectors.  Handles vectors one element at a time.  Good for
!  short, unaligned vectors but not optimized at all for long vectors.
! 
!  vector_out = vector1 + vector2
! 
!  vector1 = input vector
!  vector2 = input vector
!  vector_out = output vector
!  nelements = number of elements in the vectors.

      subroutine ef_vadd_shrt(vector1,vector2,vector_out, &
                              nelements)
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension vector1(1),vector2(1),vector_out(1)

      do 10 i=1,nelements
        vector_out(i) = vector1(i) + vector2(i)
   10 continue
      return
      end


! ======================================================================
! ======================================================================
!  ++++++++++++++++++++++++
!  SUBROUTINE: ef_vmul_const_shrt
!  ++++++++++++++++++++++++
! 
!  Multiply a vector by a constant.
!  Handles vectors one element at a time.  Good for
!  short, unaligned vectors but not optimized at all for long vectors.
! 
!  vector_out = vector1 * constant
! 
!  vector1 = input vector
!  constant = constant value
!  vector_out = output vector
!  nelements = number of elements in the vectors. Must be at least 2.

      subroutine ef_vmul_const_shrt(vector1,constant,vector_out, &
                                    nelements)
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension vector1(1),vector_out(1)

      do 10 i=1,nelements
        vector_out(i) = vector1(i) * constant
   10 continue
      return
      end

!  END_CYCLIC_ROUTINES This ends the fast cyclic reduction routines

      function prcur4(n1set,ypsi,nnn)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          prcur4 computes the plasma pressure by integration.     **
!**                                                                  **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          94/03/01..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: sifpre,bwpre,cwpre,dwpre,sfpre,sprep
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!
      if (abs(ypsi).gt.1.0) then
        prcur4=0.0
        return
      endif
!

      if (n1set.gt.0) then
       do i=2,nw-1
        siii=float(i-1)/float(nw-1)
        sifpre(i)=siii
       enddo
       sifpre(1)=0.0
       sifpre(nw)=1.0
       do i=1,nw
        sprep(i)=ppcur4(sifpre(i),kppcur)/darea
       enddo
!
       sfpre(nw)=prbdry
       delsi=sidif/float(nw-1)
       do 1000 i=1,nw-1
       sfpre(nw-i)=sfpre(nw-i+1)+0.5*(sprep(nw-i+1)+sprep(nw-i))*delsi
 1000  continue
!
       mw=nw
       call zpline(mw,sifpre,sfpre,bwpre,cwpre,dwpre)
      endif
      prcur4=seval(mw,ypsi,sifpre,sfpre,bwpre,cwpre,dwpre)
      return
      end
      function ppcur4(ypsi,nnn)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          ppcur4 computes the radial derivative of the            **
!**          pressure based on icurrt.                               **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          25/07/85..........revised                               **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!
      if (abs(ypsi).gt.1.0) then
        ppcur4=0.0
        return
      endif
      if (icurrt.eq.4) goto 2000
      if (icurrt.eq.1) goto 3000
      ppcur4=ppcurr(ypsi,nnn)
      return
!
 2000 continue
      ppcur4= ((1.-ypsi**enp)**emp*(1.-gammap)+gammap)*cratio &
                /rzero
      return
 3000 continue
      ppcur4= cratio*sbeta/srma
      return
      end
      function ppcurr(ypsi,nnn)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          ppcurr computes the radial derivative of the            **
!**          pressure.                                               **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          25/07/85..........revised                               **
!**          94/03/11..........revised                               **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
      dimension xpsii(nppcur)
!
! jm.s
      real*8, external :: linear
! jm.e
      if (abs(ypsi).gt.1.0) then
        ppcurr=0.0
        return
      endif
! jm.s
      if (npsi_ext > 0) then
!        ppcurr = linear(ypsi,psin_ext,pprime_ext,npsi_ext)
        ppcurr = seval(npsi_ext,ypsi,psin_ext,pprime_ext,bpp_ext,cpp_ext,dpp_ext)
        ppcurr = ppcurr * cratiop_ext
        return
      endif
! jm.e
      ppcurr=0.0
      call setpp(ypsi,xpsii)
      do 1400 iiij=nfcoil+1,nnn+nfcoil
        iijj=iiij-nfcoil
        ppcurr=ppcurr+brsp(iiij)*xpsii(iijj)
 1400 continue
!----------------------------------------------------------------------
!-- edge hyperbolic tangent component                                --
!----------------------------------------------------------------------
      if (kedgep.eq.0) return
      siedge=(ypsi-pe_psin)/pe_width
      p0back=pedge/pe_width/sidif
      ppcurr=ppcurr+p0back/cosh(siedge)**2
      return
!
      entry prcurr(ypsi,nnn)
      if (abs(ypsi).gt.1.0) then
        prcurr=0.0
        return
      endif
      brspp=0.0
      prcurr=0.0
      call setpr(ypsi,xpsii)
      do 1600 i=nfcoil+1,nfcoil+nnn
        nn=i-nfcoil
        prcurr=prcurr+brsp(i)*xpsii(nn)
 1600 continue
      prcurr=-sidif*prcurr/darea+prbdry
!----------------------------------------------------------------------
!-- edge hyperbolic tangent component                                --
!----------------------------------------------------------------------
      if (kedgep.eq.0) return
      siedge=(ypsi-pe_psin)/pe_width
      prcurr=prcurr+pedge/darea*(tpedge-tanh(siedge))
      return
      end
      subroutine presurw(jtime,niter)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          presurw computes the relevant parameters for rotational **
!**          pressure profile fitting.                               **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          94/03/11..........first created                         **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      common/cwork3/lkx,lky
      dimension pds(6),ybase(nwwcur),ybaseb(nwwcur)
      dimension bwork(ndata),cwork(ndata),dwork(ndata)
      data init/0/
      save init
!
      kdofit=1
      go to (1000,2000) kprfit-2
 1000 continue
      if (nmass*nomegat.gt.0) goto 1300
!------------------------------------------------------------------------
!--  specify the rotational pressure profile directly                  --
!------------------------------------------------------------------------
      do 1200 i=1,npresw
        xn=-rpresw(i)
        if (rpresw(i).le.0.0) go to 1030
        call seva2d(bkx,lkx,bky,lky,c,rpresw(i),zpresw(i),pds,ier,n333)
        xn=(simag-pds(1))/sidif
 1030   continue
        rpresws(i)=xn
        call setpw(xn,ybase)
        call setpwp(xn,ybaseb)
        xbaseb=ybaseb(kwwcur)*xn**2
        do 1150 m=1,kwwcur
          rprwpc(i,m)=-sidif*ybase(m)/darea
 1150   continue
!----------------------------------------------------------------------
!-- response for DELTAZ                                              --
!----------------------------------------------------------------------
        if (fitdelz.and.niter.ge.ndelzon) then
          if (rpresw(i).le.0.0) then
          rprwdz(i)=0.0
          else
          rprwdz(i)=pds(3)*pwpcur(xn,kwwcur)/darea
          endif
        endif
 1200 continue
      return
 1300 continue
!------------------------------------------------------------------------
!-- form rotational pressure from mass density and rotaional frequency --
!------------------------------------------------------------------------
      if (init.eq.0) then
!------------------------------------------------------------------------
!--  set up interpolation                                              --
!------------------------------------------------------------------------
        call zpline(nmass,sibeam,dmass,bwork,cwork,dwork)
        init=1
      endif
      do i=1,nomegat
        xn=-romegat(i)
        if (romegat(i).le.0.0) go to 1330
        rnow=romegat(i)
        znow=zomegat(i)
        call seva2d(bkx,lkx,bky,lky,c,rnow,znow,pds,ier,n333)
        xn=(simag-pds(1))/sidif
 1330   continue
        rpresw(i)=-xn
        rpresws(i)=xn
        dmnow=seval(nmass,xn,sibeam,dmass,bwork,cwork,dwork)
        presw(i)=dmnow*omegat(i)*rvtor**2
        sigprw(i)=abs(presw(i))*sigome(i)
        presw(i)=0.5*presw(i)*omegat(i)
        call setpw(xn,ybase)
        call setpwp(xn,ybaseb)
        xbaseb=ybaseb(kwwcur)*xn**2
        do m=1,kwwcur
          rprwpc(i,m)=-sidif*ybase(m)/darea
        enddo
!----------------------------------------------------------------------
!-- response for DELTAZ                                              --
!----------------------------------------------------------------------
        if (fitdelz.and.niter.ge.ndelzon) then
          if (romegat(i).le.0.0) then
          rprwdz(i)=0.0
          else
          rprwdz(i)=pds(3)*pwpcur(xn,kwwcur)/darea
          endif
        endif
      enddo
      npresw=nomegat
      return
!------------------------------------------------------------------------
!--  construct rotational pressure from kinetic data                   --
!------------------------------------------------------------------------
 2000 continue
      return
      end
      subroutine presur(jtime,niter,kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          presur computes the relevant parameters for pressure    **
!**          profile fitting.                                        **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          09/09/85..........first created                         **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      common/cwork3/lkx,lky
      dimension pds(6),xnsi(nppcur),xnsp(nppcur)
      character*50 edatname
      namelist/edat/npress,rpress,zpress,pressr,sigpre
! MPI >>>
      integer, intent(inout) :: kerror
      kerror = 0
! MPI <<<
      kdofit=1
      go to (1000,2000,1000) kprfit
 1000 continue
!---------------------------------------------------------------------
!--  input pressure profile                                         --
!---------------------------------------------------------------------
      do 1200 i=1,npress
        xn=-rpress(i)
        if (rpress(i).le.0.0) go to 1030
        call seva2d(bkx,lkx,bky,lky,c,rpress(i),zpress(i),pds,ier,n333)
        xn=(simag-pds(1))/sidif
 1030   continue
        call setpr(xn,xnsi)
        call setpp(xn,xnsp)
        do 1150 m=1,kppcur
          rprepc(i,m)=-sidif*xnsi(m)/darea
 1150   continue
!----------------------------------------------------------------------
!-- response for hyperbolic tangent component                        --
!----------------------------------------------------------------------
        if (kedgep.gt.0) then
          siedge=(xn-pe_psin)/pe_width
          rprepe(i)=(tpedge-tanh(siedge))/darea
        endif
!----------------------------------------------------------------------
!-- response for DELTAZ                                              --
!----------------------------------------------------------------------
        if (fitdelz.and.niter.ge.ndelzon) then
          if (rpress(i).le.0.0) then
          rpredz(i)=0.0
          else
          rpredz(i)=pds(3)*ppcurr(xn,kppcur)/darea
          endif
        endif
 1200 continue
      return
!----------------------------------------------------------------
!--  construct pressure from kinetic data                      --
!----------------------------------------------------------------
 2000 continue
!----------------------------------------------------------------
!--  use psi grid defined by Thompson data                     --
!----------------------------------------------------------------
      npress=npteth
      do 2200 i=1,npteth
        call seva2d(bkx,lkx,bky,lky,c,rteth(i),zteth(i),pds,ier,n111)
        xn=(simag-pds(1))/sidif
        if (xn.ge.1.0) then
          npress=i     -1
          go to 2220
        endif
        rpress(i)=-xn
        call setpr(xn,xnsi)
        call setpp(xn,xnsp)
        do 2150 m=1,kppcur
          rprepc(i,m)=-sidif*xnsi(m)/darea
 2150   continue
!----------------------------------------------------------------------
!-- response for hyperbolic tangent component                        --
!----------------------------------------------------------------------
        if (kedgep.gt.0) then
          siedge=(xn-pe_psin)/pe_width
          rprepe(i)=(tpedge-tanh(siedge))/darea
        endif
 2200 continue
 2220 continue
!---------------------------------------------------------------
!--  npress=0 can cause problems                              --
!---------------------------------------------------------------
      if (npress.le.0) then
        kdofit=0
        return
      endif
!----------------------------------------------------------------
!--  get ion temperature and density profile                   --
!----------------------------------------------------------------
      if (nptef.ne.0) then
! MPI >>>
        call gette(kerror)
#if defined(USEMPI)
        !if (kerror /= 0) then
        !  kerror = 1
        !  return
        !endif
#endif
! MPI <<<
      endif

! MPI >>>
      call getne(jtime,kerror)
#if defined(USEMPI)
      !if (kerror /= 0) then
      !  kerror = 1
      !  return
      !endif
#endif
! MPI <<<

! MPI >>>
      call gettion(kerror)
#if defined(USEMPI)
      !if (kerror /= 0) then
      !  kerror = 1
      !  return
      !endif
#endif
! MPI <<<
      if (nbeam.ne.0) call getbeam
!----------------------------------------------------------------
!--  construct pressure profile                                --
!----------------------------------------------------------------
      pressb=1.602e+03*(dibdry*tibdry+debdry*tebdry)+pbeamb
      prespb=1.602e+03*(dipbry*tibdry+dibdry*tipbry &
                       +depbry*tebdry+debdry*tepbry) +pbimpb
      sigpreb=(dibdry**2*stibdry**2+sdibdry**2*tibdry**2 &
              +debdry**2*stebdry**2+sdebdry**2*tebdry**2)
      sigpreb=1.602e+03*sqrt(sigpreb)
      sigppb = dibdry**2*sigtipb**2+sdibdry**2*tipbry**2 &
              +dipbry**2*stibdry**2+sigdipb**2*tibdry**2 &
              +debdry**2*sigtepb**2+sdebdry**2*tepbry**2 &
              +depbry**2*stebdry**2+sigdepb**2*tebdry**2
      sigppb =1.602e+03*sqrt(sigppb)
      do 2500 i=1,npress
        pressr(i)=1.602e+03*(dnitho(i)*tithom(i)+dnethom(i)*tethom(i)) &
                  +pbimth(i)
        sigpre(i)=(snitho(i)**2*tithom(i)**2+dnitho(i)**2*stitho(i)**2 &
                +sgneth(i)**2*tethom(i)**2+dnethom(i)**2*sgteth(i)**2)
        sigpre(i)=1.602e+03*sqrt(sigpre(i))
 2500 continue
      sgggmin=sgprmin
      if (sgprmin.lt.0.0) sgggmin=abs(sgprmin)*pressr(1)
      do 2600 i=1,npress
        sigpre(i)=max(sigpre(i),sgggmin)
 2600 continue
      if (kpressb.eq.1) then
      npress=npress+1
      pressr(npress)=pressbi
      sigpre(npress)=sigprebi
      rpress(npress)=-1.0
      endif
      if (ndokin.ge.100) then
        call getfnmu(itimeu,'k',ishot,itime,edatname)
        edatname='edat_'//edatname(2:7)// &
                       '_'//edatname(9:13)//'.pressure'
        open(unit=nin,status='old',file=edatname,err=12916)
        close(unit=nin,status='delete')
12916 continue
        open(unit=nin,status='new',file=edatname &
                                 )
        write (nin,edat)
        close(unit=nin)
      endif
      return
      end
      subroutine prtout(it)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          prtout performs printing.                               **
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
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: worka2
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      dimension xrsp(npcurn)
      dimension patmpz(magpri),xmpz(magpri),ympz(magpri),ampz(magpri)
      common/jwork4/workb(nsilop)
      character*30 sfname
      namelist/in3/mpnam2,xmp2,ymp2,amp2,smp2,rsi,zsi,wsi,hsi,as, &
        as2,lpname,rsisvs,vsname,turnfc,patmp2,racoil,zacoil, &
        wacoil,hacoil
!
      if (itek.gt.0) go to 100
!vas      write (nttyo,10000)
! MPI >>>
      if (rank == 0) then
! MPI <<<
        write (nttyo,10000) trim(ch1),trim(ch2) 
        jtime=time(it)
        saisq=tsaisq(it)
        write (nttyo,9300)
        write (nttyo,9320) ipsi(it)
        write (nttyo,9340) imag2(it)
        write (nttyo,9380) iplasm(it)
        write (nttyo,9385) idlopc(it)
        write (nttyo,10480)
        write (nttyo,10500) ishot,jtime,saisq
        write (nttyo,10520) betat(it),betap(it),ali(it)
        write (nttyo,10540) vout(it),rout(it),zout(it)
        write (nttyo,10560) eout(it),doutu(it),doutl(it)
        write (nttyo,10580) aout(it),oleft(it),oright(it)
        write (nttyo,10600) otop(it),qsta(it),rcurrt(it)
        write (nttyo,10610) zcurrt(it),bcentr(it),qout(it)
! MPI >>>
      endif
! MPI <<<
  100 continue
!
! --- delete fitout.dat if IOUT does not contain 1
!
      if (iand(iout,1).eq.0) then
         close(unit=nout,status='delete',err=300)
      else                      ! goto 300
!vas      write (nout,10000) 
      write (nout,10000) trim(ch1),trim(ch2) 
!
        jtime=time(it)
        saisq=tsaisq(it)
        write (nout,9300)
        write (nout,9320) ipsi(it)
        write (nout,9340) imag2(it)
        write (nout,9380) iplasm(it)
        write (nout,9385) idlopc(it)
        write (nout,10480)
        write (nout,10500) ishot,jtime,saisq
        write (nout,10520) betat(it),betap(it),ali(it)
        write (nout,10540) vout(it),rout(it),zout(it)
        write (nout,10560) eout(it),doutu(it),doutl(it)
        write (nout,10580) aout(it),oleft(it),oright(it)
        write (nout,10600) otop(it),qsta(it),rcurrt(it)
        write (nout,10610) zcurrt(it),bcentr(it),qout(it)
        write (nout,10620) olefs(it),orighs(it),otops(it)
        write (nout,10623) betat2
!
      if (icalbet.gt.0) &
        write (nout,10622) betat(it),vbtvac,vbtot2,vbtvac2,vbtor2,vbeta0
      if (icalbet.gt.0) &
        write (nout,10624) vbtmag,btvvac2,btvtor2,btvtot2
!
      if (scalea) then
        rowcnd=1./rowcnd
        colcnd=1./colcnd
        write (nout,11001)
        write (nout,11004) infosc,rowcnd,colcnd,arspmax
      endif
!
        write (nout,11030)
        write (nout,11020) (brsp(i)/turnfc(i),i=1,nfcoil)
        write (nout,11000)
        write (nout,11020) (brsp(i),i=1,nfcoil)
        sumif=0.0
        do 200 i=1,5
  200   sumif=sumif+brsp(i)+brsp(i+9)
        sumif=sumif+brsp(8)+brsp(17)
      sumift=0.0
      sumifs=0.0
      do 220 i=1,nfcoil
        sumift=sumift+brsp(i)
        sumifs=sumifs+brsp(i)**2
  220 continue
      sumifs=sqrt(sumifs/float(nfcoil-1))
!
        write (nout,10020) sumif,sumift,sumifs
        write (nout,11010)
        write (nout,11020) (rsisfc(i),i=1,nfcoil)
        if (ivacum.gt.0) go to 228
        if (icurrt.ne.2.and.icurrt.ne.5) go to 228
        xnorm=brsp(nfcoil+1)
        write (nout,11040) xnorm
        do 225 i=1,kwcurn
          xrsp(i)=brsp(nfcoil+i)/xnorm
  225   continue
        write (nout,11020) (xrsp(i),i=1,kwcurn)
        xnorm=darea
        write (nout,11043) xnorm
        do 226 i=1,kwcurn
          xrsp(i)=brsp(nfcoil+i)/xnorm
  226   continue
        write (nout,11020) (xrsp(i),i=1,kwcurn)
        if (keecur.gt.0) then
          write (nout,11032)
          write (nout,11020) (cerer(i),i=1,keecur)
        endif
        if (kedgep.gt.0) then
          write (nout,11034) pedge
        endif
        if (kedgef.gt.0) then
          write (nout,11036) f2edge
        endif
  228   continue
!
        write (nout,11005)
        write (nout,11020) (cecurr(i),i=1,nesum)
        write (nout,11008)
        write (nout,11020) (rsisec(i),i=1,nesum)
!
        if (ivesel.le.0) go to 300
        sumif=0.0
        do 280 i=1,nvesel
          sumif=sumif+vcurrt(i)
  280   continue
        nves2=nvesel/2
        sumift=0.0
        do 282 i=1,nves2
          sumift=sumift+vcurrt(i)
  282   continue
        sumifb=0.0
        do 284 i=nves2+1,nvesel
          sumifb=sumifb+vcurrt(i)
  284   continue
        if (ivesel.eq.11) sumif=sumvs0
        write (nout,11060) sumif,sumift,sumifb
        write (nout,11020) (vcurrt(i),i=1,nvesel)
        write (nout,11022) fzpol
        write (nout,11020) (vforcep(i),i=1,nvesel)
        write (nout,11024) fztor
        write (nout,11020) (vforcet(i),i=1,nvesel)
      endif
  300   continue
!
      if (patmp2(1).gt.0.0) go to 340
      open(unit=80,status='old', &
           file=table_di2(1:ltbdi2)//'dprobe.dat' &
                                          )
      read (80,in3)
      close(unit=80)
      if (xmp2(1).le.0.0) go to 360
      if (patmp2(1).gt.0.0) go to 340
      xmin=xmp2(1)
      xmax=xmin
      ymin=ymp2(1)
      ymax=ymin
      do 305 i=1,magpri67
        xmpz(i)=xmp2(i)
        ympz(i)=ymp2(i)
        xmin=min(xmin,xmp2(i))
        xmax=max(xmax,xmp2(i))
        ymin=min(ymin,ymp2(i))
        ymax=max(ymax,ymp2(i))
  305 continue
      xtest=(xmin+xmax)/2.
      ytest=(ymin+ymax)/2.
      nzz=0
      call packps(xmpz,ympz,magpri67,xtest,ytest,nzz)
      do 310 i=1,magpri67
        do 310 k=1,magpri67
          if ((xmp2(k).eq.xmpz(i)).and.(ymp2(k).eq.ympz(i))) &
           ampz(i)=amp2(k)
  310 continue
      do 320 i=1,magpri67
        ip1=i+1
        im1=i-1
        if (i.eq.1) im1=magpri67
        if (i.eq.magpri67) ip1=1
        dang90=abs(abs(ampz(i))-90.)
        dang270=abs(abs(ampz(i))-270.)
        danm90=abs(abs(ampz(im1))-90.)
        danm270=abs(abs(ampz(im1))-270.)
        danp90=abs(abs(ampz(ip1))-90.)
        danp270=abs(abs(ampz(ip1))-270.)
        if (dang90.lt.1.0.or.dang270.lt.1.0) then
         xxm=xmpz(i)
         xxp=xmpz(i)
         if (danm90.lt.1.0.or.danm270.lt.1.0) then
          yym=(ympz(i)+ympz(im1))/2.
         else
          sm1=tand(ampz(im1))
          yym=sm1*(xxm-xmpz(im1))+ympz(im1)
         endif
         if (danp90.lt.1.0.or.danp270.lt.1.0) then
          yyp=(ympz(i)+ympz(ip1))/2.
         else
          sm2=tand(ampz(ip1))
          yym=sm2*(xxm-xmpz(ip1))+ympz(ip1)
         endif
        else
         if (danm90.lt.1.0.or.danm270.lt.1.0) then
          xxm=xmpz(im1)
          sm2=tand(ampz(i))
          yym=sm2*(xxm-xmpz(i))+ympz(i)
         else
          dampz1=abs(ampz(im1)-ampz(i))
          dampz2=abs(dampz1-360.)
          if (dampz1.lt.1.0.or.dampz2.lt.1.0) then
           xxm=(xmpz(i)+xmpz(im1))/2.
           yym=(ympz(i)+ympz(im1))/2.
          else
           sm1=tand(ampz(im1))
           sm2=tand(ampz(i))
           xxm=(sm1*xmpz(im1)-sm2*xmpz(i)-ympz(im1)+ympz(i))/(sm1-sm2)
           yym=sm1*(xxm-xmpz(im1))+ympz(im1)
          endif
         endif
         if (danp90.lt.1.0.or.danp270.lt.1.0) then
          xxp=xmpz(ip1)
          sm1=tand(ampz(i))
          yyp=sm1*(xxp-xmpz(i))+ympz(i)
         else
          dampz1=abs(ampz(ip1)-ampz(i))
          dampz2=abs(dampz1-360.)
          if (dampz1.lt.1.0.or.dampz2.lt.1.0) then
           xxp=(xmpz(i)+xmpz(ip1))/2.
           yyp=(ympz(i)+ympz(ip1))/2.
          else
           sm1=tand(ampz(i))
           sm2=tand(ampz(ip1))
           xxp=(sm1*xmpz(i)-sm2*xmpz(ip1)-ympz(i)+ympz(ip1))/(sm1-sm2)
           yyp=sm1*(xxp-xmpz(i))+ympz(i)
          endif
         endif
        endif
        patmpz(i)=sqrt((xxp-xxm)**2+(yyp-yym)**2)
  320 continue
      do 330 i=1,magpri67
        do 330 k=1,magpri67
          if ((xmpz(k).eq.xmp2(i)).and.(ympz(k).eq.ymp2(i))) &
           patmp2(i)=patmpz(k)
  330 continue
  340 continue
      cipmp2=0.0
      do 350 i=1,magpri67
        cipmp2=cipmp2+cmpr2(i,it)*patmp2(i)
  350 continue
      cipmp2=cipmp2/tmu/twopi
  360 continue
!-----------------------------------------------------------------
!--   322 degree probes                                         --
!-----------------------------------------------------------------
      mb=magpri67+1
      mbb=magpri322
      if (xmp2(mb).le.0.0) go to 22360
      if (patmp2(mb).gt.0.0) go to 22340
      xmin=xmp2(mb)
      xmax=xmin
      ymin=ymp2(mb)
      ymax=ymin
      do 22305 i=mb,magpri67+magpri322
        xmpz(i)=xmp2(i)
        ympz(i)=ymp2(i)
        xmin=min(xmin,xmp2(i))
        xmax=max(xmax,xmp2(i))
        ymin=min(ymin,ymp2(i))
        ymax=max(ymax,ymp2(i))
22305 continue
      xtest=(xmin+xmax)/2.
      ytest=(ymin+ymax)/2.
      nzz=0
      call packps(xmpz(mb),ympz(mb),mbb,xtest,ytest,nzz)
      do 22310 i=mb,magpri67+magpri322
        do 22310 k=mb,magpri67+magpri322
          if ((xmp2(k).eq.xmpz(i)).and.(ymp2(k).eq.ympz(i))) &
           ampz(i)=amp2(k)
22310 continue
      do 22320 i=mb,magpri67+magpri322
        ip1=i+1
        im1=i-1
        if (i.eq.mb) im1=magpri67+magpri322
        if (i.eq.magpri67+magpri322) ip1=mb
        dang90=abs(abs(ampz(i))-90.)
        dang270=abs(abs(ampz(i))-270.)
        danm90=abs(abs(ampz(im1))-90.)
        danm270=abs(abs(ampz(im1))-270.)
        danp90=abs(abs(ampz(ip1))-90.)
        danp270=abs(abs(ampz(ip1))-270.)
        if (dang90.lt.1.0.or.dang270.lt.1.0) then
         xxm=xmpz(i)
         xxp=xmpz(i)
         if (danm90.lt.1.0.or.danm270.lt.1.0) then
          yym=(ympz(i)+ympz(im1))/2.
         else
          sm1=tand(ampz(im1))
          yym=sm1*(xxm-xmpz(im1))+ympz(im1)
         endif
         if (danp90.lt.1.0.or.danp270.lt.1.0) then
          yyp=(ympz(i)+ympz(ip1))/2.
         else
          sm2=tand(ampz(ip1))
          yym=sm2*(xxm-xmpz(ip1))+ympz(ip1)
         endif
        else
         if (danm90.lt.1.0.or.danm270.lt.1.0) then
          xxm=xmpz(im1)
          sm2=tand(ampz(i))
          yym=sm2*(xxm-xmpz(i))+ympz(i)
         else
          dampz1=abs(ampz(im1)-ampz(i))
          dampz2=abs(dampz1-360.)
          if (dampz1.lt.1.0.or.dampz2.lt.1.0) then
           xxm=(xmpz(i)+xmpz(im1))/2.
           yym=(ympz(i)+ympz(im1))/2.
          else
           sm1=tand(ampz(im1))
           sm2=tand(ampz(i))
           xxm=(sm1*xmpz(im1)-sm2*xmpz(i)-ympz(im1)+ympz(i))/(sm1-sm2)
           yym=sm1*(xxm-xmpz(im1))+ympz(im1)
          endif
         endif
         if (danp90.lt.1.0.or.danp270.lt.1.0) then
          xxp=xmpz(ip1)
          sm1=tand(ampz(i))
          yyp=sm1*(xxp-xmpz(i))+ympz(i)
         else
          dampz1=abs(ampz(ip1)-ampz(i))
          dampz2=abs(dampz1-360.)
          if (dampz1.lt.1.0.or.dampz2.lt.1.0) then
           xxp=(xmpz(i)+xmpz(ip1))/2.
           yyp=(ympz(i)+ympz(ip1))/2.
          else
           sm1=tand(ampz(i))
           sm2=tand(ampz(ip1))
           xxp=(sm1*xmpz(i)-sm2*xmpz(ip1)-ympz(i)+ympz(ip1))/(sm1-sm2)
           yyp=sm1*(xxp-xmpz(i))+ympz(i)
          endif
         endif
        endif
        patmpz(i)=sqrt((xxp-xxm)**2+(yyp-yym)**2)
22320 continue
      do 22330 i=mb,magpri67+magpri322
        do 22330 k=mb,magpri67+magpri322
          if ((xmpz(k).eq.xmp2(i)).and.(ympz(k).eq.ymp2(i))) &
           patmp2(i)=patmpz(k)
22330 continue
!
      open(unit=80,status='old',file='dprobe.new',err=12914)
      close(unit=80,status='delete')
12914 continue
      open(unit=80,status='new',file='dprobe.new' &
                                 )
      write (80,in3)
      close(unit=80)
22340 continue
!
      cipmp3=0.0
      do 22350 i=magpri67+1,magpri67+magpri322
        cipmp3=cipmp3+cmpr2(i,it)*patmp2(i)
22350 continue
      cipmp3=cipmp3/tmu/twopi
22360 continue
!
      if (.not.fitsiref) then
      ssiref=csilop(iabs(nslref),it)
      do 363 i=1,nsilop
        workb(i)=csilop(i,it)-ssiref
  363 continue
      else
      ssiref=0.0
      do i=1,nsilop
        workb(i)=csilop(i,it)+csiref
      enddo
      endif
!
      if (iand(iout,1).ne.0) then  ! goto 850
!
        write (nout,11100) ssiref
        write (nout,11020) (workb(i),i=1,nsilop)
        write (nout,11100) csiref
        write (nout,11020) (csilop(i,it),i=1,nsilop)
        workb(1:nsilop)=csilop(1:nsilop,it)*twopi
        write (nout,11101) csiref*twopi
        write (nout,11020) (workb(i),i=1,nsilop)
        write (nout,11102)
        write (nout,11020) (csilopv(i,it),i=1,nsilop)
        write (nout,11120) cipmp2,cipmp3
        write (nout,11020) (cmpr2(i,it),i=1,magpri)
        write (nout,11122)
        write (nout,11020) (cmpr2v(i,it),i=1,magpri)
        if (kstark.gt.0) then
          write (nout,11140)
          write (nout,11020) (cmgam(i,it),i=1,nstark)
        endif
        write (nout,11160) cpasma(it)
        write (nout,11180) cdflux(it),sbpp,delbp(it),sbppa
        if (kecebz.gt.0) write(nout,11185) cmecebz(it)
        if (kece.gt.0)  then
              write(nout,11186)
              write(nout,11020)  (cmece(m,it),m=1,nece)
        endif
!
        write (nout,11200) psiref(it)
        write (nout,11020) (silopt(it,i),i=1,nsilop)
        write (nout,11220)
        write (nout,11020) (expmpi(it,i),i=1,magpri)
        if (kstark.gt.0) then
          write (nout,11240)
          write (nout,11020) (tangam(it,i),i=1,nstark)
        endif
        write (nout,11260) pasmat(it)
        write (nout,11280) diamag(it)
        write (nout,11270) vloopt(it)
        write (nout,11292)
        write (nout,11020) (fccurt(it,i),i=1,nfcoil)
        write (nout,11294)
        write (nout,11020) (eccurt(it,i),i=1,nesum)
!
      if (abs(sigdia(it)).le.1.0e-08) go to 520
      write (nout,11300) chidlc
  520 continue
      if (iconvr.ne.3) go to 540
      write (nout,11320) emf,emp,enf,enp,betap0,rzero
      write (nout,11330) cbetap,cli,cqqxis,cbetat,ci0
  540 continue
      if (nbdry.le.0) go to 544
      write (nout,11324) erbmax,erbave
      write (nout,11326) (erbloc(i),i=1,nbdry)
  544 continue
!
      if (kvtor.gt.0) then
        write (nout,13000)
        write (nout,13020) betatw(it),betapw(it),wplasw(it)
      endif
!
      write (nout,12000)
      do 600 i=1,nitera
        write (nout,12020) i,cerror(i),csibry(i),csimag(i),cvolp(i), &
                  crmaxi(i),czmaxi(i),cemaxi(i),cqmaxi(i),cchisq(i)
  600 continue
      if ((kwripre.gt.0).and.(kwripre.le.9)) then
          call getfnmd('n',ishot,itime,sfname)
          sfname=sfname(1:13)//'_error'
          open(unit=74,status='old',file=sfname,err=12918)
          close(unit=74,status='delete')
12918     continue
          open(unit=74,status='new',file=sfname                       )
          do i=1,nitera
           write (74,*) i,cerror(i),xdum,xdum
          enddo
          close(unit=74)
      endif
      write (nout,12010)
      do 620 i=1,nitera
        write (nout,12025) i,aveerr(i),csumip(i),tratio(i),iermax(i), &
                  jermax(i)
  620 continue
!
      if ((kecebz.gt.0).or.(kece.gt.0)) then
        write (nout,12015)
       do i=1,nitera
        write (nout,12020) i,receoi(i),(recemi(i,k),k=1,nece)
       enddo
        write (nout,12016)
       do i=1,nitera
        write (nout,12020) i,(recepi(i,k),k=1,nece)
       enddo
      endif
!
      dsi=1./float(nw-1)
      write (nout,12040)
      do 700 i=1,nw
        sinow=dsi*(i-1)
        write (nout,12020) i,sinow,volp(i),pprime(i),curmid(i),ffprim(i) &
                ,pres(i),fpol(i),qpsi(i),rpres(i)
  700 continue
!
      write (nout,12043)
      do 800 i=1,nw
        sinow=dsi*(i-1)
        write (nout,12020) i,sinow,cjor(i)
  800 continue
!
      endif
  850 continue
      return
 9300 format (/,4x,16h   data used:   )
 9320 format (1x,i2,11h flux loops)
 9340 format (1x,i2,19h magnetic probes(i))
 9360 format (1x,i2,18h partial rogowskis)
 9380 format (1x,i2,14h full rogowski)
 9385 format (1x,i2,17h diamagnetic loop)
 9390 format (1x,12h bt0(t)   = ,f10.3)
!10000 format(/,6x,20('*'),' EFITD 129dx2 output ',20('*'))
10000 format(/,6x,20('*'),' EFITD',a3,' x ',a3,'  output ',20('*'))
10020 format (1x,15h  sumif(amp) = ,e10.3,15h sumift(amp) = ,e10.3, &
           15h sumifs(amp) = ,e10.3)
10480 format (1x,/)
10500 format(12h shot #   = ,i10,12h time(ms) = ,i10, &
             12h chi**2   = ,e10.3)
10520 format(12h betat(%) = ,f10.3,12h betap    = ,f10.3, &
             12h li       = ,f10.3)
10540 format(12h vol(cm3) = ,e10.3,12h rout(cm) = ,f10.3, &
             12h zout(cm) = ,f10.3)
10560 format(12h elong    = ,f10.3,12h utriang  = ,f10.3, &
             12h ltriang  = ,f10.3)
10580 format(12h a(cm)    = ,f10.3,12h lin(cm)  = ,f10.3, &
             12h lout(cm) = ,f10.3)
10600 format(12h ltop(cm) = ,f10.3,12h q*       = ,f10.3, &
             12h rc(cm)   = ,f10.3)
10610 format(12h zc(cm)   = ,f10.3,12h bt0(t)   = ,f10.3, &
             12h qout     = ,f10.3)
10620 format(12h lins(cm) = ,f10.3,12h louts(cm)= ,f10.3, &
             12h ltops(cm)= ,f10.3)
10623 format(12h beta*(%) = ,f10.3)
10622 format(//, &
      ' betat, betat-btvac, beta-total, beta-btvac2, beta-btv :',/, &
      ' betat0 :',/, &
             5(2x,e12.4,2x),/,1(2x,e12.4,2x))
10624 format(//, &
      ' betatm, betat-vbtvac2, beta-vbtor2, beta-vbtot2:',/, &
             4(2x,e12.4,2x))
11000 format(//,22x,'F-coils currents (Amp)')
11001 format(//,22x,'A matrix condition    ')
11002 format(//,22x,16hE-coils phases  )
11004 format(1x,' info = ',i4,' row = ',1pe10.3,' col = ',1pe10.3, &
                ' max = ',1pe10.3)
11005 format(//,22x,16hE-coils currents)
11008 format(//,22x,24hE-coils resistance(Ohms))
11010 format(//,22x,24hF-coils resistance(Ohms))
11017 format(//,2x,'power supply current (A) = ',e12.5, &
                5x,'phase (degree) = ',e12.5,/, &
                2x,'resistance (Ohm)         = ',e12.5, &
                5x,'inductance (H) = ',e12.5)
11020 format(4e15.6)
11022 format(/,12x,' vessel vertical forces (p, newton) ',e15.6)
11024 format(/,12x,' vessel vertical forces (t, newton) ',e15.6)
11030 format(//,22x,'F-coils currents (Amp/turn)')
11032 format(//,12x,'Electrostatic potential derivative PIEPRIM:')
11034 format(//,12x,'Hyperbolic P:',e15.6)
11036 format(//,12x,'Hyperbolic FF:',e15.6)
11040 format(//,22x,15hplasma currents,18h  normalization = ,e15.6)
11043 format(//,22x,15hplasma coeffics,18h  normalization = ,e15.6)
11060 format(//,22x,15hvessel currents,12h sum(amp) = ,e10.3, &
                12h  top   =   ,e10.3,12h  bot     = ,e10.3)
11100 format(//,16x,28hcalculated psi-loops signals, &
        18h ssiref(vs/rad) = ,e12.5)
11101 format(//,16x,28hcalculated psi-loops signals, &
        14h ssiref(vs) = ,e12.5)
11102 format(//,16x,35hcalculated vacuum psi-loops signals)
11120 format(//,4x,37hcalculated magnetic probes(i) signals, &
        12h ipmp2(A) = ,e12.5,1x,e12.5)
11122 format(//,4x,44hcalculated vacuum magnetic probes(i) signals)
11140 format(//,14x,31hcalculated polarimetry signals )
11160 format(//,14x,32hcalculated total plasma current ,/,16x,e15.6)
11180 format(//,14x,32hcalculated diamagnetic flux(vs) ,/,16x,e15.6,/, &
           '     bpdia = ',e10.3,'     delbp = ',e10.3, &
           ' app bpdia = ',e10.3)
11185 format(//,14x,32hcalculated Bz(receo,zeceo) (T)  ,/,16x,e15.6)
11186 format(//,22x,35hcalculated psi(R-)-psi(R+) (VS/rad))
11200 format(//,16x,28h  measured psi-loops signals, &
        18h psiref(vs/rad) = ,e12.5)
11220 format(//,12x,37h  measured magnetic probes(i) signals)
11240 format(//,14x,31h  measured polarimetry signals )
11260 format(//,14x,32h  measured total plasma current ,/,16x,e15.6)
11270 format(//,14x,32h  measured loop voltage (V)     ,/,16x,e15.6)
11280 format(//,14x,32h  measured diamagnetic flux     ,/,16x,e15.6)
11292 format(//,14x,32h  measured F-coil currents(A)   )
11294 format(//,14x,32h  measured E-coil currents(A)   )
11300 format(//,14x,32h    chisqr diamagnetic flux  =  ,e10.3)
11320 format (//,12h   emf  =   ,e10.3,12h   emp  =   ,e10.3, &
           12h   enf  =   ,e10.3,12h   enp  =   ,e10.3,/, &
           12h betap0 =   ,e10.3,12h rzero  =   ,e10.3)
11324 format (//,12h  erbmax =  ,e12.5,12h   erbave = ,e12.5)
11326 format (8(1x,e12.5))
11330 format (//,12h  cbetap =  ,e12.5,12h   cli    = ,e12.5, &
                 12h  cqqxis =  ,e12.5,12h   cbetat = ,e12.5,/, &
                 12h  ci0    =  ,e12.5)
12000 format (//,19h iteration summary:,/,4h  i ,2x,12h   error    ,2x, &
        12h   psibry   ,2x,12h   psimag   ,2x,12h   volume   ,2x, &
        12h   rmaxis   ,2x,12h   zmaxis   ,12h   emaxis   ,2x, &
        12h   qmaxis   ,2x,12h   chisqr   )
12010 format (//,4h  i ,2x,12h   errave   ,2x, &
        12h   current  ,2x,12h   cratio   ,2x,12h   iermax   ,2x, &
        12h   jermax   ,2x,12h            ,12h            ,2x, &
        12h            ,2x,12h            )
12015 format (//,23h ECE iteration summary:,/,4h  i ,2x,12h   receo    ,2x, &
        12h   recem    )
12016 format (//,4h  i ,2x,12h   recep    )
12020 format (i4,10(1x,e12.5))
12025 format (i4,3(2x,e12.5),2(6x,i4,4x),4(2x,e12.5))
12040 format (//,19h    plasma summary:,/,4h  i ,1x,12h   pflux    ,1x, &
        12h   vol(m3)  ,1x,12h   pprime   ,1x,12h  current   ,1x, &
        12h ffprime    ,1x,12h  pressure  ,12h   fpol     ,1x, &
        12h    q       ,1x,12h     rm     )
12043 format (//,4h  i ,1x,12h    pflux   ,1x, &
        12h   <j>      ,1x,12h            ,1x,12h            ,1x, &
        12h            ,1x,12h            ,12h            ,1x, &
        12h            ,1x,12h            )
13000 format (//,16x,'    toroidal rotation         ')
13020 format (12h  betatw =  ,1pe12.5,12h betapw   = ,1pe12.5, &
               12h  W(J)   =  ,1pe12.5)
      end
      function pwcur4(n1set,ypsi,nnn)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          pwcur4 computes the rotational pressure by integration. **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          94/02/28..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: sifprw,bwprw,cwprw,dwprw,sfprw,sprwp
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!
      if (abs(ypsi).gt.1.0) then
        pwcur4=0.0
        return
      endif
!
      if (n1set.gt.0) then
       do i=2,nw-1
        siii=float(i-1)/float(nw-1)
        sifprw(i)=siii
       enddo
       sifprw(1)=0.0
       sifprw(nw)=1.0
       do i=1,nw
        sprwp(i)=pwpcu4(sifprw(i),kwwcur)/darea
       enddo
!
       sfprw(nw)=preswb
       delsi=sidif/float(nw-1)
       do 1000 i=1,nw-1
       sfprw(nw-i)=sfprw(nw-i+1)+0.5*(sprwp(nw-i+1)+sprwp(nw-i))*delsi
 1000  continue
!
       mw=nw
       call zpline(mw,sifprw,sfprw,bwprw,cwprw,dwprw)
      endif
      pwcur4=seval(mw,ypsi,sifprw,sfprw,bwprw,cwprw,dwprw)
      return
      end
      function pwpcu4(ypsi,nnn)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          pwpcu4 computes the radial derivative of the            **
!**          rotational pressure based on icurrt.                    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          94/03/08..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!
      if (abs(ypsi).gt.1.0) then
        pwpcu4=0.0
        return
      endif
      if (icurrt.eq.4) goto 2000
      if (icurrt.eq.1) goto 3000
      pwpcu4=pwpcur(ypsi,nnn)
      return
!
 2000 continue
      pwpcu4=((1.-ypsi**enw)**emw*(1.-gammaw)+gammaw)*rbetaw &
         *cratio/rzero
      return
!
 3000 continue
      pwpcu4=sbetaw*cratio/srma
      return
      end
      function pwpcur(ypsi,nnn)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          ppcurr computes the radial derivative of the            **
!**          rotational pressure.                                    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          94/03/08..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
      dimension xpsii(nwwcur)
!
      if (abs(ypsi).gt.1.0) then
        pwpcur=0.0
        return
      endif
      pwpcur=0.0
      call setpwp(ypsi,xpsii)
      do 1400 iiij=nfnpcr+1,nnn+nfnpcr
        iijj=iiij-nfnpcr
        pwpcur=pwpcur+brsp(iiij)*xpsii(iijj)
 1400 continue
      return
!
      entry pwcurr(ypsi,nnn)
      if (abs(ypsi).gt.1.0) then
        pwcurr=0.0
        return
      endif
      pwcurr=0.0
      call setpw(ypsi,xpsii)
      do 1600 i=nfnpcr+1,nfnpcr+nnn
        nn=i-nfnpcr
        pwcurr=pwcurr+brsp(i)*xpsii(nn)
 1600 continue
      pwcurr=-sidif*pwcurr/darea+preswb
      return
      end
      subroutine residu(nx,jtime)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**         residu computes the flux variations on the r-z grid.     **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**          28/03/93..........fixe relax                            **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: psiold,psipold,psipp
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!
      if (ivacum.gt.0) return
      errold=errorm
      errave=0.0
      errorm=0.0
      do 1000 i=1,nw
      do 1000 j=1,nh
        kk=(i-1)*nh+j
        change=abs(psi(kk)-psiold(kk))
        errorm=max(errorm,change)
        errave=errave+change
        if (errorm.gt.change) go to 1000
        iermax(nx)=i
        jermax(nx)=j
 1000 continue
      errorm=errorm/abs(sidif)
!
      aveerr(nx)=errave/abs(sidif)/float(nwnh)
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
          if (nx.eq.1) write (nttyo,10017) itime, rank
          if (nsol.eq.0) then
           write (nttyo,10019)  nx,tsaisq(jtime),zmaxis,errorm,delzmm, &
             sum(chigam)
          else
           write (nttyo,10020)  nx,tsaisq(jtime),zmaxis,errorm,delzmm, &
             sum(chigam),erbmax,erbsmax           
          endif
        elseif (itell.eq.2) then
          write (nttyo,10021)  nx,ali(jtime),abs(betatn),errorm,qsiw(1)
        elseif (itell.eq.3) then
          write (nttyo,10023)  nx,difpsi,zmaxis,errorm,delzmm
        elseif (itell.eq.4) then
          if (nx.eq.1) then
            write (nttyo,10017) itime, rank
            write (nttyo,10019)  nx,tsaisq(jtime),zmaxis,errorm,delzmm, &
        sum(chigam)
          else
            write (nttyo,10025)  nx,tsaisq(jtime),zmaxis,errorm &
                                 ,delzmm,cdelz(nx-1),cdeljsum
          endif
        endif
      endif
      if (isetfb.ne.0) then
           write (4,10009) nx,tsaisq(jtime),zmaxis,errorm,delzmm &
                           ,brfb(1)
           if (isetfb.lt.0) &
           write (6,10009) nx,tsaisq(jtime),zmaxis,errorm,delzmm &
                           ,brfb(1)
        elseif (eelip.gt.2.25.and.itell.eq.0) then
           write (6,10009) nx,tsaisq(jtime),zmaxis,errorm,delzmm &
                           ,brfb(1)
      endif
      if (idebug /= 0) then
        write (nttyo,*) 'cratio,cratio_ext,cratiop_ext,cratiof_ext= ', &
           cratio,cratio_ext,cratiop_ext,cratiof_ext
        write (nttyo,*) 'scalepp_ext,scaleffp_ext= ', &
           scalepp_ext,scaleffp_ext
      endif
      return
10009   format (x,'iter',i3.3, &
        ' chsq=',1pe8.2,' zmag=',1pe9.2,' err=',1pe8.2,' dz=',1pe10.3, &
        ' Ifb=',1pe9.2)
10017   format (/,x,' ----- time =',i6,' ms ----- (',i2,')')
10019   format (x,'it=',i3, &
        ' chi2=',1pe8.2,' zm=',1pe9.2,' err=',1pe9.3, &
        ' dz=',1pe10.3,' chigam=',1pe9.2)
10020   format (x,'it=',i3, &
        ' chi2=',1pe8.2,' zm=',1pe9.2,' err=',1pe9.3, &
        ' dz=',1pe10.3,' chigam=',1pe9.2,' errb=',1pe9.2,' errbs=',1pe9.2)
10021   format (x,'it=',i3, &
        ' li=',1pe9.3,' betan=',1pe9.3,' err=',1pe9.3, &
        ' qs=',1pe9.3)
10023   format (x,'it=',i3, &
        ' dpsi=',1pe10.3,' zm=',1pe9.2,' err=',1pe9.3, &
        ' dz=',1pe10.3)
10025   format (x,'it=',i3, &
        ' chi2=',1pe8.2,' zm=',1pe9.2,' err=',1pe9.3, &
        ' dz=',1pe10.3,' delz=',1pe10.3,' dj=',1pe9.3)
      end
      subroutine setece(jtime,kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          setece obtains the response matrix                      **
!**          for ECE measurement                                     **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          11/98..........first created, Cheng Zhang               **
!**     2013/08/07..........updated for real-space Ti                **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
!
      common/cwork3/lkx,lky
      dimension rrrecep(nnece),rrrecem(nnece),iwp(nnece),iwm(nnece)
      dimension zece(nnece),pds(6),rsplt(2500),zsplt(2500),csplt(2500)
      real*8,dimension(:),allocatable :: dsidr,ddsiddr
      character*40 filenmme
      data nset/20/
      save nset
      integer, intent(inout) :: kerror
      kerror = 0
!
      if (idebug.ge.3) write (6,*) 'Enter SETECE, nw = ',nw
      ALLOCATE(dsidr(nw),ddsiddr(nw))
!
      if (idebug.ge.3) write (6,*) 'Enter SETECE, ksetece = ',ksetece
      kece=0
      kecebz=0
!
      if (ksetece.ne.0) go to 60000
!------------------------------------------------------------------
!--     zeceo from k-file, jo--zeceo on zgrid
!------------------------------------------------------------------
      zeceo=zteo
      receo=rteo
      dho=abs(zgrid(1)-zeceo)
      jo=1
      do jh=2,nh
         if (abs(zgrid(jh)-zeceo).lt.dho) then
          dho=abs(zgrid(jh)-zeceo)
          jo=jh
         endif
      enddo
      if (idebug.ge.3) write (6,*) 'SETECE, receo/zeteo = ',receo,zeceo
!-------------------------------------------------------------------
!--     kfixro=1, receo  from k-file
!------------------------------------------------------------------- 
      if (kfixro.eq.1) then
         receo=rteo
         if (kfitece.eq.1) go to 20
      endif
!------------------------------------------------------------------
!--     kfixrece=1, R+(recep) R-(recem) from k-file
!------------------------------------------------------------------
      if (kfixrece.eq.1) then
        do k=1,nnece
         if ((rtep(k).gt.1.E-5).and.(rtem(k).gt.1.E-5)) then
           recep(k)=rtep(k)
           recem(k)=rtem(k)
           nece=k
         endif
        enddo
      endif
!-----------------------------------------------------------------
!--   kfixro=0 or kfixrece=0, receo or R+ R- from getecer
!-----------------------------------------------------------------
      if (kfitece.le.2) then
!      if ((kfixro.eq.0).or.(kfixrece.eq.0)) call getecer(jtime,kerrora)
       if (kfixrece.eq.0) call getecer(jtime,kerrora)
       if ((kfixro.eq.-1).or.(kfixrece.eq.-1)) call geteceb(jtime,kerror)
      endif
      if (kfitece.eq.3) then
        call gettir(jtime,kerror)
      endif
20    continue
!----------------------------------------------------------------
!--   get iwo iwp iwm (receo, R+ R- on rgrid)
!----------------------------------------------------------------
      ddo=abs(rgrid(1)-receo)
      iwo=1
      do iw=2,nw
        if (abs(rgrid(iw)-receo).lt.ddo) then
          ddo=abs(rgrid(iw)-receo)
          iwo=iw
        endif
      enddo
      if (kfitece.eq.1) go to 29
      do 26 k=1,nece
        zece(k)=zeceo
   26 continue
      do k=1,nece
        ddp=abs(rgrid(1)-recep(k))
        ddm=abs(rgrid(1)-recem(k))
        iwp(k)=1
        iwm(k)=1
        do iw=2,nw
         if (abs(rgrid(iw)-recep(k)).lt.ddp) then
          ddp=abs(rgrid(iw)-recep(k))
          iwp(k)=iw
         endif
        enddo
        do iw=2,nw
         if (abs(rgrid(iw)-recem(k)).lt.ddm) then
          ddm=abs(rgrid(iw)-recem(k))
          iwm(k)=iw
         endif
        enddo
      enddo  
   29 continue
!-------------------------------------------------------------------
!-- try read responce function from recefile
!--------------------------------------------------------------------
      do 20000 iname=1,nset
        if (iname.le.9) then
          write(filenmme,6530) iname
        else
          write(filenmme,6540) iname
        endif
        open(unit=nffile, &
             status='old',form='unformatted', &
             file=filenmme,err=10000)
        goto 15000
10000   continue
        open(unit=nffile, &
             status='old',form='unformatted', &
             file=table_dir(1:ltbdir)//filenmme,err=30)
15000   continue
        read (nffile,err=30) nnnece
        if (nnnece.ne.nece) then
          close(unit=nffile)
          go to 20000
        endif
        read (nffile,err=30)  rrrecem
        read (nffile,err=30)  rrrecep
        read (nffile,err=30)  rrreceo
        do ii=1,nece
         if (abs(rrrecem(ii)-recem(ii)).gt.1.e-4) go to 20000
         if (abs(rrrecep(ii)-recep(ii)).gt.1.e-4) go to 20000
        enddo
        if (abs(rrreceo-receo).gt.1.e-4) go to 20000
        read (nffile,err=30) recebzfc
        read (nffile,err=30) gecebzpc
        read (nffile,err=30) recebzec
        read (nffile,err=30) recefc
        read (nffile,err=30) gecepc
        read (nffile,err=30) receec
        close(unit=nffile)
        go to 25000
20000 continue
      go to 30
25000 return
!-------------------------------------------------------------------       
!-- compute the response function about Te peak point constraint    
!--  response due to F coils --- recebzfc(ncoil)
!-------------------------------------------------------------------
30    continue
      if (idebug.ge.3) write (6,*) 'SETECE, rf/zf = ',rf(1),zf(1)
      isplit=10
      itot=isplit*isplit
      fitot=itot 
      recebzfc = 0.0
      do 60 k=1,mfcoil
      bzct=0 
      call splitc(isplit,rsplt,zsplt,csplt, &
                  rf(k),zf(k),wf(k),hf(k),af(k),af2(k),cdum)
      r1=receo
      do 50 l=1,itot
        a=rsplt(l)  
        z1=zeceo-zsplt(l)
        bzct=bzct+bz(a,r1,z1)      
   50 continue   
      kkm=fcid(k)
      recebzfc(kkm)=recebzfc(kkm)+bzct/fitot*tmu
   60 continue
!---------------------------------------------------------------------
!--   plasma response   gecebzpc(nwnh)                              --
!---------------------------------------------------------------------
      zzx=0.0
      drzg=(drgrid**2+dzgrid**2)*25.
      isplit=18
      if (receo.le.1.e-8)  go to 90
      r=receo
      do 80 ii=1,nw
        a=rgrid(ii)
        do 80 jj=1,nh
          z=zeceo-zgrid(jj)
          kk=(ii-1)*nh+jj
          dist=(a-r)**2+z**2
          if (dist.lt.drzg) then
   68       call splitc(isplit,rsplt,zsplt,csplt, &
                  rgrid(ii),zgrid(jj),drgrid,dzgrid,zzx,zzx,cdum)
            gbzt=0.0
            itot=isplit*isplit
            do 72 k=1,itot
              a=rsplt(k)
              z=zeceo-zsplt(k)
              distt=(a-r)**2+z**2
              if (distt.lt.1.e-8.and.isplit.lt.49) then
                isplit=isplit+2
                go to 68
              endif
              gbzt=gbzt+bz(a,r,z)
   72        continue    
             gecebzpc(kk)=gbzt*tmu/float(itot)
          else
             gecebzpc(kk)=bz(a,r,z)*tmu
          endif
   80 continue
   90 continue
!---------------------------------------------------------------------
!-- E coils response   recebzec                                     --
!---------------------------------------------------------------------
      if (iecurr.le.0) go to 201 
      isplit=10
      itot=isplit*isplit
      fitot=float(itot)
      do 170 i=1,nesum
         recebzec(i)=0.0
  170 continue
      if (receo.le.1.e-8) goto 201
      r1=receo
      do 200 k=1,necoil
        bzct=0.0
        call splitc(isplit,rsplt,zsplt,csplt, &
                    re(k),ze(k),we(k),he(k),zzx,zzx,cdum)
        do 180 l=1,itot
          a=rsplt(l)
          z1=zeceo-zsplt(l)      
          bzct=bzct+bz(a,r1,z1)
  180 continue  
      bzct=bzct*tmu/fitot
      kkm=ecid(k)  
      recebzec(kkm)=recebzec(kkm)+bzct
  200 continue
  201 continue
!-------------------------------------------------------------------------
!--  compute the response function about R+ R- constraint     
!--  F coils response -----recefc(nece,nfcoil)
!-------------------------------------------------------------------------
      if (kfitece.eq.1) go to 1710
      do 1500 n=1,nfcoil
        call sets2d(gridfc(1,n),c,rgrid,nw,bkx,lkx,zgrid,nh,bky, &
                       lky,wk,ier)
        do 1450 i=1,nece
        if (eceiter.eq.'pair') then
        call seva2d(bkx,lkx,bky,lky,c,recem(i),zece(i),pds,ier,n111)
        recefc(i,n)=pds(1)
        endif
        if (eceiter.eq.'flux') recefc(i,n)=0.0
 1450   continue
        do 1460 i=1,nece
         call seva2d(bkx,lkx,bky,lky,c,recep(i),zece(i),pds,ier,n111)
         recefc(i,n)=pds(1)-recefc(i,n)
 1460   continue
 1500 continue 
!------------------------------------------------------------------------------
!--    plasma response -----gecempc(nece,nwnh),geceppc(nece,nwnh)            --
!--                 R-, R+ flux from pc,gecepc=geceppc-gecempc               --
!------------------------------------------------------------------------------
      if (eceiter.eq.'pair') then
      do 72070 m=1,nece
        do 72060 i=1,nw
        do 72060 j=1,nh
          k=(i-1)*nh+j
          rdif=recem(m)-rgrid(i)
          zdif=zece(m)-zgrid(j)
          rsum=rdif**2+zdif**2
          if (rsum.gt.1.0e-08) go to 72054
          mk=(i-1)*nh+1
          gecempc(m,k)=gridpc(mk,i)
          go to 72056
72054     continue
          gecempc(m,k)=psical(recem(m),rgrid(i),zdif)*tmu
72056    continue
72060   continue
72070 continue
      endif
      if (eceiter.eq.'flux') gecempc = 0.0
      do 73070 m=1,nece
        do 73060 i=1,nw
        do 73060 j=1,nh
          k=(i-1)*nh+j
          rdif=recep(m)-rgrid(i)
          zdif=zece(m)-zgrid(j)
          rsum=rdif**2+zdif**2
          if (rsum.gt.1.0e-08) go to 73054
          mk=(i-1)*nh+1
          geceppc(m,k)=gridpc(mk,i)
          go to 73056
73054     continue
          geceppc(m,k)=psical(recep(m),rgrid(i),zdif)*tmu
73056    continue
          gecepc(m,k)=geceppc(m,k)-gecempc(m,k)
73060   continue      
73070 continue       
!-----------------------------------------------------------------------
!-- Ohmic coils   receec(nece,nesum)                                 --
!-----------------------------------------------------------------------
      if (iecurr.le.0) go to 1710
      do 1700 n=1,nesum
        call sets2d(gridec(1,n),c,rgrid,nw,bkx,lkx,zgrid,nh,bky, &
                       lky,wk,ier)
        if (eceiter.eq.'pair') then
        do 1650 i=1,nece
        call seva2d(bkx,lkx,bky,lky,c,recem(i),zece(i),pds,ier,n111)
        receec(i,n)=pds(1)
 1650   continue
        endif
        if (eceiter.eq.'flux') then
        do i=1,nece
          receec(i,n)=0.0
        enddo
        endif
        do 1660 i=1,nece
        call seva2d(bkx,lkx,bky,lky,c,recep(i),zece(i),pds,ier,n111)
        receec(i,n)=pds(1)-receec(i,n)
 1660   continue
 1700 continue
 1710 continue
!----------------------------------------------------------------------
      open(unit=nffile,status='old',form='unformatted',err=12917, &
           file='recexx.dat')
      close(unit=nffile,status='delete')
12917  continue
      open(unit=nffile,status='new',form='unformatted', &
           file='recexx.dat')
      nnnece=nece
      write (nffile) nnnece 
      do 2150 ii=1,nece
        rrrecem(ii)=recem(ii)
        rrrecep(ii)=recep(ii)
 2150 continue
      rrreceo=receo
      write (nffile) rrrecem
      write (nffile) rrrecep
      write (nffile) rrreceo
      write (nffile) recebzfc
      write (nffile) gecebzpc
      write (nffile) recebzec
      write (nffile) recefc
      write (nffile) gecepc 
      write (nffile) receec
      close(unit=nffile)
!-----------------------------------------------------------------------
!--    do every time from here
!-----------------------------------------------------------------------
60000  continue
      ksetece=ksetece+1
      if (ksetece.eq.mtxece) ksetece=0
!-------------------------------------------------------------------
!--    get d(psi)/dr and d2(psi)/dr2
!-------------------------------------------------------------------
      call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
      do 208 iw=1,nw
         rw=rgrid(iw)
         rh=zgrid(jo)
         kk=(iw-1)*nh+jo
         call seva2d(bkx,lkx,bky,lky,c,rw,rh,pds,ier,n555)
         if (ier.eq.0) go to 206
         write (nttyo,9910) ier,rw,rh
         return
 206     continue
         dsidr(iw)=pds(2)
         ddsiddr(iw)=pds(5)
 208   continue
!
      if (eceiter.eq.'flux') then
        do i=1,nnece
         rw=recem(i)
         rh=zece(i)
         call seva2d(bkx,lkx,bky,lky,c,rw,rh,pds,ier,n555)
         if (ier.eq.0) go to 211
         write (nttyo,9910) ier,rw,rh
         return
 211     continue
         brspece(jtime,i)=pds(1)
        enddo
        if (idebug.ge.3) write (6,*) 'SETECE, recem/brspece = ',(recem(ii), &
            brspece(jtime,ii),ii=1,nece)
      endif
 9910 format('  error in spline =',i4,' (r,z)= ( ',f5.2,',',f5.2,')')
!----------------------------------------------------------------------
!-- give error (ecebit,ecebzbit)                                     --
!-- if kfixro kfixrece=0, robit and ecebit  from getecer             --
!-- if kfixro kfixrece=1, robit and rmbit, rpbit from k-file         --
!--   ecebit=sqrt((dpsi/dR|m*rmbit)**2+(dpsi/dR|p*rpbit)**2)         --
!--   ecebzbit=(d2psi/dR2)*robit/receo                               --
!----------------------------------------------------------------------
      if (kfitece.eq.1) go to 350
      if (kfixrece.eq.1) then
        do k=1,nece
         kp=iwp(k)
         km=iwm(k)
         ecebit(k)=(dsidr(kp)*rpbit(k))**2
         ecebit(k)=ecebit(k)+(dsidr(km)*rmbit(k))**2
         ecebit(k)=sqrt(ecebit(k))
        enddo
      endif
350   continue
      ddsiddro=ddsiddr(iwo)
      ecebzbit=ddsiddro*robit/receo
      if (kfitece.eq.1) go to 468
      do 360 m=1,nece
!       tdata1=serror*abs(brspece(jtime,m))
        tdata1=eceerror*abs(brspece(jtime,m))
        tdata2=abs(ecebit(m))
        tdata=max(tdata1,tdata2)
        if (tdata.gt.1.0e-10) fwtece(m)=fwtece0(m)/tdata
        if (tdata.le.1.0e-10) fwtece(m)=0.0
  360 continue
      do 466 i=1,nece
        if (fwtece(i).gt.0.0) kece=kece+1
  466 continue
  468 continue
      tdata1=serror*abs(brspecebz(jtime))
      tdata2=abs(ecebzbit)
      tdata=max(tdata1,tdata2)
      if (tdata.gt.1.0e-10) fwtecebz=fwtecebz0/tdata
      if (tdata.le.1.0e-10) fwtecebz=0.0
      if (fwtecebz.gt.0.0) kecebz=kecebz+1
      receoi(nitera)=receo
      do k=1,nece
         if (fwtece(k).gt.1.e-10) then
           recemi(nitera,k)=recem(k)
           recepi(nitera,k)=recep(k)
         else
           recemi(nitera,k)=1.e-10
           recepi(nitera,k)=1.e-10
         endif
      enddo
!
      if (idebug.ge.3) then 
         write (6,*) 'SETECE, serror/eceerror = ',serror,eceerror
         write (6,*) 'SETECE, nece/fwtece = ',nece,(fwtece(i),i=1,nece)
      endif

      DEALLOCATE(dsidr,ddsiddr)
!
      return
 6530 format ('recefile_0',i1,'.ddd')
 6540 format ('recefile_',i2,'.ddd')
      end
      subroutine seter(ypsi,xpsii)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**         setff sets up the basis functions for er.                **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          97/04/23..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
      include 'basiscomdu.f90'
      dimension xpsii(1)
      do i=1,keecur
         xpsii(i) = bserel(keefnc,i,ypsi)
      enddo
      return
      end
      subroutine seterp(ypsi,xpsii)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**         seterp computes derivative of er.                        **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          97/04/25..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
      include 'basiscomdu.f90'
      dimension xpsii(1)
      do i=1,kffcur
	xpsii(i) = bserpel(keefnc,i,ypsi)
      enddo
      return
      end
      subroutine setff(ypsi,xpsii)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**         setff sets up the basis functions for ff-ff(1).          **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          94/03/09..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
      include 'basiscomdu.f90'
      dimension xpsii(1)
      do i=1,kffcur
        xpsii(i) = bsffin(kfffnc,i,ypsi)
      enddo
      return
      end
      subroutine setfp(ypsi,xpsii)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**         setfp sets up the basis functions for ffp.               **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          94/03/07..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
      include 'basiscomdu.f90'
      dimension xpsii(1)
      do i=1,kffcur
        xpsii(i) = bsffel(kfffnc,i,ypsi)
      enddo
      return
      end
      subroutine setfpp(ypsi,xpsii)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**         setfpp computes derivative of ffp.                       **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          94/03/07..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
      include 'basiscomdu.f90'
      dimension xpsii(1)
      do i=1,kffcur
        xpsii(i) = bsffpel(kfffnc,i,ypsi)
      enddo
      return
      end
      subroutine setpp(ypsi,xpsii)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**         setpp sets up the basis functions for pp.                **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          94/03/07..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
      include 'basiscomdu.f90'
      dimension xpsii(1)
      do i=1,kppcur
				 xpsii(i) = bsppel(kppfnc,i,ypsi)
      enddo
			return
      end
      subroutine setppp(ypsi,xpsii)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**         setppp computes derivative of pp.                        **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          94/03/07..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
      include 'basiscomdu.f90'
      dimension xpsii(1)
      do i=1,kppcur
				 xpsii(i) = bspppel(kppfnc,i,ypsi)
      enddo
			return
      end
      subroutine setpr(ypsi,xpsii)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**         setpr sets up the basis functions for p-p(1).            **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          94/03/09..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
      include 'basiscomdu.f90'
      dimension xpsii(1)
      do i=1,kppcur
				 xpsii(i) = bsppin(kppfnc,i,ypsi)
      enddo
			return
      end
      subroutine setpw(ypsi,xpsii)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**         setpw sets up the basis functions for pw-pw(1).          **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          94/03/09..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
      include 'basiscomdu.f90'
      dimension xpsii(1)
      do i=1,kwwcur
				 xpsii(i) = bswwin(kwwfnc,i,ypsi)
      enddo
      end
      subroutine setpwp(ypsi,xpsii)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**         setpwp sets up the basis functions for pwp.              **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          94/03/07..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
      include 'basiscomdu.f90'
      dimension xpsii(1)
      do i=1,kwwcur
				 xpsii(i) = bswwel(kwwfnc,i,ypsi)
      enddo
			return
      end
      subroutine setpwpp(ypsi,xpsii)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**         setpwpp computes derivative of pwp.                      **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          94/03/07..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
      include 'basiscomdu.f90'
      dimension xpsii(1)
      do i=1,kwwcur
				 xpsii(i) = bswwpel(kwwfnc,i,ypsi)
      enddo
			return
      end
      subroutine setstark(jtime)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          setstark obtains the response matrix                    **
!**          for polarimetry measurement                             **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          23/03/90..........first created                         **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      dimension rsplt(2500),zsplt(2500),csplt(2500) &
                ,rrgamin(nstark),zzgamin(nstark)
      character*40 filenmme
      data nset/20/
      save nset
!---------------------------------------------------------------------
!--  try read in, first locally, then efit area                     --
!---------------------------------------------------------------------
      do 20 iname=1,nset
        if (iname.le.9) then
          write(filenmme,6530) iname
        else
          write(filenmme,6540) iname
        endif
        open(unit=nffile, &
             status='old',form='unformatted', &
             file=filenmme,err=10)
        goto 15
   10   continue
        open(unit=nffile, &
             status='old',form='unformatted', &
             file=table_dir(1:ltbdir)//filenmme,err=30)
   15   continue
        read (nffile,err=30) kkstark
        if (kkstark.ne.nstark) then
          close(unit=nffile)
          go to 20
        endif
        read (nffile,err=30)  rrgamin
        read (nffile,err=30)  zzgamin
        do ii=1,nstark
          if (abs(rrgamin(ii)-rrgam(jtime,ii)).gt.1.e-4) go to 20
          if (abs(zzgamin(ii)-zzgam(jtime,ii)).gt.1.e-4) go to 20
        enddo
        read (nffile,err=30) rbrfc
        read (nffile,err=30) rbzfc
        read (nffile,err=30) gbrpc
        read (nffile,err=30) gbzpc
        read (nffile,err=30) rbrec
        read (nffile,err=30) rbzec
        close(unit=nffile)
        go to 25
   20 continue
      go to 30
   25 return
!---------------------------------------------------------------------
!-- response due to F coils                                         --
!---------------------------------------------------------------------
   30 isplit=10
      itot=isplit*isplit
      fitot=itot
      rbrfc = 0.0
      rbzfc = 0.0
      do 60 k=1,mfcoil
      call splitc(isplit,rsplt,zsplt,csplt, &
                  rf(k),zf(k),wf(k),hf(k),af(k),af2(k),cdum)
      do 60 mmm=1,nstark
        if (rrgam(jtime,mmm).le.1.e-8)  goto 60
        brct=0.0
        bzct=0.0
        r1=rrgam(jtime,mmm)
        do 50 l=1,itot
        a=rsplt(l)
        z1=zzgam(jtime,mmm)-zsplt(l)
        brct=brct+br(a,r1,z1)
        bzct=bzct+bz(a,r1,z1)
   50   continue
        kkm=fcid(k)
        rbrfc(mmm,kkm)=rbrfc(mmm,kkm)+brct/fitot*tmu
        rbzfc(mmm,kkm)=rbzfc(mmm,kkm)+bzct/fitot*tmu
   60 continue
!---------------------------------------------------------------------
!--   plasma response                                               --
!---------------------------------------------------------------------
      zzx=0.0
      drzg=(drgrid**2+dzgrid**2)*25.
      isplit=18
      kk = 1
      do 90 mmm=1,nstark
        gbrpc(mmm,kk)=0.0
        gbzpc(mmm,kk)=0.0
        if (rrgam(jtime,mmm).le.1.e-8)  goto 90
        r=rrgam(jtime,mmm)
        do 80 ii=1,nw
          a=rgrid(ii)
          do 80 jj=1,nh
            z=zzgam(jtime,mmm)-zgrid(jj)
            kk=(ii-1)*nh+jj
            dist=(a-r)**2+z**2
            if (dist.lt.drzg) then
   68        call splitc(isplit,rsplt,zsplt,csplt, &
                  rgrid(ii),zgrid(jj),drgrid,dzgrid,zzx,zzx,cdum)
             gbrt=0.0
             gbzt=0.0
             itot=isplit*isplit
             do 72 k=1,itot
               a=rsplt(k)
               z=zzgam(jtime,mmm)-zsplt(k)
               distt=(a-r)**2+z**2
               if (distt.lt.1.e-8.and.isplit.lt.49) then
                 isplit=isplit+2
                 go to 68
               endif
               gbrt=gbrt+br(a,r,z)
               gbzt=gbzt+bz(a,r,z)
   72        continue
             gbrpc(mmm,kk)=gbrt*tmu/float(itot)
             gbzpc(mmm,kk)=gbzt*tmu/float(itot)
            else
             gbrpc(mmm,kk)=br(a,r,z)*tmu
             gbzpc(mmm,kk)=bz(a,r,z)*tmu
            endif
   80 continue
   90 continue
!---------------------------------------------------------------------
!-- E coils response                                                --
!---------------------------------------------------------------------
      isplit=10
      itot=isplit*isplit
      fitot=float(itot)
      do 201 m=1,nstark
      do 170 i=1,nesum
        rbrec(m,i)=0.0
        rbzec(m,i)=0.0
  170 continue
      if (rrgam(jtime,m).le.1.e-8) goto 201
      r1=rrgam(jtime,m)
      do 200 k=1,necoil
        brct=0.0
        bzct=0.0
        call splitc(isplit,rsplt,zsplt,csplt, &
                    re(k),ze(k),we(k),he(k),zzx,zzx,cdum)
        do 180 l=1,itot
          a=rsplt(l)
          z1=zzgam(jtime,m)-zsplt(l)
          brct=brct+br(a,r1,z1)
          bzct=bzct+bz(a,r1,z1)
  180 continue
      brct=brct*tmu/fitot
      bzct=bzct*tmu/fitot
      kkm=ecid(k)
      rbrec(m,kkm)=rbrec(m,kkm)+brct
      rbzec(m,kkm)=rbzec(m,kkm)+bzct
  200 continue
  201 continue
!
! --- write out rstarkxx.dat if flag IOUT contains 8.
!
      if (iand(iout,8).ne.0) then
      open(unit=nffile,status='old',form='unformatted',err=12917, &
           file='rstarkxx.dat')
      close(unit=nffile,status='delete')
12917  continue
      open(unit=nffile,status='new',form='unformatted', &
           file='rstarkxx.dat')
      kkstark=nstark
      write (nffile) kkstark
      do 215 ii=1,nstark
        rrgamin(ii)=rrgam(jtime,ii)
        zzgamin(ii)=zzgam(jtime,ii)
  215 continue
      write (nffile) rrgamin
      write (nffile) zzgamin
      write (nffile) rbrfc
      write (nffile) rbzfc
      write (nffile) gbrpc
      write (nffile) gbzpc
      write (nffile) rbrec
      write (nffile) rbzec
      close(unit=nffile)
      endif
      return
 6530 format ('rs129129_0',i1,'.ddd')
 6540 format ('rs129129_',i2,'.ddd')
      end
      subroutine shape(iges,igmax,kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          shape finds the outermost plasma surface                **
!**          and computes various global plasma parameters.          **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          10/06/83..........first created                         **
!**          10/18/83..........modified distance calculation         **
!**          24/07/85..........revised                               **
!**          24/07/96..........revised by Q.Peng to add Rmidin/out   **
!**          30/01/98..........Q.Peng, calculate vertical stability  **
!**                            parameter n/nc outside of pltout      **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky,psiold,psipold,psipp, &
                worka,zeros,byringr,byringz,xouts,youts,bpoo,bpooz, &
                bpooc,bfpol,cfpol,dfpol,xxtra,yxtra,bpxtra,flxtra,fpxtra
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      include 'basiscomdu.f90'
      common/cwork3/lkx,lky
      common/cwork4/npxtra(nxtram),scraps(nxtram)
      common/wwork1/xlims(5),ylims(5),limtrs,xlmins
      common/adp/ringr(6),ringz(6),ringap
      dimension pds(6),amer(2,2),bmer(2),wmer(2),imer(2),temp(ntime)
      dimension rmid2(2)
      dimension sigams(nstark)
      real*8,dimension(:),allocatable :: xsisii,bpres,cpres, &
                    dpres,sjtli,sjtlir,sjtliz,rjtli,bpresw, &
                    cpresw,dpresw,copyn,cjtli,x,y
      character*30 sfname
      logical byring,double,onedone
      data floorz/-1.366/
      data psitol/1.0e-04/,idiart/1/
			data czero/0.0/
      data double/.false./, onedone/.false./
!
      ALLOCATE(xsisii(nw),bpres(nw),cpres(nw),dpres(nw), &
         sjtli(nw),sjtlir(nw),sjtliz(nw),rjtli(nw), &
         bpresw(nw),cpresw(nw),dpresw(nw),copyn(nwnh), &
         cjtli(nw),x(nw),y(nh))
!-----------------------------------------------------------------------
!-- ringr and ringz define plasma facing surfaces where strike point  --
!-- may contact surface                                               --
!-- byring is true if r>rsep, z<zsep and z>zvs                        --
!-- double refers to double-valued q-profile                          --
!-- onedone=true if the location (psi) of one q=1 surface has been    --
!-- found                                                             --
!-----------------------------------------------------------------------
      oring(iges)=999 ! initialize at ridiculous value
      ringap=999
!
      if (idebug /= 0) write (6,*) 'Enter SHAPE kerror = ', kerror
      if (kerror.gt.0) go to 1500
      if (ivacum.gt.0) go to 1500
      if (iges.gt.1) go to 100
      xguess=(rgrid(1)+rgrid(nw))/2.
      yguess=(zgrid(1)+zgrid(nh))/2.
      xlims(1)=rgrid(2)+0.2*(rgrid(2)-rgrid(1))
      xlims(2)=xlims(1)
      xlims(3)=rgrid(nw-1)-0.2*(rgrid(nw)-rgrid(nw-1))
      xlims(4)=xlims(3)
      xlims(5)=xlims(1)
      ylims(1)=zgrid(2)+0.2*(zgrid(2)-zgrid(1))
      ylims(2)=zgrid(nh-1)-0.2*(zgrid(nh)-zgrid(nh-1))
      ylims(3)=ylims(2)
      ylims(4)=ylims(1)
      ylims(5)=ylims(1)
      xlmins=xlims(1)
      call zlim(zeros,nw,nh,limtrs,xlims,ylims,rgrid,zgrid,limfag)
  100 continue
!----------------------------------------------------------------------
!--  get outermost flux surface and its shape parameters             --
!----------------------------------------------------------------------
      jges=iges
      itime=time(iges)
      sibdry(iges)=psibry
      simagx(iges)=simag
      rout(iges)=(xmin+xmax)/2.0
      zout(iges)=(ymin+ymax)/2.0
      eout(iges)=(ymax-ymin)/(xmax-xmin)
      aout(iges)=100.*(xmax-xmin)/2.0
      rexpmx=xmax+rexpan
      zexpmx=zxmax
      zzmax(nw)=ymax
      rzzmax(nw)=rymax
      xminn=xmin
      xmaxx=xmax
      xoutp=xout(2)-xout(1)
      xout(nfound+1)=xout(2)
      yout(nfound+1)=yout(2)
      do 400 i=2,nfound
        xoutm=xoutp
        xoutp=xout(i+1)-xout(i)
        if (xoutp*xoutm.ge.0.0) go to 400
        if (xoutp.gt.0.) go to 370
        if (abs(xout(i)-xmax).le.1.0e-04) go to 400
        if (abs(yout(i)-zout(iges)).ge.abs(zxmin-zout(iges))) go to &
             400
        if (xout(i).ge.rout(iges)) go to 400
        xminn=xout(i)
        go to 400
  370   if (abs(xout(i)-xmin).le.1.0e-04) go to 400
        if (abs(yout(i)-zout(iges)).ge.abs(zxmax-zout(iges))) go to 400
        if (xout(i).le.rout(iges)) go to 400
        xmaxx=xout(i)
  400 continue
      xndnt(iges)=(xminn-xmin+xmax-xmaxx)/2./aout(iges)*100.
      rout(iges)=100.*rout(iges)
      zout(iges)=100.*zout(iges)
!-----------------------------------------------------------------------
!--  the distance to the top limiter is found only for values of yout(i)
!--   which are greater than yulim below.
!-----------------------------------------------------------------------
      crymin=100.*rymin
      crymax=100.*rymax
      doutu(iges)=(rout(iges)-crymax)/aout(iges)
      doutl(iges)=(rout(iges)-crymin)/aout(iges)
!---------------------------------------------------------------------
!-- set up P' and FF', then integration                             --
!-- ffprim = (RBt) * d/dpsi(RBt)                                    --
!---------------------------------------------------------------------
      if (icurrt.ne.1) go to 7540
      pprime(1)=cratio*sbeta/darea/srma
      ffprim(1)=cratio*srma*2.*salpha/darea*twopi*tmu
      pprime(nw)=pprime(1)
      ffprim(nw)=ffprim(1)
 7540 continue
      if (icurrt.ne.2.and.icurrt.ne.5) go to 7550
      pprime(nw)=ppcurr(x111,kppcur)/darea
      ffprim(nw)=fpcurr(x111,kffcur)/darea*twopi*tmu
      pprime(1)=ppcurr(x000,kppcur)/darea
      ffprim(1)=fpcurr(x000,kffcur)/darea*twopi*tmu
      if (kfffnc.eq.8) then
         ffprec(nw)=fpecrr(x111,kffcur)/darea*twopi*tmu
         ffprec(1)=fpecrr(x000,kffcur)/darea*twopi*tmu
      else
         ffprec(nw)=0.0
         ffprec(1)=0.0
      endif
 7550 continue
      if (icurrt.ne.4) go to 7600
      call currnt(n222,iges,n222,n222,kerror)
      pprime(1)=cratio/darea/rzero
      ffprim(1)=rbetap*cratio*rzero*twopi*tmu/darea
      ffprim(nw)=ffprim(1)*gammaf
      pprime(nw)=pprime(1)*gammap
 7600 continue
!
      do i=2,nw-1
        ii=nw-i+1
        siii=1.0-1./float(nw-1)*(i-1)
        sigrid(ii)=siii
        if (icurrt.ne.2.and.icurrt.ne.5) go to 7792
        pprime(ii)=ppcurr(siii,kppcur)/darea
        ffprim(ii)=fpcurr(siii,kffcur)/darea*twopi*tmu
         if (kfffnc.eq.8) then
           ffprec(ii)=fpecrr(siii,kffcur)/darea*twopi*tmu
         else
           ffprec(ii)=0.0
         endif
 7792   continue
        if (icurrt.ne.4) go to 7794
        pprime(ii)=(1.-siii**enp)**emp*(1.-gammap)+gammap
        ffprim(ii)=ffprim(1)*pprime(ii)
        pprime(ii)=pprime(1)*pprime(ii)
 7794   continue
        if (icurrt.ne.1) go to 7796
        pprime(ii)=pprime(1)
        ffprim(ii)=ffprim(1)
 7796   continue
      enddo
!
      sigrid(1)=0.0
      sigrid(nw)=1.0
      pres(nw)=prbdry
      fpol(nw)=fbrdy*tmu
      sumf=fpol(nw)**2/2.
      delsi=-(psibry-simag)/float(nw-1)
      do i=1,nw-1
        pres(nw-i)=pres(nw-i+1)+0.5*(pprime(nw-i+1)+pprime(nw-i))*delsi
        sumf=sumf+0.5*delsi*(ffprim(nw-i+1)+ffprim(nw-i))
        if(sumf .ge. 0.0) then
          fpol(nw-i)=sqrt(2.*sumf)*fpol(nw)/abs(fpol(nw))
        else
          fpol(nw-i)=fpol(nw)
        endif
      enddo
      call zpline(nw,sigrid,fpol,bbfpol,ccfpol,ddfpol)
!
      vout(iges)=0.0
      rhovn(nw)=0.0
      xym=xout(1)*yout(1)
      yoxm=yout(1)/xout(1)
!------------------------------------------------------------------
!-- integration over z from 0 to bpolz                           --
!------------------------------------------------------------------
      xww=xout(1)
      dyww=yout(1)/float(nh-1)
      bpolzs=0.5*fpol(nw)
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
      yoxm=bpolzs/xout(1)*dyww
!
      area=0.0
      xyma=yout(1)
      zzm=zmaxis-yout(1)
      do 450 i=2,nfound
        xyp=xout(i)*yout(i)
!------------------------------------------------------------------
!-- integration over z from 0 to bpolz                           --
!------------------------------------------------------------------
        xww=xout(i)
        dyww=yout(i)/float(nh-1)
        bpolzs=0.5*fpol(nw)
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
        yoxp=bpolzs/xout(i)*dyww
!
        xypa=yout(i)
        dx=xout(i)-xout(i-1)
        vout(iges)=vout(iges)+(xyp+xym)/2.0*dx
        area=area+(xyma+xypa)/2.0*dx
        rhovn(nw)=rhovn(nw)+(yoxp+yoxm)/2.0*dx
        xym=xyp
        xyma=xypa
        yoxm=yoxp
        zzp=zmaxis-yout(i)
        if (zzp*zzm.le.0.0) then
          slope=(xout(i)-xout(i-1))/(yout(i)-yout(i-1))
          if (xout(i).lt.rmaxis) rminzm=xout(i)+zzp*slope
          if (xout(i).gt.rmaxis) rmaxzm=xout(i)+zzp*slope
        endif
        zzm=zzp
  450 continue
      vout(iges)=abs(vout(iges))*1.0e+06*twopi
      area=abs(area)
      areao(iges)=area*1.0e+04
      rmagx(iges)=rmaxis*100.0
      zmagx(iges)=zmaxis*100.0
      elongm(iges)=emaxis
      qqmagx(iges)=qmaxis
      terror(iges)=errorm
!---------------------------------------------------------------------
!--  gap calculation                                                --
!---------------------------------------------------------------------
      dleft=1.0e+10
      dright=1.0e+10
      dtop=1.0e+10
      dbott=1.0e+10
      seplim(iges)=1000.
      do 841 j=1,limitr-1
      call dslant(xout,yout,nfound,xmin,xmax,ymin,ymax, &
              xlim(j),ylim(j),xlim(j+1),ylim(j+1),disnow)
      if (xlim(j).lt.1.02.and.xlim(j+1).lt.1.02) then
        dleft = min(dleft,disnow)
      endif
      if (ylim(j).gt.1.20.and.ylim(j+1).gt.1.20) then
        dtop = min(dtop,disnow)
      endif
      if (ylim(j).lt.-1.20.and.ylim(j+1).lt.-1.20) then
        dbott = min(dbott,disnow)
      endif
      if (xlim(j).gt.1.70.and.xlim(j+1).gt.1.70) then
        dright = min(dright,disnow)
      endif
  841 continue
      dismin=min(dleft,dright,dtop,dbott)
      oleft(iges)=dleft
      oright(iges)=dright
      otop(iges)=dtop
      obott(iges)=dbott
      seplim(iges)=dismin
      ssep(iges)=40.
      if (abs(dismin).le.0.100) then
        if (dleft.eq.dismin) limloc(iges)='IN '
        if (dright.eq.dismin) limloc(iges)='OUT'
        if (dtop.eq.dismin) limloc(iges)='TOP'
        if (dbott.eq.dismin) limloc(iges)='BOT'
      else
!--------------------------------------------------------------------
!-- diverted configuration                                         --
!--------------------------------------------------------------------
        deltaa=0.0025
        if (delrmax1.lt.deltaa.and.delrmax2.lt.deltaa) then
           limloc(iges)='DN '
           ssep(iges)=max(delrmax1,delrmax2)*100.
           if (zseps(1,iges).lt.zout(iges)) then
             ssep(iges)=-ssep(iges)
           endif
        else
           if (delrmax1.lt.deltaa) then
            if (zseps(1,iges).gt.zout(iges)) then
             limloc(iges)='SNT'
             ssep(iges)=abs(delrmax2)*100.
            else
             limloc(iges)='SNB'
             ssep(iges)=-abs(delrmax2)*100.
            endif
           else
             ssep(iges)=delrmax1*100.
             limloc(iges)='MAR'
           endif
        endif
      endif
!
      xlimxs=0.0
      ylimys=0.0
      if (dismin.lt.0.500) go to 830
      xsepsl=100.
      if (zseps(1,iges).lt.0.0) xsepsl=rseps(1,iges)/100.
      if (zseps(2,iges).lt.0.0) xsepsl=rseps(2,iges)/100.
      call seva2d(bkx,lkx,bky,lky,c,xlim(1),ylim(1),pds,ier,n111)
      silimp=pds(1)-psibry
      do 820 i=2,limitr
        call seva2d(bkx,lkx,bky,lky,c,xlim(i),ylim(i),pds,ier,n111)
        silimm=silimp
        silimp=pds(1)-psibry
        if (silimp*silimm.gt.0.0) go to 820
        dsilim=silimp-silimm
        xlimxx=xlim(i-1)-(xlim(i)-xlim(i-1))/dsilim*silimm
        ylimyy=ylim(i-1)-(ylim(i)-ylim(i-1))/dsilim*silimm
        if (ylimyy.ge.0.0) go to 820
        if (xlimxx.lt.xsepsl) go to 820
        xlimxs=xlimxx
        ylimys=ylimyy
        go to 830
  820 continue
  830 continue
!----------------------------------------------------------------------
!--    find separatrix outside a limited plasma if one exists.       --
!--    if the plasma is diverted skip this calculation.              --
!--    the distances are defaulted to 99.0 cm.                       --
!--    DPSI > PSITOL  diverted plasma                                --
!--           1.e10   well diverted                                  --
!----------------------------------------------------------------------
      m20=-20
      olefs(iges)=-50.0
      orighs(iges)=-50.0
      otops(iges)=-50.0
      obots(iges)=-50.0
      if (itrace.le.0) go to 1085
      if(abs(dpsi).gt.psitol.or.dismin.gt.0.1) go to 1085
      radbou=1.10*xmin
      if (pasmat(iges).lt.-1.e3) then
         nerr=10000
      else
         nerr=0
      endif
      call bound(psi,nw,nh,nwnh,psiots,xmins,xmaxs,ymins, &
      ymaxs,zeros,rgrid,zgrid,xguess,yguess,jges,limtrs,xlims,ylims, &
      xouts,youts,nfouns,xlmins,npoint,rymins,rymaxs,dpsis, &
      zxmins,zxmaxs,nerr,ishot,itime,limfag,radbou,kbound)
      if (nerr.gt.0) then
        olefs(iges)=-89.0
        orighs(iges)=-89.0
        otops(iges)=-89.0
        obots(iges)=-89.0
        seplim(iges)=-89.0
        go to 1085
      endif
      if (abs(dpsis).le.psitol) then
        olefs(iges)=-45.0
        orighs(iges)=-45.0
        otops(iges)=-45.0
        obots(iges)=-45.0
        seplim(iges)=-45.0
        go to 1085
      endif
      call findax(nw,nh,rgrid,zgrid,rmaxis,zmaxis,simag, &
                  psiots ,rseps(1,iges),zseps(1,iges),m20, &
                  xouts,youts,nfouns,psi,xmins,xmaxs,ymins,ymaxs, &
                  zxmins,zxmaxs,rymins,rymaxs,dpsis,bpoo,bpooz, &
                  limtrs,xlims,ylims,limfag)
!---------------------------------------------------------------------
!--  gap calculation                                                --
!---------------------------------------------------------------------
      seplim(iges)=1000.
      do 861 j=1,nfouns-1
      call dslant(xout,yout,nfound,xmin,xmax,ymin,ymax, &
              xouts(j),youts(j),xouts(j+1),youts(j+1),disnow)
      seplim(iges)=-min(abs(seplim(iges)),disnow)
  861 continue
!----------------------------------------------------------------------
!--    start distance calculation                                    --
!----------------------------------------------------------------------
      zhp=zout(iges)*0.01
      radp=rout(iges)*0.01
      do 1080 j=1,nfouns-1
      if(xouts(j).gt.radp)go to 1030
      dxll=(youts(j)-zleft)*(youts(j+1)-zleft)
      if(dxll.gt.0.0)go to 1030
      if(youts(j).eq.zleft)olefs(iges)=(xouts(j)-xleft)*100.0
      if(youts(j+1).eq.zleft)olefs(iges)=(xouts(j+1)-xleft)*100.0
      if(xouts(j).eq.xouts(j+1))olefs(iges)=(xouts(j)-xleft)*100.0
      if(xouts(j).eq.xouts(j+1))go to 1030
      slope=(youts(j+1)-youts(j))/(xouts(j+1)-xouts(j))
      if(slope.eq.0.0)go to 1030
      bincp=youts(j)-slope*xouts(j)
      xl=(zleft-bincp)/slope
      olefs(iges)=(xl-xleft)*100.0
 1030 if(xouts(j).lt.radp)go to 1040
      dxrr=(youts(j+1)-zright)*(youts(j)-zright)
      if(dxrr.gt.0.0)go to 1040
      if(youts(j).eq.zright)orighs(iges)=(xright-xouts(j))*100.0
      if(youts(j+1).eq.zright)orighs(iges)=(xright-xouts(j+1))*100.0
      if(xouts(j).eq.xouts(j+1))orighs(iges)=(xright-xouts(j))*100.0
      if(xouts(j).eq.xouts(j+1))go to 1040
      slope=(youts(j+1)-youts(j))/(xouts(j+1)-xouts(j))
      if(slope.eq.0.0)go to 1040
      bincp=youts(j)-slope*xouts(j)
      xr=(zright-bincp)/slope
      orighs(iges)=(xright-xr)*100.0
 1040 if(youts(j).lt.zhp)go to 1050
      dytt=(xouts(j)-xztop)*(xouts(j+1)-xztop)
      if(dytt.gt.0.0)go to 1050
      if(xouts(j).eq.xouts(j+1))otops(iges)=(ytop-youts(j))*100.0
      if(xouts(j).eq.xouts(j+1))go to 1050
      slope=(youts(j+1)-youts(j))/(xouts(j+1)-xouts(j))
      bincp=youts(j)-slope*xouts(j)
      yt=slope*xztop+bincp
      otops(iges)=(ytop-yt)*100.0
 1050 continue
      if(youts(j).gt.zhp)go to 1080
      dybb=(xouts(j)-xzbot)*(xouts(j+1)-xzbot)
      if(dybb.gt.0.0)go to 1080
      if(xouts(j).eq.xouts(j+1))obots(iges)=(youts(j)-ybot)*100.0
      if(xouts(j).eq.xouts(j+1))go to 1080
      slope=(youts(j+1)-youts(j))/(xouts(j+1)-xouts(j))
      bincp=youts(j)-slope*xouts(j)
      yb=slope*xzbot+bincp
      obots(iges)=(yb-ybot)*100.0
 1080 continue
!
 1085 continue
      call chisqr(iges)
      nnn=1
      call betali(iges,rgrid,zgrid,nnn,kerror)
      peak(iges)=pres(1)/(.667*wplasm(iges)/(vout(iges)/1.e6))
        do 1088 i=2,nw
          if (rzzmax(i).gt.0.0) go to 1090
 1088   continue
 1090   continue
        amer(1,1)=2.*(rzzmax(1)-rzzmax(i))
        amer(1,2)=2.*(zzmax(1)-zzmax(i))
        amer(2,1)=2.*(rzzmax(1)-rzzmax(i))
        amer(2,2)=2.*(zzmax(1)-(zzmax(1)-zzmax(i)))
        bmer(1)=rzzmax(1)**2+zzmax(1)**2-rzzmax(i)**2-(zzmax(1)- &
                                                  zzmax(i))**2
        bmer(2)=rzzmax(1)**2+zzmax(1)**2-rzzmax(i)**2-(zzmax(1)- &
                                                  zzmax(i))**2
        n22=2
        x11=-1.0
        call decomp(n22,n22,amer,x11 ,imer,wmer)
        call solve(n22,n22,amer,bmer,imer)
        rmer=sqrt((bmer(1)-rzzmax(1))**2+(bmer(2)-zzmax(1))**2)
        i=i+1
        if (rzzmax(i).gt.0.0.and.rzzmax(i).lt.rzzmax(i-1)) then
        amer(1,1)=2.*(rzzmax(1)-rzzmax(i))
        amer(1,2)=2.*(zzmax(1)-zzmax(i))
        amer(2,1)=2.*(rzzmax(1)-rzzmax(i))
        amer(2,2)=2.*(zzmax(1)-(zzmax(1)-zzmax(i)))
        bmer(1)=rzzmax(1)**2+zzmax(1)**2-rzzmax(i)**2-(zzmax(1)- &
                                                  zzmax(i))**2
        bmer(2)=rzzmax(1)**2+zzmax(1)**2-rzzmax(i)**2-(zzmax(1)- &
                                                  zzmax(i))**2
        call decomp(n22,n22,amer,x11 ,imer,wmer)
        call solve(n22,n22,amer,bmer,imer)
!-----------------------------------------------------------------------
!--  need check for RMER 10/91 llao                                   --
!-----------------------------------------------------------------------
        rmer=(sqrt((bmer(1)-rzzmax(1))**2+(bmer(2)-zzmax(1))**2) &
              +rmer)/2.
        endif
        qmerci(iges)=2./(1.+elongm(iges)**2)-2.*(elongm(iges)-1.) &
                   *betap(iges)/elongm(iges)**2/(1.+elongm(iges)) &
                   +(elongm(iges)**2-1.)/(elongm(iges)**2+1.) &
                   *rmaxis/rmer
        if (qmerci(iges).gt.1.e-10) then
          qmerci(iges)=sqrt(1./qmerci(iges))
        else
          qmerci(iges)=-99.
        endif
!
      do 1102 i=1,nw-1
        xsisii(i)=float(i-1)/float(nw-1)
 1102 continue
      xsisii(nw)=1.
!-----------------------------------------------------------------------
!--  write out S(shot).(time)_X files in flux space                   --
!-----------------------------------------------------------------------
      if (kwripre.eq.2) then
        call getfnmd('s',ishot,itime,sfname)
        if (npress.gt.0) then
          sfname=sfname(1:13)//'_presd'
          open(unit=74,status='old',file=sfname,err=12918)
          close(unit=74,status='delete')
12918     continue
          open(unit=74,status='new',file=sfname                       )
          do 42272 i=1,npress
           xrpres=-rpress(i)
           write (74,*) xrpres,pressr(i),xdum,xdum
42272     continue
          close(unit=74)
        endif
        if (nbeam.gt.0) then
          sfname=sfname(1:13)//'_pbeam'
          open(unit=74,status='old',file=sfname,err=12919)
          close(unit=74,status='delete')
12919     continue
          open(unit=74,status='new',file=sfname                       )
          do 42274 i=1,nbeam
           write (74,*) sibeam(i),pbeam(i),xdum,xdum
42274     continue
          close(unit=74)
        endif
        sfname=sfname(1:13)//'_qpsi'
          open(unit=74,status='old',file=sfname,err=12920)
          close(unit=74,status='delete')
12920     continue
          open(unit=74,status='new',file=sfname                       )
        do 42276 i=1,nw
          write (74,*) xsisii(i),qpsi(i),xdum,xdum
42276   continue
        close(unit=74)
        sfname=sfname(1:13)//'_jor'
          open(unit=74,status='old',file=sfname,err=12921)
          close(unit=74,status='delete')
12921     continue
          open(unit=74,status='new',file=sfname                       )
        do 42278 i=1,nw
          cjorka=cjor(i)/1000.
          write (74,*) xsisii(i),cjorka,xdum,xdum
42278   continue
        close(unit=74)
        sfname=sfname(1:13)//'_jorec'
          open(unit=74,status='old',file=sfname,err=17921)
          close(unit=74,status='delete')
17921     continue
          open(unit=74,status='new',file=sfname                       )
        do i=1,nw
          cjorka=cjorec(i)/1000.
          write (74,*) xsisii(i),cjorka,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_cjmse'
          open(unit=74,status='old',file=sfname,err=18031)
          close(unit=74,status='delete')
18031     continue
          open(unit=74,status='new',file=sfname                       )
        if (ishot.le.97400) then
          mcentral=15
        else
          mcentral=10
        endif
        do i=1,mcentral+1
          cjorka=cjmse(i)/1000.
          write (74,*) sigam(i),cjorka,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_cjmsec'
          open(unit=74,status='old',file=sfname,err=18033)
          close(unit=74,status='delete')
18033     continue
          open(unit=74,status='new',file=sfname                       )
        if (ishot.le.97400) then
          mcentral=15
        else
          mcentral=10
        endif
        do i=1,mcentral+1
          cjorka=cjmsec(i)/1000.
          write (74,*) sigam(i),cjorka,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_bzmse'
          open(unit=74,status='old',file=sfname,err=18035)
          close(unit=74,status='delete')
18035     continue
          open(unit=74,status='new',file=sfname                       )
        do i=1,nstark
          cjorka=bzmse(i)
          write (74,*)  sigam(i),cjorka,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_bzmsec'
          open(unit=74,status='old',file=sfname,err=18037)
          close(unit=74,status='delete')
18037     continue
          open(unit=74,status='new',file=sfname                       )
        do i=1,nstark
          cjorka=bzmsec(i)
          write (74,*)  sigam(i),cjorka,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_chigam'
          open(unit=74,status='old',file=sfname,err=18039)
          close(unit=74,status='delete')
18039     continue
          open(unit=74,status='new',file=sfname                       )
        do i=1,nstark
          sigams(i)=sigam(i)
        enddo
        do j=1,nstark
        sigamnow=1.e9
        do i=1,nstark
          if (sigams(i).lt.sigamnow.and.fwtgam(i).gt.0.0) then
          if (rrgam(iges,i).gt.0.0) then
          sigamnow=sigams(i)
          cjorka=chigam(i)
          inow=i
          endif
          endif
        enddo
        write (74,*)  sigamnow,cjorka,xdum,xdum
        sigams(inow)=1.e10
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_presf'
          open(unit=74,status='old',file=sfname,err=12922)
          close(unit=74,status='delete')
12922     continue
          open(unit=74,status='new',file=sfname                       )
        do 42280 i=1,nw
          write (74,*) xsisii(i),pres(i),xdum,xdum
42280   continue
        close(unit=74)
        sfname=sfname(1:13)//'_prespf'
          open(unit=74,status='old',file=sfname,err=12923)
          close(unit=74,status='delete')
12923     continue
          open(unit=74,status='new',file=sfname                       )
        do 42282 i=1,nw
          write (74,*) xsisii(i),pprime(i),xdum,xdum
42282   continue
        close(unit=74)
        sfname=sfname(1:13)//'_prespfn'
          open(unit=74,status='old',file=sfname,err=12924)
          close(unit=74,status='delete')
12924     continue
          open(unit=74,status='new',file=sfname                       )
        do 42284 i=1,nw
          ppnow=pprime(i)/pprime(1)
          write (74,*) xsisii(i),ppnow,xdum,xdum
42284   continue
        close(unit=74)
      endif
!
      psiq1=-1000.
! if two q=1 surfaces, both psi values will be encoded in psiq1
! one to left of decimal and one to the right
!
      idoqn=1
      call zpline(nw,xsisii,qpsi,bfpol,cfpol,dfpol)
      do 1103 i=1,nw
         qpwant=speval(nw,xsisii(i),xsisii,qpsi,bfpol,cfpol,dfpol)
         if (qpwant.le.0.0) then
           idoqn=2
           go to 1107
         endif
 1103 continue
 1107 continue
!
      if (idebug >= 2) write (6,*) 'SHAPE q=1,2,3 q= ',idoqn,qpsi(1),qpsi(nw)
      nnn=1
      d11=30.
      d22=0.03
      d33=0.01
      nzz=0
      n22=2
      zxx=0.0
      do 1104 i=1,3
        pds(i)=100.
 1104 continue
      iend=3
      if (idoqn.eq.1) then
        call zpline(nw,qpsi,xsisii,bfpol,cfpol,dfpol)
      endif
      if (qpsi(nw).gt.1.0) iend=1
      if (qpsi(nw).gt.2.0) iend=2
      if (qpsi(nw).gt.3.0) iend=3
      jstart=2
      double=(idoqn.eq.2).and.(qpsi(1).gt.1.)
      do 1108 i=1,iend
        qwant=i
        if (idoqn.eq.1.and.i.ge.2) then
         if (qwant.lt.qpsi(1)+0.001) go to 1108
         siwant=seval(nw,qwant,qpsi,xsisii,bfpol,cfpol,dfpol)
        else
40009   continue
         do 40010 jjj=jstart,nw
          jj=jjj
          qppp=qwant-qpsi(jj)
          qmmm=qwant-qpsi(jj-1)
          if (qppp*qmmm.le.0.0) then
            siwant=xsisii(jj-1)+ (xsisii(jj)-xsisii(jj-1))/ &
                   (qpsi(jj)-qpsi(jj-1))*qmmm
            if(jstart.eq.2)onedone=.true.
            jstart=jj+1
            go to 1105
          endif
40010    continue
         go to 1108
        endif
 1105   continue
        if (idebug >= 2) then 
            write (6,*) 'i, iend, sw, qw = ',i, iend,siwant, qwant
        endif
        if (siwant.lt.0.0) go to 1108
        siwant=simag-siwant*(simag-psibry)
        call surfac(siwant,psi,nw,nh,rgrid,zgrid,xxtra(1,1),yxtra(1,1), &
                    nfind,npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nzz &
                    ,rmaxis,zmaxis,negcur)
        if (nfind.le.40.and.icntour.eq.0) then
        if (idebug >= 2) write (6,*) ' SHAPE/SURFAC kerror,i,nfind,qp,qm,si = ', &
                            kerror,i,nfind,qppp,qmmm,siwant
        call cntour(rmaxis,zmaxis,siwant,rqmin,rqmax,ycmin,ycmax, &
                    yxcmin,yxcmax,xycmin,xycmax,d11,drgrid,d22, &
                    d33 ,d33 ,xmin,xmax,ymin,ymax,nzz,iautoc, &
                    xxtra(1,1),yxtra(1,1),nfind,rgrid,nw,zgrid,nh, &
                    c,n22,nh2,nttyo,npoint, &
                    negcur,bkx,lkx,bky,lky,kerrorq)
        if (idebug >= 2) write (6,*) ' SHAPE/CNTOUR kerrorq,i,nfind = ', &
                            kerrorq,i,nfind
 
        else
        rqmax=xxtra(1,1)
        rqmin=rqmax
        do 1106 k=2,nfind
          rqmax=max(xxtra(k,1),rqmax)
          rqmin=min(xxtra(k,1),rqmin)
 1106   continue
        endif
        pds(i)=50.*(rqmax-rqmin)
         if (i.eq.1) then
         if(.not.double)psiq1=siwant
!    
!     This code will not compile on the HP or SuperCard.  The logic
!     needs to be fixed so that it will compile and then this code
!     can be put back in.
!    
       if(double.and.(.not.onedone))psiq1=psiq1+siwant	! second value
        endif
 1108 continue
      aaq1(iges)=pds(1)
      aaq2(iges)=pds(2)
      aaq3(iges)=pds(3)
!---------------------------------------------------------------------
!-- Get 3/2 and 2/1 surface information        2002Jan25            --
!---------------------------------------------------------------------
      if (idebug >= 2) write (6,*) 'SHAPE q=3/2, 2/1'
      iend=2
      psin32(iges)=-99.0
      rq32in(iges)=-99.0
      psin21(iges)=-99.0
      rq21top(iges)=-99.0
      do i=1,iend
       if (i.eq.1) qwant=1.5
       if (i.eq.2) qwant=2.0
        do j=1,nw-1
          jj=nw-j+1
          qppp=qwant-qpsi(jj)
          qmmm=qwant-qpsi(jj-1)
          if (qppp*qmmm.le.0.0) then
            siwant=xsisii(jj-1)+ (xsisii(jj)-xsisii(jj-1))/ &
                   (qpsi(jj)-qpsi(jj-1))*qmmm
            psiwan=simag-siwant*(simag-psibry)
        call surfac(psiwan,psi,nw,nh,rgrid,zgrid,xxtra(1,1),yxtra(1,1), &
                    nfind,npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nzz &
                    ,rmaxis,zmaxis,negcur)
            if (nfind.le.40.and.icntour.eq.0) then
        if (idebug >= 2) write (6,*) ' SHAPE/SURFAC kerror,i,nfind,qp,qm,si = ', &
                            kerror,i,nfind,qppp,qmmm,psiwan
        call cntour(rmaxis,zmaxis,psiwan,rqmin,rqmax,ycmin,ycmax, &
                    yxcmin,yxcmax,xycmin,xycmax,d11,drgrid,d22, &
                    d33 ,d33 ,xmin,xmax,ymin,ymax,nzz,iautoc, &
                    xxtra(1,1),yxtra(1,1),nfind,rgrid,nw,zgrid,nh, &
                    c,n22,nh2,nttyo,npoint, &
                    negcur,bkx,lkx,bky,lky,kerrorq)
        if (idebug >= 2) write (6,*) ' SHAPE/CNTOUR kerrorq,i,nfind = ', &
                            kerrorq,i,nfind
            endif
            if (i.eq.1) then
              psin32(iges)=siwant
              rqmin=xxtra(1,1)
              do k=2,nfind
                rqmin=min(xxtra(k,1),rqmin)
              enddo
              rq32in(iges)=100.0*rqmin
            endif
            if (i.eq.2) then
              psin21(iges)=siwant
              zqmax=yxtra(1,1)
              do k=2,nfind
                if (zqmax.lt.yxtra(k,1)) then
                  rqmax=xxtra(k,1)
                  zqmax=yxtra(k,1)
                endif
              enddo
              rq21top(iges)=100.0*rqmax
            endif
            go to 20225
          endif
        enddo
20225   continue
      enddo
!
      call zpline(nw,xsisii,qpsi,bfpol,cfpol,dfpol)
      siwant=0.95
      qpsib(iges)=seval(nw,siwant,xsisii,qpsi,bfpol,cfpol,dfpol)
      shearb(iges)=speval(nw,siwant,xsisii,qpsi,bfpol,cfpol,dfpol)
      shearb(iges)=shearb(iges)/qpsib(iges)
      siwant=0.01
      qpsic=seval(nw,siwant,xsisii,qpsi,bfpol,cfpol,dfpol)
      shearc=speval(nw,siwant,xsisii,qpsi,bfpol,cfpol,dfpol)
      shearc=shearc/qpsic
!----------------------------------------------------------------------
!--  if not diverted, get q at PSIWANT by interpolation              --
!----------------------------------------------------------------------
      if ((abs(psiwant-1.0).le.1.e-05).or.(abs(dismin).le.0.100) &
        .or.(psiwant.le.1.e-5)) then
        siwwww=psiwant
        if (siwwww.le.1.e-5) siwwww=1.e-5
        qsiwant(iges)=seval(nw,siwwww,xsisii,qpsi,bfpol,cfpol,dfpol)
      else
       siwant=simag+psiwant*(psibry-simag)
       call surfac(siwant,psi,nw,nh,rgrid,zgrid,bfpol,dfpol,nfounc &
                    ,npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nnn, &
                    rmaxis,zmaxis,negcur)
       do 51900 k=1,nfounc
        cfpol(k)=1./bfpol(k)**2
51900  continue
       call fluxav(cfpol,bfpol,dfpol,nfounc,psi,rgrid,nw,zgrid,nh, &
                  r2sdry(i),nzz ,sdlobp,sdlbp)
       do 51920 k=1,nfounc
        cfpol(k)=1./bfpol(k)
51920  continue
       call fluxav(cfpol,bfpol,dfpol,nfounc,psi,rgrid,nw,zgrid,nh, &
                  r1sdry(i),nzz ,sdlobp,sdlbp)
       if (kvtor.gt.0) then
         do k=1,nfounc
           cfpol(k)= bfpol(k)**2
         enddo
         call fluxav(cfpol,bfpol,dfpol,nfounc,psi,rgrid,nw,zgrid, &
                  nh,r2wdry,nzz ,sdlobp,sdlbp)
         do k=1,nfounc
           cfpol(k)= ((bfpol(k)/rvtor)**2-1.)**2
         enddo
         call fluxav(cfpol,bfpol,dfpol,nfounc,psi,rgrid,nw,zgrid, &
                  nh,r4wdry,nzz ,sdlobp,sdlbp)
         if (kvtor.eq.3) then
           prew0=pwcurr(psiwant,kwwcur)
           pres0=prcurr(psiwant,kppcur)
           if (abs(pres0).gt.1.e-10) then
             pwop0=prew0/pres0
           else
             pwop0=0.0
           endif
           do k=1,nfounc
             cfpol(k)= ((bfpol(k)/rvtor)**2-1.)
             cfpol(k)= cfpol(k)*exp(pwop0*cfpol(k))
           enddo
           call fluxav(cfpol,bfpol,dfpol,nfounc,psi,rgrid,nw,zgrid, &
                  nh,rp2wdry,nzz ,sdlobp,sdlbp)
           do k=1,nfounc
             cfpol(k)= ((bfpol(k)/rvtor)**2-1.)
             cfpol(k)= exp(pwop0*cfpol(k))
           enddo
           call fluxav(cfpol,bfpol,dfpol,nfounc,psi,rgrid,nw,zgrid, &
                  nh,rpwdry,nzz ,sdlobp,sdlbp)
         endif
       endif
       r2surs = r2sdry(i)*sdlobp
       fpnow = ffcurr(psiwant,kffcur)
       fpnow = fpnow*tmu
       qsiwant(iges) = abs(fpnow)/twopi*r2surs
      endif
!
      call zpline(nw,xsisii,cjor,bfpol,cfpol,dfpol)
      siwant=0.95
      cjor95(iges)=seval(nw,siwant,xsisii,cjor,bfpol,cfpol,dfpol)
!-----------------------------------------------------------------------
!--   shear=dq/d(sqrtVn)/q and J at psiwant                           --
!--   normalize all J's to I/area                         96/04/11    --
!-----------------------------------------------------------------------
      cjorsw(iges)=seval(nw,psiwant,xsisii,cjor,bfpol,cfpol,dfpol)
      if (abs(cpasma(iges)).gt.1.e-1) then
        cjor0(iges)=cjor(1)*areao(iges)/cpasma(iges)/1.0e4
        cjor95(iges)=cjor95(iges)*areao(iges)/cpasma(iges)/1.0e4
        cjorsw(iges)=cjorsw(iges)*areao(iges)/cpasma(iges)/1.0e4
        siwant=0.995
        cjor99(iges)=seval(nw,siwant,xsisii,cjor,bfpol,cfpol,dfpol)
        cjor99(iges)=cjor99(iges)*areao(iges)/cpasma(iges)/1.0e4
      endif
      do 1113 i=1,nw
          bpres(i)=sqrt(volp(i)/volp(nw))
 1113 continue
      call zpline(nw,xsisii,bpres,bfpol,cfpol,dfpol)
      siwant=0.95
      xsi95=seval(nw,siwant,xsisii,bpres,bfpol,cfpol,dfpol)
      siwant=0.01
      xsi01=seval(nw,siwant,xsisii,bpres,bfpol,cfpol,dfpol)
      xsiwww=psiwant
      if (xsiwww.le.1.e-5) xsiwww=1.e-5
      xsiww=seval(nw,xsiwww,xsisii,bpres,bfpol,cfpol,dfpol)
      call zpline(nw,bpres,qpsi,bfpol,cfpol,dfpol)
      ssiwant(iges)=speval(nw,xsiww,bpres,qpsi,bfpol,cfpol,dfpol) &
                    /qsiwant(iges)
      ssi95(iges)=speval(nw,xsi95,bpres,qpsi,bfpol,cfpol,dfpol) &
                    /qpsib(iges)
      ssi01=speval(nw,xsi01,bpres,qpsi,bfpol,cfpol,dfpol) &
                    /qpsic
!-------------------------------------------------------------------
!-- get average current density in outer at 95% flux              --
!-------------------------------------------------------------------
       siavej=0.95
       siwant=siavej
       siwant=simag+siwant*(psibry-simag)
       call surfac(siwant,psi,nw,nh,rgrid,zgrid,bfpol,dfpol,nfounc &
                    ,npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nnn, &
                    rmaxis,zmaxis,negcur)
       call fluxav(bfpol,bfpol,dfpol,nfounc,psi,rgrid,nw,zgrid,nh, &
                   rxxrry,nzz ,sdlobp,sdlbp)
      area=0.0
      xyma=dfpol(1)
      do i=2,nfounc
        xypa=dfpol(i)
        dx=bfpol(i)-bfpol(i-1)
        area=area+(xyma+xypa)/2.0*dx
        xyma=xypa
      enddo
      area=abs(area)*1.e4
      delarea=areao(iges)-area
      pasnow=sdlbp/tmu/twopi
      cj1now=(cpasma(iges)-pasnow)/delarea
      cj1ave(iges)=cj1now/cpasma(iges)*areao(iges)
      cj1now=cj1now*10.
      if (iand(iout,1).ne.0) then
      write (nout,*) cpasma(iges),pasnow,areao(iges),area,cj1ave(iges) &
                     ,cj1now,xsi95
      endif
!
      call zpline(nw,xsisii,pprime,bfpol,cfpol,dfpol)
      siwant=0.95
      pp95(iges)=seval(nw,siwant,xsisii,pprime,bfpol,cfpol,dfpol)
      pp95(iges)=pp95(iges)/pres(1)*sidif
!-----------------------------------------------------------------------
!--  switch order of separatrix, Z of first < 0                       --
!-----------------------------------------------------------------------
      if (zseps(1,iges).gt.0.0) then
        rtemp=rseps(1,iges)
        rseps(1,iges)=rseps(2,iges)
        rseps(2,iges)=rtemp
        rtemp=zseps(1,iges)
        zseps(1,iges)=zseps(2,iges)
        zseps(2,iges)=rtemp
      endif
      aspect=rout(iges)/aout(iges)
      qsta(iges)=rcentr*abs(bcentr(iges))/tmu/abs(cpasma(iges))* &
                 (eout(iges)**2+1.)/2./aspect**2
      olamda=betap(iges)+ali(iges)/2.
      olamda=1.+0.5*olamda**2
      douta=(doutu(iges)+doutl(iges))/2.

      fkdel=1.24-0.54*eout(iges)+0.3*(eout(iges)**2+douta**2) &
           +0.13*douta
      qsta(iges)=qsta(iges)*fkdel*(1.+olamda/aspect**2)
      sepexp(iges)=-100.
      if (kdata.ne.4) then
        pohm=vloopt(iges)*cpasma(iges)
        pohm=max(czero,pohm)
        pin=pohm+pbinj(iges)
        if (pin.gt.1.e-06) then
          taumhd(iges)=wplasm(iges)/pin*1000.
          taudia(iges)=wplasmd(iges)/pin*1000.
        endif
      endif
!
      call lenco2(xout,yout,nfound,iges)
!----------------------------------------------------------------------
!--  compute fast and thermal stored energy                          --
!--  using Heibrink and Duong's approximate formula   L.Lao 10/91    --
!----------------------------------------------------------------------
       wtherm(iges)=wplasm(iges)
       wfbeam(iges)=0.0
       taujd3(iges)=0.0
       tauthn(iges)=0.0
       if (dco2v(iges,2).lt.1.e+10) go to 99993
       if (pbinj(iges).le.100.) go to 99993
       if (abs(cpasma(iges)).le.10000.) go to 99993
       if (wplasm(iges).le.1000.) go to 99993
       cons1=dco2v(iges,2)*1.e-13                      !normalized density
       cons2=wplasm(iges)/cons1/vout(iges)/3.36/1.e-6  !electron temperature (ev)
       cons3=75000.                                    !w_beam (ev)
       cons4=22.*cons2                                  !w_crit (ev)
       cons5=log(1.+(75000./cons4)**1.5)
       cons6=8.37e-2*(cons2/1000)**1.5/cons1*cons5  !thermalization tau_s (s)
       cons7=sqrt(7.5e4)                               !normalized v_beam
       cons8=sqrt(cons4)                               !normalized v_crit
       cons9=cons7**2-cons7*cons8+cons8**2
       cons10=log((cons7+cons8)**2/cons9)
       cons11=atan((2*cons7-cons8)/1.732/cons8)
       cons12=2./3.*(cons8/cons7)**2
       cons13=1.+cons12*(0.5*cons10-1.732*(3.1416/6+cons11))
       cons14=cons6*3/cons5                            !spitzer tau_se (sec)
       cons15=2.75e-1*pbinj(iges)*cons14*cons13          !stored beam energy (j)
       cons16=wplasm(iges)-cons15                         !thermal energy (j)
       cons17=cons16/wplasm(iges)
       wtherm(iges)=cons16
       wfbeam(iges)=cons15
!-----------------------------------------------------------------------
!-- get thermal and jet-diiid confinement times                       --
!-----------------------------------------------------------------------
       cons21=pohm/1.e6                               !ohmic power in mw
       cons22=taumhd(iges)*cons17
       ptotal = pin/1.e6                              !toal power in mw
       if (ptotal.gt.0.01) then
       consa1=ptotal**0.49
       consa2=(abs(cpasma(iges))/1.e6)**1.08
       consa3=(rout(iges)/100.)**1.45
       taujd3(iges)=110.*consa2*consa3/consa1
       else
         taujd3(iges)=0.0
       endif
       if (taujd3(iges).gt.0.0001) then
        tauthn(iges)=cons22/taujd3(iges)
       else
        tauthn(iges)=0.0
       endif
99993  continue
!----------------------------------------------------------------------------
!-- find the vessel intersecting points if diverted plasmas                --
!----------------------------------------------------------------------------
      dolubaf(iges)=-89.0
      diludom(iges)=-89.0
      dolubafm(iges)=-89.0
      diludomm(iges)=-89.0
      dminux(iges)=-89.0
      dminlx(iges)=-89.0
      ratsol(iges)=-89.0
      rvsiu(iges)=-89.0
      zvsiu(iges)=-89.0
      rvsid(iges)=-89.0
      zvsid(iges)=-89.0
      rvsou(iges)=-89.0
      zvsou(iges)=-89.0
      rvsod(iges)=-89.0
      zvsod(iges)=-89.0
      if (dpsi.le.psitol) go to 6500
      if (limloc(iges).eq.'IN') go to 6500
!----------------------------------------------------------------------
!-- get 1 cm SOL width ratio                                         --
!----------------------------------------------------------------------
      rnow=xmax+0.01
      znow=zxmax
      call seva2d(bkx,lkx,bky,lky,c,rnow,znow,pds,ier,n111)
      psinow=pds(1)
      nmax=(xmin-rgrid(1))/(rgrid(2)-rgrid(1))+1
      znow=zxmin
      rrout=xmin
      psiout=psibry
      dsiout=psinow-psiout
      do i=1,nmax
       k= nmax-i +1
       rrin=rgrid(k)
      call seva2d(bkx,lkx,bky,lky,c,rrin,znow,pds,ier,n111)
       psiin =pds(1)
       dsiin =psinow-psiin
       if (dsiin*dsiout.le.0.0) then
         do n=1,20
          rmm=0.5*(rrin+rrout)
      call seva2d(bkx,lkx,bky,lky,c,rmm ,znow,pds,ier,n111)
          psimm=pds(1)
          dsimm=psinow-psimm
          if (dsimm*dsiin.le.0.0) then
            rrout=rmm
          else
            rrin =rmm
          endif
         enddo
         goto 66010
       endif
      enddo
66010 ratsol(iges)=100.*(xmin-rrin)
!----------------------------------------------------------------------
!-- get distance from x point to limiter                             --
!----------------------------------------------------------------------
      ycut=0.5*ymax
      if (zfsep.gt.ycut) then
      dminux(iges)=1.e10
      do i=1,limitr
       if (ylim(i).gt.ycut) then
        dminow=(rfsep-xlim(i))**2+(zfsep-ylim(i))**2
        dminux(iges)=min(dminux(iges),dminow)
       endif
      enddo
      dminux(iges)=100.0*sqrt(dminux(iges))
      else
      dminlx(iges)=1.e10
      ycut=0.5*ymin
      do i=1,limitr
       if (ylim(i).lt.ycut) then
        dminow=(rfsep-xlim(i))**2+(zfsep-ylim(i))**2
        dminlx(iges)=min(dminlx(iges),dminow)
       endif
      enddo
      dminlx(iges)=100.0*sqrt(dminlx(iges))
      endif
!---------------------------------------------------------------
!-- Get distance from second X point to limiter               --
!---------------------------------------------------------------
      ycut=0.5*ymax
      if (sissep.gt.-1.0e10) then
      if (zssep.gt.ycut) then
      dminux(iges)=1.e10
      do i=1,limitr
       if (ylim(i).gt.ycut) then
        dminow=(rssep-xlim(i))**2+(zssep-ylim(i))**2
        dminux(iges)=min(dminux(iges),dminow)
       endif
      enddo
      dminux(iges)=100.0*sqrt(dminux(iges))
      else
      dminlx(iges)=1.e10
      ycut=0.5*ymin
      do i=1,limitr
       if (ylim(i).lt.ycut) then
        dminow=(rssep-xlim(i))**2+(zssep-ylim(i))**2
        dminlx(iges)=min(dminlx(iges),dminow)
       endif
      enddo
      dminlx(iges)=100.0*sqrt(dminlx(iges))
      endif
      endif
!
      rkkk=xmax
      zkkk=zxmax
      call seva2d(bkx,lkx,bky,lky,c,rkkk,zkkk,pds,ier,n111)
      ssimax=pds(1)
      rkkk=rmaxfs
      zkkk=zmaxfs
      call seva2d(bkx,lkx,bky,lky,c,rkkk,zkkk,pds,ier,n111)
      ssimfs=pds(1)
      ixyz=-2
73723 dmaxfs=0.00001
      dmaxfs0=dmaxfs
      nmaxfs=1
77723 xxtras=rmaxfs+dmaxfs
      yxtras=zmaxfs
      rkkk=xxtras
      zkkk=yxtras
      call seva2d(bkx,lkx,bky,lky,c,rkkk,zkkk,pds,ier,n111)
      ssitra=pds(1)
      ixl=1
      if (pasmat(iges).lt.-1.e3) then
         nerr=10000
      else
         nerr=0
      endif
      call bound(psi,nw,nh,nwnh,psiout,xmin,xmax,ymin,ymax, &
          zeros,rgrid,zgrid,xxtras,yxtras,ixyz,limtrs,xlims,ylims, &
          xxtra(1,ixl),yxtra(1,ixl),npxtra(ixl),xlims(1),nxtrap, &
          rymin,rymax,dpsi,zxmin,zxmax,nerr,ishot,itime, &
          limfag,radum,kbound)
      dis2p=(xxtra(1,ixl)-xxtra(npxtra(ixl),ixl))**2
      dis2p=sqrt(dis2p+(yxtra(1,ixl)-yxtra(npxtra(ixl),ixl))**2)
      if (dis2p.lt.0.1*drgrid) then
        nmaxfs=nmaxfs+1
        dmaxfs=dmaxfs0*nmaxfs
        if (nmaxfs.le.20) goto 77723
      endif
      if (dis2p.lt.0.1*drgrid) go to 73823
      do i=1,npxtra(ixl)-1
        if ((yxtra(i,ixl)-znose)*(yxtra(i+1,ixl)-znose).le.0.0) &
             goto 6093
      enddo
      rsepnose=-999.
      goto 6097
 6093 continue
      rsepnose=xxtra(i,ixl)+(xxtra(i+1,ixl)-xxtra(i,ixl))/ &
             (yxtra(i+1,ixl)-yxtra(i,ixl))*(znose-yxtra(i,ixl))
 6097 continue
      zerovs=1.0
      do 6100 i=1,npxtra(ixl)
        zerold=zerovs
        call zlim(zerovs,nnn,nnn,limitr,xlim,ylim,xxtra(i,ixl), &
                  yxtra(i,ixl),limfag)
        if (zerold.gt.0.01.and.zerovs.lt.0.01) go to 6120
 6100 continue
 6120 continue
      if (zerold.gt.0.01.and.zerovs.lt.0.01) then
        rinvs=xxtra(i-1,ixl)
        zinvs=yxtra(i-1,ixl)
        routvs=xxtra(i,ixl)
        zoutvs=yxtra(i,ixl)
        do 6130 i=1,20
        ravs=0.5*(rinvs+routvs)
        zavs=0.5*(zinvs+zoutvs)
        call zlim(zerovs,nnn,nnn,limitr,xlim,ylim,ravs,zavs,limfag)
        if (zerovs.gt.0.01) then
          rinvs=ravs
          zinvs=zavs
        else
          routvs=ravs
          zoutvs=zavs
        endif
 6130   continue
        rvsout(iges)=ravs*100.
        zvsout(iges)=zavs*100.
      else
        rvsout(iges)=0.
        zvsout(iges)=0.
      endif
73823 if ((zavs*zfsep.le.0.0).and.(ixyz.eq.-2)) then
          ixyz=-1
          goto 73723
      endif
      if (dis2p.lt.0.1*drgrid) go to 6198
      call seva2d(bkx,lkx,bky,lky,c,ravs,zavs,pds,ier,n111)
      ssinow=pds(1)
!------------------------------------------------------------------
!-- get distance between inner leg and upper dome                --
!------------------------------------------------------------------
      ycut=ymax*0.5
      ycutm=ymin*0.5
      if (zavs.gt.ycut.and.ravs.lt.rymax) then
      rvsiu(iges)=rvsout(iges)
      zvsiu(iges)=zvsout(iges)
      if (ishot.lt.100771) go to 6198
      diludom(iges)=1.e6
      do i=1,npxtra(ixl)
        if (yxtra(i,ixl).ge.ycut) then
          dilnow=(rudom -xxtra(i,ixl))**2 + &
                          (zudom -yxtra(i,ixl))**2
          if (dilnow.lt.diludom(iges)) then
            diludom(iges)=dilnow
            zilnow=yxtra(i,ixl)
          endif
        endif
      enddo
      diludom(iges)=100.0*sqrt(diludom(iges))
      if (zilnow.lt.zudom) diludom(iges)=-diludom(iges)
!------------------------------------------------------------------
!-- get distance between outer leg and upper baffle              --
!------------------------------------------------------------------
      elseif (zavs.gt.ycut.and.ravs.gt.rymax) then
      rvsou(iges)=rvsout(iges)
      zvsou(iges)=zvsout(iges)
      if (ishot.lt.100771) go to 6198
      dolubaf(iges)=1.e6
      do i=1,npxtra(ixl)
        if (yxtra(i,ixl).ge.ycut) then
          dolnow=(rubaf -xxtra(i,ixl))**2 + &
                          (zubaf -yxtra(i,ixl))**2
          if (dolnow.lt.dolubaf(iges)) then
            dolubaf(iges)=dolnow
            rolnow=xxtra(i,ixl)
          endif
        endif
      enddo
      dolubaf(iges)=100.0*sqrt(dolubaf(iges))
      if (rolnow.gt.rubaf) dolubaf(iges)=-dolubaf(iges)
!------------------------------------------------------------------
!-- get lower strike points                                      --
!------------------------------------------------------------------
      elseif (zavs.lt.ycutm.and.ravs.gt.rymin) then
        rvsod(iges)=rvsout(iges)
        zvsod(iges)=zvsout(iges)
      elseif (zavs.lt.ycutm.and.ravs.lt.rymin) then
        rvsid(iges)=rvsout(iges)
        zvsid(iges)=zvsout(iges)
      endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ring calculation here!!!
6198  if (lring.eq.0) go to 6199
      if(zvsout(iges).eq.0.)go to 6199
      iring=1
      rxp=rseps(1,iges)/100.
      zxp=zseps(1,iges)/100.
      zrmin=floorz
      byringr(iring)=rxp  ! make first point the xpoint
      byringz(iring)=zxp
      do i=1,npxtra(ixl)
       zrmin=min(zrmin,yxtra(i,ixl))
      enddo
      do i=1,npxtra(ixl)
       byring=((xxtra(i,ixl).gt.rxp).and. &
        (yxtra(i,ixl).lt.zxp).and. &
        (yxtra(i,ixl).ge.zrmin))
       if(byring)then ! byring >>> points are close to the ring
        iring=iring+1
        byringr(iring)=xxtra(i,ixl)       ! save points in byring arrays
        byringz(iring)=yxtra(i,ixl)
       endif
      enddo
      if(iring.lt.3)go to 6199
      call order(byringr,byringz,iring) ! put point in increasing order in z
      call zpline(iring,byringz,byringr,bfpol,cfpol,dfpol)  ! spline r(z)
      zringmin=zrmin
      zringmax=max(ringz(1)+.005,zxp)
      dzring=(zringmax-zringmin)/float(nh2)
      do i=1,nh2
       zval=zringmin+(i-1)*dzring
       yxtra(i,ixl)=zval
       rval=seval(iring,zval,byringz,byringr,bfpol,cfpol,dfpol)
       xxtra(i,ixl)=rval
      enddo
      do i=1,nh2
       byringr(i)=xxtra(i,ixl)  !overwrite byring arrays with splined values
       byringz(i)=yxtra(i,ixl)
      enddo
      xrmin=1e20
      xrmax=0
      rbar=0
      zbar=0
      do i=1,6
       xrmax=max(xrmax,ringr(i))
       xrmin=min(xrmin,ringr(i))
       rbar=rbar+ringr(i)/6.
       zbar=zbar+ringz(i)/6.
      enddo
      dismin=1e20
      dsmin=1e20
      do i=1,6
       dsmin=min(dsmin,sqrt((rbar-ringr(i))**2+(zbar-ringz(i))**2))
      enddo
      call dslant(byringr,byringz,nh2,rxp,xrmax,floorz,zxp, &
        ringr(4),ringz(4),ringr(5),ringz(5),gap1) ! vertical ring face
      call dslant(byringr,byringz,nh2,rxp,xrmax,floorz,zxp, &
        ringr(5),ringz(5),ringr(6),ringz(6),gap2) ! slanted ring face
      ringap=gap2
      if(gap1.lt.gap2)ringap=gap1 ! closest approach, separatrix to ring
      do i=1,nh2
       dismin=min(dismin, &
        sqrt((rbar-byringr(i))**2+(zbar-byringz(i))**2))
      enddo
      if(ringap.eq.gap2)then
       if(dismin.lt.dsmin)gap2=-gap2 ! separatrix intersects ring
      endif
      if(ringap.eq.gap1) then	! check for intersection more carefully	
       do i=1,nh2
        if(dismin.eq.sqrt((rbar-byringr(i))**2+(zbar-byringz(i))**2)) &
        iring=i
       enddo
       do i=1,6
        if((ringz(i)-byringz(iring).ne.0.).and.(ringap.eq.gap1))then
         if(atan2(ringr(i)-byringr(iring),ringr(i)-byringz(iring)).lt. &
         0.)gap1=-gap1	 ! separatrix intersects ring
        endif
       enddo
      endif
      if (abs(gap1).gt.1.e-6.and.abs(gap1).lt.1.e6) then
      oring(iges)=gap2*gap1/abs(gap1)	! closest approach to slanted face
      endif
      iring=0 !start cleaning up byrung arrays for plots in subr expand
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ring calculation end!!!!
6199  ixyz=-2
26199 dminfs=0.0001
      dminfs0=dminfs
      nminfs=1
16199 xxtras=rminfs-dminfs
      yxtras=zminfs
      ixl=1
      if (pasmat(iges).lt.-1.e3) then
         nerr=10000
      else
         nerr=0
      endif
      call bound(psi,nw,nh,nwnh,psiout,xmin,xmax,ymin,ymax, &
          zeros,rgrid,zgrid,xxtras,yxtras,ixyz,limtrs,xlims,ylims, &
          xxtra(1,ixl),yxtra(1,ixl),npxtra(ixl),xlims(1),nxtrap, &
          rymin,rymax,dpsi,zxmin,zxmax,nerr,ishot,itime, &
          limfag,radum,kbound)
      dis2p=(xxtra(1,ixl)-xxtra(npxtra(ixl),ixl))**2
      dis2p=sqrt(dis2p+(yxtra(1,ixl)-yxtra(npxtra(ixl),ixl))**2)
      if (dis2p.lt.0.1*drgrid)then
        nminfs=nminfs+1
        dminfs=dminfs0*nminfs
        if (nminfs.le.20) goto 16199
      endif
      if (dis2p.lt.0.1*drgrid) go to 16392
      zerovs=1.0
      do 6200 i=1,npxtra(ixl)
        zerold=zerovs
        call zlim(zerovs,nnn,nnn,limitr,xlim,ylim,xxtra(i,ixl), &
                  yxtra(i,ixl),limfag)
        if (zerold.gt.0.01.and.zerovs.lt.0.01) go to 6220
 6200 continue
 6220 continue
      if (zerold.gt.0.01.and.zerovs.lt.0.01) then
        rinvs=xxtra(i-1,ixl)
        zinvs=yxtra(i-1,ixl)
        routvs=xxtra(i,ixl)
        zoutvs=yxtra(i,ixl)
        do 6230 i=1,20
        ravs=0.5*(rinvs+routvs)
        zavs=0.5*(zinvs+zoutvs)
        call zlim(zerovs,nnn,nnn,limitr,xlim,ylim,ravs,zavs,limfag)
        if (zerovs.gt.0.01) then
          rinvs=ravs
          zinvs=zavs
        else
          routvs=ravs
          zoutvs=zavs
        endif
 6230   continue
        rvsin(iges)=ravs*100.
        zvsin(iges)=zavs*100.
      else
        rvsin(iges)=0.
        zvsin(iges)=0.
      endif
16392 if ((zavs*zfsep.le.0.0).and.(ixyz.eq.-2)) then
          ixyz=-1
          goto 26199
      endif
      if (dis2p.lt.0.1*drgrid) go to 66500
      call seva2d(bkx,lkx,bky,lky,c,ravs,zavs,pds,ier,n111)
      ssinow=pds(1)
!------------------------------------------------------------------
!-- get distance between inner leg and upper dome                --
!------------------------------------------------------------------
      ycut=ymax*0.5
      ycutm=ymin*0.5
      if (zavs.gt.ycut.and.ravs.lt.rymax) then
      rvsiu(iges)=rvsin(iges)
      zvsiu(iges)=zvsin(iges)
      if (ishot.lt.100771) go to 66500
      diludom(iges)=1.e6
      do i=1,npxtra(ixl)
        if (yxtra(i,ixl).ge.ycut) then
          dilnow=(rudom -xxtra(i,ixl))**2 + &
                          (zudom -yxtra(i,ixl))**2
          if (dilnow.lt.diludom(iges)) then
            diludom(iges)=dilnow
            zilnow=yxtra(i,ixl)
          endif
        endif
      enddo
      diludom(iges)=100.0*sqrt(diludom(iges))
      if (zilnow.lt.zudom) diludom(iges)=-diludom(iges)
!------------------------------------------------------------------
!-- get distance between outer leg and upper baffle              --
!------------------------------------------------------------------
      elseif (zavs.gt.ycut.and.ravs.gt.rymax) then
      rvsou(iges)=rvsin(iges)
      zvsou(iges)=zvsin(iges)
      if (ishot.lt.100771) go to 66500
      dolubaf(iges)=1.e6
      do i=1,npxtra(ixl)
        if (yxtra(i,ixl).ge.ycut) then
          dolnow=(rubaf -xxtra(i,ixl))**2 + &
                          (zubaf -yxtra(i,ixl))**2
          if (dolnow.lt.dolubaf(iges)) then
            dolubaf(iges)=dolnow
            rolnow=xxtra(i,ixl)
          endif
        endif
      enddo
      dolubaf(iges)=100.0*sqrt(dolubaf(iges))
      if (rolnow.gt.rubaf ) dolubaf(iges)=-dolubaf(iges)
!------------------------------------------------------------------
!-- get lower strike points                                      --
!------------------------------------------------------------------
      elseif (zavs.lt.ycutm.and.ravs.gt.rymin) then
        rvsod(iges)=rvsin(iges)
        zvsod(iges)=zvsin(iges)
      elseif (zavs.lt.ycutm.and.ravs.lt.rymin) then
        rvsid(iges)=rvsin(iges)
        zvsid(iges)=zvsin(iges)
      endif
66500 if (sissep.le.-1.e10) go to 6500
!------------------------------------------------------------------
!-- second separatrix distances from outboard side               --
!------------------------------------------------------------------
      rnow=rssep
      znow=zssep
      psinow=sissep
      nmin=(xmax-rgrid(1))/(rgrid(2)-rgrid(1))+2
      iznow=(zxmax-zgrid(1)+1.0e-6)/(zgrid(2)-zgrid(1))+1
      znow=zgrid(iznow)
      rrin=xmax
      psiin=psibry
      dsiin =psinow-psiin
      do k=nmin,nw
       rrout=rgrid(k)
      call seva2d(bkx,lkx,bky,lky,c,rrout,znow,pds,ier,n111)
       psiout=pds(1)
       dsiout=psinow-psiout
       if (dsiin*dsiout.le.0.0) then
         do n=1,20
          rmm=0.5*(rrin+rrout)
      call seva2d(bkx,lkx,bky,lky,c,rmm ,znow,pds,ier,n111)
          psimm=pds(1)
          dsimm=psinow-psimm
          if (dsimm*dsiin.le.0.0) then
            rrout=rmm
          else
            rrin=rmm
          endif
         enddo
         goto 59700
       endif
      enddo
59700 rmaxss=rrout
      zmaxss=znow
      ixyz=-2
66501 xxtras=rmaxss+0.0001
      yxtras=zmaxss
      ixl=1
      if (pasmat(iges).lt.-1.e3) then
         nerr=10000
      else
         nerr=0
      endif
      call bound(psi,nw,nh,nwnh,psiout,xmin,xmax,ymin,ymax, &
          zeros,rgrid,zgrid,xxtras,yxtras,ixyz,limtrs,xlims,ylims, &
          xxtra(1,ixl),yxtra(1,ixl),npxtra(ixl),xlims(1),nxtrap, &
          rymin,rymax,dpsi,zxmin,zxmax,nerr,ishot,itime, &
          limfag,radum,kbound)
      zerovs=1.0
      do 66100 i=1,npxtra(ixl)
        zerold=zerovs
        call zlim(zerovs,nnn,nnn,limitr,xlim,ylim,xxtra(i,ixl), &
                  yxtra(i,ixl),limfag)
        if (zerold.gt.0.01.and.zerovs.lt.0.01) go to 66120
66100 continue
66120 continue
      if (zerold.gt.0.01.and.zerovs.lt.0.01) then
        rinvs=xxtra(i-1,ixl)
        zinvs=yxtra(i-1,ixl)
        routvs=xxtra(i,ixl)
        zoutvs=yxtra(i,ixl)
        do 66130 i=1,20
        ravs=0.5*(rinvs+routvs)
        zavs=0.5*(zinvs+zoutvs)
        call zlim(zerovs,nnn,nnn,limitr,xlim,ylim,ravs,zavs,limfag)
        if (zerovs.gt.0.01) then
          rinvs=ravs
          zinvs=zavs
        else
          routvs=ravs
          zoutvs=zavs
        endif
66130   continue
        if ((zavs*zssep.lt.0.0).and.(ixyz.eq.-2)) then
          ixyz=-1
          goto 66501
        endif
        rvsnow=ravs*100.
        zvsnow=zavs*100.
      else
        rvsnow=0.
        zvsnow=0.
      endif
      call seva2d(bkx,lkx,bky,lky,c,ravs,zavs,pds,ier,n111)
      ssinow=pds(1)
!------------------------------------------------------------------
!-- get distance between inner leg and upper dome                --
!------------------------------------------------------------------
      ycut=ymax*0.5
      ycutm=ymin*0.5
      rsymax=rymax
      rsymin=rymin
      if (zssep.lt.ymin) rsymin=rssep
      if (zssep.gt.ymax) rsymax=rssep
      if (zavs*zssep.lt.0.0) then
        ravssa=0.0
        zavssa=0.0
        go to 66198
      endif
      ravssa=ravs
      zavssa=zavs
      if (zavs.gt.ycut.and.ravs.lt.rsymax) then
      rvsiu(iges)=rvsnow
      zvsiu(iges)=zvsnow
      if (ishot.lt.100771) go to 6500
      diludom(iges)=1.e6
      do i=1,npxtra(ixl)
        if (yxtra(i,ixl).ge.ycut) then
          dilnow=(rudom -xxtra(i,ixl))**2 + &
                          (zudom -yxtra(i,ixl))**2
          if (dilnow.lt.diludom(iges)) then
            diludom(iges)=dilnow
            zilnow=yxtra(i,ixl)
          endif
        endif
      enddo
      diludom(iges)=100.0*sqrt(diludom(iges))
      if (zilnow.lt.zudom) diludom(iges)=-diludom(iges)
!------------------------------------------------------------------
!-- get distance between outer leg and upper baffle              --
!------------------------------------------------------------------
      elseif (zavs.gt.ycut.and.ravs.gt.rsymax) then
      rvsou(iges)=rvsnow
      zvsou(iges)=zvsnow
      if (ishot.lt.100771) go to  6500
      dolubaf(iges)=1.e6
      do i=1,npxtra(ixl)
        if (yxtra(i,ixl).ge.ycut) then
          dolnow=(rubaf -xxtra(i,ixl))**2 + &
                          (zubaf -yxtra(i,ixl))**2
          if (dolnow.lt.dolubaf(iges)) then
            dolubaf(iges)=dolnow
            rolnow=xxtra(i,ixl)
          endif
        endif
      enddo
      dolubaf(iges)=100.0*sqrt(dolubaf(iges))
      if (rolnow.gt.rubaf) dolubaf(iges)=-dolubaf(iges)
!------------------------------------------------------------------
!-- get lower strike points                                      --
!------------------------------------------------------------------
      elseif (zavs.lt.ycutm.and.ravs.gt.rsymin) then
        rvsod(iges)=rvsnow
        zvsod(iges)=zvsnow
      elseif (zavs.lt.ycutm.and.ravs.lt.rsymin) then
        rvsid(iges)=rvsnow
        zvsid(iges)=zvsnow
      endif
!------------------------------------------------------------------
!-- second separatrix distances from inboard side                --
!------------------------------------------------------------------
66198 ixyz=-1
      rnow=rssep
      znow=zssep
      psinow=sissep
      nmax=(xmin-rgrid(1))/(rgrid(2)-rgrid(1))+1
      iznow=(zxmin-zgrid(1)+1.0e-6)/(zgrid(2)-zgrid(1))+1
      znow=zgrid(iznow)
      rrout=xmin
      psiout=psibry
      dsiout=psinow-psiout
      do i=1,nmax
       k= nmax-i +1
       rrin=rgrid(k)
      call seva2d(bkx,lkx,bky,lky,c,rrin,znow,pds,ier,n111)
       psiin =pds(1)
       dsiin =psinow-psiin
       if (dsiin*dsiout.le.0.0) then
         do n=1,20
          rmm=0.5*(rrin+rrout)
      call seva2d(bkx,lkx,bky,lky,c,rmm ,znow,pds,ier,n111)
          psimm=pds(1)
          dsimm=psinow-psimm
          if (dsimm*dsiin.le.0.0) then
            rrout=rmm
          else
            rrin=rmm
          endif
         enddo
         goto 59800
       endif
      enddo
59800 rminss=rrout
      zminss=znow
66199 xxtras=rminss-0.0001
      yxtras=zminss
      ixl=1
      if (pasmat(iges).lt.-1.e3) then
         nerr=10000
      else
         nerr=0
      endif
      call bound(psi,nw,nh,nwnh,psiout,xmin,xmax,ymin,ymax, &
          zeros,rgrid,zgrid,xxtras,yxtras,ixyz,limtrs,xlims,ylims, &
          xxtra(1,ixl),yxtra(1,ixl),npxtra(ixl),xlims(1),nxtrap, &
          rymin,rymax,dpsi,zxmin,zxmax,nerr,ishot,itime, &
          limfag,radum,kbound)
      zerovs=1.0
      do 67100 i=1,npxtra(ixl)
        zerold=zerovs
        call zlim(zerovs,nnn,nnn,limitr,xlim,ylim,xxtra(i,ixl), &
                  yxtra(i,ixl),limfag)
        if (zerold.gt.0.01.and.zerovs.lt.0.01) go to 67120
67100 continue
67120 continue
      if (zerold.gt.0.01.and.zerovs.lt.0.01) then
        rinvs=xxtra(i-1,ixl)
        zinvs=yxtra(i-1,ixl)
        routvs=xxtra(i,ixl)
        zoutvs=yxtra(i,ixl)
        do 67130 i=1,20
        ravs=0.5*(rinvs+routvs)
        zavs=0.5*(zinvs+zoutvs)
        call zlim(zerovs,nnn,nnn,limitr,xlim,ylim,ravs,zavs,limfag)
        if (zerovs.gt.0.01) then
          rinvs=ravs
          zinvs=zavs
        else
          routvs=ravs
          zoutvs=zavs
        endif
67130   continue
        if (((abs(ravssa-ravs).le.1.e-04).and. &
         (abs(zavssa-zavs).le.1.e-04)).and.(ixyz.eq.-1)) then
        ixyz=-2
        goto 66199
        endif
        if ((zavs*zssep.lt.0.0).and.(ixyz.eq.-1)) then
          ixyz=-2
          goto 66199
        endif
        rvsnow=ravs*100.
        zvsnow=zavs*100.
      else
        rvsnow=0.
        zvsnow=0.
      endif
!------------------------------------------------------------------
!-- get distance between inner leg and upper dome                --
!------------------------------------------------------------------
      if (zavs*zssep.lt.0.0) go to 6500
      ycut=ymax*0.5
      ycutm=ymin*0.5
      if (zavs.gt.ycut.and.ravs.lt.rsymax) then
      rvsiu(iges)=rvsnow
      zvsiu(iges)=zvsnow
      if (ishot.lt.100771) go to 6500
      diludom(iges)=1.e6
      do i=1,npxtra(ixl)
        if (yxtra(i,ixl).ge.ycut) then
          dilnow=(rudom -xxtra(i,ixl))**2 + &
                          (zudom -yxtra(i,ixl))**2
          if (dilnow.lt.diludom(iges)) then
            diludom(iges)=dilnow
            zilnow=yxtra(i,ixl)
          endif
        endif
      enddo
      diludom(iges)=100.0*sqrt(diludom(iges))
      if (zilnow.lt.zudom) diludom(iges)=-diludom(iges)
!------------------------------------------------------------------
!-- get distance between outer leg and upper baffle              --
!------------------------------------------------------------------
      elseif (zavs.gt.ycut.and.ravs.gt.rsymax) then
      rvsou(iges)=rvsnow
      zvsou(iges)=zvsnow
      if (ishot.lt.100771) go to  6500
      dolubaf(iges)=1.e6
      do i=1,npxtra(ixl)
        if (yxtra(i,ixl).ge.ycut) then
          dolnow=(rubaf -xxtra(i,ixl))**2 + &
                          (zubaf -yxtra(i,ixl))**2
          if (dolnow.lt.dolubaf(iges)) then
            dolubaf(iges)=dolnow
            rolnow=xxtra(i,ixl)
          endif
        endif
      enddo
      dolubaf(iges)=100.0*sqrt(dolubaf(iges))
      if (rolnow.gt.rubaf) dolubaf(iges)=-dolubaf(iges)
!------------------------------------------------------------------
!-- get lower strike points                                      --
!------------------------------------------------------------------
      elseif (zavs.lt.ycutm.and.ravs.gt.rsymin) then
        rvsod(iges)=rvsnow
        zvsod(iges)=zvsnow
      elseif (zavs.lt.ycutm.and.ravs.lt.rsymin) then
        rvsid(iges)=rvsnow
        zvsid(iges)=zvsnow
      endif
!
 6500 continue
      if (rvsin(iges).gt.rvsout(iges)) then
         rtemp=rvsin(iges)
         ztemp=zvsin(iges)
         rvsin(iges)=rvsout(iges)
         zvsin(iges)=zvsout(iges)
         rvsout(iges)=rtemp
         zvsout(iges)=ztemp
      endif
!---------------------------------------------------------------------------
!-- compute distances at outboard and inboard sides                       --
!---------------------------------------------------------------------------
      if (diludom(iges).gt.0.0) then
      rnow=rudom
      znow=zudom
      call seva2d(bkx,lkx,bky,lky,c,rnow,znow,pds,ier,n111)
      psinow=pds(1)
      nmax=(xmin-rgrid(1))/(rgrid(2)-rgrid(1))+1
      znow=zxmin
      rrout=xmin
      psiout=psibry
      dsiout=psinow-psiout
      do i=1,nmax
       k= nmax-i +1
       rrin=rgrid(k)
      call seva2d(bkx,lkx,bky,lky,c,rrin,znow,pds,ier,n111)
       psiin =pds(1)
       dsiin =psinow-psiin
       if (dsiin*dsiout.le.0.0) then
         do n=1,20
          rmm=0.5*(rrin+rrout)
      call seva2d(bkx,lkx,bky,lky,c,rmm ,znow,pds,ier,n111)
          psimm=pds(1)
          dsimm=psinow-psimm
          if (dsimm*dsiin.le.0.0) then
            rrout=rmm
          else
            rrin=rmm
          endif
         enddo
         goto 66600
       endif
      enddo
66600 diludomm(iges)=100.0*(xmin-rrin)
      endif
!
      if (dolubaf(iges).gt.0.0) then
      rnow=rubaf
      znow=zubaf
      call seva2d(bkx,lkx,bky,lky,c,rnow,znow,pds,ier,n111)
      psinow=pds(1)
      nmin=(xmax-rgrid(1))/(rgrid(2)-rgrid(1))+2
      znow=zxmax
      rrin=xmax
      psiin=psibry
      dsiin =psinow-psiin
      do k=nmin,nw
       rrout=rgrid(k)
      call seva2d(bkx,lkx,bky,lky,c,rrout,znow,pds,ier,n111)
       psiout=pds(1)
       dsiout=psinow-psiout
       if (dsiin*dsiout.le.0.0) then
         do n=1,20
          rmm=0.5*(rrin+rrout)
      call seva2d(bkx,lkx,bky,lky,c,rmm ,znow,pds,ier,n111)
          psimm=pds(1)
          dsimm=psinow-psimm
          if (dsimm*dsiin.le.0.0) then
            rrout=rmm
          else
            rrin=rmm
          endif
         enddo
         goto 66700
       endif
      enddo
66700 dolubafm(iges)=100.*(rrout-xmax)
      endif
!---------------------------------------------------------------------------
!-- compute flux expansion parameter                                      --
!---------------------------------------------------------------------------
      fexpan=-10.
      fexpvs=-10.
      if (limloc(iges).ne.'DN '.and.limloc(iges).ne.'SNB'.and. &
          limloc(iges).ne.'SNT') go to 76000
      fxrmin=rexpmx-rexpan
      fxzmin=zexpmx
      if (limloc(iges).eq.'DN '.or.limloc(iges).eq.'SNB') then
        ixyz=-2
        zexpmx=zexpmx+0.00001
        fxrmax=rseps(1,iges)/100.
        fxzmax=zseps(1,iges)/100.
      else
        ixyz=-1
        zexpmx=zexpmx-0.0001
        fxrmax=rseps(2,iges)/100.
        fxzmax=zseps(2,iges)/100.
      endif
      if (pasmat(iges).lt.-1.e3) then
         nerr=10000
      else
         nerr=0
      endif
      ixl=1
      call bound(psi,nw,nh,nwnh,psiout,xmin,xmax,ymin,ymax, &
          zeros,rgrid,zgrid,rexpmx,zexpmx,ixyz,limtrs,xlims,ylims, &
          xxtra(1,ixl),yxtra(1,ixl),npxtra(ixl),xlims(1),nxtrap, &
          rymin,rymax,dpsi,zxmin,zxmax,nerr,ishot,itime, &
          limfag,radum,kbound)
!-------------------------------------------------------------------------
!-- get engineering slot parameter in cm, SEPNOSE                       --
!-------------------------------------------------------------------------
      do i=1,npxtra(ixl)-1
        if ((yxtra(i,ixl)-znose)*(yxtra(i+1,ixl)-znose).le.0.0) &
             goto 6393
      enddo
      sepnose=-999.
      goto 6397
 6393 continue
      rlinnose=xxtra(i,ixl)+(xxtra(i+1,ixl)-xxtra(i,ixl))/ &
             (yxtra(i+1,ixl)-yxtra(i,ixl))*(znose-yxtra(i,ixl))
      sepnose=(rlinnose-rsepnose)*100.
      if (rsepnose.le.0.0) sepnose=-999.
 6397 continue
      fxmin=1.e30
      fxmax=1.e30
      fvsmax=1.e30
      fsrmax=rvsout(iges)/100.
      fszmax=zvsout(iges)/100.
      do i=1,npxtra(ixl)
       fminow=sqrt((fxrmin-xxtra(i,ixl))**2+(fxzmin-yxtra(i,ixl))**2)
       fmanow=sqrt((fxrmax-xxtra(i,ixl))**2+(fxzmax-yxtra(i,ixl))**2)
       fvsnow=sqrt((fsrmax-xxtra(i,ixl))**2+(fszmax-yxtra(i,ixl))**2)
       fxmin=min(fminow,fxmin)
       fxmax=min(fmanow,fxmax)
       fvsmax=min(fvsnow,fvsmax)
      enddo
      if (fxmin.gt.1.e-30) then
        fexpan=fxmax/fxmin
        fexpvs=fvsmax/fxmin
        if (fszmax.ge.0.0) fexpvs=-999.
      endif
76000 continue
!---------------------------------------------------------------------------
!-- trace external field lines                                            --
!---------------------------------------------------------------------------
      if(nextra.eq.0)go to 900
      if(ixstrt.eq.-1)go to 832
      xxtraa=xmax
      yxtraa=zxmax
  827 continue
      go to 834
  832 continue
      xxtraa=xmin
      yxtraa=zxmin
  833 continue
  834 ixtpls=0
      ixyz=-1
  420 continue
      dxtra=scrape/iabs(nextra)
      do 40000 kkk=1,iabs(ixstrt)
      do 875 i=1,iabs(nextra)
      ixl=i+ixtpls*iabs(ixstrt)+(kkk-1)*iabs(nextra)
      scraps(ixl)=1000.*i*dxtra
      if (ixstrt.eq.2) then
      if (kkk.eq.2) then
        xxtraa=xmin
        yxtraa=zxmin
      else
        xxtraa=xmax
        yxtraa=zxmax
      endif
      endif
      xxtras=xxtraa+dxtra*i
      yxtras=yxtraa
      if (kkk.eq.2.or.ixstrt.eq.-1) then
      xxtras=xxtraa-dxtra*i
      yxtras=yxtraa
      endif
  610 continue
      if (xxtras.le.rgrid(1).or.xxtras.ge.rgrid(nw)) go to 875
      if (pasmat(iges).lt.-1.e3) then
         nerr=10000
      else
         nerr=0
      endif
      call bound(psi,nw,nh,nwnh,psiout,xmin,xmax,ymin,ymax, &
          zeros,rgrid,zgrid,xxtras,yxtras,ixyz,limtrs,xlims,ylims, &
          xxtra(1,ixl),yxtra(1,ixl),npxtra(ixl),xlims(1),nxtrap, &
          rymin,rymax,dpsi,zxmin,zxmax,nerr,ishot,itime, &
          limfag,radum,kbound)
      if ((i.gt.1).or.(ixyz.ne.-2)) go to 620
      if (xlimxs.le.0.0) go to 620
      do 615 n=1,npxtra(ixl)
        sepnow=sqrt((xlimxs-xxtra(n,ixl))**2+(ylimys- &
                     yxtra(n,ixl))**2)
        sepexp(iges)=min(abs(sepexp(iges)),sepnow)
  615 continue
      sepexp(iges)=sepexp(iges)*100.
  620 continue
      if (nextra.gt.0) go to 875
      do 623 n=1,npxtra(ixl)
        call seva2d(bkx,lkx,bky,lky,c,xxtra(n,ixl),yxtra(n,ixl), &
                    pds,ier,n333)
        bpxtra(n,ixl)=sqrt(pds(2)**2+pds(3)**2)/xxtra(n,ixl)
  623 continue
      fpxtra(1,ixl)=0.0
      flxtra(1,ixl)=0.0
      do 627 n=2,npxtra(ixl)
        nm1=n-1
        dlpol=sqrt((xxtra(n,ixl)-xxtra(nm1,ixl))**2 &
                   +(yxtra(n,ixl)-yxtra(nm1,ixl))**2)
        btnow=abs(fbrdy)*tmu/2.*(1./xxtra(n,ixl)+1./xxtra(nm1,ixl))
        bpnow=(bpxtra(n,ixl)+bpxtra(nm1,ixl))/2.
        dltol=dlpol*btnow/bpnow
        dlll=sqrt(dlpol**2+dltol**2)
        fpxtra(n,ixl)=fpxtra(nm1,ixl)+dlpol
        flxtra(n,ixl)=flxtra(nm1,ixl)+dlll
  627 continue
  875 continue
40000 continue
      if(ixyz.eq.-2)go to 900
      ixyz=-2
      ixtpls=iabs(nextra)
      go to 420
  900 continue
!-----------------------------------------------------------------------------
!--  compute the diamagnetic flux                                           --
!-----------------------------------------------------------------------------
      vbtot2=0.0
      vbtvac2=0.0
      vbtvac=0.0
      vbtor2=0.0
      volbt=0.0
      cdflux(iges)=0.0
      edflux(iges)=0.0
      btvtot2=0.0
      btvvac2=0.0
      btvtor2=0.0
      volrt=0.0
      areart=0.0
      sumbz2=0.0
      sumbp2=0.0
!
      call zpline(nw,xsisii,pres,bpres,cpres,dpres)
      if (kdopre.gt.0) then
        ii=0
        do i=1,nw,kdopre
          ii=ii+1
          rpress(ii)=-xsisii(i)
          pressr(ii)=pres(i)
          sigpre(ii)=0.1*pres(i)
        enddo
        if (abs(rpress(ii)).ne.xsisii(nw)) then
          ii=ii+1
          i=nw
          rpress(ii)=-xsisii(i)
          pressr(ii)=pres(i)
          sigpre(ii)=0.1*pres(i)
        endif
        npress=ii
      endif
!
      if (kvtor.gt.0) &
        call zpline(nw,xsisii,pressw,bpresw,cpresw,dpresw)
      delerr=0.0
      delerb=0.0
      do 950 i=1,nw
      do 950 j=1,nh
        kk=(i-1)*nh+j
        fnow=fbrdy
        presss=0.0
        prewww=0.0
        prettt=presss
        if (xpsi(kk).gt.1.0) go to 940
        fnow=seval(nw,xpsi(kk),sigrid,fpol,bbfpol,ccfpol,ddfpol)
        presss=seval(nw,xpsi(kk),xsisii,pres,bpres,cpres,dpres) &
                -prbdry
        fnow=fnow/tmu
        prettt=presss
        if (kvtor.gt.0) then
           preww0=seval(nw,xpsi(kk),xsisii,pressw,bpresw,cpresw,dpresw) &
                  -preswb
           press0=presss
           presss=press0+preww0*rgrvt(i)
           prewww=preww0*(rgrid(i)/rvtor)**2*2.
           prew0=preww0
           pres0=press0
           if (abs(pres0).gt.1.e-10) then
               pwop0=prew0/pres0
               pwp0r2=pwop0*rgrvt(i)
           else
               pwop0=0.0
               pwp0r2=0.0
           endif
           if (kvtor.eq.2) then
             presss=press0*(1.+0.5*pwp0r2**2)+preww0*rgrvt(i)
           endif
           if (kvtor.eq.3) then
             presss=press0*exp(pwp0r2)
           endif
           prettt=prewww+presss
        endif
  940   continue
        cdflux(iges)=cdflux(iges)+(fbrdy-fnow)/rgrid(i)*www(kk)
        edflux(iges)=edflux(iges)+(fbrdy**2-fnow**2)/rgrid(i)*www(kk)
        call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n666)
        bpzsq=pds(2)**2/rgrid(i)
        bpolsq=bpzsq+pds(3)**2/rgrid(i)
        sumbz2=sumbz2+bpzsq*www(kk)
        sumbp2=sumbp2+bpolsq*www(kk)
        bpolsq=bpolsq/rgrid(i)
        btvac2=(tmu*fbrdy/rgrid(i))**2
        btttt2=(tmu*fnow/rgrid(i))**2
        enrgy=(2.*twopi*tmu*prettt+bpolsq+btvac2-btttt2)*www(kk)
        volrt=volrt+rgrid(i)*enrgy
        areart=areart+enrgy
        vbtot2=vbtot2+rgrid(i)*www(kk)*(bpolsq+btttt2)
        vbtvac2=vbtvac2+rgrid(i)*www(kk)*btvac2
        vbtvac=vbtvac+rgrid(i)*www(kk)*sqrt(btvac2)
        vbtor2=vbtor2+rgrid(i)*www(kk)*btttt2
        volbt=volbt+rgrid(i)*www(kk)
        btvvac2=btvvac2+2.*twopi*tmu*presss*www(kk)/btvac2*rgrid(i)
        btvtor2=btvtor2+2.*twopi*tmu*presss*www(kk)/btttt2*rgrid(i)
        btvtot2=btvtot2+2.*twopi*tmu*presss*www(kk)/(bpolsq+btttt2)* &
                rgrid(i)
!--------------------------------------------------------------------------
!--  evaluate -del*psi/R/mu0                                             --
!--------------------------------------------------------------------------
        if (i.eq.1.or.i.eq.nw) goto 950
        if (j.eq.1.or.j.eq.nh) goto 950
        kip=i*nh+j
        kim=(i-2)*nh+j
        kjp=(i-1)*nh+j+1
        kjm=(i-1)*nh+j-1
        d2sidr2=(psipold(kip)-2.*psipold(kk)+psipold(kim))/drgrid**2
        d2sidz2=(psipold(kjp)-2.*psipold(kk)+psipold(kjm))/dzgrid**2
        dsidr=(psipold(kip)-psipold(kim))/2./drgrid
        delsbu=d2sidr2-dsidr/rgrid(i)+d2sidz2
        delsbu=-delsbu/rgrid(i)/tmu0
        pcnow=pcurrt(kk)/darea
        delerrb=abs(delsbu-pcnow)
        delerb=max(delerrb,delerb)
!--------------------------------------------------------------------------
!-- total psi for inside plasma only                                     --
!--------------------------------------------------------------------------
        if (xpsi(kk).gt.0.999) go to 950
        d2sidr2=(psiold(kip)-2.*psiold(kk)+psiold(kim))/drgrid**2
        d2sidz2=(psiold(kjp)-2.*psiold(kk)+psiold(kjm))/dzgrid**2
        dsidr=(psiold(kip)-psiold(kim))/2./drgrid
        delssi=d2sidr2-dsidr/rgrid(i)+d2sidz2
        delssi=-delssi/rgrid(i)/tmu0
        delerrx=abs(delssi-pcnow)
        delerr=max(delerrx,delerr)
  950 continue
      curnow=abs(cpasma(iges))/areao(iges)*1.e4
      delerr=delerr/curnow
      delerb=delerb/curnow
      alpha(iges)=2*sumbz2/sumbp2
      rttt(iges)=volrt/areart*100.0
      vbtot2=vbtot2/volbt
      vbtvac2=vbtvac2/volbt
      vbtor2=vbtor2/volbt
      vbtvac=vbtvac/volbt
      btuse=bcentr(iges)*rcentr*100./rout(iges)
      btuse=btuse**2
      vbtot2=betat(iges)*btuse/vbtot2
      vbtvac2=betat(iges)*btuse/vbtvac2
      vbtor2=betat(iges)*btuse/vbtor2
      vbtvac=betat(iges)*btuse/vbtvac**2
      btvvac2=btvvac2*100./volbt
      btvtor2=btvtor2*100./volbt
      btvtot2=btvtot2*100./volbt
      vbtmag=betat(iges)*(rmaxis*100./rout(iges))**2
      betped=0.0
      betnped=0.0
      if (kedgep.gt.0) then
        sssiie=pe_psin-pe_width
        presped=seval(nw,sssiie,xsisii,pres,bpres,cpres,dpres) &
                -prbdry
        betped=presped/btuse*tmu*twopi*2.*100.
        betnped=betped/betat(iges)*betatn
      endif
!--------------------------------------------------------------------
!-- change sign convention of computed diamagnetic flux            --
!-- to be consistent with DIII-D measurements.  The sign now       --
!-- is opposite to that in Nuc Fusion paper      06/23/87          --
!--------------------------------------------------------------------
      cdflux(iges)=-cdflux(iges)*darea*tmu
      edflux(iges)=-edflux(iges)*darea*tmu**2
      sbpli=(s1(iges)+s2(iges)*(1.+rttt(iges)*0.01/rcentr))/4.0
      vbpli=betap(iges)+ali(iges)/2.
      if (kvtor.gt.0) vbpli=vbpli+betapw(iges)
      dbpli(iges)=abs((sbpli-vbpli)/sbpli)
      chidlc=0.0
      if (fwtdlc.le.0.0) go to 960
      chidlc=((diamag(iges)-cdflux(iges))/sigdia(iges))**2
  960 continue
      chitot=chipre+tsaisq(iges)+chidlc
      if (idiart.gt.0) then
      bp2flx=bpolav(iges)**2
      dmui=1.0e+06*diamag(iges)*4.*pi*bcentr(iges)*rcentr &
            /bp2flx/vout(iges)
      rcurrm=rttt(iges)*0.01
      betapd(iges)=s1(iges)/2.+s2(iges)/2.*(1.-rcurrm/rcentr)-dmui
      betatd(iges)=betapd(iges)*bp2flx/bcentr(iges)**2*100.
      betatd(iges)=betatd(iges)*(rout(iges)/100./rcentr)**2
      wplasmd(iges)=1.5*betapd(iges)*bp2flx/2./tmu/2./pi*vout(iges) &
                    /1.0e+06
      if (kdata.ne.4) then
        pohm=vloopt(iges)*cpasma(iges)
        pohm=max(czero,pohm)
        pin=pohm+pbinj(iges)
        if (pin.gt.1.e-06) then
          taudia(iges)=wplasmd(iges)/pin*1000.
        endif
      endif
      endif
!
      xmui=2.*twopi*bcentr(iges)*rcentr*cdflux(iges)/vout(iges) &
         *1.0e+06/bpolav(iges)**2
      exmui=twopi*edflux(iges)/vout(iges) &
         *1.0e+06/bpolav(iges)**2
      sbppa=(s1(iges)+s2(iges)*(1.-rttt(iges)*0.01/rcentr))/2.-xmui
      sbpp=(s1(iges)+s2(iges)*(1.-rttt(iges)*0.01/rcentr))/2.-exmui
      delbp(iges)=abs((sbpp-betap(iges))/sbpp)
      if (eout(iges).gt.1.2) then
        sbli=(s1(iges)/2.+s2(iges)/2.*(1.-rttt(iges)*0.01/rcentr)- &
              s3(iges))/(alpha(iges)-1.)
        delli=abs((sbli-ali(iges))/ali(iges))
        delbp(iges)=max(delbp(iges),delli)
      endif
!
      if (nbskip.le.0) then
        nbabs=1
      else
        nbabs=max(nfound/min(mbdry,nbdrymx)+1,nbskip)
      endif
      jb=0
      iskip=nbabs-1
!-------------------------------------------------------------------
!-- skip points if needed but make sure X point and last point    --
!-- are included                                                  --
!-------------------------------------------------------------------
      do 970 i=1,nfound
        iskip=iskip+1
        dzzz1=abs(zseps(1,iges)/100.-yout(i))
        dzzz2=abs(zseps(2,iges)/100.-yout(i))
        dzzz1=min(dzzz1,dzzz2)
        if ((iskip.eq.nbabs.or.dzzz1.le.1.e-4).or.(i.eq.nfound)) &
        then
         jb=jb+1
         rbbbs(jb)=xout(i)
         zbbbs(jb)=yout(i)
         iskip=0
        endif
  970 continue
      nbbbs=jb
!--------------------------------------------------------------------
!--   rmidin and rmidout are boundary R (rbbbs) at Z (zbbbs) = 0.0 --
!--------------------------------------------------------------------
      i2=1
      do i=1,nbbbs-1
         if (zbbbs(i).eq.0.0) then 
            rmid2(i2)= rbbbs(i)
            i2=i2+1
         else if (zbbbs(i)*zbbbs(i+1).lt.0.0) then
            rmid2(i2)=(rbbbs(i)+rbbbs(i+1))/2.0
            i2=i2+1
         endif
         if (i2.gt.2) go to 971
      enddo
  971 rmidin(iges) = min(rmid2(1),rmid2(2))
      rmidout(iges)= max(rmid2(1),rmid2(2))
!
      if (kprfit.le.0) go to 980
      do 975 i=1,npress
      if (rpress(i).le.0.0) go to 975
      call seva2d(bkx,lkx,bky,lky,c,rpress(i),zpress(i),pds,ier,n111)
      rpress(i)=-(simag-pds(1))/sidif
  975 continue
  980 continue
!----------------------------------------------------------------
!-- Current at Z=Z_Libeam                                      --
!----------------------------------------------------------------
      if (kwripre.gt.0) then
      if (nstark.gt.nmtark) then
        znow=zzgam(iges,nmtark+1)
        delrnow=(xmax-rmaxis)/float(nw-1)
        do i=1,nw
           rnow=rmaxis+(i-1)*delrnow
      call seva2d(bkx,lkx,bky,lky,c,rnow,znow,pds,ier,n666)
           rjtli(i)=rnow
           sjtliz(i)=(pds(2)/rnow-pds(5))/rnow/twopi/tmu/1000.
           sjtlir(i)=-pds(6)/rnow/twopi/tmu/1000.
           sjtli(i)=sjtlir(i)+sjtliz(i)
           xpsikk=(pds(1)-simag)/(psibry-simag)
           if (xpsikk.lt.0.0) xpsikk=0.0
           cjtli(i)=rnow*ppcurr(xpsikk,kppcur) &
                      +fpcurr(xpsikk,kffcur)/rnow
           cjtli(i)=cjtli(i)/darea/1000.
        enddo
        call getfnmd('q',ishot,itime,sfname)
        sfname=sfname(1:13)//'_jtorli'
          open(unit=74,status='old',file=sfname,err=19977)
          close(unit=74,status='delete')
19977     continue
          open(unit=74,status='new',file=sfname                       )
        xdum=0.0
        do i=1,nw
          write (74,*) rjtli(i),cjtli(i),xdum,xdum
        enddo
        close(unit=74)
        call getfnmd('q',ishot,itime,sfname)
        sfname=sfname(1:13)//'_jtsli'
          open(unit=74,status='old',file=sfname,err=19979)
          close(unit=74,status='delete')
19979     continue
          open(unit=74,status='new',file=sfname                       )
        xdum=0.0
        do i=1,nw
          write (74,*) rjtli(i),sjtli(i),xdum,xdum
        enddo
        close(unit=74)
        call getfnmd('q',ishot,itime,sfname)
        sfname=sfname(1:13)//'_jtslir'
          open(unit=74,status='old',file=sfname,err=19983)
          close(unit=74,status='delete')
19983     continue
          open(unit=74,status='new',file=sfname                       )
        do i=1,nw
          write (74,*) rjtli(i),sjtlir(i),xdum,xdum
        enddo
        close(unit=74)
        call getfnmd('q',ishot,itime,sfname)
        sfname=sfname(1:13)//'_jtsliz'
          open(unit=74,status='old',file=sfname,err=19987)
          close(unit=74,status='delete')
19987     continue
          open(unit=74,status='new',file=sfname                       )
        do i=1,nw
          write (74,*) rjtli(i),sjtliz(i),xdum,xdum
        enddo
        close(unit=74)
      endif
      endif
!
      if (dco2v(iges,2).gt.1.e+10) then
      partic=dco2v(iges,2)*vout(iges)
      if (partic.gt.1.0e+10) tave(iges)=wplasm(iges)/partic/3. &
            /1.602e-16
      resist=vloopt(iges)/cpasma(iges)
      zeta=resist*areao(iges)/twopi/rout(iges)
      xlam=20.
      if (tave(iges).gt.0.001) then
         tevolt=sqrt((tave(iges)*1000.)**3)
      endif
      zeffr(iges)=zeta*tevolt/xlam/1.03e-02*2.
      if (iges.eq.igmax) then
      do 1100 m=1,igmax
        temp(m)=dco2v(m,2)
 1100 continue
         dttt=time(2)-time(1)
         if (igmax.eq.1) dttt=5.
         call getzeff(ishot,igmax,time(1),dttt,temp,tave,zeff,ier)
      endif
      endif
!------------------------------------------------------------------
!-- compute vessel forces                                        --
!------------------------------------------------------------------
      if (ifitvs.gt.0.or.icutfp.eq.2) then
        call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
        fztor=0.0
        fzpol=0.0
        do 50030 i=1,nvesel
          signz=zvs(i)/abs(zvs(i))
          signr=(rvs(i)-rcentr)/abs(rvs(i)-rcentr)
          if (avs(i).eq.0.0.and.avs2(i).eq.0.0) then
            if (wvs(i).lt.hvs(i)) then
              cosalp=0.0
              dells=hvs(i)
            else
              cosalp=1.0*signz
              dells=wvs(i)
            endif
          endif
          if (avs(i).ne.0.) then
            cosalp=abs(cosd(avs(i)))*signz
            dells=sqrt(hvs(i)**2+wvs(i)**2)
          endif
          if (avs2(i).ne.0.) then
            cosalp=abs(cosd(avs2(i)))*signz
            dells=sqrt(hvs(i)**2+wvs(i)**2)
          endif
          sinalp=signr*sqrt(1.-cosalp**2)
          rsnow0=rvs(i)-(wvs(i)/2.+wvs(i)/40.)*signz
          zsnow0=zvs(i)+(hvs(i)/2.+hvs(i)/40.)*signr
          sbrrs=0.0
          sumfzp=0.0
          do 50025 k=1,20
            rsnow=rsnow0+wvs(i)*k/20.*signz
            zsnow=zsnow0-hvs(i)*k/20.*signr
            call seva2d(bkx,lkx,bky,lky,c,rsnow,zsnow,pds,ier,n333)
            xpsivs=(simag-pds(1))/sidif
            sbrrs=sbrrs-pds(3)
            if (pds(1).lt.psibry.and. &
               (abs(zsnow).lt.abs(yvs2).or.zsnow*yvs2.lt.0.0)) then
               fpnow=ffcurr(xpsivs,kffcur)
            else
               fpnow=fbrdy
            endif
            delfp=fpnow -fbrdy
            btnow=fpnow*tmu/rsnow
            sumfzp=sumfzp+btnow*delfp
50025     continue
          vforcet(i)=-sbrrs/20.*twopi*vcurrt(i)
          fztor=fztor+vforcet(i)
          vforcep(i)= cosalp*sumfzp*dells/20.
          fzpol=fzpol+vforcep(i)
50030   continue
      endif
      if (icutfp.eq.2) then
        xxxx=1./xpsimin
        fpolvs=ffcurr(xxxx  ,kffcur)-ffcurr(x111,kffcur)
      endif
!
! --- prepare for vertical stability parameters
! --- calculation is done after pltout calls
!
      if (ivacum.le.0) then
         pleng=0.0
         do i=1,nfound-1
            ip1=i+1
            dli=sqrt((xout(ip1)-xout(i))**2+(yout(ip1)-yout(i))**2)
            pleng=pleng+dli
         enddo
         abar=100.*pleng/2./pi
      endif
!
      if (idebug /= 0) write (6,*) 'Call SHAPE/PLTOUT kerror = ', kerror
      if (itek.gt.0) then
	  if (idplace.ne.0) then
           call altplt(xout,yout,nfound,iges,nnn, &
            xmin,xmax,ymin,ymax,igmax)
          else
           if (idebug /= 0) write (6,*) 'Before SHAPE/PLTOUT'
           call pltout(xout,yout,nfound,iges,nnn, &
            xmin,xmax,ymin,ymax,igmax)
           if (idebug /= 0) write (6,*) 'After SHAPE/PLTOUT'
          endif
      endif
      if (idebug /= 0) write (6,*) 'exit SHAPE/PLTOUT kerror = ', kerror
      if (itrace.le.1) go to 1120
      if (abs(dpsi).gt.psitol) go to 1120
      if (itek.gt.0) then
         if (idplace.ne.0) then
           call altplt(xouts,youts,nfouns,jges,n22, &
            xmins,xmaxs,ymins,ymaxs,igmax)
         else
           call pltout(xouts,youts,nfouns,jges,n22, &
            xmins,xmaxs,ymins,ymaxs,igmax)
         endif
      endif
 1120 continue
      if ((itek.eq.2).and.(iges.eq.igmax)) call donepl
      if ((itek.eq.4).and.(iges.eq.igmax)) call donepl
      if ((itek.ge.5).and.(iges.eq.igmax)) call closepl
      if ((ilaser.gt.0).and.(iges.eq.igmax)) call donepl
      goto 1900
!
 1500 continue
      if (kerror.le.0) call chisqr(iges)
      if (itek.gt.0) then
        call pltout(xout,yout,nzz,iges,nnn,zxx,zxx,zxx,zxx,igmax)
      endif
      if ((itek.eq.2).and.(iges.eq.igmax)) call donepl
      if ((itek.eq.4).and.(iges.eq.igmax)) call donepl
      if ((itek.ge.5).and.(iges.eq.igmax)) call closepl
      if ((ilaser.gt.0).and.(iges.eq.igmax)) call donepl
!
 1900 continue
!
!-----------------------------------------------------------------------
!-- vertical stability parameter,  reference Nuc Fusion  18(1978)1331 --
!-- move out of pltout so that vertn, xnnc are indepent of itek value --
!-----------------------------------------------------------------------
      do 550 i=1,nw
         do 550 j=1,nh
            kk=(i-1)*nh+j
            copy(i,j)=0.0
            do m=1,nfcoil
               copy(i,j)=copy(i,j)+gridfc(kk,m)*brsp(m)
            enddo
            if (ivesel.le.0) go to 526
            do m=1,nvesel
               copy(i,j)=copy(i,j)+gridvs(kk,m)*vcurrt(m)
            enddo
  526       continue
            if (iecurr.le.0) go to 530
            do m=1,nesum
               copy(i,j)=copy(i,j)+gridec(kk,m)*ecurrt(m)
            enddo
  530       continue
  550 continue
      do i=1,nw
      do j=1,nh
        k=(i-1)*nh+j
        copyn(k)=copy(i,j)
      enddo
      enddo
      call sets2d(copyn,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
      rcur=rcurrt(jges)/100! m
      zcur=zcurrt(jges)/100! m
      call seva2d(bkx,lkx,bky,lky,c,rcur,zcur,pds,ier,n555)
      vertn(jges)=1.-pds(5)/pds(2)*rcur
      if (ivacum.le.0) then
         rx=rmagx(jges)/100.
         f_0=log(8*rout(jges)/abar)-2+betap(jges)+ali(jges)/2+.5
         delr=rout(jges)/100.-1.67
!-----------------------------------------------------------------------
!-- metal wall                                                        --
!-----------------------------------------------------------------------
         xnnc(jges)=vertn(jges)/((10.77*delr**2+8.08*delr+2.54)/f_0)
      endif
!
!------------------------------------------------------------------
!-- compute shearing rate eshear                                 --
!------------------------------------------------------------------
      if (keecur.le.0) goto 1990
      do i=1,nw
        eshear(i)=esradial(sipmid(i),keecur,rpmid(i),zmaxis)
      enddo
 1990 continue
!-----------------------------------------------------------------------
!--  write out r(shot).(time)_X files in rho space                    --
!-----------------------------------------------------------------------
      if (kwripre.eq.3) then
        call getfnmd('r',ishot,itime,sfname)
        sfname=sfname(1:13)//'_qpsi'
          open(unit=74,status='old',file=sfname,err=19920)
          close(unit=74,status='delete')
19920     continue
          open(unit=74,status='new',file=sfname                       )
        do i=1,nw
          write (74,*) rhovn(i),qpsi(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_jor'
          open(unit=74,status='old',file=sfname,err=19921)
          close(unit=74,status='delete')
19921     continue
          open(unit=74,status='new',file=sfname                       )
        do i=1,nw
          cjorka=cjor(i)/1000.
          write (74,*) rhovn(i),cjorka,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_jorec'
          open(unit=74,status='old',file=sfname,err=18921)
          close(unit=74,status='delete')
18921     continue
          open(unit=74,status='new',file=sfname                       )
        do i=1,nw
          cjorka=cjorec(i)/1000.
          write (74,*) rhovn(i),cjorka,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_cjmse'
          open(unit=74,status='old',file=sfname,err=18021)
          close(unit=74,status='delete')
18021     continue
          open(unit=74,status='new',file=sfname                       )
        if (ishot.le.97400) then
          mcentral=15
        else
          mcentral=10
        endif
        do i=1,mcentral+1
          cjorka=cjmse(i)/1000.
          write (74,*) rhogam(i),cjorka,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_cjmsec'
          open(unit=74,status='old',file=sfname,err=18023)
          close(unit=74,status='delete')
18023     continue
          open(unit=74,status='new',file=sfname                       )
        do i=1,mcentral+1
          cjorka=cjmsec(i)/1000.
          write (74,*) rhogam(i),cjorka,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_bzmse'
          open(unit=74,status='old',file=sfname,err=18025)
          close(unit=74,status='delete')
18025     continue
          open(unit=74,status='new',file=sfname                       )
        do i=1,nstark
          write (74,*) rhogam(i),bzmse(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_bzmsec'
          open(unit=74,status='old',file=sfname,err=18027)
          close(unit=74,status='delete')
18027     continue
          open(unit=74,status='new',file=sfname                       )
        do i=1,nstark
          write (74,*) rhogam(i),bzmsec(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_presf'
          open(unit=74,status='old',file=sfname,err=17922)
          close(unit=74,status='delete')
17922     continue
          open(unit=74,status='new',file=sfname                       )
        do i=1,nw
            write (74,*) rhovn(i),pres(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_presd'
          open(unit=74,status='old',file=sfname,err=18922)
          close(unit=74,status='delete')
18922     continue
          open(unit=74,status='new',file=sfname                       )
        do i=1,npress
            xmnow=-rpress(i)
            ymnow=seval(nw,xmnow,sigrid,rhovn,brhovn,crhovn,drhovn)
            write (74,*) ymnow,pressr(i),xdum,sigpre(i)
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_pbeam'
          open(unit=74,status='old',file=sfname,err=18932)
          close(unit=74,status='delete')
18932     continue
          open(unit=74,status='new',file=sfname                       )
        do i=1,nbeam
            xmnow=sibeam(i)
            ymnow=seval(nw,xmnow,sigrid,rhovn,brhovn,crhovn,drhovn)
            write (74,*) ymnow,pbeam(i),xdum,xdum
        enddo
        close(unit=74)
       if (keecur.gt.0) then
        sfname=sfname(1:13)//'_wexb'
          open(unit=74,status='old',file=sfname,err=19932)
          close(unit=74,status='delete')
19932     continue
          open(unit=74,status='new',file=sfname                       )
        do i=1,nw
          if (rpmid(i).ge.rmaxis) then
            write (74,*) rhopmid(i),eshear(i),xdum,xdum
          endif
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_errho'
          open(unit=74,status='old',file=sfname,err=19934)
          close(unit=74,status='delete')
19934     continue
          open(unit=74,status='new',file=sfname                       )
        do i=1,nw
          if (rpmid(i).ge.rmaxis) then
            write (74,*) rhopmid(i),ermid(i),xdum,xdum
          endif
        enddo
        close(unit=74)
       endif
      endif
!
      DEALLOCATE(xsisii,bpres,cpres,dpres,sjtli,sjtlir,sjtliz, &
            rjtli,bpresw,cpresw,dpresw,copyn,cjtli,x,y)
!
      return
 2000 format (1x,'NBDRY = ',i5,/,1x,'RBDRY = ',/)
 2010 format (1x,'ZBDRY = ',/)
 2020 format (6(1x,e12.5))
 2980 format (1x,i6)
 3000 format (1x,6e12.5)
 3020 format (1x,i6,2e12.5,2i7)
      end
      subroutine splitc(is,rs,zs,cs,rc,zc,wc,hc,ac,ac2,cc)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**                                                                  **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      data pi/3.1415926535897932/
      dimension rs(1),zs(1),cs(1)
!
      frd=pi/180.
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
      side=tan(frd*ac)*wc
      hdelt=hc/is
      wdelt=wc/is
      zdelt=tan(frd*ac)*wdelt
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
      side=hc/tan(frd*ac2)
      hdelt=hc/is
      wdelt=wc/is
      zstrt=zc-hc/2.+hdelt/2.
      rdelt=hdelt/tan(frd*ac2)
      rstrt=rc-side/2.-wc/2.+rdelt/2.+wdelt/2.
      side=hc/tan(frd*ac2)
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
      end
      subroutine steps(ix,ixt,ixout,jtime,kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          step computes the dimensionless poloidal fluxes for     **
!**          the r-z grid.                                           **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky,zeros,xouts,youts, &
           rsplt,zsplt,csplt
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      common/cwork3/lkx,lky
      common/wwork1/xlims(5),ylims(5),limtrs,xlmins
      dimension pds(6)
      integer iii
      real :: zmaxis_last = 0.0
      data isplit/8/,psitol/1.0e-04/
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
      call findax(nw,nh,rgrid,zgrid,rmaxis,zmaxis,simag, &
                  psibry,rseps(1,jtime),zseps(1,jtime),m10, &
                  xout,yout,nfound,psi,xmin,xmax,ymin,ymax, &
                  zxmin,zxmax,rymin,rymax,dpsi,bpol,bpolz, &
                  limitr,xlim,ylim,limfag)
      if (nsol.gt.0) then
        if (idebug >= 2) then
          write (6,*) 'STEPS R,Z,Si,Err = ', rsol(1),zsol(1),wsisol,ier
          call seva2d(bkx,lkx,bky,lky,c,rbdry(1),zbdry(1),pds,ier,n111)
          write (6,*) 'STEPS R,Z,Si,Err = ', rbdry(1),zbdry(1),pds(1),ier
          call seva2d(bkx,lkx,bky,lky,c,rbdry(nbdry),zbdry(nbdry),pds,ier &
             ,n111)
          write (6,*) 'STEPS R,Z,Si,Err = ', rbdry(nbdry),zbdry(nbdry) &
            ,pds(1),ier
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
      if (pasmat(jtime).lt.-1.e3) then
        nnerr=10000
      else
        nnerr=0
      endif
      call bound(psi,nw,nh,nwnh,psibry,xmin,xmax,ymin,ymax, &
                 zero,rgrid,zgrid,xguess,yguess,ixt,limitr,xlim,ylim, &
                 xout,yout,nfound,xltrac,npoint,rymin,rymax,dpsi, &
                 zxmin,zxmax,nnerr,ishot,itime, &
                 limfag,radbou,kbound)
      if (nnerr.ge.1) then
        kerror=1
        write(nttyo,2100) ishot,itime
        open(unit=40,file='errfil.out',status='unknown',access='append' &
                                 )
        write(40,2100) ishot,itime
        close(unit=40)
        return
      endif
  155 continue
!----------------------------------------------------------------------
!--  find magnetic axis and poloidal flux at axis simag              --
!----------------------------------------------------------------------
      m20=20
      call findax(nw,nh,rgrid,zgrid,rmaxis,zmaxis,simag, &
                  psibry ,rseps(1,jtime),zseps(1,jtime),m20, &
                  xout,yout,nfound,psi,xmin,xmax,ymin,ymax, &
                  zxmin,zxmax,rymin,rymax,dpsi,bpol,bpolz, &
                  limitr,xlim,ylim,limfag)
      sidif=simag-psibry
      eouter=(ymax-ymin)/(xmax-xmin)
      zplasm=(ymin+ymax)/2.
      aouter=(xmax-xmin)/2.
!-----------------------------------------------------------------------
!--   force free current in the scrape-off layer                      --
!-----------------------------------------------------------------------
      if (icutfp.eq.2) then
        xvsmaxo=xvsmax
        xvsmin=1.e10
        xvsmax=-1.e10
        if (limvs.eq.0) then
        itot=isplit*isplit
        do 51000 k=1,nvesel
          call splitc(isplit,rsplt,zsplt,csplt, &
                  rvs(k),zvs(k),wvs(k),hvs(k),avs(k),avs2(k),cdum)
          do 50900 kk=2,itot
            call seva2d(bkx,lkx,bky,lky,c,rsplt(kk),zsplt(kk), &
                        pds,ier,n111)
            xvsmin=min(xvsmin,pds(1))
            xvsmax=max(xvsmax,pds(1))
50900     continue
51000   continue
        else
        do 51009 k=1,limitr-1
           delx=xlim(k+1)-xlim(k)
           dely=ylim(k+1)-ylim(k)
           dels=sqrt(delx**2+dely**2)
           nn=dels/0.002
           nn=max(5,nn)
           delx=delx/float(nn-1)
           dely=dely/float(nn-1)
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
          if (zero(kk).gt.0.0005.and.www(kk).lt.0.1) then
           call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n333)
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
        if (bpmin.le.0.10*avebp.and.relsi.gt.0.005) then
!------------------------------------------------------------------
!-- find second separatrix                                       --
!------------------------------------------------------------------
          do 9300 j=1,40
          call seva2d(bkx,lkx,bky,lky,c,xs,ys,pds,ier,n666)
          det=pds(5)*pds(6)-pds(4)*pds(4)
          if (abs(det).lt.1.0e-15) go to 9305
          xerr=(-pds(2)*pds(6)+pds(4)*pds(3))/det
          yerr=(-pds(5)*pds(3)+pds(2)*pds(4))/det
          xs=xs+xerr
          ys=ys+yerr
          if (xerr*xerr+yerr*yerr.lt.1.0e-12) go to 9310
 9300     continue
 9305     continue
          epssep=xerr*xerr+yerr*yerr
          write (nttyo,11001) epssep,ixt
          if (iand(iout,1).ne.0) write (nout,11001) epssep,ixt
          if (epssep.lt.1.0e-10) go to 9310
          go to 9320
 9310     continue
         sibpmin=pds(1)
         yvs2=ys
         rsepex=xs
         relsi=abs((sibpmin-psibry)/sidif)
         if (relsi.gt.0.005) then
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
        xpsi(kk)=1.1
        if ((rgrid(i).lt.xmin).or.(rgrid(i).gt.xmax)) go to 1000
        if ((zgrid(j).lt.ymin).or.(zgrid(j).gt.ymax)) go to 1000
        xpsi(kk)=(simag-psi(kk))/sidif
        else
        if (zero(kk).gt.0.0005) then
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
        write (6,*) 'STEPS lkx,lky,c = ',bkx(33),bky(33),c(1,33,1,33)

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
        if (rrgam(jtime,k).le.0.0) goto 50299
        call seva2d(bkx,lkx,bky,lky,c,rrgam(jtime,k) &
                    ,zzgam(jtime,k),pds,ier,n111)
        sistark(k)=pds(1)
        sisinow=(simag-pds(1))/sidif
        sigam(k)=sisinow
        fpnow=ffcurr(sisinow,kffcur)
        btgam(k)=fpnow*tmu/rrgam(jtime,k)
50299 continue
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
      if (abs(sizeroj(1)-1.0).le.1.e-05.and.kzeroj.eq.1) goto 51977
      if (kzeroj.gt.0) then
       do i=1,kzeroj
       if (sizeroj(i).ge.1.0) sizeroj(i)=0.99999
       siwant=simag+sizeroj(i)*(psibry-simag)
       call surfac(siwant,psi,nw,nh,rgrid,zgrid,rsplt,zsplt,nfounc &
                    ,npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nnn, &
                    rmaxis,zmaxis,negcur)
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
       call surfac(siwant,psi,nw,nh,rgrid,zgrid,rsplt,zsplt,nfounc &
                    ,npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nnn, &
                    rmaxis,zmaxis,negcur)
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
      cvolp(ixt)=abs(cvolp(ixt))*1.0e+06
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
	elseif (eelip.gt.2.25.and.itell.eq.0) then
           delzmm = zmaxis - zmaxis_last
           zmaxis_last=zmaxis
      endif
!----------------------------------------------------------------------
!--  magnetic axis parameters if needed                              --
!----------------------------------------------------------------------
      if (icinit.gt.0) then
        if ((iconvr.ne.3).and.(ixout.le.1)) go to 1580
      endif
      go to (1570,1590,1595,1580,1590) icurrt
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
      if (errorm.lt.0.1000.and.icurrt.eq.4) nqend=nqiter
      do 1585 i=1,nqend
      if (i.gt.1) call currnt(n22,jtime,n22,n22,kerror)
! MPI >>>
#if defined(USEMPI)
      !if (kerror /= 0) then
      !  kerror = 1
      !  return
      !endif
#endif
! MPI <<<
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
 1587 continue
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
 2100 format(/,1x,4hshot,i6,4h at ,i6,4h ms ,'** Problem in BOUND **')
11001 format(/,1x,'** 2nd seperatrix **',2x,e10.3,2x,i4)
      end
      subroutine tsorder(mbox,zprof,nemprof,temprof,nerprof,terprof)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          this subroutine reorders the z profile data to be       **
!**          in ascending order and sets the ne and te data to       **
!**          correspond to the new order.                            **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          27/11/90..........first created, T. Carlstrom           **
!**          08/07/91..........revised for EFIT, L. Lao              **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      implicit integer*4 (i-n), real*8 (a-h, o-z)
	dimension ztemp(40),temtemp(40),tertemp(40),zinc(40), &
                  zprof(1), temprof(1), terprof(1)
	real*8 nemtemp(40),nertemp(40),nemprof(1),nerprof(1)
!---------------------------------------------------------------------
!--	copy zprof to ztemp (temporary work space)                  --
!---------------------------------------------------------------------
	do 100 j=1,mbox
          ztemp(j)=zprof(j)
  100   continue
!---------------------------------------------------------------------
!--	find min z in ztemp                                         --
!---------------------------------------------------------------------
	do 1000 j=1,mbox
		zmin=999.
			do 800 i=1,mbox-j+1
				if(ztemp(i).lt.zmin) then
                                  zmin=ztemp(i)
                                  kmin=i
                                end if
  800			continue
!---------------------------------------------------------------------
!--	put zmin into new vectors                                   --
!---------------------------------------------------------------------
		zinc(j)=zmin
                nemtemp(j)=nemprof(kmin)
                temtemp(j)=temprof(kmin)
                nertemp(j)=nerprof(kmin)
                tertemp(j)=terprof(kmin)
!---------------------------------------------------------------------
!--	create new ztemp with remaining data                        --
!---------------------------------------------------------------------
		k=0
		do 900 i=1,mbox-j+1
			if(zmin.ne.ztemp(i))then
				k=k+1
				ztemp(k)=ztemp(i)
			end if
  900 		continue
 1000 	continue
!---------------------------------------------------------------------
!--	rewrite new vectors                                         --
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
	end
      subroutine vescur(jtime)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          vescur computes the currents induced in the vessel      **
!**          segments due to E coils and F coils.                    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/11/85..........first created                         **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
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
      end
      subroutine weight(x,y)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          weight computes the weighting function w.               **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
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
      go to (10,11,12,13,14),in
   10 www(kk) = 1.
      go to 20
   11 xxw = (p1+p2+p3+p4)/4.
      yyw = (a1+a2+a3+a4)
      yyw = (yyw/(xxw-yyw))**2
      www(kk) = 1.-0.5*yyw
      go to 20
   12 xxw = (p1+p2+p3+p4)
      yyw = (a1+a2+a3+a4)
      www(kk) = xxw/(xxw-yyw)
      go to 20
   13 xxw = (p1+p2+p3+p4)
      yyw = (a1+a2+a3+a4)/4.
      xxw = (xxw/(xxw-yyw))**2
      www(kk) = 0.5*xxw
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
      end
      subroutine old_new(psi,nw,nh,nwh,psivl,xmin,xmax,ymin,ymax, &
           zero,x,y,xctr,yctr,ix,limitr,xlim,ylim,xcontr,ycontr, &
           ncontr,xlmin,npoint,rymin,rymax,dpsi,zxmin,zxmax,nerr, &
           ishot,itime,limfag,radold,kbound)
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      return
      end
      block data efit_bdata
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          block data routine to hold all data statements for      **
!**          variables that are in common blocks.                    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**                91..........first created, John Ferron            **
!**          11/11/93..........revised                               **
!**          10/01/97..........revised by Q.Peng,.ddd specific data  **
!**                            are moved to a different file         **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: zeros,xouts,youts,bpoo,bpooz,bpooc
      include 'eparmdud129.f90'
      include 'modules1.f90'
      include 'env2d.inc'
      common/wwork1/xlims(5),ylims(5),limtrs,xlmins
      data limtrs/5/
      data iunit/35/, m_write/1/, m_read/1/
      data out2d /'curve2d.dat'/, out2d_bin/'curve2d_bin.dat'/
!
      end
!----------------------------------------------------------------------
!--  Subroutine for writing K file in K file mode Qilong Ren         -- 
!----------------------------------------------------------------------
      subroutine write_K(ksstime,kerror)
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
      include 'basiscomdu.f90'
      
! MPI >>>
#if defined(USEMPI)
      include 'mpif.h'
#endif
! MPI <<<
      
      logical lopened
      character filenm*15,ishotime*12,news*72,table_s*72, &
                eqdsk*30,comfile*15,prefix1*1,header*42,fit_type*3
      integer*4 ltbdis
      dimension coils(nsilop),expmp2(magpri), &
                denr(nco2r),denv(nco2v), &
                tgamma(nmtark),sgamma(nmtark),rrrgam(nmtark), &
                zzzgam(nmtark),aa1gam(nmtark),aa2gam(nmtark), &
                aa3gam(nmtark),aa4gam(nmtark),aa5gam(nmtark), &
                aa6gam(nmtark)
      dimension tlibim(libim),slibim(libim),rrrlib(libim)
      character*82 snap_ext
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
           ,iplcout,imagsigma,errmag,saimin,errmagb,fitfcsum,fwtfcsum,efitversion
      namelist/inwant/psiwant,vzeroj,nccoil,currc79,currc139,rexpan, &
           znose,sizeroj,fitdelz,relaxdz,errdelz,oldccomp,nicoil, &
           oldcomp,currc199,curriu30,curriu90, &
           curriu150,curril30,curril90,curril150,ifitdelz,scaledz
      namelist/ink/isetfb,ioffr,ioffz,ishiftz,gain,gainp,idplace &
           ,symmetrize,backaverage
      namelist/ins/tgamma,sgamma,fwtgam,rrrgam,zzzgam,aa1gam,aa2gam, &
                   aa3gam,aa4gam,aa5gam,aa6gam,aa7gam,msebkp, &
            msefitfun,mse_quiet,mse_spave_on,kwaitmse, &
            dtmsefull,mse_strict,t_max_beam_off,ok_30rt,ok_210lt
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
           ,ok_30rt,ok_210lt,vbit,nbdrymx
      namelist/efitink/isetfb,ioffr,ioffz,ishiftz,gain,gainp,idplace &
           ,symmetrize,backaverage,lring
      data mcontr/35/,lfile/36/,ifpsi/0/
      data currn1/0.0/,currc79/0.0/,currc139/0.0/,currc199/0.0/ &
                      ,curriu30/0.0/,curriu90/0.0/,curriu150/0.0/ &
                      ,curril30/0.0/,curril90/0.0/,curril150/0.0/
      integer, intent(out) :: ksstime
      integer, intent(inout) :: kerror
      kerror = 0
!---------------------------------------------------------------------
!--  generate input files and command file for running EFITD        --
!---------------------------------------------------------------------
!      if (kdata .ge. 5 .and. kdata .lt. 7) then
      if ((kdata .ge. 5 .and. kdata .lt. 7).or.(kdata.eq.8)) then
 3020 continue
       if (rank == 0) then
        if (use_opt_input .eqv. .false.) then
          write (nttyo,6610)
          read (ntty,6620) comfile
          write (nttyo,6600)
          read (ntty,6620) filenm
          if (filenm.eq.'0') then
            write (nttyo,6040)
            read (ntty,*) lshot,timeb,dtime,ktime
          endif
        else
          comfile = cmdfile_in
          filenm = shotfile_in
          if (filenm .eq. '0') then
            lshot = shot_in
            timeb = starttime_in
            dtime = deltatime_in
            ktime = steps_in
          endif
        endif
       endif
! MPI >>>
#if defined(USEMPI)
       if (nproc > 1) then
        call MPI_BCAST(comfile,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(filenm,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
        if (filenm == '0') then
          call MPI_BCAST(lshot,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(ktime,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(timeb,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(dtime,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        endif
       endif
#endif
! MPI <<<
       if (comfile.ne.'0') then
         open(unit=nffile,status='old',file=comfile,err=12928)
         close(unit=nffile,status='delete')
12928    continue
         open(unit=nffile,status='new',file=comfile)
         write (nffile,4958)
       endif
       if (filenm.ne.'0') &
          open(unit=nout,status='old',file=filenm &
                                  ,err=3025         )
       go to 3028
 3025  open(unit=nout,status='old', &
            file='phys_data:[d3phys.diiid.gsl]'//filenm &
           ,err=3020)
      endif
 3028 continue
      open(unit=neqdsk,status='old', &
           file='efit_snap.dat',err=3030)
      snapfile='efit_snap.dat'
      go to 3032
 3030 continue
      open(unit=neqdsk,status='old', &
           file= input_dir(1:lindir)//'efit_snap.dat'         )
      snapfile=input_dir(1:lindir)//'efit_snap.dat'
 3032 continue
      do 3034 i=1,nfcoil
        fwtfc(i)=0.
 3034 continue
      do i=1,nesum
        fwtec(i)=0.
      enddo
      do i=1,mbdry
       fwtbdry(i)=1.0
       fwtsol(i)=1.0
       sigrbd(i)=1.e10
       sigzbd(i)=1.e10
      enddo
      do 3036 i=1,magpri
        fwtmp2(i)=0.
        bitmpi(i)=0.0
 3036 continue
      do 3037 i=1,nstark
        fwtgam(i)=0.0
 3037 continue
      do 3038 i=1,nsilop
        fwtsi(i)=0.
        psibit(i)=0.0
 3038 continue
      error=1.e-03
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
      scrape=0.030
      qvfit=0.95
      serror=0.03
      xltype=0.
      xltype_180=0.0
!
      read (neqdsk,efitin,end=111)
 111  continue
      read (neqdsk,efitink,err=3039,end=112)
 112  continue
 3039 close(unit=neqdsk)
      if (imagsigma.eq.1) then
         do_spline_fit=.false.
         saimin=150.
      endif
!----------------------------------------------------------------------
!-- recalculate length of default directories in case any change     --
!----------------------------------------------------------------------
      ltbdir=0
      lindir=0
      lstdir=0
      do i=1,len(table_dir)
         if (table_dir(i:i).ne.' ') ltbdir=ltbdir+1
         if (input_dir(i:i).ne.' ') lindir=lindir+1
         if (store_dir(i:i).ne.' ') lstdir=lstdir+1
      enddo
      if (lshot.ge.112000.and.jtime.le.1) then
        if (lshot.lt.156000) then
          table_di2 = table_dir(1:ltbdir)//'112000/'
        else
          if (kdata.ne.2) then
            table_di2 = table_dir(1:ltbdir)//'156014/'
          else
            if (efitversion <= 20140331) then
               table_di2 = table_dir(1:ltbdir)//'112000/'
            else
               table_di2 = table_dir(1:ltbdir)//'156014/'
            endif
          endif
        endif
        ltbdi2=ltbdir+7
       endif
       ksstime = ktime
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
        calpa(1,1)=0.1
        calpa(2,1)=0.1
        calpa(3,1)=0.1
        xalpa(1)=0.0
        kcgama=1
        cgama(1,1)=0.1
        cgama(2,1)=0.1
        cgama(3,1)=0.1
        xgama(1)=0.0
      endif
!
      limitr=-limid
      do 3040 i=1,nfcoil
        swtfc(i)=fwtfc(i)
 3040 continue
      do i=1,nesum
        swtec(i)=fwtec(i)
      enddo
      do 3042 i=1,magpri
        swtmp2(i)=fwtmp2(i)
 3042 continue
      do 3044 i=1,nsilop
        swtsi(i)=fwtsi(i)
 3044 continue
      swtcur=fwtcur
      do 3045 i=1,nstark
        swtgam(i)=fwtgam(i)
 3045 continue
 3046 continue
!
      if (filenm.ne.'0') then
        read (nout,4970,end=3200) ishotime
        read (ishotime,fmt='(i6,1x,i5)',err=3046) ishot,itime
        times=itime/1000.
        delt=0.002
        ktime=1
        timeb=itime
        dtime=0.
      else
        times=timeb/1000.
        delt=dtime/1000.
        ishot=lshot
      endif
      if (ifitvs.gt.0) then
          istop=-1
      else
          istop=0
      endif
      if (ishot.gt.108281) n1coil = 0
!----------------------------------------------------------------------
!-- Fetch data                                                       --
!----------------------------------------------------------------------
! MPI >>>
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
          if (filenm.ne.'0') then
            go to 3046
          else
! MPI >>>
#if defined(USEMPI)
            call mpi_stop
#else
            stop
#endif
! MPI <<<
          endif
      endif
!
      mmstark=0
      do 3056 i=1,nstark
         if (swtgam(i).gt.1.e-06) mmstark=mmstark+1
 3056 continue
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
!----------------------------------------------------------------------
! --- get xltype, xltype_180 without reading limiters                --
!----------------------------------------------------------------------
      call getlim(0,xltype,xltype_180)
!----------------------------------------------------------------------
!-- adjust fitting weights based on FITWEIGHT.DAT                    --
!----------------------------------------------------------------------
      if (lookfw.ge.0) then
        do 3048 i=1,magpri
        rwtmp2(i)=0.0
 3048   continue
        do i=1,nsilop
         rwtsi(i)=0.0
        enddo
        open(unit=neqdsk,status='old', &
             file=table_di2(1:ltbdi2)//'fitweight.dat'         )
 3050   read (neqdsk,*,end=3052) irshot
        if (irshot.gt.ishot) go to 3052
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
 3052   continue
        close(unit=neqdsk)
      endif
!
      do 3060 i=1,ktime
        time(i)=time(i)*1000.
 3060 continue
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
!-- Get edge pedestal tanh paramters                                  --
!-----------------------------------------------------------------------
      if (fitzts.eq.'te') then
        call gettanh(ishot,fitzts,ktime,time,ztssym,ztswid, &
                          ptssym,ztserr)
      endif
!
      do 3190 jtime=1,ktime
        itime=time(jtime)
        timems=itime
        timeus=(time(jtime)-timems)*1000.
        itimeu=timeus
!-----------------------------------------------------------------------
!--  correction for truncation                                        --
!-----------------------------------------------------------------------
        if (itimeu.ge.990) then
          itime=itime+1
          itimeu=0
          time(jtime)=itime
        endif
        do 3068 i=1,nsilop
          coils(i)=silopt(jtime,i)
          fwtsi(i)=swtsi(i)
          if (lookfw.gt.0.and.fwtsi(i).gt.0.0) fwtsi(i)=rwtsi(i)
          if (ierpsi(i).ne.0) fwtsi(i)=0.0
 3068   continue
        do 3070 i=1,magpri
          expmp2(i)=expmpi(jtime,i)
          fwtmp2(i)=swtmp2(i)
          if (lookfw.gt.0.and.fwtmp2(i).gt.0.0) fwtmp2(i)=rwtmp2(i)
          if (iermpi(i).ne.0) fwtmp2(i)=0.0
 3070   continue
        do 3072 i=1,nfcoil
          brsp(i)=fccurt(jtime,i)
          fwtfc(i)=swtfc(i)
          if (ierfc(i).ne.0) fwtfc(i)=0.0
 3072   continue
        do 3074 i=1,nesum
          ecurrt(i)=eccurt(jtime,i)
          fwtec(i)=swtec(i)
          if (ierec(i).ne.0) fwtec(i)=0.0
 3074   continue
        do 3080 i=1,nco2r
          denr(i)=denrt(jtime,i)
 3080   continue
        do 3085 i=1,nco2v
          denv(i)=denvt(jtime,i)
 3085   continue
        if (mmstark.gt.0) then
         do 3090 i=1,nmtark
          tgamma(i)=tangam(jtime,i)
          sgamma(i)=siggam(jtime,i)
          rrrgam(i)=rrgam(jtime,i)
          zzzgam(i)=zzgam(jtime,i)
          aa1gam(i)=a1gam(jtime,i)
          aa2gam(i)=a2gam(jtime,i)
          aa3gam(i)=a3gam(jtime,i)
          aa4gam(i)=a4gam(jtime,i)
          aa5gam(i)=a5gam(jtime,i)
          aa6gam(i)=a6gam(jtime,i)
          fwtgam(i)=swtgam(i)
          if (iergam(i).ne.0) fwtgam(i)=0.0
 3090    continue
        endif
        fwtcur=swtcur
        if (ierpla.ne.0) fwtcur=0.0
        btor=bcentr(jtime)
        plasma=pasmat(jtime)
        siref=psiref(jtime)
        vloop=vloopt(jtime)
        dflux=1.0e+03*diamag(jtime)
        sigdlc=1.0e+03*sigdia(jtime)
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
!-- Set edge pedestal tanh paramters                                  --
!-----------------------------------------------------------------------
        if (fitzts.eq.'te'.and.ztserr(jtime)) then
          nbdry=1
          rbdry(1)=1.94
          zbdry(1)=ztssym(jtime)+0.5*ztswid(jtime)
        endif
 3180   continue

        if (kdata .ge. 5 .and. kdata .lt. 7) then
          call getfnmu(itimeu,'k',ishot,itime,eqdsk)
        endif
        if (kdata.eq.8) then
          call getfnmu_pefit(itimeu,'k',ishot,itime,eqdsk)
        endif
        open(unit=neqdsk,status='old',file=eqdsk,err=12915, &
           recl=72,delim='APOSTROPHE')
        close(unit=neqdsk,status='delete')
12915   continue
!-------------------------------------------------------------------------------
!--  Set bit noise for ishot > 152000                                         --
!-------------------------------------------------------------------------------
        if (ishot.gt.152000) vbit = 80
        open(unit=neqdsk,                       file=eqdsk,status='new', &
           recl=72,delim='APOSTROPHE')
        write (neqdsk,in1)
        write (neqdsk,inwant)
        if (isetfb.ne.0) write (neqdsk,ink)
        if (mmstark.gt.0) write (neqdsk,ins)
        if (kwaitmse.ne.0) write (neqdsk,ina)
        if (kfitece.gt.0) write (neqdsk,inece)
        if (keecur.gt.0) write (neqdsk,iner)
!-----------------------------------------------------------------------
!--  fitting type flag                                                --
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
!---------------------------------------------------------------------
!-- Append SNAP file                                                --
!---------------------------------------------------------------------
        if (appendsnap.eq.'K'.or.appendsnap.eq.'KG') then
        if (snapfile/='none') then
          open(unit=nsnapf,status='old', &
            file=snapfile,err=9981)
        do i=1,1000000
          read (nsnapf,9991,end=9981) tmpdata
          if (INDEX(tmpdata,'&efitin')/=0) go to 381
        enddo
 381    continue
        do i=1,1000000
           write (neqdsk,9991) tmpdata
           read (nsnapf,9991,end=9981) tmpdata
           if (INDEX(tmpdata,'/')/=0) then
             write (neqdsk,9991) tmpdata
             go to 9979
           endif
        enddo
 9979   close (unit=nsnapf)
 9981   continue
 9991   format (a)
        endif
        endif
!
        close(unit=neqdsk)
        if (filenm.eq.'0') then
          read (eqdsk,6700) prefix1,ishotime
        endif
        if (comfile.ne.'0') &
          write (nffile,4960) ishotime
 3190 continue
      if (filenm.ne.'0') go to 3046
 3200 continue
      if (comfile.ne.'0') then
        write (nffile,4962)
        close(unit=nffile)
      endif
      if (filenm.ne.'0') close(unit=nout)
 4042 format (1x,a42,1x,a3)
 4958 format ('#!/bin/csh -f')
 4960 format ('      runefit.sc k',a12)
 4962 format ('#',/,'exit')
 4970 format (2x,a,1x,a,1x,a,1x,a)
 4980 format (i5)
 5000 format (2e12.6)
!vas 5500 format (/,10x,'EFITD 129dx2 Version  ',2a5,/)
 5500 format (/,10x,'EFITD ',a3,' x ',a3,' Version  ',2a5,/)
 6000 format (/,1x,'type mode (2=file, 3=snap, 4=time', &
               ', 5=input, 6=com file, 7=snap_ext):')
 6040 format (/,1x,'type shot #, start time(ms), time step(ms), steps' &
        ,'(<401):')
 6080 format (/,1x,'type limiter position (cm, 0=ptdata):')
 6090 format(' enter number of extra field lines to trace:')
 6091 format(' enter scrape off depth(m),'/ &
       '       sense of tracing (+1 for down, -1 for up),'/ &
       '       ixstrt (+1 for start on outside, -1' &
       ' for start inside):')
 6100 format(/,1x,48htype plot mode (0=none, 1=tektronix, 2=versatec, &
           ,17h 3=qms, -=x ray):)
 6200 format (/,1x,22hnumber of time slices?)
 6220 format (/,1x,22htype input file names:)
 6230 format (1x,1h#)
 6240 format (a)
 6600 format (/,1x,'good shot list file name ( 0=tty) ?')
 6610 format (/,1x,'command file name ( 0=none) ?')
 6617 format (/,1x,'type snap file extension (def for default):')
 6620 format (a)
 6700 format (a1,a12)
20000 format (/,1x,'shot data not on disk')
30000 format (i9)
30200 format (10f3.0)
      end

!----------------------------------------------------------------------
!--  Subroutine for writing K file in SNAP mode 2014/04/16           -- 
!----------------------------------------------------------------------
      subroutine write_K2(jtime,kerror)
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
      include 'basiscomdu.f90'
      
! MPI >>>
#if defined(USEMPI)
      include 'mpif.h'
#endif
! MPI <<<
      
      logical lopened
      character filenm*15,ishotime*12,news*72,table_s*72, &
                eqdsk*20,comfile*15,prefix1*1,header*42,fit_type*3
      integer*4 ltbdis,limitrss
      dimension coils(nsilop),expmp2(magpri), &
                denr(nco2r),denv(nco2v), &
                tgamma(nmtark),sgamma(nmtark),rrrgam(nmtark), &
                zzzgam(nmtark),aa1gam(nmtark),aa2gam(nmtark), &
                aa3gam(nmtark),aa4gam(nmtark),aa5gam(nmtark), &
                aa6gam(nmtark)
      dimension tlibim(libim),slibim(libim),rrrlib(libim)
      character*82 snap_ext
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
           ,iplcout,imagsigma,errmag,saimin,errmagb,fitfcsum,fwtfcsum,efitversion
      namelist/inwant/psiwant,vzeroj,nccoil,currc79,currc139,rexpan, &
           znose,sizeroj,fitdelz,relaxdz,errdelz,oldccomp,nicoil, &
           oldcomp,currc199,curriu30,curriu90, &
           curriu150,curril30,curril90,curril150,ifitdelz,scaledz
      namelist/ink/isetfb,ioffr,ioffz,ishiftz,gain,gainp,idplace &
           ,symmetrize,backaverage
      namelist/ins/tgamma,sgamma,fwtgam,rrrgam,zzzgam,aa1gam,aa2gam, &
                   aa3gam,aa4gam,aa5gam,aa6gam,aa7gam,msebkp, &
            msefitfun,mse_quiet,mse_spave_on,kwaitmse, &
            dtmsefull,mse_strict,t_max_beam_off,ok_30rt,ok_210lt
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
           ,ok_30rt,ok_210lt,vbit,nbdrymx
      namelist/efitink/isetfb,ioffr,ioffz,ishiftz,gain,gainp,idplace &
           ,symmetrize,backaverage,lring
      data mcontr/35/,lfile/36/,ifpsi/0/
      data currn1/0.0/,currc79/0.0/,currc139/0.0/,currc199/0.0/ &
                      ,curriu30/0.0/,curriu90/0.0/,curriu150/0.0/ &
                      ,curril30/0.0/,curril90/0.0/,curril150/0.0/
      integer, intent(in) :: jtime
      integer, intent(inout) :: kerror
      kerror = 0
!----------------------------------------------------------------
!-- recover the value of table_dir for mode 3 or 7             --
!----------------------------------------------------------------
      if (ishot.ge.112000) then
        ltbdis =ltbdir
        ltbdir = ltbdir-7
        table_s = table_dir
        table_dir = table_dir(1:ltbdir)
      endif
 3028 continue
!
      itime=time(jtime)
      timems=itime
      timeus=(time(jtime)-timems)*1000.
      itimeu=timeus
      siref=psirefs(jtime)
      do 3068 i=1,nsilop
        coils(i)=silopt(jtime,i)-siref
        fwtsi(i)=swtsi(i)
 3068 continue
      do 3070 i=1,magpri
        expmp2(i)=expmpi(jtime,i)
        fwtmp2(i)=swtmp2(i)
 3070 continue
      brspss=brsp
      brsp=0.0
      do 3072 i=1,nfcoil
        brsp(i)=fccurt(jtime,i)
        fwtfc(i)=swtfc(i)
 3072 continue
      do 3074 i=1,nesum
        ecurrt(i)=eccurt(jtime,i)
        fwtec(i)=swtec(i)
 3074 continue
      do 3080 i=1,nco2r
        denr(i)=denrt(jtime,i)
 3080 continue
      do 3085 i=1,nco2v
        denv(i)=denvt(jtime,i)
 3085 continue
      if (mmstark.gt.0) then
      do 3090 i=1,nmtark
        tgamma(i)=tangam(jtime,i)
        sgamma(i)=siggam(jtime,i)
        rrrgam(i)=rrgam(jtime,i)
        zzzgam(i)=zzgam(jtime,i)
        aa1gam(i)=a1gam(jtime,i)
        aa2gam(i)=a2gam(jtime,i)
        aa3gam(i)=a3gam(jtime,i)
        aa4gam(i)=a4gam(jtime,i)
        aa5gam(i)=a5gam(jtime,i)
        aa6gam(i)=a6gam(jtime,i)
        fwtgam(i)=swtgam(i)
 3090 continue
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
        rbdry(1)=1.94
        zbdry(1)=ztssym(jtime)+0.5*ztswid(jtime)
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
      dflux=1.0e+03*diamag(jtime)
      sigdlc=1.0e+03*sigdia(jtime)
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
 3180 continue
      call getfnmu(itimeu,'k',ishot,itime,eqdsk)
      open(unit=neqdsk,status='old',file=eqdsk,err=12915, &
           recl=72,delim='APOSTROPHE')
      close(unit=neqdsk,status='delete')
12915 continue
!-------------------------------------------------------------------------------
!--  Write K file                                                             --
!-------------------------------------------------------------------------------
      open(unit=neqdsk,                       file=eqdsk,status='new', &
           recl=72,delim='APOSTROPHE')
      write (neqdsk,in1)
      write (neqdsk,inwant)
      if (isetfb.ne.0) write (neqdsk,ink)
      if (mmstark.gt.0) write (neqdsk,ins)
      if (kwaitmse.ne.0) write (neqdsk,ina)
      if (kfitece.gt.0) write (neqdsk,inece)
      if (keecur.gt.0) write (neqdsk,iner)
!-----------------------------------------------------------------------
!--  fitting type flag                                                --
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
!-- Restore variables                                                --
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
!-- Append SNAP file                                                --
!---------------------------------------------------------------------
      if (kdata .eq.  7) then
         open(unit=nsnapf,status='old', &
           file= snap_file,err=81)
         snapfile=snap_file
         go to 3032
 81    continue
         open(unit=nsnapf,status='old', &
           file= input_dir(1:lindir)//snap_file,err=83  )
         snapfile=input_dir(1:lindir)//snap_file
         go to 3032
  83   continue
         open(unit=nsnapf,status='old', &
           file= snapextin       )
         snapfile=snapextin
      else
      open(unit=nsnapf,status='old', &
           file='efit_snap.dat',err=3030)
      snapfile='efit_snap.dat'
      go to 3032
 3030 continue
      open(unit=nsnapf,status='old', &
           file= input_dir(1:lindir)//'efit_snap.dat'         )
      snapfile=input_dir(1:lindir)//'efit_snap.dat'
      endif
 3032 continue
      if (appendsnap.eq.'K'.or.appendsnap.eq.'KG') then
      do i=1,1000000
         read (nsnapf,9991,end=9981) tmpdata
         if (INDEX(tmpdata,'&efitin')/=0) go to 381
      enddo
 381  continue
      do i=1,1000000
         write (neqdsk,9991) tmpdata
         read (nsnapf,9991,end=9981) tmpdata
         if (INDEX(tmpdata,'/')/=0) then
            write (neqdsk,9991) tmpdata
            go to 9979
         endif
      enddo
 9979 close (unit=nsnapf)
 9981 continue
 9991 format (a)
      endif
!
      close(unit=neqdsk)
      ltbdir=ltbdis
      table_dir = table_s(1:ltbdis)
 4042 format (1x,a42,1x,a3)
 4958 format ('#!/bin/csh -f')
 4960 format ('      runefit.sc k',a12)
 4962 format ('#',/,'exit')
 4970 format (2x,a,1x,a,1x,a,1x,a)
 4980 format (i5)
 5000 format (2e12.6)
 5500 format (/,10x,'EFITD ',a3,' x ',a3,' Version  ',2a5,/)
 6000 format (/,1x,'type mode (2=file, 3=snap, 4=time', &
               ', 5=input, 6=com file, 7=snap_ext):')
 6040 format (/,1x,'type shot #, start time(ms), time step(ms), steps' &
        ,'(<401):')
 6080 format (/,1x,'type limiter position (cm, 0=ptdata):')
 6090 format(' enter number of extra field lines to trace:')
 6091 format(' enter scrape off depth(m),'/ &
       '       sense of tracing (+1 for down, -1 for up),'/ &
       '       ixstrt (+1 for start on outside, -1' &
       ' for start inside):')
 6100 format(/,1x,48htype plot mode (0=none, 1=tektronix, 2=versatec, &
           ,17h 3=qms, -=x ray):)
 6200 format (/,1x,22hnumber of time slices?)
 6220 format (/,1x,22htype input file names:)
 6230 format (1x,1h#)
 6240 format (a)
 6600 format (/,1x,'good shot list file name ( 0=tty) ?')
 6610 format (/,1x,'command file name ( 0=none) ?')
 6617 format (/,1x,'type snap file extension (def for default):')
 6620 format (a)
 6700 format (a1,a12)
20000 format (/,1x,'shot data not on disk')
30000 format (i9)
30200 format (10f3.0)
      end



!
!   This routine is required if the CVS revision numbers are to 
!   survive an optimization.
!
!
!   1997/05/23 22:53:27 peng
!
      subroutine efitdx_rev(i)
      CHARACTER*100 opt
      character*10 s 
      if( i .eq. 0) s =  &
      '@(#)$RCSFILE: efitdx.for,v $ 4.61\000'
      return
      end

! jm.s
real*8 function linear(x,xa,ya,n)

      implicit none
      
      real*8, intent(in) :: x
      real*8, dimension(65), intent(inout) :: xa, ya
      integer, intent(in) :: n
      
      integer :: klo, khi, k
      real*8 :: h, b, a
      
      klo=1
      khi=n
      
      do 
         if (khi-klo <= 1) exit

         k=(khi+klo)/2
         if (xa(k) > x) then
            khi=k
         else
            klo=k
         endif
      enddo

      h=xa(khi)-xa(klo)
      
      if (h == 0.0) then
         print *, 'Bad xa input to routine linear'
! MPI >>>
         ! NOTE : Do NOT need to replace STOP command since function currently unused
         stop
! MPI <<<
      endif

      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
    
      linear = a*ya(klo)+b*ya(khi)

      end function linear
! jm.e

! MPI >>>
#if defined(USEMPI)
    ! Shutdown MPI before calling STOP to terminate program
    subroutine mpi_stop()
    
      include 'modules1.f90'
      include 'mpif.h'
      
      if (allocated(dist_data)) deallocate(dist_data)
      if (allocated(dist_data_displs)) deallocate(dist_data_displs)
      if (allocated(fwtgam_mpi)) deallocate(fwtgam_mpi)
      
      if (rank == 0) then
        print *, 'STOPPING MPI'
      endif
      call MPI_FINALIZE(ierr)
      STOP
    
    end subroutine mpi_stop
#endif
! MPI <<<
