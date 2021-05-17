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
     program efitd
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


