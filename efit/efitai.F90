!**********************************************************************
!>  
!!    efitai is the main driver for equilibrium analysis
!!    in the EFIT-AI project.  It is meant to be simpler than
!!    the full EFIT for more rapid development
!!
!!
!!     REFERENCES:\n                                                
!!          (1) L.L. Lao, H. St. John, R.D. Stambaugh, \
!!              A.G. Kellman, and W. Pfeiffer, Nuclear Fusion
!!             25 (1985) 1611.\n
!!          (2) L.L. Lao, H. St. John, R.D. Stambaugh, and 
!!               W. Pfeiffer, Nuclear Fusion 25 (1985) 1421.\n
!!          (3) L.L. Lao, J.R. Ferron, R.J. Groebner, W. Howl.
!!              H. St. John, E.J. Strait, and T.S. Taylor \n
!!              Nuclear Fusion 30 (1990) 1035.\n
!!          (4) L.L. Lao and T.H. Jensen Nuclear Fusion
!!              31 (1991) 1909.\n
!!          (5) L.L. Lao, H. St. John, et al, Fusion Sci. Technol. 
!!              48 (2005) 968.                                     
!!                                                                 
!**********************************************************************
      program efitai
      use commonblocks
      use set_kinds
      use mpi_efit
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      parameter (krord=4,kzord=4)
      character inp1*4,inp2*4
      integer*4 :: nargs, iargc, finfo, kerror, terr

      integer*4 :: iend1, iend2

      kerror = 0
!------------------------------------------------------------------------------
!     Initialize MPI environment
!------------------------------------------------------------------------------
      call MPI_INIT_THREAD(MPI_THREAD_SINGLE,terr,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

      ! Set global constants for each rank
      call set_constants()
!------------------------------------------------------------------------------
!--   Set external variables and print build info
!------------------------------------------------------------------------------ 
      call set_extvars
!----------------------------------------------------------------------
!-- Read in grid size from command line and set global variables     --
!-- ONLY root process reads command-line arguments                   --
!----------------------------------------------------------------------
      nw = 0
      nh = 0
      !TODO Need to process args here to get grid size
      if (rank == 0) then
       nargs = iargc()
! Using mpirun command so will have different number of arguments than serial case
       if (nh == 0) nh = nw

      endif

      ! Distribute command-line information (meaning grid dimensions) to all
      ! processes if necessary
      if (nproc > 1) then
       call MPI_BCAST(nw,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(nh,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      endif
      if (nw == 0 .or. nh == 0) then
        if (rank == 0) then
          call errctrl_msg('efitai', &
                           'Must specify grid dimensions as arguments')
        endif
        deallocate(dist_data,dist_data_displs)
        call mpi_finalize(ierr)
        STOP
      endif

      !TODO;  This stuff is mysterious and should be put somewhere
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
      
      call read_optin
      ! Create character versions of nw, nh for labelling
      call inp_file_ch(nw,nh,ch1,ch2)

      ! This chooses the EFIT mode, but efitai driver just does file input.
      !call get_opt_input(ktime)
      kdata=2
      ntime = ktime

      ! Black voodoo magic of setting poorly defined variables
      call set_eparm_defaults()
      call read_machinein(ifname(1))!this assume machine is always the same
      call set_eparm_dependents()

!----------------------------------------------------------------------
!-- Global Allocations                                               --
!   Now that we have the black voodoo magic of poorly documented integers
!   done, let's use that to set poorly documented variables
!----------------------------------------------------------------------
      include 'global_allocs.f90'
      call set_ecom_mod1_arrays
      call set_ecom_mod2_arrays

!----------------------------------------------------------------------
!-- get data                                                         --
!----------------------------------------------------------------------
      allocate(dist_data(nproc),dist_data_displs(nproc))
      !call set_table_dir
      !call read_tables
      !TODO: SEK: ZZ: Stopping here for now
      print *, 'Entering setup_data_fetch'
  20  call setup_data_fetch(ktime,kerror)
      print * ,'exiting setup_data_fetch'

#if defined(USEMPI)
      if (nproc > 1) then
        call MPI_ALLREDUCE(kerror,MPI_IN_PLACE,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
      endif
      if (kerror.gt.0) then
        call errctrl_msg('efitai', &
          'Aborting due to fatal error in setup_data_fetch')
        call mpi_abort(MPI_COMM_WORLD,ierr) ! kill all processes, something is wrong with the setup.
      endif
#else
      if (kerror.gt.0) then
        stop
      end if
#endif

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
        
        if(idebug>=2) write(6,*) ' Entering print_header subroutine'
        call print_header()
        if(idebug>=2) write(6,*) ' Entering data_input subroutine'

        call data_input(ks,iconvr,ktime,kerror)

        if(idebug>=2) write(6,*) ' Entering errctrl_setstate'
        call errctrl_setstate(rank,time(ks))
        if (kerror.gt.0) go to 500
        if (iconvr.lt.0) go to 500
        if (kautoknt .eq. 1) then
           call autoknot(ks,iconvr,ktime,kerror)
        else
!----------------------------------------------------------------------
!--  initialize current profile                                      --
!----------------------------------------------------------------------

           if(idebug>=2) write(6,*) 'Entering set_init subroutine'
           call set_init(ks)
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
        if (idebug /= 0) write (6,*) 'Main/PRINT_STATS ks/kerror = ', ks, kerror
        call print_stats(ks)
        if ((kwaitmse.ne.0).and.(kmtark.gt.0)) call fixstark(-ks,kerror)
!----------------------------------------------------------------------
!--  write A and G EQDSKs                                            --
!----------------------------------------------------------------------
        if (idebug /= 0) write (6,*) 'Main/write_g ks/kerror = ', ks, kerror
        call write_g(ks)
        if (iconvr.ge.0) then
           call write_a(ktime,ks)
!DEPRECATED           call wtear(mtear,ks)
        endif
        if (idebug /= 0) write (6,*) 'Main/write_m ks/kerror = ', ks, kerror
        call write_m(ktime,ks,ks,1)
!----------------------------------------------------------------------
! -- write Kfile if needed                                           --
!----------------------------------------------------------------------
  500 if (k.lt.ktime) then
        kerrot(ks)=kerror
        go to 100
      endif
      call write_m(ktime,1,ktime,2)
      call write_ot(ktime)

#if defined(USEMPI)
      ! Finalize MPI
      if (allocated(dist_data)) deallocate(dist_data)
      if (allocated(dist_data_displs)) deallocate(dist_data_displs)
      call errctrl_msg('efitai','Done processing',3)
      call mpi_finalize(ierr)
#endif
      stop
      end program
