#include "config.f"
!**********************************************************************
!>  
!!    efit is the main driver for equilibrium analysis.
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
      program efit
      use commonblocks
      use set_kinds
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
#if defined(USEMPI)
      include 'mpif.h'
#endif
      character inp1*4,inp2*4
      integer*4 :: k,krord,kzord,nargs,iargc,finfo,kerror,terr,ioerr, &
                   nwrk,mfila,ktime,mtear,ks

      integer*4 :: iend1,iend2
      character*80 :: cmdline
      parameter (krord=4,kzord=4)

      kerror = 0
      kwake = 0
!------------------------------------------------------------------------------
!--   MPI if we have it
!------------------------------------------------------------------------------ 
#if defined(USEMPI)
! Initialize MPI environment
      call MPI_INIT_THREAD(MPI_THREAD_SINGLE,terr,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
! Arrays can only be allocated after MPI has been initialized because dimension is # of processes
      allocate(dist_data(nproc),dist_data_displs(nproc),fwtgam_mpi(nstark,nproc))
#ifdef DEBUG_PLTS
      if (nproc.gt.1) then
        call errctrl_msg('efit', &
               'Surface files debugging/plotting is serial only')
        stop
      endif
#endif
#else
      rank  = 0
      nproc = 0
#endif

      ! Set global constants for each rank
      call set_constants()
!------------------------------------------------------------------------------
!--   Set external variables and print build info
!------------------------------------------------------------------------------ 
      call set_exvars
      ! If environment variable exists, then override values from set_exvars
      call getenv("link_efit",link_efitx)
      call getenv("link_store",link_storex)
      if (link_efitx(1:1).ne.' ') then
        table_dir=trim(link_efitx)//'green/'
        input_dir=trim(link_efitx)
      endif
      if (link_storex(1:1).ne.' ')  store_dir=trim(link_storex)
!----------------------------------------------------------------------
!--   Read in grid size from command line and set global variables   --
!--   ONLY root process reads command-line arguments                 --
!----------------------------------------------------------------------
      nw = 0
      nh = 0
      if (rank == 0) then
        nargs = iargc()
! WARNING: is this true?? if so LF95 -> USEMPI ...
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
        read (inp1,'(i4)',iostat=ioerr) nw
        if(ioerr.ne.0) read (inp2,'(i4)') nh

! Ensure grid size is defined
        if (nw == 0) then
          call errctrl_msg('efit', &
                 'Must specify grid dimensions as arguments')
          stop
        endif
        if(nh == 0) nh = nw
      endif

#if defined(USEMPI)
! Distribute command-line information (meaning grid dimensions) to all processes if necessary
      if (nproc > 1) then
        call MPI_BCAST(nw,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(nh,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      endif
#endif

      if (nw .le. 129) then
        npoint=800
      elseif (nw .le. 516) then
        npoint=3200
      else
        npoint=12800
      endif
      nwnh=nw*nh
      nh2=2*nh
      nwrk=2*(nw+1)*nh
!     nwwf=2*nw
      nwwf=3*nw
      nwf=nwwf
      kubicx=4
      kubicy=4
      lubicx=nw-kubicx+1
      lubicy=nh-kubicy+1
      kujunk=kubicx*kubicy*lubicx*lubicy
      boundary_count=2*nh+2*(nw-2)
      lr0=nw-krord+1
      lz0=nh-kzord+1
      nxtrap=npoint
      mfila=10
      
      call read_optin
      call inp_file_ch(nw,nh,ch1,ch2)

      call get_opt_input(ktime)
      call get_eparmdud_defaults()
      ntime = ktime

      select case (kdata)
      case (1)
        call read_omas_in1(ifname(1))     !this assumes machine is always the same
      case (2)
        call read_dirs_shot(ifname(1))     !this assumes machine is always the same
      case (4)
        call read_dirs_shot('efit_time.dat')
      case (5)
        if (snapextin.ne.'none') then ! could come from efit.input
          call read_dirs_shot('efit_snap.dat_'//adjustl(snapextin))     !this assumes machine is always the same
        else
          call read_dirs_shot('efit_snap.dat')
        endif
      case (7)
        call read_dirs_shot('efit_snap.dat_'//adjustl(snapextin))     !this assumes machine is always the same
      case default
        call read_dirs_shot('efit_snap.dat')
      end select

      call set_table_dir
      call read_eparmdud
      call get_eparmdud_dependents
!----------------------------------------------------------------------
!--   Global Allocations                                             --
!----------------------------------------------------------------------
      include 'global_allocs.f90'
      call set_mod_arrays
      call set_ecom_mod1_arrays
      call set_ecom_mod2_arrays
!----------------------------------------------------------------------
!--   get data                                                       --
!----------------------------------------------------------------------
      call efit_read_tables
  20  call getsets(ktime,mtear,kerror)
#if defined(USEMPI)
      if (nproc > 1) &
        call MPI_ALLREDUCE(kerror,MPI_IN_PLACE,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
      if (kerror.gt.0) then
        call errctrl_msg('efit','Aborting due to fatal error in getsets')
        call mpi_abort(MPI_COMM_WORLD,ierr) ! kill all processes, something is wrong with the setup.
      endif
#else
      if(kerror.gt.0) stop
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
!--   start simulation for KTIME time slices per rank                --
!----------------------------------------------------------------------
      do k=1,ktime
        ks=k ! ks=1,2,3... in serial, but ks=1,1,1,... in parallel
!----------------------------------------------------------------------
!--     set up data                                                  --
!----------------------------------------------------------------------
        
#ifdef DEBUG_LEVEL2
        write(6,*) ' Entering prtoutheader subroutine'
#endif
        call prtoutheader()
#ifdef DEBUG_LEVEL2
        write(6,*) ' Entering data_input subroutine'
#endif

        call data_input(ks,iconvr,ktime,mtear,kerror)

! Check that the grid sizes are compatible with Buneman's algorithm
! (this is not relevant to pefit - not yet available)
        if (ibunmn.ne.0) then
          if (nh .ne. 0) then
            select case (nh)
            case (3,5,9,17,33,65,129,257,513,1025,2049)
              ! all good
            case default
              call errctrl_msg('efit', &
                   'Chosen grid dimensions cannot be run')
              stop
            end select
          endif
          select case (nw)
          case (3,5,9,17,33,65,129,257,513,1025,2049)
            ! all good
          case default
            call errctrl_msg('efit', &
                 'Chosen grid dimensions cannot be run')
            stop
          end select
        endif

#ifdef DEBUG_LEVEL2
        write(6,*) 'Entering errctrl_setstate'
#endif
        call errctrl_setstate(rank,time(ks))
        if ((kerror.gt.0).or.(iconvr.lt.0)) then
          if (k.lt.ktime) then
            kerrot(ks)=kerror
            cycle
          else
            exit
          endif
        endif
        if (kautoknt .eq. 1) then
          call autoknot(ks,iconvr,ktime,mtear,kerror)
        else
!----------------------------------------------------------------------
!--       initialize current profile                                 --
!----------------------------------------------------------------------
#ifdef DEBUG_LEVEL2
          write(6,*) 'Entering inicur subroutine'
#endif
          call inicur(ks)
!----------------------------------------------------------------------
!--       get equilibrium                                            --
!----------------------------------------------------------------------
#ifdef DEBUG_LEVEL2
          write(6,*) 'Entering fit subroutine'
#endif
          call fit(ks,kerror)
          if (kerror.gt.0) then
            if (k.lt.ktime) then
              kerrot(ks)=kerror
              cycle
            else
              exit
            endif
          endif
        endif
!----------------------------------------------------------------------
!--     post processing for graphic and text outputs                 --
!----------------------------------------------------------------------
#ifdef DEBUG_LEVEL2
        write(6,*) 'Entering shapesurf'
#endif
        call shapesurf(ks,ktime,kerror)
        if (kerror.gt.0) then
          if (k.lt.ktime) then
            kerrot(ks)=kerror
            cycle
          else
            exit
          endif
        endif
!DEPRECATED        if (mtear.ne.0) call tearing(ks,mtear,kerror)
!DEPRECATED        if (kerror.gt.0) then
!DEPRECATED          if (k.lt.ktime) then
!DEPRECATED            kerrot(ks)=kerror
!DEPRECATED            cycle
!DEPRECATED          else
!DEPRECATED            exit
!DEPRECATED          endif
!DEPRECATED        endif
#ifdef DEBUG_LEVEL1
        write (6,*) 'Main/PRTOUT ks/kerror = ', ks, kerror
#endif
        call prtout(ks)
        if((kwaitmse.ne.0).and.(kmtark.gt.0)) call fixstark(-ks,kerror)
!----------------------------------------------------------------------
!--     write A and G EQDSKs                                         --
!----------------------------------------------------------------------
#ifdef DEBUG_LEVEL1
        write (6,*) 'Main/WQEDSK ks/kerror = ', ks, kerror
#endif
        call weqdsk(ks)
        if (iconvr.ge.0) then
           call shipit(ktime,ks,ks)
!DEPRECATED           call wtear(mtear,ks)
        endif
#ifdef USE_NETCDF
#ifdef DEBUG_LEVEL2
        write (6,*) 'Main/wmeasure ks/kerror = ', ks, kerror
#endif
        call wmeasure(ktime,ks,ks,1)
#endif
!----------------------------------------------------------------------
! --    write Kfile if needed                                        --
!----------------------------------------------------------------------
#if defined(USE_SNAP)
        if (kdata.eq.3 .or. kdata.eq.7) then
          if(write_Kfile) call write_K2(ks,kerror)
        endif
#endif
        if(k.lt.ktime) kerrot(ks)=kerror
      enddo

      if(kwake.ne.0) go to 20
#ifdef USE_NETCDF
      call wmeasure(ktime,1,ktime,2)
#else
      if (.not.((iand(iout,2).eq.0).and.(iand(iout,4).eq.0))) &
        write(nttyo,*) 'netcdf needs to be linked to write m-files'
#endif
      call wtime(ktime)

#if defined(USEMPI)
      ! Finalize MPI
      if(allocated(dist_data)) deallocate(dist_data)
      if(allocated(dist_data_displs)) deallocate(dist_data_displs)
      if(allocated(fwtgam_mpi)) deallocate(fwtgam_mpi)
      call errctrl_msg('efit','Done processing',3)
      call mpi_finalize(ierr)
#endif
      stop
      end program efit
