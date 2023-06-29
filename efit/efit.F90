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
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'

      implicit none
#if defined(USEMPI)
      include 'mpif.h'
#endif
      integer*4 k,krord,kzord,finfo,kerror,terr,ioerr, &
                   nwrk,ktime,ks
      integer*4 i,iend1,iend2
      character inp1*4,inp2*4
      character*80 :: cmdline
      parameter (krord=4,kzord=4)

      kerror = 0
!------------------------------------------------------------------------------
!--   MPI if we have it
!------------------------------------------------------------------------------ 
#if defined(USEMPI)
! Initialize MPI environment
      call MPI_INIT_THREAD(MPI_THREAD_SINGLE,terr,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
! Arrays can only be allocated after MPI has been initialized because dimension is # of processes
      allocate(dist_data(nproc),dist_data_displs(nproc))
#ifdef DEBUG_PLTS
      if (nproc.gt.1) then
        call errctrl_msg('efit', &
               'Surface files debugging/plotting is serial only')
        stop
      endif
#endif
#else
      rank  = 0
      nproc = 1
#endif

      ! Set global constants for each rank
      pi = 4.0_dp*atan(1.0_dp) ! calculate pi to machine precision
      twopi = 2.0*pi
      radeg = pi/180.0
      tmu0=twopi*tmu
      tmu02=tmu0*2.0
!------------------------------------------------------------------------------
!--   Set external variables and print build info
!------------------------------------------------------------------------------ 
      call set_extvars
!----------------------------------------------------------------------
!--   Read in grid size from command line and set global variables   --
!--   ONLY root process reads command-line arguments                 --
!----------------------------------------------------------------------
      nw = 0
      nh = 0
      if (rank == 0) then
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
        if(ioerr.eq.0) read (inp2,'(i4)') nh

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
      
      call read_optin
      call table_name_ch(nw,nh,ch1,ch2)

      call get_opt_input(ktime)
      call set_eparm_defaults()
      ntime = ktime

      select case (kdata)
      case (1)
        call read_dirs_shot_imas(ifname(1))
      case (2)
        call read_dirs_shot(ifname(1))
      case (4)
        call read_dirs_shot('efit_time.dat')
      case default
        if (snapextin.eq.'none') then
          snap_file = 'efit_snap.dat'
        else
          snap_file = 'efit_snap.dat_'//adjustl(snapextin)
        endif
      end select

      call set_table_dir
      call read_machinein
      call set_eparm_dependents
!----------------------------------------------------------------------
!--   Global Allocations                                             --
!----------------------------------------------------------------------
      include 'global_allocs.f90' ! this prevents changing the machine during a run
!----------------------------------------------------------------------
!--   Read in Green Function Tables
!----------------------------------------------------------------------
      call read_tables
!----------------------------------------------------------------------
!--   K-file from snap mode                                          --
!----------------------------------------------------------------------
#if defined(USE_SNAP)
      if (kdata.eq.5.or.kdata.eq.6) then
#ifdef DEBUG_LEVEL2
        write(6,*) ' Entering write_k subroutine'
#endif
        call write_k(ktime,kerror)
        if (kerror.gt.0) then
          call errctrl_msg('efit','Error in write k-files',3)
        else
          call errctrl_msg('efit','Done writing k-files',3)
        endif
#if defined(USEMPI)
        call mpi_finalize(ierr)
#endif
        stop
      endif
#endif
!----------------------------------------------------------------------
!--   get measurement data
!----------------------------------------------------------------------
#ifdef DEBUG_LEVEL2
      write(6,*) ' Entering setup_data_fetch subroutine'
#endif
  20  call setup_data_fetch(ktime,kerror)
#if defined(USEMPI)
      if (nproc > 1) &
        call MPI_ALLREDUCE(kerror,MPI_IN_PLACE,1,MPI_INTEGER,MPI_MAX, &
                           MPI_COMM_WORLD,ierr)
      if (kerror.gt.0) then
        call errctrl_msg('efit', &
          'Aborting due to fatal error in setup_data_fetch')
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
        write(6,*) ' Entering print_header subroutine'
#endif
        call print_header()

#ifdef DEBUG_LEVEL2
        write(6,*) ' Entering data_input subroutine'
#endif
        call data_input(ks,iconvr,ktime,kerror)
        if ((kerror.gt.0).or.(iconvr.lt.0)) then
          if(k.lt.ktime) kerrot(ks)=kerror
          cycle
        endif
        if (req_valid) then
          ! don't solve times without plasma
          if (ivacum.eq.1) then
            call errctrl_msg('efit', &
            'Solution does contain any plasma, no outputs are generated')
            cycle
          endif
          ! don't solve times without mse
          if (sum(abs(swtgam)).gt.nstark*1.e-6_dp .and. kstark.eq.0) then
            call errctrl_msg('efit', &
              'No MSE data found, not solving for equilibrium')
            cycle
          endif
          ! don't write times without cer
          if (mse_usecer.ne.0 .and. &
              maxval(abs(tangam(ks,1:nmselp) &
                        -tangam_uncor(ks,1:nmselp))).le.1.e-10_dp) then
            call errctrl_msg('efit', &
              'No CER correction used, not solving for equilibrium')
            cycle
          endif
        endif

        if (kautoknt .eq. 1) then
#ifdef DEBUG_LEVEL2
          write(6,*) ' Entering autoknot subroutine'
#endif
          call autoknot(ks,iconvr,ktime,kerror)
        else
!----------------------------------------------------------------------
!--       initialize current profile                                 --
!----------------------------------------------------------------------
#ifdef DEBUG_LEVEL2
          write(6,*) 'Entering set_init subroutine'
#endif
          call set_init(ks)
          ! don't solve times without plasma solution
          if (req_valid) then
            if (icinit.eq.1) then
              if (abs(sum(pcurrt)).lt.1.e-3_dp) then
                call errctrl_msg('efit', &
                  'Insufficient current to create plasma, not solving')
                cycle
              endif
            else
              do i=1,nwnh
                if((xpsi(i).ge.0.0).and.(xpsi(i).le.1.0)) exit
              enddo
              if (i.ge.nwnh) then
                call errctrl_msg('efit', &
                  'Insufficient current to create plasma, not solving')
                cycle
              endif
            endif
          endif
!----------------------------------------------------------------------
!--       get equilibrium                                            --
!----------------------------------------------------------------------
#ifdef DEBUG_LEVEL2
          write(6,*) 'Entering fit subroutine'
#endif
          call fit(ks,kerror)
          if (kerror.gt.0) then
            if(k.lt.ktime) kerrot(ks)=kerror
            cycle
          endif
        endif
        ! for use in optimization loops
        if (ierchk.lt.0 .and. lflag.gt.0) then
          call write_a(ktime,ks)
          cycle
        endif
!----------------------------------------------------------------------
!--     post processing for graphic and text outputs                 --
!----------------------------------------------------------------------
#ifdef DEBUG_LEVEL2
        write(6,*) 'Entering shapesurf'
#endif
        call shapesurf(ks,ktime,kerror)
        if (kerror.gt.0) then
          if(k.lt.ktime) kerrot(ks)=kerror
          cycle
        endif
#ifdef DEBUG_LEVEL2
        write (6,*) 'Main/print_stats ks/kerror = ', ks, kerror
#endif
        call print_stats(ks)
        if((kwaitmse.ne.0).and.(kmtark.gt.0)) call fixstark(-ks,kerror)
!----------------------------------------------------------------------
!--     write A, G, and M EQDSKs
!----------------------------------------------------------------------
        if (iconvr.ge.0) then
#ifdef DEBUG_LEVEL2
          write (6,*) 'Main/write_a ks/kerror = ', ks, kerror
#endif
          call write_a(ktime,ks)
          if(abs(ierchk).gt.1 .and. lflag.gt.0) cycle
        endif
#ifdef DEBUG_LEVEL2
        write (6,*) 'Main/write_g ks/kerror = ', ks, kerror
#endif
        call write_g(ks)
#ifdef USE_NETCDF
#ifdef DEBUG_LEVEL2
        write (6,*) 'Main/write_m ks/kerror = ', ks, kerror
#endif
        call write_m(ktime,ks,ks,1)
#endif
!----------------------------------------------------------------------
! --    write k-file if needed                                        --
!----------------------------------------------------------------------
#if defined(USE_SNAP)
        if (kdata.eq.3 .or. kdata.eq.7) then
#ifdef DEBUG_LEVEL2
          write (6,*) 'Main/write_k2 ks/kerror = ', ks, kerror
#endif
          if(write_kfile) call write_k2(ks,kerror)
        endif
#endif
        if(k.lt.ktime) kerrot(ks)=kerror
      enddo

#ifdef USE_NETCDF
      call write_m(ktime,1,ktime,2)
#else
      if (.not.((iand(iout,2).eq.0).and.(iand(iout,4).eq.0))) &
        write(nttyo,*) 'netcdf needs to be linked to write m-files'
#endif
      call write_ot(ktime)

      call errctrl_msg('efit','Done processing',3)
#if defined(USEMPI)
      ! Finalize MPI
      if(allocated(dist_data)) deallocate(dist_data)
      if(allocated(dist_data_displs)) deallocate(dist_data_displs)
      call mpi_finalize(ierr)
#endif
      stop
      end program efit
