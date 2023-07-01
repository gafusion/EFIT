#include "config.f"
!**********************************************************************
!>
!!    get_opt_input performs command line inputing\n
!!
!!    kdata:
!!      1: mimics option 2 but reads input from an hdf5
!!           file that has the IMAS format
!!      2: produces g-files (and others) from k-files
!!      3-7: query databases for diagnostic inputs
!!        3,7: produces g-files (and others)
!!        5: produces k-files from a snap file
!!        4,6: variations of snap files??
!!      8: indended for CUDA parallel execution which has not been setup
!!      -#: behaves the same as #, but sets ilaser=1
!!
!**********************************************************************
      subroutine get_opt_input(ktime)
      use opt_input
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(out) :: ktime
      integer*4 i,j
      logical file_stat
#if defined(USEMPI)
      include 'mpif.h'
#endif
      character*86 snap_ext     
      snap_ext='none'
      snapextin='none'

      ! ONLY root process allowed to interface with terminal
      if (rank == 0) then
        if (use_opt_input .eqv. .false.) then
          write (nttyo,6000)
          read (ntty,*) kdata
        else
          kdata = mode_in
        endif
      endif
#if defined(USEMPI)
      if(nproc > 1) &
        call MPI_BCAST(kdata,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
!----------------------------------------------------------------------
!--   Turn on ilaser                                                 --
!----------------------------------------------------------------------
      if (kdata.lt.0) then
        kdata=-kdata
        ilaser=1
      endif

      if (kdata.eq.4) return
!----------------------------------------------------------------------
!--   Check that input option is valid                               --
!----------------------------------------------------------------------
#if defined(USE_HDF5)
#if defined(USE_SNAP)
      if (kdata.lt.1 .or. kdata.gt.7) then
        call errctrl_msg('get_opt_input', 'kdata run type is not available')
        stop
      endif
#else
      if (kdata.lt.1 .or. kdata.gt.2) then
        call errctrl_msg('get_opt_input', 'kdata run type is not available')
        stop
      endif
#endif
#else
#if defined(USE_SNAP)
      if (kdata.lt.2 .or. kdata.gt.7) then
        call errctrl_msg('get_opt_input', 'kdata run type is not available')
        stop
      endif
#else
      if (kdata.ne.2) then
        call errctrl_msg('get_opt_input', 'kdata run type is not available')
        stop
      endif
#endif
#endif
!      TODO: pefit has not been setup
!      if (kdata.lt.1 .or. kdata.gt.8) then
!        call errctrl_msg('get_opt_input', 'kdata run type is not available')
!        stop
!      endif
!----------------------------------------------------------------------
!--   Get command-line input for each mode                           --
!----------------------------------------------------------------------
      if(kdata.eq.1) ALLOCATE(ifname(1))
      rank0: if (rank == 0) then
        select case (kdata)
        case (2) 
          if (use_opt_input .eqv. .false.) then
            write (nttyo,6200)
            read (ntty,*) ktime
            ALLOCATE(ifname(ktime))
            write (nttyo,6220)
            do i=1,ktime
              write (nttyo,6230)
              read (ntty,6240) ifname(i)
            enddo
          else
            ktime = steps_in
            ALLOCATE(ifname(ktime))
            ifname(1:ktime) = inpfile(1:ktime)
          endif

        case (1)
          if (use_opt_input .eqv. .false.) then
            write (nttyo,6220)
            read (ntty,6240) ifname(1)
          else
            ifname(1) = inpfile(1)
          endif
#if defined(USE_HDF5)
          inquire(file=trim(ifname(1)),exist=file_stat)
          if (.not. file_stat) then
            call errctrl_msg('get_opt_input',trim(ifname(1))//' not found')
            stop
          endif
          call fch5init
          call open_oldh5file(trim(ifname(1)),fileid,rootgid,h5in,h5err)
          call test_group(rootgid,"equilibrium",file_stat,h5err)
          if (.not. file_stat) then
            call errctrl_msg('get_opt_input','equilibrium group not found')
            stop
          endif
          call open_group(rootgid,"equilibrium",eqid,h5err)
          call test_group(eqid,"code",file_stat,h5err)
          if (.not. file_stat) then
            call errctrl_msg('data_input','code group not found')
            stop
          endif
          call open_group(eqid,"code",cid,h5err)
          call test_group(cid,"parameters",file_stat,h5err)
          if (.not. file_stat) then
            call errctrl_msg('data_input','parameters group not found')
            stop
          endif
          call open_group(cid,"parameters",pid,h5err)
          call test_group(pid,"time_slice",file_stat,h5err)
          if (.not. file_stat) then
            call errctrl_msg('data_input','time_slice group not found')
            stop
          endif
          call get_nmembers(pid,"time_slice",ktime,h5err)
          call close_group("parameters",pid,h5err)
          call close_group("code",cid,h5err)
          call close_group("equilibrium",eqid,h5err)
          call close_h5file(fileid,rootgid,h5err)
#else
          ! this code should not be reachable
          call errctrl_msg('get_opt_input','HDF5 needs to be linked')
          stop
#endif

        case (3,7)
          ! Snap-Extension mode
          if (kdata.eq.7) then
            if (use_opt_input .eqv. .false.) then
              write (nttyo,6617)
              read (ntty,6620) snap_ext
            else
              snap_ext = snapext_in
            endif
          endif
          if (use_opt_input .eqv. .false.) then
            write (nttyo,6040)
            read (ntty,*) ishot,timeb,dtime,ktime
            shot_in = ishot
            starttime_in = timeb
            deltatime_in = dtime
            steps_in= ktime 
          else
            ishot = shot_in
            timeb = starttime_in
            dtime = deltatime_in
            ktime = steps_in
          endif

        case (5,6)
          if (use_opt_input .eqv. .false.) then
            write (nttyo,6610)
            read (ntty,6620) cmdfile_in
            write (nttyo,6600)
            read (ntty,6620) shotfile_in
            if (shotfile_in.eq.'0' .or. shotfile_in.eq.'') then
              write (nttyo,6040)
              read (ntty,*) ishot,timeb,dtime,ktime
            endif
          else
            if (shotfile_in.eq.'0' .or. shotfile_in.eq.'') then
              ishot = shot_in
              timeb = starttime_in
              dtime = deltatime_in
              ktime = steps_in
              snapextin = snapext_in
            endif
          endif

        end select
      endif rank0

      select case (kdata)
      case (3,7) 
#if defined(USEMPI)
        if (nproc > 1) then
          call MPI_BCAST(ishot,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(timeb,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(dtime,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(snap_ext,82,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
        endif
#endif 
        snapextin = snap_ext
      case (5,6)
        ! only single process is used in this mode
#if defined(USEMPI)
        if(nproc > 1) &
            write(nttyo,*) 'Warning: only 1 processor is active'
#endif 
        return
      end select
     
#if defined(USEMPI)
      if (nproc > 1) then
        if (rank == 0) then
          ! Ensure there are not more processors than time slices
          if (nproc > ktime) then
            call errctrl_msg('get_opt_input', &
                             'MPI processes have nothing to do')
          endif
          ! Warn if the time slice distribution is not balanced
          if (mod(ktime,nproc) .ne. 0) &
            write(nttyo,*) 'Warning: time slices are not balanced across processors'
        endif

        ! Broadcast inputs 
        call MPI_BCAST(ktime,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        ! if (nproc > ktime) stop ! optional?...

        if (kdata.eq. 1 .or. kdata.eq.2) then
          ! Distribute steps among ALL processes
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
              ! Input filenames are up to 80 characters and displacements given as number of bytes
              dist_data_displs(i) = sum(dist_data(1:i-1))*80
            enddo
          endif
          ! Explicitly synchronize processes
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          ! Distribute time step and filename information to ALL processes
          call MPI_SCATTER(dist_data,1,MPI_INTEGER,ktime,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          if (kdata.eq.2) then
            ! Recall each filename 80 characters
            if (rank == 0) then
              dist_data(:) = dist_data(:)*80
              call MPI_SCATTERV(ifname,dist_data,dist_data_displs,MPI_CHARACTER, &
                     MPI_IN_PLACE,dist_data(rank+1),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
            else
              ALLOCATE(ifname(ktime))
              call MPI_SCATTERV(ifname,dist_data,dist_data_displs,MPI_CHARACTER, &
                                ifname,ktime*80,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
            endif
          else
            call MPI_BCAST(ifname(1),80,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
          endif
        endif
      endif
#endif

#if defined(USE_HDF5)
#if defined(USE_SNAP)
 6000 format (/,1x,'type mode (1=imas, 2=file, 3=snap, 4=time', &
               ', 5=input, 6=com file, 7=snap_ext):')
#else
 6000 format (/,1x,'type mode (1=imas, 2=file):')
#endif
#else
#if defined(USE_SNAP)
 6000 format (/,1x,'type mode (2=file, 3=snap, 4=time', &
               ', 5=input, 6=com file, 7=snap_ext):')
#else
 6000 format (/,1x,'type mode (2=file):')
#endif
#endif
! TODO: pefit has not been setup
! 6000 format (/,1x,'type mode (1=imas, 2=file, 3=snap, 4=time', &
!               ', 5=input, 6=com file, 7=snap_ext, 8=pefit):')
 6040 format (/,1x,'type shot #, start time(ms), time step(ms), steps:')
 6200 format (/,1x,'number of time slices?')
 6220 format (/,1x,'type input file names:')
 6230 format (1x,'#')
 6240 format (a)
 6600 format (/,1x,'good shot list file name (0=tty)?')
 6610 format (/,1x,'command file name (0=none)?')
 6617 format (/,1x,'type snap file extension (no extension if left blank):')
 6620 format (a)
     return
end subroutine get_opt_input
