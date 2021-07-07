subroutine get_opt_input(ktime)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          getsets performs inputing, pulled out of getsets for    **
!**          generalization of EFIT                                  **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**********************************************************************
      use set_kinds
      include 'eparmdud129.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      USE mpi_efit
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      character*82 snap_ext
     

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
      if (nproc > 1) then
        call MPI_BCAST(kdata,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      endif

      if (kdata.eq.7) then
         if (rank == 0) then
           if (use_opt_input .eqv. .false.) then
             write (nttyo,6617)
             read (ntty,6620) snap_ext
           else
             snap_ext = snapext_in
           endif
         endif
         snapextin=snap_ext
         if (nproc > 1) then
           call MPI_BCAST(snap_ext,82,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
         endif
      endif


      if (kdata.eq.2) then 
       if (rank == 0) then
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
          do i=1,ktime
             ifname(i) = inpfile_in(i)
          enddo
        endif
       endif

      elseif (kdata.eq.3 .or. kdata.eq.7) then
    
       if (kwake.eq.0) then
        if (rank == 0) then
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
        endif
       endif
    else
     return
    endif
     
#if defined(USEMPI)
        if (nproc > 1) then
          call MPI_BCAST(ishot,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(ktime,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(timeb,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(dtime,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        endif
#endif 

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

      3000 continue


 5500 format (/,10x,'EFITD Version  ',2a5,/)
 6000 format (/,1x,'type mode (2=file, 3=snap, 4=time', &
               ', 5=input, 6=com file, 7=snap_ext,',    &
               ' 8=pefit):')
 6040 format (/,1x,'type shot #, start time(ms), time step(ms), steps' &
        ,'(<1001):')
 6200 format (/,1x,'number of time slices?')
 6220 format (/,1x,'type input file names:')
 6230 format (1x,'#')
 6240 format (a)
 6617 format (/,1x,'type snap file extension (def for default):')
 6620 format (a)
     return
end subroutine
