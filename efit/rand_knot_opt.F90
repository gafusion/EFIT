#include "config.f"
!*******************************************************************
!!
!!    rand_knot varies knot positions randomly until convergence
!!      is achieved
!!
!!    @param ks: time index
!!    @param ktime : number of time slices
!!    @param kerror: error flag
!!
!*******************************************************************
      subroutine rand_knot(ks,ktime,kerror)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
#if defined(USEMPI)
      include 'mpif.h'
#endif

      integer*4, intent(in) :: ks,ktime
      integer*4, intent(out) :: kerror
      integer*4 i,j,k,kloop,loc(1),nvary,saveiter
      real*8 berror
      real*8, dimension(nproc) :: errall
      real*8, dimension(:), allocatable :: lbnd,krange,kpos,kbest

      kerror = 0

      ! Run with input settings first
      if(iconvr.lt.0) return
      call set_init(ks)
      call fit(ktime,ks,kerror)
      if(rank==0) write(6,*) 'Initial error is: ',terror(ks)
      if ((kerror == 0) .and. (terror(ks).le.error)) then
        if(rank==0) write(6,*) 'Input settings converged'
        return
      endif
      berror=terror(ks)

      ! Setup bounds for varying each knot
      nvary=0
      if(kppfnc .eq. 6) nvary=nvary+kppknt-2
      if(kfffnc .eq. 6) nvary=nvary+kffknt-2
      if(kwwfnc .eq. 6) nvary=nvary+kwwknt-2
      if(keefnc .eq. 6) nvary=nvary+keeknt-2
      if (nvary.le.0) then
        if(rank==0) write(6,*) 'No knots to vary'
        return
      endif
      allocate(lbnd(nvary),krange(nvary),kpos(nvary),kbest(nvary))
      j=0
      if (kppfnc .eq. 6) then
        do i=2,kppknt-1
          j=j+1
          lbnd(j) = ppknt(i)-appdf(i)*(ppknt(i)-ppknt(i-1))
          krange(j) = appdf(i)*(ppknt(i+1)-ppknt(i-1))
          kbest(j) = ppknt(i)
        enddo
      endif
      if (kfffnc .eq. 6) then
        do i=2,kffknt-1
          j=j+1
          lbnd(j) = ffknt(i)-affdf(i)*(ffknt(i)-ffknt(i-1))
          krange(j) = affdf(i)*(ffknt(i+1)-ffknt(i-1))
          kbest(j) = ffknt(i)
        enddo
      endif
      if (kwwfnc .eq. 6) then
        do i=2,kwwknt-1
          j=j+1
          lbnd(j) = wwknt(i)-awwdf(i)*(wwknt(i)-wwknt(i-1))
          krange(j) = awwdf(i)*(wwknt(i+1)-wwknt(i-1))
          kbest(j) = wwknt(i)
        enddo
      endif
      if (keefnc .eq. 6) then
        do i=2,keeknt-1
          j=j+1
          lbnd(j) = eeknt(i)-aeedf(i)*(eeknt(i)-eeknt(i-1))
          krange(j) = aeedf(i)*(eeknt(i+1)-eeknt(i-1))
          kbest(j) = eeknt(i)
        enddo
      endif

#if defined(USEMPI)
      if (nproc > 1) then
        if (rank == 0) then
          ! TODO: it should be possible to distribute processors effectively for multiple times and 
          !       leave others for knot optimization, but that logic hasn't been setup yet
          !       (only can parallelize over one or the other and times takes precedence)
          if (nproc > kakloop) then
            ! Ensure there are not more processors than time slices
            call errctrl_msg('rand_knot', &
                             'MPI processes have nothing to do')
          endif
          ! Warn if the time slice distribution is not balanced
          if (mod(kakloop,nproc) .ne. 0) &
            write(nttyo,*) 'Warning: knot varitaion loops are not balanced across processors'
        endif
        ! This could be avoided if MPI is handled differently
        if(nproc > kakloop) stop
        kloop=kakloop/nproc
        do k=1,mod(kakloop,nproc)
          if(rank == k) kloop=kloop+1
        enddo
      else
        kloop=kakloop
      endif
#else
      kloop=kakloop
#endif 

      ! Vary knots randomly until convergence or max iterations is reached
      do k=1,kloop
        kerror = 0

        ! With a bit of careful debugging we should be able to avoid this...
        call data_input(ks,ktime,kerror)
        if (kerror.gt.0 .or. iconvr.lt.0) then
          write(6,*) 'Error re-reading input, aborting'
          return
        endif

        ! Set random knot positions
        call random_number(kpos)
        kpos=lbnd+kpos*krange
        j=0
        if (kppfnc .eq. 6) then
          do i=2,kppknt-1
            j=j+1
            ppknt(i) = kpos(j)
          enddo
        endif
        if (kfffnc .eq. 6) then
          do i=2,kffknt-1
            j=j+1
            ffknt(i) = kpos(j)
          enddo
        endif
        if (kwwfnc .eq. 6) then
          do i=2,kwwknt-1
            j=j+1
            wwknt(i) = kpos(j)
          enddo
        endif
        if (keefnc .eq. 6) then
          do i=2,keeknt-1
            j=j+1
            eeknt(i) = kpos(j)
          enddo
        endif

        ! Run fit
        mxiter = kakiter
        call set_init(ks)
        call fit(ktime,ks,kerror)
        if(kerror .ne. 0) terror(ks)=999.
        write(6,*) 'Iteration ',(k-1)*nproc+rank+1,' has error: ',terror(ks)

#if defined(USEMPI)
        ! Get best error and corresponding knot positions from all processors
        if (nproc > 1) then
          errall=999.
          call MPI_GATHER(terror(ks),1,MPI_DOUBLE_PRECISION,errall,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
          if(rank==0) loc=minloc(errall)
          call MPI_BCAST(loc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(terror(ks),1,MPI_DOUBLE_PRECISION,loc-1,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(kpos,nvary,MPI_DOUBLE_PRECISION,loc-1,MPI_COMM_WORLD,ierr)
        endif
#endif

        ! Check convergence and save best error and knots
        if (terror(ks).le.error) then
          if (rank == 0) then
            write(6,*) 'New knot locations converged'
            return
          else
            ! Finished with parallel processing
#if defined(USEMPI)
            call mpi_finalize(ierr)
#endif
            stop
          endif
        endif
        if (rank.eq.0 .and. terror(ks).le.berror) then
          berror=terror(ks)
          kbest=kpos
        endif
      enddo

#if defined(USEMPI)
      ! Finished with parallel processing
      if (rank.gt.0) then
        call mpi_finalize(ierr)
        stop
      endif
#endif

      ! Re-run the best case to use as the final output
      kerror = 0
      ! With a bit of careful debugging we should be able to avoid this...
      call data_input(ks,ktime,kerror)
      if (kerror.gt.0 .or. iconvr.lt.0) then
        write(6,*) 'Error re-reading input, aborting'
        return
      endif
      j=0
      if (kppfnc .eq. 6) then
        do i=2,kppknt-1
          j=j+1
          ppknt(i) = kbest(j)
        enddo
      endif
      if (kfffnc .eq. 6) then
        do i=2,kffknt-1
          j=j+1
          ffknt(i) = kbest(j)
        enddo
      endif
      if (kwwfnc .eq. 6) then
        do i=2,kwwknt-1
          j=j+1
          wwknt(i) = kbest(j)
        enddo
      endif
      if (keefnc .eq. 6) then
        do i=2,keeknt-1
          j=j+1
          eeknt(i) = kbest(j)
        enddo
      endif
      call set_init(ks)
      call fit(ktime,ks,kerror)
      write(6,*) 'Lowest error found is ',terror(ks)

      return
      end subroutine rand_knot
