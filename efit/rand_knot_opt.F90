#include "config.f"
!*******************************************************************
!>  
!!    rand_knot varies knot positions randomly until convergence
!!      is achieved
!!
!!    @param ks : Number of time slices
!!    @param ktime : Time index
!!    @param kerror: Error Flag
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
      integer*4 i,j,k,kloop,loc(1),nvary,npll,ndims,ncoef!,lead_rank
      real*8 berror
      !real*8, dimension(nproc) :: errall
      character filename*300
      integer*4, dimension(:), allocatable :: knot_sizes,coef_sizes
      real*8, dimension(:), allocatable :: lbnd,krange,kpos,kbest,errall
      real*8, dimension(:,:), allocatable :: kall,coefall

      kerror=0
      loc=1

      ! Run with input settings first
      if(iconvr.lt.0) return
      call set_init(ks)
      call fit(ks,kerror)

      ! Setup bounds for varying each knot
      nvary=0
      ndims=0
      if (kppfnc.eq.6) then
        nvary=nvary+kppknt-2
        ndims=ndims+1
      endif
      if (kfffnc.eq.6) then
        nvary=nvary+kffknt-2
        ndims=ndims+1
      endif
      if (kwwfnc.eq.6) then
        nvary=nvary+kwwknt-2
        ndims=ndims+1
      endif
      if (keefnc.eq.6) then
        nvary=nvary+keeknt-2
        ndims=ndims+1
      endif
      if (nvary.le.0) then
        if(rank.eq.0 .or. ttime.gt.1) write(6,*) 'No knots to vary'
        return
      endif
      allocate(lbnd(nvary),krange(nvary),kpos(nvary),kbest(nvary))
      allocate(knot_sizes(ndims),coef_sizes(ndims))
      j=0
      if (kppfnc.eq.6) then
        do i=2,kppknt-1
          j=j+1
          lbnd(j)=ppknt(i)-appdf(i)*(ppknt(i)-ppknt(i-1))
          krange(j)=appdf(i)*(ppknt(i+1)-ppknt(i-1))
          kbest(j)=ppknt(i)
        enddo
      endif
      if (kfffnc.eq.6) then
        do i=2,kffknt-1
          j=j+1
          lbnd(j)=ffknt(i)-affdf(i)*(ffknt(i)-ffknt(i-1))
          krange(j)=affdf(i)*(ffknt(i+1)-ffknt(i-1))
          kbest(j)=ffknt(i)
        enddo
      endif
      if (kwwfnc.eq.6) then
        do i=2,kwwknt-1
          j=j+1
          lbnd(j)=wwknt(i)-awwdf(i)*(wwknt(i)-wwknt(i-1))
          krange(j)=awwdf(i)*(wwknt(i+1)-wwknt(i-1))
          kbest(j)=wwknt(i)
        enddo
      endif
      if (keefnc.eq.6) then
        do i=2,keeknt-1
          j=j+1
          lbnd(j)=eeknt(i)-aeedf(i)*(eeknt(i)-eeknt(i-1))
          krange(j)=aeedf(i)*(eeknt(i+1)-eeknt(i-1))
          kbest(j)=eeknt(i)
        enddo
      endif

      ! Define additional output sizes
      i=0
      if (kppfnc.eq.6) then
        i=i+1
        knot_sizes(i)=kppknt-2
        coef_sizes(i)=kppcur
      endif
      if (kfffnc.eq.6) then
        i=i+1
        knot_sizes(i)=kffknt-2
        coef_sizes(i)=kffcur
      endif
      if (kwwfnc.eq.6) then
        i=i+1
        knot_sizes(i)=kwwknt-2
        coef_sizes(i)=kwwcur
      endif
      if (keefnc.eq.6) then
        i=i+1
        knot_sizes(i)=keeknt-2
        coef_sizes(i)=keecur
      endif
      ncoef=0
      do i=1,ndims
        ncoef=ncoef+coef_sizes(i)
      enddo

      ! Open file to store points tested :(ascii for now)
      if (rank.eq.0 .or. ttime.gt.1) then 
        call setfnm('n',ishot,int(time(ktime)),itimeu,'',filename)
        call open_new(neqdsk,filename,'formatted','quote')
        write(neqdsk,*) 1,(knot_sizes(i),i=1,ndims), &
                        (coef_sizes(i),i=1,ndims)
        write(neqdsk,*) terror(ks),kbest,(brsp(i),i=nfsum+1,nfsum+ncoef)
        write(6,*) 'Initial error is: ',terror(ks)
      endif

      ! Stop if converged
      if ((kerror.eq.0) .and. (terror(ks).le.error)) then
        if (rank.eq.0 .or. ttime.gt.1) then
          write(6,*) 'Input settings converged'
          close(unit=neqdsk)
          return
        else
          ! No need for parallel processing
#if defined(USEMPI)
          call mpi_finalize(ierr)
#endif
          stop
        endif
      endif
      berror=terror(ks)

      if (nproc.gt.1 .and. ttime.eq.1) then
#if defined(USEMPI)
        if (rank.eq.0) then
          ! TODO: it should be possible to distribute processors effectively for multiple times and 
          !       leave others for knot optimization, but that logic hasn't been setup yet
          !       (only can parallelize over one or the other and times takes precedence)
          if (nproc .gt. kakloop) then
            ! Ensure there are not more processors than knot attempts
            call errctrl_msg('rand_knot', &
                             'MPI processes have nothing to do')
          endif
          ! Warn that more knot positions will be attempted to equalize work across all processors
          if(mod(kakloop,nproc) .ne. 0) &
            write(nttyo,*) &
              'Warning: more knots attempted to balanced processors'
        endif
        ! This could be avoided if MPI is handled differently
        if(nproc.gt.kakloop) stop
        kloop=kakloop/nproc
        if(kakloop.gt.kloop*nproc) kloop=kloop+1
        npll=nproc
#endif 
      else
        kloop=kakloop
        npll=1
      endif
      allocate(errall(npll),kall(nvary,npll),coefall(ncoef,npll))

      ! Vary knots randomly until convergence or max iterations is reached
      do k=1,kloop
        kerror=0

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
        if (kppfnc.eq.6) then
          do i=2,kppknt-1
            j=j+1
            ppknt(i)=kpos(j)
          enddo
        endif
        if (kfffnc.eq.6) then
          do i=2,kffknt-1
            j=j+1
            ffknt(i)=kpos(j)
          enddo
        endif
        if (kwwfnc.eq.6) then
          do i=2,kwwknt-1
            j=j+1
            wwknt(i)=kpos(j)
          enddo
        endif
        if (keefnc.eq.6) then
          do i=2,keeknt-1
            j=j+1
            eeknt(i)=kpos(j)
          enddo
        endif

        ! Run fit
        mxiter=kakiter
        call set_init(ks)
        call fit(ks,kerror)
        if(kerror.ne.0) terror(ks)=999.
        if (ttime.eq.1) then
          write(6,*) 'Iteration ',(k-1)*nproc+rank+1,' has error: ',terror(ks)
        else
          write(6,*) 'Iteration ',k,' has error: ',terror(ks)
        endif

        ! Gather the knot positions tested and the errors from all processors
        errall=999.
        kall=999.
        coefall=999.
        if (nproc.gt.1 .and. ttime.eq.1) then
#if defined(USEMPI)
          call MPI_GATHER(terror(ks),1,MPI_DOUBLE_PRECISION,errall,1, &
                          MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
          call MPI_GATHER(kpos,nvary,MPI_DOUBLE_PRECISION,kall,nvary, &
                          MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
          call MPI_GATHER(brsp(nfsum+1:nfsum+ncoef),ncoef, &
                          MPI_DOUBLE_PRECISION,coefall,ncoef, &
                          MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
          loc=minloc(errall)
          terror(ks)=minval(errall)
          kpos=kall(:,loc(1))

          ! Broadcast the best error so that processors know to stop execution
          call MPI_BCAST(terror(ks),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!          ! Get best error and corresponding knot positions from all processors
!          errall=999.
!          call MPI_GATHER(terror(ks),1,MPI_DOUBLE_PRECISION,errall,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!          if(rank.eq.0) loc=minloc(errall)
!          call MPI_BCAST(loc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!          call MPI_BCAST(terror(ks),1,MPI_DOUBLE_PRECISION,loc-1,MPI_COMM_WORLD,ierr)
!          call MPI_BCAST(kpos,nvary,MPI_DOUBLE_PRECISION,loc-1,MPI_COMM_WORLD,ierr)
#endif
        else
          errall(1)=terror(ks)
          kall(:,1)=kpos
          coefall(:,1)=brsp(nfsum+1:nfsum+ncoef)
        endif

        ! Write knot positions attempted and error to file
        if (rank.eq.0 .or. ttime.gt.1) then
          do j=1,npll
            write(neqdsk,*) errall(j),kall(:,j),(coefall(i,j),i=1,ncoef)
          enddo
        endif

        ! Check convergence and save best error and knots
        if (terror(ks).le.error) then
          !lead_rank=loc(1)-1
          !if (rank.eq.lead_rank .or. ttime.gt.1) then
          if (rank.eq.0 .or. ttime.gt.1) then
            write(6,*) 'New knot locations converged'
            close(unit=neqdsk)
            return
          else
            ! Finished with parallel processing
#if defined(USEMPI)
            call mpi_finalize(ierr)
#endif
            stop
          endif
        endif
        if ((rank.eq.0 .or. ttime.gt.1) .and. terror(ks).le.berror) then
          berror=terror(ks)
          kbest=kpos
        endif
      enddo

#if defined(USEMPI)
      ! Finished with parallel processing
      if (rank.gt.0 .and. ttime.eq.1) then
        call mpi_finalize(ierr)
        stop
      endif
#endif
      close(unit=neqdsk)

      ! Re-run the best case to use as the final output
      kerror = 0
      ! With a bit of careful debugging we should be able to avoid this...
      call data_input(ks,ktime,kerror)
      if (kerror.gt.0 .or. iconvr.lt.0) then
        write(6,*) 'Error re-reading input, aborting'
        return
      endif
      j=0
      if (kppfnc.eq.6) then
        do i=2,kppknt-1
          j=j+1
          ppknt(i)=kbest(j)
        enddo
      endif
      if (kfffnc.eq.6) then
        do i=2,kffknt-1
          j=j+1
          ffknt(i)=kbest(j)
        enddo
      endif
      if (kwwfnc.eq.6) then
        do i=2,kwwknt-1
          j=j+1
          wwknt(i)=kbest(j)
        enddo
      endif
      if (keefnc.eq.6) then
        do i=2,keeknt-1
          j=j+1
          eeknt(i)=kbest(j)
        enddo
      endif
      call set_init(ks)
      call fit(ks,kerror)
      write(6,*) 'Lowest error found is ',terror(ks)

      return
      end subroutine rand_knot
