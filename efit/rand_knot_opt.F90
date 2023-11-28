!*******************************************************************
!!
!!    rand_knot varies knot positions randomly until convergence
!!      is achieved
!!
!!    @param ks: time index
!!    @param lconvr: convergence flag
!!    @param ktime : number of time slices
!!    @param kerror: error flag
!!
!*******************************************************************
      subroutine rand_knot(ks,lconvr,ktime,kerror)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none

      integer*4, intent(in) :: ks,ktime
      integer*4, intent(inout) :: lconvr
      integer*4, intent(out) :: kerror
      integer*4 i,j,k,nvary,saveiter
      real*8 berror
      real*8, dimension(:), allocatable :: lbnd,krange,kpos,kbest

      kerror = 0

      ! Run with input settings first
      if(lconvr.lt.0) return
      call set_init(ks)
      call fit(ks,kerror)
      write(6,*) 'Initial error is: ',terror(ks)
      if ((kerror == 0) .and. (terror(ks).le.error)) then
        write(6,*) 'Input settings converged'
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
        write(6,*) 'No knots to vary'
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

      ! Vary knots randomly until convergence or max iterations is reached
      do k=1,kakloop
        kerror = 0
        ! With a bit of careful debugging we should be able to avoid this...
        call data_input(ks,lconvr,ktime,kerror)
        if (kerror.gt.0 .or. lconvr.lt.0) then
          write(6,*) 'Error re-reading input, aborting'
          return
        endif
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
        mxiter = kakiter
        call set_init(ks)
        call fit(ks,kerror)
        write(6,*) 'Iteration ',k,' has error: ',terror(ks)
        if ((kerror == 0) .and. (terror(ks).le.error)) then
          write(6,*) 'New knot locations converged'
          return
        endif
        if ((kerror == 0) .and. (terror(ks).le.berror)) then
          berror=terror(ks)
          kbest=kpos
        endif
      enddo

      ! re-run the best case to use as the final output
      kerror = 0
      ! With a bit of careful debugging we should be able to avoid this...
      call data_input(ks,lconvr,ktime,kerror)
      if (kerror.gt.0 .or. lconvr.lt.0) then
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
      call fit(ks,kerror)
      write(6,*) 'Lowest error found is ',terror(ks)

      return
      end subroutine rand_knot
