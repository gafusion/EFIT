!*******************************************************************
!!    autoknot minimizes chi-squared and grad-shfranov error
!!      as a function of knot position
!!
!!    @param ks: time index
!!    @param lconvr: convergence flag
!!    @param ktime : number of time slices
!!    @param kerror: error flag
!!
!*******************************************************************
      subroutine autoknot(ks,lconvr,ktime,kerror)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none

      real*8 fmin
      real*8, external :: ppakfunc,ffakfunc,wwakfunc,eeakfunc
      integer*4, intent(in) :: ks,ktime
      integer*4, intent(inout) :: lconvr
      integer*4, intent(out) :: kerror
      integer*4 i,j,kloop,saveiter
      real*8 lbnd,rbnd,prevknt
      if(aktol.le.0.0) &
       aktol=0.1/max(kppknt,max(kffknt,max(kwwknt,keeknt)))*akrange

      kerror = 0
!
!     Store the values away so the functions can restore them 
!     after re-reading the k file. 
!
      ks_a = ks
      lconvr_a = lconvr
      ktime_a = ktime
      kerror_a = kerror
      saveiter = mxiter
      call store_autoknotvals
      mxiter_a = kakiter
!
!     Minimize chi^2 and error for pp knot locations
!
      if (kppfnc .eq. 6) then
        kloop = kakloop
        if(kppknt .le. 3) kloop = 1
        do j=1,kloop
          do i=2,kppknt-1
            kadknt = i
            lbnd = appknt(i)-akrange*(appknt(i)-appknt(i-1))
            rbnd = appknt(i)+akrange*(appknt(i+1)-appknt(i))
            prevknt = appknt(kadknt)
            appknt(kadknt) = fmin(lbnd,rbnd,ppakfunc,aktol)
            write(6,*) 'pp knot ',kadknt,' set to ',appknt(kadknt)
          enddo
          if (abs(prevknt - appknt(kadknt)) .le. aktol) then
            write(6,*) 'Last PPknot changed by less that use tolerance'
            exit
          endif
        enddo
      endif
!
!     Minimize chi^2 and error for ff knot locations
!
      if (kfffnc .eq. 6) then
        kloop = kakloop
        if(kffknt .le. 3) kloop = 1
        do j=1,kloop
          do i=2,kffknt-1
            kadknt = i
            lbnd = affknt(i)-akrange*(affknt(i)-affknt(i-1))
            rbnd = affknt(i)+akrange*(affknt(i+1)-affknt(i))
            prevknt = affknt(kadknt)
            affknt(kadknt) = fmin(lbnd,rbnd,ffakfunc,aktol)
            write(6,*) 'ff knot ',kadknt,' set to ',affknt(kadknt)
          enddo
          if (abs(prevknt - affknt(kadknt)) .le. aktol) then
            write(6,*) 'Last FFknot changed by less that use tolerance'
            exit
          endif
        enddo
      endif
!
!     Minimize chi^2 and error for ww knot locations
!
      if (kwwfnc .eq. 6) then
        kloop = kakloop
        if(kwwknt .le. 3) kloop = 1
        do j=1,kloop
          do i=2,kwwknt-1
            kadknt = i
            lbnd = awwknt(i)-akrange*(awwknt(i)-awwknt(i-1))
            rbnd = awwknt(i)+akrange*(awwknt(i+1)-awwknt(i))
            prevknt = awwknt(kadknt)
            awwknt(kadknt) = fmin(lbnd,rbnd,wwakfunc,aktol)
            write(6,*) 'ww knot ',kadknt,' set to ',awwknt(kadknt)
          enddo
          if (abs(prevknt - awwknt(kadknt)) .le. aktol) then
            write(6,*) 'Last WWknot changed by less that use tolerance'
            exit
          endif
        enddo
      endif
!
!     Minimize chi^2 and error for ee knot locations
!
      if (keefnc .eq. 6) then
        kloop = kakloop
        if(keeknt .le. 3) kloop = 1
        do j=1,kloop
          do i=2,keeknt-1
            kadknt = i
            lbnd = aeeknt(i)-akrange*(aeeknt(i)-aeeknt(i-1))
            rbnd = aeeknt(i)+akrange*(aeeknt(i+1)-aeeknt(i))
            prevknt = aeeknt(kadknt)
            aeeknt(kadknt) = fmin(lbnd,rbnd,eeakfunc,aktol)
            write(6,*) 'ee knot ',kadknt,' set to ',aeeknt(kadknt)
          enddo
          if (abs(prevknt - aeeknt(kadknt)) .le. aktol) then
            write(6,*) 'Last EEknot changed by less that use tolerance'
            exit
          endif
        enddo
      endif
!
!     Now do the final fit with adjusted knots and full iterations
!
      call data_input(ks,lconvr,ktime,kerror)
      if(kerror.gt.0) return
      if(lconvr.lt.0) return
      mxiter_a = saveiter
      call restore_autoknotvals
      call set_init(ks)
      call fit(ks,kerror)

      return
      end subroutine autoknot

!*******************************************************************
!!    knot_opt varies knot positions to target reduce grad-shfranov
!!      error untile convergence is achieved
!!
!!    @param ks: time index
!!    @param lconvr: convergence flag
!!    @param ktime : number of time slices
!!    @param kerror: error flag
!!
!*******************************************************************
      subroutine knot_opt(ks,lconvr,ktime,kerror)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none

      real*8 fmin
      real*8, external :: ppakfunc,ffakfunc,wwakfunc,eeakfunc
      integer*4, intent(in) :: ks,ktime
      integer*4, intent(inout) :: lconvr
      integer*4, intent(out) :: kerror
      integer*4 i,j,kloop,saveiter
      real*8 lbnd,rbnd

      kerror = 0

      ! Store the values away so the functions can restore them 
      !  after re-reading the k file. 
      ks_a = ks
      lconvr_a = lconvr
      ktime_a = ktime
      kerror_a = kerror
      saveiter = mxiter
      call store_autoknotvals
      mxiter_a = kakiter
      if(aktol.le.0.0) &
       aktol=0.1/max(kppknt,max(kffknt,max(kwwknt,keeknt)))*akrange

      ! Run with input settings first
      if(lconvr.lt.0) return
      call set_init(ks)
      call fit(ks,kerror)
      if ((kerror == 0) .and. (terror(ks).le.error)) then
        write(6,*) 'Input settings converged'
        return
      endif

      ! Minimize chi^2 and error for pp knot locations
      if (kppfnc .eq. 6) then
        kloop = kakloop
        if(kppknt .le. 3) kloop = 1
        do j=1,kloop
          do i=2,kppknt-1
            kadknt = i
            lbnd = appknt(i)-akrange*(appknt(i)-appknt(i-1))
            rbnd = appknt(i)+akrange*(appknt(i+1)-appknt(i))
            appknt(i) = fmin(lbnd,rbnd,ppakfunc,aktol)
            write(6,*) 'pp knot ',i,' set to ',appknt(i)
            if ((kerror == 0) .and. (terror(ks).le.error)) then
              write(6,*) 'New pp knot location converged'
              return
            endif
          enddo
        enddo
      endif

      ! Minimize chi^2 and error for ff knot locations
      if (kfffnc .eq. 6) then
        kloop = kakloop
        if(kffknt .le. 3) kloop = 1
        do j=1,kloop
          do i=2,kffknt-1
            kadknt = i
            lbnd = affknt(i)-akrange*(affknt(i)-affknt(i-1))
            rbnd = affknt(i)+akrange*(affknt(i+1)-affknt(i))
            affknt(i) = fmin(lbnd,rbnd,ffakfunc,aktol)
            write(6,*) 'ff knot ',i,' set to ',affknt(i)
            if ((kerror == 0) .and. (terror(ks).le.error)) then
              write(6,*) 'New ff knot location converged'
              return
            endif
          enddo
        enddo
      endif

      ! Minimize chi^2 and error for ww knot locations
      if (kwwfnc .eq. 6) then
        kloop = kakloop
        if(kwwknt .le. 3) kloop = 1
        do j=1,kloop
          do i=2,kwwknt-1
            kadknt = i
            lbnd = awwknt(i)-akrange*(awwknt(i)-awwknt(i-1))
            rbnd = awwknt(i)+akrange*(awwknt(i+1)-awwknt(i))
            awwknt(i) = fmin(lbnd,rbnd,wwakfunc,aktol)
            write(6,*) 'ww knot ',i,' set to ',awwknt(i)
            if ((kerror == 0) .and. (terror(ks).le.error)) then
              write(6,*) 'New ww knot location converged'
              return
            endif
          enddo
        enddo
      endif

      ! Minimize chi^2 and error for ee knot locations
      if (keefnc .eq. 6) then
        kloop = kakloop
        if(keeknt .le. 3) kloop = 1
        do j=1,kloop
          do i=2,keeknt-1
            kadknt = i
            lbnd = aeeknt(i)-akrange*(aeeknt(i)-aeeknt(i-1))
            rbnd = aeeknt(i)+akrange*(aeeknt(i+1)-aeeknt(i))
            aeeknt(i) = fmin(lbnd,rbnd,eeakfunc,aktol)
            write(6,*) 'ee knot ',i,' set to ',aeeknt(i)
            if ((kerror == 0) .and. (terror(ks).le.error)) then
              write(6,*) 'New ee knot location converged'
              return
            endif
          enddo
        enddo
      endif

      return
      end subroutine knot_opt
!
!    store values read from k file into autoknot variables
!
      subroutine restore_autoknotvals()
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none

      ppknt = appknt
      ffknt = affknt
      wwknt = awwknt
      eeknt = aeeknt
      mxiter = mxiter_a
      return
      end subroutine restore_autoknotvals
!
!     store autoknot variables into standard efit names
!     for example knot locations
!
      subroutine store_autoknotvals()
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none

      appknt = ppknt
      affknt = ffknt
      awwknt = wwknt
      aeeknt = eeknt
      mxiter_a = mxiter
      return
      end subroutine store_autoknotvals
!
! used by autoknot which passes the routine to a minimization subroutine
! which calls it to evaulate the function being minimized
!
      real*8 function ppakfunc(xknot)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      real*8, intent(in) :: xknot
      integer*4 kerror

      kerror = 0
      ppakfunc = 1000.0
      write(6,*)
      write(6,*) ' trying pp knot ',kadknt,' at location ',xknot, &
                 ' out of ',kppknt,' knots'
      call data_input(ks_a,lconvr_a,ktime_a,kerror)
      if(kerror.gt.0) return
      if(lconvr_a.lt.0) return
      call restore_autoknotvals
      ppknt(kadknt) = xknot
      call set_init(ks_a)
      call fit(ks_a,kerror_a)
      if(kerror_a .gt. 0) return
      ppakfunc = akchiwt * chisq(ks_a) + akerrwt * errorm &
               + akgamwt * chigamt + akprewt * chipre
      return
      end function ppakfunc
!
! used by autoknot which passes the routine to a minimization subroutine
! which calls it to evaulate the function being minimized
!
      real*8 function ffakfunc(xknot)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      real*8, intent(in) :: xknot
      integer*4 kerror

      kerror = 0
      ffakfunc = 1000.0
      write(6,*)
      write(6,*) ' trying ff knot ',kadknt,' at location ',xknot, &
                 ' out of ',kffknt,' knots'
      call data_input(ks_a,lconvr_a,ktime_a,kerror)
      if(kerror.gt.0) return
      if(lconvr_a.lt.0) return
      call restore_autoknotvals
      ffknt(kadknt) = xknot
      call set_init(ks_a)
      call fit(ks_a,kerror_a)
      if(kerror_a .gt. 0) return
      ffakfunc = akchiwt * chisq(ks_a) + akerrwt * errorm &
               + akgamwt * chigamt + akprewt * chipre
      return
      end function ffakfunc
!
! used by autoknot which passes the routine to a minimization subroutine
! which calls it to evaulate the function being minimized
!
      real*8 function wwakfunc(xknot)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      real*8, intent(in) :: xknot
      integer*4 kerror

      kerror = 0
      wwakfunc = 1000.0
      write(6,*)
      write(6,*) ' trying ww knot ',kadknt,' at location ',xknot, &
                 ' out of ',kwwknt,' knots'
      call data_input(ks_a,lconvr_a,ktime_a,kerror)
      if(kerror.gt.0) return
      if(lconvr_a.lt.0) return
      call restore_autoknotvals
      wwknt(kadknt) = xknot
      call set_init(ks_a)
      call fit(ks_a,kerror_a)
      if(kerror_a .gt. 0) return
      wwakfunc = akchiwt * chisq(ks_a) + akerrwt * errorm &
               + akgamwt * chigamt + akprewt * chipre
      return
      end function wwakfunc
!
! used by autoknot which passes the routine to a minimization subroutine
! which calls it to evaulate the function being minimized
!
      real*8 function eeakfunc(xknot)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      real*8, intent(in) :: xknot
      integer*4 kerror

      kerror = 0
      eeakfunc = 1000.0
      write(6,*)
      write(6,*) ' trying ee knot ',kadknt,' at location ',xknot, &
                 ' out of ',keeknt,' knots'
      call data_input(ks_a,lconvr_a,ktime_a,kerror)
      if(kerror.gt.0) return
      if(lconvr_a.lt.0) return
      call restore_autoknotvals
      eeknt(kadknt) = xknot
      call set_init(ks_a)
      call fit(ks_a,kerror_a)
      if(kerror_a .gt. 0) return
      eeakfunc = akchiwt * chisq(ks_a) + akerrwt * errorm &
               + akgamwt * chigamt + akprewt * chipre
      return
      end function eeakfunc
