!*******************************************************************
!>     SUBPROGRAM DESCRIPTION:
!!          autoknot minimizes chi-squared as a function of knot
!!          location
!!
!!
!!     @param ks:
!!
!!     @param lconvr:
!!
!!     @param ktime: 
!!
!!     @param kerror: error flag
!!
!*******************************************************************
      subroutine autoknot(ks,lconvr,ktime,kerror)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      external ppakfunc,ffakfunc,wwakfunc,eeakfunc
      integer*4, intent(inout) :: kerror
      kerror = 0
!
!    Store the values away so the functions can restore them 
!    after reading the k file. 
!
      ks_a = ks
      lconvr_a = lconvr
      ktime_a = ktime
      kerror_a = kerror_a
      tol = aktol
      saveiter = mxiter
      call store_autoknotvals
      mxiter_a = kakiter
!
!     Minimize chi^2 and error for pp knot location
!
      if (kppfnc .eq. 6) then
        kakpploop = kakloop
        if(kppknt .le. 3) kakpploop = 1
        do j=1,kakpploop
          do i=2,kppknt-1
            kappknt = kppknt
            kadknt = i
            dzero = appknt(i-1)
            done = appknt(i+1)
            prevknt = appknt(kadknt)
            appknt(kadknt) = fmin(dzero,done,ppakfunc,tol)
            write(6,*) 'pp knot ',kadknt,' set to ',appknt(kadknt)
          enddo
          if (abs(prevknt - appknt(kadknt)) .le. aktol) then
            write(6,*) 'Last PPknot changed by less that use tolerance'
            exit
          endif
        enddo
      endif
!
!     Minimize chi^2 and error for ff knot location
!
      if (kfffnc .eq. 6) then
        kakffloop = kakloop
        if(kffknt .le. 3) kakffloop = 1
        do j=1,kakffloop
          do i=2,kffknt-1
            kaffknt = kffknt
            kadknt = i
            dzero = affknt(i-1)
            done = affknt(i+1)
            prevknt = affknt(kadknt)
            affknt(kadknt) = fmin(dzero,done,ffakfunc,tol)
            write(6,*) 'ff knot ',kadknt,' set to ',affknt(kadknt)
          enddo
           if (abs(prevknt - affknt(kadknt)) .le. aktol) then
             write(6,*) 'Last FFknot changed by less that use tolerance'
             exit
           endif
        enddo
      endif
!
!     Minimize chi^2 and error for ww knot location
!
      if (kwwfnc .eq. 6) then
        kakwwloop = kakloop
        if(kwwknt .le. 3) kakwwloop = 1
        do j=1,kakwwloop
          do i=2,kwwknt-1
            kawwknt = kwwknt
            kadknt = i
            dzero = awwknt(i-1)
            done = awwknt(i+1)
            prevknt = awwknt(kadknt)
            awwknt(kadknt) = fmin(dzero,done,wwakfunc,tol)
            write(6,*) 'ww knot ',kadknt,' set to ',awwknt(kadknt)
          enddo
          if (abs(prevknt - awwknt(kadknt)) .le. aktol) then
            write(6,*) 'Last WWknot changed by less that use tolerance'
            exit
          endif
        enddo
      endif
!
!     Minimize chi^2 and error for ee knot location
!
      if (keefnc .eq. 6) then
        kakeeloop = kakloop
        if(keeknt .le. 3) kakeeloop = 1
        do j=1,kakeeloop
          do i=2,keeknt-1
            kaeeknt = keeknt
            kadknt = i
            dzero = aeeknt(i-1)
            done = aeeknt(i+1)
            prevknt = aeeknt(kadknt)
            aeeknt(kadknt) = fmin(dzero,done,eeakfunc,tol)
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
      call data_input(ks_a,lconvr_a,ktime_a,kerror)
      if(kerror.gt.0) return
      if(lconvr_a.lt.0) return
      mxiter_a = saveiter
      call restore_autoknotvals
      call set_init(ks_a)
      call fit(ks_a,kerror_a)
      if (kerror_a /= 0) then
        kerror = 1
        return
      endif

      return
      end subroutine autoknot
!
!    store values read from k file into autoknot variables
!
      subroutine restore_autoknotvals()
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      ppknt(1:npcurn) = appknt(1:npcurn)
      ffknt(1:npcurn) = affknt(1:npcurn)
      wwknt(1:npcurn) = awwknt(1:npcurn)
      eeknt(1:npcurn) = aeeknt(1:npcurn)
      kppknt = kappknt
      kffknt = kaffknt
      kwwknt = kawwknt
      keeknt = kaeeknt
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
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      appknt(1:npcurn) = ppknt(1:npcurn)
      affknt(1:npcurn) = ffknt(1:npcurn)
      awwknt(1:npcurn) = wwknt(1:npcurn)
      aeeknt(1:npcurn) = eeknt(1:npcurn)
      kappknt = kppknt
      kaffknt = kffknt
      kawwknt = kwwknt
      kaeeknt = keeknt
      mxiter_a = mxiter
      return
      end subroutine store_autoknotvals
!
! used by autoknot which passes the routine to a minimization subroutine
! which calls it to evaulate the function being minimized
!
      function ppakfunc(xknot)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      kerror = 0
      ppakfunc = 1000.0
      write(6,*)
      write(6,*) ' trying pp knot ',kadknt,' at location ',xknot, &
                 ' out of ',kappknt,' knots'
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
      function ffakfunc(xknot)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      kerror = 0
      ffakfunc = 1000.0
      write(6,*)
      write(6,*) ' trying ff knot ',kadknt,' at location ',xknot, &
                 ' out of ',kaffknt,' knots'
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
      function wwakfunc(xknot)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      kerror = 0
      wwakfunc = 1000.0
      write(6,*)
      write(6,*) ' trying ww knot ',kadknt,' at location ',xknot, &
                 ' out of ',kawwknt,' knots'
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
      function eeakfunc(xknot)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      kerror = 0
      eeakfunc = 1000.0
      write(6,*)
      write(6,*) ' trying ee knot ',kadknt,' at location ',xknot, &
                 ' out of ',kaeeknt,' knots'
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
