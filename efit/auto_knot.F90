!*******************************************************************
!>  
!!    autoknot minimizes chi-squared and grad-shfranov error
!!      as a function of knot position
!!
!!    @param ktime : Number of time slices
!!    @param jtime : Time index
!!    @param kerror: Error Flag
!!
!*******************************************************************
      subroutine autoknot(ks,ktime,kerror)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none

      real*8 fmin
      real*8, external :: ppakfunc,ffakfunc,wwakfunc,eeakfunc
      integer*4, intent(in) :: ks,ktime
      integer*4, intent(out) :: kerror
      integer*4 i,j,kloop,saveiter
      real*8 lbnd,rbnd,prevknt
      if(aktol.le.0.0) &
       aktol=0.1/max(kppknt,max(kffknt,max(kwwknt,keeknt))) &
                *min(minval(appdf),min(minval(affdf),min(minval(awwdf), &
                                                         minval(aeedf))))

      kerror = 0
!
!     Store the values away so the functions can restore them 
!     after re-reading the k file. 
      ks_a = ks
      ktime_a = ktime
      kerror_a = kerror
      saveiter = mxiter
      call store_autoknotvals
      mxiter_a = kakiter
!
!     Minimize chi^2 and error for pp knot locations
      if (kppfnc .eq. 6) then
        kloop = kakloop
        if(kppknt .le. 3) kloop = 1
        do j=1,kloop
          do i=2,kppknt-1
            kadknt = i
            lbnd = appknt(i)-appdf(i)*(appknt(i)-appknt(i-1))
            rbnd = appknt(i)+appdf(i)*(appknt(i+1)-appknt(i))
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
      if (kfffnc .eq. 6) then
        kloop = kakloop
        if(kffknt .le. 3) kloop = 1
        do j=1,kloop
          do i=2,kffknt-1
            kadknt = i
            lbnd = affknt(i)-affdf(i)*(affknt(i)-affknt(i-1))
            rbnd = affknt(i)+affdf(i)*(affknt(i+1)-affknt(i))
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
      if (kwwfnc .eq. 6) then
        kloop = kakloop
        if(kwwknt .le. 3) kloop = 1
        do j=1,kloop
          do i=2,kwwknt-1
            kadknt = i
            lbnd = awwknt(i)-awwdf(i)*(awwknt(i)-awwknt(i-1))
            rbnd = awwknt(i)+awwdf(i)*(awwknt(i+1)-awwknt(i))
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
      if (keefnc .eq. 6) then
        kloop = kakloop
        if(keeknt .le. 3) kloop = 1
        do j=1,kloop
          do i=2,keeknt-1
            kadknt = i
            lbnd = aeeknt(i)-aeedf(i)*(aeeknt(i)-aeeknt(i-1))
            rbnd = aeeknt(i)+aeedf(i)*(aeeknt(i+1)-aeeknt(i))
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
      call data_input(ks,ktime,kerror)
      if(kerror.gt.0) return
      if(iconvr.lt.0) return
      mxiter_a = saveiter
      call restore_autoknotvals
      call set_init(ks)
      call fit(ks,kerror)

      return
      end subroutine autoknot

!*******************************************************************
!>  
!!    knot_opt varies knot positions to target reduce grad-shfranov
!!      error until convergence is achieved
!!
!!    @param ktime : Number of time slices
!!    @param jtime : Time index
!!    @param kerror: Error Flag
!!
!*******************************************************************
      subroutine knot_opt(ks,ktime,kerror)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none

      real*8 fmin_opt
      real*8, external :: ppakfunc,ffakfunc,wwakfunc,eeakfunc
      integer*4, intent(in) :: ks,ktime
      integer*4, intent(out) :: kerror
      integer*4 i,j,kloop,saveiter
      real*8 lbnd,rbnd

      kerror = 0

      ! Store the values away so the functions can restore them 
      !  after re-reading the k file. 
      ks_a = ks
      ktime_a = ktime
      kerror_a = kerror
      saveiter = mxiter
      call store_autoknotvals
      mxiter_a = kakiter
      if(aktol.le.0.0) &
       aktol=0.1/max(kppknt,max(kffknt,max(kwwknt,keeknt))) &
                *min(minval(appdf),min(minval(affdf),min(minval(awwdf), &
                                                         minval(aeedf))))

      ! Run with input settings first
      if(iconvr.lt.0) return
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
            lbnd = appknt(i)-appdf(i)*(appknt(i)-appknt(i-1))
            rbnd = appknt(i)+appdf(i)*(appknt(i+1)-appknt(i))
            appknt(i) = fmin_opt(lbnd,rbnd,ppakfunc,aktol,error)
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
            lbnd = affknt(i)-affdf(i)*(affknt(i)-affknt(i-1))
            rbnd = affknt(i)+affdf(i)*(affknt(i+1)-affknt(i))
            affknt(i) = fmin_opt(lbnd,rbnd,ffakfunc,aktol,error)
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
            lbnd = awwknt(i)-awwdf(i)*(awwknt(i)-awwknt(i-1))
            rbnd = awwknt(i)+awwdf(i)*(awwknt(i+1)-awwknt(i))
            awwknt(i) = fmin_opt(lbnd,rbnd,wwakfunc,aktol,error)
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
            lbnd = aeeknt(i)-aeedf(i)*(aeeknt(i)-aeeknt(i-1))
            rbnd = aeeknt(i)+aeedf(i)*(aeeknt(i+1)-aeeknt(i))
            aeeknt(i) = fmin_opt(lbnd,rbnd,eeakfunc,aktol,error)
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

!*******************************************************************
!
!    this subroutine stores the knot locations and iteration
!       limit read from a k file into autoknot variables
!
!*******************************************************************
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

!*******************************************************************
!
!     this subroutine replaces the efit knot locations and iteration
!       limit with autoknot values
!
!*******************************************************************
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

!*******************************************************************
!
!     function used by autoknot to update a P' knot postion, call
!       fit, and evaluate the quality function within the fmin calls
!
!*******************************************************************
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
      call data_input(ks_a,ktime_a,kerror)
      if(kerror.gt.0) return
      if(iconvr.lt.0) return
      call restore_autoknotvals
      ppknt(kadknt) = xknot
      call set_init(ks_a)
      call fit(ks_a,kerror_a)
      if(kerror_a .gt. 0) return
      ppakfunc = akchiwt * chisq(ks_a) + akerrwt * errorm &
               + akgamwt * chigamt + akprewt * chipre
      return
      end function ppakfunc

!*******************************************************************
!
!     function used by autoknot to update an FF' knot postion, call
!       fit, and evaluate the quality function within the fmin calls
!
!*******************************************************************
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
      call data_input(ks_a,ktime_a,kerror)
      if(kerror.gt.0) return
      if(iconvr.lt.0) return
      call restore_autoknotvals
      ffknt(kadknt) = xknot
      call set_init(ks_a)
      call fit(ks_a,kerror_a)
      if(kerror_a .gt. 0) return
      ffakfunc = akchiwt * chisq(ks_a) + akerrwt * errorm &
               + akgamwt * chigamt + akprewt * chipre
      return
      end function ffakfunc

!*******************************************************************
!
!     function used by autoknot to update an omega knot postion, call
!       fit, and evaluate the quality function within the fmin calls
!
!*******************************************************************
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
      call data_input(ks_a,ktime_a,kerror)
      if(kerror.gt.0) return
      if(iconvr.lt.0) return
      call restore_autoknotvals
      wwknt(kadknt) = xknot
      call set_init(ks_a)
      call fit(ks_a,kerror_a)
      if(kerror_a .gt. 0) return
      wwakfunc = akchiwt * chisq(ks_a) + akerrwt * errorm &
               + akgamwt * chigamt + akprewt * chipre
      return
      end function wwakfunc

!*******************************************************************
!
!     function used by autoknot to update an Er knot postion, call
!       fit, and evaluate the quality function within the fmin calls
!
!*******************************************************************
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
      call data_input(ks_a,ktime_a,kerror)
      if(kerror.gt.0) return
      if(iconvr.lt.0) return
      call restore_autoknotvals
      eeknt(kadknt) = xknot
      call set_init(ks_a)
      call fit(ks_a,kerror_a)
      if(kerror_a .gt. 0) return
      eeakfunc = akchiwt * chisq(ks_a) + akerrwt * errorm &
               + akgamwt * chigamt + akprewt * chipre
      return
      end function eeakfunc

!*******************************************************************
!!
!!    Optimized fmin function with early stopping criteria.
!!    This is based (copied) from fmin in netlib_lite/r8slatec/fmin.f90
!!      (see there for a more detailed description of the algorithm)
!!    An approximation x to the point where f attains a minimum on
!!      the interval (ax,bx) is determined.
!!
!!    @param ax: left endpoint of initial interval
!!    @param bx: right endpoint of initial interval
!!    @param f: function subprogram which evaluates f(x) for any x
!!              in the interval (ax,bx)
!!    @param tol: desired length of the interval of uncertainty of
!!                the final result (.ge. 0.0d0)
!!    @param f_tol: threshold value of f 
!!    @param fmin_opt: abcissa approximating the point where f
!!                     attains a minimum 
!!
!*******************************************************************
      double precision function fmin_opt(ax,bx,f,tol,f_tol)
      double precision ax,bx,f,tol,f_tol
      double precision a,b,d,e,eps,xm,p,q,r,tol1,tol2,u,v,w
      double precision fu,fv,fw,fx,x
      double precision dabs,dsqrt,dsign
!
!     c is the squared inverse of the golden ratio
      double precision, parameter :: c = 0.5d0*(3. - dsqrt(5.0d0))
!
!     eps is approximately the square root of the relative machine
!       precision.
      eps = 1.0d00
      eps = eps/2.0d00
      tol1 = 1.0d0 + eps
      do while (tol1 .gt. 1.0d00)
        eps = eps/2.0d00
        tol1 = 1.0d0 + eps
      enddo
      eps = dsqrt(eps)
!
!     initialization
      a = ax
      b = bx
      v = a + c*(b - a)
      w = v
      x = v
      e = 0.0d0
      fx = f(x)
      fv = fx
      fw = fx
!
!     exit if f(x) is below the threshold
      if (fx .lt. f_tol) then
        fmin_opt = x
        return
      endif
!
!     setup  stopping criterion
      xm = 0.5d0*(a + b)
      tol1 = eps*dabs(x) + tol/3.0d0
      tol2 = 2.0d0*tol1
!
!     main loop starts here
      do while (dabs(x - xm) .gt. (tol2 - 0.5d0*(b - a)))
        xm = 0.5d0*(a + b)
        tol1 = eps*dabs(x) + tol/3.0d0
        tol2 = 2.0d0*tol1
!
!       is parabola fit necessary
        if (dabs(e) .gt. tol1) then
!
!         fit parabola
          r = (x - w)*(fx - fv)
          q = (x - v)*(fx - fw)
          p = (x - v)*q - (x - w)*r
          q = 2.0d00*(q - r)
          if(q .gt. 0.0d0) p = -p
          q =  dabs(q)
          r = e
          e = d
!
!         is parabola acceptable
          if ((dabs(p) .ge. dabs(0.5d0*q*r)) .or. (p .le. q*(a - x)) &
              .or. (p .ge. q*(b - x))) then
!
!           a golden-section step
            if(x .ge. xm) e = a - x
            if(x .lt. xm) e = b - x
            d = c*e
          else
!
!           a parabolic interpolation step
            d = p/q
            u = x + d
!
!           f must not be evaluated too close to ax or bx
            if((u - a) .lt. tol2) d = dsign(tol1, xm - x)
            if((b - u) .lt. tol2) d = dsign(tol1, xm - x)
          endif
        else
!
!         a golden-section step
          if(x .ge. xm) e = a - x
          if(x .lt. xm) e = b - x
          d = c*e
        endif
!
!       f must not be evaluated too close to x
        if(dabs(d) .ge. tol1) u = x + d
        if(dabs(d) .lt. tol1) u = x + dsign(tol1, d)
        fu = f(u)
!
!       exit if f(u) is below the threshold
        if (fu .le. f_tol) then
          fmin_opt = u
          return
        endif
!
!       update  a, b, v, w, and x
        if (fu .le. fx) then
          if(u .ge. x) a = x
          if(u .lt. x) b = x
          v = w
          fv = fw
          w = x
          fw = fx
          x = u
          fx = fu
          cycle
        endif
        if(u .lt. x) a = u
        if(u .ge. x) b = u
        if ((fu .le. fw) .or. (w .eq. x)) then
          v = w
          fv = fw
          w = u
          fw = fu
          cycle
        endif
        if ((fu .le. fv) .or. (v .eq. x) .or. (v .eq. w)) then
          v = u
          fv = fu
        endif
      enddo
!
!     end of main loop
      fmin_opt = x
      return
      end function fmin_opt
