!----------------------------------------------------------------------
!>
!!    the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!!    for a cubic interpolating spline
!!    
!!    s(x) = y(i) + b(i)(x-x(i)) + c(i)(x-x(i))**2 + d(i)(x-x(i))**3
!!    
!!    for  x(i) .le. x .le. x(i+1)\n
!!
!!    y(i) = s(x(i))\n
!!    b(i) = s'(x(i))\n
!!    c(i) = s''(x(i))/2\n
!!    d(i) = s'''(x(i))/6  (derivative from the right)
!!
!!    @param n : the number of data points or knots (n.ge.2)
!!    @param x : the abscissas of the knots in strictly increasing order
!!    @param y : the ordinates of the knots
!!    @param b : array of spline coefficients as defined above
!!    @param c : array of spline coefficients as defined above
!!    @param d : array of spline coefficients as defined above
!!
!----------------------------------------------------------------------
 subroutine zpline(n,x,y,b,c,d)
   implicit none
   integer*4, intent(in) :: n
   real*8, intent(in) :: x(n),y(n)
   real*8, intent(out) :: b(n),c(n),d(n)
   integer*4 nm1,ib,i
   real*8 t
 
   nm1 = n-1
   if(n .lt. 2) return
   if (n .lt. 3) then
     b(1) = (y(2)-y(1))/(x(2)-x(1))
     c(1) = 0.
     d(1) = 0.
     b(2) = b(1)
     c(2) = 0.
     d(2) = 0.
     return
   endif
!
!  set up tridiagonal system
!
!  b = diagonal, d = offdiagonal, c = right hand side
!
   d(1) = x(2) - x(1)
   c(2) = (y(2) - y(1))/d(1)
   do i=2,nm1
     d(i) = x(i+1) - x(i)
     b(i) = 2.*(d(i-1) + d(i))
     c(i+1) = (y(i+1) - y(i))/d(i)
     c(i) = c(i+1) - c(i)
   enddo
!
!  end conditions.  third derivatives at  x(1)  and  x(n)
!  obtained from divided differences
!
   b(1) = -d(1)
   b(n) = -d(n-1)
   c(1) = 0.
   c(n) = 0.
   if (n .ne. 3) then
     c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
     c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
     c(1) = c(1)*d(1)**2/(x(4)-x(1))
     c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
   endif
!
!  forward elimination
!
   do i=2,n
     t = d(i-1)/b(i-1)
     b(i) = b(i) - t*d(i-1)
     c(i) = c(i) - t*c(i-1)
   enddo
!
!  back substitution
!
   c(n) = c(n)/b(n)
   do ib=1,nm1
     i = n-ib
     c(i) = (c(i) - d(i)*c(i+1))/b(i)
   enddo
!
!  c(i) is now the sigma(i) of the text
!
!  compute polynomial coefficients
!
   b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
   do i = 1, nm1
     b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
     d(i) = (c(i+1) - c(i))/d(i)
     c(i) = 3.*c(i)
   enddo
   c(n) = 3.*c(n)
   d(n) = d(n-1)
   return
 end subroutine zpline

!**********************************************************************
!>
!!    the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!!    for a cubic interpolating spline *for which s'(xn)=0
!!    
!!    s(x) = y(i) + b(i)(x-x(i)) + c(i)(x-x(i))**2 + d(i)(x-x(i))**3
!!    
!!    for  x(i) .le. x .le. x(i+1)\n
!!
!!    y(i) = s(x(i))\n
!!    b(i) = s'(x(i))\n
!!    c(i) = s''(x(i))/2\n\n
!!    d(i) = s'''(x(i))/6  (derivative from the right)
!!    
!!    @param n : the number of data points or knots (n.ge.2)
!!    @param x : the abscissas of the knots in strictly increasing order
!!    @param y : the ordinates of the knots
!!    @param b : arrays of spline coefficients as defined above
!!    @param c : arrays of spline coefficients as defined above
!!    @param d : arrays of spline coefficients as defined above
!!
!**********************************************************************
 subroutine spleen(n,x,y,b,c,d)
   implicit none
   integer*4, intent(in) :: n
   real*8, intent(in) :: x(n),y(n)
   real*8, intent(out) :: b(n),c(n),d(n)
   integer*4 nm1,ib,i
   real*8 t

   nm1 = n-1
   if(n .lt. 2) return
   if (n .lt. 3) then
     b(1) = (y(2)-y(1))/(x(2)-x(1))
     c(1) = 0.
     d(1) = 0.
     b(2) = b(1)
     c(2) = 0.
     d(2) = 0.
     return
   endif
!
!  set up tridiagonal system
!
!  b = diagonal, d = offdiagonal, c = right hand side
!
   d(1) = x(2) - x(1)
   c(2) = (y(2) - y(1))/d(1)
   do i=2,nm1
     d(i) = x(i+1) - x(i)
     b(i) = 2.*(d(i-1) + d(i))
     c(i+1) = (y(i+1) - y(i))/d(i)
     c(i) = c(i+1) - c(i)
   enddo
!
!  end conditions.  third derivatives at  x(1)  and  x(n)
!  obtained from divided differences
!
!  *difference from zpline is here
!
   b(1) = -d(1)
   b(n) = 2.*d(n-1)
   c(1) = 0.
   c(n) = 0.
   if (n .eq. 3) then
     c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
     c(1) = c(1)*d(1)**2/(x(4)-x(1))
     c(n)= -(y(n)-y(n-1))/d(n-1)
   endif
!
!  forward elimination
!
   do i=2,n
     t = d(i-1)/b(i-1)
     b(i) = b(i) - t*d(i-1)
     c(i) = c(i) - t*c(i-1)
   enddo
!
!  back substitution
!
   c(n) = c(n)/b(n)
   do ib=1,nm1
     i = n-ib
     c(i) = (c(i) - d(i)*c(i+1))/b(i)
   enddo
!
!  c(i) is now the sigma(i) of the text
!
!  compute polynomial coefficients
!
   b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
   do i=1,nm1
     b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
     d(i) = (c(i+1) - c(i))/d(i)
     c(i) = 3.*c(i)
   enddo
   c(n) = 3.*c(n)
   d(n) = d(n-1)
   return
 end subroutine spleen

!**********************************************************************
!>
!!   the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!!   for a cubic interpolating spline *for which s'(x1)=0
!!    
!!    s(x) = y(i) + b(i)(x-x(i)) + c(i)(x-x(i))2 + d(i)(x-x(i))^3
!!    
!!    for  x(i) .le. x .le. x(i+1)\n
!!
!!    y(i) = s(x(i))\n
!!    b(i) = s'(x(i))\n
!!    c(i) = s''(x(i))/2\n\n
!!    d(i) = s'''(x(i))/6  (derivative from the right)
!!
!!    @param n : the number of data points or knots (n.ge.2)
!!    @param x : the abscissas of the knots in strictly increasing order
!!    @param y : the ordinates of the knots
!!    @param b : arrays of spline coefficients as defined above.
!!    @param c : arrays of spline coefficients as defined above.
!!    @param d : arrays of spline coefficients as defined above.
!!
!**********************************************************************
 subroutine splaan(n,x,y,b,c,d)
   implicit none
   integer*4, intent(in) :: n
   real*8, intent(in) :: x(n),y(n)
   real*8, intent(out) :: b(n),c(n),d(n)
   integer*4 nm1,ib,i
   real*8 t

   nm1 = n-1
   if(n .lt. 2) return
   if (n .lt. 3) then
     b(1) = (y(2)-y(1))/(x(2)-x(1))
     c(1) = 0.
     d(1) = 0.
     b(2) = b(1)
     c(2) = 0.
     d(2) = 0.
     return
   endif
!
!  set up tridiagonal system
!
!  b = diagonal, d = offdiagonal, c = right hand side.
!
   d(1) = x(2) - x(1)
   c(2) = (y(2) - y(1))/d(1)
   do i=2,nm1
     d(i) = x(i+1) - x(i)
     b(i) = 2.*(d(i-1) + d(i))
     c(i+1) = (y(i+1) - y(i))/d(i)
     c(i) = c(i+1) - c(i)
   enddo
!
!  end conditions.  third derivatives at  x(1)  and  x(n)
!  obtained from divided differences
!
!  *difference from zpline is here
!
   b(1)=2.*d(1)
   b(n) = -d(n-1)
   c(1) = 0.
   c(n) = 0.
   if (n .eq. 3) then
     c(1) = (y(2)-y(1))/d(1)
     c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
     c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
   endif
!
!  forward elimination
!
   do i=2,n
     t = d(i-1)/b(i-1)
     b(i) = b(i) - t*d(i-1)
     c(i) = c(i) - t*c(i-1)
   enddo
!
!  back substitution
!
   c(n) = c(n)/b(n)
   do ib=1,nm1
     i = n-ib
     c(i) = (c(i) - d(i)*c(i+1))/b(i)
   enddo
!
!  c(i) is now the sigma(i) of the text
!
!  compute polynomial coefficients
!
   b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
   do i=1,nm1
     b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
     d(i) = (c(i+1) - c(i))/d(i)
     c(i) = 3.*c(i)
   enddo
   c(n) = 3.*c(n)
   d(n) = d(n-1)
   return
 end subroutine splaan

!**********************************************************************
!>
!!    this subroutine evaluates the cubic spline function
!!
!!    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3\n
!!
!!    where  x(i) .lt. u .lt. x(i+1), using horner's rule\n
!! 
!!    if  u .lt. x(1) then  i = 1  is used.\n
!!    if  u .ge. x(n) then  i = n  is used.\n
!!
!!    if  u  is not in the same interval as the previous call, then a
!!      binary search is performed to determine the proper interval.
!!
!!    y(i) = s(x(i))\n\n
!!    b(i) = s'(x(i))\n\n
!!    c(i) = s''(x(i))/2\n\n\n
!!    d(i) = s'''(x(i))/6 (derivative from the right)
!!
!!    @param n : the number of data points or knots (n.ge.2)
!!    @param u : function to be splined
!!    @param x : the abscissas of the knots in strictly increasing order
!!    @param y : the ordinates of the knots
!!    @param b : arrays of spline coefficients as defined above.
!!    @param c : rrays of spline coefficients as defined above.
!!    @param d : rrays of spline coefficients as defined above.
!!
!**********************************************************************
 real*8 function seval(n,u,x,y,b,c,d)
   implicit none
   integer*4, intent(in) :: n
   real*8, intent(in) :: u,x(n),y(n),b(n),c(n),d(n)
   integer*4 i,j,k
   real*8 dx
   data i/1/

   if(i .ge. n) i = 1
!
!  binary search
!
   if ((u .lt. x(i)) .or. (u .gt. x(i+1))) then
     i = 1
     j = n+1
     k = (i+j)/2
     if(u .lt. x(k)) j = k
     if(u .ge. x(k)) i = k
     do while (j .gt. i+1)
       k = (i+j)/2
       if(u .lt. x(k)) j = k
       if(u .ge. x(k)) i = k
     enddo
   endif
!
!  evaluate spline
!
   dx = u - x(i)
   seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
   return
 end function seval

!**********************************************************************
!>
!!    this subroutine evaluates the *derivative of the cubic spline function
!!
!!    where  x(i) .lt. u .lt. x(i+1), using horner's rule\n
!!
!!    if  u .lt. x(1) then  i = 1  is used.\n
!!    if  u .ge. x(n) then  i = n  is used.\n
!!
!!    @param n : the number of data points or knots (n.ge.2)
!!    @param u : the function to be splined
!!    @param x : the abscissas of the knots in strictly increasing order
!!    @param y : the ordinates of the knots
!!    @param b : arrays of spline coefficients as defined above.
!!    @param c : rrays of spline coefficients as defined above.
!!    @param d : rrays of spline coefficients as defined above.
!!
!**********************************************************************
 real*8 function speval(n,u,x,y,b,c,d)
   implicit none
   integer*4, intent(in) :: n
   real*8, intent(in) :: u,x(n),y(n),b(n),c(n),d(n)
   integer*4 i,j,k
   real*8 dx
   data i/1/

   if(i .ge. n) i = 1
!
!  binary search
!
   if ((u .lt. x(i)) .or. (u .gt. x(i+1))) then
     i = 1
     j = n+1
     k = (i+j)/2
     if(u .lt. x(k)) j = k
     if(u .ge. x(k)) i = k
     do while (j .gt. i+1)
       k = (i+j)/2
       if(u .lt. x(k)) j = k
       if(u .ge. x(k)) i = k
     enddo
   endif
!
!  evaluate spline *derivative
!
   dx = u - x(i)
   speval = b(i) + dx*(2.*c(i) + 3.*dx*d(i))
   return
 end function speval
