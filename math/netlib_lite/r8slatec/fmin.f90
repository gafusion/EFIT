      double precision function fmin(ax,bx,f,tol)
      double precision ax,bx,f,tol
!
!      an approximation  x  to the point where  f  attains a minimum  on
!  the interval  (ax,bx)  is determined.
!
!
!  input..
!
!  ax    left endpoint of initial interval
!  bx    right endpoint of initial interval
!  f     function subprogram which evaluates  f(x)  for any  x
!        in the interval  (ax,bx)
!  tol   desired length of the interval of uncertainty of the final
!        result ( .ge. 0.0d0)
!
!
!  output..
!
!  fmin  abcissa approximating the point where  f  attains a minimum
!
!
!      the method used is a combination of  golden  section  search  and
!  successive parabolic interpolation.  convergence is never much slower
!  than  that  for  a  fibonacci search.  if  f  has a continuous second
!  derivative which is positive at the minimum (which is not  at  ax  or
!  bx),  then  convergence  is  superlinear, and usually of the order of
!  about  1.324....
!      the function  f  is never evaluated at two points closer together
!  than  eps*abs(fmin) + (tol/3), where eps is  approximately the square
!  root  of  the  relative  machine  precision.   if   f   is a unimodal
!  function and the computed values of   f   are  always  unimodal  when
!  separated by at least  eps*abs(x) + (tol/3), then  fmin  approximates
!  the abcissa of the global minimum of  f  on the interval  ax,bx  with
!  an error less than  3*eps*abs(fmin) + tol.  if   f   is not unimodal,
!  then fmin may approximate a local, but perhaps non-global, minimum to
!  the same accuracy.
!      this function subprogram is a slightly modified  version  of  the
!  algol  60 procedure  localmin  given in richard brent, algorithms for
!  minimization without derivatives, prentice - hall, inc. (1973).
!
!
      double precision  a,b,c,d,e,eps,xm,p,q,r,tol1,tol2,u,v,w
      double precision  fu,fv,fw,fx,x
      double precision  dabs,dsqrt,dsign
!
!  c is the squared inverse of the golden ratio
!
      c = 0.5d0*(3. - dsqrt(5.0d0))
!
!  eps is approximately the square root of the relative machine
!  precision.
!
      eps = 1.0d00
   10 eps = eps/2.0d00
      tol1 = 1.0d0 + eps
      if (tol1 .gt. 1.0d00) go to 10
      eps = dsqrt(eps)
!
!  initialization
!
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
!  main loop starts here
!
   20 xm = 0.5d0*(a + b)
      tol1 = eps*dabs(x) + tol/3.0d0
      tol2 = 2.0d0*tol1
!
!  check stopping criterion
!
      if (dabs(x - xm) .le. (tol2 - 0.5d0*(b - a))) go to 90
!
! is golden-section necessary
!
      if (dabs(e) .le. tol1) go to 40
!
!  fit parabola
!
      r = (x - w)*(fx - fv)
      q = (x - v)*(fx - fw)
      p = (x - v)*q - (x - w)*r
      q = 2.0d00*(q - r)
      if (q .gt. 0.0d0) p = -p
      q =  dabs(q)
      r = e
      e = d
!
!  is parabola acceptable
!
   30 if (dabs(p) .ge. dabs(0.5d0*q*r)) go to 40
      if (p .le. q*(a - x)) go to 40
      if (p .ge. q*(b - x)) go to 40
!
!  a parabolic interpolation step
!
      d = p/q
      u = x + d
!
!  f must not be evaluated too close to ax or bx
!
      if ((u - a) .lt. tol2) d = dsign(tol1, xm - x)
      if ((b - u) .lt. tol2) d = dsign(tol1, xm - x)
      go to 50
!
!  a golden-section step
!
   40 if (x .ge. xm) e = a - x
      if (x .lt. xm) e = b - x
      d = c*e
!
!  f must not be evaluated too close to x
!
   50 if (dabs(d) .ge. tol1) u = x + d
      if (dabs(d) .lt. tol1) u = x + dsign(tol1, d)
      fu = f(u)
!
!  update  a, b, v, w, and x
!
      if (fu .gt. fx) go to 60
      if (u .ge. x) a = x
      if (u .lt. x) b = x
      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
      go to 20
   60 if (u .lt. x) a = u
      if (u .ge. x) b = u
      if (fu .le. fw) go to 70
      if (w .eq. x) go to 70
      if (fu .le. fv) go to 80
      if (v .eq. x) go to 80
      if (v .eq. w) go to 80
      go to 20
   70 v = w
      fv = fw
      w = u
      fw = fu
      go to 20
   80 v = u
      fv = fu
      go to 20
!
!  end of main loop
!
   90 fmin = x
      return
      end
