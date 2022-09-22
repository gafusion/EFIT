!**********************************************************************
!>
!!    Bicubic spline routines.
!!      Put together with routines from E.Solano.
!!
!!
!!    @param bkx : interval coefficients of length lubicx+1 from
!!                 sets2d.
!!
!!    @param lx : number of terms in bkx  from sets2d.
!!
!!    @param bky : interval coefficients of length lubicy+1 from
!!                 sets2d.
!!
!!    @param ly : number of terms in bky  from sets2d.
!!
!!    @param cs : array of spline coefficients of dimension (kubicx,
!!                lubicx,kubicy,lubicy) from sets2d.
!!
!!    @param xl : the point at which interpolations are desired.
!!
!!    @param yl : the point at which interpolations are desired.
!!
!!    @param fs : vector containing results depending on icalc\n
!!                  icalc              fs\n
!!                    1                f\n
!!                    2                fx\n
!!                    3                fy\n
!!                    4                fxy\n
!!                    5                fxx\n
!!                    6                fyy
!!
!!    @param ier : error flag
!!
!!    @param icalc : flag for return variable fs
!!
!**********************************************************************
      subroutine seva2d(bkx,lx,bky,ly,cs,xl,yl,fs,ier,icalc)
      use error_control, only: errctrl_msg
      include 'eparm.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!
      INTEGER*4 ier,lx,ly
      REAL*8, intent(in) :: xl,yl
      REAL*8 cs(kubicx,lubicx,kubicy,lubicy),fs(6), &
             bkx(lubicx+1),bky(lubicy+1)
!
!     Local Variable Specifications:
!
      dimension work0(4),work1(4),work2(4)
      data n00/0/,n11/1/,n22/2/
!
!     Check for consistent variable size inputs
!
      ier = 0
      if ((lx.gt.lubicx+1).or.(ly.gt.lubicy+1)) then
        ier = 1
        call errctrl_msg('seva2d', &
                         'input variables are not consistent')
        stop
      endif
!
!     Evaluate function and its partial derivatives at (XL, YL):
!
!     First do all the lookup and interpolation stuff.
!     This is the most time consuming part of the evaluation, so
!     don't do more than needed.
!
      call interv(bky(1:ly),ly,yl,lef,mflag)
      call interv(bkx(1:lx),lx,xl,ibk,ndummy)
      h = xl - bkx(ibk)
      do jj=1,4
         work0(jj) = ppvalw(cs(1,ibk,jj,lef),h,n00)
         if (icalc.eq.1) cycle
         work1(jj) = ppvalw(cs(1,ibk,jj,lef),h,n11)
         if (icalc.le.4) cycle
         work2(jj) = ppvalw(cs(1,ibk,jj,lef),h,n22)
      enddo
      h = yl - bky(lef)
      fs(1) = ppvalw(work0,h,n00)
      if(icalc.eq.1) return
      fs(2) = ppvalw(work1,h,n00)
      if(icalc.eq.2) return
      fs(3) = ppvalw(work0,h,n11)
      if(icalc.eq.3) return
      fs(4) = ppvalw(work1,h,n11)
      if(icalc.eq.4) return
      fs(5) = ppvalw(work2,h,n00)
      if(icalc.eq.5) return
      fs(6) = ppvalw(work0,h,n22)
!
      return
      end subroutine seva2d

!**********************************************************************
!>
!!    Bicubic spline routines.
!!      Put together with routines from E.Solano.
!!
!!    Inputs:
!!
!!    @param nx : Size of first direction
!!
!!    @param ny : Size of second direction
!!
!!    @param s : nx by ny array containing the function values at (x,y).
!!                This is a 1-d array, k=k=(i-1)*ny+j.
!! 
!!    @param x : (x,y) location, arrays of length nx and ny.
!!
!!    @param y : (x,y) location, arrays of length nx and ny.
!!
!!    Outputs:
!!
!!    @param cs : array of spline coefficients of dimension (kubicx,
!!                lubicx,kubicy,lubicy).
!!
!!    @param bkx : interval coefficients of length lubicx+1.
!!
!!    @param bky : interval coefficients of length lubicy+1.
!!
!!    @param lx : number of terms in bkx.
!!
!!    @param ly : number of terms in bky.
!!
!!    @param ier : error parameter.
!!
!!    Work arrays:
!!
!!    @param wk : of dimension at least nx by ny.
!!
!**********************************************************************
      subroutine sets2d(s,cs,x,nx,bkx,lx,y,ny,bky,ly,wk,ier)
!!      use commonblocks,only: bkx,bky
      include 'eparm.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!
      parameter (krord=4,kzord=4)
      dimension s(nx*ny), x(nx), y(ny), wk(nx,ny), &
                cs(kubicx, lubicx, kubicy, lubicy), &
                bkx(lubicx+1), bky(lubicy+1)
      real*8,allocatable :: xknot(:),yknot(:),rknot(:), &
           rgrid(:),zgrid(:),zknot(:),copynew(:,:)

      ier = 0
      ALLOCATE(xknot(kubicx + nw),yknot(kubicy + nh), &
               rknot(nw+krord),rgrid(nw),zgrid(nh), &
               zknot(nh+kzord),copynew(nw,nh))
!
!     Set up knots:
!     
      call eknot(nx, x, kubicx, xknot)
      call eknot(ny, y, kubicy, yknot)
!
!     Save the original, use the work array
!
      do i=1,nx
         do j=1,ny
            k=(i-1)*ny+j
            wk(i,j) = s(k)
         enddo
      enddo
!
!     Calculate spline coefficients:
!

      call spl2bc(x, y, xknot, yknot, wk)
!
!     Coefficients stored in bkx, bky, and c:

      call spl2pp(xknot, yknot, wk, bkx, lx, bky, ly, cs)
!
      DEALLOCATE(xknot,yknot,rknot,rgrid,zgrid,zknot,copynew)
!
      return
      end subroutine sets2d

!**********************************************************************
!>
!!    calculates the b-spline coeficients
!!    
!!
!!    @param rgrid :
!!
!!    @param zgrid :
!!
!!    @param rknot :
!!
!!    @param zknot :
!!
!!    @param copynew :
!!
!**********************************************************************
      subroutine spl2bc(rgrid,zgrid,rknot,zknot,copynew)
      use eparm,only:nw,nh
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      parameter (krord=4,kzord=4)
      dimension rgrid(*),zgrid(*),rknot(*),zknot(*),copynew(*)
      real*8,allocatable :: work1(:,:),work2(:),work3(:,:)
!------------------------------------------------------------------
!--   change dimension of work2 and work3 from nw to nh          --
!--   to ensure the cases when nh > nw                           --
!------------------------------------------------------------------
!
      ALLOCATE(work1(nw,nh),work2(nh),work3(nh,2*krord-1))
!
      call spli2d(rgrid,copynew,rknot,nw,krord,nh,work2,work3,work1,iflag)
      if (iflag.ne.1) write(*,*) ' error in first spli2d, iflag=',iflag
      call spli2d(zgrid,work1,zknot,nh,kzord,nw,work2,work3,copynew,iflag)
      if (iflag.ne.1) write(*,*) ' error in second spli2d, iflag=',iflag
!
      DEALLOCATE(work1,work2,work3)
!
      return
      end subroutine spl2bc

!**********************************************************************
!>
!!    translates to pp representation
!!    
!!
!!    @param rknot :
!!
!!    @param zknot :
!!
!!    @param copy :
!!
!!    @param breakr :
!!
!!    @param lr :
!!
!!    @param breakz :
!!
!!    @param lz :
!!
!!    @param coef :
!!
!**********************************************************************
      subroutine spl2pp(rknot,zknot,copy,breakr,lr,breakz,lz,coef)
      use eparm,only:nw,nh,lr0,lz0
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      parameter (krord=4,kzord=4)
      dimension rknot(*),zknot(*),copy(*),coef(*),breakr(*),breakz(*)
      real*8,allocatable :: work4(:,:,:),work5(:,:,:), &
              work6(:,:,:,:)

      ALLOCATE(work4(krord,nw,nh),work5(nh,krord,lr0), &
         work6(kzord,kzord,nw,krord))
!
      call bspp2d(rknot,copy,nw,krord,nh,work4,breakr,work5,lr)
      ndum=lr*krord
      call bspp2d(zknot,work5,nh,kzord,ndum,work6,breakz,coef,lz)
! 
      DEALLOCATE(work4,work5,work6) 
      return
      end subroutine spl2pp

!**********************************************************************
!>
!!    given the ordered data points x(1)<...<x(n), this subroutine generates
!!    a knot sequence with not-a-knot end conditions (like BSNAK from IMSL)
!!    Some of this is discussed in de Boor(1978), page 211.
!!    
!!
!!    @param n :
!!
!!    @param x :
!!
!!    @param k :
!!
!!    @param xk :
!!
!**********************************************************************
      subroutine eknot(n,x,k,xk)
      use set_kinds
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension x(n),xk(n+k)
      INTEGER*4 kh
!
      do i=1,k
        xk(i)=x(1)
        ii=i+n
        xk(ii)= x(n)+1.e-5_dp
      enddo
      kh=k/2
      k2=kh+kh
      if (k2.eq.k) then
!       even k, place knots at data points
        xk(k+1:n)=x(k+1-kh:n-kh)
      else
!       odd k, place knots in between data points
        xk(k+1:n)=.5_dp*(x(k+1-kh:n-kh)+x(k-kh:n-1-kh))
      end if
      return
      end subroutine eknot

!**********************************************************************
!>
!!    alls bsplvb, banfac/slv
!!    this is an extended version of  splint , for the use in tensor prod-
!!    uct interpolation.
!!    
!!    spli2d  produces the b-spline coeff.s  bcoef(j,.)  of the spline of
!!    order  k  with knots  t (i), i=1,..., n + k , which takes on the
!!    value  gtau (i,j)  at  tau (i), i=1,..., n ; j=1,..., m .
!!    
!!    i n p u t
!!    tau   array of length  n , containing data point abscissae.
!!    a s s u m p t i o n . . .  tau  is strictly increasing
!!    gtau(.,j)  corresponding array of length  n , containing data point
!!    ordinates, j=1,...,m
!!    t     knot sequence, of length  n+k
!!    n     number of data points and dimension of spline space  s(k,t)
!!    k     order of spline
!!    m     number of data sets
!!    
!!    w o r k   a r e a
!!    work  a vector of length  n
!!    
!!    o u t p u t
!!    q     array of order  (n,2k-1), containing the triangular factoriz-
!!    ation of the coefficient matrix of the linear system for the b-
!!    coefficients of the spline interpolant.
!!    the b-coeffs for the interpolant of an additional data set
!!    (tau(i),htau(i)), i=1,...,n  with the same data abscissae can
!!    be obtained without going through all the calculations in this
!!    routine, simply by loading  htau  into  bcoef  and then execut-
!!    ing the    call banslv ( q, n, n, 2k-1, k, bcoef )
!!    bcoef the b-coefficients of the interpolant, of length  n
!!    iflag an INTEGER4 indicating success (= 1)  or failure (= 2)
!!    the linear system to be solved is (theoretically) invertible if
!!    and only if
!!    t(i) .lt. tau(i) .lt. tau(i+k),    all i.
!!    violation of this condition is certain to lead to  iflag = 2 .
!!    
!!    m e t h o d
!!    the i-th equation of the linear system  abcoef = b  for the b-co-
!!    effs of the interpolant enforces interpolation at  tau(i), i=1,...,n.
!!    hence,  b(i) = gtau(i), all i, and  a  is a band matrix with  2k-1
!!    bands (if it is invertible).
!!    the matrix  a  is generated row by row and stored, diagonal by di-
!!    agonal, in the  c o l u m n s  of the array  q , with the main diag-
!!    onal going into column  k .  see comments in the program below.
!!    the banded system is then solved by a call to  banfac (which con-
!!    structs the triangular factorization for  a  and stores it again in
!!    q ), followed by a call to  banslv (which then obtains the solution
!!    bcoef  by substitution).
!!    banfac  does no pivoting, since the total positivity of the matrix
!!    a  makes this unnecessary.
!!    
!!    INTEGER4 iflag,k,m,n,i,ilp1mx,j,jj,kpkm1,left,np1
!!    REAL8 bcoef(m,n),gtau(n,m),q(n,7),t(n+k),tau(n),work(n),taui
!!    
!**********************************************************************
      subroutine spli2d ( tau, gtau, t, n, k, m, work, q, bcoef, iflag )
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension bcoef(m,n),gtau(n,m),q(n,2*k-1),t(n+k),tau(n),work(n)
!
      nnn=1
      np1 = n + 1
      kpkm1 = 2*k - 1
      left = k
!
!     ***   loop over  i  to construct the  n  interpolation equations
      do i=1,n
         iindex=i
         taui = tau(iindex)
         ilp1mx = min(iindex+k,np1)
!        *** zero out all entries in row  i  of  a (in the 2k-1 bands)
         q(iindex,1:kpkm1) = 0.
!        *** find  left  in the closed interval (i,i+k-1) such that
!                t(left) .le. tau(i) .lt. t(left+1)
!        matrix is singular if this is not possible
         left = max(left,i)
         if (taui .lt. t(left)) then
            iflag = 2
            write(*,*) ' linear system in  splint  not invertible'
            return
         endif

   15    if (taui .ge. t(left+1)) then
            left = left + 1
            if(left .lt. ilp1mx) go to 15
            left = left - 1
            if (taui .gt. t(left+1)) then
               iflag = 2
               write(*,*) ' linear system in  splint  not invertible'
               return
            endif
         endif
!        *** the i-th equation enforces interpolation at taui, hence
!        a(i,j) = b(j,k,t)(taui), all j. only the  k  entries with  j =
!        left-k+1,...,left actually might be nonzero. these  k  numbers
!        are returned, in  work  (used for temp.storage here), by the
!        following
         call bsplvb( t, k, nnn, taui, left, work )
!        we therefore want  work(j) = b(left-k+j)(taui) to go into
!        a(i,left-k+j), i.e., into  q(i,left-i+j), since the i-th row of
!        a  is so stored in the i-th row of  q  that the (i,i)-entry of
!        a  goes into the  k-th  entry of  q.
         jj = left - iindex
         do j=1,k
            jj = jj+1
            q(iindex,jj) = work(j)
         enddo
      enddo
!
!     ***obtain factorization of  a  , stored again in  q.
      call banfac( q, n, n, kpkm1, k, iflag )
      if (iflag.eq.2) then
         write(*,*) ' linear system in  splint  not invertible'
         return
      endif
!     *** solve  a*bcoef = gtau  by backsubstitution
      do j=1,m
         work(1:n) = gtau(1:n,j)
         call banslv( q, n, n, kpkm1, k, work )
         bcoef(j,1:n) = work(1:n)
      enddo
      return
      end subroutine spli2d

!**********************************************************************
!>
!!    calls  bsplvb
!!    this is an extended version of  bsplpp  for use with tensor products
!!    
!!    onverts the b-representation  t, bcoef(.,j), n, k  of some spline into
!!    its pp-representation  breakpts, coef(j,.,.), l, k ; j=1, ..., m  .
!!    
!!    i n p u t
!!    t     knot sequence, of length  n+k
!!    bcoef(.,j) b-spline coefficient sequence, of length  n ;j=1,...,m
!!    n     length of  bcoef  and  dimension of spline space  s(k,t)
!!    k     order of the spline
!!    m     number of data sets
!!    
!!    w o r k   a r e a
!!    scrtch   of size  (k,k,m), needed to contain bcoeffs of a piece of
!!    the spline and its  k-1  derivatives   for each of the m sets
!!    
!!    o u t p u t
!!    breakpts breakpoint sequence, of length  l+1, contains (in increasing
!!    order) the distinct points in the sequence  t(k), ..., t(n+1)
!!    coef(mm,.,.)  array of size (k,n), with  coef(mm,i,j) = (i-1)st der-
!!    ivative of  mm-th  spline at breakpts(j) from the right, mm=1,.,m
!!    l     number of polynomial pieces which make up the spline in the
!!    interval  (t(k), t(n+1))
!!    
!!    m e t h o d
!!    for each breakpoint interval, the  k  relevant b-coeffs of the
!!    spline are found and then differenced repeatedly to get the b-coeffs
!!    of all the derivatives of the spline on that interval. the spline and
!!    its first  k-1  derivatives are then evaluated at the left end
!!    point of that interval, using  bsplvb  repeatedly to obtain the val-
!!    ues of all b-splines of the appropriate order at that point.
!!
!**********************************************************************
      subroutine bspp2d ( t, bcoef, n, k, m, scrtch, breakpts, coef, l )
      implicit integer*4 (i-n), real*8 (a-h, o-z)
        parameter (kmax=4)
      INTEGER*4 k,l,m,n,i,j,jp1,kmj,left
      dimension bcoef(n,m),breakpts(*),coef(m,k,*),scrtch(k,k,m),t(*), &
           biatx(kmax)
      REAL*8 diff,fkmj,sm
!
      n11=1
      n22=2
      l = 0
      breakpts(1) = t(k)
      do left=k,n
!        find the next nontrivial knot interval.
         if(t(left+1) .eq. t(left)) cycle
         l = l + 1
         breakpts(l+1) = t(left+1)
         if (k .le. 1) then
            coef(1:m,1,l) = bcoef(left,1:m)
            cycle
         endif
!        store the k b-spline coeff.s relevant to current knot interval
!        in  scrtch(.,1) .
         do i=1,k
            scrtch(i,1,1:m) = bcoef(left-k+i,1:m)
         enddo
!        for j=1,...,k-1, compute the k-j b-spline coeff.s relevant to
!        current knot interval for the j-th derivative by differencing
!        those for the (j-1)st derivative, and store in scrtch(.,j+1) .
         do jp1=2,k
            j = jp1 - 1
            kmj = k - j
            fkmj = kmj
            do i=1,kmj
               diff = (t(left+i) - t(left+i - kmj))/fkmj
               if(diff .le. 0.) cycle
               scrtch(i,jp1,1:m) = &
                (scrtch(i+1,j,1:m) - scrtch(i,j,1:m))/diff
            enddo
         enddo
!        starting with the one b-spline of order 1 not zero at t(left),
!        find the values at t(left) of the j+1 b-splines of order j+1
!        not identically zero there from those of order j, then combine
!        with the b-spline coeff.s found earlier to compute the (k-j)-
!        th derivative at t(left) of the given spline.
         call bsplvb( t, n11, n11, t(left), left, biatx )
         coef(1:m,k,l) = scrtch(1,k,1:m)
         do jp1=2,k
            call bsplvb( t, jp1, n22, t(left), left, biatx )
            kmj = k+1 - jp1
            do mm=1,m
               coef(mm,kmj,l) = sum(biatx(1:jp1)*scrtch(1:jp1,kmj,mm))
           enddo
         enddo
      enddo
      return
      end subroutine bspp2d

!**********************************************************************
!>
!!    Calculates the value of all possibly nonzero b-splines at  x  of order
!!    
!!    jout  =  max( jhigh , (j+1)(index-1) )
!!    
!!    with knot sequence  t .
!!    
!!    i n p u t
!!    t.....knot sequence, of length  left + jout  , assumed to be nonde-
!!    creasing.  a s s u m p t i o n . . . .
!!    t(left)  .lt.  t(left + 1)   .
!!    d i v i s i o n  b y  z e r o  will result if  t(left) = t(left+1)
!!    jhigh,
!!    index.....integers which determine the order  jout = max(jhigh,
!!    (j+1)(index-1))  of the b-splines whose values at  x  are to
!!    be returned.  index  is used to avoid recalculations when seve-
!!    ral columns of the triangular array of b-spline values are nee-
!!    ded (e.g., in  bvalue  or in  bsplvd ). precisely,
!!    if  index = 1 ,
!!    the calculation starts from scratch and the entire triangular
!!    array of b-spline values of orders 1,2,...,jhigh  is generated
!!    order by order , i.e., column by column .
!!    if  index = 2 ,
!!    only the b-spline values of order  j+1, j+2, ..., jout  are ge-
!!    nerated, the assumption being that  biatx , j , deltal , deltar
!!    are, on entry, as they were on exit at the previous call.
!!    in particular, if  jhigh = 0, then  jout = j+1, i.e., just
!!    the next column of b-spline values is generated.
!!    
!!    w a r n i n g . . .  the restriction   jout .le. jmax (= 20)  is im-
!!    posed arbitrarily by the dimension statement for  deltal  and
!!    deltar  below, but is  n o w h e r e  c h e c k e d  for .
!!    
!!    x.....the point at which the b-splines are to be evaluated.
!!    left.....an INTEGER4 chosen (usually) so that
!!    t(left) .le. x .le. t(left+1)  .
!!    
!!    o u t p u t
!!    biatx.....array of length  jout , with  biatx(i)  containing the val-
!!    ue at  x  of the polynomial of order  jout  which agrees with
!!    the b-spline  b(left-jout+i,jout,t)  on the interval (t(left),
!!    t(left+1)) .
!!    
!!
!**********************************************************************
      subroutine bsplvb( t, jhigh, index, x, left, biatx )
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      parameter(jmax = 4)
      INTEGER*4 index,jhigh,left,   i,j,jp1
      REAL*8 x,saved,term
!      REAL*8 biatx(jhigh),t(1),x,
      dimension deltal(jmax),deltar(jmax)
      dimension biatx(jhigh), t(left+jhigh)
!urrent fortran standard makes it impossible to specify the length of
!  t  and of  biatx  precisely without the introduction of otherwise
!  superfluous additional arguments.
      data j/1/
      save deltal,deltar  ! (valid in fortran 77)
!
      if (index .ne. 2) then
         j = 1
         biatx(1) = 1.
         if(j .ge. jhigh) return
      endif
!
   20 continue 
      jp1 = j + 1
      deltar(j) = t(left+j) - x
      deltal(j) = x - t(left+1-j)
      saved = 0.
      do i=1,j
         term = biatx(i)/(deltar(i) + deltal(jp1-i))
         biatx(i) = saved + deltar(i)*term
         saved = deltal(jp1-i)*term
      enddo
      biatx(jp1) = saved
      j = jp1
      if(j .lt. jhigh) go to 20
      return
      end subroutine bsplvb

!**********************************************************************
!>
!!    Modified for optimization by S.J. Thompson, 30-Aug-1993
!!    Revised to eliminate call to interv by S.M.Wolfe, 17-Dec-1993
!!    and to use ASF's for evaluation
!!    This routine performs only the innermost guts of the spline evaluation
!!    Assumes k=4 (cubic spline only). No other cases considered.
!!    does not call  interv
!!    alculates value at  x  of  jd-th derivative of pp fct from pp-repr
!!    
!!    i n p u t   to PPVALU, on which this is based.
!!    break, coef, l, k.....forms the pp-representation of the function  f
!!    to be evaluated. specifically, the j-th derivative of  f  is
!!    given by
!!    
!!    (dj)f(x) = coef(j+1,i) + h(coef(j+2,i) + h( ... (coef(k-1,i) +
!!    + hcoef(k,i)/(k-j-1))/(k-j-2) ... )/2)/1
!!    
!!    with  h = x - break(i),  and
!!    
!!    i  =  max( 1 , max( j ;  break(j) .le. x , 1 .le. j .le. l ) ).
!!    
!!    x.....the point at which to evaluate.
!!    as used here, x is the distance from the break, not the absolute
!!    position.
!!    jd.....integer4 giving the order of the derivative to be evaluat-
!!    ed.  a s s u m e d  to be zero or positive.
!!    
!!    o u t p u t
!!    ppvalw.....the value of the (jd)-th derivative of  f  at  x.
!!    
!!    m e t h o d
!!    the interval index  i , appropriate for  x , is found through a
!!    call to  interv . the formula above for the  jd-th derivative
!!    of  f  is then evaluated (by nested multipication).
!!
!**********************************************************************
      function ppvalw(coef, x, jd )
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      REAL*8 ppvalw,x
      dimension coef(4)
!----------------------------------------------------------------------
! ASF's may be slightly more efficient than the alternative
!----------------------------------------------------------------------
      d2(xx) = coef(4)*xx + coef(3)
      d1(xx) = (coef(4)*xx/2. + coef(3))*xx + coef(2)
      d0(xx) = ((coef(4)*xx/3. + coef(3))*xx/2. +  &
                 coef(2))*xx + coef(1)
!-----------------------------------------------------------------------        
!   Derivatives of order k or higher are identically zero.
!-----------------------------------------------------------------------        
!   Evaluate jd-th derivative of i-th polynomial piece at x .
!-----------------------------------------------------------------------        
      select case (jd+1)
        case (1)
          ppvalw = d0(x) ! k = 4 , jd = 0
          return
        case (2)
          ppvalw = d1(x) ! k = 4 , jd = 1
          return
        case (3)
          ppvalw = d2(x) ! k = 4 , jd = 2
          return
      end select
      ppvalw = 0.
      write(*,*) 'Error (ppvalw): JD must be 0, 1, or 2.'
      write(*,*) 'Execution terminated.'
      return
      end function ppvalw

!**********************************************************************
!>
!!    This subroutine does ?
!!    
!!
!!    @param a :
!!
!!    @param nrow :
!!
!!    @param n :
!!
!!    @param ndiag :
!!
!!    @param middle :
!!
!!    @param b :
!!
!**********************************************************************        
      subroutine banslv( a, nrow, n, ndiag, middle, b )
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension a(nrow,ndiag),b(n)
      if (n .ne. 1) then
         ilo = middle - 1
         if (ilo .ge. 1) then
            do i=2,n
               jmax = min(i-1,ilo)
               do j=1,jmax
                  b(i) = b(i) - b(i-j)*a(i,middle-j)
               enddo
            enddo
         endif
      endif
!
      ihi = ndiag-middle
      do i=n,1,-1
         jmax = min(n-i,ihi)
         if (jmax .ge. 1) then
             do j=1,jmax
                b(i) = b(i) - b(i+j)*a(i,middle+j)
             enddo
         endif
         b(i) = b(i)/a(i,middle)
      enddo
      return
      end subroutine banslv

!**********************************************************************
!>
!!    This subroutine does ?
!!    
!!
!!    @param a :
!!
!!    @param nrow :
!!
!!    @param n :
!!
!!    @param ndiag :
!!
!!    @param middle :
!!
!!    @param iflag :
!!
!**********************************************************************    
      subroutine banfac( a, nrow, n, ndiag, middle, iflag )
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension a(nrow,ndiag)
      iflag = 1
      ilo = middle - 1
      if (ilo.lt.0) then
         iflag=2
         return
      elseif (ilo.eq.0) then
         do i=1,n
            if (a(i,1) .eq. 0.)  then
               iflag = 2
               return
            endif
         enddo
         return
      else
         ihi = ndiag - middle
         if (ihi.lt.0) then
            iflag = 2
            return
         elseif (ihi.eq.0) then
            do i=1,n
               if (a(i,middle) .eq. 0.) then
                  iflag = 2
                  return
               endif
               jmax = min(ilo,n-i)
               if (jmax .lt. 1) cycle
               do j=1,jmax
                  a(i+j,middle-j) = a(i+j,middle-j)/a(i,middle)
               enddo
            enddo
            return
         else
            do i=1,n
               diag = a(i,middle)
               if (diag .eq. 0.) then
                  iflag = 2
                  return
               endif
               jmax = min(ilo,n-i)
               if (jmax .lt. 1) cycle
               kmax = min(ihi,n-i)
               do j=1,jmax
                  mmj = middle-j
                  a(i+j,mmj) = a(i+j,mmj)/diag
                  do k=1,kmax
                     a(i+j,mmj+k) = a(i+j,mmj+k) &
                                  - a(i+j,mmj)*a(i,middle+k)
                  enddo
               enddo
            enddo
            return
         endif
      endif
      end subroutine banfac

!**********************************************************************
!>
!!    Computes  left = max( i ; 1 .le. i .le. lxt  .and.  xt(i) .le. x )  .
!!    
!!    i n p u t
!!    @param xt :  a REAL8 sequence, of length  lxt , assumed to be nondecreasing
!!    @param lxt : number of terms in the sequence  xt .
!!    @param x : the point whose location with respect to the sequence  xt  is
!!    to be determined.
!!    @param left :  
!!    @param mflag :  
!!    
!!    o u t p u t
!!    left, mflag.....both integers, whose value is
!!    
!!    1     -1      if               x .lt.  xt(1)
!!    i      0      if   xt(i)  .le. x .lt. xt(i+1)
!!    lxt     1      if  xt(lxt) .le. x
!!    
!!    in particular,  mflag = 0 is the 'usual' case.  mflag .ne. 0
!!    indicates that  x  lies outside the halfopen interval
!!    xt(1) .le. y .lt. xt(lxt) . the asymmetric treatment of the
!!    interval is due to the decision to make all pp functions cont-
!!    inuous from the right.
!!    
!!    m e t h o d
!!    the program is designed to be efficient in the common situation that
!!    it is called repeatedly, with  x  taken from an increasing or decrea-
!!    sing sequence. this will happen, e.g., when a pp function is to be
!!    graphed. the first guess for  left  is therefore taken to be the val-
!!    ue returned at the previous call and stored in the  l o c a l  varia-
!!    ble  ilo . a first check ascertains that  ilo .lt. lxt (this is nec-
!!    essary since the present call may have nothing to do with the previ-
!!    ous call). then, if  xt(ilo) .le. x .lt. xt(ilo+1), we set  left =
!!    ilo  and are done after just three comparisons.
!!    otherwise, we repeatedly double the difference  istep = ihi - ilo
!!    while also moving  ilo  and  ihi  in the direction of  x , until
!!    xt(ilo) .le. x .lt. xt(ihi) ,
!!    after which we use bisection to get, in addition, ilo+1 = ihi .
!!    left = ilo  is then returned.
!!    
!**********************************************************************      
      subroutine interv ( xt, lxt, x, left, mflag )
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      INTEGER*4 left,lxt,mflag,ihi,ilo,istep,middle
      REAL*8, intent(in) :: x
      dimension xt(lxt)
      data ilo /1/
      ihi = ilo + 1
      if (ihi .ge. lxt) then
         if (x .ge. xt(lxt)) then
            mflag = 1
            left = lxt
            return
         endif
         if (lxt .le. 1) then
            mflag = -1
            left = 1
            return
         endif
         ilo = lxt - 1
         ihi = lxt
      endif
!
      if (x .lt. xt(ihi)) then
         if (x .ge. xt(ilo)) then
            mflag = 0
            left = ilo
            return
         endif
!
!        **** now x .lt. xt(ilo) . decrease  ilo  to capture  x .
         istep = 1
   31    ihi = ilo
         ilo = ihi - istep
         if (ilo .gt. 1) then
            if(x .ge. xt(ilo)) go to 50
            istep = istep*2
            go to 31
         endif
         ilo = 1
         if (x .lt. xt(1)) then
            mflag = -1
            left = 1
            return
         endif
         go to 50
!        **** now x .ge. xt(ihi) . increase  ihi  to capture  x .
      endif
      istep = 1
   41 ilo = ihi
      ihi = ilo + istep
      if (ihi .lt. lxt) then
         if(x .lt. xt(ihi)) go to 50
         istep = istep*2
         go to 41
      endif
      if (x .ge. xt(lxt)) then
         mflag = 1
         left = lxt
         return
      endif
      ihi = lxt
!
!           **** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval.
   50 middle = (ilo + ihi)/2
      if (middle .eq. ilo) then
         mflag = 0
         left = ilo
         return
      endif
!     note. it is assumed that middle = ilo in case ihi = ilo+1 .
      if (x .lt. xt(middle)) then
         ihi = middle
         go to 50
      endif
      ilo = middle
      go to 50
      end subroutine interv

!**********************************************************************
!>
!!    These are interface routines to bridge from IMSL on the VAX to LINPACK on
!!    the Multiflow.
!!    
!!    LINV1F
!!    This routine is an interface from the IMSL call used on the VAX to the
!!    LINPACK call used on the Multiflow.  This routine performs the inversion
!!    of an N x N matrix.  This is not a general purpose routine and is specific to
!!    the EFIT code and should only be used for that code.  This routine has local
!!    arrays that are sized to correspond to the largest size of an EFIT call.  If
!!    the parameters which define dimensions for EFIT arrays should change, then
!!    the sizes of these arrays may need to change as well.
!!    
!!    Correspondence of the variables between the IMSL LINV1F routine and the
!!    LINPACK DGEFA and DGEDI routines.  See the IMSL and LINPACK documentation for
!!    further information.
!!    
!!    A contains the N x N matrix that is to be transposed.  This is the
!!    same matrix for both routines, however the returned contents are not
!!    the same.  This is not a problem as this matrix is not used again after
!!    the call is made.  This routine calls DGEFA with the input matrix A and
!!    it returns results which are input to DGEDI.
!!    N is the row dimension of A where number rows equal number columns.
!!    IA is the actual leading storage dimension of the matrix A and AINV.
!!    AINV is the resultant transposed N x N matrix.
!!    IDGT is an accuracy option used by the IMSL routine.  However there is no
!!    corresponding option for the LINPACK routines, so this is not used.
!!    WK is a scratch array which can be used by both calls.
!!    IER 129 is an error return for LINV1F, non-zero is an error for DGEFA.
!!    If DGEFA gets any error, then the error return is set to 129.
!!    
!!    IPVT is a vector of pivot indices returned by DGEFA and used by DGEDI.
!!    DET is a determinant of the original matrix but is not optioned or used.
!!
!**********************************************************************   
      SUBROUTINE LINV1F(A,N,IA,AINV,IDGT,WK,IER)
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      !
      REAL*8   A(IA,IA),AINV(IA,IA),WK(N)
      REAL*8   DET(2)
      INTEGER*4 IPVT(9)
      !
      CALL DGEFA(A,IA,N,IPVT,IER)
      IF (IER .NE. 0) THEN                   ! return if error
        IER = 129
        RETURN
      ENDIF
      nnn=1
      CALL DGEDI(A,IA,N,IPVT,DET,WK,nnn)
      ! move result to output array
      AINV(1:N,1:N) = A(1:N,1:N)
      RETURN
      END SUBROUTINE LINV1F
