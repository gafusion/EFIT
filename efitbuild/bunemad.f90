      subroutine buneto(psi,nwb,nhb,sia)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          buneto sets up the appropriate arrays for the           **
!**          Buneman's solver.                                       **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          05 03/85..........first created                         **
!**                                                                  **
!**                                                                  **
!**********************************************************************
!vas  f90 modifi.
      use var_bunemn

      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension   psi(*), sia(*)
!vas  f90 modifi
!vas      common/bunemn/mno,m,n,s,shift,dr,dz
      mno = nbmdim
      m = nww
      n = nhh
      s = drdz2
      shift = rgrid1
      dr = delrgrid
      dz = delz


!----------------------------------------------------------------------
! copy psi into sia rowwise                                          --
!----------------------------------------------------------------------
      do 2 i = 1,nwb
      ii = (i-1)*nhb + 1
      do 2 j = 1, nhb
      sia(i+j*nwb-nwb) = psi(ii-1+j)
    2 continue
      ia = nwb+nwb
      ju = n*nwb
!-----------------------------
! set up for rzpois
!-----------------------------
      do 3 i = ia,ju,nwb
      sia(i-m+1) = sia(i-m+1)+(.5+.25/(1.+shift/dr))*sia(i-m)/s
      sia(i-1) = sia(i-1)+(.5-.25/(m-1+shift/dr))*sia(i)/s
    3 continue
      call rzpois(sia)
      nwhbb = nwb*nhb
!------------------------------------------------------------------------------
! copy sia back into psi columnwise                                          --
!------------------------------------------------------------------------------
      do 6 i = 2,n
      ii = (i-1)*nwb + 1
      do 6 j = 2,m
      psi(i+j*nhb-nhb) = sia(ii-1+j)
    6 continue
      return
      end

      subroutine rzpois(q)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          rzpois solves for the poloidal flux using the           **
!**          Buneman's method.                                       **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          05 03/85..........first created                         **
!**                                                                  **
!**                                                                  **
!**********************************************************************
!vas  f90 modifi
      use set_kinds
      use var_bunemn

!vas      implicit integer*4 (i-n), real*8 (a-h, o-z)
!vas      common/bunemn/mno,m,n,s,shift,dr,dz
!vas      dimension   g(mno),p(mno),c(mno),d(mno),temp(mno)
!vas      dimension  q(1)
!vas  f90 modifi
      implicit none
      real(rprec) flag
      real(rprec) shftdr,a,pi,as
      integer(iprec) i,l,lo,ju,n222,id,jd,ko,k4,li,jh,jt,ji,jo, &
      j2,iu,j,ii,io,iallocate_stat

      real*8 q(1)
      real(rprec),dimension(:), allocatable :: g
      real(rprec),dimension(:), allocatable :: p
      real(rprec),dimension(:), allocatable :: c
      real(rprec),dimension(:), allocatable :: d
      real(rprec),dimension(:), allocatable :: temp

      flag = 0.0
!vas  f90 modifi
      allocate(g(mno),stat=iallocate_stat)
      if(iallocate_stat/=0) stop "*** Not enough space for g ***"
      allocate(p(mno),stat=iallocate_stat)
      if(iallocate_stat/=0) stop "*** Not enough space for p ***"
      allocate(c(mno),stat=iallocate_stat)
      if(iallocate_stat/=0) stop "*** Not enough space for c ***"
      allocate(d(mno),stat=iallocate_stat)
      if(iallocate_stat/=0) stop "*** Not enough space for d ***"
      allocate(temp(mno),stat=iallocate_stat)
      if(iallocate_stat/=0) stop "*** Not enough space for temp ***"

!vas  f90 modifi
!vas      do 50 i = 1,mno
!vas      g(i) = 0.
!vas      p(i) = 0.
!vas      d(i) = 0.
!vas   50 continue
      g = 0.
      p = 0.
      d = 0.

      if (flag .ne. 0.) go to 200
!vas  f90 modifi
!vas      do 60 i = 1,mno
!vas      temp(i) = 0.
!vas      c(i) = 0.
!vas   60 continue
      temp = 0.
      c = 0.

      shftdr = shift/dr
      do 40 i = 2,m
      temp(i) = 1. - .5/(i+shftdr-1.)
   40 continue
      ju = (n-1)*(m+1)
      n222 = n/2
      c(n222) = 0.
      lo = n/2
   1  l = lo/2
      c(l) = sqrt(2.+c(lo))
      lo = l
   2  c(n-l) = -c(l)
      l = l+2*lo
      if((2*l/n)*(2*lo-3)< 0) then
        goto 4
      else if((2*l/n)*(2*lo-3) == 0) then 
        goto 3
      else 
        goto 1
      end if
   3  c(l) = (c(l+lo)+c(l-lo))/c(lo)
      go to 2
   4  do 5 l = 2,n
   5  c(l-1) = 1./(2.+s*(2.-c(l-1)))
      flag = 1.
  200 lo = n/2
      ko = 2
      id = 1
  15  li = 2*lo
      k4 = 2*ko-li/n
      jd = (m+1)*n/li
      jh = (m+1)*(n/(2*li))
      jt = jd+jh
      ji = 2*jd
      jo = jd*ko
      do 11 j = jo,ju,ji
      j2 = j+2
      iu = j+m
      select case (k4)
        case (1)
          go to 20
        case (2)
          go to 24
        case (3)
          go to 26
        case (4)
          go to 28
      end select
  28  do 29 i = j2,iu
      pi = q(i)-q(i+jt)-q(i-jt)
      q(i) = q(i)-q(i+jh)-q(i-jh)+q(i+jd)+q(i-jd)
  29  p(i-j) = pi+q(i)
      go to 10
  26  do 27 i = j2,iu
      p(i-j) = 2.*q(i)
  27  q(i) = q(i+jd)+q(i-jd)
      go to 10
  24  do 25 i = j2,iu
      p(i-j) = 2.*q(i)+q(i+jd)+q(i-jd)
  25  q(i) = q(i)-q(i+jh)-q(i-jh)
      go to 10
  20  do 23 i = j2,iu
      p(i-j) = 2.*q(i)+q(i+jd)+q(i-jd)
  23  q(i) = 0.
  10  do 22 l = lo,n,li
      a = c(l)
      as = a*s
      do 18 i = 2,m
      p(i) = as*p(i)
      d(i) = a*temp(i)
  18  g(i) = 2*a-d(i)
      g(2) = 0.
      d(m) = 0.
  19  ii = 2*id
      io = ii+1
      do 21 i = io,m,ii
      a = 1./(1.-d(i)*g(i+id)-g(i)*d(i-id))
      p(i) = a*(p(i)+d(i)*p(i+id)+g(i)*p(i-id))
      d(i) = d(i)*d(i+id)*a
  21  g(i) = g(i)*g(i-id)*a
      id = ii
      if (id-m/2.lt.0)go to 19
  16  id = ii/2
      io = id+1
      do 17 i = io,m,ii
  17  p(i) = p(i)+d(i)*p(i+id)+g(i)*p(i-id)
      ii = id
      if (id.gt.1)go to 16
  22  continue
      do 11 i = j2,iu
  11  q(i) = q(i)+p(i-j)
      select case (ko)
        case (1)
          go to 13
        case (2)
          go to 12
      end select
  12  lo = lo/2
      if (lo.eq.1)ko = 1
      go to 15
  13  lo = 2*lo
      if (lo.lt.n)go to 15
!vas  f90 modifi.
      deallocate(g)
      deallocate(p)
      deallocate(c)
      deallocate(d)
      deallocate(temp)

      return
      end
!
!   This routine is required if the CVS revision numbers are to
!   survive an optimization.
!
!
!   1998/02/03 23:54:21 meyer
!
      subroutine bunema_rev(i)
      CHARACTER*100 opt
      character*10 s
      if( i .eq. 0) s =  &
      '@(#)bunema.for,v 4.14\000'
      return
      end
