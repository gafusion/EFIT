!********************************************************************** 
!!
!>    buneto sets up the appropriate arrays for the
!!     Buneman's solver.
!!
!!    @param psi :
!!    @param nwb :
!!    @param nhb :
!!    @param sia :
!!********************************************************************* 
      subroutine buneto(psi,nwb,nhb,sia) 
      use set_kinds 
      use var_buneman
      use eparm, only: nw, nh
 
      implicit integer*4 (i-n), real*8 (a-h, o-z) 
      real(dp), dimension(nwb*nhb), intent(inout) :: psi, sia
      integer(i4), intent(in) :: nwb, nhb
      mno = max(nw,nh)+1 
      m = nw-1
      n = nh-1
      s = drdz2 
      shift = rgrid1 
      dr = delrgrid 
      dz = delz 
!---------------------------------------------------------------------- 
!     copy psi into sia rowwise                                      -- 
!---------------------------------------------------------------------- 
      do i = 1,nwb 
        ii = (i-1)*nhb + 1 
        do j = 1, nhb 
          sia(i+j*nwb-nwb) = psi(ii-1+j) 
        enddo
      enddo
      ia = nwb+nwb 
      ju = n*nwb 
!----------------------------- 
!     set up for rzpois 
!----------------------------- 
      do i = ia,ju,nwb 
        sia(i-m+1) = sia(i-m+1)+(.5_dp+.25_dp/(1.0+shift/dr))*sia(i-m)/s 
        sia(i-1) = sia(i-1)+(.5_dp-.25_dp/(m-1+shift/dr))*sia(i)/s 
      enddo
      nwhbb = nwb*nhb 
      call rzpois(nwhbb,sia) 
!------------------------------------------------------------------------------ 
!     copy sia back into psi columnwise                                      -- 
!------------------------------------------------------------------------------ 
      do i = 2,n 
        ii = (i-1)*nwb + 1 
        do j = 2,m 
          psi(i+j*nhb-nhb) = sia(ii-1+j) 
        enddo
      enddo
      return 
      end subroutine buneto
 
!********************************************************************** 
!!
!>    rzpois solves for the poloidal flux using the
!!    Buneman's method.
!! 
!!    @param nq :
!!    @param q :
!!
!!********************************************************************* 
      subroutine rzpois(nq,q) 
      use set_kinds 
      use var_buneman 
 
      implicit none 
      integer*4, intent(in) :: nq
      real(dp), dimension(nq), intent(inout) :: q
      real(dp) flag 
      real(dp) shftdr,a,pitmp,as 
      integer*4 i,l,lo,ju,n222,id,jd,ko,k4,li,jh,jt,ji,jo
      integer*4 j2,iu,j,ii,io,iallocate_stat 
 
      real(dp),dimension(:), allocatable :: g 
      real(dp),dimension(:), allocatable :: p 
      real(dp),dimension(:), allocatable :: c 
      real(dp),dimension(:), allocatable :: d 
      real(dp),dimension(:), allocatable :: temp 
 
      flag = 0.0 
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
 
      g = 0. 
      p = 0. 
      d = 0. 
 
      if (flag .eq. 0.) then 
        temp = 0. 
        c = 0. 
 
        shftdr = shift/dr 
        do i = 2,m 
          temp(i) = 1. - .5_dp/(i+shftdr-1.) 
        enddo
        ju = (n-1)*(m+1) 
        n222 = n/2 
        c(n222) = 0. 
        lo = n/2 
   1    l = lo/2 
        c(l) = sqrt(2.+c(lo)) 
        lo = l 
   2    c(n-l) = -c(l) 
        l = l+2*lo 
        if ((2*l/n)*(2*lo-3)< 0) then
          goto 4
        elseif ((2*l/n)*(2*lo-3) == 0) then 
          goto 3
        else 
          goto 1
        end if
   3    c(l) = (c(l+lo)+c(l-lo))/c(lo) 
        go to 2 
   4    do l = 2,n 
          c(l-1) = 1./(2.+s*(2.-c(l-1))) 
        enddo
        flag = 1.
      endif
      lo = n/2 
      ko = 2 
      id = 1 
  15  li = 2*lo 
      k4 = 2*ko-li/n 
      jd = (m+1)*n/li 
      jh = (m+1)*(n/(2*li)) 
      jt = jd+jh 
      ji = 2*jd 
      jo = jd*ko 
      do j = jo,ju,ji 
        j2 = j+2 
        iu = j+m 
        select case (k4)
        case (1)
          do i = j2,iu 
            p(i-j) = 2.*q(i)+q(i+jd)+q(i-jd) 
            q(i) = 0.
          enddo 
        case (2)
          do i = j2,iu 
            p(i-j) = 2.*q(i)+q(i+jd)+q(i-jd) 
            q(i) = q(i)-q(i+jh)-q(i-jh) 
          enddo
        case (3)
          do i = j2,iu 
            p(i-j) = 2.*q(i) 
            q(i) = q(i+jd)+q(i-jd)
          enddo 
        case default
          do i = j2,iu 
            pitmp = q(i)-q(i+jt)-q(i-jt) 
            q(i) = q(i)-q(i+jh)-q(i-jh)+q(i+jd)+q(i-jd) 
            p(i-j) = pitmp+q(i)
          enddo 
        end select
        do l = lo,n,li 
          a = c(l) 
          as = a*s 
          do i = 2,m 
            p(i) = as*p(i) 
            d(i) = a*temp(i) 
            g(i) = 2*a-d(i)
          enddo
          g(2) = 0. 
          d(m) = 0. 
  19      ii = 2*id 
          io = ii+1 
          do i = io,m,ii 
            a = 1./(1.-d(i)*g(i+id)-g(i)*d(i-id)) 
            p(i) = a*(p(i)+d(i)*p(i+id)+g(i)*p(i-id)) 
            d(i) = d(i)*d(i+id)*a 
            g(i) = g(i)*g(i-id)*a
          enddo
          id = ii 
          if (id-m/2.lt.0) go to 19 
  16      id = ii/2 
          io = id+1 
          do i = io,m,ii 
            p(i) = p(i)+d(i)*p(i+id)+g(i)*p(i-id)
          enddo
          ii = id 
          if (id.gt.1) go to 16
        enddo
        do i = j2,iu 
          q(i) = q(i)+p(i-j)
        enddo
      enddo
      select case (ko)
      case (1)
        lo = 2*lo 
        if (lo.lt.n) go to 15 
      case default
        lo = lo/2 
        if (lo.eq.1) ko = 1 
        go to 15 
      end select
      deallocate(g) 
      deallocate(p) 
      deallocate(c) 
      deallocate(d) 
      deallocate(temp) 
 
      return 
      end subroutine rzpois
