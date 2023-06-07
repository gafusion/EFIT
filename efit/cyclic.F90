!> The following set of routines are specific to the fast cyclic solver
!! originated by Holger St. John, optimised to realtime format by John Ferron
!! and then modified for use as an option here by Dylan Brennan.
!********************************************************************** 
!!
!> ef_init_cycred_data
!! Initialize a set of precomputed data for the cyclic reduction algorithm.
!! These data depend on the grid size and spacings.
!!
!**********************************************************************  
      function ef_init_cycred_data()
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      use var_buneman, only: delrgrid, delz
      implicit none
      integer*4 ef_init_cycred_data
      real*8 cosdii, sindi, denom, dr, dumy, dumy1
      real*8 dzsq, dz, dzdrsq, dzdr
      integer*4 i, j, k, l, nred, m1,jpow, jpowm1, jstep, jend, jstart, ind

      ef_init_cycred_data = 0

! ----------------------------------------------------------------------
!     Power of 2 that specifies the grid height.
      i=1
      nhpwr = -1
      do j=1,11
        i = i*2
        if (i.eq.(nh - 1)) then
          nhpwr = j
          exit
        endif
      enddo

      if (nhpwr.eq.-1) then
         nhpwr = 5
         ef_init_cycred_data=1
         call errctrl_msg('ef_init_cycred_data',&
                          'grid height must be a power of 2 + 1')
         return 
      endif  

! ----------------------------------------------------------------------
!     Elements of the tridiagonal matrix.

      dr = delrgrid
      dz = delz

      dzsq = dz*dz
      dzdrsq = (dz/dr) * (dz/dr)
      dumy = dzsq/(2.0*dr)
      dumy1 = 2.0 * (1.0 + dzdrsq)

!     All elements of the matrix diagonal have the same value
      diag = dumy1 
      do i=0,nw-3
         denom = dumy/rgrid(i+2)
         diagl(i+1) = -1.0 * dzdrsq - denom
         diagu(i+1) = -1.0 * dzdrsq + denom
      enddo

! ----------------------------------------------------------------------
!     Misc. constants used in computing the rhs vector.
!
      rhs_a_dumy = dzdrsq + dzsq/(2.0 * rgrid(2) * dr)
      rhs_b_dumy = dzdrsq - dzsq/(2.0 * rgrid(nw-1) * dr)

! ----------------------------------------------------------------------
!     Constants used during the forward reduction procedure.

      ind = 0
!     nred steps are required
      nred = nhpwr - 1    
      jpow = 1   
!     k is the reduction step index
      do k=1,nred   
!                2**(k-1)
        jpowm1 = jpow   
!                2**k
        jpow   = jpow * 2 
!                2**k + 1
        jstart = jpow + 1   
!                nh - 2**k
        jend   = nh - jpow 
!                2**k 
        jstep  = jpow    

        do j=jstart,jend,jstep
          m1 = -1.0
          do l=1,jpowm1
            m1 = -1.0 * m1
 
            ind = ind + 1
            if (ind.gt.icycred_loopmax) then
              ef_init_cycred_data=1
              call errctrl_msg('ef_init_cycred_data', &
                   'constant data index in forward reduction is > max')
              return
            endif

              cosdii = -2.0 * cos( (2.0 * l - 1.0) * pi/jpow)
              sindi = sin( (2.0 * l - 1.0) * pi/jpow)
              alphab(ind) = m1 * sindi/jpowm1
              diag1(ind) = diag + cosdii
          enddo 

          if (ind.gt.icycred_loopmax) exit
              
        enddo
      enddo

! ----------------------------------------------------------------------
!     Constants used during the back solution procedure.

!            nhpwr
      nred = nred + 1 
!            2**nhpwr
      jpow = 2**nred 

      do k=nred,1,-1
!               2**k
        jstep = jpow
!              2**(k-1)
        jpow = jpow/2   
        jstart = jpow + 1
        jend = nh - jpow
        do j=jstart,jend,jstep
          m1 = -1.0
          do l=1,jpow
            m1 = -1.0 * m1

            ind = ind + 1
            if (ind.gt.icycred_loopmax) then
              ef_init_cycred_data=1
              call errctrl_msg('ef_init_cycred_data', &
                   'constant data index in back solution is > max')
              return
            endif 

            cosdii = -2.0 * cos( (2.0 * l - 1.0) * pi/(2**k)) 
            sindi = sin( (2.0 * l - 1.0) * pi/(2**k))
            alphab(ind) = m1 * sindi/jpow
            diag1(ind) = diag + cosdii
          enddo 

          if(ind.gt.icycred_loopmax) exit
             
        enddo 
      enddo 

! ----------------------------------------------------------------------
!     Vectors of precalculated values used by the tridiag routine.
!     These vectors are various combinations of the diagonals of the
!     matrices that are solved.

!     At this point "index" holds the number of reduction loops that are used.
      do i=1,ind
        beti(i,1) = 1.0/diag1(i)
!       not actually used
        abeti(i,1) = 0.0 
!       not actually used
        wk1(i,1) = 0.0   

        do j=2,nw-2
          wk1(i,j) = diagu(j-1) * beti(i,j-1)
          beti(i,j) =  &
            1.0/( diag1(i) - diagl(j) * wk1(i,j))
          abeti(i,j) =  &
            diagl(j) * beti(i,j)
        enddo
      enddo  
      end function ef_init_cycred_data


!********************************************************************** 
!! 
!>    cyclic_reduction
!!    The core routine to compute the flux on the grid using Holger StJohn's
!!    single cyclic reduction algorithm.
!!    
!!    @param f :
!!
!**********************************************************************
      subroutine cyclic_reduction(f)

      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      real*8, intent(inout) :: f(nwnh)
      integer*4 nv,ind,nred,jpow,i,k,jpowm1,jstart,jend,jstep,j,jd,jdp,jdm,l

! ----------------------------------------------------------------------
!     Forward reduction.

      nv = nw - 2

      ind = 0 
!     nred steps are required.
      nred = nhpwr - 1
      jpow = 1   
      do i=0,nv
        wk2(i+1) = 0.0
      enddo

!     k is the reduction step index. 
      do k=1,nred   
!                2**(k-1)
        jpowm1 = jpow 
!                2**k
        jpow   = jpow * 2  
!                2**k + 1 
        jstart = jpow + 1
!                nh - 2**k 
        jend   = nh - jpow 
!                2**k     
        jstep  = jpow  

        do j=jstart,jend,jstep
       
!         Index of the first element of jth row 
          jd = (j-1) * nw + 1
!         Next row up
          jdp = jd + jpowm1*nw  
!         Next row down
          jdm = jd - jpowm1*nw 

          phi(1:nv)=f(jdm+1:jdm+nv)+f(jdp+1:jdp+nv)
          
          do l=1,jpowm1
            ind = ind + 1
            call ef_tridiag2(f(jd+1),nv,ind)
          enddo

!         use the declared var constant, instead of 0.5 directly because
!         of change in alignment for different compilations

          f(jd+1:jd+nv)=f(jd+1:jd+nv)*0.5_dp

        enddo
      enddo

! ----------------------------------------------------------------------
!     Back solution.

!            nhpwr
      nred = nred + 1
!            2**nhpwr 
      jpow = 2**nred 

      do k=nred,1,-1
!               2**k
        jstep = jpow    
!              2**(k-1) 
        jpow = jpow/2  
        jstart = jpow + 1
        jend = nh - jpow

        do j=jstart,jend,jstep
!         Index of the first element of jth row 
          jd = (j-1) * nw   
!         Next row up 
          jdp = jd + jpow*nw 
!         Next row down 
          jdm = jd - jpow*nw  
 
          phi(1:nv)=f(jdm+2:jdm+1+nv)+f(jdp+2:jdp+1+nv)
          do l=1,jpow
            ind = ind + 1

            if (l.eq.1) then
              call ef_tridiag1(f(jd+2),nv,ind)
            else
              call ef_tridiag2(f(jd+2),nv,ind)
            endif
          enddo
        enddo
      enddo
      return
      end subroutine cyclic_reduction


!**********************************************************************
!!
!>  pflux_cycred use e the single cyclic reduction method of Holger StJohn
!! and John Ferron 
!! to get the flux on the grid resulting from the plasma current.
!!    
!!    @param psigrid :
!!    @param kerror :
!!
!**********************************************************************
      subroutine pflux_cycred(psigrid,kerror)

      include 'eparm.inc'
      include 'modules1.inc'
      implicit none
      real*8, intent(inout) :: psigrid(nwnh)
      integer*4, intent(out) :: kerror
      integer*4 i,ii,j,is,nn
      real*8 sia(nwnh)

      kerror = 0

! ----------------------------------------------------------------------
!     transpose the matrix, solver needs transposition (same as buneto)

      do i = 1,nw
        ii = (i-1)*nh 
        do j = 1,nh
          sia(i+(j-1)*nw) = psigrid(ii+j)
        enddo
      enddo

! ----------------------------------------------------------------------
!     Add the finite difference expressions for the second and nw-1 grid
!     columns to the rhs. These terms are known boundary terms and
!     hence do not form part of the tridiagonal matrix.

      is = nw
      nn = nh-2
      call vsma(sia, nw+1, nw+2, is, rhs_a_dumy, nn)
      call vsma(sia, 2*nw, 2*nw-1, is, rhs_b_dumy, nn)
! ----------------------------------------------------------------------
!     Do the cyclic reduction to get psigrid.

      call cyclic_reduction(sia)

! ----------------------------------------------------------------------
!     transpose the matrix back

      do i = 1,nh
        ii = (i-1)*nw
        do j = 1,nw
          psigrid(i+(j-1)*nh) = sia(ii+j)
        enddo
      enddo
      return
      end subroutine pflux_cycred


!**********************************************************************
!>
!!    Add the finite difference expressions for the boundary terms
!!    
!!    @param mat : matrix
!!    @param iu : starting index
!!    @param iv : starting index
!!    @param is : step size
!!    @param b : right hand side constant
!!    @param n : number of steps
!!
!**********************************************************************
      subroutine vsma(mat, iu, iv, is, b, n)

      use eparm, only: nwnh
      implicit none
      integer*4, intent(in) :: iu,iv,is,n
      real*8, intent(inout) :: mat(nwnh)
      real*8, intent(in) :: b
      integer*4 ii,i

      if(n.le.0) return

      ii = 0
      do i=1,n 
        mat(iv+ii) = mat(iu+ii)*b + mat(iv+ii)
        ii = ii + is
      enddo
      return
      end subroutine vsma


!**********************************************************************
!!
!>  ef_tridiag2 solves system of equations defined by a tridiagonal matrix.
!!  This is a version of the routine in Numerical Recipies which
!!  uses some values precalculated from the diagonals of the matrix.
!! 
!!  This routine also does a bunch of preparation to compute the rhs vector
!!  and accumulate the result.
!! 
!!    @param f : the vector in which the result is being accumulated.
!!    @param n : number of elements in the vector. 
!!    @param ind : index
!!
!!    Other useful variables:
!!     con : array of structures giving precalculated values
!!     wk1 : vector of precalculated values. 
!!     wk2,phi : vectors that when combined yield the right hand side vector
!!     alphab : a constant that would be used in computing the right hand side
!!             vector.
!!     v : temporary storage area for the result vector
!!
!**********************************************************************
      subroutine ef_tridiag2(f,n,ind)

      include 'eparm.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: n,ind
      real*8, intent(inout) :: f(nwnh)
      integer*4 j
      real*8 usave

      v(1)=(wk2(1)+alphab(ind)*phi(1))* &
           beti(ind,1)
      do j=2,n
        v(j)=wk2(j)*beti(ind,j)+  &
          alphab(ind)*phi(j)*beti(ind,j)-  &
          v(j-1) * abeti(ind,j)
      enddo

      usave = v(n)
      f(n) = usave + f(n)

      do j=n-1,1,-1
        usave = v(j)-wk1(ind,j+1)*usave
        f(j) = usave + f(j)
      enddo
      return
      end subroutine ef_tridiag2


!*********************************************************************
!!
!>  ef_tridiag1 solves system of equations defined by a tridiagonal matrix.
!!  This is a version of the routine in Numerical Recipies which
!!  uses some values precalculated from the diagonals of the matrix.
!! 
!!  This routine also does a bunch of preparation to compute the rhs vector
!!  and accumulate the result.
!! 
!!  Basically the same as tridiag2 except:
!!  f is the input instead of wk2.  At the end of the routine, f as it
!!  was at the start of the routine has been copied into wk2.
!! 
!!  The result is written directly into f rather than adding the result
!!  to the input value of f.
!! 
!!  @param con : array of structures giving precalculated values.
!!  @param wk1 : vector of precalculated values. 
!!  @param wk2,phi : vectors that when combined yield the right hand side vector
!!  @param alphab : a constant that would be used in computing the right hand side
!!         vector.
!!  @param f : the vector in which the result is being accumulated.
!!  @param v : temporary storage area for the result vector
!!  @param n : number of elements in the vector. 
!!
!!********************************************************************
      subroutine ef_tridiag1(f,n,ind)

      include 'eparm.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: n,ind
      real*8, intent(inout) :: f(nwnh)
      integer*4 j
      real*8 usave
 
      v(1)=(f(1)+alphab(ind)*phi(1))* &
                    beti(ind,1)
 
      do j=2,n
       v(j) = f(j)*beti(ind,j)+  &
         alphab(ind)*phi(j)*beti(ind,j)-  &
         v(j-1)*abeti(ind,j)
      enddo
 
      wk2(n) = f(n)
      f(n) = v(n)
 
      do j=n-1,1,-1
        wk2(j) = f(j)
        v(j)=v(j)-wk1(ind,j+1)*v(j+1)
        f(j) = v(j)
      enddo
      return
      end subroutine ef_tridiag1
