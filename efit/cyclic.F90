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
      real*8 u0, cosdii, sindi, denom, dr, dumy, dumy1
      real*8 dzsq, dz, dzdrsq, dzdr
      integer*4 i,j,k,l, nred, m1,jpow, jpowm1, jstep, jend, jstart, index

      ef_init_cycred_data = 0

      u0 = 4.0 * pi * 1.0e-7_dp

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

      index = 0
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
 
            index = index + 1
            if (index.gt.icycred_loopmax) then
              ef_init_cycred_data=1
              call errctrl_msg('ef_init_cycred_data', &
                   'constant data index in forward reduction is > max')
              return
            endif

              cosdii = -2.0 * cos( (2.0 * l - 1.0) * pi/jpow)
              sindi = sin( (2.0 * l - 1.0) * pi/jpow)
              alphab(index) = m1 * sindi/jpowm1
              diag1(index) = diag + cosdii
          enddo 

          if (index.gt.icycred_loopmax) exit
              
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

            index = index + 1
            if (index.gt.icycred_loopmax) then
              ef_init_cycred_data=1
              call errctrl_msg('ef_init_cycred_data', &
                   'constant data index in back solution is > max')
              return
            endif 

            cosdii = -2.0 * cos( (2.0 * l - 1.0) * pi/(2**k)) 
            sindi = sin( (2.0 * l - 1.0) * pi/(2**k))
            alphab(index) = m1 * sindi/jpow
            diag1(index) = diag + cosdii
          enddo 

          if(index.gt.icycred_loopmax) exit
             
        enddo 
      enddo 

! ----------------------------------------------------------------------
!     Vectors of precalculated values used by the tridiag routine.
!     These vectors are various combinations of the diagonals of the
!     matrices that are solved.

!     At this point "index" holds the number of reduction loops that are used.
      do i=1,index
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
!**********************************************************************
      subroutine cyclic_reduction(f)

      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      dimension f(*)

      constant=0.5_dp
! ----------------------------------------------------------------------
!     Forward reduction.

      nv = nw - 2

      index = 0 
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

          call ef_vadd_shrt(f(jdm+1),f(jdp+1),phi,nv)
          
          do l=1,jpowm1
            index = index + 1
            call ef_tridiag2(f(jd+1),nv,index)
          enddo

!         use the declared var constant, instead of 0.5 directly because
!         of change in alignment for different compilations

          call ef_vmul_const_shrt(f(jd+1),constant,f(jd+1),nv)

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
 
          call ef_vadd_shrt(f(jdm+2),f(jdp+2),phi,nv)
          do l=1,jpow
            index = index + 1

            if (l.eq.1) then
              call ef_tridiag1(f(jd+2),nv,index)
            else
              call ef_tridiag2(f(jd+2),nv,index)
            endif
          enddo
        enddo
      enddo

! All done.

      return
      end subroutine cyclic_reduction


!**********************************************************************
!!
!>  pflux_cycred use e the single cyclic reduction method of Holger StJohn
!! and John Ferron 
!! to get the flux on the grid resulting from the plasma current.
!!
!**********************************************************************
      subroutine pflux_cycred(psigrid,sia,kerror)

      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      dimension psigrid(*),sia(*)
      integer*4 initresult, ef_init_cycred_data

      kerror = 0

! ----------------------------------------------------------------------
!     Call the initialization and check the result
      initresult=ef_init_cycred_data()
      if (initresult.eq.1) then
          kerror = 1
          return
      endif

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

      ius = nw
      ivs = nw
      iys = nw
      nn = nh-2
      call vsma_(sia(nw+1), ius, &
          rhs_a_dumy, &
          sia(nw+2), ivs, &
          sia(nw+2), iys, &
          nn)
      call vsma_(sia(2*nw), ius, &
          rhs_b_dumy, &
          sia(2*nw-1), ivs, &
          sia(2*nw-1), iys, &
          nn)
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

! All done.

      return
      end subroutine pflux_cycred

! ======================================================================
! ======================================================================
!  +++++++++++++++++++++++
!  SUBROUTINE: vsma_
!  +++++++++++++++++++++++
      subroutine vsma_(a, ia, b, c, ic, d, id, n)
      implicit integer*4 (i-n), real*8 (a-h, o-z)

      dimension a(1),c(1),d(1) ! used as pointer inside array

      if (n.le.0) then
        return
      endif

      iai = 1
      ici = 1
      idi = 1

      do i=1,n 
        d(idi) = a(iai)*b + c(ici)
        iai = iai + ia
        ici = ici + ic
        idi = idi + id
      enddo
      return
      end subroutine vsma_


!**********************************************************************
!>
!!    this subroutine multiplies two vectors together
!!    
!!
!!    @param vin1 :
!!
!!    @param vin2 :
!!
!!    @param out :
!!
!!    @param nelements :
!!
!**********************************************************************
      subroutine ef_vvmul(vin1,vin2,out,nelements)
      implicit integer*4 (i-n), real*8 (a-h, o-z)

      dimension vin1(1),vin2(1),out(1) ! Could be used as pointer inside array

      do i=1,nelements
        out(i) = vin1(i) * vin2(i)
      enddo
      return
      end subroutine ef_vvmul

!**********************************************************************
!!
!>  ef_tridiag2 solves system of equations defined by a tridiagonal matrix.
!!  This is a version of the routine in Numerical Recipies which
!!  uses some values precalculated from the diagonals of the matrix.
!! 
!!  This routine also does a bunch of preparation to compute the rhs vector
!!  and accumulate the result.
!! 
!!    @param con : array of structures giving precalculated values
!!    @param wk1 : vector of precalculated values. 
!!    @param wk2,phi : vectors that when combined yield the right hand side vector
!!    @param alphab : a constant that would be used in computing the right hand side
!!         vector.
!!    @param f : the vector in which the result is being accumulated.
!!    @param v : temporary storage area for the result vector
!!    @param  n : number of elements in the vector. 
!!
!**********************************************************************
      subroutine ef_tridiag2(f,n,index)

      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      dimension f(n)

      v(1)=(wk2(1)+alphab(index)*phi(1))* &
           beti(index,1)
      do j=2,n
        v(j)=wk2(j)*beti(index,j)+  &
          alphab(index)*phi(j)*beti(index,j)-  &
          v(j-1) * abeti(index,j)
      enddo

      usave = v(n)
      f(n) = usave + f(n)

      do j=n-1,1,-1
        usave = v(j)-wk1(index,j+1)*usave
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
      subroutine ef_tridiag1(f,n,index)

      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      dimension f(n)
 
      v(1)=(f(1)+alphab(index)*phi(1))* &
                    beti(index,1)
 
      do j=2,n
       v(j) = f(j)*beti(index,j)+  &
         alphab(index)*phi(j)*beti(index,j)-  &
         v(j-1)*abeti(index,j)
      enddo
 
      wk2(n) = f(n)
      f(n) = v(n)
 
      do j=n-1,1,-1
        wk2(j) = f(j)
        v(j)=v(j)-wk1(index,j+1)*v(j+1)
        f(j) = v(j)
      enddo
      return
      end subroutine ef_tridiag1


!*********************************************************************
!!
!>  ef_vadd_shrt adds two vectors.  Handles vectors one 
!!  element at a time.  Good for
!!  short, unaligned vectors but not optimized at all for long vectors.
!! 
!!  @param vector_out : vector1 + vector2
!! 
!!  @param mvector1 : input vector
!!  @param vector2 : input vector
!!  @param vector_out : output vector
!!  @param nelements : number of elements in the vectors.
!!
!***********************************************************************
      subroutine ef_vadd_shrt(vector1,vector2,vector_out, &
                              nelements)
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension vector1(1),vector2(1),vector_out(1)

      do i=1,nelements
        vector_out(i) = vector1(i) + vector2(i)
      enddo
      return
      end subroutine ef_vadd_shrt

      
!**************************************************************************
!>    ef_vmul_const_shrt Multiply a vector by a constant.
!!    Handles vectors one element at a time.  Good for short, unaligned vectors
!!    but not optimized at all for long vectors.\n
!! 
!!    vector_out = vector1 * constant
!! 
!!    @param vector1 : input vector
!!    @param constant : constant value
!!    @param vector_out : output vector
!!    @param nelements : number of elements in the vectors. Must be at least 2.
!**************************************************************************
      subroutine ef_vmul_const_shrt(vector1,constant,vector_out,nelements)
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension vector1(1),vector_out(1)

      do i=1,nelements
        vector_out(i) = vector1(i) * constant
      enddo
      return
      end subroutine ef_vmul_const_shrt
