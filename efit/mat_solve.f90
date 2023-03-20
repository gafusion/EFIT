!**********************************************************************
!>
!!    decomp decomposes a matrix by gaussian elimination.
!!
!!    @param ndim : declared row dimension of a 
!!
!!    @param n : order of a     
!!
!!    @param a : matrix to be triangularized
!!
!!    @param cond : an estimate of a condition, for the linear
!!                       system a*x=b changes in a and b may cause 
!!                       cond times as large in x.  estimate only
!!                       if cond > 0   
!!
!!    @param ipvt : pivot vector, ipvt(k) index of the k pivot 
!!                       row, ipvt(n)=(-1)**(number of inter- 
!!                       changes), det(a)=ipvt(n)*a(1,1)*...*  
!!                       a(n,n)  
!!
!!    @param work : work space array  
!!
!**********************************************************************
      subroutine decomp(ndim,n,a,cond,ipvt,work)
      integer*4 ndim,n
      real*8 a(ndim,n),cond,work(n)
      integer*4 ipvt(n)
      real*8 ek, t, anorm, ynorm, znorm
      integer*4 nm1, i, j, k, kp1, kb, km1, m
      ipvt(n) = 1
      if (n .eq. 1) go to 80
      nm1 = n - 1
!---------------------------------------------------------------------- 
!--   compute 1-norm of a                                            -- 
!---------------------------------------------------------------------- 
      anorm = 0.0                                                       
      do 10 j = 1, n                                                    
         t = 0.0                                                        
         do 5 i = 1, n                                                  
            t = t + abs(a(i,j))                                         
    5    continue                                                       
         if (t .gt. anorm) anorm = t                                    
   10 continue                                                          
!---------------------------------------------------------------------- 
!--   gaussian elimination with partial pivoting                     -- 
!---------------------------------------------------------------------- 
      do 35 k = 1,nm1                                                   
         kp1= k+1                                                       
!---------------------------------------------------------------------- 
!--      find pivot                                                  -- 
!---------------------------------------------------------------------- 
         m = k                                                          
         do 15 i = kp1,n                                                
            if (abs(a(i,k)) .gt. abs(a(m,k))) m = i                     
   15    continue                                                       
         ipvt(k) = m                                                    
         if (m .ne. k) ipvt(n) = -ipvt(n)                               
         t = a(m,k)                                                     
         a(m,k) = a(k,k)                                                
         a(k,k) = t                                                     
!---------------------------------------------------------------------- 
!--      skip step if pivot is zero                                  -- 
!---------------------------------------------------------------------- 
         if (t .eq. 0.0) go to 35                                       
!---------------------------------------------------------------------- 
!--      compute multipliers                                         -- 
!---------------------------------------------------------------------- 
         do 20 i = kp1,n                                                
             a(i,k) = -a(i,k)/t                                         
   20    continue                                                       
!---------------------------------------------------------------------- 
!--      interchange and eliminate by columns                        -- 
!---------------------------------------------------------------------- 
         do 30 j = kp1,n                                                
             t = a(m,j)                                                 
             a(m,j) = a(k,j)                                            
             a(k,j) = t                                                 
             if (t .eq. 0.0) go to 30                                   
             do 25 i = kp1,n                                            
                a(i,j) = a(i,j) + a(i,k)*t                              
   25        continue                                                   
   30    continue                                                       
   35 continue                                                          
!---------------------------------------------------------------------- 
!--   cond = (1-norm of a)*(an estimate of 1-norm of a-inverse)      -- 
!--   estimate obtained by one step of inverse iteration for the     -- 
!--   small singular vector.  this involves solving two systems      -- 
!--   of equations, (a-transpose)*y = e  and  a*z = y  where  e      -- 
!--   is a vector of +1 or -1 chosen to cause growth in y.           -- 
!--   estimate = (1-norm of z)/(1-norm of y)                         -- 
!--                                                                  -- 
!--   solve (a-transpose)*y = e                                      -- 
!---------------------------------------------------------------------- 
      if (cond.le.0.0) return
      do 50 k = 1, n                                                    
         t = 0.0                                                        
         if (k .eq. 1) go to 45                                         
         km1 = k-1                                                      
         do 40 i = 1, km1                                               
            t = t + a(i,k)*work(i)                                      
   40    continue                                                       
   45    ek = 1.0                                                       
         if (t .lt. 0.0) ek = -1.0                                      
         if (a(k,k) .eq. 0.0) go to 90                                  
         work(k) = -(ek + t)/a(k,k)                                     
   50 continue                                                          
      do 60 kb = 1, nm1                                                 
         k = n - kb                                                     
         t = 0.0                                                        
         kp1 = k+1                                                      
         do 55 i = kp1, n                                               
            t = t + a(i,k)*work(k)                                      
   55    continue                                                       
         work(k) = t                                                    
         m = ipvt(k)                                                    
         if (m .eq. k) go to 60                                         
         t = work(m)                                                    
         work(m) = work(k)                                              
         work(k) = t                                                    
   60 continue                                                          
!                                                                       
      ynorm = 0.0                                                       
      do 65 i = 1, n                                                    
         ynorm = ynorm + abs(work(i))                                   
   65 continue                                                          
!---------------------------------------------------------------------- 
!--   solve a*z = y                                                  -- 
!---------------------------------------------------------------------- 
      call solve(ndim, n, a, work, ipvt)                                
!                                                                       
      znorm = 0.0                                                       
      do 70 i = 1, n                                                    
         znorm = znorm + abs(work(i))                                   
   70 continue                                                          
!---------------------------------------------------------------------- 
!--   estimate condition                                             -- 
!---------------------------------------------------------------------- 
      cond = anorm*znorm/ynorm                                          
      if (cond .lt. 1.0) cond = 1.0                                     
      return                                                            
!---------------------------------------------------------------------- 
!--   1-by-1                                                         -- 
!---------------------------------------------------------------------- 
   80 cond = 1.0                                                        
      if (a(1,1) .ne. 0.0) return                                       
!---------------------------------------------------------------------- 
!--   exact singularity                                              -- 
!---------------------------------------------------------------------- 
   90 cond = 1.0e+32   
      return                                                            
      end subroutine decomp 

!**********************************************************************
!>
!!    solve finds the solution of a linear set of
!!    algebric equations.  a call to decomp is required.
!!
!!    @param ndim : declared row dimension of a 
!!
!!    @param n : rder of a  
!!
!!    @param a : triangularized matrix from decomp
!!
!!    @param b : right hand side vector, contain solution
!!                       on return
!!
!!    @param ipvt : pivot vector from decomp
!!
!**********************************************************************   
      subroutine solve(ndim, n, a, b, ipvt)
      integer*4 ndim, n, ipvt(n)
      real*8 a(ndim,n),b(n)
      integer*4 kb, km1, nm1, kp1, i, k, m
      real*8 t
!---------------------------------------------------------------------- 
!--   forward elimination                                            -- 
!---------------------------------------------------------------------- 
      if (n .eq. 1) go to 50                                            
      nm1 = n-1                                                         
      do 20 k = 1, nm1                                                  
         kp1 = k+1                                                      
         m = ipvt(k)                                                    
         t = b(m)                                                       
         b(m) = b(k)                                                    
         b(k) = t                                                       
         do 10 i = kp1, n                                               
             b(i) = b(i) + a(i,k)*t                                     
   10    continue                                                       
   20 continue                                                          
!---------------------------------------------------------------------- 
!--   back substitution                                              -- 
!---------------------------------------------------------------------- 
      do 40 kb = 1,nm1                                                  
         km1 = n-kb                                                     
         k = km1+1                                                      
         b(k) = b(k)/a(k,k)                                             
         t = -b(k)                                                      
         do 30 i = 1, km1                                               
             b(i) = b(i) + a(i,k)*t                                     
   30    continue                                                       
   40 continue                                                          
   50 b(1) = b(1)/a(1,1)  
      return                                                            
      end subroutine solve 

!**********************************************************************
!>
!!    The singular value decomposition A = UQVT of an m-by-n matrix A
!!    can be computed using the lapack subroutine DGESVD:
!!    
!!    CALL DGESVD(JOBU,JOBVT,M,N,A,IA,S,U,MM,V,NN,WK3,MM3,INFO)
!!    
!!    where
!!    JOBU - Specifies options for computing all or part of the matrix U:
!!         = 'A':  all M columns of U are returned in array U:
!!         = 'S':  the first min(m,n) columns of U (the left singular
!!                 vectors) are returned in the array U;
!!         = 'O':  the first min(m,n) columns of U (the left singular
!!                 vectors) are overwritten on the array A;
!!         = 'N':  no columns of U (no left singular vectors) are
!!                 computed.
!!
!!    JOBVT - Specifies options for computing all or part of the matrix V**T:
!!          = 'A':  all N rows of V**T are returned in the array VT;
!!          = 'S':  the first min(m,n) rows of V**T (the right singular
!!                  vectors) are returned in the array VT;
!!          = 'O':  the first min(m,n) rows of V**T (the right singular
!!                  vectors) are overwritten on the array A;
!!          = 'N':  no rows of V**T (no right singular vectors) are
!!                  computed.
!!    
!!    M    - Row dimension of the matrix A.
!!    
!!    N    - Column dimension of matrix A.
!!    
!!    A    - REAL8 M-by-N matrix A  [A(IA,N)].
!!    
!!    IA   - Row dimension of the array A.
!!    
!!    S    - Returns the singular values in descending order [S(N)].
!!    
!!    U    - Returns the orthogonal matrix U if requested [U(MM,M)].
!!    
!!    MM   - Row dimension of the array U.
!!    
!!    V    - Returns the orthogonal matrix V if requested [V(NN,N)].
!!    
!!    NN   - Row dimension of the array V.
!!    
!!    WK3 - Work vector [WORK(3*MM)].
!!    
!!    INFO = 0 implies SVD is found
!!    INFO > 0 implies SVD is did not converge
!!    INFO < 0 implies illegal input
!!    
!!    Remarks:
!!
!!    1. JOBVT and JOBU cannot both be 'O'.
!!    
!!    2.  Anything less than the full SVD allows for substantial overwriting
!!    by identifying the arrays U and/or V with the array A in the calling
!!    sequence.
!!
!!    3.  Any option that returns the matrix V in array A requires M >= N.
!!    
!!    4.  A leading application of the SVD is the linear least squares
!!    problem, especially when the data matrix is numerically rank
!!    deficient.
!!    
!**********************************************************************
      subroutine sdecm(a, ia, m, n, b, ib, nb, s, wk, ier)
      use set_kinds
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      parameter (nn = 100, mm = 3000)
      parameter (mm3 = 3*mm)
      dimension a(ia,n), s(n), wk(n), b(ib,nb) 
      dimension bb(mm), u(mm,mm), v(nn,nn)
      dimension wk3(mm3)
      character*1 jobu, jobvt
!
      if ((nn .lt. n) .or. (mm .lt. m)) then 
         write (6,100) n,m
         ier = 129
         return
      endif
      n11=11
      jobu='A'
      jobvt='A'
      info=0
!---------------------------------------------------------------------
!--   Changed to LAPACK routine                                     --
!---------------------------------------------------------------------
      call dgesvd(jobu,jobvt,int(m,8),int(n,8),a,int(ia,8),s,u, &
                  int(mm,8),v,int(nn,8),wk3,int(mm3,8),int(info,8))
!
      if (info.gt.0) ier = 129                   ! Error message
      if (s(1).le.0) ier = 129                   ! This better not happen
      if (s(n)/s(1) .le. 1.e-7_dp) ier = 33      ! Ill-conditioned or rank-deficient matrix
!
      do i=1,n              ! Copy V into A
        a(i,1:n) = v(1:n,i)
      enddo
      if (nb.eq.1) then     ! B is a vector
        do i=1,m            ! Compute U**T * B
          bb(I) = sum(u(1:m,i)*B(1:m,1))
        enddo
        b(1:m,1) = bb(1:m)  ! Replace B with U**T * B 
      else                  ! B is a matrix
        do i=1,nb           ! Replace B with U**T
          b(1:m,i) = u(i,1:m)
        enddo
      endif
!
 100  format(' Array dimensions', 2i6 ,' too large. Recompile.')
      return
      end subroutine sdecm
