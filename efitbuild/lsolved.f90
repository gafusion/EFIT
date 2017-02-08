      subroutine decomp(ndim,n,a,cond,ipvt,work)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          decomp decomposes a matrix by gaussian elimination.     **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       ndim............declared row dimension of a                **
!**       n...............order of a                                 **
!**       a...............matrix to be triangularized                **
!**       cond............an estimate of a condition, for the linear **
!**                       system a*x=b changes in a and b may cause  **
!**                       cond times as large in x.  estimate only   **
!**                       if cond > 0                                **
!**       ipvt............pivot vector, ipvt(k) index of the k pivot **
!**                       row, ipvt(n)=(-1)**(number of inter-       **
!**                       changes), det(a)=ipvt(n)*a(1,1)*...*       **
!**                       a(n,n)                                     **
!**       work............work space array                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
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
      end                                                               
      subroutine solve(ndim, n, a, b, ipvt)                             
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          solve finds the solution of a linear set of             **
!**          algebric equations.  a call to decomp is required.      **
!**                                                                  **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       ndim............declared row dimension of a                **
!**       n...............order of a                                 **
!**       a...............triangularized matrix from decomp          **
!**       b...............right hand side vector, contain solution   **
!**                       on return                                  **
!**       ipvt............pivot vector from decomp                   **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
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
      end                                                               
!
!   This routine is required if the CVS revision numbers are to 
!   survive an optimization.
!
!
!   $Date: 2008/07/29 23:42:35 $ $Author: radhakri $
!
      subroutine lsolve_rev(i)
      CHARACTER*100 opt
      character*10 s 
      if( i .eq. 0) s =  &
      '@(#)$RCSfile: lsolved.f90,v $ $Revision: 1.1.2.1 $\000'
      return
      end
