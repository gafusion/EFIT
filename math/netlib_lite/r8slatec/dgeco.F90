      SUBROUTINE DGECO(A,LDA,N,IPVT,RCOND,Z)                            
      INTEGER LDA,N,IPVT(1)                                             
      DOUBLE PRECISION A(LDA,1),Z(1)                                    
      DOUBLE PRECISION RCOND                                            
!                                                                       
!     DGECO FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION   
!     AND ESTIMATES THE CONDITION OF THE MATRIX.                        
!                                                                       
!     IF  RCOND  IS NOT NEEDED, DGEFA IS SLIGHTLY FASTER.               
!     TO SOLVE  A*X = B , FOLLOW DGECO BY DGESL.                        
!     TO COMPUTE  INVERSE(A)*C , FOLLOW DGECO BY DGESL.                 
!     TO COMPUTE  DETERMINANT(A) , FOLLOW DGECO BY DGEDI.               
!     TO COMPUTE  INVERSE(A) , FOLLOW DGECO BY DGEDI.                   
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        A       DOUBLE PRECISION(LDA, N)                               
!                THE MATRIX TO BE FACTORED.                             
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  A .                
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS         
!                WHICH WERE USED TO OBTAIN IT.                          
!                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE       
!                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER          
!                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.       
!                                                                       
!        IPVT    INTEGER(N)                                             
!                AN INTEGER VECTOR OF PIVOT INDICES.                    
!                                                                       
!        RCOND   DOUBLE PRECISION                                       
!                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .        
!                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS       
!                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE             
!                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . 
!                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION     
!                           1.0 + RCOND .EQ. 1.0                        
!                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING           
!                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF         
!                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE          
!                UNDERFLOWS.                                            
!                                                                       
!        Z       DOUBLE PRECISION(N)                                    
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.  
!                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS      
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT           
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .                    
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     LINPACK DGEFA                                                     
!     BLAS DAXPY,DDOT,DSCAL,DASUM                                       
!     FORTRAN DABS,DMAX1,DSIGN                                          
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION DDOT,EK,T,WK,WKM                                 
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM                           
      INTEGER INFO,J,K,KB,KP1,L                                         
!                                                                       
!                                                                       
!     COMPUTE 1-NORM OF A                                               
!                                                                       
      ANORM = 0.0D0                                                     
      DO 10 J = 1, N                                                    
         ANORM = DMAX1(ANORM,DASUM(N,A(1,J),1))                         
   10 CONTINUE                                                          
!                                                                       
!     FACTOR                                                            
!                                                                       
      CALL DGEFA(A,LDA,N,IPVT,INFO)                                     
!                                                                       
!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .              
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .  
!     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE      
!     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE  
!     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID    
!     OVERFLOW.                                                         
!                                                                       
!     SOLVE TRANS(U)*W = E                                              
!                                                                       
      EK = 1.0D0                                                        
      DO 20 J = 1, N                                                    
         Z(J) = 0.0D0                                                   
   20 CONTINUE                                                          
      DO 100 K = 1, N                                                   
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))                      
         IF (DABS(EK-Z(K)) .LE. DABS(A(K,K))) GO TO 30                  
            S = DABS(A(K,K))/DABS(EK-Z(K))                              
            CALL DSCAL(N,S,Z,1)                                         
            EK = S*EK                                                   
   30    CONTINUE                                                       
         WK = EK - Z(K)                                                 
         WKM = -EK - Z(K)                                               
         S = DABS(WK)                                                   
         SM = DABS(WKM)                                                 
         IF (A(K,K) .EQ. 0.0D0) GO TO 40                                
            WK = WK/A(K,K)                                              
            WKM = WKM/A(K,K)                                            
         GO TO 50                                                       
   40    CONTINUE                                                       
            WK = 1.0D0                                                  
            WKM = 1.0D0                                                 
   50    CONTINUE                                                       
         KP1 = K + 1                                                    
         IF (KP1 .GT. N) GO TO 90                                       
            DO 60 J = KP1, N                                            
               SM = SM + DABS(Z(J)+WKM*A(K,J))                          
               Z(J) = Z(J) + WK*A(K,J)                                  
               S = S + DABS(Z(J))                                       
   60       CONTINUE                                                    
            IF (S .GE. SM) GO TO 80                                     
               T = WKM - WK                                             
               WK = WKM                                                 
               DO 70 J = KP1, N                                         
                  Z(J) = Z(J) + T*A(K,J)                                
   70          CONTINUE                                                 
   80       CONTINUE                                                    
   90    CONTINUE                                                       
         Z(K) = WK                                                      
  100 CONTINUE                                                          
      S = 1.0D0/DASUM(N,Z,1)                                            
      CALL DSCAL(N,S,Z,1)                                               
!                                                                       
!     SOLVE TRANS(L)*Y = W                                              
!                                                                       
      DO 120 KB = 1, N                                                  
         K = N + 1 - KB                                                 
         IF (K .LT. N) Z(K) = Z(K) + DDOT(N-K,A(K+1,K),1,Z(K+1),1)      
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 110                           
            S = 1.0D0/DABS(Z(K))                                        
            CALL DSCAL(N,S,Z,1)                                         
  110    CONTINUE                                                       
         L = IPVT(K)                                                    
         T = Z(L)                                                       
         Z(L) = Z(K)                                                    
         Z(K) = T                                                       
  120 CONTINUE                                                          
      S = 1.0D0/DASUM(N,Z,1)                                            
      CALL DSCAL(N,S,Z,1)                                               
!                                                                       
      YNORM = 1.0D0                                                     
!                                                                       
!     SOLVE L*V = Y                                                     
!                                                                       
      DO 140 K = 1, N                                                   
         L = IPVT(K)                                                    
         T = Z(L)                                                       
         Z(L) = Z(K)                                                    
         Z(K) = T                                                       
         IF (K .LT. N) CALL DAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)            
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 130                           
            S = 1.0D0/DABS(Z(K))                                        
            CALL DSCAL(N,S,Z,1)                                         
            YNORM = S*YNORM                                             
  130    CONTINUE                                                       
  140 CONTINUE                                                          
      S = 1.0D0/DASUM(N,Z,1)                                            
      CALL DSCAL(N,S,Z,1)                                               
      YNORM = S*YNORM                                                   
!                                                                       
!     SOLVE  U*Z = V                                                    
!                                                                       
      DO 160 KB = 1, N                                                  
         K = N + 1 - KB                                                 
         IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 150                    
            S = DABS(A(K,K))/DABS(Z(K))                                 
            CALL DSCAL(N,S,Z,1)                                         
            YNORM = S*YNORM                                             
  150    CONTINUE                                                       
         IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)                      
         IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0                            
         T = -Z(K)                                                      
         CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)                              
  160 CONTINUE                                                          
!     MAKE ZNORM = 1.0                                                  
      S = 1.0D0/DASUM(N,Z,1)                                            
      CALL DSCAL(N,S,Z,1)                                               
      YNORM = S*YNORM                                                   
!                                                                       
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM                         
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0                               
      RETURN                                                            
      END                                                               
