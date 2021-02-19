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
      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)                               
      INTEGER LDA,N,IPVT(1),INFO                                        
      DOUBLE PRECISION A(LDA,1)                                         
!                                                                       
!     DGEFA FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION.  
!                                                                       
!     DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED            
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.          
!     (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA) .                   
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
!        INFO    INTEGER                                                
!                = 0  NORMAL VALUE.                                     
!                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR       
!                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES        
!                     INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO  
!                     IF CALLED.  USE  RCOND  IN DGECO FOR A RELIABLE   
!                     INDICATION OF SINGULARITY.                        
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DAXPY,DSCAL,IDAMAX                                           
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION T                                                
      INTEGER IDAMAX,J,K,KP1,L,NM1                                      
!                                                                       
!                                                                       
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING                        
!                                                                       
      INFO = 0                                                          
      NM1 = N - 1                                                       
      IF (NM1 .LT. 1) GO TO 70                                          
      DO 60 K = 1, NM1                                                  
         KP1 = K + 1                                                    
!                                                                       
!        FIND L = PIVOT INDEX                                           
!                                                                       
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1                             
         IPVT(K) = L                                                    
!                                                                       
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED          
!                                                                       
         IF (A(L,K) .EQ. 0.0D0) GO TO 40                                
!                                                                       
!           INTERCHANGE IF NECESSARY                                    
!                                                                       
            IF (L .EQ. K) GO TO 10                                      
               T = A(L,K)                                               
               A(L,K) = A(K,K)                                          
               A(K,K) = T                                               
   10       CONTINUE                                                    
!                                                                       
!           COMPUTE MULTIPLIERS                                         
!                                                                       
            T = -1.0D0/A(K,K)                                           
            CALL DSCAL(N-K,T,A(K+1,K),1)                                
!                                                                       
!           ROW ELIMINATION WITH COLUMN INDEXING                        
!                                                                       
            DO 30 J = KP1, N                                            
               T = A(L,J)                                               
               IF (L .EQ. K) GO TO 20                                   
                  A(L,J) = A(K,J)                                       
                  A(K,J) = T                                            
   20          CONTINUE                                                 
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)                  
   30       CONTINUE                                                    
         GO TO 50                                                       
   40    CONTINUE                                                       
            INFO = K                                                    
   50    CONTINUE                                                       
   60 CONTINUE                                                          
   70 CONTINUE                                                          
      IPVT(N) = N                                                       
      IF (A(N,N) .EQ. 0.0D0) INFO = N                                   
      RETURN                                                            
      END                                                               
      SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)                              
      INTEGER LDA,N,IPVT(1),JOB                                         
      DOUBLE PRECISION A(LDA,1),B(1)                                    
!                                                                       
!     DGESL SOLVES THE DOUBLE PRECISION SYSTEM                          
!     A * X = B  OR  TRANS(A) * X = B                                   
!     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.                     
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        A       DOUBLE PRECISION(LDA, N)                               
!                THE OUTPUT FROM DGECO OR DGEFA.                        
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  A .                
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!        IPVT    INTEGER(N)                                             
!                THE PIVOT VECTOR FROM DGECO OR DGEFA.                  
!                                                                       
!        B       DOUBLE PRECISION(N)                                    
!                THE RIGHT HAND SIDE VECTOR.                            
!                                                                       
!        JOB     INTEGER                                                
!                = 0         TO SOLVE  A*X = B ,                        
!                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE            
!                            TRANS(A)  IS THE TRANSPOSE.                
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        B       THE SOLUTION VECTOR  X .                               
!                                                                       
!     ERROR CONDITION                                                   
!                                                                       
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A   
!        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY  
!        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER       
!        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE     
!        CALLED CORRECTLY AND IF DGECO HAS SET RCOND .GT. 0.0           
!        OR DGEFA HAS SET INFO .EQ. 0 .                                 
!                                                                       
!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 
!     WITH  P  COLUMNS                                                  
!           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)                            
!           IF (RCOND IS TOO SMALL) GO TO ...                           
!           DO 10 J = 1, P                                              
!              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)                        
!        10 CONTINUE                                                    
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DAXPY,DDOT                                                   
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION DDOT,T                                           
      INTEGER K,KB,L,NM1                                                
!                                                                       
      NM1 = N - 1                                                       
      IF (JOB .NE. 0) GO TO 50                                          
!                                                                       
!        JOB = 0 , SOLVE  A * X = B                                     
!        FIRST SOLVE  L*Y = B                                           
!                                                                       
         IF (NM1 .LT. 1) GO TO 30                                       
         DO 20 K = 1, NM1                                               
            L = IPVT(K)                                                 
            T = B(L)                                                    
            IF (L .EQ. K) GO TO 10                                      
               B(L) = B(K)                                              
               B(K) = T                                                 
   10       CONTINUE                                                    
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)                       
   20    CONTINUE                                                       
   30    CONTINUE                                                       
!                                                                       
!        NOW SOLVE  U*X = Y                                             
!                                                                       
         DO 40 KB = 1, N                                                
            K = N + 1 - KB                                              
            B(K) = B(K)/A(K,K)                                          
            T = -B(K)                                                   
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)                           
   40    CONTINUE                                                       
      GO TO 100                                                         
   50 CONTINUE                                                          
!                                                                       
!        JOB = NONZERO, SOLVE  TRANS(A) * X = B                         
!        FIRST SOLVE  TRANS(U)*Y = B                                    
!                                                                       
         DO 60 K = 1, N                                                 
            T = DDOT(K-1,A(1,K),1,B(1),1)                               
            B(K) = (B(K) - T)/A(K,K)                                    
   60    CONTINUE                                                       
!                                                                       
!        NOW SOLVE TRANS(L)*X = Y                                       
!                                                                       
         IF (NM1 .LT. 1) GO TO 90                                       
         DO 80 KB = 1, NM1                                              
            K = N - KB                                                  
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)                 
            L = IPVT(K)                                                 
            IF (L .EQ. K) GO TO 70                                      
               T = B(L)                                                 
               B(L) = B(K)                                              
               B(K) = T                                                 
   70       CONTINUE                                                    
   80    CONTINUE                                                       
   90    CONTINUE                                                       
  100 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DGEDI(A,LDA,N,IPVT,DET,WORK,JOB)                       
      INTEGER LDA,N,IPVT(1),JOB                                         
      DOUBLE PRECISION A(LDA,1),DET(2),WORK(1)                          
!                                                                       
!     DGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX            
!     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.                     
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        A       DOUBLE PRECISION(LDA, N)                               
!                THE OUTPUT FROM DGECO OR DGEFA.                        
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  A .                
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!        IPVT    INTEGER(N)                                             
!                THE PIVOT VECTOR FROM DGECO OR DGEFA.                  
!                                                                       
!        WORK    DOUBLE PRECISION(N)                                    
!                WORK VECTOR.  CONTENTS DESTROYED.                      
!                                                                       
!        JOB     INTEGER                                                
!                = 11   BOTH DETERMINANT AND INVERSE.                   
!                = 01   INVERSE ONLY.                                   
!                = 10   DETERMINANT ONLY.                               
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.               
!                OTHERWISE UNCHANGED.                                   
!                                                                       
!        DET     DOUBLE PRECISION(2)                                    
!                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.           
!                OTHERWISE NOT REFERENCED.                              
!                DETERMINANT = DET(1) * 10.0**DET(2)                    
!                WITH  1.0 .LE. DABS(DET(1)) .LT. 10.0                  
!                OR  DET(1) .EQ. 0.0 .                                  
!                                                                       
!     ERROR CONDITION                                                   
!                                                                       
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS     
!        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.           
!        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY      
!        AND IF DGECO HAS SET RCOND .GT. 0.0 OR DGEFA HAS SET           
!        INFO .EQ. 0 .                                                  
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DAXPY,DSCAL,DSWAP                                            
!     FORTRAN DABS,MOD                                                  
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION T                                                
      DOUBLE PRECISION TEN                                              
      INTEGER I,J,K,KB,KP1,L,NM1                                        
!                                                                       
!                                                                       
!     COMPUTE DETERMINANT                                               
!                                                                       
      IF (JOB/10 .EQ. 0) GO TO 70                                       
         DET(1) = 1.0D0                                                 
         DET(2) = 0.0D0                                                 
         TEN = 10.0D0                                                   
         DO 50 I = 1, N                                                 
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)                        
            DET(1) = A(I,I)*DET(1)                                      
!        ...EXIT                                                        
            IF (DET(1) .EQ. 0.0D0) GO TO 60                             
   10       IF (DABS(DET(1)) .GE. 1.0D0) GO TO 20                       
               DET(1) = TEN*DET(1)                                      
               DET(2) = DET(2) - 1.0D0                                  
            GO TO 10                                                    
   20       CONTINUE                                                    
   30       IF (DABS(DET(1)) .LT. TEN) GO TO 40                         
               DET(1) = DET(1)/TEN                                      
               DET(2) = DET(2) + 1.0D0                                  
            GO TO 30                                                    
   40       CONTINUE                                                    
   50    CONTINUE                                                       
   60    CONTINUE                                                       
   70 CONTINUE                                                          
!                                                                       
!     COMPUTE INVERSE(U)                                                
!                                                                       
      IF (MOD(JOB,10) .EQ. 0) GO TO 150                                 
         DO 100 K = 1, N                                                
            A(K,K) = 1.0D0/A(K,K)                                       
            T = -A(K,K)                                                 
            CALL DSCAL(K-1,T,A(1,K),1)                                  
            KP1 = K + 1                                                 
            IF (N .LT. KP1) GO TO 90                                    
            DO 80 J = KP1, N                                            
               T = A(K,J)                                               
               A(K,J) = 0.0D0                                           
               CALL DAXPY(K,T,A(1,K),1,A(1,J),1)                        
   80       CONTINUE                                                    
   90       CONTINUE                                                    
  100    CONTINUE                                                       
!                                                                       
!        FORM INVERSE(U)*INVERSE(L)                                     
!                                                                       
         NM1 = N - 1                                                    
         IF (NM1 .LT. 1) GO TO 140                                      
         DO 130 KB = 1, NM1                                             
            K = N - KB                                                  
            KP1 = K + 1                                                 
            DO 110 I = KP1, N                                           
               WORK(I) = A(I,K)                                         
               A(I,K) = 0.0D0                                           
  110       CONTINUE                                                    
            DO 120 J = KP1, N                                           
               T = WORK(J)                                              
               CALL DAXPY(N,T,A(1,J),1,A(1,K),1)                        
  120       CONTINUE                                                    
            L = IPVT(K)                                                 
            IF (L .NE. K) CALL DSWAP(N,A(1,K),1,A(1,L),1)               
  130    CONTINUE                                                       
  140    CONTINUE                                                       
  150 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)                    
      INTEGER LDA,N,ML,MU,IPVT(1)                                       
      DOUBLE PRECISION ABD(LDA,1),Z(1)                                  
      DOUBLE PRECISION RCOND                                            
!                                                                       
!     DGBCO FACTORS A DOUBLE PRECISION BAND MATRIX BY GAUSSIAN          
!     ELIMINATION AND ESTIMATES THE CONDITION OF THE MATRIX.            
!                                                                       
!     IF  RCOND  IS NOT NEEDED, DGBFA IS SLIGHTLY FASTER.               
!     TO SOLVE  A*X = B , FOLLOW DGBCO BY DGBSL.                        
!     TO COMPUTE  INVERSE(A)*C , FOLLOW DGBCO BY DGBSL.                 
!     TO COMPUTE  DETERMINANT(A) , FOLLOW DGBCO BY DGBDI.               
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        ABD     DOUBLE PRECISION(LDA, N)                               
!                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS      
!                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND   
!                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS         
!                ML+1 THROUGH 2*ML+MU+1 OF  ABD .                       
!                SEE THE COMMENTS BELOW FOR DETAILS.                    
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  ABD .              
!                LDA MUST BE .GE. 2*ML + MU + 1 .                       
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE ORIGINAL MATRIX.                      
!                                                                       
!        ML      INTEGER                                                
!                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.           
!                0 .LE. ML .LT. N .                                     
!                                                                       
!        MU      INTEGER                                                
!                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.           
!                0 .LE. MU .LT. N .                                     
!                MORE EFFICIENT IF  ML .LE. MU .                        
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND         
!                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.          
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
!     BAND STORAGE                                                      
!                                                                       
!           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT      
!           WILL SET UP THE INPUT.                                      
!                                                                       
!                   ML = (BAND WIDTH BELOW THE DIAGONAL)                
!                   MU = (BAND WIDTH ABOVE THE DIAGONAL)                
!                   M = ML + MU + 1                                     
!                   DO 20 J = 1, N                                      
!                      I1 = MAX0(1, J-MU)                               
!                      I2 = MIN0(N, J+ML)                               
!                      DO 10 I = I1, I2                                 
!                         K = I - J + M                                 
!                         ABD(K,J) = A(I,J)                             
!                10    CONTINUE                                         
!                20 CONTINUE                                            
!                                                                       
!           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .         
!           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR      
!           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.            
!           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .    
!           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE            
!           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.          
!                                                                       
!     EXAMPLE..  IF THE ORIGINAL MATRIX IS                              
!                                                                       
!           11 12 13  0  0  0                                           
!           21 22 23 24  0  0                                           
!            0 32 33 34 35  0                                           
!            0  0 43 44 45 46                                           
!            0  0  0 54 55 56                                           
!            0  0  0  0 65 66                                           
!                                                                       
!      THEN  N = 6, ML = 1, MU = 2, LDA .GE. 5  AND ABD SHOULD CONTAIN  
!                                                                       
!            *  *  *  +  +  +  , * = NOT USED                           
!            *  * 13 24 35 46  , + = USED FOR PIVOTING                  
!            * 12 23 34 45 56                                           
!           11 22 33 44 55 66                                           
!           21 32 43 54 65  *                                           
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     LINPACK DGBFA                                                     
!     BLAS DAXPY,DDOT,DSCAL,DASUM                                       
!     FORTRAN DABS,DMAX1,MAX0,MIN0,DSIGN                                
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION DDOT,EK,T,WK,WKM                                 
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM                           
      INTEGER IS,INFO,J,JU,K,KB,KP1,L,LA,LM,LZ,M,MM                     
!                                                                       
!                                                                       
!     COMPUTE 1-NORM OF A                                               
!                                                                       
      ANORM = 0.0D0                                                     
      L = ML + 1                                                        
      IS = L + MU                                                       
      DO 10 J = 1, N                                                    
         ANORM = DMAX1(ANORM,DASUM(L,ABD(IS,J),1))                      
         IF (IS .GT. ML + 1) IS = IS - 1                                
         IF (J .LE. MU) L = L + 1                                       
         IF (J .GE. N - ML) L = L - 1                                   
   10 CONTINUE                                                          
!                                                                       
!     FACTOR                                                            
!                                                                       
      CALL DGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)                             
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
      M = ML + MU + 1                                                   
      JU = 0                                                            
      DO 100 K = 1, N                                                   
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))                      
         IF (DABS(EK-Z(K)) .LE. DABS(ABD(M,K))) GO TO 30                
            S = DABS(ABD(M,K))/DABS(EK-Z(K))                            
            CALL DSCAL(N,S,Z,1)                                         
            EK = S*EK                                                   
   30    CONTINUE                                                       
         WK = EK - Z(K)                                                 
         WKM = -EK - Z(K)                                               
         S = DABS(WK)                                                   
         SM = DABS(WKM)                                                 
         IF (ABD(M,K) .EQ. 0.0D0) GO TO 40                              
            WK = WK/ABD(M,K)                                            
            WKM = WKM/ABD(M,K)                                          
         GO TO 50                                                       
   40    CONTINUE                                                       
            WK = 1.0D0                                                  
            WKM = 1.0D0                                                 
   50    CONTINUE                                                       
         KP1 = K + 1                                                    
         JU = MIN0(MAX0(JU,MU+IPVT(K)),N)                               
         MM = M                                                         
         IF (KP1 .GT. JU) GO TO 90                                      
            DO 60 J = KP1, JU                                           
               MM = MM - 1                                              
               SM = SM + DABS(Z(J)+WKM*ABD(MM,J))                       
               Z(J) = Z(J) + WK*ABD(MM,J)                               
               S = S + DABS(Z(J))                                       
   60       CONTINUE                                                    
            IF (S .GE. SM) GO TO 80                                     
               T = WKM - WK                                             
               WK = WKM                                                 
               MM = M                                                   
               DO 70 J = KP1, JU                                        
                  MM = MM - 1                                           
                  Z(J) = Z(J) + T*ABD(MM,J)                             
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
         LM = MIN0(ML,N-K)                                              
         IF (K .LT. N) Z(K) = Z(K) + DDOT(LM,ABD(M+1,K),1,Z(K+1),1)     
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
         LM = MIN0(ML,N-K)                                              
         IF (K .LT. N) CALL DAXPY(LM,T,ABD(M+1,K),1,Z(K+1),1)           
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
!     SOLVE  U*Z = W                                                    
!                                                                       
      DO 160 KB = 1, N                                                  
         K = N + 1 - KB                                                 
         IF (DABS(Z(K)) .LE. DABS(ABD(M,K))) GO TO 150                  
            S = DABS(ABD(M,K))/DABS(Z(K))                               
            CALL DSCAL(N,S,Z,1)                                         
            YNORM = S*YNORM                                             
  150    CONTINUE                                                       
         IF (ABD(M,K) .NE. 0.0D0) Z(K) = Z(K)/ABD(M,K)                  
         IF (ABD(M,K) .EQ. 0.0D0) Z(K) = 1.0D0                          
         LM = MIN0(K,M) - 1                                             
         LA = M - LM                                                    
         LZ = K - LM                                                    
         T = -Z(K)                                                      
         CALL DAXPY(LM,T,ABD(LA,K),1,Z(LZ),1)                           
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
      SUBROUTINE DGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)                       
      INTEGER LDA,N,ML,MU,IPVT(1),INFO                                  
      DOUBLE PRECISION ABD(LDA,1)                                       
!                                                                       
!     DGBFA FACTORS A DOUBLE PRECISION BAND MATRIX BY ELIMINATION.      
!                                                                       
!     DGBFA IS USUALLY CALLED BY DGBCO, BUT IT CAN BE CALLED            
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.          
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        ABD     DOUBLE PRECISION(LDA, N)                               
!                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS      
!                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND   
!                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS         
!                ML+1 THROUGH 2*ML+MU+1 OF  ABD .                       
!                SEE THE COMMENTS BELOW FOR DETAILS.                    
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  ABD .              
!                LDA MUST BE .GE. 2*ML + MU + 1 .                       
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE ORIGINAL MATRIX.                      
!                                                                       
!        ML      INTEGER                                                
!                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.           
!                0 .LE. ML .LT. N .                                     
!                                                                       
!        MU      INTEGER                                                
!                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.           
!                0 .LE. MU .LT. N .                                     
!                MORE EFFICIENT IF  ML .LE. MU .                        
!     ON RETURN                                                         
!                                                                       
!        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND         
!                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.          
!                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE       
!                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER          
!                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.       
!                                                                       
!        IPVT    INTEGER(N)                                             
!                AN INTEGER VECTOR OF PIVOT INDICES.                    
!                                                                       
!        INFO    INTEGER                                                
!                = 0  NORMAL VALUE.                                     
!                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR       
!                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES        
!                     INDICATE THAT DGBSL WILL DIVIDE BY ZERO IF        
!                     CALLED.  USE  RCOND  IN DGBCO FOR A RELIABLE      
!                     INDICATION OF SINGULARITY.                        
!                                                                       
!     BAND STORAGE                                                      
!                                                                       
!           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT      
!           WILL SET UP THE INPUT.                                      
!                                                                       
!                   ML = (BAND WIDTH BELOW THE DIAGONAL)                
!                   MU = (BAND WIDTH ABOVE THE DIAGONAL)                
!                   M = ML + MU + 1                                     
!                   DO 20 J = 1, N                                      
!                      I1 = MAX0(1, J-MU)                               
!                      I2 = MIN0(N, J+ML)                               
!                      DO 10 I = I1, I2                                 
!                         K = I - J + M                                 
!                         ABD(K,J) = A(I,J)                             
!                10    CONTINUE                                         
!                20 CONTINUE                                            
!                                                                       
!           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .         
!           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR      
!           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.            
!           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .    
!           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE            
!           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.          
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DAXPY,DSCAL,IDAMAX                                           
!     FORTRAN MAX0,MIN0                                                 
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION T                                                
      INTEGER I,IDAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1             
!                                                                       
!                                                                       
      M = ML + MU + 1                                                   
      INFO = 0                                                          
!                                                                       
!     ZERO INITIAL FILL-IN COLUMNS                                      
!                                                                       
      J0 = MU + 2                                                       
      J1 = MIN0(N,M) - 1                                                
      IF (J1 .LT. J0) GO TO 30                                          
      DO 20 JZ = J0, J1                                                 
         I0 = M + 1 - JZ                                                
         DO 10 I = I0, ML                                               
            ABD(I,JZ) = 0.0D0                                           
   10    CONTINUE                                                       
   20 CONTINUE                                                          
   30 CONTINUE                                                          
      JZ = J1                                                           
      JU = 0                                                            
!                                                                       
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING                        
!                                                                       
      NM1 = N - 1                                                       
      IF (NM1 .LT. 1) GO TO 130                                         
      DO 120 K = 1, NM1                                                 
         KP1 = K + 1                                                    
!                                                                       
!        ZERO NEXT FILL-IN COLUMN                                       
!                                                                       
         JZ = JZ + 1                                                    
         IF (JZ .GT. N) GO TO 50                                        
         IF (ML .LT. 1) GO TO 50                                        
            DO 40 I = 1, ML                                             
               ABD(I,JZ) = 0.0D0                                        
   40       CONTINUE                                                    
   50    CONTINUE                                                       
!                                                                       
!        FIND L = PIVOT INDEX                                           
!                                                                       
         LM = MIN0(ML,N-K)                                              
         L = IDAMAX(LM+1,ABD(M,K),1) + M - 1                            
         IPVT(K) = L + K - M                                            
!                                                                       
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED          
!                                                                       
         IF (ABD(L,K) .EQ. 0.0D0) GO TO 100                             
!                                                                       
!           INTERCHANGE IF NECESSARY                                    
!                                                                       
            IF (L .EQ. M) GO TO 60                                      
               T = ABD(L,K)                                             
               ABD(L,K) = ABD(M,K)                                      
               ABD(M,K) = T                                             
   60       CONTINUE                                                    
!                                                                       
!           COMPUTE MULTIPLIERS                                         
!                                                                       
            T = -1.0D0/ABD(M,K)                                         
            CALL DSCAL(LM,T,ABD(M+1,K),1)                               
!                                                                       
!           ROW ELIMINATION WITH COLUMN INDEXING                        
!                                                                       
            JU = MIN0(MAX0(JU,MU+IPVT(K)),N)                            
            MM = M                                                      
            IF (JU .LT. KP1) GO TO 90                                   
            DO 80 J = KP1, JU                                           
               L = L - 1                                                
               MM = MM - 1                                              
               T = ABD(L,J)                                             
               IF (L .EQ. MM) GO TO 70                                  
                  ABD(L,J) = ABD(MM,J)                                  
                  ABD(MM,J) = T                                         
   70          CONTINUE                                                 
               CALL DAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)              
   80       CONTINUE                                                    
   90       CONTINUE                                                    
         GO TO 110                                                      
  100    CONTINUE                                                       
            INFO = K                                                    
  110    CONTINUE                                                       
  120 CONTINUE                                                          
  130 CONTINUE                                                          
      IPVT(N) = N                                                       
      IF (ABD(M,N) .EQ. 0.0D0) INFO = N                                 
      RETURN                                                            
      END                                                               
      SUBROUTINE DGBSL(ABD,LDA,N,ML,MU,IPVT,B,JOB)                      
      INTEGER LDA,N,ML,MU,IPVT(1),JOB                                   
      DOUBLE PRECISION ABD(LDA,1),B(1)                                  
!                                                                       
!     DGBSL SOLVES THE DOUBLE PRECISION BAND SYSTEM                     
!     A * X = B  OR  TRANS(A) * X = B                                   
!     USING THE FACTORS COMPUTED BY DGBCO OR DGBFA.                     
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        ABD     DOUBLE PRECISION(LDA, N)                               
!                THE OUTPUT FROM DGBCO OR DGBFA.                        
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  ABD .              
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE ORIGINAL MATRIX.                      
!                                                                       
!        ML      INTEGER                                                
!                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.           
!                                                                       
!        MU      INTEGER                                                
!                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.           
!                                                                       
!        IPVT    INTEGER(N)                                             
!                THE PIVOT VECTOR FROM DGBCO OR DGBFA.                  
!                                                                       
!        B       DOUBLE PRECISION(N)                                    
!                THE RIGHT HAND SIDE VECTOR.                            
!                                                                       
!        JOB     INTEGER                                                
!                = 0         TO SOLVE  A*X = B ,                        
!                = NONZERO   TO SOLVE  TRANS(A)*X = B , WHERE           
!                            TRANS(A)  IS THE TRANSPOSE.                
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        B       THE SOLUTION VECTOR  X .                               
!                                                                       
!     ERROR CONDITION                                                   
!                                                                       
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A   
!        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY  
!        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER       
!        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE     
!        CALLED CORRECTLY AND IF DGBCO HAS SET RCOND .GT. 0.0           
!        OR DGBFA HAS SET INFO .EQ. 0 .                                 
!                                                                       
!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 
!     WITH  P  COLUMNS                                                  
!           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)                    
!           IF (RCOND IS TOO SMALL) GO TO ...                           
!           DO 10 J = 1, P                                              
!              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)                
!        10 CONTINUE                                                    
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DAXPY,DDOT                                                   
!     FORTRAN MIN0                                                      
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION DDOT,T                                           
      INTEGER K,KB,L,LA,LB,LM,M,NM1                                     
!                                                                       
      M = MU + ML + 1                                                   
      NM1 = N - 1                                                       
      IF (JOB .NE. 0) GO TO 50                                          
!                                                                       
!        JOB = 0 , SOLVE  A * X = B                                     
!        FIRST SOLVE L*Y = B                                            
!                                                                       
         IF (ML .EQ. 0) GO TO 30                                        
         IF (NM1 .LT. 1) GO TO 30                                       
            DO 20 K = 1, NM1                                            
               LM = MIN0(ML,N-K)                                        
               L = IPVT(K)                                              
               T = B(L)                                                 
               IF (L .EQ. K) GO TO 10                                   
                  B(L) = B(K)                                           
                  B(K) = T                                              
   10          CONTINUE                                                 
               CALL DAXPY(LM,T,ABD(M+1,K),1,B(K+1),1)                   
   20       CONTINUE                                                    
   30    CONTINUE                                                       
!                                                                       
!        NOW SOLVE  U*X = Y                                             
!                                                                       
         DO 40 KB = 1, N                                                
            K = N + 1 - KB                                              
            B(K) = B(K)/ABD(M,K)                                        
            LM = MIN0(K,M) - 1                                          
            LA = M - LM                                                 
            LB = K - LM                                                 
            T = -B(K)                                                   
            CALL DAXPY(LM,T,ABD(LA,K),1,B(LB),1)                        
   40    CONTINUE                                                       
      GO TO 100                                                         
   50 CONTINUE                                                          
!                                                                       
!        JOB = NONZERO, SOLVE  TRANS(A) * X = B                         
!        FIRST SOLVE  TRANS(U)*Y = B                                    
!                                                                       
         DO 60 K = 1, N                                                 
            LM = MIN0(K,M) - 1                                          
            LA = M - LM                                                 
            LB = K - LM                                                 
            T = DDOT(LM,ABD(LA,K),1,B(LB),1)                            
            B(K) = (B(K) - T)/ABD(M,K)                                  
   60    CONTINUE                                                       
!                                                                       
!        NOW SOLVE TRANS(L)*X = Y                                       
!                                                                       
         IF (ML .EQ. 0) GO TO 90                                        
         IF (NM1 .LT. 1) GO TO 90                                       
            DO 80 KB = 1, NM1                                           
               K = N - KB                                               
               LM = MIN0(ML,N-K)                                        
               B(K) = B(K) + DDOT(LM,ABD(M+1,K),1,B(K+1),1)             
               L = IPVT(K)                                              
               IF (L .EQ. K) GO TO 70                                   
                  T = B(L)                                              
                  B(L) = B(K)                                           
                  B(K) = T                                              
   70          CONTINUE                                                 
   80       CONTINUE                                                    
   90    CONTINUE                                                       
  100 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DGBDI(ABD,LDA,N,ML,MU,IPVT,DET)                        
      INTEGER LDA,N,ML,MU,IPVT(1)                                       
      DOUBLE PRECISION ABD(LDA,1),DET(2)                                
!                                                                       
!     DGBDI COMPUTES THE DETERMINANT OF A BAND MATRIX                   
!     USING THE FACTORS COMPUTED BY DGBCO OR DGBFA.                     
!     IF THE INVERSE IS NEEDED, USE DGBSL  N  TIMES.                    
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        ABD     DOUBLE PRECISION(LDA, N)                               
!                THE OUTPUT FROM DGBCO OR DGBFA.                        
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  ABD .              
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE ORIGINAL MATRIX.                      
!                                                                       
!        ML      INTEGER                                                
!                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.           
!                                                                       
!        MU      INTEGER                                                
!                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.           
!                                                                       
!        IPVT    INTEGER(N)                                             
!                THE PIVOT VECTOR FROM DGBCO OR DGBFA.                  
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        DET     DOUBLE PRECISION(2)                                    
!                DETERMINANT OF ORIGINAL MATRIX.                        
!                DETERMINANT = DET(1) * 10.0**DET(2)                    
!                WITH  1.0 .LE. DABS(DET(1)) .LT. 10.0                  
!                OR  DET(1) = 0.0 .                                     
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     FORTRAN DABS                                                      
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION TEN                                              
      INTEGER I,M                                                       
!                                                                       
!                                                                       
      M = ML + MU + 1                                                   
      DET(1) = 1.0D0                                                    
      DET(2) = 0.0D0                                                    
      TEN = 10.0D0                                                      
      DO 50 I = 1, N                                                    
         IF (IPVT(I) .NE. I) DET(1) = -DET(1)                           
         DET(1) = ABD(M,I)*DET(1)                                       
!     ...EXIT                                                           
         IF (DET(1) .EQ. 0.0D0) GO TO 60                                
   10    IF (DABS(DET(1)) .GE. 1.0D0) GO TO 20                          
            DET(1) = TEN*DET(1)                                         
            DET(2) = DET(2) - 1.0D0                                     
         GO TO 10                                                       
   20    CONTINUE                                                       
   30    IF (DABS(DET(1)) .LT. TEN) GO TO 40                            
            DET(1) = DET(1)/TEN                                         
            DET(2) = DET(2) + 1.0D0                                     
         GO TO 30                                                       
   40    CONTINUE                                                       
   50 CONTINUE                                                          
   60 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DPOCO(A,LDA,N,RCOND,Z,INFO)                            
      INTEGER LDA,N,INFO                                                
      DOUBLE PRECISION A(LDA,1),Z(1)                                    
      DOUBLE PRECISION RCOND                                            
!                                                                       
!     DPOCO FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE      
!     MATRIX AND ESTIMATES THE CONDITION OF THE MATRIX.                 
!                                                                       
!     IF  RCOND  IS NOT NEEDED, DPOFA IS SLIGHTLY FASTER.               
!     TO SOLVE  A*X = B , FOLLOW DPOCO BY DPOSL.                        
!     TO COMPUTE  INVERSE(A)*C , FOLLOW DPOCO BY DPOSL.                 
!     TO COMPUTE  DETERMINANT(A) , FOLLOW DPOCO BY DPODI.               
!     TO COMPUTE  INVERSE(A) , FOLLOW DPOCO BY DPODI.                   
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        A       DOUBLE PRECISION(LDA, N)                               
!                THE SYMMETRIC MATRIX TO BE FACTORED.  ONLY THE         
!                DIAGONAL AND UPPER TRIANGLE ARE USED.                  
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  A .                
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        A       AN UPPER TRIANGULAR MATRIX  R  SO THAT  A = TRANS(R)*R 
!                WHERE  TRANS(R)  IS THE TRANSPOSE.                     
!                THE STRICT LOWER TRIANGLE IS UNALTERED.                
!                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.   
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
!                UNDERFLOWS.  IF INFO .NE. 0 , RCOND IS UNCHANGED.      
!                                                                       
!        Z       DOUBLE PRECISION(N)                                    
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.  
!                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS      
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT           
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .                    
!                IF  INFO .NE. 0 , Z  IS UNCHANGED.                     
!                                                                       
!        INFO    INTEGER                                                
!                = 0  FOR NORMAL RETURN.                                
!                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR    
!                     OF ORDER  K  IS NOT POSITIVE DEFINITE.            
!                                                                       
!     LINPACK.  THIS VERSION DATED 08/14/78 .                           
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     LINPACK DPOFA                                                     
!     BLAS DAXPY,DDOT,DSCAL,DASUM                                       
!     FORTRAN DABS,DMAX1,DREAL,DSIGN                                    
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION DDOT,EK,T,WK,WKM                                 
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM                           
      INTEGER I,J,JM1,K,KB,KP1                                          
!                                                                       
!                                                                       
!     FIND NORM OF A USING ONLY UPPER HALF                              
!                                                                       
      DO 30 J = 1, N                                                    
         Z(J) = DASUM(J,A(1,J),1)                                       
         JM1 = J - 1                                                    
         IF (JM1 .LT. 1) GO TO 20                                       
         DO 10 I = 1, JM1                                               
            Z(I) = Z(I) + DABS(A(I,J))                                  
   10    CONTINUE                                                       
   20    CONTINUE                                                       
   30 CONTINUE                                                          
      ANORM = 0.0D0                                                     
      DO 40 J = 1, N                                                    
         ANORM = DMAX1(ANORM,Z(J))                                      
   40 CONTINUE                                                          
!                                                                       
!     FACTOR                                                            
!                                                                       
      CALL DPOFA(A,LDA,N,INFO)                                          
      IF (INFO .NE. 0) GO TO 180                                        
!                                                                       
!        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .           
!        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .      
!        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL        
!        GROWTH IN THE ELEMENTS OF W  WHERE  TRANS(R)*W = E .           
!        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.         
!                                                                       
!        SOLVE TRANS(R)*W = E                                           
!                                                                       
         EK = 1.0D0                                                     
         DO 50 J = 1, N                                                 
            Z(J) = 0.0D0                                                
   50    CONTINUE                                                       
         DO 110 K = 1, N                                                
            IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))                   
            IF (DABS(EK-Z(K)) .LE. A(K,K)) GO TO 60                     
               S = A(K,K)/DABS(EK-Z(K))                                 
               CALL DSCAL(N,S,Z,1)                                      
               EK = S*EK                                                
   60       CONTINUE                                                    
            WK = EK - Z(K)                                              
            WKM = -EK - Z(K)                                            
            S = DABS(WK)                                                
            SM = DABS(WKM)                                              
            WK = WK/A(K,K)                                              
            WKM = WKM/A(K,K)                                            
            KP1 = K + 1                                                 
            IF (KP1 .GT. N) GO TO 100                                   
               DO 70 J = KP1, N                                         
                  SM = SM + DABS(Z(J)+WKM*A(K,J))                       
                  Z(J) = Z(J) + WK*A(K,J)                               
                  S = S + DABS(Z(J))                                    
   70          CONTINUE                                                 
               IF (S .GE. SM) GO TO 90                                  
                  T = WKM - WK                                          
                  WK = WKM                                              
                  DO 80 J = KP1, N                                      
                     Z(J) = Z(J) + T*A(K,J)                             
   80             CONTINUE                                              
   90          CONTINUE                                                 
  100       CONTINUE                                                    
            Z(K) = WK                                                   
  110    CONTINUE                                                       
         S = 1.0D0/DASUM(N,Z,1)                                         
         CALL DSCAL(N,S,Z,1)                                            
!                                                                       
!        SOLVE R*Y = W                                                  
!                                                                       
         DO 130 KB = 1, N                                               
            K = N + 1 - KB                                              
            IF (DABS(Z(K)) .LE. A(K,K)) GO TO 120                       
               S = A(K,K)/DABS(Z(K))                                    
               CALL DSCAL(N,S,Z,1)                                      
  120       CONTINUE                                                    
            Z(K) = Z(K)/A(K,K)                                          
            T = -Z(K)                                                   
            CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)                           
  130    CONTINUE                                                       
         S = 1.0D0/DASUM(N,Z,1)                                         
         CALL DSCAL(N,S,Z,1)                                            
!                                                                       
         YNORM = 1.0D0                                                  
!                                                                       
!        SOLVE TRANS(R)*V = Y                                           
!                                                                       
         DO 150 K = 1, N                                                
            Z(K) = Z(K) - DDOT(K-1,A(1,K),1,Z(1),1)                     
            IF (DABS(Z(K)) .LE. A(K,K)) GO TO 140                       
               S = A(K,K)/DABS(Z(K))                                    
               CALL DSCAL(N,S,Z,1)                                      
               YNORM = S*YNORM                                          
  140       CONTINUE                                                    
            Z(K) = Z(K)/A(K,K)                                          
  150    CONTINUE                                                       
         S = 1.0D0/DASUM(N,Z,1)                                         
         CALL DSCAL(N,S,Z,1)                                            
         YNORM = S*YNORM                                                
!                                                                       
!        SOLVE R*Z = V                                                  
!                                                                       
         DO 170 KB = 1, N                                               
            K = N + 1 - KB                                              
            IF (DABS(Z(K)) .LE. A(K,K)) GO TO 160                       
               S = A(K,K)/DABS(Z(K))                                    
               CALL DSCAL(N,S,Z,1)                                      
               YNORM = S*YNORM                                          
  160       CONTINUE                                                    
            Z(K) = Z(K)/A(K,K)                                          
            T = -Z(K)                                                   
            CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)                           
  170    CONTINUE                                                       
!        MAKE ZNORM = 1.0                                               
         S = 1.0D0/DASUM(N,Z,1)                                         
         CALL DSCAL(N,S,Z,1)                                            
         YNORM = S*YNORM                                                
!                                                                       
         IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM                      
         IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0                            
  180 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DPOFA(A,LDA,N,INFO)                                    
      INTEGER LDA,N,INFO                                                
      DOUBLE PRECISION A(LDA,1)                                         
!                                                                       
!     DPOFA FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE      
!     MATRIX.                                                           
!                                                                       
!     DPOFA IS USUALLY CALLED BY DPOCO, BUT IT CAN BE CALLED            
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.          
!     (TIME FOR DPOCO) = (1 + 18/N)*(TIME FOR DPOFA) .                  
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        A       DOUBLE PRECISION(LDA, N)                               
!                THE SYMMETRIC MATRIX TO BE FACTORED.  ONLY THE         
!                DIAGONAL AND UPPER TRIANGLE ARE USED.                  
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  A .                
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        A       AN UPPER TRIANGULAR MATRIX  R  SO THAT  A = TRANS(R)*R 
!                WHERE  TRANS(R)  IS THE TRANSPOSE.                     
!                THE STRICT LOWER TRIANGLE IS UNALTERED.                
!                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.   
!                                                                       
!        INFO    INTEGER                                                
!                = 0  FOR NORMAL RETURN.                                
!                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR    
!                     OF ORDER  K  IS NOT POSITIVE DEFINITE.            
!                                                                       
!     LINPACK.  THIS VERSION DATED 08/14/78 .                           
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DDOT                                                         
!     FORTRAN DSQRT                                                     
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION DDOT,T                                           
      DOUBLE PRECISION S                                                
      INTEGER J,JM1,K                                                   
!     BEGIN BLOCK WITH ...EXITS TO 40                                   
!                                                                       
!                                                                       
         DO 30 J = 1, N                                                 
            INFO = J                                                    
            S = 0.0D0                                                   
            JM1 = J - 1                                                 
            IF (JM1 .LT. 1) GO TO 20                                    
            DO 10 K = 1, JM1                                            
               T = A(K,J) - DDOT(K-1,A(1,K),1,A(1,J),1)                 
               T = T/A(K,K)                                             
               A(K,J) = T                                               
               S = S + T*T                                              
   10       CONTINUE                                                    
   20       CONTINUE                                                    
            S = A(J,J) - S                                              
!     ......EXIT                                                        
            IF (S .LE. 0.0D0) GO TO 40                                  
            A(J,J) = DSQRT(S)                                           
   30    CONTINUE                                                       
         INFO = 0                                                       
   40 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DPOSL(A,LDA,N,B)                                       
      INTEGER LDA,N                                                     
      DOUBLE PRECISION A(LDA,1),B(1)                                    
!                                                                       
!     DPOSL SOLVES THE DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE     
!     SYSTEM A * X = B                                                  
!     USING THE FACTORS COMPUTED BY DPOCO OR DPOFA.                     
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        A       DOUBLE PRECISION(LDA, N)                               
!                THE OUTPUT FROM DPOCO OR DPOFA.                        
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  A .                
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!        B       DOUBLE PRECISION(N)                                    
!                THE RIGHT HAND SIDE VECTOR.                            
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        B       THE SOLUTION VECTOR  X .                               
!                                                                       
!     ERROR CONDITION                                                   
!                                                                       
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS     
!        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES            
!        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE    
!        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED    
!        CORRECTLY AND  INFO .EQ. 0 .                                   
!                                                                       
!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 
!     WITH  P  COLUMNS                                                  
!           CALL DPOCO(A,LDA,N,RCOND,Z,INFO)                            
!           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ...          
!           DO 10 J = 1, P                                              
!              CALL DPOSL(A,LDA,N,C(1,J))                               
!        10 CONTINUE                                                    
!                                                                       
!     LINPACK.  THIS VERSION DATED 08/14/78 .                           
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DAXPY,DDOT                                                   
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION DDOT,T                                           
      INTEGER K,KB                                                      
!                                                                       
!     SOLVE TRANS(R)*Y = B                                              
!                                                                       
      DO 10 K = 1, N                                                    
         T = DDOT(K-1,A(1,K),1,B(1),1)                                  
         B(K) = (B(K) - T)/A(K,K)                                       
   10 CONTINUE                                                          
!                                                                       
!     SOLVE R*X = Y                                                     
!                                                                       
      DO 20 KB = 1, N                                                   
         K = N + 1 - KB                                                 
         B(K) = B(K)/A(K,K)                                             
         T = -B(K)                                                      
         CALL DAXPY(K-1,T,A(1,K),1,B(1),1)                              
   20 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DPODI(A,LDA,N,DET,JOB)                                 
      INTEGER LDA,N,JOB                                                 
      DOUBLE PRECISION A(LDA,1)                                         
      DOUBLE PRECISION DET(2)                                           
!                                                                       
!     DPODI COMPUTES THE DETERMINANT AND INVERSE OF A CERTAIN           
!     DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE MATRIX (SEE BELOW)   
!     USING THE FACTORS COMPUTED BY DPOCO, DPOFA OR DQRDC.              
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        A       DOUBLE PRECISION(LDA, N)                               
!                THE OUTPUT  A  FROM DPOCO OR DPOFA                     
!                OR THE OUTPUT  X  FROM DQRDC.                          
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  A .                
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!        JOB     INTEGER                                                
!                = 11   BOTH DETERMINANT AND INVERSE.                   
!                = 01   INVERSE ONLY.                                   
!                = 10   DETERMINANT ONLY.                               
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        A       IF DPOCO OR DPOFA WAS USED TO FACTOR  A  THEN          
!                DPODI PRODUCES THE UPPER HALF OF INVERSE(A) .          
!                IF DQRDC WAS USED TO DECOMPOSE  X  THEN                
!                DPODI PRODUCES THE UPPER HALF OF INVERSE(TRANS(X)*X)   
!                WHERE TRANS(X) IS THE TRANSPOSE.                       
!                ELEMENTS OF  A  BELOW THE DIAGONAL ARE UNCHANGED.      
!                IF THE UNITS DIGIT OF JOB IS ZERO,  A  IS UNCHANGED.   
!                                                                       
!        DET     DOUBLE PRECISION(2)                                    
!                DETERMINANT OF  A  OR OF  TRANS(X)*X  IF REQUESTED.    
!                OTHERWISE NOT REFERENCED.                              
!                DETERMINANT = DET(1) * 10.0**DET(2)                    
!                WITH  1.0 .LE. DET(1) .LT. 10.0                        
!                OR  DET(1) .EQ. 0.0 .                                  
!                                                                       
!     ERROR CONDITION                                                   
!                                                                       
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS     
!        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.           
!        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY      
!        AND IF DPOCO OR DPOFA HAS SET INFO .EQ. 0 .                    
!                                                                       
!     LINPACK.  THIS VERSION DATED 08/14/78 .                           
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DAXPY,DSCAL                                                  
!     FORTRAN MOD                                                       
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION T                                                
      DOUBLE PRECISION S                                                
      INTEGER I,J,JM1,K,KP1                                             
!                                                                       
!     COMPUTE DETERMINANT                                               
!                                                                       
      IF (JOB/10 .EQ. 0) GO TO 70                                       
         DET(1) = 1.0D0                                                 
         DET(2) = 0.0D0                                                 
         S = 10.0D0                                                     
         DO 50 I = 1, N                                                 
            DET(1) = A(I,I)**2*DET(1)                                   
!        ...EXIT                                                        
            IF (DET(1) .EQ. 0.0D0) GO TO 60                             
   10       IF (DET(1) .GE. 1.0D0) GO TO 20                             
               DET(1) = S*DET(1)                                        
               DET(2) = DET(2) - 1.0D0                                  
            GO TO 10                                                    
   20       CONTINUE                                                    
   30       IF (DET(1) .LT. S) GO TO 40                                 
               DET(1) = DET(1)/S                                        
               DET(2) = DET(2) + 1.0D0                                  
            GO TO 30                                                    
   40       CONTINUE                                                    
   50    CONTINUE                                                       
   60    CONTINUE                                                       
   70 CONTINUE                                                          
!                                                                       
!     COMPUTE INVERSE(R)                                                
!                                                                       
      IF (MOD(JOB,10) .EQ. 0) GO TO 140                                 
         DO 100 K = 1, N                                                
            A(K,K) = 1.0D0/A(K,K)                                       
            T = -A(K,K)                                                 
            CALL DSCAL(K-1,T,A(1,K),1)                                  
            KP1 = K + 1                                                 
            IF (N .LT. KP1) GO TO 90                                    
            DO 80 J = KP1, N                                            
               T = A(K,J)                                               
               A(K,J) = 0.0D0                                           
               CALL DAXPY(K,T,A(1,K),1,A(1,J),1)                        
   80       CONTINUE                                                    
   90       CONTINUE                                                    
  100    CONTINUE                                                       
!                                                                       
!        FORM  INVERSE(R) * TRANS(INVERSE(R))                           
!                                                                       
         DO 130 J = 1, N                                                
            JM1 = J - 1                                                 
            IF (JM1 .LT. 1) GO TO 120                                   
            DO 110 K = 1, JM1                                           
               T = A(K,J)                                               
               CALL DAXPY(K,T,A(1,J),1,A(1,K),1)                        
  110       CONTINUE                                                    
  120       CONTINUE                                                    
            T = A(J,J)                                                  
            CALL DSCAL(J,T,A(1,J),1)                                    
  130    CONTINUE                                                       
  140 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DPPCO(AP,N,RCOND,Z,INFO)                               
      INTEGER N,INFO                                                    
      DOUBLE PRECISION AP(1),Z(1)                                       
      DOUBLE PRECISION RCOND                                            
!                                                                       
!     DPPCO FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE      
!     MATRIX STORED IN PACKED FORM                                      
!     AND ESTIMATES THE CONDITION OF THE MATRIX.                        
!                                                                       
!     IF  RCOND  IS NOT NEEDED, DPPFA IS SLIGHTLY FASTER.               
!     TO SOLVE  A*X = B , FOLLOW DPPCO BY DPPSL.                        
!     TO COMPUTE  INVERSE(A)*C , FOLLOW DPPCO BY DPPSL.                 
!     TO COMPUTE  DETERMINANT(A) , FOLLOW DPPCO BY DPPDI.               
!     TO COMPUTE  INVERSE(A) , FOLLOW DPPCO BY DPPDI.                   
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        AP      DOUBLE PRECISION (N*(N+1)/2)                           
!                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE        
!                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY  
!                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .      
!                SEE COMMENTS BELOW FOR DETAILS.                        
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        AP      AN UPPER TRIANGULAR MATRIX  R , STORED IN PACKED       
!                FORM, SO THAT  A = TRANS(R)*R .                        
!                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.   
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
!                UNDERFLOWS.  IF INFO .NE. 0 , RCOND IS UNCHANGED.      
!                                                                       
!        Z       DOUBLE PRECISION(N)                                    
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.  
!                IF  A  IS SINGULAR TO WORKING PRECISION, THEN  Z  IS   
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT           
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .                    
!                IF  INFO .NE. 0 , Z  IS UNCHANGED.                     
!                                                                       
!        INFO    INTEGER                                                
!                = 0  FOR NORMAL RETURN.                                
!                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR    
!                     OF ORDER  K  IS NOT POSITIVE DEFINITE.            
!                                                                       
!     PACKED STORAGE                                                    
!                                                                       
!          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER            
!          TRIANGLE OF A SYMMETRIC MATRIX.                              
!                                                                       
!                K = 0                                                  
!                DO 20 J = 1, N                                         
!                   DO 10 I = 1, J                                      
!                      K = K + 1                                        
!                      AP(K) = A(I,J)                                   
!             10    CONTINUE                                            
!             20 CONTINUE                                               
!                                                                       
!     LINPACK.  THIS VERSION DATED 08/14/78 .                           
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     LINPACK DPPFA                                                     
!     BLAS DAXPY,DDOT,DSCAL,DASUM                                       
!     FORTRAN DABS,DMAX1,DREAL,DSIGN                                    
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION DDOT,EK,T,WK,WKM                                 
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM                           
      INTEGER I,IJ,J,JM1,J1,K,KB,KJ,KK,KP1                              
!                                                                       
!                                                                       
!     FIND NORM OF A                                                    
!                                                                       
      J1 = 1                                                            
      DO 30 J = 1, N                                                    
         Z(J) = DASUM(J,AP(J1),1)                                       
         IJ = J1                                                        
         J1 = J1 + J                                                    
         JM1 = J - 1                                                    
         IF (JM1 .LT. 1) GO TO 20                                       
         DO 10 I = 1, JM1                                               
            Z(I) = Z(I) + DABS(AP(IJ))                                  
            IJ = IJ + 1                                                 
   10    CONTINUE                                                       
   20    CONTINUE                                                       
   30 CONTINUE                                                          
      ANORM = 0.0D0                                                     
      DO 40 J = 1, N                                                    
         ANORM = DMAX1(ANORM,Z(J))                                      
   40 CONTINUE                                                          
!                                                                       
!     FACTOR                                                            
!                                                                       
      CALL DPPFA(AP,N,INFO)                                             
      IF (INFO .NE. 0) GO TO 180                                        
!                                                                       
!        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .           
!        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .      
!        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL        
!        GROWTH IN THE ELEMENTS OF W  WHERE  TRANS(R)*W = E .           
!        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.         
!                                                                       
!        SOLVE TRANS(R)*W = E                                           
!                                                                       
         EK = 1.0D0                                                     
         DO 50 J = 1, N                                                 
            Z(J) = 0.0D0                                                
   50    CONTINUE                                                       
         KK = 0                                                         
         DO 110 K = 1, N                                                
            KK = KK + K                                                 
            IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))                   
            IF (DABS(EK-Z(K)) .LE. AP(KK)) GO TO 60                     
               S = AP(KK)/DABS(EK-Z(K))                                 
               CALL DSCAL(N,S,Z,1)                                      
               EK = S*EK                                                
   60       CONTINUE                                                    
            WK = EK - Z(K)                                              
            WKM = -EK - Z(K)                                            
            S = DABS(WK)                                                
            SM = DABS(WKM)                                              
            WK = WK/AP(KK)                                              
            WKM = WKM/AP(KK)                                            
            KP1 = K + 1                                                 
            KJ = KK + K                                                 
            IF (KP1 .GT. N) GO TO 100                                   
               DO 70 J = KP1, N                                         
                  SM = SM + DABS(Z(J)+WKM*AP(KJ))                       
                  Z(J) = Z(J) + WK*AP(KJ)                               
                  S = S + DABS(Z(J))                                    
                  KJ = KJ + J                                           
   70          CONTINUE                                                 
               IF (S .GE. SM) GO TO 90                                  
                  T = WKM - WK                                          
                  WK = WKM                                              
                  KJ = KK + K                                           
                  DO 80 J = KP1, N                                      
                     Z(J) = Z(J) + T*AP(KJ)                             
                     KJ = KJ + J                                        
   80             CONTINUE                                              
   90          CONTINUE                                                 
  100       CONTINUE                                                    
            Z(K) = WK                                                   
  110    CONTINUE                                                       
         S = 1.0D0/DASUM(N,Z,1)                                         
         CALL DSCAL(N,S,Z,1)                                            
!                                                                       
!        SOLVE R*Y = W                                                  
!                                                                       
         DO 130 KB = 1, N                                               
            K = N + 1 - KB                                              
            IF (DABS(Z(K)) .LE. AP(KK)) GO TO 120                       
               S = AP(KK)/DABS(Z(K))                                    
               CALL DSCAL(N,S,Z,1)                                      
  120       CONTINUE                                                    
            Z(K) = Z(K)/AP(KK)                                          
            KK = KK - K                                                 
            T = -Z(K)                                                   
            CALL DAXPY(K-1,T,AP(KK+1),1,Z(1),1)                         
  130    CONTINUE                                                       
         S = 1.0D0/DASUM(N,Z,1)                                         
         CALL DSCAL(N,S,Z,1)                                            
!                                                                       
         YNORM = 1.0D0                                                  
!                                                                       
!        SOLVE TRANS(R)*V = Y                                           
!                                                                       
         DO 150 K = 1, N                                                
            Z(K) = Z(K) - DDOT(K-1,AP(KK+1),1,Z(1),1)                   
            KK = KK + K                                                 
            IF (DABS(Z(K)) .LE. AP(KK)) GO TO 140                       
               S = AP(KK)/DABS(Z(K))                                    
               CALL DSCAL(N,S,Z,1)                                      
               YNORM = S*YNORM                                          
  140       CONTINUE                                                    
            Z(K) = Z(K)/AP(KK)                                          
  150    CONTINUE                                                       
         S = 1.0D0/DASUM(N,Z,1)                                         
         CALL DSCAL(N,S,Z,1)                                            
         YNORM = S*YNORM                                                
!                                                                       
!        SOLVE R*Z = V                                                  
!                                                                       
         DO 170 KB = 1, N                                               
            K = N + 1 - KB                                              
            IF (DABS(Z(K)) .LE. AP(KK)) GO TO 160                       
               S = AP(KK)/DABS(Z(K))                                    
               CALL DSCAL(N,S,Z,1)                                      
               YNORM = S*YNORM                                          
  160       CONTINUE                                                    
            Z(K) = Z(K)/AP(KK)                                          
            KK = KK - K                                                 
            T = -Z(K)                                                   
            CALL DAXPY(K-1,T,AP(KK+1),1,Z(1),1)                         
  170    CONTINUE                                                       
!        MAKE ZNORM = 1.0                                               
         S = 1.0D0/DASUM(N,Z,1)                                         
         CALL DSCAL(N,S,Z,1)                                            
         YNORM = S*YNORM                                                
!                                                                       
         IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM                      
         IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0                            
  180 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DPPFA(AP,N,INFO)                                       
      INTEGER N,INFO                                                    
      DOUBLE PRECISION AP(1)                                            
!                                                                       
!     DPPFA FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE      
!     MATRIX STORED IN PACKED FORM.                                     
!                                                                       
!     DPPFA IS USUALLY CALLED BY DPPCO, BUT IT CAN BE CALLED            
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.          
!     (TIME FOR DPPCO) = (1 + 18/N)*(TIME FOR DPPFA) .                  
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        AP      DOUBLE PRECISION (N*(N+1)/2)                           
!                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE        
!                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY  
!                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .      
!                SEE COMMENTS BELOW FOR DETAILS.                        
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        AP      AN UPPER TRIANGULAR MATRIX  R , STORED IN PACKED       
!                FORM, SO THAT  A = TRANS(R)*R .                        
!                                                                       
!        INFO    INTEGER                                                
!                = 0  FOR NORMAL RETURN.                                
!                = K  IF THE LEADING MINOR OF ORDER  K  IS NOT          
!                     POSITIVE DEFINITE.                                
!                                                                       
!                                                                       
!     PACKED STORAGE                                                    
!                                                                       
!          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER            
!          TRIANGLE OF A SYMMETRIC MATRIX.                              
!                                                                       
!                K = 0                                                  
!                DO 20 J = 1, N                                         
!                   DO 10 I = 1, J                                      
!                      K = K + 1                                        
!                      AP(K) = A(I,J)                                   
!             10    CONTINUE                                            
!             20 CONTINUE                                               
!                                                                       
!     LINPACK.  THIS VERSION DATED 08/14/78 .                           
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DDOT                                                         
!     FORTRAN DSQRT                                                     
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION DDOT,T                                           
      DOUBLE PRECISION S                                                
      INTEGER J,JJ,JM1,K,KJ,KK                                          
!     BEGIN BLOCK WITH ...EXITS TO 40                                   
!                                                                       
!                                                                       
         JJ = 0                                                         
         DO 30 J = 1, N                                                 
            INFO = J                                                    
            S = 0.0D0                                                   
            JM1 = J - 1                                                 
            KJ = JJ                                                     
            KK = 0                                                      
            IF (JM1 .LT. 1) GO TO 20                                    
            DO 10 K = 1, JM1                                            
               KJ = KJ + 1                                              
               T = AP(KJ) - DDOT(K-1,AP(KK+1),1,AP(JJ+1),1)             
               KK = KK + K                                              
               T = T/AP(KK)                                             
               AP(KJ) = T                                               
               S = S + T*T                                              
   10       CONTINUE                                                    
   20       CONTINUE                                                    
            JJ = JJ + J                                                 
            S = AP(JJ) - S                                              
!     ......EXIT                                                        
            IF (S .LE. 0.0D0) GO TO 40                                  
            AP(JJ) = DSQRT(S)                                           
   30    CONTINUE                                                       
         INFO = 0                                                       
   40 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DPPSL(AP,N,B)                                          
      INTEGER N                                                         
      DOUBLE PRECISION AP(1),B(1)                                       
!                                                                       
!     DPPSL SOLVES THE DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE     
!     SYSTEM A * X = B                                                  
!     USING THE FACTORS COMPUTED BY DPPCO OR DPPFA.                     
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        AP      DOUBLE PRECISION (N*(N+1)/2)                           
!                THE OUTPUT FROM DPPCO OR DPPFA.                        
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!        B       DOUBLE PRECISION(N)                                    
!                THE RIGHT HAND SIDE VECTOR.                            
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        B       THE SOLUTION VECTOR  X .                               
!                                                                       
!     ERROR CONDITION                                                   
!                                                                       
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS     
!        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES            
!        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE    
!        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED    
!        CORRECTLY AND  INFO .EQ. 0 .                                   
!                                                                       
!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 
!     WITH  P  COLUMNS                                                  
!           CALL DPPCO(AP,N,RCOND,Z,INFO)                               
!           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ...          
!           DO 10 J = 1, P                                              
!              CALL DPPSL(AP,N,C(1,J))                                  
!        10 CONTINUE                                                    
!                                                                       
!     LINPACK.  THIS VERSION DATED 08/14/78 .                           
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DAXPY,DDOT                                                   
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION DDOT,T                                           
      INTEGER K,KB,KK                                                   
!                                                                       
      KK = 0                                                            
      DO 10 K = 1, N                                                    
         T = DDOT(K-1,AP(KK+1),1,B(1),1)                                
         KK = KK + K                                                    
         B(K) = (B(K) - T)/AP(KK)                                       
   10 CONTINUE                                                          
      DO 20 KB = 1, N                                                   
         K = N + 1 - KB                                                 
         B(K) = B(K)/AP(KK)                                             
         KK = KK - K                                                    
         T = -B(K)                                                      
         CALL DAXPY(K-1,T,AP(KK+1),1,B(1),1)                            
   20 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DPPDI(AP,N,DET,JOB)                                    
      INTEGER N,JOB                                                     
      DOUBLE PRECISION AP(1)                                            
      DOUBLE PRECISION DET(2)                                           
!                                                                       
!     DPPDI COMPUTES THE DETERMINANT AND INVERSE                        
!     OF A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE MATRIX          
!     USING THE FACTORS COMPUTED BY DPPCO OR DPPFA .                    
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        AP      DOUBLE PRECISION (N*(N+1)/2)                           
!                THE OUTPUT FROM DPPCO OR DPPFA.                        
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!        JOB     INTEGER                                                
!                = 11   BOTH DETERMINANT AND INVERSE.                   
!                = 01   INVERSE ONLY.                                   
!                = 10   DETERMINANT ONLY.                               
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        AP      THE UPPER TRIANGULAR HALF OF THE INVERSE .             
!                THE STRICT LOWER TRIANGLE IS UNALTERED.                
!                                                                       
!        DET     DOUBLE PRECISION(2)                                    
!                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.           
!                OTHERWISE NOT REFERENCED.                              
!                DETERMINANT = DET(1) * 10.0**DET(2)                    
!                WITH  1.0 .LE. DET(1) .LT. 10.0                        
!                OR  DET(1) .EQ. 0.0 .                                  
!                                                                       
!     ERROR CONDITION                                                   
!                                                                       
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS     
!        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.           
!        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY      
!        AND IF DPOCO OR DPOFA HAS SET INFO .EQ. 0 .                    
!                                                                       
!     LINPACK.  THIS VERSION DATED 08/14/78 .                           
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DAXPY,DSCAL                                                  
!     FORTRAN MOD                                                       
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION T                                                
      DOUBLE PRECISION S                                                
      INTEGER I,II,J,JJ,JM1,J1,K,KJ,KK,KP1,K1                           
!                                                                       
!     COMPUTE DETERMINANT                                               
!                                                                       
      IF (JOB/10 .EQ. 0) GO TO 70                                       
         DET(1) = 1.0D0                                                 
         DET(2) = 0.0D0                                                 
         S = 10.0D0                                                     
         II = 0                                                         
         DO 50 I = 1, N                                                 
            II = II + I                                                 
            DET(1) = AP(II)**2*DET(1)                                   
!        ...EXIT                                                        
            IF (DET(1) .EQ. 0.0D0) GO TO 60                             
   10       IF (DET(1) .GE. 1.0D0) GO TO 20                             
               DET(1) = S*DET(1)                                        
               DET(2) = DET(2) - 1.0D0                                  
            GO TO 10                                                    
   20       CONTINUE                                                    
   30       IF (DET(1) .LT. S) GO TO 40                                 
               DET(1) = DET(1)/S                                        
               DET(2) = DET(2) + 1.0D0                                  
            GO TO 30                                                    
   40       CONTINUE                                                    
   50    CONTINUE                                                       
   60    CONTINUE                                                       
   70 CONTINUE                                                          
!                                                                       
!     COMPUTE INVERSE(R)                                                
!                                                                       
      IF (MOD(JOB,10) .EQ. 0) GO TO 140                                 
         KK = 0                                                         
         DO 100 K = 1, N                                                
            K1 = KK + 1                                                 
            KK = KK + K                                                 
            AP(KK) = 1.0D0/AP(KK)                                       
            T = -AP(KK)                                                 
            CALL DSCAL(K-1,T,AP(K1),1)                                  
            KP1 = K + 1                                                 
            J1 = KK + 1                                                 
            KJ = KK + K                                                 
            IF (N .LT. KP1) GO TO 90                                    
            DO 80 J = KP1, N                                            
               T = AP(KJ)                                               
               AP(KJ) = 0.0D0                                           
               CALL DAXPY(K,T,AP(K1),1,AP(J1),1)                        
               J1 = J1 + J                                              
               KJ = KJ + J                                              
   80       CONTINUE                                                    
   90       CONTINUE                                                    
  100    CONTINUE                                                       
!                                                                       
!        FORM  INVERSE(R) * TRANS(INVERSE(R))                           
!                                                                       
         JJ = 0                                                         
         DO 130 J = 1, N                                                
            J1 = JJ + 1                                                 
            JJ = JJ + J                                                 
            JM1 = J - 1                                                 
            K1 = 1                                                      
            KJ = J1                                                     
            IF (JM1 .LT. 1) GO TO 120                                   
            DO 110 K = 1, JM1                                           
               T = AP(KJ)                                               
               CALL DAXPY(K,T,AP(J1),1,AP(K1),1)                        
               K1 = K1 + K                                              
               KJ = KJ + 1                                              
  110       CONTINUE                                                    
  120       CONTINUE                                                    
            T = AP(JJ)                                                  
            CALL DSCAL(J,T,AP(J1),1)                                    
  130    CONTINUE                                                       
  140 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DPBCO(ABD,LDA,N,M,RCOND,Z,INFO)                        
      INTEGER LDA,N,M,INFO                                              
      DOUBLE PRECISION ABD(LDA,1),Z(1)                                  
      DOUBLE PRECISION RCOND                                            
!                                                                       
!     DPBCO FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE      
!     MATRIX STORED IN BAND FORM AND ESTIMATES THE CONDITION OF THE     
!     MATRIX.                                                           
!                                                                       
!     IF  RCOND  IS NOT NEEDED, DPBFA IS SLIGHTLY FASTER.               
!     TO SOLVE  A*X = B , FOLLOW DPBCO BY DPBSL.                        
!     TO COMPUTE  INVERSE(A)*C , FOLLOW DPBCO BY DPBSL.                 
!     TO COMPUTE  DETERMINANT(A) , FOLLOW DPBCO BY DPBDI.               
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        ABD     DOUBLE PRECISION(LDA, N)                               
!                THE MATRIX TO BE FACTORED.  THE COLUMNS OF THE UPPER   
!                TRIANGLE ARE STORED IN THE COLUMNS OF ABD AND THE      
!                DIAGONALS OF THE UPPER TRIANGLE ARE STORED IN THE      
!                ROWS OF ABD .  SEE THE COMMENTS BELOW FOR DETAILS.     
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  ABD .              
!                LDA MUST BE .GE. M + 1 .                               
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!        M       INTEGER                                                
!                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.       
!                0 .LE. M .LT. N .                                      
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        ABD     AN UPPER TRIANGULAR MATRIX  R , STORED IN BAND         
!                FORM, SO THAT  A = TRANS(R)*R .                        
!                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.   
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
!                UNDERFLOWS.  IF INFO .NE. 0 , RCOND IS UNCHANGED.      
!                                                                       
!        Z       DOUBLE PRECISION(N)                                    
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.  
!                IF  A  IS SINGULAR TO WORKING PRECISION, THEN  Z  IS   
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT           
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .                    
!                IF  INFO .NE. 0 , Z  IS UNCHANGED.                     
!                                                                       
!        INFO    INTEGER                                                
!                = 0  FOR NORMAL RETURN.                                
!                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR    
!                     OF ORDER  K  IS NOT POSITIVE DEFINITE.            
!                                                                       
!     BAND STORAGE                                                      
!                                                                       
!           IF  A  IS A SYMMETRIC POSITIVE DEFINITE BAND MATRIX,        
!           THE FOLLOWING PROGRAM SEGMENT WILL SET UP THE INPUT.        
!                                                                       
!                   M = (BAND WIDTH ABOVE DIAGONAL)                     
!                   DO 20 J = 1, N                                      
!                      I1 = MAX0(1, J-M)                                
!                      DO 10 I = I1, J                                  
!                         K = I-J+M+1                                   
!                         ABD(K,J) = A(I,J)                             
!                10    CONTINUE                                         
!                20 CONTINUE                                            
!                                                                       
!           THIS USES  M + 1  ROWS OF  A , EXCEPT FOR THE  M BY M       
!           UPPER LEFT TRIANGLE, WHICH IS IGNORED.                      
!                                                                       
!     EXAMPLE..  IF THE ORIGINAL MATRIX IS                              
!                                                                       
!           11 12 13  0  0  0                                           
!           12 22 23 24  0  0                                           
!           13 23 33 34 35  0                                           
!            0 24 34 44 45 46                                           
!            0  0 35 45 55 56                                           
!            0  0  0 46 56 66                                           
!                                                                       
!     THEN  N = 6 , M = 2  AND  ABD  SHOULD CONTAIN                     
!                                                                       
!            *  * 13 24 35 46                                           
!            * 12 23 34 45 56                                           
!           11 22 33 44 55 66                                           
!                                                                       
!     LINPACK.  THIS VERSION DATED 08/14/78 .                           
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     LINPACK DPBFA                                                     
!     BLAS DAXPY,DDOT,DSCAL,DASUM                                       
!     FORTRAN DABS,DMAX1,MAX0,MIN0,DREAL,DSIGN                          
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION DDOT,EK,T,WK,WKM                                 
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM                           
      INTEGER I,J,J2,K,KB,KP1,L,LA,LB,LM,MU                             
!                                                                       
!                                                                       
!     FIND NORM OF A                                                    
!                                                                       
      DO 30 J = 1, N                                                    
         L = MIN0(J,M+1)                                                
         MU = MAX0(M+2-J,1)                                             
         Z(J) = DASUM(L,ABD(MU,J),1)                                    
         K = J - L                                                      
         IF (M .LT. MU) GO TO 20                                        
         DO 10 I = MU, M                                                
            K = K + 1                                                   
            Z(K) = Z(K) + DABS(ABD(I,J))                                
   10    CONTINUE                                                       
   20    CONTINUE                                                       
   30 CONTINUE                                                          
      ANORM = 0.0D0                                                     
      DO 40 J = 1, N                                                    
         ANORM = DMAX1(ANORM,Z(J))                                      
   40 CONTINUE                                                          
!                                                                       
!     FACTOR                                                            
!                                                                       
      CALL DPBFA(ABD,LDA,N,M,INFO)                                      
      IF (INFO .NE. 0) GO TO 180                                        
!                                                                       
!        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .           
!        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .      
!        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL        
!        GROWTH IN THE ELEMENTS OF W  WHERE  TRANS(R)*W = E .           
!        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.         
!                                                                       
!        SOLVE TRANS(R)*W = E                                           
!                                                                       
         EK = 1.0D0                                                     
         DO 50 J = 1, N                                                 
            Z(J) = 0.0D0                                                
   50    CONTINUE                                                       
         DO 110 K = 1, N                                                
            IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))                   
            IF (DABS(EK-Z(K)) .LE. ABD(M+1,K)) GO TO 60                 
               S = ABD(M+1,K)/DABS(EK-Z(K))                             
               CALL DSCAL(N,S,Z,1)                                      
               EK = S*EK                                                
   60       CONTINUE                                                    
            WK = EK - Z(K)                                              
            WKM = -EK - Z(K)                                            
            S = DABS(WK)                                                
            SM = DABS(WKM)                                              
            WK = WK/ABD(M+1,K)                                          
            WKM = WKM/ABD(M+1,K)                                        
            KP1 = K + 1                                                 
            J2 = MIN0(K+M,N)                                            
            I = M + 1                                                   
            IF (KP1 .GT. J2) GO TO 100                                  
               DO 70 J = KP1, J2                                        
                  I = I - 1                                             
                  SM = SM + DABS(Z(J)+WKM*ABD(I,J))                     
                  Z(J) = Z(J) + WK*ABD(I,J)                             
                  S = S + DABS(Z(J))                                    
   70          CONTINUE                                                 
               IF (S .GE. SM) GO TO 90                                  
                  T = WKM - WK                                          
                  WK = WKM                                              
                  I = M + 1                                             
                  DO 80 J = KP1, J2                                     
                     I = I - 1                                          
                     Z(J) = Z(J) + T*ABD(I,J)                           
   80             CONTINUE                                              
   90          CONTINUE                                                 
  100       CONTINUE                                                    
            Z(K) = WK                                                   
  110    CONTINUE                                                       
         S = 1.0D0/DASUM(N,Z,1)                                         
         CALL DSCAL(N,S,Z,1)                                            
!                                                                       
!        SOLVE  R*Y = W                                                 
!                                                                       
         DO 130 KB = 1, N                                               
            K = N + 1 - KB                                              
            IF (DABS(Z(K)) .LE. ABD(M+1,K)) GO TO 120                   
               S = ABD(M+1,K)/DABS(Z(K))                                
               CALL DSCAL(N,S,Z,1)                                      
  120       CONTINUE                                                    
            Z(K) = Z(K)/ABD(M+1,K)                                      
            LM = MIN0(K-1,M)                                            
            LA = M + 1 - LM                                             
            LB = K - LM                                                 
            T = -Z(K)                                                   
            CALL DAXPY(LM,T,ABD(LA,K),1,Z(LB),1)                        
  130    CONTINUE                                                       
         S = 1.0D0/DASUM(N,Z,1)                                         
         CALL DSCAL(N,S,Z,1)                                            
!                                                                       
         YNORM = 1.0D0                                                  
!                                                                       
!        SOLVE TRANS(R)*V = Y                                           
!                                                                       
         DO 150 K = 1, N                                                
            LM = MIN0(K-1,M)                                            
            LA = M + 1 - LM                                             
            LB = K - LM                                                 
            Z(K) = Z(K) - DDOT(LM,ABD(LA,K),1,Z(LB),1)                  
            IF (DABS(Z(K)) .LE. ABD(M+1,K)) GO TO 140                   
               S = ABD(M+1,K)/DABS(Z(K))                                
               CALL DSCAL(N,S,Z,1)                                      
               YNORM = S*YNORM                                          
  140       CONTINUE                                                    
            Z(K) = Z(K)/ABD(M+1,K)                                      
  150    CONTINUE                                                       
         S = 1.0D0/DASUM(N,Z,1)                                         
         CALL DSCAL(N,S,Z,1)                                            
         YNORM = S*YNORM                                                
!                                                                       
!        SOLVE  R*Z = W                                                 
!                                                                       
         DO 170 KB = 1, N                                               
            K = N + 1 - KB                                              
            IF (DABS(Z(K)) .LE. ABD(M+1,K)) GO TO 160                   
               S = ABD(M+1,K)/DABS(Z(K))                                
               CALL DSCAL(N,S,Z,1)                                      
               YNORM = S*YNORM                                          
  160       CONTINUE                                                    
            Z(K) = Z(K)/ABD(M+1,K)                                      
            LM = MIN0(K-1,M)                                            
            LA = M + 1 - LM                                             
            LB = K - LM                                                 
            T = -Z(K)                                                   
            CALL DAXPY(LM,T,ABD(LA,K),1,Z(LB),1)                        
  170    CONTINUE                                                       
!        MAKE ZNORM = 1.0                                               
         S = 1.0D0/DASUM(N,Z,1)                                         
         CALL DSCAL(N,S,Z,1)                                            
         YNORM = S*YNORM                                                
!                                                                       
         IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM                      
         IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0                            
  180 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DPBFA(ABD,LDA,N,M,INFO)                                
      INTEGER LDA,N,M,INFO                                              
      DOUBLE PRECISION ABD(LDA,1)                                       
!                                                                       
!     DPBFA FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE      
!     MATRIX STORED IN BAND FORM.                                       
!                                                                       
!     DPBFA IS USUALLY CALLED BY DPBCO, BUT IT CAN BE CALLED            
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.          
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        ABD     DOUBLE PRECISION(LDA, N)                               
!                THE MATRIX TO BE FACTORED.  THE COLUMNS OF THE UPPER   
!                TRIANGLE ARE STORED IN THE COLUMNS OF ABD AND THE      
!                DIAGONALS OF THE UPPER TRIANGLE ARE STORED IN THE      
!                ROWS OF ABD .  SEE THE COMMENTS BELOW FOR DETAILS.     
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  ABD .              
!                LDA MUST BE .GE. M + 1 .                               
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!        M       INTEGER                                                
!                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.       
!                0 .LE. M .LT. N .                                      
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        ABD     AN UPPER TRIANGULAR MATRIX  R , STORED IN BAND         
!                FORM, SO THAT  A = TRANS(R)*R .                        
!                                                                       
!        INFO    INTEGER                                                
!                = 0  FOR NORMAL RETURN.                                
!                = K  IF THE LEADING MINOR OF ORDER  K  IS NOT          
!                     POSITIVE DEFINITE.                                
!                                                                       
!     BAND STORAGE                                                      
!                                                                       
!           IF  A  IS A SYMMETRIC POSITIVE DEFINITE BAND MATRIX,        
!           THE FOLLOWING PROGRAM SEGMENT WILL SET UP THE INPUT.        
!                                                                       
!                   M = (BAND WIDTH ABOVE DIAGONAL)                     
!                   DO 20 J = 1, N                                      
!                      I1 = MAX0(1, J-M)                                
!                      DO 10 I = I1, J                                  
!                         K = I-J+M+1                                   
!                         ABD(K,J) = A(I,J)                             
!                10    CONTINUE                                         
!                20 CONTINUE                                            
!                                                                       
!     LINPACK.  THIS VERSION DATED 08/14/78 .                           
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DDOT                                                         
!     FORTRAN MAX0,DSQRT                                                
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION DDOT,T                                           
      DOUBLE PRECISION S                                                
      INTEGER IK,J,JK,K,MU                                              
!     BEGIN BLOCK WITH ...EXITS TO 40                                   
!                                                                       
!                                                                       
         DO 30 J = 1, N                                                 
            INFO = J                                                    
            S = 0.0D0                                                   
            IK = M + 1                                                  
            JK = MAX0(J-M,1)                                            
            MU = MAX0(M+2-J,1)                                          
            IF (M .LT. MU) GO TO 20                                     
            DO 10 K = MU, M                                             
               T = ABD(K,J) - DDOT(K-MU,ABD(IK,JK),1,ABD(MU,J),1)       
               T = T/ABD(M+1,JK)                                        
               ABD(K,J) = T                                             
               S = S + T*T                                              
               IK = IK - 1                                              
               JK = JK + 1                                              
   10       CONTINUE                                                    
   20       CONTINUE                                                    
            S = ABD(M+1,J) - S                                          
!     ......EXIT                                                        
            IF (S .LE. 0.0D0) GO TO 40                                  
            ABD(M+1,J) = DSQRT(S)                                       
   30    CONTINUE                                                       
         INFO = 0                                                       
   40 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DPBSL(ABD,LDA,N,M,B)                                   
      INTEGER LDA,N,M                                                   
      DOUBLE PRECISION ABD(LDA,1),B(1)                                  
!                                                                       
!     DPBSL SOLVES THE DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE     
!     BAND SYSTEM  A*X = B                                              
!     USING THE FACTORS COMPUTED BY DPBCO OR DPBFA.                     
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        ABD     DOUBLE PRECISION(LDA, N)                               
!                THE OUTPUT FROM DPBCO OR DPBFA.                        
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  ABD .              
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!        M       INTEGER                                                
!                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.       
!                                                                       
!        B       DOUBLE PRECISION(N)                                    
!                THE RIGHT HAND SIDE VECTOR.                            
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        B       THE SOLUTION VECTOR  X .                               
!                                                                       
!     ERROR CONDITION                                                   
!                                                                       
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS     
!        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES            
!        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE    
!        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED    
!        CORRECTLY AND  INFO .EQ. 0 .                                   
!                                                                       
!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 
!     WITH  P  COLUMNS                                                  
!           CALL DPBCO(ABD,LDA,N,RCOND,Z,INFO)                          
!           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ...          
!           DO 10 J = 1, P                                              
!              CALL DPBSL(ABD,LDA,N,C(1,J))                             
!        10 CONTINUE                                                    
!                                                                       
!     LINPACK.  THIS VERSION DATED 08/14/78 .                           
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DAXPY,DDOT                                                   
!     FORTRAN MIN0                                                      
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION DDOT,T                                           
      INTEGER K,KB,LA,LB,LM                                             
!                                                                       
!     SOLVE TRANS(R)*Y = B                                              
!                                                                       
      DO 10 K = 1, N                                                    
         LM = MIN0(K-1,M)                                               
         LA = M + 1 - LM                                                
         LB = K - LM                                                    
         T = DDOT(LM,ABD(LA,K),1,B(LB),1)                               
         B(K) = (B(K) - T)/ABD(M+1,K)                                   
   10 CONTINUE                                                          
!                                                                       
!     SOLVE R*X = Y                                                     
!                                                                       
      DO 20 KB = 1, N                                                   
         K = N + 1 - KB                                                 
         LM = MIN0(K-1,M)                                               
         LA = M + 1 - LM                                                
         LB = K - LM                                                    
         B(K) = B(K)/ABD(M+1,K)                                         
         T = -B(K)                                                      
         CALL DAXPY(LM,T,ABD(LA,K),1,B(LB),1)                           
   20 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DPBDI(ABD,LDA,N,M,DET)                                 
      INTEGER LDA,N,M                                                   
      DOUBLE PRECISION ABD(LDA,1)                                       
      DOUBLE PRECISION DET(2)                                           
!                                                                       
!     DPBDI COMPUTES THE DETERMINANT                                    
!     OF A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE BAND MATRIX     
!     USING THE FACTORS COMPUTED BY DPBCO OR DPBFA.                     
!     IF THE INVERSE IS NEEDED, USE DPBSL  N  TIMES.                    
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        ABD     DOUBLE PRECISION(LDA, N)                               
!                THE OUTPUT FROM DPBCO OR DPBFA.                        
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  ABD .              
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!        M       INTEGER                                                
!                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.       
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        DET     DOUBLE PRECISION(2)                                    
!                DETERMINANT OF ORIGINAL MATRIX IN THE FORM             
!                DETERMINANT = DET(1) * 10.0**DET(2)                    
!                WITH  1.0 .LE. DET(1) .LT. 10.0                        
!                OR  DET(1) .EQ. 0.0 .                                  
!                                                                       
!     LINPACK.  THIS VERSION DATED 08/14/78 .                           
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION S                                                
      INTEGER I                                                         
!                                                                       
!     COMPUTE DETERMINANT                                               
!                                                                       
      DET(1) = 1.0D0                                                    
      DET(2) = 0.0D0                                                    
      S = 10.0D0                                                        
      DO 50 I = 1, N                                                    
         DET(1) = ABD(M+1,I)**2*DET(1)                                  
!     ...EXIT                                                           
         IF (DET(1) .EQ. 0.0D0) GO TO 60                                
   10    IF (DET(1) .GE. 1.0D0) GO TO 20                                
            DET(1) = S*DET(1)                                           
            DET(2) = DET(2) - 1.0D0                                     
         GO TO 10                                                       
   20    CONTINUE                                                       
   30    IF (DET(1) .LT. S) GO TO 40                                    
            DET(1) = DET(1)/S                                           
            DET(2) = DET(2) + 1.0D0                                     
         GO TO 30                                                       
   40    CONTINUE                                                       
   50 CONTINUE                                                          
   60 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DSICO(A,LDA,N,KPVT,RCOND,Z)                            
      INTEGER LDA,N,KPVT(1)                                             
      DOUBLE PRECISION A(LDA,1),Z(1)                                    
      DOUBLE PRECISION RCOND                                            
!                                                                       
!     DSICO FACTORS A DOUBLE PRECISION SYMMETRIC MATRIX BY ELIMINATION  
!     WITH SYMMETRIC PIVOTING AND ESTIMATES THE CONDITION OF THE        
!     MATRIX.                                                           
!                                                                       
!     IF  RCOND  IS NOT NEEDED, DSIFA IS SLIGHTLY FASTER.               
!     TO SOLVE  A*X = B , FOLLOW DSICO BY DSISL.                        
!     TO COMPUTE  INVERSE(A)*C , FOLLOW DSICO BY DSISL.                 
!     TO COMPUTE  INVERSE(A) , FOLLOW DSICO BY DSIDI.                   
!     TO COMPUTE  DETERMINANT(A) , FOLLOW DSICO BY DSIDI.               
!     TO COMPUTE  INERTIA(A), FOLLOW DSICO BY DSIDI.                    
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        A       DOUBLE PRECISION(LDA, N)                               
!                THE SYMMETRIC MATRIX TO BE FACTORED.                   
!                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED.         
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  A .                
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!     OUTPUT                                                            
!                                                                       
!        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH      
!                WERE USED TO OBTAIN IT.                                
!                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)     
!                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT         
!                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE            
!                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL            
!                WITH 1 BY 1 AND 2 BY 2 BLOCKS.                         
!                                                                       
!        KPVT    INTEGER(N)                                             
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
!     LINPACK DSIFA                                                     
!     BLAS DAXPY,DDOT,DSCAL,DASUM                                       
!     FORTRAN DABS,DMAX1,IABS,DSIGN                                     
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DDOT,DENOM,EK,T                  
      DOUBLE PRECISION ANORM,S,DASUM,YNORM                              
      INTEGER I,INFO,J,JM1,K,KP,KPS,KS                                  
!                                                                       
!                                                                       
!     FIND NORM OF A USING ONLY UPPER HALF                              
!                                                                       
      DO 30 J = 1, N                                                    
         Z(J) = DASUM(J,A(1,J),1)                                       
         JM1 = J - 1                                                    
         IF (JM1 .LT. 1) GO TO 20                                       
         DO 10 I = 1, JM1                                               
            Z(I) = Z(I) + DABS(A(I,J))                                  
   10    CONTINUE                                                       
   20    CONTINUE                                                       
   30 CONTINUE                                                          
      ANORM = 0.0D0                                                     
      DO 40 J = 1, N                                                    
         ANORM = DMAX1(ANORM,Z(J))                                      
   40 CONTINUE                                                          
!                                                                       
!     FACTOR                                                            
!                                                                       
      CALL DSIFA(A,LDA,N,KPVT,INFO)                                     
!                                                                       
!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .              
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .         
!     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL           
!     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E .                   
!     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.            
!                                                                       
!     SOLVE U*D*W = E                                                   
!                                                                       
      EK = 1.0D0                                                        
      DO 50 J = 1, N                                                    
         Z(J) = 0.0D0                                                   
   50 CONTINUE                                                          
      K = N                                                             
   60 IF (K .EQ. 0) GO TO 120                                           
         KS = 1                                                         
         IF (KPVT(K) .LT. 0) KS = 2                                     
         KP = IABS(KPVT(K))                                             
         KPS = K + 1 - KS                                               
         IF (KP .EQ. KPS) GO TO 70                                      
            T = Z(KPS)                                                  
            Z(KPS) = Z(KP)                                              
            Z(KP) = T                                                   
   70    CONTINUE                                                       
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,Z(K))                       
         Z(K) = Z(K) + EK                                               
         CALL DAXPY(K-KS,Z(K),A(1,K),1,Z(1),1)                          
         IF (KS .EQ. 1) GO TO 80                                        
            IF (Z(K-1) .NE. 0.0D0) EK = DSIGN(EK,Z(K-1))                
            Z(K-1) = Z(K-1) + EK                                        
            CALL DAXPY(K-KS,Z(K-1),A(1,K-1),1,Z(1),1)                   
   80    CONTINUE                                                       
         IF (KS .EQ. 2) GO TO 100                                       
            IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 90                  
               S = DABS(A(K,K))/DABS(Z(K))                              
               CALL DSCAL(N,S,Z,1)                                      
               EK = S*EK                                                
   90       CONTINUE                                                    
            IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)                   
            IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0                         
         GO TO 110                                                      
  100    CONTINUE                                                       
            AK = A(K,K)/A(K-1,K)                                        
            AKM1 = A(K-1,K-1)/A(K-1,K)                                  
            BK = Z(K)/A(K-1,K)                                          
            BKM1 = Z(K-1)/A(K-1,K)                                      
            DENOM = AK*AKM1 - 1.0D0                                     
            Z(K) = (AKM1*BK - BKM1)/DENOM                               
            Z(K-1) = (AK*BKM1 - BK)/DENOM                               
  110    CONTINUE                                                       
         K = K - KS                                                     
      GO TO 60                                                          
  120 CONTINUE                                                          
      S = 1.0D0/DASUM(N,Z,1)                                            
      CALL DSCAL(N,S,Z,1)                                               
!                                                                       
!     SOLVE TRANS(U)*Y = W                                              
!                                                                       
      K = 1                                                             
  130 IF (K .GT. N) GO TO 160                                           
         KS = 1                                                         
         IF (KPVT(K) .LT. 0) KS = 2                                     
         IF (K .EQ. 1) GO TO 150                                        
            Z(K) = Z(K) + DDOT(K-1,A(1,K),1,Z(1),1)                     
            IF (KS .EQ. 2)                                               &
               Z(K+1) = Z(K+1) + DDOT(K-1,A(1,K+1),1,Z(1),1)            
            KP = IABS(KPVT(K))                                          
            IF (KP .EQ. K) GO TO 140                                    
               T = Z(K)                                                 
               Z(K) = Z(KP)                                             
               Z(KP) = T                                                
  140       CONTINUE                                                    
  150    CONTINUE                                                       
         K = K + KS                                                     
      GO TO 130                                                         
  160 CONTINUE                                                          
      S = 1.0D0/DASUM(N,Z,1)                                            
      CALL DSCAL(N,S,Z,1)                                               
!                                                                       
      YNORM = 1.0D0                                                     
!                                                                       
!     SOLVE U*D*V = Y                                                   
!                                                                       
      K = N                                                             
  170 IF (K .EQ. 0) GO TO 230                                           
         KS = 1                                                         
         IF (KPVT(K) .LT. 0) KS = 2                                     
         IF (K .EQ. KS) GO TO 190                                       
            KP = IABS(KPVT(K))                                          
            KPS = K + 1 - KS                                            
            IF (KP .EQ. KPS) GO TO 180                                  
               T = Z(KPS)                                               
               Z(KPS) = Z(KP)                                           
               Z(KP) = T                                                
  180       CONTINUE                                                    
            CALL DAXPY(K-KS,Z(K),A(1,K),1,Z(1),1)                       
            IF (KS .EQ. 2) CALL DAXPY(K-KS,Z(K-1),A(1,K-1),1,Z(1),1)    
  190    CONTINUE                                                       
         IF (KS .EQ. 2) GO TO 210                                       
            IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 200                 
               S = DABS(A(K,K))/DABS(Z(K))                              
               CALL DSCAL(N,S,Z,1)                                      
               YNORM = S*YNORM                                          
  200       CONTINUE                                                    
            IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)                   
            IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0                         
         GO TO 220                                                      
  210    CONTINUE                                                       
            AK = A(K,K)/A(K-1,K)                                        
            AKM1 = A(K-1,K-1)/A(K-1,K)                                  
            BK = Z(K)/A(K-1,K)                                          
            BKM1 = Z(K-1)/A(K-1,K)                                      
            DENOM = AK*AKM1 - 1.0D0                                     
            Z(K) = (AKM1*BK - BKM1)/DENOM                               
            Z(K-1) = (AK*BKM1 - BK)/DENOM                               
  220    CONTINUE                                                       
         K = K - KS                                                     
      GO TO 170                                                         
  230 CONTINUE                                                          
      S = 1.0D0/DASUM(N,Z,1)                                            
      CALL DSCAL(N,S,Z,1)                                               
      YNORM = S*YNORM                                                   
!                                                                       
!     SOLVE TRANS(U)*Z = V                                              
!                                                                       
      K = 1                                                             
  240 IF (K .GT. N) GO TO 270                                           
         KS = 1                                                         
         IF (KPVT(K) .LT. 0) KS = 2                                     
         IF (K .EQ. 1) GO TO 260                                        
            Z(K) = Z(K) + DDOT(K-1,A(1,K),1,Z(1),1)                     
            IF (KS .EQ. 2)                                               &
               Z(K+1) = Z(K+1) + DDOT(K-1,A(1,K+1),1,Z(1),1)            
            KP = IABS(KPVT(K))                                          
            IF (KP .EQ. K) GO TO 250                                    
               T = Z(K)                                                 
               Z(K) = Z(KP)                                             
               Z(KP) = T                                                
  250       CONTINUE                                                    
  260    CONTINUE                                                       
         K = K + KS                                                     
      GO TO 240                                                         
  270 CONTINUE                                                          
!     MAKE ZNORM = 1.0                                                  
      S = 1.0D0/DASUM(N,Z,1)                                            
      CALL DSCAL(N,S,Z,1)                                               
      YNORM = S*YNORM                                                   
!                                                                       
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM                         
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0                               
      RETURN                                                            
      END                                                               
      SUBROUTINE DSIFA(A,LDA,N,KPVT,INFO)                               
      INTEGER LDA,N,KPVT(1),INFO                                        
      DOUBLE PRECISION A(LDA,1)                                         
!                                                                       
!     DSIFA FACTORS A DOUBLE PRECISION SYMMETRIC MATRIX BY ELIMINATION  
!     WITH SYMMETRIC PIVOTING.                                          
!                                                                       
!     TO SOLVE  A*X = B , FOLLOW DSIFA BY DSISL.                        
!     TO COMPUTE  INVERSE(A)*C , FOLLOW DSIFA BY DSISL.                 
!     TO COMPUTE  DETERMINANT(A) , FOLLOW DSIFA BY DSIDI.               
!     TO COMPUTE  INERTIA(A) , FOLLOW DSIFA BY DSIDI.                   
!     TO COMPUTE  INVERSE(A) , FOLLOW DSIFA BY DSIDI.                   
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        A       DOUBLE PRECISION(LDA,N)                                
!                THE SYMMETRIC MATRIX TO BE FACTORED.                   
!                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED.         
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  A .                
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH      
!                WERE USED TO OBTAIN IT.                                
!                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)     
!                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT         
!                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE            
!                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL            
!                WITH 1 BY 1 AND 2 BY 2 BLOCKS.                         
!                                                                       
!        KPVT    INTEGER(N)                                             
!                AN INTEGER VECTOR OF PIVOT INDICES.                    
!                                                                       
!        INFO    INTEGER                                                
!                = 0  NORMAL VALUE.                                     
!                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS      
!                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE,       
!                     BUT IT DOES INDICATE THAT DSISL OR DSIDI MAY      
!                     DIVIDE BY ZERO IF CALLED.                         
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.            
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DAXPY,DSWAP,IDAMAX                                           
!     FORTRAN DABS,DMAX1,DSQRT                                          
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T              
      DOUBLE PRECISION ABSAKK,ALPHA,COLMAX,ROWMAX                       
      INTEGER IMAX,IMAXP1,J,JJ,JMAX,K,KM1,KM2,KSTEP,IDAMAX              
      LOGICAL SWAP                                                      
!                                                                       
!                                                                       
!     INITIALIZE                                                        
!                                                                       
!     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.                       
      ALPHA = (1.0D0 + DSQRT(17.0D0))/8.0D0                             
!                                                                       
      INFO = 0                                                          
!                                                                       
!     MAIN LOOP ON K, WHICH GOES FROM N TO 1.                           
!                                                                       
      K = N                                                             
   10 CONTINUE                                                          
!                                                                       
!        LEAVE THE LOOP IF K=0 OR K=1.                                  
!                                                                       
!     ...EXIT                                                           
         IF (K .EQ. 0) GO TO 200                                        
         IF (K .GT. 1) GO TO 20                                         
            KPVT(1) = 1                                                 
            IF (A(1,1) .EQ. 0.0D0) INFO = 1                             
!     ......EXIT                                                        
            GO TO 200                                                   
   20    CONTINUE                                                       
!                                                                       
!        THIS SECTION OF CODE DETERMINES THE KIND OF                    
!        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,            
!        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND          
!        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS                
!        REQUIRED.                                                      
!                                                                       
         KM1 = K - 1                                                    
         ABSAKK = DABS(A(K,K))                                          
!                                                                       
!        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN                  
!        COLUMN K.                                                      
!                                                                       
         IMAX = IDAMAX(K-1,A(1,K),1)                                    
         COLMAX = DABS(A(IMAX,K))                                       
         IF (ABSAKK .LT. ALPHA*COLMAX) GO TO 30                         
            KSTEP = 1                                                   
            SWAP = .FALSE.                                              
         GO TO 90                                                       
   30    CONTINUE                                                       
!                                                                       
!           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN               
!           ROW IMAX.                                                   
!                                                                       
            ROWMAX = 0.0D0                                              
            IMAXP1 = IMAX + 1                                           
            DO 40 J = IMAXP1, K                                         
               ROWMAX = DMAX1(ROWMAX,DABS(A(IMAX,J)))                   
   40       CONTINUE                                                    
            IF (IMAX .EQ. 1) GO TO 50                                   
               JMAX = IDAMAX(IMAX-1,A(1,IMAX),1)                        
               ROWMAX = DMAX1(ROWMAX,DABS(A(JMAX,IMAX)))                
   50       CONTINUE                                                    
            IF (DABS(A(IMAX,IMAX)) .LT. ALPHA*ROWMAX) GO TO 60          
               KSTEP = 1                                                
               SWAP = .TRUE.                                            
            GO TO 80                                                    
   60       CONTINUE                                                    
            IF (ABSAKK .LT. ALPHA*COLMAX*(COLMAX/ROWMAX)) GO TO 70      
               KSTEP = 1                                                
               SWAP = .FALSE.                                           
            GO TO 80                                                    
   70       CONTINUE                                                    
               KSTEP = 2                                                
               SWAP = IMAX .NE. KM1                                     
   80       CONTINUE                                                    
   90    CONTINUE                                                       
         IF (DMAX1(ABSAKK,COLMAX) .NE. 0.0D0) GO TO 100                 
!                                                                       
!           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.           
!                                                                       
            KPVT(K) = K                                                 
            INFO = K                                                    
         GO TO 190                                                      
  100    CONTINUE                                                       
         IF (KSTEP .EQ. 2) GO TO 140                                    
!                                                                       
!           1 X 1 PIVOT BLOCK.                                          
!                                                                       
            IF (.NOT.SWAP) GO TO 120                                    
!                                                                       
!              PERFORM AN INTERCHANGE.                                  
!                                                                       
               CALL DSWAP(IMAX,A(1,IMAX),1,A(1,K),1)                    
               DO 110 JJ = IMAX, K                                      
                  J = K + IMAX - JJ                                     
                  T = A(J,K)                                            
                  A(J,K) = A(IMAX,J)                                    
                  A(IMAX,J) = T                                         
  110          CONTINUE                                                 
  120       CONTINUE                                                    
!                                                                       
!           PERFORM THE ELIMINATION.                                    
!                                                                       
            DO 130 JJ = 1, KM1                                          
               J = K - JJ                                               
               MULK = -A(J,K)/A(K,K)                                    
               T = MULK                                                 
               CALL DAXPY(J,T,A(1,K),1,A(1,J),1)                        
               A(J,K) = MULK                                            
  130       CONTINUE                                                    
!                                                                       
!           SET THE PIVOT ARRAY.                                        
!                                                                       
            KPVT(K) = K                                                 
            IF (SWAP) KPVT(K) = IMAX                                    
         GO TO 190                                                      
  140    CONTINUE                                                       
!                                                                       
!           2 X 2 PIVOT BLOCK.                                          
!                                                                       
            IF (.NOT.SWAP) GO TO 160                                    
!                                                                       
!              PERFORM AN INTERCHANGE.                                  
!                                                                       
               CALL DSWAP(IMAX,A(1,IMAX),1,A(1,K-1),1)                  
               DO 150 JJ = IMAX, KM1                                    
                  J = KM1 + IMAX - JJ                                   
                  T = A(J,K-1)                                          
                  A(J,K-1) = A(IMAX,J)                                  
                  A(IMAX,J) = T                                         
  150          CONTINUE                                                 
               T = A(K-1,K)                                             
               A(K-1,K) = A(IMAX,K)                                     
               A(IMAX,K) = T                                            
  160       CONTINUE                                                    
!                                                                       
!           PERFORM THE ELIMINATION.                                    
!                                                                       
            KM2 = K - 2                                                 
            IF (KM2 .EQ. 0) GO TO 180                                   
               AK = A(K,K)/A(K-1,K)                                     
               AKM1 = A(K-1,K-1)/A(K-1,K)                               
               DENOM = 1.0D0 - AK*AKM1                                  
               DO 170 JJ = 1, KM2                                       
                  J = KM1 - JJ                                          
                  BK = A(J,K)/A(K-1,K)                                  
                  BKM1 = A(J,K-1)/A(K-1,K)                              
                  MULK = (AKM1*BK - BKM1)/DENOM                         
                  MULKM1 = (AK*BKM1 - BK)/DENOM                         
                  T = MULK                                              
                  CALL DAXPY(J,T,A(1,K),1,A(1,J),1)                     
                  T = MULKM1                                            
                  CALL DAXPY(J,T,A(1,K-1),1,A(1,J),1)                   
                  A(J,K) = MULK                                         
                  A(J,K-1) = MULKM1                                     
  170          CONTINUE                                                 
  180       CONTINUE                                                    
!                                                                       
!           SET THE PIVOT ARRAY.                                        
!                                                                       
            KPVT(K) = 1 - K                                             
            IF (SWAP) KPVT(K) = -IMAX                                   
            KPVT(K-1) = KPVT(K)                                         
  190    CONTINUE                                                       
         K = K - KSTEP                                                  
      GO TO 10                                                          
  200 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DSISL(A,LDA,N,KPVT,B)                                  
      INTEGER LDA,N,KPVT(1)                                             
      DOUBLE PRECISION A(LDA,1),B(1)                                    
!                                                                       
!     DSISL SOLVES THE DOUBLE PRECISION SYMMETRIC SYSTEM                
!     A * X = B                                                         
!     USING THE FACTORS COMPUTED BY DSIFA.                              
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        A       DOUBLE PRECISION(LDA,N)                                
!                THE OUTPUT FROM DSIFA.                                 
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  A .                
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!        KPVT    INTEGER(N)                                             
!                THE PIVOT VECTOR FROM DSIFA.                           
!                                                                       
!        B       DOUBLE PRECISION(N)                                    
!                THE RIGHT HAND SIDE VECTOR.                            
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        B       THE SOLUTION VECTOR  X .                               
!                                                                       
!     ERROR CONDITION                                                   
!                                                                       
!        A DIVISION BY ZERO MAY OCCUR IF  DSICO  HAS SET RCOND .EQ. 0.0 
!        OR  DSIFA  HAS SET INFO .NE. 0  .                              
!                                                                       
!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 
!     WITH  P  COLUMNS                                                  
!           CALL DSIFA(A,LDA,N,KPVT,INFO)                               
!           IF (INFO .NE. 0) GO TO ...                                  
!           DO 10 J = 1, P                                              
!              CALL DSISL(A,LDA,N,KPVT,C(1,J))                          
!        10 CONTINUE                                                    
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.            
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DAXPY,DDOT                                                   
!     FORTRAN IABS                                                      
!                                                                       
!     INTERNAL VARIABLES.                                               
!                                                                       
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DDOT,DENOM,TEMP                  
      INTEGER K,KP                                                      
!                                                                       
!     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND                    
!     D INVERSE TO B.                                                   
!                                                                       
      K = N                                                             
   10 IF (K .EQ. 0) GO TO 80                                            
         IF (KPVT(K) .LT. 0) GO TO 40                                   
!                                                                       
!           1 X 1 PIVOT BLOCK.                                          
!                                                                       
            IF (K .EQ. 1) GO TO 30                                      
               KP = KPVT(K)                                             
               IF (KP .EQ. K) GO TO 20                                  
!                                                                       
!                 INTERCHANGE.                                          
!                                                                       
                  TEMP = B(K)                                           
                  B(K) = B(KP)                                          
                  B(KP) = TEMP                                          
   20          CONTINUE                                                 
!                                                                       
!              APPLY THE TRANSFORMATION.                                
!                                                                       
               CALL DAXPY(K-1,B(K),A(1,K),1,B(1),1)                     
   30       CONTINUE                                                    
!                                                                       
!           APPLY D INVERSE.                                            
!                                                                       
            B(K) = B(K)/A(K,K)                                          
            K = K - 1                                                   
         GO TO 70                                                       
   40    CONTINUE                                                       
!                                                                       
!           2 X 2 PIVOT BLOCK.                                          
!                                                                       
            IF (K .EQ. 2) GO TO 60                                      
               KP = IABS(KPVT(K))                                       
               IF (KP .EQ. K - 1) GO TO 50                              
!                                                                       
!                 INTERCHANGE.                                          
!                                                                       
                  TEMP = B(K-1)                                         
                  B(K-1) = B(KP)                                        
                  B(KP) = TEMP                                          
   50          CONTINUE                                                 
!                                                                       
!              APPLY THE TRANSFORMATION.                                
!                                                                       
               CALL DAXPY(K-2,B(K),A(1,K),1,B(1),1)                     
               CALL DAXPY(K-2,B(K-1),A(1,K-1),1,B(1),1)                 
   60       CONTINUE                                                    
!                                                                       
!           APPLY D INVERSE.                                            
!                                                                       
            AK = A(K,K)/A(K-1,K)                                        
            AKM1 = A(K-1,K-1)/A(K-1,K)                                  
            BK = B(K)/A(K-1,K)                                          
            BKM1 = B(K-1)/A(K-1,K)                                      
            DENOM = AK*AKM1 - 1.0D0                                     
            B(K) = (AKM1*BK - BKM1)/DENOM                               
            B(K-1) = (AK*BKM1 - BK)/DENOM                               
            K = K - 2                                                   
   70    CONTINUE                                                       
      GO TO 10                                                          
   80 CONTINUE                                                          
!                                                                       
!     LOOP FORWARD APPLYING THE TRANSFORMATIONS.                        
!                                                                       
      K = 1                                                             
   90 IF (K .GT. N) GO TO 160                                           
         IF (KPVT(K) .LT. 0) GO TO 120                                  
!                                                                       
!           1 X 1 PIVOT BLOCK.                                          
!                                                                       
            IF (K .EQ. 1) GO TO 110                                     
!                                                                       
!              APPLY THE TRANSFORMATION.                                
!                                                                       
               B(K) = B(K) + DDOT(K-1,A(1,K),1,B(1),1)                  
               KP = KPVT(K)                                             
               IF (KP .EQ. K) GO TO 100                                 
!                                                                       
!                 INTERCHANGE.                                          
!                                                                       
                  TEMP = B(K)                                           
                  B(K) = B(KP)                                          
                  B(KP) = TEMP                                          
  100          CONTINUE                                                 
  110       CONTINUE                                                    
            K = K + 1                                                   
         GO TO 150                                                      
  120    CONTINUE                                                       
!                                                                       
!           2 X 2 PIVOT BLOCK.                                          
!                                                                       
            IF (K .EQ. 1) GO TO 140                                     
!                                                                       
!              APPLY THE TRANSFORMATION.                                
!                                                                       
               B(K) = B(K) + DDOT(K-1,A(1,K),1,B(1),1)                  
               B(K+1) = B(K+1) + DDOT(K-1,A(1,K+1),1,B(1),1)            
               KP = IABS(KPVT(K))                                       
               IF (KP .EQ. K) GO TO 130                                 
!                                                                       
!                 INTERCHANGE.                                          
!                                                                       
                  TEMP = B(K)                                           
                  B(K) = B(KP)                                          
                  B(KP) = TEMP                                          
  130          CONTINUE                                                 
  140       CONTINUE                                                    
            K = K + 2                                                   
  150    CONTINUE                                                       
      GO TO 90                                                          
  160 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DSIDI(A,LDA,N,KPVT,DET,INERT,WORK,JOB)                 
      INTEGER LDA,N,JOB                                                 
      DOUBLE PRECISION A(LDA,1),WORK(1)                                 
      DOUBLE PRECISION DET(2)                                           
      INTEGER KPVT(1),INERT(3)                                          
!                                                                       
!     DSIDI COMPUTES THE DETERMINANT, INERTIA AND INVERSE               
!     OF A DOUBLE PRECISION SYMMETRIC MATRIX USING THE FACTORS FROM     
!     DSIFA.                                                            
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        A       DOUBLE PRECISION(LDA,N)                                
!                THE OUTPUT FROM DSIFA.                                 
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY A.                  
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX A.                             
!                                                                       
!        KPVT    INTEGER(N)                                             
!                THE PIVOT VECTOR FROM DSIFA.                           
!                                                                       
!        WORK    DOUBLE PRECISION(N)                                    
!                WORK VECTOR.  CONTENTS DESTROYED.                      
!                                                                       
!        JOB     INTEGER                                                
!                JOB HAS THE DECIMAL EXPANSION  ABC  WHERE              
!                   IF  C .NE. 0, THE INVERSE IS COMPUTED,              
!                   IF  B .NE. 0, THE DETERMINANT IS COMPUTED,          
!                   IF  A .NE. 0, THE INERTIA IS COMPUTED.              
!                                                                       
!                FOR EXAMPLE, JOB = 111  GIVES ALL THREE.               
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        VARIABLES NOT REQUESTED BY JOB ARE NOT USED.                   
!                                                                       
!        A      CONTAINS THE UPPER TRIANGLE OF THE INVERSE OF           
!               THE ORIGINAL MATRIX.  THE STRICT LOWER TRIANGLE         
!               IS NEVER REFERENCED.                                    
!                                                                       
!        DET    DOUBLE PRECISION(2)                                     
!               DETERMINANT OF ORIGINAL MATRIX.                         
!               DETERMINANT = DET(1) * 10.0**DET(2)                     
!               WITH 1.0 .LE. DABS(DET(1)) .LT. 10.0                    
!               OR DET(1) = 0.0.                                        
!                                                                       
!        INERT  INTEGER(3)                                              
!               THE INERTIA OF THE ORIGINAL MATRIX.                     
!               INERT(1)  =  NUMBER OF POSITIVE EIGENVALUES.            
!               INERT(2)  =  NUMBER OF NEGATIVE EIGENVALUES.            
!               INERT(3)  =  NUMBER OF ZERO EIGENVALUES.                
!                                                                       
!     ERROR CONDITION                                                   
!                                                                       
!        A DIVISION BY ZERO MAY OCCUR IF THE INVERSE IS REQUESTED       
!        AND  DSICO  HAS SET RCOND .EQ. 0.0                             
!        OR  DSIFA  HAS SET  INFO .NE. 0 .                              
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB             
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DAXPY,DCOPY,DDOT,DSWAP                                       
!     FORTRAN DABS,IABS,MOD                                             
!                                                                       
!     INTERNAL VARIABLES.                                               
!                                                                       
      DOUBLE PRECISION AKKP1,DDOT,TEMP                                  
      DOUBLE PRECISION TEN,D,T,AK,AKP1                                  
      INTEGER J,JB,K,KM1,KS,KSTEP                                       
      LOGICAL NOINV,NODET,NOERT                                         
!                                                                       
      NOINV = MOD(JOB,10) .EQ. 0                                        
      NODET = MOD(JOB,100)/10 .EQ. 0                                    
      NOERT = MOD(JOB,1000)/100 .EQ. 0                                  
!                                                                       
      IF (NODET .AND. NOERT) GO TO 140                                  
         IF (NOERT) GO TO 10                                            
            INERT(1) = 0                                                
            INERT(2) = 0                                                
            INERT(3) = 0                                                
   10    CONTINUE                                                       
         IF (NODET) GO TO 20                                            
            DET(1) = 1.0D0                                              
            DET(2) = 0.0D0                                              
            TEN = 10.0D0                                                
   20    CONTINUE                                                       
         T = 0.0D0                                                      
         DO 130 K = 1, N                                                
            D = A(K,K)                                                  
!                                                                       
!           CHECK IF 1 BY 1                                             
!                                                                       
            IF (KPVT(K) .GT. 0) GO TO 50                                
!                                                                       
!              2 BY 2 BLOCK                                             
!              USE DET (D  S)  =  (D/T * C - T) * T  ,  T = DABS(S)     
!                      (S  C)                                           
!              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.                    
!              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.      
!                                                                       
               IF (T .NE. 0.0D0) GO TO 30                               
                  T = DABS(A(K,K+1))                                    
                  D = (D/T)*A(K+1,K+1) - T                              
               GO TO 40                                                 
   30          CONTINUE                                                 
                  D = T                                                 
                  T = 0.0D0                                             
   40          CONTINUE                                                 
   50       CONTINUE                                                    
!                                                                       
            IF (NOERT) GO TO 60                                         
               IF (D .GT. 0.0D0) INERT(1) = INERT(1) + 1                
               IF (D .LT. 0.0D0) INERT(2) = INERT(2) + 1                
               IF (D .EQ. 0.0D0) INERT(3) = INERT(3) + 1                
   60       CONTINUE                                                    
!                                                                       
            IF (NODET) GO TO 120                                        
               DET(1) = D*DET(1)                                        
               IF (DET(1) .EQ. 0.0D0) GO TO 110                         
   70             IF (DABS(DET(1)) .GE. 1.0D0) GO TO 80                 
                     DET(1) = TEN*DET(1)                                
                     DET(2) = DET(2) - 1.0D0                            
                  GO TO 70                                              
   80             CONTINUE                                              
   90             IF (DABS(DET(1)) .LT. TEN) GO TO 100                  
                     DET(1) = DET(1)/TEN                                
                     DET(2) = DET(2) + 1.0D0                            
                  GO TO 90                                              
  100             CONTINUE                                              
  110          CONTINUE                                                 
  120       CONTINUE                                                    
  130    CONTINUE                                                       
  140 CONTINUE                                                          
!                                                                       
!     COMPUTE INVERSE(A)                                                
!                                                                       
      IF (NOINV) GO TO 270                                              
         K = 1                                                          
  150    IF (K .GT. N) GO TO 260                                        
            KM1 = K - 1                                                 
            IF (KPVT(K) .LT. 0) GO TO 180                               
!                                                                       
!              1 BY 1                                                   
!                                                                       
               A(K,K) = 1.0D0/A(K,K)                                    
               IF (KM1 .LT. 1) GO TO 170                                
                  CALL DCOPY(KM1,A(1,K),1,WORK,1)                       
                  DO 160 J = 1, KM1                                     
                     A(J,K) = DDOT(J,A(1,J),1,WORK,1)                   
                     CALL DAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)          
  160             CONTINUE                                              
                  A(K,K) = A(K,K) + DDOT(KM1,WORK,1,A(1,K),1)           
  170          CONTINUE                                                 
               KSTEP = 1                                                
            GO TO 220                                                   
  180       CONTINUE                                                    
!                                                                       
!              2 BY 2                                                   
!                                                                       
               T = DABS(A(K,K+1))                                       
               AK = A(K,K)/T                                            
               AKP1 = A(K+1,K+1)/T                                      
               AKKP1 = A(K,K+1)/T                                       
               D = T*(AK*AKP1 - 1.0D0)                                  
               A(K,K) = AKP1/D                                          
               A(K+1,K+1) = AK/D                                        
               A(K,K+1) = -AKKP1/D                                      
               IF (KM1 .LT. 1) GO TO 210                                
                  CALL DCOPY(KM1,A(1,K+1),1,WORK,1)                     
                  DO 190 J = 1, KM1                                     
                     A(J,K+1) = DDOT(J,A(1,J),1,WORK,1)                 
                     CALL DAXPY(J-1,WORK(J),A(1,J),1,A(1,K+1),1)        
  190             CONTINUE                                              
                  A(K+1,K+1) = A(K+1,K+1) + DDOT(KM1,WORK,1,A(1,K+1),1) 
                  A(K,K+1) = A(K,K+1) + DDOT(KM1,A(1,K),1,A(1,K+1),1)   
                  CALL DCOPY(KM1,A(1,K),1,WORK,1)                       
                  DO 200 J = 1, KM1                                     
                     A(J,K) = DDOT(J,A(1,J),1,WORK,1)                   
                     CALL DAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)          
  200             CONTINUE                                              
                  A(K,K) = A(K,K) + DDOT(KM1,WORK,1,A(1,K),1)           
  210          CONTINUE                                                 
               KSTEP = 2                                                
  220       CONTINUE                                                    
!                                                                       
!           SWAP                                                        
!                                                                       
            KS = IABS(KPVT(K))                                          
            IF (KS .EQ. K) GO TO 250                                    
               CALL DSWAP(KS,A(1,KS),1,A(1,K),1)                        
               DO 230 JB = KS, K                                        
                  J = K + KS - JB                                       
                  TEMP = A(J,K)                                         
                  A(J,K) = A(KS,J)                                      
                  A(KS,J) = TEMP                                        
  230          CONTINUE                                                 
               IF (KSTEP .EQ. 1) GO TO 240                              
                  TEMP = A(KS,K+1)                                      
                  A(KS,K+1) = A(K,K+1)                                  
                  A(K,K+1) = TEMP                                       
  240          CONTINUE                                                 
  250       CONTINUE                                                    
            K = K + KSTEP                                               
         GO TO 150                                                      
  260    CONTINUE                                                       
  270 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DSPCO(AP,N,KPVT,RCOND,Z)                               
      INTEGER N,KPVT(1)                                                 
      DOUBLE PRECISION AP(1),Z(1)                                       
      DOUBLE PRECISION RCOND                                            
!                                                                       
!     DSPCO FACTORS A DOUBLE PRECISION SYMMETRIC MATRIX STORED IN       
!     PACKED FORM BY ELIMINATION WITH SYMMETRIC PIVOTING AND ESTIMATES  
!     THE CONDITION OF THE MATRIX.                                      
!                                                                       
!     IF  RCOND  IS NOT NEEDED, DSPFA IS SLIGHTLY FASTER.               
!     TO SOLVE  A*X = B , FOLLOW DSPCO BY DSPSL.                        
!     TO COMPUTE  INVERSE(A)*C , FOLLOW DSPCO BY DSPSL.                 
!     TO COMPUTE  INVERSE(A) , FOLLOW DSPCO BY DSPDI.                   
!     TO COMPUTE  DETERMINANT(A) , FOLLOW DSPCO BY DSPDI.               
!     TO COMPUTE  INERTIA(A), FOLLOW DSPCO BY DSPDI.                    
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        AP      DOUBLE PRECISION (N*(N+1)/2)                           
!                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE        
!                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY  
!                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .      
!                SEE COMMENTS BELOW FOR DETAILS.                        
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!     OUTPUT                                                            
!                                                                       
!        AP      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH      
!                WERE USED TO OBTAIN IT STORED IN PACKED FORM.          
!                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)     
!                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT         
!                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE            
!                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL            
!                WITH 1 BY 1 AND 2 BY 2 BLOCKS.                         
!                                                                       
!        KPVT    INTEGER(N)                                             
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
!     PACKED STORAGE                                                    
!                                                                       
!          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER            
!          TRIANGLE OF A SYMMETRIC MATRIX.                              
!                                                                       
!                K = 0                                                  
!                DO 20 J = 1, N                                         
!                   DO 10 I = 1, J                                      
!                      K = K + 1                                        
!                      AP(K) = A(I,J)                                   
!             10    CONTINUE                                            
!             20 CONTINUE                                               
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     LINPACK DSPFA                                                     
!     BLAS DAXPY,DDOT,DSCAL,DASUM                                       
!     FORTRAN DABS,DMAX1,IABS,DSIGN                                     
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DDOT,DENOM,EK,T                  
      DOUBLE PRECISION ANORM,S,DASUM,YNORM                              
      INTEGER I,IJ,IK,IKM1,IKP1,INFO,J,JM1,J1                           
      INTEGER K,KK,KM1K,KM1KM1,KP,KPS,KS                                
!                                                                       
!                                                                       
!     FIND NORM OF A USING ONLY UPPER HALF                              
!                                                                       
      J1 = 1                                                            
      DO 30 J = 1, N                                                    
         Z(J) = DASUM(J,AP(J1),1)                                       
         IJ = J1                                                        
         J1 = J1 + J                                                    
         JM1 = J - 1                                                    
         IF (JM1 .LT. 1) GO TO 20                                       
         DO 10 I = 1, JM1                                               
            Z(I) = Z(I) + DABS(AP(IJ))                                  
            IJ = IJ + 1                                                 
   10    CONTINUE                                                       
   20    CONTINUE                                                       
   30 CONTINUE                                                          
      ANORM = 0.0D0                                                     
      DO 40 J = 1, N                                                    
         ANORM = DMAX1(ANORM,Z(J))                                      
   40 CONTINUE                                                          
!                                                                       
!     FACTOR                                                            
!                                                                       
      CALL DSPFA(AP,N,KPVT,INFO)                                        
!                                                                       
!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .              
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .         
!     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL           
!     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E .                   
!     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.            
!                                                                       
!     SOLVE U*D*W = E                                                   
!                                                                       
      EK = 1.0D0                                                        
      DO 50 J = 1, N                                                    
         Z(J) = 0.0D0                                                   
   50 CONTINUE                                                          
      K = N                                                             
      IK = (N*(N - 1))/2                                                
   60 IF (K .EQ. 0) GO TO 120                                           
         KK = IK + K                                                    
         IKM1 = IK - (K - 1)                                            
         KS = 1                                                         
         IF (KPVT(K) .LT. 0) KS = 2                                     
         KP = IABS(KPVT(K))                                             
         KPS = K + 1 - KS                                               
         IF (KP .EQ. KPS) GO TO 70                                      
            T = Z(KPS)                                                  
            Z(KPS) = Z(KP)                                              
            Z(KP) = T                                                   
   70    CONTINUE                                                       
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,Z(K))                       
         Z(K) = Z(K) + EK                                               
         CALL DAXPY(K-KS,Z(K),AP(IK+1),1,Z(1),1)                        
         IF (KS .EQ. 1) GO TO 80                                        
            IF (Z(K-1) .NE. 0.0D0) EK = DSIGN(EK,Z(K-1))                
            Z(K-1) = Z(K-1) + EK                                        
            CALL DAXPY(K-KS,Z(K-1),AP(IKM1+1),1,Z(1),1)                 
   80    CONTINUE                                                       
         IF (KS .EQ. 2) GO TO 100                                       
            IF (DABS(Z(K)) .LE. DABS(AP(KK))) GO TO 90                  
               S = DABS(AP(KK))/DABS(Z(K))                              
               CALL DSCAL(N,S,Z,1)                                      
               EK = S*EK                                                
   90       CONTINUE                                                    
            IF (AP(KK) .NE. 0.0D0) Z(K) = Z(K)/AP(KK)                   
            IF (AP(KK) .EQ. 0.0D0) Z(K) = 1.0D0                         
         GO TO 110                                                      
  100    CONTINUE                                                       
            KM1K = IK + K - 1                                           
            KM1KM1 = IKM1 + K - 1                                       
            AK = AP(KK)/AP(KM1K)                                        
            AKM1 = AP(KM1KM1)/AP(KM1K)                                  
            BK = Z(K)/AP(KM1K)                                          
            BKM1 = Z(K-1)/AP(KM1K)                                      
            DENOM = AK*AKM1 - 1.0D0                                     
            Z(K) = (AKM1*BK - BKM1)/DENOM                               
            Z(K-1) = (AK*BKM1 - BK)/DENOM                               
  110    CONTINUE                                                       
         K = K - KS                                                     
         IK = IK - K                                                    
         IF (KS .EQ. 2) IK = IK - (K + 1)                               
      GO TO 60                                                          
  120 CONTINUE                                                          
      S = 1.0D0/DASUM(N,Z,1)                                            
      CALL DSCAL(N,S,Z,1)                                               
!                                                                       
!     SOLVE TRANS(U)*Y = W                                              
!                                                                       
      K = 1                                                             
      IK = 0                                                            
  130 IF (K .GT. N) GO TO 160                                           
         KS = 1                                                         
         IF (KPVT(K) .LT. 0) KS = 2                                     
         IF (K .EQ. 1) GO TO 150                                        
            Z(K) = Z(K) + DDOT(K-1,AP(IK+1),1,Z(1),1)                   
            IKP1 = IK + K                                               
            IF (KS .EQ. 2)                                               &
               Z(K+1) = Z(K+1) + DDOT(K-1,AP(IKP1+1),1,Z(1),1)          
            KP = IABS(KPVT(K))                                          
            IF (KP .EQ. K) GO TO 140                                    
               T = Z(K)                                                 
               Z(K) = Z(KP)                                             
               Z(KP) = T                                                
  140       CONTINUE                                                    
  150    CONTINUE                                                       
         IK = IK + K                                                    
         IF (KS .EQ. 2) IK = IK + (K + 1)                               
         K = K + KS                                                     
      GO TO 130                                                         
  160 CONTINUE                                                          
      S = 1.0D0/DASUM(N,Z,1)                                            
      CALL DSCAL(N,S,Z,1)                                               
!                                                                       
      YNORM = 1.0D0                                                     
!                                                                       
!     SOLVE U*D*V = Y                                                   
!                                                                       
      K = N                                                             
      IK = N*(N - 1)/2                                                  
  170 IF (K .EQ. 0) GO TO 230                                           
         KK = IK + K                                                    
         IKM1 = IK - (K - 1)                                            
         KS = 1                                                         
         IF (KPVT(K) .LT. 0) KS = 2                                     
         IF (K .EQ. KS) GO TO 190                                       
            KP = IABS(KPVT(K))                                          
            KPS = K + 1 - KS                                            
            IF (KP .EQ. KPS) GO TO 180                                  
               T = Z(KPS)                                               
               Z(KPS) = Z(KP)                                           
               Z(KP) = T                                                
  180       CONTINUE                                                    
            CALL DAXPY(K-KS,Z(K),AP(IK+1),1,Z(1),1)                     
            IF (KS .EQ. 2) CALL DAXPY(K-KS,Z(K-1),AP(IKM1+1),1,Z(1),1)  
  190    CONTINUE                                                       
         IF (KS .EQ. 2) GO TO 210                                       
            IF (DABS(Z(K)) .LE. DABS(AP(KK))) GO TO 200                 
               S = DABS(AP(KK))/DABS(Z(K))                              
               CALL DSCAL(N,S,Z,1)                                      
               YNORM = S*YNORM                                          
  200       CONTINUE                                                    
            IF (AP(KK) .NE. 0.0D0) Z(K) = Z(K)/AP(KK)                   
            IF (AP(KK) .EQ. 0.0D0) Z(K) = 1.0D0                         
         GO TO 220                                                      
  210    CONTINUE                                                       
            KM1K = IK + K - 1                                           
            KM1KM1 = IKM1 + K - 1                                       
            AK = AP(KK)/AP(KM1K)                                        
            AKM1 = AP(KM1KM1)/AP(KM1K)                                  
            BK = Z(K)/AP(KM1K)                                          
            BKM1 = Z(K-1)/AP(KM1K)                                      
            DENOM = AK*AKM1 - 1.0D0                                     
            Z(K) = (AKM1*BK - BKM1)/DENOM                               
            Z(K-1) = (AK*BKM1 - BK)/DENOM                               
  220    CONTINUE                                                       
         K = K - KS                                                     
         IK = IK - K                                                    
         IF (KS .EQ. 2) IK = IK - (K + 1)                               
      GO TO 170                                                         
  230 CONTINUE                                                          
      S = 1.0D0/DASUM(N,Z,1)                                            
      CALL DSCAL(N,S,Z,1)                                               
      YNORM = S*YNORM                                                   
!                                                                       
!     SOLVE TRANS(U)*Z = V                                              
!                                                                       
      K = 1                                                             
      IK = 0                                                            
  240 IF (K .GT. N) GO TO 270                                           
         KS = 1                                                         
         IF (KPVT(K) .LT. 0) KS = 2                                     
         IF (K .EQ. 1) GO TO 260                                        
            Z(K) = Z(K) + DDOT(K-1,AP(IK+1),1,Z(1),1)                   
            IKP1 = IK + K                                               
            IF (KS .EQ. 2)                                               &
               Z(K+1) = Z(K+1) + DDOT(K-1,AP(IKP1+1),1,Z(1),1)          
            KP = IABS(KPVT(K))                                          
            IF (KP .EQ. K) GO TO 250                                    
               T = Z(K)                                                 
               Z(K) = Z(KP)                                             
               Z(KP) = T                                                
  250       CONTINUE                                                    
  260    CONTINUE                                                       
         IK = IK + K                                                    
         IF (KS .EQ. 2) IK = IK + (K + 1)                               
         K = K + KS                                                     
      GO TO 240                                                         
  270 CONTINUE                                                          
!     MAKE ZNORM = 1.0                                                  
      S = 1.0D0/DASUM(N,Z,1)                                            
      CALL DSCAL(N,S,Z,1)                                               
      YNORM = S*YNORM                                                   
!                                                                       
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM                         
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0                               
      RETURN                                                            
      END                                                               
      SUBROUTINE DSPFA(AP,N,KPVT,INFO)                                  
      INTEGER N,KPVT(1),INFO                                            
      DOUBLE PRECISION AP(1)                                            
!                                                                       
!     DSPFA FACTORS A DOUBLE PRECISION SYMMETRIC MATRIX STORED IN       
!     PACKED FORM BY ELIMINATION WITH SYMMETRIC PIVOTING.               
!                                                                       
!     TO SOLVE  A*X = B , FOLLOW DSPFA BY DSPSL.                        
!     TO COMPUTE  INVERSE(A)*C , FOLLOW DSPFA BY DSPSL.                 
!     TO COMPUTE  DETERMINANT(A) , FOLLOW DSPFA BY DSPDI.               
!     TO COMPUTE  INERTIA(A) , FOLLOW DSPFA BY DSPDI.                   
!     TO COMPUTE  INVERSE(A) , FOLLOW DSPFA BY DSPDI.                   
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        AP      DOUBLE PRECISION (N*(N+1)/2)                           
!                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE        
!                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY  
!                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .      
!                SEE COMMENTS BELOW FOR DETAILS.                        
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!     OUTPUT                                                            
!                                                                       
!        AP      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH      
!                WERE USED TO OBTAIN IT STORED IN PACKED FORM.          
!                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)     
!                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT         
!                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE            
!                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL            
!                WITH 1 BY 1 AND 2 BY 2 BLOCKS.                         
!                                                                       
!        KPVT    INTEGER(N)                                             
!                AN INTEGER VECTOR OF PIVOT INDICES.                    
!                                                                       
!        INFO    INTEGER                                                
!                = 0  NORMAL VALUE.                                     
!                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS      
!                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE,       
!                     BUT IT DOES INDICATE THAT DSPSL OR DSPDI MAY      
!                     DIVIDE BY ZERO IF CALLED.                         
!                                                                       
!     PACKED STORAGE                                                    
!                                                                       
!          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER            
!          TRIANGLE OF A SYMMETRIC MATRIX.                              
!                                                                       
!                K = 0                                                  
!                DO 20 J = 1, N                                         
!                   DO 10 I = 1, J                                      
!                      K = K + 1                                        
!                      AP(K)  = A(I,J)                                  
!             10    CONTINUE                                            
!             20 CONTINUE                                               
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.            
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DAXPY,DSWAP,IDAMAX                                           
!     FORTRAN DABS,DMAX1,DSQRT                                          
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T              
      DOUBLE PRECISION ABSAKK,ALPHA,COLMAX,ROWMAX                       
      INTEGER IDAMAX,IJ,IJJ,IK,IKM1,IM,IMAX,IMAXP1,IMIM,IMJ,IMK         
      INTEGER J,JJ,JK,JKM1,JMAX,JMIM,K,KK,KM1,KM1K,KM1KM1,KM2,KSTEP     
      LOGICAL SWAP                                                      
!                                                                       
!                                                                       
!     INITIALIZE                                                        
!                                                                       
!     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.                       
      ALPHA = (1.0D0 + DSQRT(17.0D0))/8.0D0                             
!                                                                       
      INFO = 0                                                          
!                                                                       
!     MAIN LOOP ON K, WHICH GOES FROM N TO 1.                           
!                                                                       
      K = N                                                             
      IK = (N*(N - 1))/2                                                
   10 CONTINUE                                                          
!                                                                       
!        LEAVE THE LOOP IF K=0 OR K=1.                                  
!                                                                       
!     ...EXIT                                                           
         IF (K .EQ. 0) GO TO 200                                        
         IF (K .GT. 1) GO TO 20                                         
            KPVT(1) = 1                                                 
            IF (AP(1) .EQ. 0.0D0) INFO = 1                              
!     ......EXIT                                                        
            GO TO 200                                                   
   20    CONTINUE                                                       
!                                                                       
!        THIS SECTION OF CODE DETERMINES THE KIND OF                    
!        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,            
!        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND          
!        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS                
!        REQUIRED.                                                      
!                                                                       
         KM1 = K - 1                                                    
         KK = IK + K                                                    
         ABSAKK = DABS(AP(KK))                                          
!                                                                       
!        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN                  
!        COLUMN K.                                                      
!                                                                       
         IMAX = IDAMAX(K-1,AP(IK+1),1)                                  
         IMK = IK + IMAX                                                
         COLMAX = DABS(AP(IMK))                                         
         IF (ABSAKK .LT. ALPHA*COLMAX) GO TO 30                         
            KSTEP = 1                                                   
            SWAP = .FALSE.                                              
         GO TO 90                                                       
   30    CONTINUE                                                       
!                                                                       
!           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN               
!           ROW IMAX.                                                   
!                                                                       
            ROWMAX = 0.0D0                                              
            IMAXP1 = IMAX + 1                                           
            IM = IMAX*(IMAX - 1)/2                                      
            IMJ = IM + 2*IMAX                                           
            DO 40 J = IMAXP1, K                                         
               ROWMAX = DMAX1(ROWMAX,DABS(AP(IMJ)))                     
               IMJ = IMJ + J                                            
   40       CONTINUE                                                    
            IF (IMAX .EQ. 1) GO TO 50                                   
               JMAX = IDAMAX(IMAX-1,AP(IM+1),1)                         
               JMIM = JMAX + IM                                         
               ROWMAX = DMAX1(ROWMAX,DABS(AP(JMIM)))                    
   50       CONTINUE                                                    
            IMIM = IMAX + IM                                            
            IF (DABS(AP(IMIM)) .LT. ALPHA*ROWMAX) GO TO 60              
               KSTEP = 1                                                
               SWAP = .TRUE.                                            
            GO TO 80                                                    
   60       CONTINUE                                                    
            IF (ABSAKK .LT. ALPHA*COLMAX*(COLMAX/ROWMAX)) GO TO 70      
               KSTEP = 1                                                
               SWAP = .FALSE.                                           
            GO TO 80                                                    
   70       CONTINUE                                                    
               KSTEP = 2                                                
               SWAP = IMAX .NE. KM1                                     
   80       CONTINUE                                                    
   90    CONTINUE                                                       
         IF (DMAX1(ABSAKK,COLMAX) .NE. 0.0D0) GO TO 100                 
!                                                                       
!           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.           
!                                                                       
            KPVT(K) = K                                                 
            INFO = K                                                    
         GO TO 190                                                      
  100    CONTINUE                                                       
         IF (KSTEP .EQ. 2) GO TO 140                                    
!                                                                       
!           1 X 1 PIVOT BLOCK.                                          
!                                                                       
            IF (.NOT.SWAP) GO TO 120                                    
!                                                                       
!              PERFORM AN INTERCHANGE.                                  
!                                                                       
               CALL DSWAP(IMAX,AP(IM+1),1,AP(IK+1),1)                   
               IMJ = IK + IMAX                                          
               DO 110 JJ = IMAX, K                                      
                  J = K + IMAX - JJ                                     
                  JK = IK + J                                           
                  T = AP(JK)                                            
                  AP(JK) = AP(IMJ)                                      
                  AP(IMJ) = T                                           
                  IMJ = IMJ - (J - 1)                                   
  110          CONTINUE                                                 
  120       CONTINUE                                                    
!                                                                       
!           PERFORM THE ELIMINATION.                                    
!                                                                       
            IJ = IK - (K - 1)                                           
            DO 130 JJ = 1, KM1                                          
               J = K - JJ                                               
               JK = IK + J                                              
               MULK = -AP(JK)/AP(KK)                                    
               T = MULK                                                 
               CALL DAXPY(J,T,AP(IK+1),1,AP(IJ+1),1)                    
               IJJ = IJ + J                                             
               AP(JK) = MULK                                            
               IJ = IJ - (J - 1)                                        
  130       CONTINUE                                                    
!                                                                       
!           SET THE PIVOT ARRAY.                                        
!                                                                       
            KPVT(K) = K                                                 
            IF (SWAP) KPVT(K) = IMAX                                    
         GO TO 190                                                      
  140    CONTINUE                                                       
!                                                                       
!           2 X 2 PIVOT BLOCK.                                          
!                                                                       
            KM1K = IK + K - 1                                           
            IKM1 = IK - (K - 1)                                         
            IF (.NOT.SWAP) GO TO 160                                    
!                                                                       
!              PERFORM AN INTERCHANGE.                                  
!                                                                       
               CALL DSWAP(IMAX,AP(IM+1),1,AP(IKM1+1),1)                 
               IMJ = IKM1 + IMAX                                        
               DO 150 JJ = IMAX, KM1                                    
                  J = KM1 + IMAX - JJ                                   
                  JKM1 = IKM1 + J                                       
                  T = AP(JKM1)                                          
                  AP(JKM1) = AP(IMJ)                                    
                  AP(IMJ) = T                                           
                  IMJ = IMJ - (J - 1)                                   
  150          CONTINUE                                                 
               T = AP(KM1K)                                             
               AP(KM1K) = AP(IMK)                                       
               AP(IMK) = T                                              
  160       CONTINUE                                                    
!                                                                       
!           PERFORM THE ELIMINATION.                                    
!                                                                       
            KM2 = K - 2                                                 
            IF (KM2 .EQ. 0) GO TO 180                                   
               AK = AP(KK)/AP(KM1K)                                     
               KM1KM1 = IKM1 + K - 1                                    
               AKM1 = AP(KM1KM1)/AP(KM1K)                               
               DENOM = 1.0D0 - AK*AKM1                                  
               IJ = IK - (K - 1) - (K - 2)                              
               DO 170 JJ = 1, KM2                                       
                  J = KM1 - JJ                                          
                  JK = IK + J                                           
                  BK = AP(JK)/AP(KM1K)                                  
                  JKM1 = IKM1 + J                                       
                  BKM1 = AP(JKM1)/AP(KM1K)                              
                  MULK = (AKM1*BK - BKM1)/DENOM                         
                  MULKM1 = (AK*BKM1 - BK)/DENOM                         
                  T = MULK                                              
                  CALL DAXPY(J,T,AP(IK+1),1,AP(IJ+1),1)                 
                  T = MULKM1                                            
                  CALL DAXPY(J,T,AP(IKM1+1),1,AP(IJ+1),1)               
                  AP(JK) = MULK                                         
                  AP(JKM1) = MULKM1                                     
                  IJJ = IJ + J                                          
                  IJ = IJ - (J - 1)                                     
  170          CONTINUE                                                 
  180       CONTINUE                                                    
!                                                                       
!           SET THE PIVOT ARRAY.                                        
!                                                                       
            KPVT(K) = 1 - K                                             
            IF (SWAP) KPVT(K) = -IMAX                                   
            KPVT(K-1) = KPVT(K)                                         
  190    CONTINUE                                                       
         IK = IK - (K - 1)                                              
         IF (KSTEP .EQ. 2) IK = IK - (K - 2)                            
         K = K - KSTEP                                                  
      GO TO 10                                                          
  200 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DSPSL(AP,N,KPVT,B)                                     
      INTEGER N,KPVT(1)                                                 
      DOUBLE PRECISION AP(1),B(1)                                       
!                                                                       
!     DSISL SOLVES THE DOUBLE PRECISION SYMMETRIC SYSTEM                
!     A * X = B                                                         
!     USING THE FACTORS COMPUTED BY DSPFA.                              
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        AP      DOUBLE PRECISION(N*(N+1)/2)                            
!                THE OUTPUT FROM DSPFA.                                 
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!        KPVT    INTEGER(N)                                             
!                THE PIVOT VECTOR FROM DSPFA.                           
!                                                                       
!        B       DOUBLE PRECISION(N)                                    
!                THE RIGHT HAND SIDE VECTOR.                            
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        B       THE SOLUTION VECTOR  X .                               
!                                                                       
!     ERROR CONDITION                                                   
!                                                                       
!        A DIVISION BY ZERO MAY OCCUR IF  DSPCO  HAS SET RCOND .EQ. 0.0 
!        OR  DSPFA  HAS SET INFO .NE. 0  .                              
!                                                                       
!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 
!     WITH  P  COLUMNS                                                  
!           CALL DSPFA(AP,N,KPVT,INFO)                                  
!           IF (INFO .NE. 0) GO TO ...                                  
!           DO 10 J = 1, P                                              
!              CALL DSPSL(AP,N,KPVT,C(1,J))                             
!        10 CONTINUE                                                    
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.            
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DAXPY,DDOT                                                   
!     FORTRAN IABS                                                      
!                                                                       
!     INTERNAL VARIABLES.                                               
!                                                                       
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DDOT,DENOM,TEMP                  
      INTEGER IK,IKM1,IKP1,K,KK,KM1K,KM1KM1,KP                          
!                                                                       
!     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND                    
!     D INVERSE TO B.                                                   
!                                                                       
      K = N                                                             
      IK = (N*(N - 1))/2                                                
   10 IF (K .EQ. 0) GO TO 80                                            
         KK = IK + K                                                    
         IF (KPVT(K) .LT. 0) GO TO 40                                   
!                                                                       
!           1 X 1 PIVOT BLOCK.                                          
!                                                                       
            IF (K .EQ. 1) GO TO 30                                      
               KP = KPVT(K)                                             
               IF (KP .EQ. K) GO TO 20                                  
!                                                                       
!                 INTERCHANGE.                                          
!                                                                       
                  TEMP = B(K)                                           
                  B(K) = B(KP)                                          
                  B(KP) = TEMP                                          
   20          CONTINUE                                                 
!                                                                       
!              APPLY THE TRANSFORMATION.                                
!                                                                       
               CALL DAXPY(K-1,B(K),AP(IK+1),1,B(1),1)                   
   30       CONTINUE                                                    
!                                                                       
!           APPLY D INVERSE.                                            
!                                                                       
            B(K) = B(K)/AP(KK)                                          
            K = K - 1                                                   
            IK = IK - K                                                 
         GO TO 70                                                       
   40    CONTINUE                                                       
!                                                                       
!           2 X 2 PIVOT BLOCK.                                          
!                                                                       
            IKM1 = IK - (K - 1)                                         
            IF (K .EQ. 2) GO TO 60                                      
               KP = IABS(KPVT(K))                                       
               IF (KP .EQ. K - 1) GO TO 50                              
!                                                                       
!                 INTERCHANGE.                                          
!                                                                       
                  TEMP = B(K-1)                                         
                  B(K-1) = B(KP)                                        
                  B(KP) = TEMP                                          
   50          CONTINUE                                                 
!                                                                       
!              APPLY THE TRANSFORMATION.                                
!                                                                       
               CALL DAXPY(K-2,B(K),AP(IK+1),1,B(1),1)                   
               CALL DAXPY(K-2,B(K-1),AP(IKM1+1),1,B(1),1)               
   60       CONTINUE                                                    
!                                                                       
!           APPLY D INVERSE.                                            
!                                                                       
            KM1K = IK + K - 1                                           
            KK = IK + K                                                 
            AK = AP(KK)/AP(KM1K)                                        
            KM1KM1 = IKM1 + K - 1                                       
            AKM1 = AP(KM1KM1)/AP(KM1K)                                  
            BK = B(K)/AP(KM1K)                                          
            BKM1 = B(K-1)/AP(KM1K)                                      
            DENOM = AK*AKM1 - 1.0D0                                     
            B(K) = (AKM1*BK - BKM1)/DENOM                               
            B(K-1) = (AK*BKM1 - BK)/DENOM                               
            K = K - 2                                                   
            IK = IK - (K + 1) - K                                       
   70    CONTINUE                                                       
      GO TO 10                                                          
   80 CONTINUE                                                          
!                                                                       
!     LOOP FORWARD APPLYING THE TRANSFORMATIONS.                        
!                                                                       
      K = 1                                                             
      IK = 0                                                            
   90 IF (K .GT. N) GO TO 160                                           
         IF (KPVT(K) .LT. 0) GO TO 120                                  
!                                                                       
!           1 X 1 PIVOT BLOCK.                                          
!                                                                       
            IF (K .EQ. 1) GO TO 110                                     
!                                                                       
!              APPLY THE TRANSFORMATION.                                
!                                                                       
               B(K) = B(K) + DDOT(K-1,AP(IK+1),1,B(1),1)                
               KP = KPVT(K)                                             
               IF (KP .EQ. K) GO TO 100                                 
!                                                                       
!                 INTERCHANGE.                                          
!                                                                       
                  TEMP = B(K)                                           
                  B(K) = B(KP)                                          
                  B(KP) = TEMP                                          
  100          CONTINUE                                                 
  110       CONTINUE                                                    
            IK = IK + K                                                 
            K = K + 1                                                   
         GO TO 150                                                      
  120    CONTINUE                                                       
!                                                                       
!           2 X 2 PIVOT BLOCK.                                          
!                                                                       
            IF (K .EQ. 1) GO TO 140                                     
!                                                                       
!              APPLY THE TRANSFORMATION.                                
!                                                                       
               B(K) = B(K) + DDOT(K-1,AP(IK+1),1,B(1),1)                
               IKP1 = IK + K                                            
               B(K+1) = B(K+1) + DDOT(K-1,AP(IKP1+1),1,B(1),1)          
               KP = IABS(KPVT(K))                                       
               IF (KP .EQ. K) GO TO 130                                 
!                                                                       
!                 INTERCHANGE.                                          
!                                                                       
                  TEMP = B(K)                                           
                  B(K) = B(KP)                                          
                  B(KP) = TEMP                                          
  130          CONTINUE                                                 
  140       CONTINUE                                                    
            IK = IK + K + K + 1                                         
            K = K + 2                                                   
  150    CONTINUE                                                       
      GO TO 90                                                          
  160 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DSPDI(AP,N,KPVT,DET,INERT,WORK,JOB)                    
      INTEGER N,JOB                                                     
      DOUBLE PRECISION AP(1),WORK(1)                                    
      DOUBLE PRECISION DET(2)                                           
      INTEGER KPVT(1),INERT(3)                                          
!                                                                       
!     DSPDI COMPUTES THE DETERMINANT, INERTIA AND INVERSE               
!     OF A DOUBLE PRECISION SYMMETRIC MATRIX USING THE FACTORS FROM     
!     DSPFA, WHERE THE MATRIX IS STORED IN PACKED FORM.                 
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        AP      DOUBLE PRECISION (N*(N+1)/2)                           
!                THE OUTPUT FROM DSPFA.                                 
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX A.                             
!                                                                       
!        KPVT    INTEGER(N)                                             
!                THE PIVOT VECTOR FROM DSPFA.                           
!                                                                       
!        WORK    DOUBLE PRECISION(N)                                    
!                WORK VECTOR.  CONTENTS IGNORED.                        
!                                                                       
!        JOB     INTEGER                                                
!                JOB HAS THE DECIMAL EXPANSION  ABC  WHERE              
!                   IF  C .NE. 0, THE INVERSE IS COMPUTED,              
!                   IF  B .NE. 0, THE DETERMINANT IS COMPUTED,          
!                   IF  A .NE. 0, THE INERTIA IS COMPUTED.              
!                                                                       
!                FOR EXAMPLE, JOB = 111  GIVES ALL THREE.               
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        VARIABLES NOT REQUESTED BY JOB ARE NOT USED.                   
!                                                                       
!        AP     CONTAINS THE UPPER TRIANGLE OF THE INVERSE OF           
!               THE ORIGINAL MATRIX, STORED IN PACKED FORM.             
!               THE COLUMNS OF THE UPPER TRIANGLE ARE STORED            
!               SEQUENTIALLY IN A ONE-DIMENSIONAL ARRAY.                
!                                                                       
!        DET    DOUBLE PRECISION(2)                                     
!               DETERMINANT OF ORIGINAL MATRIX.                         
!               DETERMINANT = DET(1) * 10.0**DET(2)                     
!               WITH 1.0 .LE. DABS(DET(1)) .LT. 10.0                    
!               OR DET(1) = 0.0.                                        
!                                                                       
!        INERT  INTEGER(3)                                              
!               THE INERTIA OF THE ORIGINAL MATRIX.                     
!               INERT(1)  =  NUMBER OF POSITIVE EIGENVALUES.            
!               INERT(2)  =  NUMBER OF NEGATIVE EIGENVALUES.            
!               INERT(3)  =  NUMBER OF ZERO EIGENVALUES.                
!                                                                       
!     ERROR CONDITION                                                   
!                                                                       
!        A DIVISION BY ZERO WILL OCCUR IF THE INVERSE IS REQUESTED      
!        AND  DSPCO  HAS SET RCOND .EQ. 0.0                             
!        OR  DSPFA  HAS SET  INFO .NE. 0 .                              
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.            
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DAXPY,DCOPY,DDOT,DSWAP                                       
!     FORTRAN DABS,IABS,MOD                                             
!                                                                       
!     INTERNAL VARIABLES.                                               
!                                                                       
      DOUBLE PRECISION AKKP1,DDOT,TEMP                                  
      DOUBLE PRECISION TEN,D,T,AK,AKP1                                  
      INTEGER IJ,IK,IKP1,IKS,J,JB,JK,JKP1                               
      INTEGER K,KK,KKP1,KM1,KS,KSJ,KSKP1,KSTEP                          
      LOGICAL NOINV,NODET,NOERT                                         
!                                                                       
      NOINV = MOD(JOB,10) .EQ. 0                                        
      NODET = MOD(JOB,100)/10 .EQ. 0                                    
      NOERT = MOD(JOB,1000)/100 .EQ. 0                                  
!                                                                       
      IF (NODET .AND. NOERT) GO TO 140                                  
         IF (NOERT) GO TO 10                                            
            INERT(1) = 0                                                
            INERT(2) = 0                                                
            INERT(3) = 0                                                
   10    CONTINUE                                                       
         IF (NODET) GO TO 20                                            
            DET(1) = 1.0D0                                              
            DET(2) = 0.0D0                                              
            TEN = 10.0D0                                                
   20    CONTINUE                                                       
         T = 0.0D0                                                      
         IK = 0                                                         
         DO 130 K = 1, N                                                
            KK = IK + K                                                 
            D = AP(KK)                                                  
!                                                                       
!           CHECK IF 1 BY 1                                             
!                                                                       
            IF (KPVT(K) .GT. 0) GO TO 50                                
!                                                                       
!              2 BY 2 BLOCK                                             
!              USE DET (D  S)  =  (D/T * C - T) * T  ,  T = DABS(S)     
!                      (S  C)                                           
!              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.                    
!              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.      
!                                                                       
               IF (T .NE. 0.0D0) GO TO 30                               
                  IKP1 = IK + K                                         
                  KKP1 = IKP1 + K                                       
                  T = DABS(AP(KKP1))                                    
                  D = (D/T)*AP(KKP1+1) - T                              
               GO TO 40                                                 
   30          CONTINUE                                                 
                  D = T                                                 
                  T = 0.0D0                                             
   40          CONTINUE                                                 
   50       CONTINUE                                                    
!                                                                       
            IF (NOERT) GO TO 60                                         
               IF (D .GT. 0.0D0) INERT(1) = INERT(1) + 1                
               IF (D .LT. 0.0D0) INERT(2) = INERT(2) + 1                
               IF (D .EQ. 0.0D0) INERT(3) = INERT(3) + 1                
   60       CONTINUE                                                    
!                                                                       
            IF (NODET) GO TO 120                                        
               DET(1) = D*DET(1)                                        
               IF (DET(1) .EQ. 0.0D0) GO TO 110                         
   70             IF (DABS(DET(1)) .GE. 1.0D0) GO TO 80                 
                     DET(1) = TEN*DET(1)                                
                     DET(2) = DET(2) - 1.0D0                            
                  GO TO 70                                              
   80             CONTINUE                                              
   90             IF (DABS(DET(1)) .LT. TEN) GO TO 100                  
                     DET(1) = DET(1)/TEN                                
                     DET(2) = DET(2) + 1.0D0                            
                  GO TO 90                                              
  100             CONTINUE                                              
  110          CONTINUE                                                 
  120       CONTINUE                                                    
            IK = IK + K                                                 
  130    CONTINUE                                                       
  140 CONTINUE                                                          
!                                                                       
!     COMPUTE INVERSE(A)                                                
!                                                                       
      IF (NOINV) GO TO 270                                              
         K = 1                                                          
         IK = 0                                                         
  150    IF (K .GT. N) GO TO 260                                        
            KM1 = K - 1                                                 
            KK = IK + K                                                 
            IKP1 = IK + K                                               
            KKP1 = IKP1 + K                                             
            IF (KPVT(K) .LT. 0) GO TO 180                               
!                                                                       
!              1 BY 1                                                   
!                                                                       
               AP(KK) = 1.0D0/AP(KK)                                    
               IF (KM1 .LT. 1) GO TO 170                                
                  CALL DCOPY(KM1,AP(IK+1),1,WORK,1)                     
                  IJ = 0                                                
                  DO 160 J = 1, KM1                                     
                     JK = IK + J                                        
                     AP(JK) = DDOT(J,AP(IJ+1),1,WORK,1)                 
                     CALL DAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IK+1),1)      
                     IJ = IJ + J                                        
  160             CONTINUE                                              
                  AP(KK) = AP(KK) + DDOT(KM1,WORK,1,AP(IK+1),1)         
  170          CONTINUE                                                 
               KSTEP = 1                                                
            GO TO 220                                                   
  180       CONTINUE                                                    
!                                                                       
!              2 BY 2                                                   
!                                                                       
               T = DABS(AP(KKP1))                                       
               AK = AP(KK)/T                                            
               AKP1 = AP(KKP1+1)/T                                      
               AKKP1 = AP(KKP1)/T                                       
               D = T*(AK*AKP1 - 1.0D0)                                  
               AP(KK) = AKP1/D                                          
               AP(KKP1+1) = AK/D                                        
               AP(KKP1) = -AKKP1/D                                      
               IF (KM1 .LT. 1) GO TO 210                                
                  CALL DCOPY(KM1,AP(IKP1+1),1,WORK,1)                   
                  IJ = 0                                                
                  DO 190 J = 1, KM1                                     
                     JKP1 = IKP1 + J                                    
                     AP(JKP1) = DDOT(J,AP(IJ+1),1,WORK,1)               
                     CALL DAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IKP1+1),1)    
                     IJ = IJ + J                                        
  190             CONTINUE                                              
                  AP(KKP1+1) = AP(KKP1+1)                                &
                               + DDOT(KM1,WORK,1,AP(IKP1+1),1)          
                  AP(KKP1) = AP(KKP1)                                    &
                             + DDOT(KM1,AP(IK+1),1,AP(IKP1+1),1)        
                  CALL DCOPY(KM1,AP(IK+1),1,WORK,1)                     
                  IJ = 0                                                
                  DO 200 J = 1, KM1                                     
                     JK = IK + J                                        
                     AP(JK) = DDOT(J,AP(IJ+1),1,WORK,1)                 
                     CALL DAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IK+1),1)      
                     IJ = IJ + J                                        
  200             CONTINUE                                              
                  AP(KK) = AP(KK) + DDOT(KM1,WORK,1,AP(IK+1),1)         
  210          CONTINUE                                                 
               KSTEP = 2                                                
  220       CONTINUE                                                    
!                                                                       
!           SWAP                                                        
!                                                                       
            KS = IABS(KPVT(K))                                          
            IF (KS .EQ. K) GO TO 250                                    
               IKS = (KS*(KS - 1))/2                                    
               CALL DSWAP(KS,AP(IKS+1),1,AP(IK+1),1)                    
               KSJ = IK + KS                                            
               DO 230 JB = KS, K                                        
                  J = K + KS - JB                                       
                  JK = IK + J                                           
                  TEMP = AP(JK)                                         
                  AP(JK) = AP(KSJ)                                      
                  AP(KSJ) = TEMP                                        
                  KSJ = KSJ - (J - 1)                                   
  230          CONTINUE                                                 
               IF (KSTEP .EQ. 1) GO TO 240                              
                  KSKP1 = IKP1 + KS                                     
                  TEMP = AP(KSKP1)                                      
                  AP(KSKP1) = AP(KKP1)                                  
                  AP(KKP1) = TEMP                                       
  240          CONTINUE                                                 
  250       CONTINUE                                                    
            IK = IK + K                                                 
            IF (KSTEP .EQ. 2) IK = IK + K + 1                           
            K = K + KSTEP                                               
         GO TO 150                                                      
  260    CONTINUE                                                       
  270 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DTRCO(T,LDT,N,RCOND,Z,JOB)                             
      INTEGER LDT,N,JOB                                                 
      DOUBLE PRECISION T(LDT,1),Z(1)                                    
      DOUBLE PRECISION RCOND                                            
!                                                                       
!     DTRCO ESTIMATES THE CONDITION OF A DOUBLE PRECISION TRIANGULAR    
!     MATRIX.                                                           
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        T       DOUBLE PRECISION(LDT,N)                                
!                T CONTAINS THE TRIANGULAR MATRIX. THE ZERO             
!                ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND         
!                THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE         
!                USED TO STORE OTHER INFORMATION.                       
!                                                                       
!        LDT     INTEGER                                                
!                LDT IS THE LEADING DIMENSION OF THE ARRAY T.           
!                                                                       
!        N       INTEGER                                                
!                N IS THE ORDER OF THE SYSTEM.                          
!                                                                       
!        JOB     INTEGER                                                
!                = 0         T  IS LOWER TRIANGULAR.                    
!                = NONZERO   T  IS UPPER TRIANGULAR.                    
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        RCOND   DOUBLE PRECISION                                       
!                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  T .        
!                FOR THE SYSTEM  T*X = B , RELATIVE PERTURBATIONS       
!                IN  T  AND  B  OF SIZE  EPSILON  MAY CAUSE             
!                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . 
!                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION     
!                           1.0 + RCOND .EQ. 1.0                        
!                IS TRUE, THEN  T  MAY BE SINGULAR TO WORKING           
!                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF         
!                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE          
!                UNDERFLOWS.                                            
!                                                                       
!        Z       DOUBLE PRECISION(N)                                    
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.  
!                IF  T  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS      
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT           
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .                    
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DAXPY,DSCAL,DASUM                                            
!     FORTRAN DABS,DMAX1,DSIGN                                          
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION W,WK,WKM,EK                                      
      DOUBLE PRECISION TNORM,YNORM,S,SM,DASUM                           
      INTEGER I1,J,J1,J2,K,KK,L                                         
      LOGICAL LOWER                                                     
!                                                                       
      LOWER = JOB .EQ. 0                                                
!                                                                       
!     COMPUTE 1-NORM OF T                                               
!                                                                       
      TNORM = 0.0D0                                                     
      DO 10 J = 1, N                                                    
         L = J                                                          
         IF (LOWER) L = N + 1 - J                                       
         I1 = 1                                                         
         IF (LOWER) I1 = J                                              
         TNORM = DMAX1(TNORM,DASUM(L,T(I1,J),1))                        
   10 CONTINUE                                                          
!                                                                       
!     RCOND = 1/(NORM(T)*(ESTIMATE OF NORM(INVERSE(T)))) .              
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  T*Z = Y  AND  TRANS(T)*Y = E .  
!     TRANS(T)  IS THE TRANSPOSE OF T .                                 
!     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL           
!     GROWTH IN THE ELEMENTS OF Y .                                     
!     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.            
!                                                                       
!     SOLVE TRANS(T)*Y = E                                              
!                                                                       
      EK = 1.0D0                                                        
      DO 20 J = 1, N                                                    
         Z(J) = 0.0D0                                                   
   20 CONTINUE                                                          
      DO 100 KK = 1, N                                                  
         K = KK                                                         
         IF (LOWER) K = N + 1 - KK                                      
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))                      
         IF (DABS(EK-Z(K)) .LE. DABS(T(K,K))) GO TO 30                  
            S = DABS(T(K,K))/DABS(EK-Z(K))                              
            CALL DSCAL(N,S,Z,1)                                         
            EK = S*EK                                                   
   30    CONTINUE                                                       
         WK = EK - Z(K)                                                 
         WKM = -EK - Z(K)                                               
         S = DABS(WK)                                                   
         SM = DABS(WKM)                                                 
         IF (T(K,K) .EQ. 0.0D0) GO TO 40                                
            WK = WK/T(K,K)                                              
            WKM = WKM/T(K,K)                                            
         GO TO 50                                                       
   40    CONTINUE                                                       
            WK = 1.0D0                                                  
            WKM = 1.0D0                                                 
   50    CONTINUE                                                       
         IF (KK .EQ. N) GO TO 90                                        
            J1 = K + 1                                                  
            IF (LOWER) J1 = 1                                           
            J2 = N                                                      
            IF (LOWER) J2 = K - 1                                       
            DO 60 J = J1, J2                                            
               SM = SM + DABS(Z(J)+WKM*T(K,J))                          
               Z(J) = Z(J) + WK*T(K,J)                                  
               S = S + DABS(Z(J))                                       
   60       CONTINUE                                                    
            IF (S .GE. SM) GO TO 80                                     
               W = WKM - WK                                             
               WK = WKM                                                 
               DO 70 J = J1, J2                                         
                  Z(J) = Z(J) + W*T(K,J)                                
   70          CONTINUE                                                 
   80       CONTINUE                                                    
   90    CONTINUE                                                       
         Z(K) = WK                                                      
  100 CONTINUE                                                          
      S = 1.0D0/DASUM(N,Z,1)                                            
      CALL DSCAL(N,S,Z,1)                                               
!                                                                       
      YNORM = 1.0D0                                                     
!                                                                       
!     SOLVE T*Z = Y                                                     
!                                                                       
      DO 130 KK = 1, N                                                  
         K = N + 1 - KK                                                 
         IF (LOWER) K = KK                                              
         IF (DABS(Z(K)) .LE. DABS(T(K,K))) GO TO 110                    
            S = DABS(T(K,K))/DABS(Z(K))                                 
            CALL DSCAL(N,S,Z,1)                                         
            YNORM = S*YNORM                                             
  110    CONTINUE                                                       
         IF (T(K,K) .NE. 0.0D0) Z(K) = Z(K)/T(K,K)                      
         IF (T(K,K) .EQ. 0.0D0) Z(K) = 1.0D0                            
         I1 = 1                                                         
         IF (LOWER) I1 = K + 1                                          
         IF (KK .GE. N) GO TO 120                                       
            W = -Z(K)                                                   
            CALL DAXPY(N-KK,W,T(I1,K),1,Z(I1),1)                        
  120    CONTINUE                                                       
  130 CONTINUE                                                          
!     MAKE ZNORM = 1.0                                                  
      S = 1.0D0/DASUM(N,Z,1)                                            
      CALL DSCAL(N,S,Z,1)                                               
      YNORM = S*YNORM                                                   
!                                                                       
      IF (TNORM .NE. 0.0D0) RCOND = YNORM/TNORM                         
      IF (TNORM .EQ. 0.0D0) RCOND = 0.0D0                               
      RETURN                                                            
      END                                                               
      SUBROUTINE DTRSL(T,LDT,N,B,JOB,INFO)                              
      INTEGER LDT,N,JOB,INFO                                            
      DOUBLE PRECISION T(LDT,1),B(1)                                    
!                                                                       
!                                                                       
!     DTRSL SOLVES SYSTEMS OF THE FORM                                  
!                                                                       
!                   T * X = B                                           
!     OR                                                                
!                   TRANS(T) * X = B                                    
!                                                                       
!     WHERE T IS A TRIANGULAR MATRIX OF ORDER N. HERE TRANS(T)          
!     DENOTES THE TRANSPOSE OF THE MATRIX T.                            
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!         T         DOUBLE PRECISION(LDT,N)                             
!                   T CONTAINS THE MATRIX OF THE SYSTEM. THE ZERO       
!                   ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND      
!                   THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE      
!                   USED TO STORE OTHER INFORMATION.                    
!                                                                       
!         LDT       INTEGER                                             
!                   LDT IS THE LEADING DIMENSION OF THE ARRAY T.        
!                                                                       
!         N         INTEGER                                             
!                   N IS THE ORDER OF THE SYSTEM.                       
!                                                                       
!         B         DOUBLE PRECISION(N).                                
!                   B CONTAINS THE RIGHT HAND SIDE OF THE SYSTEM.       
!                                                                       
!         JOB       INTEGER                                             
!                   JOB SPECIFIES WHAT KIND OF SYSTEM IS TO BE SOLVED.  
!                   IF JOB IS                                           
!                                                                       
!                        00   SOLVE T*X=B, T LOWER TRIANGULAR,          
!                        01   SOLVE T*X=B, T UPPER TRIANGULAR,          
!                        10   SOLVE TRANS(T)*X=B, T LOWER TRIANGULAR,   
!                        11   SOLVE TRANS(T)*X=B, T UPPER TRIANGULAR.   
!                                                                       
!     ON RETURN                                                         
!                                                                       
!         B         B CONTAINS THE SOLUTION, IF INFO .EQ. 0.            
!                   OTHERWISE B IS UNALTERED.                           
!                                                                       
!         INFO      INTEGER                                             
!                   INFO CONTAINS ZERO IF THE SYSTEM IS NONSINGULAR.    
!                   OTHERWISE INFO CONTAINS THE INDEX OF                
!                   THE FIRST ZERO DIAGONAL ELEMENT OF T.               
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     G. W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DAXPY,DDOT                                                   
!     FORTRAN MOD                                                       
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION DDOT,TEMP                                        
      INTEGER CASE,J,JJ                                                 
!                                                                       
!     BEGIN BLOCK PERMITTING ...EXITS TO 150                            
!                                                                       
!        CHECK FOR ZERO DIAGONAL ELEMENTS.                              
!                                                                       
         DO 10 INFO = 1, N                                              
!     ......EXIT                                                        
            IF (T(INFO,INFO) .EQ. 0.0D0) GO TO 150                      
   10    CONTINUE                                                       
         INFO = 0                                                       
!                                                                       
!        DETERMINE THE TASK AND GO TO IT.                               
!                                                                       
         CASE = 1                                                       
         IF (MOD(JOB,10) .NE. 0) CASE = 2                               
         IF (MOD(JOB,100)/10 .NE. 0) CASE = CASE + 2                    
         select case (CASE)
         case (1)
           go to 20
         case (2)
           go to 50
         case (3)
           go to 80
         case (4)
           go to 110
         end select
!                                                                       
!        SOLVE T*X=B FOR T LOWER TRIANGULAR                             
!                                                                       
   20    CONTINUE                                                       
            B(1) = B(1)/T(1,1)                                          
            IF (N .LT. 2) GO TO 40                                      
            DO 30 J = 2, N                                              
               TEMP = -B(J-1)                                           
               CALL DAXPY(N-J+1,TEMP,T(J,J-1),1,B(J),1)                 
               B(J) = B(J)/T(J,J)                                       
   30       CONTINUE                                                    
   40       CONTINUE                                                    
         GO TO 140                                                      
!                                                                       
!        SOLVE T*X=B FOR T UPPER TRIANGULAR.                            
!                                                                       
   50    CONTINUE                                                       
            B(N) = B(N)/T(N,N)                                          
            IF (N .LT. 2) GO TO 70                                      
            DO 60 JJ = 2, N                                             
               J = N - JJ + 1                                           
               TEMP = -B(J+1)                                           
               CALL DAXPY(J,TEMP,T(1,J+1),1,B(1),1)                     
               B(J) = B(J)/T(J,J)                                       
   60       CONTINUE                                                    
   70       CONTINUE                                                    
         GO TO 140                                                      
!                                                                       
!        SOLVE TRANS(T)*X=B FOR T LOWER TRIANGULAR.                     
!                                                                       
   80    CONTINUE                                                       
            B(N) = B(N)/T(N,N)                                          
            IF (N .LT. 2) GO TO 100                                     
            DO 90 JJ = 2, N                                             
               J = N - JJ + 1                                           
               B(J) = B(J) - DDOT(JJ-1,T(J+1,J),1,B(J+1),1)             
               B(J) = B(J)/T(J,J)                                       
   90       CONTINUE                                                    
  100       CONTINUE                                                    
         GO TO 140                                                      
!                                                                       
!        SOLVE TRANS(T)*X=B FOR T UPPER TRIANGULAR.                     
!                                                                       
  110    CONTINUE                                                       
            B(1) = B(1)/T(1,1)                                          
            IF (N .LT. 2) GO TO 130                                     
            DO 120 J = 2, N                                             
               B(J) = B(J) - DDOT(J-1,T(1,J),1,B(1),1)                  
               B(J) = B(J)/T(J,J)                                       
  120       CONTINUE                                                    
  130       CONTINUE                                                    
  140    CONTINUE                                                       
  150 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DTRDI(T,LDT,N,DET,JOB,INFO)                            
      INTEGER LDT,N,JOB,INFO                                            
      DOUBLE PRECISION T(LDT,1),DET(2)                                  
!                                                                       
!     DTRDI COMPUTES THE DETERMINANT AND INVERSE OF A DOUBLE PRECISION  
!     TRIANGULAR MATRIX.                                                
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        T       DOUBLE PRECISION(LDT,N)                                
!                T CONTAINS THE TRIANGULAR MATRIX. THE ZERO             
!                ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND         
!                THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE         
!                USED TO STORE OTHER INFORMATION.                       
!                                                                       
!        LDT     INTEGER                                                
!                LDT IS THE LEADING DIMENSION OF THE ARRAY T.           
!                                                                       
!        N       INTEGER                                                
!                N IS THE ORDER OF THE SYSTEM.                          
!                                                                       
!        JOB     INTEGER                                                
!                = 010       NO DET, INVERSE OF LOWER TRIANGULAR.       
!                = 011       NO DET, INVERSE OF UPPER TRIANGULAR.       
!                = 100       DET, NO INVERSE.                           
!                = 110       DET, INVERSE OF LOWER TRIANGULAR.          
!                = 111       DET, INVERSE OF UPPER TRIANGULAR.          
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        T       INVERSE OF ORIGINAL MATRIX IF REQUESTED.               
!                OTHERWISE UNCHANGED.                                   
!                                                                       
!        DET     DOUBLE PRECISION(2)                                    
!                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.           
!                OTHERWISE NOT REFERENCED.                              
!                DETERMINANT = DET(1) * 10.0**DET(2)                    
!                WITH  1.0 .LE. DABS(DET(1)) .LT. 10.0                  
!                OR  DET(1) .EQ. 0.0 .                                  
!                                                                       
!        INFO    INTEGER                                                
!                INFO CONTAINS ZERO IF THE SYSTEM IS NONSINGULAR        
!                AND THE INVERSE IS REQUESTED.                          
!                OTHERWISE INFO CONTAINS THE INDEX OF                   
!                A ZERO DIAGONAL ELEMENT OF T.                          
!                                                                       
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DAXPY,DSCAL                                                  
!     FORTRAN DABS,MOD                                                  
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      DOUBLE PRECISION TEMP                                             
      DOUBLE PRECISION TEN                                              
      INTEGER I,J,K,KB,KM1,KP1                                          
!                                                                       
!     BEGIN BLOCK PERMITTING ...EXITS TO 180                            
!                                                                       
!        COMPUTE DETERMINANT                                            
!                                                                       
         IF (JOB/100 .EQ. 0) GO TO 70                                   
            DET(1) = 1.0D0                                              
            DET(2) = 0.0D0                                              
            TEN = 10.0D0                                                
            DO 50 I = 1, N                                              
               DET(1) = T(I,I)*DET(1)                                   
!           ...EXIT                                                     
               IF (DET(1) .EQ. 0.0D0) GO TO 60                          
   10          IF (DABS(DET(1)) .GE. 1.0D0) GO TO 20                    
                  DET(1) = TEN*DET(1)                                   
                  DET(2) = DET(2) - 1.0D0                               
               GO TO 10                                                 
   20          CONTINUE                                                 
   30          IF (DABS(DET(1)) .LT. TEN) GO TO 40                      
                  DET(1) = DET(1)/TEN                                   
                  DET(2) = DET(2) + 1.0D0                               
               GO TO 30                                                 
   40          CONTINUE                                                 
   50       CONTINUE                                                    
   60       CONTINUE                                                    
   70    CONTINUE                                                       
!                                                                       
!        COMPUTE INVERSE OF UPPER TRIANGULAR                            
!                                                                       
         IF (MOD(JOB/10,10) .EQ. 0) GO TO 170                           
            IF (MOD(JOB,10) .EQ. 0) GO TO 120                           
!              BEGIN BLOCK PERMITTING ...EXITS TO 110                   
                  DO 100 K = 1, N                                       
                     INFO = K                                           
!              ......EXIT                                               
                     IF (T(K,K) .EQ. 0.0D0) GO TO 110                   
                     T(K,K) = 1.0D0/T(K,K)                              
                     TEMP = -T(K,K)                                     
                     CALL DSCAL(K-1,TEMP,T(1,K),1)                      
                     KP1 = K + 1                                        
                     IF (N .LT. KP1) GO TO 90                           
                     DO 80 J = KP1, N                                   
                        TEMP = T(K,J)                                   
                        T(K,J) = 0.0D0                                  
                        CALL DAXPY(K,TEMP,T(1,K),1,T(1,J),1)            
   80                CONTINUE                                           
   90                CONTINUE                                           
  100             CONTINUE                                              
                  INFO = 0                                              
  110          CONTINUE                                                 
            GO TO 160                                                   
  120       CONTINUE                                                    
!                                                                       
!              COMPUTE INVERSE OF LOWER TRIANGULAR                      
!                                                                       
               DO 150 KB = 1, N                                         
                  K = N + 1 - KB                                        
                  INFO = K                                              
!     ............EXIT                                                  
                  IF (T(K,K) .EQ. 0.0D0) GO TO 180                      
                  T(K,K) = 1.0D0/T(K,K)                                 
                  TEMP = -T(K,K)                                        
                  IF (K .NE. N) CALL DSCAL(N-K,TEMP,T(K+1,K),1)         
                  KM1 = K - 1                                           
                  IF (KM1 .LT. 1) GO TO 140                             
                  DO 130 J = 1, KM1                                     
                     TEMP = T(K,J)                                      
                     T(K,J) = 0.0D0                                     
                     CALL DAXPY(N-K+1,TEMP,T(K,K),1,T(K,J),1)           
  130             CONTINUE                                              
  140             CONTINUE                                              
  150          CONTINUE                                                 
               INFO = 0                                                 
  160       CONTINUE                                                    
  170    CONTINUE                                                       
  180 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DGTSL(N,C,D,E,B,INFO)                                  
      INTEGER N,INFO                                                    
      DOUBLE PRECISION C(1),D(1),E(1),B(1)                              
!                                                                       
!     DGTSL GIVEN A GENERAL TRIDIAGONAL MATRIX AND A RIGHT HAND         
!     SIDE WILL FIND THE SOLUTION.                                      
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        N       INTEGER                                                
!                IS THE ORDER OF THE TRIDIAGONAL MATRIX.                
!                                                                       
!        C       DOUBLE PRECISION(N)                                    
!                IS THE SUBDIAGONAL OF THE TRIDIAGONAL MATRIX.          
!                C(2) THROUGH C(N) SHOULD CONTAIN THE SUBDIAGONAL.      
!                ON OUTPUT C IS DESTROYED.                              
!                                                                       
!        D       DOUBLE PRECISION(N)                                    
!                IS THE DIAGONAL OF THE TRIDIAGONAL MATRIX.             
!                ON OUTPUT D IS DESTROYED.                              
!                                                                       
!        E       DOUBLE PRECISION(N)                                    
!                IS THE SUPERDIAGONAL OF THE TRIDIAGONAL MATRIX.        
!                E(1) THROUGH E(N-1) SHOULD CONTAIN THE SUPERDIAGONAL.  
!                ON OUTPUT E IS DESTROYED.                              
!                                                                       
!        B       DOUBLE PRECISION(N)                                    
!                IS THE RIGHT HAND SIDE VECTOR.                         
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        B       IS THE SOLUTION VECTOR.                                
!                                                                       
!        INFO    INTEGER                                                
!                = 0 NORMAL VALUE.                                      
!                = K IF THE K-TH ELEMENT OF THE DIAGONAL BECOMES        
!                    EXACTLY ZERO.  THE SUBROUTINE RETURNS WHEN         
!                    THIS IS DETECTED.                                  
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     JACK DONGARRA, ARGONNE NATIONAL LABORATORY.                       
!                                                                       
!     NO EXTERNALS                                                      
!     FORTRAN DABS                                                      
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      INTEGER K,KB,KP1,NM1,NM2                                          
      DOUBLE PRECISION T                                                
!     BEGIN BLOCK PERMITTING ...EXITS TO 100                            
!                                                                       
         INFO = 0                                                       
         C(1) = D(1)                                                    
         NM1 = N - 1                                                    
         IF (NM1 .LT. 1) GO TO 40                                       
            D(1) = E(1)                                                 
            E(1) = 0.0D0                                                
            E(N) = 0.0D0                                                
!                                                                       
            DO 30 K = 1, NM1                                            
               KP1 = K + 1                                              
!                                                                       
!              FIND THE LARGEST OF THE TWO ROWS                         
!                                                                       
               IF (DABS(C(KP1)) .LT. DABS(C(K))) GO TO 10               
!                                                                       
!                 INTERCHANGE ROW                                       
!                                                                       
                  T = C(KP1)                                            
                  C(KP1) = C(K)                                         
                  C(K) = T                                              
                  T = D(KP1)                                            
                  D(KP1) = D(K)                                         
                  D(K) = T                                              
                  T = E(KP1)                                            
                  E(KP1) = E(K)                                         
                  E(K) = T                                              
                  T = B(KP1)                                            
                  B(KP1) = B(K)                                         
                  B(K) = T                                              
   10          CONTINUE                                                 
!                                                                       
!              ZERO ELEMENTS                                            
!                                                                       
               IF (C(K) .NE. 0.0D0) GO TO 20                            
                  INFO = K                                              
!     ............EXIT                                                  
                  GO TO 100                                             
   20          CONTINUE                                                 
               T = -C(KP1)/C(K)                                         
               C(KP1) = D(KP1) + T*D(K)                                 
               D(KP1) = E(KP1) + T*E(K)                                 
               E(KP1) = 0.0D0                                           
               B(KP1) = B(KP1) + T*B(K)                                 
   30       CONTINUE                                                    
   40    CONTINUE                                                       
         IF (C(N) .NE. 0.0D0) GO TO 50                                  
            INFO = N                                                    
         GO TO 90                                                       
   50    CONTINUE                                                       
!                                                                       
!           BACK SOLVE                                                  
!                                                                       
            NM2 = N - 2                                                 
            B(N) = B(N)/C(N)                                            
            IF (N .EQ. 1) GO TO 80                                      
               B(NM1) = (B(NM1) - D(NM1)*B(N))/C(NM1)                   
               IF (NM2 .LT. 1) GO TO 70                                 
               DO 60 KB = 1, NM2                                        
                  K = NM2 - KB + 1                                      
                  B(K) = (B(K) - D(K)*B(K+1) - E(K)*B(K+2))/C(K)        
   60          CONTINUE                                                 
   70          CONTINUE                                                 
   80       CONTINUE                                                    
   90    CONTINUE                                                       
  100 CONTINUE                                                          
!                                                                       
      RETURN                                                            
      END                                                               
      SUBROUTINE DPTSL(N,D,E,B)                                         
      INTEGER N                                                         
      DOUBLE PRECISION D(1),E(1),B(1)                                   
!                                                                       
!     DPTSL GIVEN A POSITIVE DEFINITE TRIDIAGONAL MATRIX AND A RIGHT    
!     HAND SIDE WILL FIND THE SOLUTION.                                 
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        N        INTEGER                                               
!                 IS THE ORDER OF THE TRIDIAGONAL MATRIX.               
!                                                                       
!        D        DOUBLE PRECISION(N)                                   
!                 IS THE DIAGONAL OF THE TRIDIAGONAL MATRIX.            
!                 ON OUTPUT D IS DESTROYED.                             
!                                                                       
!        E        DOUBLE PRECISION(N)                                   
!                 IS THE OFFDIAGONAL OF THE TRIDIAGONAL MATRIX.         
!                 E(1) THROUGH E(N-1) SHOULD CONTAIN THE                
!                 OFFDIAGONAL.                                          
!                                                                       
!        B        DOUBLE PRECISION(N)                                   
!                 IS THE RIGHT HAND SIDE VECTOR.                        
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        B        CONTAINS THE SOULTION.                                
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     JACK DONGARRA, ARGONNE NATIONAL LABORATORY.                       
!                                                                       
!     NO EXTERNALS                                                      
!     FORTRAN MOD                                                       
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      INTEGER K,KBM1,KE,KF,KP1,NM1,NM1D2                                
      DOUBLE PRECISION T1,T2                                            
!                                                                       
!     CHECK FOR 1 X 1 CASE                                              
!                                                                       
      IF (N .NE. 1) GO TO 10                                            
         B(1) = B(1)/D(1)                                               
      GO TO 70                                                          
   10 CONTINUE                                                          
         NM1 = N - 1                                                    
         NM1D2 = NM1/2                                                  
         IF (N .EQ. 2) GO TO 30                                         
            KBM1 = N - 1                                                
!                                                                       
!           ZERO TOP HALF OF SUBDIAGONAL AND BOTTOM HALF OF             
!           SUPERDIAGONAL                                               
!                                                                       
            DO 20 K = 1, NM1D2                                          
               T1 = E(K)/D(K)                                           
               D(K+1) = D(K+1) - T1*E(K)                                
               B(K+1) = B(K+1) - T1*B(K)                                
               T2 = E(KBM1)/D(KBM1+1)                                   
               D(KBM1) = D(KBM1) - T2*E(KBM1)                           
               B(KBM1) = B(KBM1) - T2*B(KBM1+1)                         
               KBM1 = KBM1 - 1                                          
   20       CONTINUE                                                    
   30    CONTINUE                                                       
         KP1 = NM1D2 + 1                                                
!                                                                       
!        CLEAN UP FOR POSSIBLE 2 X 2 BLOCK AT CENTER                    
!                                                                       
         IF (MOD(N,2) .NE. 0) GO TO 40                                  
            T1 = E(KP1)/D(KP1)                                          
            D(KP1+1) = D(KP1+1) - T1*E(KP1)                             
            B(KP1+1) = B(KP1+1) - T1*B(KP1)                             
            KP1 = KP1 + 1                                               
   40    CONTINUE                                                       
!                                                                       
!        BACK SOLVE STARTING AT THE CENTER, GOING TOWARDS THE TOP       
!        AND BOTTOM                                                     
!                                                                       
         B(KP1) = B(KP1)/D(KP1)                                         
         IF (N .EQ. 2) GO TO 60                                         
            K = KP1 - 1                                                 
            KE = KP1 + NM1D2 - 1                                        
            DO 50 KF = KP1, KE                                          
               B(K) = (B(K) - E(K)*B(K+1))/D(K)                         
               B(KF+1) = (B(KF+1) - E(KF)*B(KF))/D(KF+1)                
               K = K - 1                                                
   50       CONTINUE                                                    
   60    CONTINUE                                                       
         IF (MOD(N,2) .EQ. 0) B(1) = (B(1) - E(1)*B(2))/D(1)            
   70 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DCHDC(A,LDA,P,WORK,JPVT,JOB,INFO)                      
      INTEGER LDA,P,JPVT(1),JOB,INFO                                    
      DOUBLE PRECISION A(LDA,1),WORK(1)                                 
!                                                                       
!     DCHDC COMPUTES THE CHOLESKY DECOMPOSITION OF A POSITIVE DEFINITE  
!     MATRIX.  A PIVOTING OPTION ALLOWS THE USER TO ESTIMATE THE        
!     CONDITION OF A POSITIVE DEFINITE MATRIX OR DETERMINE THE RANK     
!     OF A POSITIVE SEMIDEFINITE MATRIX.                                
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!         A      DOUBLE PRECISION(LDA,P).                               
!                A CONTAINS THE MATRIX WHOSE DECOMPOSITION IS TO        
!                BE COMPUTED.  ONLT THE UPPER HALF OF A NEED BE STORED. 
!                THE LOWER PART OF THE ARRAY A IS NOT REFERENCED.       
!                                                                       
!         LDA    INTEGER.                                               
!                LDA IS THE LEADING DIMENSION OF THE ARRAY A.           
!                                                                       
!         P      INTEGER.                                               
!                P IS THE ORDER OF THE MATRIX.                          
!                                                                       
!         WORK   DOUBLE PRECISION.                                      
!                WORK IS A WORK ARRAY.                                  
!                                                                       
!         JPVT   INTEGER(P).                                            
!                JPVT CONTAINS INTEGERS THAT CONTROL THE SELECTION      
!                OF THE PIVOT ELEMENTS, IF PIVOTING HAS BEEN REQUESTED. 
!                EACH DIAGONAL ELEMENT A(K,K)                           
!                IS PLACED IN ONE OF THREE CLASSES ACCORDING TO THE     
!                VALUE OF JPVT(K).                                      
!                                                                       
!                   IF JPVT(K) .GT. 0, THEN X(K) IS AN INITIAL          
!                                      ELEMENT.                         
!                                                                       
!                   IF JPVT(K) .EQ. 0, THEN X(K) IS A FREE ELEMENT.     
!                                                                       
!                   IF JPVT(K) .LT. 0, THEN X(K) IS A FINAL ELEMENT.    
!                                                                       
!                BEFORE THE DECOMPOSITION IS COMPUTED, INITIAL ELEMENTS 
!                ARE MOVED BY SYMMETRIC ROW AND COLUMN INTERCHANGES TO  
!                THE BEGINNING OF THE ARRAY A AND FINAL                 
!                ELEMENTS TO THE END.  BOTH INITIAL AND FINAL ELEMENTS  
!                ARE FROZEN IN PLACE DURING THE COMPUTATION AND ONLY    
!                FREE ELEMENTS ARE MOVED.  AT THE K-TH STAGE OF THE     
!                REDUCTION, IF A(K,K) IS OCCUPIED BY A FREE ELEMENT     
!                IT IS INTERCHANGED WITH THE LARGEST FREE ELEMENT       
!                A(L,L) WITH L .GE. K.  JPVT IS NOT REFERENCED IF       
!                JOB .EQ. 0.                                            
!                                                                       
!        JOB     INTEGER.                                               
!                JOB IS AN INTEGER THAT INITIATES COLUMN PIVOTING.      
!                IF JOB .EQ. 0, NO PIVOTING IS DONE.                    
!                IF JOB .NE. 0, PIVOTING IS DONE.                       
!                                                                       
!     ON RETURN                                                         
!                                                                       
!         A      A CONTAINS IN ITS UPPER HALF THE CHOLESKY FACTOR       
!                OF THE MATRIX A AS IT HAS BEEN PERMUTED BY PIVOTING.   
!                                                                       
!         JPVT   JPVT(J) CONTAINS THE INDEX OF THE DIAGONAL ELEMENT     
!                OF A THAT WAS MOVED INTO THE J-TH POSITION,            
!                PROVIDED PIVOTING WAS REQUESTED.                       
!                                                                       
!         INFO   CONTAINS THE INDEX OF THE LAST POSITIVE DIAGONAL       
!                ELEMENT OF THE CHOLESKY FACTOR.                        
!                                                                       
!     FOR POSITIVE DEFINITE MATRICES INFO = P IS THE NORMAL RETURN.     
!     FOR PIVOTING WITH POSITIVE SEMIDEFINITE MATRICES INFO WILL        
!     IN GENERAL BE LESS THAN P.  HOWEVER, INFO MAY BE GREATER THAN     
!     THE RANK OF A, SINCE ROUNDING ERROR CAN CAUSE AN OTHERWISE ZERO   
!     ELEMENT TO BE POSITIVE. INDEFINITE SYSTEMS WILL ALWAYS CAUSE      
!     INFO TO BE LESS THAN P.                                           
!                                                                       
!     LINPACK. THIS VERSION DATED 03/19/79 .                            
!     J.J. DONGARRA AND G.W. STEWART, ARGONNE NATIONAL LABORATORY AND   
!     UNIVERSITY OF MARYLAND.                                           
!                                                                       
!                                                                       
!     BLAS DAXPY,DSWAP                                                  
!     FORTRAN DSQRT                                                     
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      INTEGER PU,PL,PLP1,I,J,JP,JT,K,KB,KM1,KP1,L,MAXL                  
      DOUBLE PRECISION TEMP                                             
      DOUBLE PRECISION MAXDIA                                           
      LOGICAL SWAPK,NEGK                                                
!                                                                       
      PL = 1                                                            
      PU = 0                                                            
      INFO = P                                                          
      IF (JOB .EQ. 0) GO TO 160                                         
!                                                                       
!        PIVOTING HAS BEEN REQUESTED. REARRANGE THE                     
!        THE ELEMENTS ACCORDING TO JPVT.                                
!                                                                       
         DO 70 K = 1, P                                                 
            SWAPK = JPVT(K) .GT. 0                                      
            NEGK = JPVT(K) .LT. 0                                       
            JPVT(K) = K                                                 
            IF (NEGK) JPVT(K) = -JPVT(K)                                
            IF (.NOT.SWAPK) GO TO 60                                    
               IF (K .EQ. PL) GO TO 50                                  
                  CALL DSWAP(PL-1,A(1,K),1,A(1,PL),1)                   
                  TEMP = A(K,K)                                         
                  A(K,K) = A(PL,PL)                                     
                  A(PL,PL) = TEMP                                       
                  PLP1 = PL + 1                                         
                  IF (P .LT. PLP1) GO TO 40                             
                  DO 30 J = PLP1, P                                     
                     IF (J .GE. K) GO TO 10                             
                        TEMP = A(PL,J)                                  
                        A(PL,J) = A(J,K)                                
                        A(J,K) = TEMP                                   
                     GO TO 20                                           
   10                CONTINUE                                           
                     IF (J .EQ. K) GO TO 20                             
                        TEMP = A(K,J)                                   
                        A(K,J) = A(PL,J)                                
                        A(PL,J) = TEMP                                  
   20                CONTINUE                                           
   30             CONTINUE                                              
   40             CONTINUE                                              
                  JPVT(K) = JPVT(PL)                                    
                  JPVT(PL) = K                                          
   50          CONTINUE                                                 
               PL = PL + 1                                              
   60       CONTINUE                                                    
   70    CONTINUE                                                       
         PU = P                                                         
         IF (P .LT. PL) GO TO 150                                       
         DO 140 KB = PL, P                                              
            K = P - KB + PL                                             
            IF (JPVT(K) .GE. 0) GO TO 130                               
               JPVT(K) = -JPVT(K)                                       
               IF (PU .EQ. K) GO TO 120                                 
                  CALL DSWAP(K-1,A(1,K),1,A(1,PU),1)                    
                  TEMP = A(K,K)                                         
                  A(K,K) = A(PU,PU)                                     
                  A(PU,PU) = TEMP                                       
                  KP1 = K + 1                                           
                  IF (P .LT. KP1) GO TO 110                             
                  DO 100 J = KP1, P                                     
                     IF (J .GE. PU) GO TO 80                            
                        TEMP = A(K,J)                                   
                        A(K,J) = A(J,PU)                                
                        A(J,PU) = TEMP                                  
                     GO TO 90                                           
   80                CONTINUE                                           
                     IF (J .EQ. PU) GO TO 90                            
                        TEMP = A(K,J)                                   
                        A(K,J) = A(PU,J)                                
                        A(PU,J) = TEMP                                  
   90                CONTINUE                                           
  100             CONTINUE                                              
  110             CONTINUE                                              
                  JT = JPVT(K)                                          
                  JPVT(K) = JPVT(PU)                                    
                  JPVT(PU) = JT                                         
  120          CONTINUE                                                 
               PU = PU - 1                                              
  130       CONTINUE                                                    
  140    CONTINUE                                                       
  150    CONTINUE                                                       
  160 CONTINUE                                                          
      DO 270 K = 1, P                                                   
!                                                                       
!        REDUCTION LOOP.                                                
!                                                                       
         MAXDIA = A(K,K)                                                
         KP1 = K + 1                                                    
         MAXL = K                                                       
!                                                                       
!        DETERMINE THE PIVOT ELEMENT.                                   
!                                                                       
         IF (K .LT. PL .OR. K .GE. PU) GO TO 190                        
            DO 180 L = KP1, PU                                          
               IF (A(L,L) .LE. MAXDIA) GO TO 170                        
                  MAXDIA = A(L,L)                                       
                  MAXL = L                                              
  170          CONTINUE                                                 
  180       CONTINUE                                                    
  190    CONTINUE                                                       
!                                                                       
!        QUIT IF THE PIVOT ELEMENT IS NOT POSITIVE.                     
!                                                                       
         IF (MAXDIA .GT. 0.0D0) GO TO 200                               
            INFO = K - 1                                                
!     ......EXIT                                                        
            GO TO 280                                                   
  200    CONTINUE                                                       
         IF (K .EQ. MAXL) GO TO 210                                     
!                                                                       
!           START THE PIVOTING AND UPDATE JPVT.                         
!                                                                       
            KM1 = K - 1                                                 
            CALL DSWAP(KM1,A(1,K),1,A(1,MAXL),1)                        
            A(MAXL,MAXL) = A(K,K)                                       
            A(K,K) = MAXDIA                                             
            JP = JPVT(MAXL)                                             
            JPVT(MAXL) = JPVT(K)                                        
            JPVT(K) = JP                                                
  210    CONTINUE                                                       
!                                                                       
!        REDUCTION STEP. PIVOTING IS CONTAINED ACROSS THE ROWS.         
!                                                                       
         WORK(K) = DSQRT(A(K,K))                                        
         A(K,K) = WORK(K)                                               
         IF (P .LT. KP1) GO TO 260                                      
         DO 250 J = KP1, P                                              
            IF (K .EQ. MAXL) GO TO 240                                  
               IF (J .GE. MAXL) GO TO 220                               
                  TEMP = A(K,J)                                         
                  A(K,J) = A(J,MAXL)                                    
                  A(J,MAXL) = TEMP                                      
               GO TO 230                                                
  220          CONTINUE                                                 
               IF (J .EQ. MAXL) GO TO 230                               
                  TEMP = A(K,J)                                         
                  A(K,J) = A(MAXL,J)                                    
                  A(MAXL,J) = TEMP                                      
  230          CONTINUE                                                 
  240       CONTINUE                                                    
            A(K,J) = A(K,J)/WORK(K)                                     
            WORK(J) = A(K,J)                                            
            TEMP = -A(K,J)                                              
            CALL DAXPY(J-K,TEMP,WORK(KP1),1,A(KP1,J),1)                 
  250    CONTINUE                                                       
  260    CONTINUE                                                       
  270 CONTINUE                                                          
  280 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DCHUD(R,LDR,P,X,Z,LDZ,NZ,Y,RHO,C,S)                    
      INTEGER LDR,P,LDZ,NZ                                              
      DOUBLE PRECISION RHO(1),C(1)                                      
      DOUBLE PRECISION R(LDR,1),X(1),Z(LDZ,1),Y(1),S(1)                 
!                                                                       
!     DCHUD UPDATES AN AUGMENTED CHOLESKY DECOMPOSITION OF THE          
!     TRIANGULAR PART OF AN AUGMENTED QR DECOMPOSITION.  SPECIFICALLY,  
!     GIVEN AN UPPER TRIANGULAR MATRIX R OF ORDER P, A ROW VECTOR       
!     X, A COLUMN VECTOR Z, AND A SCALAR Y, DCHUD DETERMINES A          
!     UNTIARY MATRIX U AND A SCALAR ZETA SUCH THAT                      
!                                                                       
!                                                                       
!                              (R  Z)     (RR   ZZ )                    
!                         U  * (    )  =  (        ) ,                  
!                              (X  Y)     ( 0  ZETA)                    
!                                                                       
!     WHERE RR IS UPPER TRIANGULAR.  IF R AND Z HAVE BEEN               
!     OBTAINED FROM THE FACTORIZATION OF A LEAST SQUARES                
!     PROBLEM, THEN RR AND ZZ ARE THE FACTORS CORRESPONDING TO          
!     THE PROBLEM WITH THE OBSERVATION (X,Y) APPENDED.  IN THIS         
!     CASE, IF RHO IS THE NORM OF THE RESIDUAL VECTOR, THEN THE         
!     NORM OF THE RESIDUAL VECTOR OF THE UPDATED PROBLEM IS             
!     DSQRT(RHO**2 + ZETA**2).  DCHUD WILL SIMULTANEOUSLY UPDATE        
!     SEVERAL TRIPLETS (Z,Y,RHO).                                       
!     FOR A LESS TERSE DESCRIPTION OF WHAT DCHUD DOES AND HOW           
!     IT MAY BE APPLIED, SEE THE LINPACK GUIDE.                         
!                                                                       
!     THE MATRIX U IS DETERMINED AS THE PRODUCT U(P)*...*U(1),          
!     WHERE U(I) IS A ROTATION IN THE (I,P+1) PLANE OF THE              
!     FORM                                                              
!                                                                       
!                       (     C(I)      S(I) )                          
!                       (                    ) .                        
!                       (    -S(I)      C(I) )                          
!                                                                       
!     THE ROTATIONS ARE CHOSEN SO THAT C(I) IS DOUBLE PRECISION.        
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!         R      DOUBLE PRECISION(LDR,P), WHERE LDR .GE. P.             
!                R CONTAINS THE UPPER TRIANGULAR MATRIX                 
!                THAT IS TO BE UPDATED.  THE PART OF R                  
!                BELOW THE DIAGONAL IS NOT REFERENCED.                  
!                                                                       
!         LDR    INTEGER.                                               
!                LDR IS THE LEADING DIMENSION OF THE ARRAY R.           
!                                                                       
!         P      INTEGER.                                               
!                P IS THE ORDER OF THE MATRIX R.                        
!                                                                       
!         X      DOUBLE PRECISION(P).                                   
!                X CONTAINS THE ROW TO BE ADDED TO R.  X IS             
!                NOT ALTERED BY DCHUD.                                  
!                                                                       
!         Z      DOUBLE PRECISION(LDZ,NZ), WHERE LDZ .GE. P.            
!                Z IS AN ARRAY CONTAINING NZ P-VECTORS TO               
!                BE UPDATED WITH R.                                     
!                                                                       
!         LDZ    INTEGER.                                               
!                LDZ IS THE LEADING DIMENSION OF THE ARRAY Z.           
!                                                                       
!         NZ     INTEGER.                                               
!                NZ IS THE NUMBER OF VECTORS TO BE UPDATED              
!                NZ MAY BE ZERO, IN WHICH CASE Z, Y, AND RHO            
!                ARE NOT REFERENCED.                                    
!                                                                       
!         Y      DOUBLE PRECISION(NZ).                                  
!                Y CONTAINS THE SCALARS FOR UPDATING THE VECTORS        
!                Z.  Y IS NOT ALTERED BY DCHUD.                         
!                                                                       
!         RHO    DOUBLE PRECISION(NZ).                                  
!                RHO CONTAINS THE NORMS OF THE RESIDUAL                 
!                VECTORS THAT ARE TO BE UPDATED.  IF RHO(J)             
!                IS NEGATIVE, IT IS LEFT UNALTERED.                     
!                                                                       
!     ON RETURN                                                         
!                                                                       
!         RC                                                            
!         RHO    CONTAIN THE UPDATED QUANTITIES.                        
!         Z                                                             
!                                                                       
!         C      DOUBLE PRECISION(P).                                   
!                C CONTAINS THE COSINES OF THE TRANSFORMING             
!                ROTATIONS.                                             
!                                                                       
!         S      DOUBLE PRECISION(P).                                   
!                S CONTAINS THE SINES OF THE TRANSFORMING               
!                ROTATIONS.                                             
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.       
!                                                                       
!     DCHUD USES THE FOLLOWING FUNCTIONS AND SUBROUTINES.               
!                                                                       
!     EXTENDED BLAS DROTG                                               
!     FORTRAN DSQRT                                                     
!                                                                       
      INTEGER I,J,JM1                                                   
      DOUBLE PRECISION AZETA,SCALE                                      
      DOUBLE PRECISION T,XJ,ZETA                                        
!                                                                       
!     UPDATE R.                                                         
!                                                                       
      DO 30 J = 1, P                                                    
         XJ = X(J)                                                      
!                                                                       
!        APPLY THE PREVIOUS ROTATIONS.                                  
!                                                                       
         JM1 = J - 1                                                    
         IF (JM1 .LT. 1) GO TO 20                                       
         DO 10 I = 1, JM1                                               
            T = C(I)*R(I,J) + S(I)*XJ                                   
            XJ = C(I)*XJ - S(I)*R(I,J)                                  
            R(I,J) = T                                                  
   10    CONTINUE                                                       
   20    CONTINUE                                                       
!                                                                       
!        COMPUTE THE NEXT ROTATION.                                     
!                                                                       
         CALL DROTG(R(J,J),XJ,C(J),S(J))                                
   30 CONTINUE                                                          
!                                                                       
!     IF REQUIRED, UPDATE Z AND RHO.                                    
!                                                                       
      IF (NZ .LT. 1) GO TO 70                                           
      DO 60 J = 1, NZ                                                   
         ZETA = Y(J)                                                    
         DO 40 I = 1, P                                                 
            T = C(I)*Z(I,J) + S(I)*ZETA                                 
            ZETA = C(I)*ZETA - S(I)*Z(I,J)                              
            Z(I,J) = T                                                  
   40    CONTINUE                                                       
         AZETA = DABS(ZETA)                                             
         IF (AZETA .EQ. 0.0D0 .OR. RHO(J) .LT. 0.0D0) GO TO 50          
            SCALE = AZETA + RHO(J)                                      
            RHO(J) = SCALE*DSQRT((AZETA/SCALE)**2+(RHO(J)/SCALE)**2)    
   50    CONTINUE                                                       
   60 CONTINUE                                                          
   70 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DCHDD(R,LDR,P,X,Z,LDZ,NZ,Y,RHO,C,S,INFO)               
      INTEGER LDR,P,LDZ,NZ,INFO                                         
      DOUBLE PRECISION R(LDR,1),X(1),Z(LDZ,1),Y(1),S(1)                 
      DOUBLE PRECISION RHO(1),C(1)                                      
!                                                                       
!     DCHDD DOWNDATES AN AUGMENTED CHOLESKY DECOMPOSITION OR THE        
!     TRIANGULAR FACTOR OF AN AUGMENTED QR DECOMPOSITION.               
!     SPECIFICALLY, GIVEN AN UPPER TRIANGULAR MATRIX R OF ORDER P,  A   
!     ROW VECTOR X, A COLUMN VECTOR Z, AND A SCALAR Y, DCHDD            
!     DETERMINEDS A ORTHOGONAL MATRIX U AND A SCALAR ZETA SUCH THAT     
!                                                                       
!                        (R   Z )     (RR  ZZ)                          
!                    U * (      )  =  (      ) ,                        
!                        (0 ZETA)     ( X   Y)                          
!                                                                       
!     WHERE RR IS UPPER TRIANGULAR.  IF R AND Z HAVE BEEN OBTAINED      
!     FROM THE FACTORIZATION OF A LEAST SQUARES PROBLEM, THEN           
!     RR AND ZZ ARE THE FACTORS CORRESPONDING TO THE PROBLEM            
!     WITH THE OBSERVATION (X,Y) REMOVED.  IN THIS CASE, IF RHO         
!     IS THE NORM OF THE RESIDUAL VECTOR, THEN THE NORM OF              
!     THE RESIDUAL VECTOR OF THE DOWNDATED PROBLEM IS                   
!     DSQRT(RHO**2 - ZETA**2). DCHDD WILL SIMULTANEOUSLY DOWNDATE       
!     SEVERAL TRIPLETS (Z,Y,RHO) ALONG WITH R.                          
!     FOR A LESS TERSE DESCRIPTION OF WHAT DCHDD DOES AND HOW           
!     IT MAY BE APPLIED, SEE THE LINPACK GUIDE.                         
!                                                                       
!     THE MATRIX U IS DETERMINED AS THE PRODUCT U(1)*...*U(P)           
!     WHERE U(I) IS A ROTATION IN THE (P+1,I)-PLANE OF THE              
!     FORM                                                              
!                                                                       
!                       ( C(I)     -S(I)     )                          
!                       (                    ) .                        
!                       ( S(I)       C(I)    )                          
!                                                                       
!     THE ROTATIONS ARE CHOSEN SO THAT C(I) IS DOUBLE PRECISION.        
!                                                                       
!     THE USER IS WARNED THAT A GIVEN DOWNDATING PROBLEM MAY            
!     BE IMPOSSIBLE TO ACCOMPLISH OR MAY PRODUCE                        
!     INACCURATE RESULTS.  FOR EXAMPLE, THIS CAN HAPPEN                 
!     IF X IS NEAR A VECTOR WHOSE REMOVAL WILL REDUCE THE               
!     RANK OF R.  BEWARE.                                               
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!         R      DOUBLE PRECISION(LDR,P), WHERE LDR .GE. P.             
!                R CONTAINS THE UPPER TRIANGULAR MATRIX                 
!                THAT IS TO BE DOWNDATED.  THE PART OF  R               
!                BELOW THE DIAGONAL IS NOT REFERENCED.                  
!                                                                       
!         LDR    INTEGER.                                               
!                LDR IS THE LEADING DIMENSION FO THE ARRAY R.           
!                                                                       
!         P      INTEGER.                                               
!                P IS THE ORDER OF THE MATRIX R.                        
!                                                                       
!         X      DOUBLE PRECISION(P).                                   
!                X CONTAINS THE ROW VECTOR THAT IS TO                   
!                BE REMOVED FROM R.  X IS NOT ALTERED BY DCHDD.         
!                                                                       
!         Z      DOUBLE PRECISION(LDZ,NZ), WHERE LDZ .GE. P.            
!                Z IS AN ARRAY OF NZ P-VECTORS WHICH                    
!                ARE TO BE DOWNDATED ALONG WITH R.                      
!                                                                       
!         LDZ    INTEGER.                                               
!                LDZ IS THE LEADING DIMENSION OF THE ARRAY Z.           
!                                                                       
!         NZ     INTEGER.                                               
!                NZ IS THE NUMBER OF VECTORS TO BE DOWNDATED            
!                NZ MAY BE ZERO, IN WHICH CASE Z, Y, AND RHO            
!                ARE NOT REFERENCED.                                    
!                                                                       
!         Y      DOUBLE PRECISION(NZ).                                  
!                Y CONTAINS THE SCALARS FOR THE DOWNDATING              
!                OF THE VECTORS Z.  Y IS NOT ALTERED BY DCHDD.          
!                                                                       
!         RHO    DOUBLE PRECISION(NZ).                                  
!                RHO CONTAINS THE NORMS OF THE RESIDUAL                 
!                VECTORS THAT ARE TO BE DOWNDATED.                      
!                                                                       
!     ON RETURN                                                         
!                                                                       
!         R                                                             
!         Z      CONTAIN THE DOWNDATED QUANTITIES.                      
!         RHO                                                           
!                                                                       
!         C      DOUBLE PRECISION(P).                                   
!                C CONTAINS THE COSINES OF THE TRANSFORMING             
!                ROTATIONS.                                             
!                                                                       
!         S      DOUBLE PRECISION(P).                                   
!                S CONTAINS THE SINES OF THE TRANSFORMING               
!                ROTATIONS.                                             
!                                                                       
!         INFO   INTEGER.                                               
!                INFO IS SET AS FOLLOWS.                                
!                                                                       
!                   INFO = 0  IF THE ENTIRE DOWNDATING                  
!                             WAS SUCCESSFUL.                           
!                                                                       
!                   INFO =-1  IF R COULD NOT BE DOWNDATED.              
!                             IN THIS CASE, ALL QUANTITIES              
!                             ARE LEFT UNALTERED.                       
!                                                                       
!                   INFO = 1  IF SOME RHO COULD NOT BE                  
!                             DOWNDATED.  THE OFFENDING RHOS ARE        
!                             SET TO -1.                                
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.       
!                                                                       
!     DCHDD USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.               
!                                                                       
!     FORTRAN DABS                                                      
!     BLAS DDOT, DNRM2                                                  
!                                                                       
      INTEGER I,II,J                                                    
      DOUBLE PRECISION A,ALPHA,AZETA,NORM,DNRM2                         
      DOUBLE PRECISION DDOT,T,ZETA,B,XX                                 
!                                                                       
!     SOLVE THE SYSTEM TRANS(R)*A = X, PLACING THE RESULT               
!     IN THE ARRAY S.                                                   
!                                                                       
      INFO = 0                                                          
      S(1) = X(1)/R(1,1)                                                
      IF (P .LT. 2) GO TO 20                                            
      DO 10 J = 2, P                                                    
         S(J) = X(J) - DDOT(J-1,R(1,J),1,S,1)                           
         S(J) = S(J)/R(J,J)                                             
   10 CONTINUE                                                          
   20 CONTINUE                                                          
      NORM = DNRM2(P,S,1)                                               
      IF (NORM .LT. 1.0D0) GO TO 30                                     
         INFO = -1                                                      
      GO TO 120                                                         
   30 CONTINUE                                                          
         ALPHA = DSQRT(1.0D0-NORM**2)                                   
!                                                                       
!        DETERMINE THE TRANSFORMATIONS.                                 
!                                                                       
         DO 40 II = 1, P                                                
            I = P - II + 1                                              
            SCALE = ALPHA + DABS(S(I))                                  
            A = ALPHA/SCALE                                             
            B = S(I)/SCALE                                              
            NORM = DSQRT(A**2+B**2+0.0D0**2)                            
            C(I) = A/NORM                                               
            S(I) = B/NORM                                               
            ALPHA = SCALE*NORM                                          
   40    CONTINUE                                                       
!                                                                       
!        APPLY THE TRANSFORMATIONS TO R.                                
!                                                                       
         DO 60 J = 1, P                                                 
            XX = 0.0D0                                                  
            DO 50 II = 1, J                                             
               I = J - II + 1                                           
               T = C(I)*XX + S(I)*R(I,J)                                
               R(I,J) = C(I)*R(I,J) - S(I)*XX                           
               XX = T                                                   
   50       CONTINUE                                                    
   60    CONTINUE                                                       
!                                                                       
!        IF REQUIRED, DOWNDATE Z AND RHO.                               
!                                                                       
         IF (NZ .LT. 1) GO TO 110                                       
         DO 100 J = 1, NZ                                               
            ZETA = Y(J)                                                 
            DO 70 I = 1, P                                              
               Z(I,J) = (Z(I,J) - S(I)*ZETA)/C(I)                       
               ZETA = C(I)*ZETA - S(I)*Z(I,J)                           
   70       CONTINUE                                                    
            AZETA = DABS(ZETA)                                          
            IF (AZETA .LE. RHO(J)) GO TO 80                             
               INFO = 1                                                 
               RHO(J) = -1.0D0                                          
            GO TO 90                                                    
   80       CONTINUE                                                    
               RHO(J) = RHO(J)*DSQRT(1.0D0-(AZETA/RHO(J))**2)           
   90       CONTINUE                                                    
  100    CONTINUE                                                       
  110    CONTINUE                                                       
  120 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DCHEX(R,LDR,P,K,L,Z,LDZ,NZ,C,S,JOB)                    
      INTEGER LDR,P,K,L,LDZ,NZ,JOB                                      
      DOUBLE PRECISION R(LDR,1),Z(LDZ,1),S(1)                           
      DOUBLE PRECISION C(1)                                             
!                                                                       
!     DCHEX UPDATES THE CHOLESKY FACTORIZATION                          
!                                                                       
!                   A = TRANS(R)*R                                      
!                                                                       
!     OF A POSITIVE DEFINITE MATRIX A OF ORDER P UNDER DIAGONAL         
!     PERMUTATIONS OF THE FORM                                          
!                                                                       
!                   TRANS(E)*A*E                                        
!                                                                       
!     WHERE E IS A PERMUTATION MATRIX.  SPECIFICALLY, GIVEN             
!     AN UPPER TRIANGULAR MATRIX R AND A PERMUTATION MATRIX             
!     E (WHICH IS SPECIFIED BY K, L, AND JOB), DCHEX DETERMINES         
!     A ORTHOGONAL MATRIX U SUCH THAT                                   
!                                                                       
!                           U*R*E = RR,                                 
!                                                                       
!     WHERE RR IS UPPER TRIANGULAR.  AT THE USERS OPTION, THE           
!     TRANSFORMATION U WILL BE MULTIPLIED INTO THE ARRAY Z.             
!     IF A = TRANS(X)*X, SO THAT R IS THE TRIANGULAR PART OF THE        
!     QR FACTORIZATION OF X, THEN RR IS THE TRIANGULAR PART OF THE      
!     QR FACTORIZATION OF X*E, I.E. X WITH ITS COLUMNS PERMUTED.        
!     FOR A LESS TERSE DESCRIPTION OF WHAT DCHEX DOES AND HOW           
!     IT MAY BE APPLIED, SEE THE LINPACK GUIDE.                         
!                                                                       
!     THE MATRIX Q IS DETERMINED AS THE PRODUCT U(L-K)*...*U(1)         
!     OF PLANE ROTATIONS OF THE FORM                                    
!                                                                       
!                           (    C(I)       S(I) )                      
!                           (                    ) ,                    
!                           (    -S(I)      C(I) )                      
!                                                                       
!     WHERE C(I) IS DOUBLE PRECISION, THE ROWS THESE ROTATIONS OPERATE  
!     ON ARE DESCRIBED BELOW.                                           
!                                                                       
!     THERE ARE TWO TYPES OF PERMUTATIONS, WHICH ARE DETERMINED         
!     BY THE VALUE OF JOB.                                              
!                                                                       
!     1. RIGHT CIRCULAR SHIFT (JOB = 1).                                
!                                                                       
!         THE COLUMNS ARE REARRANGED IN THE FOLLOWING ORDER.            
!                                                                       
!                1,...,K-1,L,K,K+1,...,L-1,L+1,...,P.                   
!                                                                       
!         U IS THE PRODUCT OF L-K ROTATIONS U(I), WHERE U(I)            
!         ACTS IN THE (L-I,L-I+1)-PLANE.                                
!                                                                       
!     2. LEFT CIRCULAR SHIFT (JOB = 2).                                 
!         THE COLUMNS ARE REARRANGED IN THE FOLLOWING ORDER             
!                                                                       
!                1,...,K-1,K+1,K+2,...,L,K,L+1,...,P.                   
!                                                                       
!         U IS THE PRODUCT OF L-K ROTATIONS U(I), WHERE U(I)            
!         ACTS IN THE (K+I-1,K+I)-PLANE.                                
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!         R      DOUBLE PRECISION(LDR,P), WHERE LDR.GE.P.               
!                R CONTAINS THE UPPER TRIANGULAR FACTOR                 
!                THAT IS TO BE UPDATED.  ELEMENTS OF R                  
!                BELOW THE DIAGONAL ARE NOT REFERENCED.                 
!                                                                       
!         LDR    INTEGER.                                               
!                LDR IS THE LEADING DIMENSION OF THE ARRAY R.           
!                                                                       
!         P      INTEGER.                                               
!                P IS THE ORDER OF THE MATRIX R.                        
!                                                                       
!         K      INTEGER.                                               
!                K IS THE FIRST COLUMN TO BE PERMUTED.                  
!                                                                       
!         L      INTEGER.                                               
!                L IS THE LAST COLUMN TO BE PERMUTED.                   
!                L MUST BE STRICTLY GREATER THAN K.                     
!                                                                       
!         Z      DOUBLE PRECISION(LDZ,NZ), WHERE LDZ.GE.P.              
!                Z IS AN ARRAY OF NZ P-VECTORS INTO WHICH THE           
!                TRANSFORMATION U IS MULTIPLIED.  Z IS                  
!                NOT REFERENCED IF NZ = 0.                              
!                                                                       
!         LDZ    INTEGER.                                               
!                LDZ IS THE LEADING DIMENSION OF THE ARRAY Z.           
!                                                                       
!         NZ     INTEGER.                                               
!                NZ IS THE NUMBER OF COLUMNS OF THE MATRIX Z.           
!                                                                       
!         JOB    INTEGER.                                               
!                JOB DETERMINES THE TYPE OF PERMUTATION.                
!                       JOB = 1  RIGHT CIRCULAR SHIFT.                  
!                       JOB = 2  LEFT CIRCULAR SHIFT.                   
!                                                                       
!     ON RETURN                                                         
!                                                                       
!         R      CONTAINS THE UPDATED FACTOR.                           
!                                                                       
!         Z      CONTAINS THE UPDATED MATRIX Z.                         
!                                                                       
!         C      DOUBLE PRECISION(P).                                   
!                C CONTAINS THE COSINES OF THE TRANSFORMING ROTATIONS.  
!                                                                       
!         S      DOUBLE PRECISION(P).                                   
!                S CONTAINS THE SINES OF THE TRANSFORMING ROTATIONS.    
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.       
!                                                                       
!     DCHEX USES THE FOLLOWING FUNCTIONS AND SUBROUTINES.               
!                                                                       
!     BLAS DROTG                                                        
!     FORTRAN MIN0                                                      
!                                                                       
      INTEGER I,II,IL,IU,J,JJ,KM1,KP1,LMK,LM1                           
      DOUBLE PRECISION RJP1J,T                                          
!                                                                       
!     INITIALIZE                                                        
!                                                                       
      KM1 = K - 1                                                       
      KP1 = K + 1                                                       
      LMK = L - K                                                       
      LM1 = L - 1                                                       
      !                                                                       
      !     PERFORM THE APPROPRIATE TASK.                                     
      !                                                                       
      select case (JOB)
      case (1)
        go to 10
      case (2)
        go to 130
      end select
!                                                                       
!     RIGHT CIRCULAR SHIFT.                                             
!                                                                       
   10 CONTINUE                                                          
!                                                                       
!        REORDER THE COLUMNS.                                           
!                                                                       
         DO 20 I = 1, L                                                 
            II = L - I + 1                                              
            S(I) = R(II,L)                                              
   20    CONTINUE                                                       
         DO 40 JJ = K, LM1                                              
            J = LM1 - JJ + K                                            
            DO 30 I = 1, J                                              
               R(I,J+1) = R(I,J)                                        
   30       CONTINUE                                                    
            R(J+1,J+1) = 0.0D0                                          
   40    CONTINUE                                                       
         IF (K .EQ. 1) GO TO 60                                         
            DO 50 I = 1, KM1                                            
               II = L - I + 1                                           
               R(I,K) = S(II)                                           
   50       CONTINUE                                                    
   60    CONTINUE                                                       
!                                                                       
!        CALCULATE THE ROTATIONS.                                       
!                                                                       
         T = S(1)                                                       
         DO 70 I = 1, LMK                                               
            CALL DROTG(S(I+1),T,C(I),S(I))                              
            T = S(I+1)                                                  
   70    CONTINUE                                                       
         R(K,K) = T                                                     
         DO 90 J = KP1, P                                               
            IL = MAX0(1,L-J+1)                                          
            DO 80 II = IL, LMK                                          
               I = L - II                                               
               T = C(II)*R(I,J) + S(II)*R(I+1,J)                        
               R(I+1,J) = C(II)*R(I+1,J) - S(II)*R(I,J)                 
               R(I,J) = T                                               
   80       CONTINUE                                                    
   90    CONTINUE                                                       
!                                                                       
!        IF REQUIRED, APPLY THE TRANSFORMATIONS TO Z.                   
!                                                                       
         IF (NZ .LT. 1) GO TO 120                                       
         DO 110 J = 1, NZ                                               
            DO 100 II = 1, LMK                                          
               I = L - II                                               
               T = C(II)*Z(I,J) + S(II)*Z(I+1,J)                        
               Z(I+1,J) = C(II)*Z(I+1,J) - S(II)*Z(I,J)                 
               Z(I,J) = T                                               
  100       CONTINUE                                                    
  110    CONTINUE                                                       
  120    CONTINUE                                                       
      GO TO 260                                                         
!                                                                       
!     LEFT CIRCULAR SHIFT                                               
!                                                                       
  130 CONTINUE                                                          
!                                                                       
!        REORDER THE COLUMNS                                            
!                                                                       
         DO 140 I = 1, K                                                
            II = LMK + I                                                
            S(II) = R(I,K)                                              
  140    CONTINUE                                                       
         DO 160 J = K, LM1                                              
            DO 150 I = 1, J                                             
               R(I,J) = R(I,J+1)                                        
  150       CONTINUE                                                    
            JJ = J - KM1                                                
            S(JJ) = R(J+1,J+1)                                          
  160    CONTINUE                                                       
         DO 170 I = 1, K                                                
            II = LMK + I                                                
            R(I,L) = S(II)                                              
  170    CONTINUE                                                       
         DO 180 I = KP1, L                                              
            R(I,L) = 0.0D0                                              
  180    CONTINUE                                                       
!                                                                       
!        REDUCTION LOOP.                                                
!                                                                       
         DO 220 J = K, P                                                
            IF (J .EQ. K) GO TO 200                                     
!                                                                       
!              APPLY THE ROTATIONS.                                     
!                                                                       
               IU = MIN0(J-1,L-1)                                       
               DO 190 I = K, IU                                         
                  II = I - K + 1                                        
                  T = C(II)*R(I,J) + S(II)*R(I+1,J)                     
                  R(I+1,J) = C(II)*R(I+1,J) - S(II)*R(I,J)              
                  R(I,J) = T                                            
  190          CONTINUE                                                 
  200       CONTINUE                                                    
            IF (J .GE. L) GO TO 210                                     
               JJ = J - K + 1                                           
               T = S(JJ)                                                
               CALL DROTG(R(J,J),T,C(JJ),S(JJ))                         
  210       CONTINUE                                                    
  220    CONTINUE                                                       
!                                                                       
!        APPLY THE ROTATIONS TO Z.                                      
!                                                                       
         IF (NZ .LT. 1) GO TO 250                                       
         DO 240 J = 1, NZ                                               
            DO 230 I = K, LM1                                           
               II = I - KM1                                             
               T = C(II)*Z(I,J) + S(II)*Z(I+1,J)                        
               Z(I+1,J) = C(II)*Z(I+1,J) - S(II)*Z(I,J)                 
               Z(I,J) = T                                               
  230       CONTINUE                                                    
  240    CONTINUE                                                       
  250    CONTINUE                                                       
  260 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DQRDC(X,LDX,N,P,QRAUX,JPVT,WORK,JOB)                   
      INTEGER LDX,N,P,JOB                                               
      INTEGER JPVT(1)                                                   
      DOUBLE PRECISION X(LDX,1),QRAUX(1),WORK(1)                        
!                                                                       
!     DQRDC USES HOUSEHOLDER TRANSFORMATIONS TO COMPUTE THE QR          
!     FACTORIZATION OF AN N BY P MATRIX X.  COLUMN PIVOTING             
!     BASED ON THE 2-NORMS OF THE REDUCED COLUMNS MAY BE                
!     PERFORMED AT THE USERS OPTION.                                    
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        X       DOUBLE PRECISION(LDX,P), WHERE LDX .GE. N.             
!                X CONTAINS THE MATRIX WHOSE DECOMPOSITION IS TO BE     
!                COMPUTED.                                              
!                                                                       
!        LDX     INTEGER.                                               
!                LDX IS THE LEADING DIMENSION OF THE ARRAY X.           
!                                                                       
!        N       INTEGER.                                               
!                N IS THE NUMBER OF ROWS OF THE MATRIX X.               
!                                                                       
!        P       INTEGER.                                               
!                P IS THE NUMBER OF COLUMNS OF THE MATRIX X.            
!                                                                       
!        JPVT    INTEGER(P).                                            
!                JPVT CONTAINS INTEGERS THAT CONTROL THE SELECTION      
!                OF THE PIVOT COLUMNS.  THE K-TH COLUMN X(K) OF X       
!                IS PLACED IN ONE OF THREE CLASSES ACCORDING TO THE     
!                VALUE OF JPVT(K).                                      
!                                                                       
!                   IF JPVT(K) .GT. 0, THEN X(K) IS AN INITIAL          
!                                      COLUMN.                          
!                                                                       
!                   IF JPVT(K) .EQ. 0, THEN X(K) IS A FREE COLUMN.      
!                                                                       
!                   IF JPVT(K) .LT. 0, THEN X(K) IS A FINAL COLUMN.     
!                                                                       
!                BEFORE THE DECOMPOSITION IS COMPUTED, INITIAL COLUMNS  
!                ARE MOVED TO THE BEGINNING OF THE ARRAY X AND FINAL    
!                COLUMNS TO THE END.  BOTH INITIAL AND FINAL COLUMNS    
!                ARE FROZEN IN PLACE DURING THE COMPUTATION AND ONLY    
!                FREE COLUMNS ARE MOVED.  AT THE K-TH STAGE OF THE      
!                REDUCTION, IF X(K) IS OCCUPIED BY A FREE COLUMN        
!                IT IS INTERCHANGED WITH THE FREE COLUMN OF LARGEST     
!                REDUCED NORM.  JPVT IS NOT REFERENCED IF               
!                JOB .EQ. 0.                                            
!                                                                       
!        WORK    DOUBLE PRECISION(P).                                   
!                WORK IS A WORK ARRAY.  WORK IS NOT REFERENCED IF       
!                JOB .EQ. 0.                                            
!                                                                       
!        JOB     INTEGER.                                               
!                JOB IS AN INTEGER THAT INITIATES COLUMN PIVOTING.      
!                IF JOB .EQ. 0, NO PIVOTING IS DONE.                    
!                IF JOB .NE. 0, PIVOTING IS DONE.                       
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        X       X CONTAINS IN ITS UPPER TRIANGLE THE UPPER             
!                TRIANGULAR MATRIX R OF THE QR FACTORIZATION.           
!                BELOW ITS DIAGONAL X CONTAINS INFORMATION FROM         
!                WHICH THE ORTHOGONAL PART OF THE DECOMPOSITION         
!                CAN BE RECOVERED.  NOTE THAT IF PIVOTING HAS           
!                BEEN REQUESTED, THE DECOMPOSITION IS NOT THAT          
!                OF THE ORIGINAL MATRIX X BUT THAT OF X                 
!                WITH ITS COLUMNS PERMUTED AS DESCRIBED BY JPVT.        
!                                                                       
!        QRAUX   DOUBLE PRECISION(P).                                   
!                QRAUX CONTAINS FURTHER INFORMATION REQUIRED TO RECOVER 
!                THE ORTHOGONAL PART OF THE DECOMPOSITION.              
!                                                                       
!        JPVT    JPVT(K) CONTAINS THE INDEX OF THE COLUMN OF THE        
!                ORIGINAL MATRIX THAT HAS BEEN INTERCHANGED INTO        
!                THE K-TH COLUMN, IF PIVOTING WAS REQUESTED.            
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.       
!                                                                       
!     DQRDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.               
!                                                                       
!     BLAS DAXPY,DDOT,DSCAL,DSWAP,DNRM2                                 
!     FORTRAN DABS,DMAX1,MIN0,DSQRT                                     
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      INTEGER J,JP,L,LP1,LUP,MAXJ,PL,PU                                 
      DOUBLE PRECISION MAXNRM,DNRM2,TT                                  
      DOUBLE PRECISION DDOT,NRMXL,T                                     
      LOGICAL NEGJ,SWAPJ                                                
!                                                                       
!                                                                       
      PL = 1                                                            
      PU = 0                                                            
      IF (JOB .EQ. 0) GO TO 60                                          
!                                                                       
!        PIVOTING HAS BEEN REQUESTED.  REARRANGE THE COLUMNS            
!        ACCORDING TO JPVT.                                             
!                                                                       
         DO 20 J = 1, P                                                 
            SWAPJ = JPVT(J) .GT. 0                                      
            NEGJ = JPVT(J) .LT. 0                                       
            JPVT(J) = J                                                 
            IF (NEGJ) JPVT(J) = -J                                      
            IF (.NOT.SWAPJ) GO TO 10                                    
               IF (J .NE. PL) CALL DSWAP(N,X(1,PL),1,X(1,J),1)          
               JPVT(J) = JPVT(PL)                                       
               JPVT(PL) = J                                             
               PL = PL + 1                                              
   10       CONTINUE                                                    
   20    CONTINUE                                                       
         PU = P                                                         
         DO 50 JJ = 1, P                                                
            J = P - JJ + 1                                              
            IF (JPVT(J) .GE. 0) GO TO 40                                
               JPVT(J) = -JPVT(J)                                       
               IF (J .EQ. PU) GO TO 30                                  
                  CALL DSWAP(N,X(1,PU),1,X(1,J),1)                      
                  JP = JPVT(PU)                                         
                  JPVT(PU) = JPVT(J)                                    
                  JPVT(J) = JP                                          
   30          CONTINUE                                                 
               PU = PU - 1                                              
   40       CONTINUE                                                    
   50    CONTINUE                                                       
   60 CONTINUE                                                          
!                                                                       
!     COMPUTE THE NORMS OF THE FREE COLUMNS.                            
!                                                                       
      IF (PU .LT. PL) GO TO 80                                          
      DO 70 J = PL, PU                                                  
         QRAUX(J) = DNRM2(N,X(1,J),1)                                   
         WORK(J) = QRAUX(J)                                             
   70 CONTINUE                                                          
   80 CONTINUE                                                          
!                                                                       
!     PERFORM THE HOUSEHOLDER REDUCTION OF X.                           
!                                                                       
      LUP = MIN0(N,P)                                                   
      DO 200 L = 1, LUP                                                 
         IF (L .LT. PL .OR. L .GE. PU) GO TO 120                        
!                                                                       
!           LOCATE THE COLUMN OF LARGEST NORM AND BRING IT              
!           INTO THE PIVOT POSITION.                                    
!                                                                       
            MAXNRM = 0.0D0                                              
            MAXJ = L                                                    
            DO 100 J = L, PU                                            
               IF (QRAUX(J) .LE. MAXNRM) GO TO 90                       
                  MAXNRM = QRAUX(J)                                     
                  MAXJ = J                                              
   90          CONTINUE                                                 
  100       CONTINUE                                                    
            IF (MAXJ .EQ. L) GO TO 110                                  
               CALL DSWAP(N,X(1,L),1,X(1,MAXJ),1)                       
               QRAUX(MAXJ) = QRAUX(L)                                   
               WORK(MAXJ) = WORK(L)                                     
               JP = JPVT(MAXJ)                                          
               JPVT(MAXJ) = JPVT(L)                                     
               JPVT(L) = JP                                             
  110       CONTINUE                                                    
  120    CONTINUE                                                       
         QRAUX(L) = 0.0D0                                               
         IF (L .EQ. N) GO TO 190                                        
!                                                                       
!           COMPUTE THE HOUSEHOLDER TRANSFORMATION FOR COLUMN L.        
!                                                                       
            NRMXL = DNRM2(N-L+1,X(L,L),1)                               
            IF (NRMXL .EQ. 0.0D0) GO TO 180                             
               IF (X(L,L) .NE. 0.0D0) NRMXL = DSIGN(NRMXL,X(L,L))       
               CALL DSCAL(N-L+1,1.0D0/NRMXL,X(L,L),1)                   
               X(L,L) = 1.0D0 + X(L,L)                                  
!                                                                       
!              APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS,       
!              UPDATING THE NORMS.                                      
!                                                                       
               LP1 = L + 1                                              
               IF (P .LT. LP1) GO TO 170                                
               DO 160 J = LP1, P                                        
                  T = -DDOT(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)             
                  CALL DAXPY(N-L+1,T,X(L,L),1,X(L,J),1)                 
                  IF (J .LT. PL .OR. J .GT. PU) GO TO 150               
                  IF (QRAUX(J) .EQ. 0.0D0) GO TO 150                    
                     TT = 1.0D0 - (DABS(X(L,J))/QRAUX(J))**2            
                     TT = DMAX1(TT,0.0D0)                               
                     T = TT                                             
                     TT = 1.0D0 + 0.05D0*TT*(QRAUX(J)/WORK(J))**2       
                     IF (TT .EQ. 1.0D0) GO TO 130                       
                        QRAUX(J) = QRAUX(J)*DSQRT(T)                    
                     GO TO 140                                          
  130                CONTINUE                                           
                        QRAUX(J) = DNRM2(N-L,X(L+1,J),1)                
                        WORK(J) = QRAUX(J)                              
  140                CONTINUE                                           
  150             CONTINUE                                              
  160          CONTINUE                                                 
  170          CONTINUE                                                 
!                                                                       
!              SAVE THE TRANSFORMATION.                                 
!                                                                       
               QRAUX(L) = X(L,L)                                        
               X(L,L) = -NRMXL                                          
  180       CONTINUE                                                    
  190    CONTINUE                                                       
  200 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DQRSL(X,LDX,N,K,QRAUX,Y,QY,QTY,B,RSD,XB,JOB,INFO)      
      INTEGER LDX,N,K,JOB,INFO                                          
      DOUBLE PRECISION X(LDX,1),QRAUX(1),Y(1),QY(1),QTY(1),B(1),RSD(1),  &
                       XB(1)                                            
!                                                                       
!     DQRSL APPLIES THE OUTPUT OF DQRDC TO COMPUTE COORDINATE           
!     TRANSFORMATIONS, PROJECTIONS, AND LEAST SQUARES SOLUTIONS.        
!     FOR K .LE. MIN(N,P), LET XK BE THE MATRIX                         
!                                                                       
!            XK = (X(JPVT(1)),X(JPVT(2)), ... ,X(JPVT(K)))              
!                                                                       
!     FORMED FROM COLUMNNS JPVT(1), ... ,JPVT(K) OF THE ORIGINAL        
!     N X P MATRIX X THAT WAS INPUT TO DQRDC (IF NO PIVOTING WAS        
!     DONE, XK CONSISTS OF THE FIRST K COLUMNS OF X IN THEIR            
!     ORIGINAL ORDER).  DQRDC PRODUCES A FACTORED ORTHOGONAL MATRIX Q   
!     AND AN UPPER TRIANGULAR MATRIX R SUCH THAT                        
!                                                                       
!              XK = Q * (R)                                             
!                       (0)                                             
!                                                                       
!     THIS INFORMATION IS CONTAINED IN CODED FORM IN THE ARRAYS         
!     X AND QRAUX.                                                      
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        X      DOUBLE PRECISION(LDX,P).                                
!               X CONTAINS THE OUTPUT OF DQRDC.                         
!                                                                       
!        LDX    INTEGER.                                                
!               LDX IS THE LEADING DIMENSION OF THE ARRAY X.            
!                                                                       
!        N      INTEGER.                                                
!               N IS THE NUMBER OF ROWS OF THE MATRIX XK.  IT MUST      
!               HAVE THE SAME VALUE AS N IN DQRDC.                      
!                                                                       
!        K      INTEGER.                                                
!               K IS THE NUMBER OF COLUMNS OF THE MATRIX XK.  K         
!               MUST NNOT BE GREATER THAN MIN(N,P), WHERE P IS THE      
!               SAME AS IN THE CALLING SEQUENCE TO DQRDC.               
!                                                                       
!        QRAUX  DOUBLE PRECISION(P).                                    
!               QRAUX CONTAINS THE AUXILIARY OUTPUT FROM DQRDC.         
!                                                                       
!        Y      DOUBLE PRECISION(N)                                     
!               Y CONTAINS AN N-VECTOR THAT IS TO BE MANIPULATED        
!               BY DQRSL.                                               
!                                                                       
!        JOB    INTEGER.                                                
!               JOB SPECIFIES WHAT IS TO BE COMPUTED.  JOB HAS          
!               THE DECIMAL EXPANSION ABCDE, WITH THE FOLLOWING         
!               MEANING.                                                
!                                                                       
!                    IF A.NE.0, COMPUTE QY.                             
!                    IF B,C,D, OR E .NE. 0, COMPUTE QTY.                
!                    IF C.NE.0, COMPUTE B.                              
!                    IF D.NE.0, COMPUTE RSD.                            
!                    IF E.NE.0, COMPUTE XB.                             
!                                                                       
!               NOTE THAT A REQUEST TO COMPUTE B, RSD, OR XB            
!               AUTOMATICALLY TRIGGERS THE COMPUTATION OF QTY, FOR      
!               WHICH AN ARRAY MUST BE PROVIDED IN THE CALLING          
!               SEQUENCE.                                               
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        QY     DOUBLE PRECISION(N).                                    
!               QY CONNTAINS Q*Y, IF ITS COMPUTATION HAS BEEN           
!               REQUESTED.                                              
!                                                                       
!        QTY    DOUBLE PRECISION(N).                                    
!               QTY CONTAINS TRANS(Q)*Y, IF ITS COMPUTATION HAS         
!               BEEN REQUESTED.  HERE TRANS(Q) IS THE                   
!               TRANSPOSE OF THE MATRIX Q.                              
!                                                                       
!        B      DOUBLE PRECISION(K)                                     
!               B CONTAINS THE SOLUTION OF THE LEAST SQUARES PROBLEM    
!                                                                       
!                    MINIMIZE NORM2(Y - XK*B),                          
!                                                                       
!               IF ITS COMPUTATION HAS BEEN REQUESTED.  (NOTE THAT      
!               IF PIVOTING WAS REQUESTED IN DQRDC, THE J-TH            
!               COMPONENT OF B WILL BE ASSOCIATED WITH COLUMN JPVT(J)   
!               OF THE ORIGINAL MATRIX X THAT WAS INPUT INTO DQRDC.)    
!                                                                       
!        RSD    DOUBLE PRECISION(N).                                    
!               RSD CONTAINS THE LEAST SQUARES RESIDUAL Y - XK*B,       
!               IF ITS COMPUTATION HAS BEEN REQUESTED.  RSD IS          
!               ALSO THE ORTHOGONAL PROJECTION OF Y ONTO THE            
!               ORTHOGONAL COMPLEMENT OF THE COLUMN SPACE OF XK.        
!                                                                       
!        XB     DOUBLE PRECISION(N).                                    
!               XB CONTAINS THE LEAST SQUARES APPROXIMATION XK*B,       
!               IF ITS COMPUTATION HAS BEEN REQUESTED.  XB IS ALSO      
!               THE ORTHOGONAL PROJECTION OF Y ONTO THE COLUMN SPACE    
!               OF X.                                                   
!                                                                       
!        INFO   INTEGER.                                                
!               INFO IS ZERO UNLESS THE COMPUTATION OF B HAS            
!               BEEN REQUESTED AND R IS EXACTLY SINGULAR.  IN           
!               THIS CASE, INFO IS THE INDEX OF THE FIRST ZERO          
!               DIAGONAL ELEMENT OF R AND B IS LEFT UNALTERED.          
!                                                                       
!     THE PARAMETERS QY, QTY, B, RSD, AND XB ARE NOT REFERENCED         
!     IF THEIR COMPUTATION IS NOT REQUESTED AND IN THIS CASE            
!     CAN BE REPLACED BY DUMMY VARIABLES IN THE CALLING PROGRAM.        
!     TO SAVE STORAGE, THE USER MAY IN SOME CASES USE THE SAME          
!     ARRAY FOR DIFFERENT PARAMETERS IN THE CALLING SEQUENCE.  A        
!     FREQUENTLY OCCURING EXAMPLE IS WHEN ONE WISHES TO COMPUTE         
!     ANY OF B, RSD, OR XB AND DOES NOT NEED Y OR QTY.  IN THIS         
!     CASE ONE MAY IDENTIFY Y, QTY, AND ONE OF B, RSD, OR XB, WHILE     
!     PROVIDING SEPARATE ARRAYS FOR ANYTHING ELSE THAT IS TO BE         
!     COMPUTED.  THUS THE CALLING SEQUENCE                              
!                                                                       
!          CALL DQRSL(X,LDX,N,K,QRAUX,Y,DUM,Y,B,Y,DUM,110,INFO)         
!                                                                       
!     WILL RESULT IN THE COMPUTATION OF B AND RSD, WITH RSD             
!     OVERWRITING Y.  MORE GENERALLY, EACH ITEM IN THE FOLLOWING        
!     LIST CONTAINS GROUPS OF PERMISSIBLE IDENTIFICATIONS FOR           
!     A SINGLE CALLINNG SEQUENCE.                                       
!                                                                       
!          1. (Y,QTY,B) (RSD) (XB) (QY)                                 
!                                                                       
!          2. (Y,QTY,RSD) (B) (XB) (QY)                                 
!                                                                       
!          3. (Y,QTY,XB) (B) (RSD) (QY)                                 
!                                                                       
!          4. (Y,QY) (QTY,B) (RSD) (XB)                                 
!                                                                       
!          5. (Y,QY) (QTY,RSD) (B) (XB)                                 
!                                                                       
!          6. (Y,QY) (QTY,XB) (B) (RSD)                                 
!                                                                       
!     IN ANY GROUP THE VALUE RETURNED IN THE ARRAY ALLOCATED TO         
!     THE GROUP CORRESPONDS TO THE LAST MEMBER OF THE GROUP.            
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.       
!                                                                       
!     DQRSL USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.               
!                                                                       
!     BLAS DAXPY,DCOPY,DDOT                                             
!     FORTRAN DABS,MIN0,MOD                                             
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      INTEGER I,J,JJ,JU,KP1                                             
      DOUBLE PRECISION DDOT,T,TEMP                                      
      LOGICAL CB,CQY,CQTY,CR,CXB                                        
!                                                                       
!                                                                       
!     SET INFO FLAG.                                                    
!                                                                       
      INFO = 0                                                          
!                                                                       
!     DETERMINE WHAT IS TO BE COMPUTED.                                 
!                                                                       
      CQY = JOB/10000 .NE. 0                                            
      CQTY = MOD(JOB,10000) .NE. 0                                      
      CB = MOD(JOB,1000)/100 .NE. 0                                     
      CR = MOD(JOB,100)/10 .NE. 0                                       
      CXB = MOD(JOB,10) .NE. 0                                          
      JU = MIN0(K,N-1)                                                  
!                                                                       
!     SPECIAL ACTION WHEN N=1.                                          
!                                                                       
      IF (JU .NE. 0) GO TO 40                                           
         IF (CQY) QY(1) = Y(1)                                          
         IF (CQTY) QTY(1) = Y(1)                                        
         IF (CXB) XB(1) = Y(1)                                          
         IF (.NOT.CB) GO TO 30                                          
            IF (X(1,1) .NE. 0.0D0) GO TO 10                             
               INFO = 1                                                 
            GO TO 20                                                    
   10       CONTINUE                                                    
               B(1) = Y(1)/X(1,1)                                       
   20       CONTINUE                                                    
   30    CONTINUE                                                       
         IF (CR) RSD(1) = 0.0D0                                         
      GO TO 250                                                         
   40 CONTINUE                                                          
!                                                                       
!        SET UP TO COMPUTE QY OR QTY.                                   
!                                                                       
         IF (CQY) CALL DCOPY(N,Y,1,QY,1)                                
         IF (CQTY) CALL DCOPY(N,Y,1,QTY,1)                              
         IF (.NOT.CQY) GO TO 70                                         
!                                                                       
!           COMPUTE QY.                                                 
!                                                                       
            DO 60 JJ = 1, JU                                            
               J = JU - JJ + 1                                          
               IF (QRAUX(J) .EQ. 0.0D0) GO TO 50                        
                  TEMP = X(J,J)                                         
                  X(J,J) = QRAUX(J)                                     
                  T = -DDOT(N-J+1,X(J,J),1,QY(J),1)/X(J,J)              
                  CALL DAXPY(N-J+1,T,X(J,J),1,QY(J),1)                  
                  X(J,J) = TEMP                                         
   50          CONTINUE                                                 
   60       CONTINUE                                                    
   70    CONTINUE                                                       
         IF (.NOT.CQTY) GO TO 100                                       
!                                                                       
!           COMPUTE TRANS(Q)*Y.                                         
!                                                                       
            DO 90 J = 1, JU                                             
               IF (QRAUX(J) .EQ. 0.0D0) GO TO 80                        
                  TEMP = X(J,J)                                         
                  X(J,J) = QRAUX(J)                                     
                  T = -DDOT(N-J+1,X(J,J),1,QTY(J),1)/X(J,J)             
                  CALL DAXPY(N-J+1,T,X(J,J),1,QTY(J),1)                 
                  X(J,J) = TEMP                                         
   80          CONTINUE                                                 
   90       CONTINUE                                                    
  100    CONTINUE                                                       
!                                                                       
!        SET UP TO COMPUTE B, RSD, OR XB.                               
!                                                                       
         IF (CB) CALL DCOPY(K,QTY,1,B,1)                                
         KP1 = K + 1                                                    
         IF (CXB) CALL DCOPY(K,QTY,1,XB,1)                              
         IF (CR .AND. K .LT. N) CALL DCOPY(N-K,QTY(KP1),1,RSD(KP1),1)   
         IF (.NOT.CXB .OR. KP1 .GT. N) GO TO 120                        
            DO 110 I = KP1, N                                           
               XB(I) = 0.0D0                                            
  110       CONTINUE                                                    
  120    CONTINUE                                                       
         IF (.NOT.CR) GO TO 140                                         
            DO 130 I = 1, K                                             
               RSD(I) = 0.0D0                                           
  130       CONTINUE                                                    
  140    CONTINUE                                                       
         IF (.NOT.CB) GO TO 190                                         
!                                                                       
!           COMPUTE B.                                                  
!                                                                       
            DO 170 JJ = 1, K                                            
               J = K - JJ + 1                                           
               IF (X(J,J) .NE. 0.0D0) GO TO 150                         
                  INFO = J                                              
!           ......EXIT                                                  
                  GO TO 180                                             
  150          CONTINUE                                                 
               B(J) = B(J)/X(J,J)                                       
               IF (J .EQ. 1) GO TO 160                                  
                  T = -B(J)                                             
                  CALL DAXPY(J-1,T,X(1,J),1,B,1)                        
  160          CONTINUE                                                 
  170       CONTINUE                                                    
  180       CONTINUE                                                    
  190    CONTINUE                                                       
         IF (.NOT.CR .AND. .NOT.CXB) GO TO 240                          
!                                                                       
!           COMPUTE RSD OR XB AS REQUIRED.                              
!                                                                       
            DO 230 JJ = 1, JU                                           
               J = JU - JJ + 1                                          
               IF (QRAUX(J) .EQ. 0.0D0) GO TO 220                       
                  TEMP = X(J,J)                                         
                  X(J,J) = QRAUX(J)                                     
                  IF (.NOT.CR) GO TO 200                                
                     T = -DDOT(N-J+1,X(J,J),1,RSD(J),1)/X(J,J)          
                     CALL DAXPY(N-J+1,T,X(J,J),1,RSD(J),1)              
  200             CONTINUE                                              
                  IF (.NOT.CXB) GO TO 210                               
                     T = -DDOT(N-J+1,X(J,J),1,XB(J),1)/X(J,J)           
                     CALL DAXPY(N-J+1,T,X(J,J),1,XB(J),1)               
  210             CONTINUE                                              
                  X(J,J) = TEMP                                         
  220          CONTINUE                                                 
  230       CONTINUE                                                    
  240    CONTINUE                                                       
  250 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DSVDC(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,JOB,INFO)         
      INTEGER LDX,N,P,LDU,LDV,JOB,INFO                                  
      DOUBLE PRECISION X(LDX,1),S(1),E(1),U(LDU,1),V(LDV,1),WORK(1)     
!                                                                       
!                                                                       
!     DSVDC IS A SUBROUTINE TO REDUCE A DOUBLE PRECISION NXP MATRIX X   
!     BY ORTHOGONAL TRANSFORMATIONS U AND V TO DIAGONAL FORM.  THE      
!     DIAGONAL ELEMENTS S(I) ARE THE SINGULAR VALUES OF X.  THE         
!     COLUMNS OF U ARE THE CORRESPONDING LEFT SINGULAR VECTORS,         
!     AND THE COLUMNS OF V THE RIGHT SINGULAR VECTORS.                  
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!         X         DOUBLE PRECISION(LDX,P), WHERE LDX.GE.N.            
!                   X CONTAINS THE MATRIX WHOSE SINGULAR VALUE          
!                   DECOMPOSITION IS TO BE COMPUTED.  X IS              
!                   DESTROYED BY DSVDC.                                 
!                                                                       
!         LDX       INTEGER.                                            
!                   LDX IS THE LEADING DIMENSION OF THE ARRAY X.        
!                                                                       
!         N         INTEGER.                                            
!                   N IS THE NUMBER OF COLUMNS OF THE MATRIX X.         
!                                                                       
!         P         INTEGER.                                            
!                   P IS THE NUMBER OF ROWS OF THE MATRIX X.            
!                                                                       
!         LDU       INTEGER.                                            
!                   LDU IS THE LEADING DIMENSION OF THE ARRAY U.        
!                   (SEE BELOW).                                        
!                                                                       
!         LDV       INTEGER.                                            
!                   LDV IS THE LEADING DIMENSION OF THE ARRAY V.        
!                   (SEE BELOW).                                        
!                                                                       
!         WORK      DOUBLE PRECISION(N).                                
!                   WORK IS A SCRATCH ARRAY.                            
!                                                                       
!         JOB       INTEGER.                                            
!                   JOB CONTROLS THE COMPUTATION OF THE SINGULAR        
!                   VECTORS.  IT HAS THE DECIMAL EXPANSION AB           
!                   WITH THE FOLLOWING MEANING                          
!                                                                       
!                        A.EQ.0    DO NOT COMPUTE THE LEFT SINGULAR     
!                                  VECTORS.                             
!                        A.EQ.1    RETURN THE N LEFT SINGULAR VECTORS   
!                                  IN U.                                
!                        A.GE.2    RETURN THE FIRST MIN(N,P) SINGULAR   
!                                  VECTORS IN U.                        
!                        B.EQ.0    DO NOT COMPUTE THE RIGHT SINGULAR    
!                                  VECTORS.                             
!                        B.EQ.1    RETURN THE RIGHT SINGULAR VECTORS    
!                                  IN V.                                
!                                                                       
!     ON RETURN                                                         
!                                                                       
!         S         DOUBLE PRECISION(MM), WHERE MM=MIN(N+1,P).          
!                   THE FIRST MIN(N,P) ENTRIES OF S CONTAIN THE         
!                   SINGULAR VALUES OF X ARRANGED IN DESCENDING         
!                   ORDER OF MAGNITUDE.                                 
!                                                                       
!         E         DOUBLE PRECISION(P).                                
!                   E ORDINARILY CONTAINS ZEROS.  HOWEVER SEE THE       
!                   DISCUSSION OF INFO FOR EXCEPTIONS.                  
!                                                                       
!         U         DOUBLE PRECISION(LDU,K), WHERE LDU.GE.N.  IF        
!                                   JOBA.EQ.1 THEN K.EQ.N, IF JOBA.GE.2 
!                                   THEN K.EQ.MIN(N,P).                 
!                   U CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS.    
!                   U IS NOT REFERENCED IF JOBA.EQ.0.  IF N.LE.P        
!                   OR IF JOBA.EQ.2, THEN U MAY BE IDENTIFIED WITH X    
!                   IN THE SUBROUTINE CALL.                             
!                                                                       
!         V         DOUBLE PRECISION(LDV,P), WHERE LDV.GE.P.            
!                   V CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS.    
!                   V IS NOT REFERENCED IF JOB.EQ.0.  IF P.LE.N,        
!                   THEN V MAY BE IDENTIFIED WITH X IN THE              
!                   SUBROUTINE CALL.                                    
!                                                                       
!         INFO      INTEGER.                                            
!                   THE SINGULAR VALUES (AND THEIR CORRESPONDING        
!                   SINGULAR VECTORS) S(INFO+1),S(INFO+2),...,S(M)      
!                   ARE CORRECT (HERE M=MIN(N,P)).  THUS IF             
!                   INFO.EQ.0, ALL THE SINGULAR VALUES AND THEIR        
!                   VECTORS ARE CORRECT.  IN ANY EVENT, THE MATRIX      
!                   B = TRANS(U)*X*V IS THE BIDIAGONAL MATRIX           
!                   WITH THE ELEMENTS OF S ON ITS DIAGONAL AND THE      
!                   ELEMENTS OF E ON ITS SUPER-DIAGONAL (TRANS(U)       
!                   IS THE TRANSPOSE OF U).  THUS THE SINGULAR          
!                   VALUES OF X AND B ARE THE SAME.                     
!                                                                       
!     LINPACK. THIS VERSION DATED 03/19/79 .                            
!     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.       
!                                                                       
!     DSVDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.               
!                                                                       
!     EXTERNAL DROT                                                     
!     BLAS DAXPY,DDOT,DSCAL,DSWAP,DNRM2,DROTG                           
!     FORTRAN DABS,DMAX1,MAX0,MIN0,MOD,DSQRT                            
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      INTEGER I,ITER,J,JOBU,K,KASE,KK,L,LL,LLS,LM1,LP1,LS,LU,M,MAXIT,    &
              MM,MM1,MP1,NCT,NCTP1,NCU,NRT,NRTP1                        
      DOUBLE PRECISION DDOT,T,R                                         
      DOUBLE PRECISION B,C,CS,EL,EMM1,F,G,DNRM2,SCALE,SHIFT,SL,SM,SN,    &
                       SMM1,T1,TEST,ZTEST                               
      LOGICAL WANTU,WANTV                                               
!                                                                       
!                                                                       
!     SET THE MAXIMUM NUMBER OF ITERATIONS.                             
!                                                                       
      MAXIT = 30                                                        
!                                                                       
!     DETERMINE WHAT IS TO BE COMPUTED.                                 
!                                                                       
      WANTU = .FALSE.                                                   
      WANTV = .FALSE.                                                   
      JOBU = MOD(JOB,100)/10                                            
      NCU = N                                                           
      IF (JOBU .GT. 1) NCU = MIN0(N,P)                                  
      IF (JOBU .NE. 0) WANTU = .TRUE.                                   
      IF (MOD(JOB,10) .NE. 0) WANTV = .TRUE.                            
!                                                                       
!     REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS        
!     IN S AND THE SUPER-DIAGONAL ELEMENTS IN E.                        
!                                                                       
      INFO = 0                                                          
      NCT = MIN0(N-1,P)                                                 
      NRT = MAX0(0,MIN0(P-2,N))                                         
      LU = MAX0(NCT,NRT)                                                
      IF (LU .LT. 1) GO TO 170                                          
      DO 160 L = 1, LU                                                  
         LP1 = L + 1                                                    
         IF (L .GT. NCT) GO TO 20                                       
!                                                                       
!           COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND          
!           PLACE THE L-TH DIAGONAL IN S(L).                            
!                                                                       
            S(L) = DNRM2(N-L+1,X(L,L),1)                                
            IF (S(L) .EQ. 0.0D0) GO TO 10                               
               IF (X(L,L) .NE. 0.0D0) S(L) = DSIGN(S(L),X(L,L))         
               CALL DSCAL(N-L+1,1.0D0/S(L),X(L,L),1)                    
               X(L,L) = 1.0D0 + X(L,L)                                  
   10       CONTINUE                                                    
            S(L) = -S(L)                                                
   20    CONTINUE                                                       
         IF (P .LT. LP1) GO TO 50                                       
         DO 40 J = LP1, P                                               
            IF (L .GT. NCT) GO TO 30                                    
            IF (S(L) .EQ. 0.0D0) GO TO 30                               
!                                                                       
!              APPLY THE TRANSFORMATION.                                
!                                                                       
               T = -DDOT(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)                
               CALL DAXPY(N-L+1,T,X(L,L),1,X(L,J),1)                    
   30       CONTINUE                                                    
!                                                                       
!           PLACE THE L-TH ROW OF X INTO  E FOR THE                     
!           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION.           
!                                                                       
            E(J) = X(L,J)                                               
   40    CONTINUE                                                       
   50    CONTINUE                                                       
         IF (.NOT.WANTU .OR. L .GT. NCT) GO TO 70                       
!                                                                       
!           PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK           
!           MULTIPLICATION.                                             
!                                                                       
            DO 60 I = L, N                                              
               U(I,L) = X(I,L)                                          
   60       CONTINUE                                                    
   70    CONTINUE                                                       
         IF (L .GT. NRT) GO TO 150                                      
!                                                                       
!           COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE           
!           L-TH SUPER-DIAGONAL IN E(L).                                
!                                                                       
            E(L) = DNRM2(P-L,E(LP1),1)                                  
            IF (E(L) .EQ. 0.0D0) GO TO 80                               
               IF (E(LP1) .NE. 0.0D0) E(L) = DSIGN(E(L),E(LP1))         
               CALL DSCAL(P-L,1.0D0/E(L),E(LP1),1)                      
               E(LP1) = 1.0D0 + E(LP1)                                  
   80       CONTINUE                                                    
            E(L) = -E(L)                                                
            IF (LP1 .GT. N .OR. E(L) .EQ. 0.0D0) GO TO 120              
!                                                                       
!              APPLY THE TRANSFORMATION.                                
!                                                                       
               DO 90 I = LP1, N                                         
                  WORK(I) = 0.0D0                                       
   90          CONTINUE                                                 
               DO 100 J = LP1, P                                        
                  CALL DAXPY(N-L,E(J),X(LP1,J),1,WORK(LP1),1)           
  100          CONTINUE                                                 
               DO 110 J = LP1, P                                        
                  CALL DAXPY(N-L,-E(J)/E(LP1),WORK(LP1),1,X(LP1,J),1)   
  110          CONTINUE                                                 
  120       CONTINUE                                                    
            IF (.NOT.WANTV) GO TO 140                                   
!                                                                       
!              PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT             
!              BACK MULTIPLICATION.                                     
!                                                                       
               DO 130 I = LP1, P                                        
                  V(I,L) = E(I)                                         
  130          CONTINUE                                                 
  140       CONTINUE                                                    
  150    CONTINUE                                                       
  160 CONTINUE                                                          
  170 CONTINUE                                                          
!                                                                       
!     SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M.                    
!                                                                       
      M = MIN0(P,N+1)                                                   
      NCTP1 = NCT + 1                                                   
      NRTP1 = NRT + 1                                                   
      IF (NCT .LT. P) S(NCTP1) = X(NCTP1,NCTP1)                         
      IF (N .LT. M) S(M) = 0.0D0                                        
      IF (NRTP1 .LT. M) E(NRTP1) = X(NRTP1,M)                           
      E(M) = 0.0D0                                                      
!                                                                       
!     IF REQUIRED, GENERATE U.                                          
!                                                                       
      IF (.NOT.WANTU) GO TO 300                                         
         IF (NCU .LT. NCTP1) GO TO 200                                  
         DO 190 J = NCTP1, NCU                                          
            DO 180 I = 1, N                                             
               U(I,J) = 0.0D0                                           
  180       CONTINUE                                                    
            U(J,J) = 1.0D0                                              
  190    CONTINUE                                                       
  200    CONTINUE                                                       
         IF (NCT .LT. 1) GO TO 290                                      
         DO 280 LL = 1, NCT                                             
            L = NCT - LL + 1                                            
            IF (S(L) .EQ. 0.0D0) GO TO 250                              
               LP1 = L + 1                                              
               IF (NCU .LT. LP1) GO TO 220                              
               DO 210 J = LP1, NCU                                      
                  T = -DDOT(N-L+1,U(L,L),1,U(L,J),1)/U(L,L)             
                  CALL DAXPY(N-L+1,T,U(L,L),1,U(L,J),1)                 
  210          CONTINUE                                                 
  220          CONTINUE                                                 
               CALL DSCAL(N-L+1,-1.0D0,U(L,L),1)                        
               U(L,L) = 1.0D0 + U(L,L)                                  
               LM1 = L - 1                                              
               IF (LM1 .LT. 1) GO TO 240                                
               DO 230 I = 1, LM1                                        
                  U(I,L) = 0.0D0                                        
  230          CONTINUE                                                 
  240          CONTINUE                                                 
            GO TO 270                                                   
  250       CONTINUE                                                    
               DO 260 I = 1, N                                          
                  U(I,L) = 0.0D0                                        
  260          CONTINUE                                                 
               U(L,L) = 1.0D0                                           
  270       CONTINUE                                                    
  280    CONTINUE                                                       
  290    CONTINUE                                                       
  300 CONTINUE                                                          
!                                                                       
!     IF IT IS REQUIRED, GENERATE V.                                    
!                                                                       
      IF (.NOT.WANTV) GO TO 350                                         
         DO 340 LL = 1, P                                               
            L = P - LL + 1                                              
            LP1 = L + 1                                                 
            IF (L .GT. NRT) GO TO 320                                   
            IF (E(L) .EQ. 0.0D0) GO TO 320                              
               DO 310 J = LP1, P                                        
                  T = -DDOT(P-L,V(LP1,L),1,V(LP1,J),1)/V(LP1,L)         
                  CALL DAXPY(P-L,T,V(LP1,L),1,V(LP1,J),1)               
  310          CONTINUE                                                 
  320       CONTINUE                                                    
            DO 330 I = 1, P                                             
               V(I,L) = 0.0D0                                           
  330       CONTINUE                                                    
            V(L,L) = 1.0D0                                              
  340    CONTINUE                                                       
  350 CONTINUE                                                          
!                                                                       
!     MAIN ITERATION LOOP FOR THE SINGULAR VALUES.                      
!                                                                       
      MM = M                                                            
      ITER = 0                                                          
  360 CONTINUE                                                          
!                                                                       
!        QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND.               
!                                                                       
!     ...EXIT                                                           
         IF (M .EQ. 0) GO TO 620                                        
!                                                                       
!        IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET                
!        FLAG AND RETURN.                                               
!                                                                       
         IF (ITER .LT. MAXIT) GO TO 370                                 
            INFO = M                                                    
!     ......EXIT                                                        
            GO TO 620                                                   
  370    CONTINUE                                                       
!                                                                       
!        THIS SECTION OF THE PROGRAM INSPECTS FOR                       
!        NEGLIGIBLE ELEMENTS IN THE S AND E ARRAYS.  ON                 
!        COMPLETION THE VARIABLES KASE AND L ARE SET AS FOLLOWS.        
!                                                                       
!           KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M   
!           KASE = 2     IF S(L) IS NEGLIGIBLE AND L.LT.M               
!           KASE = 3     IF E(L-1) IS NEGLIGIBLE, L.LT.M, AND           
!                        S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR STEP).  
!           KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE).         
!                                                                       
         DO 390 LL = 1, M                                               
            L = M - LL                                                  
!        ...EXIT                                                        
            IF (L .EQ. 0) GO TO 400                                     
            TEST = DABS(S(L)) + DABS(S(L+1))                            
            ZTEST = TEST + DABS(E(L))                                   
            IF (ZTEST .NE. TEST) GO TO 380                              
               E(L) = 0.0D0                                             
!        ......EXIT                                                     
               GO TO 400                                                
  380       CONTINUE                                                    
  390    CONTINUE                                                       
  400    CONTINUE                                                       
         IF (L .NE. M - 1) GO TO 410                                    
            KASE = 4                                                    
         GO TO 480                                                      
  410    CONTINUE                                                       
            LP1 = L + 1                                                 
            MP1 = M + 1                                                 
            DO 430 LLS = LP1, MP1                                       
               LS = M - LLS + LP1                                       
!           ...EXIT                                                     
               IF (LS .EQ. L) GO TO 440                                 
               TEST = 0.0D0                                             
               IF (LS .NE. M) TEST = TEST + DABS(E(LS))                 
               IF (LS .NE. L + 1) TEST = TEST + DABS(E(LS-1))           
               ZTEST = TEST + DABS(S(LS))                               
               IF (ZTEST .NE. TEST) GO TO 420                           
                  S(LS) = 0.0D0                                         
!           ......EXIT                                                  
                  GO TO 440                                             
  420          CONTINUE                                                 
  430       CONTINUE                                                    
  440       CONTINUE                                                    
            IF (LS .NE. L) GO TO 450                                    
               KASE = 3                                                 
            GO TO 470                                                   
  450       CONTINUE                                                    
            IF (LS .NE. M) GO TO 460                                    
               KASE = 1                                                 
            GO TO 470                                                   
  460       CONTINUE                                                    
               KASE = 2                                                 
               L = LS                                                   
  470       CONTINUE                                                    
  480    CONTINUE                                                       
         L = L + 1                                                      
         !                                                                       
         !        PERFORM THE TASK INDICATED BY KASE.                            
         !                                                                       
         select case (KASE)
         case (1)
           go to 490
         case (2)
           go to 520
         case (3)
           go to 540
         case (4)
           go to 570
         end select
!                                                                       
!        DEFLATE NEGLIGIBLE S(M).                                       
!                                                                       
  490    CONTINUE                                                       
            MM1 = M - 1                                                 
            F = E(M-1)                                                  
            E(M-1) = 0.0D0                                              
            DO 510 KK = L, MM1                                          
               K = MM1 - KK + L                                         
               T1 = S(K)                                                
               CALL DROTG(T1,F,CS,SN)                                   
               S(K) = T1                                                
               IF (K .EQ. L) GO TO 500                                  
                  F = -SN*E(K-1)                                        
                  E(K-1) = CS*E(K-1)                                    
  500          CONTINUE                                                 
               IF (WANTV) CALL DROT(P,V(1,K),1,V(1,M),1,CS,SN)          
  510       CONTINUE                                                    
         GO TO 610                                                      
!                                                                       
!        SPLIT AT NEGLIGIBLE S(L).                                      
!                                                                       
  520    CONTINUE                                                       
            F = E(L-1)                                                  
            E(L-1) = 0.0D0                                              
            DO 530 K = L, M                                             
               T1 = S(K)                                                
               CALL DROTG(T1,F,CS,SN)                                   
               S(K) = T1                                                
               F = -SN*E(K)                                             
               E(K) = CS*E(K)                                           
               IF (WANTU) CALL DROT(N,U(1,K),1,U(1,L-1),1,CS,SN)        
  530       CONTINUE                                                    
         GO TO 610                                                      
!                                                                       
!        PERFORM ONE QR STEP.                                           
!                                                                       
  540    CONTINUE                                                       
!                                                                       
!           CALCULATE THE SHIFT.                                        
!                                                                       
            SCALE = DMAX1(DABS(S(M)),DABS(S(M-1)),DABS(E(M-1)),          &
                          DABS(S(L)),DABS(E(L)))                        
            SM = S(M)/SCALE                                             
            SMM1 = S(M-1)/SCALE                                         
            EMM1 = E(M-1)/SCALE                                         
            SL = S(L)/SCALE                                             
            EL = E(L)/SCALE                                             
            B = ((SMM1 + SM)*(SMM1 - SM) + EMM1**2)/2.0D0               
            C = (SM*EMM1)**2                                            
            SHIFT = 0.0D0                                               
            IF (B .EQ. 0.0D0 .AND. C .EQ. 0.0D0) GO TO 550              
               SHIFT = DSQRT(B**2+C)                                    
               IF (B .LT. 0.0D0) SHIFT = -SHIFT                         
               SHIFT = C/(B + SHIFT)                                    
  550       CONTINUE                                                    
            F = (SL + SM)*(SL - SM) - SHIFT                             
            G = SL*EL                                                   
!                                                                       
!           CHASE ZEROS.                                                
!                                                                       
            MM1 = M - 1                                                 
            DO 560 K = L, MM1                                           
               CALL DROTG(F,G,CS,SN)                                    
               IF (K .NE. L) E(K-1) = F                                 
               F = CS*S(K) + SN*E(K)                                    
               E(K) = CS*E(K) - SN*S(K)                                 
               G = SN*S(K+1)                                            
               S(K+1) = CS*S(K+1)                                       
               IF (WANTV) CALL DROT(P,V(1,K),1,V(1,K+1),1,CS,SN)        
               CALL DROTG(F,G,CS,SN)                                    
               S(K) = F                                                 
               F = CS*E(K) + SN*S(K+1)                                  
               S(K+1) = -SN*E(K) + CS*S(K+1)                            
               G = SN*E(K+1)                                            
               E(K+1) = CS*E(K+1)                                       
               IF (WANTU .AND. K .LT. N)                                 &
                  CALL DROT(N,U(1,K),1,U(1,K+1),1,CS,SN)                
  560       CONTINUE                                                    
            E(M-1) = F                                                  
            ITER = ITER + 1                                             
         GO TO 610                                                      
!                                                                       
!        CONVERGENCE.                                                   
!                                                                       
  570    CONTINUE                                                       
!                                                                       
!           MAKE THE SINGULAR VALUE  POSITIVE.                          
!                                                                       
            IF (S(L) .GE. 0.0D0) GO TO 580                              
               S(L) = -S(L)                                             
               IF (WANTV) CALL DSCAL(P,-1.0D0,V(1,L),1)                 
  580       CONTINUE                                                    
!                                                                       
!           ORDER THE SINGULAR VALUE.                                   
!                                                                       
  590       IF (L .EQ. MM) GO TO 600                                    
!           ...EXIT                                                     
               IF (S(L) .GE. S(L+1)) GO TO 600                          
               T = S(L)                                                 
               S(L) = S(L+1)                                            
               S(L+1) = T                                               
               IF (WANTV .AND. L .LT. P)                                 &
                  CALL DSWAP(P,V(1,L),1,V(1,L+1),1)                     
               IF (WANTU .AND. L .LT. N)                                 &
                  CALL DSWAP(N,U(1,L),1,U(1,L+1),1)                     
               L = L + 1                                                
            GO TO 590                                                   
  600       CONTINUE                                                    
            ITER = 0                                                    
            M = M - 1                                                   
  610    CONTINUE                                                       
      GO TO 360                                                         
  620 CONTINUE                                                          
      RETURN                                                            
      END                                                               
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
