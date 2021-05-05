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
