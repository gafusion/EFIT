      SUBROUTINE r8poisd2 (MR,NR,ISTAG,BA,BB,BC,Q,IDIMQ,B,W,D,TCOS,P)
C
C     SUBROUTINE TO SOLVE POISSON'S EQUATION FOR DIRICHLET BOUNDARY
C     CONDITIONS.
C
C     ISTAG = 1 IF THE LAST DIAGONAL BLOCK IS THE MATRIX A.
C     ISTAG = 2 IF THE LAST DIAGONAL BLOCK IS THE MATRIX A+I.
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER nr,istag,idimq,mr,m,n,jsh,ip,ipstor,kr,irreg,jstsav,i
      INTEGER lr,nun,jst,jsp,l,nodd,j,jm1,jp1,jm2,jp2,jm3,jp3
      INTEGER noddpr,ip1,krpi,ideg,jdeg
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 ba,bb,bc,q,b,w,d,tcos,p,fi,t
!============
      DIMENSION       Q(IDIMQ,1) ,BA(1)      ,BB(1)      ,BC(1)      ,
     1                TCOS(1)    ,B(1)       ,D(1)       ,W(1)       ,
     2                P(1)
      M = MR
      N = NR
      JSH = 0
      FI = 1._r8/real(ISTAG,r8)
      IP = -M
      IPSTOR = 0
      GO TO (101,102),ISTAG
  101 KR = 0
      IRREG = 1
      IF (N .GT. 1) GO TO 106
      TCOS(1) = 0._r8
      GO TO 103
  102 KR = 1
      JSTSAV = 1
      IRREG = 2
      IF (N .GT. 1) GO TO 106
      TCOS(1) = -1._r8
  103 DO 104 I=1,M
         B(I) = Q(I,1)
  104 CONTINUE
      CALL r8trix (1,0,M,BA,BB,BC,B,TCOS,D,W)
      DO 105 I=1,M
         Q(I,1) = B(I)
  105 CONTINUE
      GO TO 183
  106 LR = 0
      DO 107 I=1,M
         P(I) = 0._r8
  107 CONTINUE
      NUN = N
      JST = 1
      JSP = N
C
C     IRREG = 1 WHEN NO IRREGULARITIES HAVE OCCURRED, OTHERWISE IT IS 2.
C
  108 L = 2*JST
      NODD = 2-2*((NUN+1)/2)+NUN
C
C     NODD = 1 WHEN NUN IS ODD, OTHERWISE IT IS 2.
C
      GO TO (110,109),NODD
  109 JSP = JSP-L
      GO TO 111
  110 JSP = JSP-JST
      IF (IRREG .NE. 1) JSP = JSP-L
  111 CONTINUE
C
C     REGULAR REDUCTION
C
      CALL r8cosgen (JST,1,0.5_r8,0.0_r8,TCOS)
      IF (L .GT. JSP) GO TO 118
      DO 117 J=L,JSP,L
         JM1 = J-JSH
         JP1 = J+JSH
         JM2 = J-JST
         JP2 = J+JST
         JM3 = JM2-JSH
         JP3 = JP2+JSH
         IF (JST .NE. 1) GO TO 113
         DO 112 I=1,M
            B(I) = 2._r8*Q(I,J)
            Q(I,J) = Q(I,JM2)+Q(I,JP2)
  112    CONTINUE
         GO TO 115
  113    DO 114 I=1,M
            T = Q(I,J)-Q(I,JM1)-Q(I,JP1)+Q(I,JM2)+Q(I,JP2)
            B(I) = T+Q(I,J)-Q(I,JM3)-Q(I,JP3)
            Q(I,J) = T
  114    CONTINUE
  115    CONTINUE
         CALL r8trix (JST,0,M,BA,BB,BC,B,TCOS,D,W)
         DO 116 I=1,M
            Q(I,J) = Q(I,J)+B(I)
  116    CONTINUE
  117 CONTINUE
C
C     REDUCTION FOR LAST UNKNOWN
C
  118 GO TO (119,136),NODD
  119 GO TO (152,120),IRREG
C
C     ODD NUMBER OF UNKNOWNS
C
  120 JSP = JSP+L
      J = JSP
      JM1 = J-JSH
      JP1 = J+JSH
      JM2 = J-JST
      JP2 = J+JST
      JM3 = JM2-JSH
      GO TO (123,121),ISTAG
  121 CONTINUE
      IF (JST .NE. 1) GO TO 123
      DO 122 I=1,M
         B(I) = Q(I,J)
         Q(I,J) = 0._r8
  122 CONTINUE
      GO TO 130
  123 GO TO (124,126),NODDPR
  124 DO 125 I=1,M
         IP1 = IP+I
         B(I) = .5_r8*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))+P(IP1)+Q(I,J)
  125 CONTINUE
      GO TO 128
  126 DO 127 I=1,M
         B(I) = .5_r8*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
     &        +Q(I,JP2)-Q(I,JP1)+Q(I,J)
  127 CONTINUE
  128 DO 129 I=1,M
         Q(I,J) = .5_r8*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
  129 CONTINUE
  130 CALL r8trix (JST,0,M,BA,BB,BC,B,TCOS,D,W)
      IP = IP+M
      IPSTOR = max(IPSTOR,IP+M)
      DO 131 I=1,M
         IP1 = IP+I
         P(IP1) = Q(I,J)+B(I)
         B(I) = Q(I,JP2)+P(IP1)
  131 CONTINUE
      IF (LR .NE. 0) GO TO 133
      DO 132 I=1,JST
         KRPI = KR+I
         TCOS(KRPI) = TCOS(I)
  132 CONTINUE
      GO TO 134
  133 CONTINUE
      CALL r8cosgen (LR,JSTSAV,0._r8,FI,TCOS(JST+1))
      CALL r8merge (TCOS,0,JST,JST,LR,KR)
  134 CONTINUE
      CALL r8cosgen (KR,JSTSAV,0.0_r8,FI,TCOS)
      CALL r8trix (KR,KR,M,BA,BB,BC,B,TCOS,D,W)
      DO 135 I=1,M
         IP1 = IP+I
         Q(I,J) = Q(I,JM2)+B(I)+P(IP1)
  135 CONTINUE
      LR = KR
      KR = KR+L
      GO TO 152
C
C     EVEN NUMBER OF UNKNOWNS
C
  136 JSP = JSP+L
      J = JSP
      JM1 = J-JSH
      JP1 = J+JSH
      JM2 = J-JST
      JP2 = J+JST
      JM3 = JM2-JSH
      GO TO (137,138),IRREG
  137 CONTINUE
      JSTSAV = JST
      IDEG = JST
      KR = L
      GO TO 139
  138 CALL r8cosgen (KR,JSTSAV,0.0_r8,FI,TCOS)
      CALL r8cosgen (LR,JSTSAV,0.0_r8,FI,TCOS(KR+1))
      IDEG = KR
      KR = KR+JST
  139 IF (JST .NE. 1) GO TO 141
      IRREG = 2
      DO 140 I=1,M
         B(I) = Q(I,J)
         Q(I,J) = Q(I,JM2)
  140 CONTINUE
      GO TO 150
  141 DO 142 I=1,M
         B(I) = Q(I,J)+.5_r8*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  142 CONTINUE
      GO TO (143,145),IRREG
  143 DO 144 I=1,M
         Q(I,J) = Q(I,JM2)+.5_r8*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
  144 CONTINUE
      IRREG = 2
      GO TO 150
  145 CONTINUE
      GO TO (146,148),NODDPR
  146 DO 147 I=1,M
         IP1 = IP+I
         Q(I,J) = Q(I,JM2)+P(IP1)
  147 CONTINUE
      IP = IP-M
      GO TO 150
  148 DO 149 I=1,M
         Q(I,J) = Q(I,JM2)+Q(I,J)-Q(I,JM1)
  149 CONTINUE
  150 CALL r8trix (IDEG,LR,M,BA,BB,BC,B,TCOS,D,W)
      DO 151 I=1,M
         Q(I,J) = Q(I,J)+B(I)
  151 CONTINUE
  152 NUN = NUN/2
      NODDPR = NODD
      JSH = JST
      JST = 2*JST
      IF (NUN .GE. 2) GO TO 108
C
C     START SOLUTION.
C
      J = JSP
      DO 153 I=1,M
         B(I) = Q(I,J)
  153 CONTINUE
      GO TO (154,155),IRREG
  154 CONTINUE
      CALL r8cosgen (JST,1,0.5_r8,0.0_r8,TCOS)
      IDEG = JST
      GO TO 156
  155 KR = LR+JST
      CALL r8cosgen (KR,JSTSAV,0.0_r8,FI,TCOS)
      CALL r8cosgen (LR,JSTSAV,0.0_r8,FI,TCOS(KR+1))
      IDEG = KR
  156 CONTINUE
      CALL r8trix (IDEG,LR,M,BA,BB,BC,B,TCOS,D,W)
      JM1 = J-JSH
      JP1 = J+JSH
      GO TO (157,159),IRREG
  157 DO 158 I=1,M
         Q(I,J) = .5_r8*(Q(I,J)-Q(I,JM1)-Q(I,JP1))+B(I)
  158 CONTINUE
      GO TO 164
  159 GO TO (160,162),NODDPR
  160 DO 161 I=1,M
         IP1 = IP+I
         Q(I,J) = P(IP1)+B(I)
  161 CONTINUE
      IP = IP-M
      GO TO 164
  162 DO 163 I=1,M
         Q(I,J) = Q(I,J)-Q(I,JM1)+B(I)
  163 CONTINUE
  164 CONTINUE
C
C     START BACK SUBSTITUTION.
C
      JST = JST/2
      JSH = JST/2
      NUN = 2*NUN
      IF (NUN .GT. N) GO TO 183
      DO 182 J=JST,N,L
         JM1 = J-JSH
         JP1 = J+JSH
         JM2 = J-JST
         JP2 = J+JST
         IF (J .GT. JST) GO TO 166
         DO 165 I=1,M
            B(I) = Q(I,J)+Q(I,JP2)
  165    CONTINUE
         GO TO 170
  166    IF (JP2 .LE. N) GO TO 168
         DO 167 I=1,M
            B(I) = Q(I,J)+Q(I,JM2)
  167    CONTINUE
         IF (JST .LT. JSTSAV) IRREG = 1
         GO TO (170,171),IRREG
  168    DO 169 I=1,M
            B(I) = Q(I,J)+Q(I,JM2)+Q(I,JP2)
  169    CONTINUE
  170    CONTINUE
         CALL r8cosgen (JST,1,0.5_r8,0.0_r8,TCOS)
         IDEG = JST
         JDEG = 0
         GO TO 172
  171    IF (J+L .GT. N) LR = LR-JST
         KR = JST+LR
         CALL r8cosgen (KR,JSTSAV,0.0_r8,FI,TCOS)
         CALL r8cosgen (LR,JSTSAV,0.0_r8,FI,TCOS(KR+1))
         IDEG = KR
         JDEG = LR
  172    CONTINUE
         CALL r8trix (IDEG,JDEG,M,BA,BB,BC,B,TCOS,D,W)
         IF (JST .GT. 1) GO TO 174
         DO 173 I=1,M
            Q(I,J) = B(I)
  173    CONTINUE
         GO TO 182
  174    IF (JP2 .GT. N) GO TO 177
  175    DO 176 I=1,M
            Q(I,J) = .5_r8*(Q(I,J)-Q(I,JM1)-Q(I,JP1))+B(I)
  176    CONTINUE
         GO TO 182
  177    GO TO (175,178),IRREG
  178    IF (J+JSH .GT. N) GO TO 180
         DO 179 I=1,M
            IP1 = IP+I
            Q(I,J) = B(I)+P(IP1)
  179    CONTINUE
         IP = IP-M
         GO TO 182
  180    DO 181 I=1,M
            Q(I,J) = B(I)+Q(I,J)-Q(I,JM1)
  181    CONTINUE
  182 CONTINUE
      L = L/2
      GO TO 164
  183 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR P VECTORS.
C
      W(1) = IPSTOR
      RETURN
      END
! 06Jul2005 fgtok -s r8_precision.sub fish_sub.txt "r8con.csh conversion"
! 06Jul2005 fgtok
