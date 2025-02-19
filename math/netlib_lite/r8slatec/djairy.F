*DECK DJAIRY
      SUBROUTINE DJAIRY (X, RX, C, AI, DAI)
C***BEGIN PROLOGUE  DJAIRY
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBESJ and DBESY
C***LIBRARY   SLATEC
C***TYPE      REAL*8 (JAIRY-S, DJAIRY-D)
C***AUTHOR  Amos, D. E., (SNLA)
C           Daniel, S. L., (SNLA)
C           Weston, M. K., (SNLA)
C***DESCRIPTION
C
C                  DJAIRY computes the Airy function AI(X)
C                   and its derivative DAI(X) for DASYJY
C
C                                   INPUT
C
C         X - Argument, computed by DASYJY, X unrestricted
C        RX - RX=SQRT(ABS(X)), computed by DASYJY
C         C - C=2.*(ABS(X)**1.5)/3., computed by DASYJY
C
C                                  OUTPUT
C
C        AI - Value of function AI(X)
C       DAI - Value of the derivative DAI(X)
C
C***SEE ALSO  DBESJ, DBESY
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   750101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891009  Removed unreferenced variable.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910408  Updated the AUTHOR section.  (WRB)
C***END PROLOGUE  DJAIRY
C
      IMPLICIT NONE
      integer, parameter :: r8 = selected_real_kind(12,100)
      INTEGER I, J, M1, M1D, M2, M2D, M3, M3D, M4, M4D, N1, N1D, N2,
     1 N2D, N3, N3D, N4, N4D
      REAL*8 A,AI,AJN,AJP,AK1,AK2,AK3,B,C,CCV,CON2,
     1 CON3, CON4, CON5, CV, DA, DAI, DAJN, DAJP, DAK1, DAK2, DAK3,
     2 DB, EC, E1, E2, FPI12, F1, F2, RTRX, RX, SCV, T, TEMP1, TEMP2,
     3 TT, X
      DIMENSION AJP(19), AJN(19), A(15), B(15)
      DIMENSION AK1(14), AK2(23), AK3(14)
      DIMENSION DAJP(19), DAJN(19), DA(15), DB(15)
      DIMENSION DAK1(14), DAK2(24), DAK3(14)
      SAVE N1, N2, N3, N4, M1, M2, M3, M4, FPI12, CON2, CON3,
     1 CON4, CON5, AK1, AK2, AK3, AJP, AJN, A, B,
     2 N1D, N2D, N3D, N4D, M1D, M2D, M3D, M4D, DAK1, DAK2, DAK3,
     3 DAJP, DAJN, DA, DB
      DATA N1,N2,N3,N4/14,23,19,15/
      DATA M1,M2,M3,M4/12,21,17,13/
      DATA FPI12,CON2,CON3,CON4,CON5/
     1 1.30899693899575E00_r8, 5.03154716196777E00_r8,
     >      3.80004589867293E-01_r8,
     2 8.33333333333333E-01_r8, 8.66025403784439E-01_r8/
      DATA AK1(1), AK1(2), AK1(3), AK1(4), AK1(5), AK1(6), AK1(7),
     1     AK1(8), AK1(9), AK1(10),AK1(11),AK1(12),AK1(13),
     2     AK1(14)         / 2.20423090987793E-01_r8,-
     >     1.25290242787700E-01_r8,
     3 1.03881163359194E-02_r8, 8.22844152006343E-04_r8,-
     >     2.34614345891226E-04_r8,
     4 1.63824280172116E-05_r8, 3.06902589573189E-07_r8,-
     >     1.29621999359332E-07_r8,
     5 8.22908158823668E-09_r8, 1.53963968623298E-11_r8,-
     >     3.39165465615682E-11_r8,
     6 2.03253257423626E-12_r8,-1.10679546097884E-14_r8,-
     >     5.16169497785080E-15_r8/
      DATA AK2(1), AK2(2), AK2(3), AK2(4), AK2(5), AK2(6), AK2(7),
     1     AK2(8), AK2(9), AK2(10),AK2(11),AK2(12),AK2(13),AK2(14),
     2     AK2(15),AK2(16),AK2(17),AK2(18),AK2(19),AK2(20),AK2(21),
     3     AK2(22),AK2(23) / 2.74366150869598E-01_r8,
     >      5.39790969736903E-03_r8,
     4-1.57339220621190E-03_r8, 4.27427528248750E-04_r8,-
     >     1.12124917399925E-04_r8,
     5 2.88763171318904E-05_r8,-7.36804225370554E-06_r8,
     >      1.87290209741024E-06_r8,
     6-4.75892793962291E-07_r8, 1.21130416955909E-07_r8,-
     >     3.09245374270614E-08_r8,
     7 7.92454705282654E-09_r8,-2.03902447167914E-09_r8,
     >      5.26863056595742E-10_r8,
     8-1.36704767639569E-10_r8, 3.56141039013708E-11_r8,-
     >     9.31388296548430E-12_r8,
     9 2.44464450473635E-12_r8,-6.43840261990955E-13_r8,
     >      1.70106030559349E-13_r8,
     1-4.50760104503281E-14_r8, 1.19774799164811E-14_r8,-
     >     3.19077040865066E-15_r8/
      DATA AK3(1), AK3(2), AK3(3), AK3(4), AK3(5), AK3(6), AK3(7),
     1     AK3(8), AK3(9), AK3(10),AK3(11),AK3(12),AK3(13),
     2     AK3(14)         / 2.80271447340791E-01_r8,-
     >     1.78127042844379E-03_r8,
     3 4.03422579628999E-05_r8,-1.63249965269003E-06_r8,
     >      9.21181482476768E-08_r8,
     4-6.52294330229155E-09_r8, 5.47138404576546E-10_r8,-
     >     5.24408251800260E-11_r8,
     5 5.60477904117209E-12_r8,-6.56375244639313E-13_r8,
     >      8.31285761966247E-14_r8,
     6-1.12705134691063E-14_r8, 1.62267976598129E-15_r8,-
     >     2.46480324312426E-16_r8/
      DATA AJP(1), AJP(2), AJP(3), AJP(4), AJP(5), AJP(6), AJP(7),
     1     AJP(8), AJP(9), AJP(10),AJP(11),AJP(12),AJP(13),AJP(14),
     2     AJP(15),AJP(16),AJP(17),AJP(18),
     3     AJP(19)         / 7.78952966437581E-02_r8,-
     >     1.84356363456801E-01_r8,
     4 3.01412605216174E-02_r8, 3.05342724277608E-02_r8,-
     >     4.95424702513079E-03_r8,
     5-1.72749552563952E-03_r8, 2.43137637839190E-04_r8,
     >      5.04564777517082E-05_r8,
     6-6.16316582695208E-06_r8,-9.03986745510768E-07_r8,
     >      9.70243778355884E-08_r8,
     7 1.09639453305205E-08_r8,-1.04716330588766E-09_r8,-
     >     9.60359441344646E-11_r8,
     8 8.25358789454134E-12_r8, 6.36123439018768E-13_r8,-
     >     4.96629614116015E-14_r8,
     9-3.29810288929615E-15_r8, 2.35798252031104E-16_r8/
      DATA AJN(1), AJN(2), AJN(3), AJN(4), AJN(5), AJN(6), AJN(7),
     1     AJN(8), AJN(9), AJN(10),AJN(11),AJN(12),AJN(13),AJN(14),
     2     AJN(15),AJN(16),AJN(17),AJN(18),
     3     AJN(19)         / 3.80497887617242E-02_r8,-
     >     2.45319541845546E-01_r8,
     4 1.65820623702696E-01_r8, 7.49330045818789E-02_r8,-
     >     2.63476288106641E-02_r8,
     5-5.92535597304981E-03_r8, 1.44744409589804E-03_r8,
     >      2.18311831322215E-04_r8,
     6-4.10662077680304E-05_r8,-4.66874994171766E-06_r8,
     >      7.15218807277160E-07_r8,
     7 6.52964770854633E-08_r8,-8.44284027565946E-09_r8,-
     >     6.44186158976978E-10_r8,
     8 7.20802286505285E-11_r8, 4.72465431717846E-12_r8,-
     >     4.66022632547045E-13_r8,
     9-2.67762710389189E-14_r8, 2.36161316570019E-15_r8/
      DATA A(1),   A(2),   A(3),   A(4),   A(5),   A(6),   A(7),
     1     A(8),   A(9),   A(10),  A(11),  A(12),  A(13),  A(14),
     2     A(15)           / 4.90275424742791E-01_r8,
     >      1.57647277946204E-03_r8,
     3-9.66195963140306E-05_r8, 1.35916080268815E-07_r8,
     >      2.98157342654859E-07_r8,
     4-1.86824767559979E-08_r8,-1.03685737667141E-09_r8,
     >      3.28660818434328E-10_r8,
     5-2.57091410632780E-11_r8,-2.32357655300677E-12_r8,
     >      9.57523279048255E-13_r8,
     6-1.20340828049719E-13_r8,-2.90907716770715E-15_r8,
     >      4.55656454580149E-15_r8,
     7-9.99003874810259E-16_r8/
      DATA B(1),   B(2),   B(3),   B(4),   B(5),   B(6),   B(7),
     1     B(8),   B(9),   B(10),  B(11),  B(12),  B(13),  B(14),
     2     B(15)           / 2.78593552803079E-01_r8,-
     >     3.52915691882584E-03_r8,
     3-2.31149677384994E-05_r8, 4.71317842263560E-06_r8,-
     >     1.12415907931333E-07_r8,
     4-2.00100301184339E-08_r8, 2.60948075302193E-09_r8,-
     >     3.55098136101216E-11_r8,
     5-3.50849978423875E-11_r8, 5.83007187954202E-12_r8,-
     >     2.04644828753326E-13_r8,
     6-1.10529179476742E-13_r8, 2.87724778038775E-14_r8,-
     >     2.88205111009939E-15_r8,
     7-3.32656311696166E-16_r8/
      DATA N1D,N2D,N3D,N4D/14,24,19,15/
      DATA M1D,M2D,M3D,M4D/12,22,17,13/
      DATA DAK1(1), DAK1(2), DAK1(3), DAK1(4), DAK1(5), DAK1(6),
     1     DAK1(7), DAK1(8), DAK1(9), DAK1(10),DAK1(11),DAK1(12),
     2    DAK1(13),DAK1(14)/ 2.04567842307887E-01_r8,-
     >     6.61322739905664E-02_r8,
     3-8.49845800989287E-03_r8, 3.12183491556289E-03_r8,-
     >     2.70016489829432E-04_r8,
     4-6.35636298679387E-06_r8, 3.02397712409509E-06_r8,-
     >     2.18311195330088E-07_r8,
     5-5.36194289332826E-10_r8, 1.13098035622310E-09_r8,-
     >     7.43023834629073E-11_r8,
     6 4.28804170826891E-13_r8, 2.23810925754539E-13_r8,-
     >     1.39140135641182E-14_r8/
      DATA DAK2(1), DAK2(2), DAK2(3), DAK2(4), DAK2(5), DAK2(6),
     1     DAK2(7), DAK2(8), DAK2(9), DAK2(10),DAK2(11),DAK2(12),
     2     DAK2(13),DAK2(14),DAK2(15),DAK2(16),DAK2(17),DAK2(18),
     3     DAK2(19),DAK2(20),DAK2(21),DAK2(22),DAK2(23),
     4     DAK2(24)        / 2.93332343883230E-01_r8,-
     >     8.06196784743112E-03_r8,
     5 2.42540172333140E-03_r8,-6.82297548850235E-04_r8,
     >      1.85786427751181E-04_r8,
     6-4.97457447684059E-05_r8, 1.32090681239497E-05_r8,-
     >     3.49528240444943E-06_r8,
     7 9.24362451078835E-07_r8,-2.44732671521867E-07_r8,
     >      6.49307837648910E-08_r8,
     8-1.72717621501538E-08_r8, 4.60725763604656E-09_r8,-
     >     1.23249055291550E-09_r8,
     9 3.30620409488102E-10_r8,-8.89252099772401E-11_r8,
     >      2.39773319878298E-11_r8,
     1-6.48013921153450E-12_r8, 1.75510132023731E-12_r8,-
     >     4.76303829833637E-13_r8,
     2 1.29498241100810E-13_r8,-3.52679622210430E-14_r8,
     >      9.62005151585923E-15_r8,
     3-2.62786914342292E-15_r8/
      DATA DAK3(1), DAK3(2), DAK3(3), DAK3(4), DAK3(5), DAK3(6),
     1     DAK3(7), DAK3(8), DAK3(9), DAK3(10),DAK3(11),DAK3(12),
     2    DAK3(13),DAK3(14)/ 2.84675828811349E-01_r8,
     >      2.53073072619080E-03_r8,
     3-4.83481130337976E-05_r8, 1.84907283946343E-06_r8,-
     >     1.01418491178576E-07_r8,
     4 7.05925634457153E-09_r8,-5.85325291400382E-10_r8,
     >      5.56357688831339E-11_r8,
     5-5.90889094779500E-12_r8, 6.88574353784436E-13_r8,-
     >     8.68588256452194E-14_r8,
     6 1.17374762617213E-14_r8,-1.68523146510923E-15_r8,
     >      2.55374773097056E-16_r8/
      DATA DAJP(1), DAJP(2), DAJP(3), DAJP(4), DAJP(5), DAJP(6),
     1     DAJP(7), DAJP(8), DAJP(9), DAJP(10),DAJP(11),DAJP(12),
     2     DAJP(13),DAJP(14),DAJP(15),DAJP(16),DAJP(17),DAJP(18),
     3     DAJP(19)        / 6.53219131311457E-02_r8,-
     >     1.20262933688823E-01_r8,
     4 9.78010236263823E-03_r8, 1.67948429230505E-02_r8,-
     >     1.97146140182132E-03_r8,
     5-8.45560295098867E-04_r8, 9.42889620701976E-05_r8,
     >      2.25827860945475E-05_r8,
     6-2.29067870915987E-06_r8,-3.76343991136919E-07_r8,
     >      3.45663933559565E-08_r8,
     7 4.29611332003007E-09_r8,-3.58673691214989E-10_r8,-
     >     3.57245881361895E-11_r8,
     8 2.72696091066336E-12_r8, 2.26120653095771E-13_r8,-
     >     1.58763205238303E-14_r8,
     9-1.12604374485125E-15_r8, 7.31327529515367E-17_r8/
      DATA DAJN(1), DAJN(2), DAJN(3), DAJN(4), DAJN(5), DAJN(6),
     1     DAJN(7), DAJN(8), DAJN(9), DAJN(10),DAJN(11),DAJN(12),
     2     DAJN(13),DAJN(14),DAJN(15),DAJN(16),DAJN(17),DAJN(18),
     3     DAJN(19)        / 1.08594539632967E-02_r8,
     >      8.53313194857091E-02_r8,
     4-3.15277068113058E-01_r8,-8.78420725294257E-02_r8,
     >      5.53251906976048E-02_r8,
     5 9.41674060503241E-03_r8,-3.32187026018996E-03_r8,-
     >     4.11157343156826E-04_r8,
     6 1.01297326891346E-04_r8, 9.87633682208396E-06_r8,-
     >     1.87312969812393E-06_r8,
     7-1.50798500131468E-07_r8, 2.32687669525394E-08_r8,
     >      1.59599917419225E-09_r8,
     8-2.07665922668385E-10_r8,-1.24103350500302E-11_r8,
     >      1.39631765331043E-12_r8,
     9 7.39400971155740E-14_r8,-7.32887475627500E-15_r8/
      DATA DA(1),  DA(2),  DA(3),  DA(4),  DA(5),  DA(6),  DA(7),
     1     DA(8),  DA(9),  DA(10), DA(11), DA(12), DA(13), DA(14),
     2     DA(15)          / 4.91627321104601E-01_r8,
     >      3.11164930427489E-03_r8,
     3 8.23140762854081E-05_r8,-4.61769776172142E-06_r8,-
     >     6.13158880534626E-08_r8,
     4 2.87295804656520E-08_r8,-1.81959715372117E-09_r8,-
     >     1.44752826642035E-10_r8,
     5 4.53724043420422E-11_r8,-3.99655065847223E-12_r8,-
     >     3.24089119830323E-13_r8,
     6 1.62098952568741E-13_r8,-2.40765247974057E-14_r8,
     >      1.69384811284491E-16_r8,
     7 8.17900786477396E-16_r8/
      DATA DB(1),  DB(2),  DB(3),  DB(4),  DB(5),  DB(6),  DB(7),
     1     DB(8),  DB(9),  DB(10), DB(11), DB(12), DB(13), DB(14),
     2     DB(15)          /-2.77571356944231E-01_r8,
     >      4.44212833419920E-03_r8,
     3-8.42328522190089E-05_r8,-2.58040318418710E-06_r8,
     >      3.42389720217621E-07_r8,
     4-6.24286894709776E-09_r8,-2.36377836844577E-09_r8,
     >      3.16991042656673E-10_r8,
     5-4.40995691658191E-12_r8,-5.18674221093575E-12_r8,
     >      9.64874015137022E-13_r8,
     6-4.90190576608710E-14_r8,-1.77253430678112E-14_r8,
     >      5.55950610442662E-15_r8,
     7-7.11793337579530E-16_r8/
C***FIRST EXECUTABLE STATEMENT  DJAIRY
      IF (X.LT.0.0_r8) GO TO 90
      IF (C.GT.5.0_r8) GO TO 60
      IF (X.GT.1.20_r8) GO TO 30
      T = (X+X-1.2_r8)*CON4
      TT = T + T
      J = N1
      F1 = AK1(J)
      F2 = 0.0_r8
      DO 10 I=1,M1
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + AK1(J)
        F2 = TEMP1
   10 CONTINUE
      AI = T*F1 - F2 + AK1(1)
C
      J = N1D
      F1 = DAK1(J)
      F2 = 0.0_r8
      DO 20 I=1,M1D
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + DAK1(J)
        F2 = TEMP1
   20 CONTINUE
      DAI = -(T*F1-F2+DAK1(1))
      RETURN
C
   30 CONTINUE
      T = (X+X-CON2)*CON3
      TT = T + T
      J = N2
      F1 = AK2(J)
      F2 = 0.0_r8
      DO 40 I=1,M2
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + AK2(J)
        F2 = TEMP1
   40 CONTINUE
      RTRX = SQRT(RX)
      EC = EXP(-C)
      AI = EC*(T*F1-F2+AK2(1))/RTRX
      J = N2D
      F1 = DAK2(J)
      F2 = 0.0_r8
      DO 50 I=1,M2D
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + DAK2(J)
        F2 = TEMP1
   50 CONTINUE
      DAI = -EC*(T*F1-F2+DAK2(1))*RTRX
      RETURN
C
   60 CONTINUE
      T = 10.0_r8/C - 1.0_r8
      TT = T + T
      J = N1
      F1 = AK3(J)
      F2 = 0.0_r8
      DO 70 I=1,M1
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + AK3(J)
        F2 = TEMP1
   70 CONTINUE
      RTRX = SQRT(RX)
      EC = EXP(-C)
      AI = EC*(T*F1-F2+AK3(1))/RTRX
      J = N1D
      F1 = DAK3(J)
      F2 = 0.0_r8
      DO 80 I=1,M1D
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + DAK3(J)
        F2 = TEMP1
   80 CONTINUE
      DAI = -RTRX*EC*(T*F1-F2+DAK3(1))
      RETURN
C
   90 CONTINUE
      IF (C.GT.5.0_r8) GO TO 120
      T = 0.4_r8*C - 1.0_r8
      TT = T + T
      J = N3
      F1 = AJP(J)
      E1 = AJN(J)
      F2 = 0.0_r8
      E2 = 0.0_r8
      DO 100 I=1,M3
        J = J - 1
        TEMP1 = F1
        TEMP2 = E1
        F1 = TT*F1 - F2 + AJP(J)
        E1 = TT*E1 - E2 + AJN(J)
        F2 = TEMP1
        E2 = TEMP2
  100 CONTINUE
      AI = (T*E1-E2+AJN(1)) - X*(T*F1-F2+AJP(1))
      J = N3D
      F1 = DAJP(J)
      E1 = DAJN(J)
      F2 = 0.0_r8
      E2 = 0.0_r8
      DO 110 I=1,M3D
        J = J - 1
        TEMP1 = F1
        TEMP2 = E1
        F1 = TT*F1 - F2 + DAJP(J)
        E1 = TT*E1 - E2 + DAJN(J)
        F2 = TEMP1
        E2 = TEMP2
  110 CONTINUE
      DAI = X*X*(T*F1-F2+DAJP(1)) + (T*E1-E2+DAJN(1))
      RETURN
C
  120 CONTINUE
      T = 10.0_r8/C - 1.0_r8
      TT = T + T
      J = N4
      F1 = A(J)
      E1 = B(J)
      F2 = 0.0_r8
      E2 = 0.0_r8
      DO 130 I=1,M4
        J = J - 1
        TEMP1 = F1
        TEMP2 = E1
        F1 = TT*F1 - F2 + A(J)
        E1 = TT*E1 - E2 + B(J)
        F2 = TEMP1
        E2 = TEMP2
  130 CONTINUE
      TEMP1 = T*F1 - F2 + A(1)
      TEMP2 = T*E1 - E2 + B(1)
      RTRX = SQRT(RX)
      CV = C - FPI12
      CCV = COS(CV)
      SCV = SIN(CV)
      AI = (TEMP1*CCV-TEMP2*SCV)/RTRX
      J = N4D
      F1 = DA(J)
      E1 = DB(J)
      F2 = 0.0_r8
      E2 = 0.0_r8
      DO 140 I=1,M4D
        J = J - 1
        TEMP1 = F1
        TEMP2 = E1
        F1 = TT*F1 - F2 + DA(J)
        E1 = TT*E1 - E2 + DB(J)
        F2 = TEMP1
        E2 = TEMP2
  140 CONTINUE
      TEMP1 = T*F1 - F2 + DA(1)
      TEMP2 = T*E1 - E2 + DB(1)
      E1 = CCV*CON5 + 0.5_r8*SCV
      E2 = SCV*CON5 - 0.5_r8*CCV
      DAI = (TEMP1*E1-TEMP2*E2)*RTRX
      RETURN
      END
! 22Jun2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 22Jun2000 fgtok
