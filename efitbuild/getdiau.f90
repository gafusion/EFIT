        subroutine getdia(nshot,tim,npts,tavg,ierr,phidia,sigphi)

! **********************************************************************
!   This routine returns the average of the compensated diamagnetic    *
!   fluxes and error in the diamagnetic fluxes.                        * 
!       IERR - an error array of dimension 3. Each element of the      *
!	       array contains the PTDATA error for one of the three    *
!              diamagnetic flux signals. If ierr(j) is different from  *
!              zero, that flux is not included in the average value    *
!              returned.                                               *
!       TAVG - this is the time that the data is averaged              *
!              over.                                                   *
!       PHIDIA - this array contains the conpensated diamagnetic       *
!                flux (mv-sec).                                        *
!       SIGPHI - this array contains the error in the flux value       *
!                (mv-sec).                                             *
!                                                                      *
!   written by A.Kellman  7/23/86                                      *  
!                                                                      *
!       MODIFICATIONS:                                                 *
!	  7/23/86 - DIAMAG01 is not included in the average because    *
!	            BTDOT has not yet been removed so the data is no   *
!	            good yet.                                          * 
!         8/30/88 - Modified DLCOMP so that if the DCOEF file is       *
!                   not found in the DIAMAG_DIR area, then it looks    *
!                   for it in [kellman.dia] area on the USC vax.       *
!         6/19/92 - diamag2 data not used because the diagnostic       *
!                   signal is not reliable.                            *
!        11/21/94 - New calibration factors. New updated linear comp   *
!                   terms including three c-coil terms in addition to  *
!                   previous n=1 coil term, and new nonlinear term     *
!                   for DIA3 only for ECOIL times product of ECOIL     *
!                   and BT. Now looks in SYS$USER:[LAHAYE.DIA]. (Rob   *
!                   LaHaye)                                            *
!       2013/06/25  Increase data arrays dimension, new compensation   *                                           *
!***********************************************************************   
        parameter (ntims=8192)
        dimension diamag(ntims,3),diamagc(ntims,3), &
                  sig(ntims,3),tim(ntims)
        dimension phidia(1),sigphi(1)
        dimension ierr(3),iwght(3)


        ndia=3
        call dlcomp(tim,diamag,diamagc,sig,nshot,npts,idlc,ierr, &
            tavg)
        ierr(1)=1
	ierr(2)=1	!do not use diamag2
        do 50 j=1,3 
        iwght(j)=1
        if (ierr(j).ne.0) then
            iwght(j)=0
            ndia=ndia-1	
        endif
  50    continue

        if (ndia.eq.0) then
            do 70 i=1,npts
                phidia(i)=0.0
                sigphi(i)=0.0
  70        continue
        else

            do 100 i=1,npts
                phidia(i)=0.0
                sigphi(i)=0.0
                do 80 j=1,3
                    phidia(i)=phidia(i)+iwght(j)*diamagc(i,j)
                    sigphi(i)=sigphi(i)+iwght(j)*sig(i,j)
  80            continue
                phidia(i)=phidia(i)/ndia
                sigphi(i)=sigphi(i)/ndia
 100        continue
        endif
        return
        end

!********************************************************************
        SUBROUTINE DLCOMP(TIM,DIAMAG,DIAMAGC,SIG,NSHOT,NPTS,IDLC &
      		,IERR,TAVG)
!
!   MODIFIED 9/5/87 - An incorrect RC/G was found for all three loops,
!       A.Kellman     so the correction was put in for all shots before
!                     the correct RC/G was put into the database. This
!                     was done using the variable FIX(3).
!   Modified 6/1/89 - An additional term has been added to compensate
!                     for the pickup from the n=1 coil. The coefficient
!                     was determined from shot XXXXX.
!      Revised: R. La Haye Nov. 21, 1994.
!        8/9/2004 - Remove n1coil for shot > 108281
!   9/21/2020 - R.S. Changed ifix to int
!       
!vas-oct3,08        include 'expath.inc'
        use expath
        PARAMETER (NTIMS=8192)
        DIMENSION DIAMAG(NTIMS,3),DIAMAGC(NTIMS,3),SIG(NTIMS,3)
        DIMENSION TIM(NTIMS),TTEMP(NTIMS)
        DIMENSION IDAT(NTIMS),YDAT(NTIMS),YTEMP(NTIMS),BT(NTIMS)
        DIMENSION IAR(128),RAR(8400),IERR(3),IADJ(3),BTCOMP(3),EBCOUP(3)
        DIMENSION FIX(3)
        DIMENSION ASCII(11),INT16(11),INT32(11),REAL32(11)

        DIMENSION COUP(26,3)
        INTEGER ASCII
        REAL SOURCE
        CHARACTER*100 FILIN
        CHARACTER*10 POINT(26),IPOINT,DNAME(3),DNAMES(3)
            integer*2 iipoint(5)
            equivalence (iipoint,ipoint)
        DATA SOURCE/'.PLA'/
        DATA DNAME/'DIAMAG1   ','DIAMAG2   ','DIAMAG3   '/
        DATA DNAMES/'DIAMAG1S  ','DIAMAG2S  ','DIAMAG3S  '/
        DATA IADJ/-1,1,1/
        DATA BTCOMP/-0.3,0.4,-0.4/
        DATA EBCOUP/0.00,0.00,-2.0/

        EQUIVALENCE (RAR(21),TTEMP(1))
        if (nshot.lt.85885) then
    		filin = 'dcoef.dat'
            EBCOUP(3)=-2.0
        else if (nshot.lt.108282) then
                filin = 'dcoef95.dat'
            EBCOUP(3)=0.0
        else if (nshot.lt.152400) then
                filin = 'dcoef04.dat'
            EBCOUP(3)=0.0
        else
                filin = 'dcoef13.dat'
            EBCOUP(3)=0.0
        endif
!
        filin=input_dir(1:lindir)//filin
!vas
!        print*,' vasan filr : ',filin

        ITYP=11
        IAR(1)=NPTS
        RAR(1)=-0.050000
        RAR(2)=0.0020005
        NAVG=int((TAVG+.5)/1000./RAR(2))

        NIN=21
        IER=1
        ASCII(1)=0
        INT16(1)=0
        INT32(1)=0
        REAL32(1)=0
        IDLC=0
        SIGMAD=0.01
        SIGMAF=0.1

!
!   GET COMPENSATION COEFFICIENTS
!
        OPEN(UNIT=NIN,ACCESS='SEQUENTIAL', &
        STATUS='OLD',FILE=FILIN,ERR=10)
        GOTO 15

   10   FILIN='VUSC::SYS$USER:[LAHAYE.DIA]dcoef.sav'
        OPEN(UNIT=NIN,ACCESS='SEQUENTIAL', &
        STATUS='OLD',FILE=FILIN,ERR=2000)
 
!
!
   15   DO 21 I=1,26
          READ(NIN,1001,err=2000) POINT(I),(COUP(I,J),J=1,3)
!         WRITE(6,1002) I,POINT(I),(COUP(I,J),J=1,3)
   21   CONTINUE
        CLOSE(UNIT=NIN)

!
!  CALCULATE COEFFICIENTS FOR BT COMPENSATION. THIS IS DONE BY
!  COMPUTING THE RATIO OF THE RESIDUAL TOROIDAL FIELD PICKUP ON THE
!  SLOW DIAMAGNETIC SIGNALS TO THE CURRENT IN THE BCOIL AT THE
!  TIME RIGHT BEFORE THE ECOIL TURNS ON I.E., HEX CODE = 69
!

!   GET TIME OF HEX CODE 69. N69 WILL BE THE INDEX INTO THE ARRAY
!      FOR THE TIME 69.
!   GET BCOIL CURRENT
!   GET SLOW DIAMAGNETIC SIGNALS
!	DO 100 J=1,3
!	IPOINT=DNAMES(J)
!	CALL PTDATA(ITYP,NSHOT,SOURCE,%REF(IPOINT),IDAT,IER,IAR,RAR,
!     *		ASCII,INT16,INT32,REAL32)
!	IF (IER.NE.0) THEN
!		WRITE (6,1012)
!		WRITE (6,1013) IPOINT,NSHOT,IER,RAR(7),RAR(3),IAR(1)
!		TYPE *,' NO BT COMPENSATION DONE FOR ',DNAME(J)
!	ENDIF
!	DO 40 N=1,NPTS
!	DIAMAG(N,J)=(RAR(6)-IDAT(N))*RAR(5)*RAR(4)
!   40	CONTINUE
!	BTCOEF(J)=DIAMAG(N69,J)/BCOIL(N69)
!  100	CONTINUE

!
!   adjust for an error in the rc/g for all three loops. The correct value was
!   measured on 9/4/87 and is correct for all shots after 56600.
!
        if (nshot.lt.56600) then
            fix(1)=1.017
            fix(2)=1.019
            fix(3)=1.035
        else
            fix(1)=1.000
            fix(2)=1.000
            fix(3)=1.000
        endif             
!
!   GET DIAMAGNETIC FLUXES
!
        DO 130 J=1,3
          IERR(J)=0
          IPOINT=DNAME(J)
          CALL PTDATA(ITYP,NSHOT,SOURCE,iIPOINT, IDAT,IER,IAR,RAR, &
              ASCII,INT16,INT32,REAL32)
          IF (IER.NE.0 .AND. IER.NE.2 .AND. IER.NE.4) THEN
              WRITE (6,1012)
              WRITE (6,1013) IPOINT,NSHOT,IER,RAR(1),RAR(7), &
                             TTEMP(NPTS)
              IERR(J)=IER
              GO TO 130
          ENDIF
          IF (IER.EQ.4) THEN
              IF (IAR(2).LT.NPTS) NPTS=IAR(2)
              IF (NPTS.LE.0) GO TO 130
          ENDIF

          DO 125 N=1,NPTS
              DIAMAG(N,J)=(RAR(6)-IDAT(N))*RAR(5)*RAR(4)*IADJ(J)*fix(j)
!        
!   COMPENSATE FOR BT PICKUP (0.1 MV-SEC IS UNCERTAINTY IN BT PICKUP)
!   NOV 21, 1994: NEW CALIB INCREASES BY 5%.
!
              DIAMAGC(N,J)=1.05*DIAMAG(N,J)-BTCOMP(J)
              SIG(N,J)=(DIAMAG(N,J)*SIGMAD)**2+0.1
  125     CONTINUE
        if (navg.gt.1) then
            CALL LOWPASS (DIAMAG(1,J),DIAMAG(1,J),NAVG,NPTS)
            CALL LOWPASS (DIAMAGC(1,J),DIAMAGC(1,J),NAVG,NPTS)
        endif
  130   CONTINUE

!   CREATE TIME ARRAY BASED ON THE DIAMAGNETIC LOOP TIME ARRAY

        IF (NPTS.GT.0) THEN
        DO N=1,NPTS
            TIM(N)=TTEMP(N)
        ENDDO
        ENDIF
!
!   COMPENSATE FOR E AND F-COIL PICKUP
!
        DO 200 I=1,26
        IPOINT=POINT(I)
        IREPLACE=0
        ier=0
        IF (IPOINT.EQ.'ECOILB    ' .or. IPOINT.EQ.'VLOOPBS   ') THEN
            IAR(1)=8180
        ELSE
            IAR(1)=8160
        ENDIF
        IERSAV=0
        IF (IPOINT.EQ.'VLOOPBS   ' .AND. NSHOT.EQ.52400) THEN
            IPOINT='VLOOPB    '
            IREPLACE=1
        ENDIF
        IF ((ISHOT .GE. 83350) .AND. (I .GE. 24)) THEN 
  145   CALL PTDATA(ITYP,NSHOT,SOURCE,iIPOINT ,IDAT,IER,IAR,RAR, &
        ASCII,INT16,INT32,REAL32)
        ELSEIF ((ISHOT .LT. 83350) .AND. (I .LT. 24)) THEN 
  146   CALL PTDATA(ITYP,NSHOT,SOURCE,iIPOINT ,IDAT,IER,IAR,RAR, &
        ASCII,INT16,INT32,REAL32)
        ENDIF
        if (ier.ne.0 .and. ier.ne.2 .and. ier.ne.4) then
            write (6,1012)
            write (6,1013) POINT(I),NSHOT,IER,RAR(1),TTEMP(1)
            idlc=idlc+1
            goto 200
        endif       

   
        IF (IPOINT.EQ.'VLOOPBS   ' .OR. IREPLACE.EQ.1) THEN
            NPTS2=IAR(2)
            TAU=0.125
            SUM=0.0
            DO N=1,NPTS2-1
                IF (TTEMP(N).GT.TIM(NPTS)) GOTO 150
                YTEMP(N)=(RAR(6)-IDAT(N))*RAR(5)*RAR(4)
                SUM=SUM+EXP(TTEMP(N)/TAU)*YTEMP(N)*(TTEMP(N+1)-TTEMP(N))
                YTEMP(N)=YTEMP(N)-EXP(-1.*TTEMP(N)/TAU)/TAU*SUM
            ENDDO

  150       call interp(ttemp,ytemp,N,tim,ydat,npts)

        ELSE IF (IPOINT.EQ.'ECOILB    ' .OR. (ABS(RAR(1)-TTEMP(1)).GT. &
               RAR(3) .OR. ABS(TTEMP(NPTS)-TIM(NPTS)).GT.1.0E-6)) then
            npts2=iar(2)
            do n=1,npts2
                ytemp(n)=(rar(6)-idat(n))*rar(5)*rar(4)
            enddo
            call interp(ttemp,ytemp,iar(2),tim,ydat,npts)
        else
            do n=1,npts
                ydat(n)=(rar(6)-idat(n))*rar(5)*rar(4)
            enddo
        endif 

        IF (ipoint.eq.'ECOILB    ') then
            do n=1,npts
                prod=((ydat(n)/6.6e4)**2)*bt(n)/2.13
                do j=1,3
                    DIAMAGC(N,J)=DIAMAGC(N,J)-EBCOUP(J)*PROD
                    sig(n,j)=sig(n,j)+(ebcoup(j)*prod)**2
                enddo
            enddo
        endif


        DO 160 J=1,3
            const=coup(i,j)*iadj(j)
            DO 155 N=1,NPTS
                DIAMAGC(N,J)=DIAMAGC(N,J)-YDAT(N)*const
                SIG(N,J)=SIG(N,J)+(YDAT(N)*CONST*SIGMAF)**2
  155       CONTINUE
  160   CONTINUE

  200	CONTINUE



        DO 475 J=1,3
        IF (IERR(J).NE.0) THEN
            DO 450 N=1,NPTS
                DIAMAG(N,J)=0.0
                DIAMAGC(N,J)=0.0
  450       CONTINUE
        ENDIF
  475   CONTINUE  

        DO 550 J=1,3
            DO 500 N=1,NPTS
                SIG(N,J)=SQRT(SIG(N,J))
  500       CONTINUE
  550   CONTINUE

 1001   FORMAT(A10,3E12.4)
 1002   FORMAT(' ',I3,2X,A10,3E12.4)
 1004   FORMAT(A10)
 1012   FORMAT(' *********************  ERROR IN PTDATA CALL ********')
 1013   FORMAT(' ',A10,2I6,4F9.4)

        RETURN

 2000   do 2010 i=1,3
        ierr(i)=1
 2010   continue
        write (6,2020)
 2020   format(' error in open or read for file DCOEF.DAT')
        return
        END

!********************************************************************
        SUBROUTINE LOWPASS(XIN,XOUT,NAVG,NPTS)
        PARAMETER (NTIMS=8192)
        REAL XIN(1),XOUT(1),XTMP(NTIMS)

        DO 300 J=1,NPTS
        IND=(NAVG-1)/2
        IF (J.LE.IND) IND=J-1
        IF ((J+IND).GT.NPTS) IND=NPTS-J
        ICNT=0
        XTMP(J)=0

        DO 200 I=J-IND,J+IND
        XTMP(J)=XTMP(J)+XIN(I)
        ICNT=ICNT+1
  200   CONTINUE

        XTMP(J)=XTMP(J)/ICNT

  300   CONTINUE

        DO 400 J=1,NPTS
        XOUT(J)=XTMP(J)
  400   CONTINUE

        RETURN
        END

!*************************************************************
        SUBROUTINE INTERP(XIN,YIN,NIN,XOUT,YOUT,NOUT)
        PARAMETER (NTIMS=8192)
        REAL XIN(NTIMS),YIN(NTIMS),XOUT(NTIMS),YOUT(NTIMS)
        INTEGER NIN,NOUT

        DO 5 I=1,NOUT
        IF (XOUT(I).GE.XIN(1)) GOTO 20
        YOUT(I)=YIN(1)
    5   CONTINUE
        GOTO 500

   20   JLAST=1
   25   DO 100 J=JLAST,NIN-1
        IF ((XOUT(I).GE.XIN(J)).AND.(XOUT(I).LT.XIN(J+1))) GOTO 200
  100   CONTINUE

        GOTO 300
  200   YOUT(I)=YIN(J)+(YIN(J)-YIN(J+1))*(XOUT(I)-XIN(J))/(XIN(J)- &
        XIN(J+1))
        I=I+1
        IF (I.GT.NOUT) GOTO 500
        JLAST=J
        GOTO 25
  300   CONTINUE
        DO 400 K=I,NOUT
        YOUT(K)=YIN(NIN)
  400   CONTINUE
  500   RETURN
        END

!
!   This routine is required if the CVS revision numbers are to 
!   survive an optimization.
!
!
!   1998/02/04 15:26:30 meyer
!
      subroutine getdiax_rev(i)
      CHARACTER*100 opt
      character*10 s 
      if( i .eq. 0) s =  &
      '@(#)getdiax.for,v 4.17\000'
      return
      end
