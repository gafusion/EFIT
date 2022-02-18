!**********************************************************************
!>
!!    This routine returns the average of the compensated diamagnetic
!!    fluxes and error in the diamagnetic fluxes.
!!
!!    @param nshot : shot number
!!
!!    @param tim :
!!
!!    @param npts :
!!
!!    @param tavg : this is the time that the data is averaged 
!              over.
!!
!!    @param ierr : an error array of dimension 3. Each element of the
!        array contains the PTDATA error for one of the three 
!              diamagnetic flux signals. If ierr(j) is different from 
!              zero, that flux is not included in the average value  
!              returne
!!
!!    @param phidia : this array contains the conpensated diamagnetic
!                flux (mv-sec)
!!
!!    @param sigphi : this array contains the error in the flux value
!                (mv-sec).
!!
!********************************************************************** 
        subroutine getdia(nshot,tim,npts,tavg,ierr,phidia,sigphi)

        use var_gggttt, only: ntims
        dimension diamag(ntims,3),diamagc(ntims,3), &
                  sig(ntims,3),tim(ntims)
        dimension phidia(ntims),sigphi(ntims)
        dimension ierr(3),iwght(3)


        ndia=3
        call dlcomp(tim,diamag,diamagc,sig,nshot,npts,idlc,ierr, &
                    tavg)
        ierr(1)=1
        ierr(2)=1 !do not use diamag2
        do j=1,3 
            iwght(j)=1
            if (ierr(j).ne.0) then
                iwght(j)=0
                ndia=ndia-1
            endif
        enddo

        if (ndia.eq.0) then
            do i=1,npts
                phidia(i)=0.0
                sigphi(i)=0.0
            enddo
        else

            do i=1,npts
                phidia(i)=0.0
                sigphi(i)=0.0
                do j=1,3
                    phidia(i)=phidia(i)+iwght(j)*diamagc(i,j)
                    sigphi(i)=sigphi(i)+iwght(j)*sig(i,j)
                enddo
                phidia(i)=phidia(i)/ndia
                sigphi(i)=sigphi(i)/ndia
            enddo
        endif
        return
        end


!**********************************************************************
!>
!!    This subroutine applies diamagnetic compensations
!!    
!!
!**********************************************************************
        SUBROUTINE DLCOMP(TIM,DIAMAG,DIAMAGC,SIG,NSHOT,NPTS,IDLC &
        ,IERR,TAVG)
        use exvars
        use var_exdata, only: ishot
        use var_gggttt, only: ntims
        PARAMETER (NPOINT=42)
        REAL*4 RAR
        INTEGER*2 IIPOINT(5)
        INTEGER*4 IDAT,IAR,ASCII,INT16,INT32,REAL32
        DIMENSION DIAMAG(NTIMS,3),DIAMAGC(NTIMS,3),SIG(NTIMS,3)
        DIMENSION TIM(NTIMS),TTEMP(NTIMS)
        DIMENSION IDAT(NTIMS),YDAT(NTIMS),YTEMP(NTIMS),BT(NTIMS)
        DIMENSION IAR(128),RAR(8400),IERR(3),IADJ(3),BTCOMP(3),EBCOUP(3)
        DIMENSION FIX(3)
        DIMENSION ASCII(11),INT16(11),INT32(11),REAL32(11)

        DIMENSION COUP(NPOINT,3)
        INTEGER IDIASHOT
        CHARACTER*4 SOURCE
        CHARACTER*100 FILIN
        CHARACTER*10 POINT(NPOINT),IPOINT,DNAME(3),DNAMES(3)

        EQUIVALENCE (IIPOINT,IPOINT)
        DATA SOURCE/'.PLA'/
        DATA DNAME/'DIAMAG1   ','DIAMAG2   ','DIAMAG3   '/
        DATA DNAMES/'DIAMAG1S  ','DIAMAG2S  ','DIAMAG3S  '/
        DATA IADJ/-1,1,1/
        DATA BTCOMP/-0.3,0.4,-0.4/
        DATA EBCOUP/0.00,0.00,0.0/

        EQUIVALENCE (RAR(21),TTEMP(1))

        filin = 'dcoef.dat'
        if (nshot.lt.85885) then
            EBCOUP(3)=-2.0
        else
            EBCOUP(3)=0.0
        endif
!
        filin=input_dir(1:lindir)//filin

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
!       GET COMPENSATION COEFFICIENTS
!
        OPEN(UNIT=NIN,ACCESS='SEQUENTIAL', &
             STATUS='OLD',FILE=FILIN,iostat=ioerr)
        if (ioerr.ne.0) then
          do i=1,3
            ierr(i)=1
          enddo
          write(6,2020)
          return
        endif

41000   READ (NIN,*,iostat=ioerr) idiashot
        if (.not.is_iostat_end(ioerr)) then
          if (ioerr.ne.0) then
            do j=1,3
              ierr(j)=1
            enddo
            write(6,2020)
            return
          endif
          DO I=1,NPOINT
            READ(NIN,1001,iostat=ioerr) POINT(I),(COUP(I,J),J=1,3)
            if (ioerr.ne.0) then
              do j=1,3
                ierr(j)=1
              enddo
              write(6,2020)
              return
            endif
          ENDDO
          IF (nshot.lt.idiashot) then
            GOTO 41000
          ENDIF
        endif
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
!   DO 100 J=1,3
!   IPOINT=DNAMES(J)
!   CALL PTDATA(ITYP,NSHOT,SOURCE,%REF(IPOINT),IDAT,IER,IAR,RAR,
!     *  ASCII,INT16,INT32,REAL32)
!   IF (IER.NE.0) THEN
!     WRITE (6,1012)
!     WRITE (6,1013) IPOINT,NSHOT,IER,RAR(7),RAR(3),IAR(1)
!     TYPE *,' NO BT COMPENSATION DONE FOR ',DNAME(J)
!   ENDIF
!   DO 40 N=1,NPTS
!     DIAMAG(N,J)=(RAR(6)-IDAT(N))*RAR(5)*RAR(4)
!40 CONTINUE
!   BTCOEF(J)=DIAMAG(N69,J)/BCOIL(N69)
! 100 CONTINUE

!
!       adjust for an error in the rc/g for all three loops. The correct value was
!       measured on 9/4/87 and is correct for all shots after 56600.
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
!       GET DIAMAGNETIC FLUXES
!
        DO J=1,3
          IERR(J)=0
          IPOINT=DNAME(J)
          CALL PTDATA(ITYP,NSHOT,SOURCE,IIPOINT,IDAT,IER,IAR,RAR, &
                      ASCII,INT16,INT32,REAL32)
          IF (IER.NE.0 .AND. IER.NE.2 .AND. IER.NE.4) THEN
            WRITE (6,1012)
            WRITE (6,1013) IPOINT,NSHOT,IER,RAR(1),RAR(7), &
                           TTEMP(NPTS)
            IERR(J)=IER
            CYCLE
          ENDIF
          IF (IER.EQ.4) THEN
            IF (IAR(2).LT.NPTS) NPTS=IAR(2)
            IF (NPTS.LE.0) CYCLE
          ENDIF

          DO N=1,NPTS
            DIAMAG(N,J)=(RAR(6)-IDAT(N))*RAR(5)*RAR(4)*IADJ(J)*fix(j)
!        
!           COMPENSATE FOR BT PICKUP (0.1 MV-SEC IS UNCERTAINTY IN BT PICKUP)
!           NOV 21, 1994: NEW CALIB INCREASES BY 5%.
!
            DIAMAGC(N,J)=1.05*DIAMAG(N,J)-BTCOMP(J)
            SIG(N,J)=(DIAMAG(N,J)*SIGMAD)**2+0.1
          ENDDO
          if (navg.gt.1) then
            CALL LOWPASS (DIAMAG(1,J),DIAMAG(1,J),NAVG,NPTS)
            CALL LOWPASS (DIAMAGC(1,J),DIAMAGC(1,J),NAVG,NPTS)
          endif
        ENDDO

!       CREATE TIME ARRAY BASED ON THE DIAMAGNETIC LOOP TIME ARRAY

        IF (NPTS.GT.0) THEN
          DO N=1,NPTS
            TIM(N)=TTEMP(N)
          ENDDO
        ENDIF
!
!       COMPENSATE FOR E AND F-COIL PICKUP
!
        DO I=1,26
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
            CALL PTDATA(ITYP,NSHOT,SOURCE,IIPOINT,IDAT,IER,IAR,RAR, &
                        ASCII,INT16,INT32,REAL32)
          ELSEIF ((ISHOT .LT. 83350) .AND. (I .LT. 24)) THEN 
            CALL PTDATA(ITYP,NSHOT,SOURCE,IIPOINT,IDAT,IER,IAR,RAR, &
                        ASCII,INT16,INT32,REAL32)
          ENDIF
          if (ier.ne.0 .and. ier.ne.2 .and. ier.ne.4) then
            write (6,1012)
            write (6,1013) POINT(I),NSHOT,IER,RAR(1),TTEMP(1)
            idlc=idlc+1
            CYCLE
          endif       

   
          IF (IPOINT.EQ.'VLOOPBS   ' .OR. IREPLACE.EQ.1) THEN
            NPTS2=IAR(2)
            TAU=0.125
            SUM=0.0
            DO N=1,NPTS2-1
              IF (TTEMP(N).GT.TIM(NPTS)) EXIT
              YTEMP(N)=(RAR(6)-IDAT(N))*RAR(5)*RAR(4)
              SUM=SUM+EXP(TTEMP(N)/TAU)*YTEMP(N)*(TTEMP(N+1)-TTEMP(N))
              YTEMP(N)=YTEMP(N)-EXP(-1.*TTEMP(N)/TAU)/TAU*SUM
            ENDDO

            call interp(ttemp,ytemp,N,tim,ydat,npts)

          ELSE IF (IPOINT.EQ.'ECOILB    ' .OR. (ABS(RAR(1)-TTEMP(1)).GT. &
              RAR(3) .OR. ABS(TTEMP(NPTS)-TIM(NPTS)).GT.1.0E-6)) then
            npts2=iar(2)
            do n=1,npts2
              ytemp(n)=(rar(6)-idat(n))*rar(5)*rar(4)
            enddo
            call interp(ttemp,ytemp,iar(2),tim,ydat,npts)
          ELSE
            do n=1,npts
              ydat(n)=(rar(6)-idat(n))*rar(5)*rar(4)
            enddo
          ENDIF 

          IF (ipoint.eq.'ECOILB    ') then
            do n=1,npts
              ! TODO: bt is undefined here...
!              prod=((ydat(n)/6.6e4)**2)*bt(n)/2.13
              prod=((ydat(n)/6.6e4)**2)/2.13
              do j=1,3
                DIAMAGC(N,J)=DIAMAGC(N,J)-EBCOUP(J)*PROD
                sig(n,j)=sig(n,j)+(ebcoup(j)*prod)**2
              enddo
            enddo
          ENDIF


          DO J=1,3
            const=coup(i,j)*iadj(j)
            DO N=1,NPTS
              DIAMAGC(N,J)=DIAMAGC(N,J)-YDAT(N)*const
              SIG(N,J)=SIG(N,J)+(YDAT(N)*CONST*SIGMAF)**2
            ENDDO
          ENDDO

        ENDDO



        DO J=1,3
          IF (IERR(J).NE.0) THEN
            DO N=1,NPTS
              DIAMAG(N,J)=0.0
              DIAMAGC(N,J)=0.0
            ENDDO
          ENDIF
        ENDDO

        DO J=1,3
          DO N=1,NPTS
            SIG(N,J)=SQRT(SIG(N,J))
          ENDDO
        ENDDO

 1001   FORMAT(A10,3E12.4)
! 1002   FORMAT(' ',I3,2X,A10,3E12.4)
 1012   FORMAT(' *********************  ERROR IN PTDATA CALL ********')
 1013   FORMAT(' ',A10,2I6,4F9.4)

        RETURN

 2020   FORMAT(' error in open or read for file DCOEF.DAT')
        return
        END


!**********************************************************************
!>
!!    This subroutine does
!!    
!!
!**********************************************************************
        SUBROUTINE LOWPASS(XIN,XOUT,NAVG,NPTS)
        USE VAR_GGGTTT, ONLY: NTIMS
        REAL XIN(1),XOUT(1),XTMP(NTIMS)

        DO J=1,NPTS
          IND=(NAVG-1)/2
          IF (J.LE.IND) IND=J-1
          IF ((J+IND).GT.NPTS) IND=NPTS-J
          ICNT=0
          XTMP(J)=0

          DO I=J-IND,J+IND
            XTMP(J)=XTMP(J)+XIN(I)
            ICNT=ICNT+1
          ENDDO

          XTMP(J)=XTMP(J)/ICNT

        ENDDO

        DO J=1,NPTS
          XOUT(J)=XTMP(J)
        ENDDO

        RETURN
        END

!**********************************************************************
!>
!!    This subroutine does
!!    
!!
!**********************************************************************
        SUBROUTINE INTERP(XIN,YIN,NIN,XOUT,YOUT,NOUT)
        USE VAR_GGGTTT, ONLY: NTIMS
        REAL XIN(NTIMS),YIN(NTIMS),XOUT(NTIMS),YOUT(NTIMS)
        INTEGER NIN,NOUT

        IF (NIN.LT.1) RETURN

        DO I=1,NOUT
          IF (XOUT(I).GE.XIN(1)) go to 20
          YOUT(I)=YIN(1)
        ENDDO
        RETURN

   20   JLAST=1
   25   DO J=JLAST,NIN-1
          IF ((XOUT(I).GE.XIN(J)).AND.(XOUT(I).LT.XIN(J+1))) go to 200
        ENDDO

        go to 300
  200   YOUT(I)=YIN(J)+(YIN(J)-YIN(J+1))*(XOUT(I)-XIN(J))/(XIN(J)- &
        XIN(J+1))
        I=I+1
        IF (I.GT.NOUT) RETURN
        JLAST=J
        go to 25
  300   CONTINUE
        DO K=I,NOUT
          YOUT(K)=YIN(NIN)
        ENDDO
        RETURN
        END
