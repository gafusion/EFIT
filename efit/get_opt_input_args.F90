!**********************************************************************
!>
!!    Subroutine for converting command line arguments into the logical
!!    parameters for control
!!
!**********************************************************************
SUBROUTINE process_arguments()
  CHARACTER(128) :: argument 
  INTEGER(KIND=4) :: iarg, nargs,iloop,kfindex,ipos1,ipos2,ilen,istart
  INTEGER(KIND=4) :: it,iscan,ifarg 
  nargs=0 
  CALL get_arg_count(nargs) 
  IF (nargs == 0) THEN
    CALL write_usage 
    STOP('please correct')
  ENDIF
                                                                 
  iarg=0 
  loop: DO iloop=1,nargs 
    iarg=iarg+1 
    CALL get_arg(iarg,argument) 
    SELECT CASE(argument) 
    CASE('-h','--help')
      CALL write_usage 
      STOP('normal termination')
    ! Example boolean
    CASE('-bool') 
      boolvar=.TRUE. 
      CYCLE 
    ! Example of additional argument
    CASE('-subarg') 
      iarg=iarg+1 
      CALL get_arg(iarg,argument) 
      CYCLE 
    CASE DEFAULT 
      ! Rest are something else
      nfiles=nargs-iarg+1 
      IF (nfiles==0) THEN
        CALL write_usage 
        CALL nim_stop('please specify dump file(s)')
      ENDIF
      DO ifarg=iarg,nargs 
        CALL get_arg(ifarg,dump_file) 
        INQUIRE(FILE=TRIM(dump_file),EXIST=file_stat) 
        IF (.NOT. file_stat) THEN 
          msg = 'File does not exist: '//TRIM(dump_file) 
          CALL nim_stop(msg) 
        ENDIF 
        ifile=nargs-ifarg+1 
        file_list(ifile)=TRIM(dump_file) 
      ENDDO 
      EXIT loop 
    END SELECT 
  ENDDO loop 
END SUBROUTINE process_arguments 



!**********************************************************************
!>     subprogram write_usage
!!
!**********************************************************************
SUBROUTINE write_usage()
  WRITE(nim_wr,*)'' 
  msg='Usage: efit <options> other' 
  WRITE(6,*) msg 
  WRITE(6,*)'' 
  WRITE(6,*)'where the  <options> are' 
  WRITE(6,*)'-h ' 
  WRITE(6,*)'  Write this message'
  WRITE(6,*)'-bool ' 
  WRITE(6,*)'  example boolean'
END SUBROUTINE write_usage 

