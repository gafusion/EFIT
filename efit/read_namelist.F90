subroutine read_efitin
     use commonblocks
     use set_kinds
     use mpi_efit
     include 'eparmdud129.inc'
     include 'modules2.inc'
     include 'modules1.inc'
     implicit integer*4 (i-n), real*8 (a-h,o-z)
     character(80),dimension(1001) :: inpfile
     logical input_flag
     integer*4 mode, shot, steps
     character cmdfile*15, shotfile*15, snapext*82
     real*8 starttime, deltatime

     namelist/optin/mode,cmdfile,shotfile,shot,starttime,deltatime,steps,snapext,inpfile
! OPT_INPUT >>>
      use_opt_input = .false.
      if (rank == 0) then
        inquire(file='efit.input',exist=input_flag)
      endif
      if (nproc > 1) then
        call MPI_BCAST(input_flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      endif
      if (input_flag) then
       if (rank == 0) then
          write(*,*) ' Using efit.input file...'
       endif
       ! ALL processes open and read efit.input file
       open(unit=nin,status='old',file='efit.input')
       read(nin,optin)
       close(nin)
       use_opt_input = .true.
       mode_in = mode
       cmdfile_in = cmdfile
       shotfile_in = shotfile
       shot_in = shot
       starttime_in = starttime
       deltatime_in = deltatime
       steps_in = steps
       snapext_in = snapext
       inpfile_in = inpfile(1:SIZE(inpfile_in))
     endif
     return

end subroutine read_efitin


subroutine read_eparmdud
  use eparmdud129
  use var_nio
  use expath
  implicit none
  integer ::istatus

  NAMELIST/machinein/nsilds,nsilol,nfcoil,nrogow,nacoil,mfcoil,necoil,nvesel, &
      mpress,nesum,magpri67,magpri322,magprirdp,magudom,maglds,mse315,mse45, &
      mse15,mse1h,mse315_2,mse210,libim,nmsels,nnece,nnecein,neceo,nnnte, &
      ngam_vars,ngam_u,ngam_w,nlimit,nlimbd,nangle,ntangle,nfbcoil,mccoil, &
      micoil,ndata,nwwcur,nffcur,nppcur,nercur,ntime,ndim,kxiter,mqwant, &
      mbdry,mbdry1,nxtram,nxtlim,nco2v,nco2r,modef,modep,modew,kubics, &
      icycred_loopmax,nfourier

  open(unit=nin,status='old',file=table_di2(1:ltbdi2)//'dprobe.dat')
  read (nin,machinein,iostat=istatus)
  close(unit=nin)

  return
end subroutine read_eparmdud

