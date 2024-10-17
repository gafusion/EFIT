!
   program driver

   integer :: ishot           ! Desired shot number
   character(len=2) :: cdata  ! Desired data, i.e.,
                                          ! 'NE', 'TE' or 'PE'
                                               ! which data are desired
   integer,parameter :: ntime=5           ! Number of times
   real*4, dimension(ntime) :: time =(/1001.0, 1010.14, 2033.0, 2040.0, 3000.0/)
   real*4, dimension(ntime) :: sym  ! Symmetry points (m)
                                               ! evaluated on time array
   real*4, dimension(ntime) :: wid  ! Full widths (m)
                                               ! evaluated on time array
   real*4, dimension(ntime) :: ped  ! Pedestal values (units
                                            ! vary on time array
   logical, dimension(ntime) :: err            ! false means no error

   print *,' starting out '

   cdata = 'te'
   ishot = 105548
   call  get_mtanh_ts (ishot, cdata, ntime, time, sym, wid, ped, err) 
   print *,' sym = ', sym
   print *,' wid = ', wid 
   print *,' ped = ', ped 
   print *,' err = ', err
   print *,' ending '
       
   end program driver
!
