#include "config.f"
!**********************************************************************
!>
!!    This subroutine set the defines file name suffixes to look for
!!    (e.g. if running nw=129 ,nw=129, efit will look for Green's function
!!     tables with ec129129, re129129, etc.)
!!    
!!
!!    @param nw : number of grid point (width)
!!
!!    @param nh : number of grid point (height)
!!
!!    @param ch1 : character of number of grid point (width)
!!
!!    @param ch2 : character of number of grid point (height)
!!
!**********************************************************************
      subroutine table_name_ch(nw,nh,ch1,ch2)

      implicit none
      integer*4,   intent(in)  :: nw, nh 
      character*4, intent(out) :: ch1, ch2
      
      write(ch1,'(i4)') nw
      ch1 = adjustl(ch1)
      write(ch2,'(i4)') nh
      ch2 = adjustl(ch2)
      
      return
      end subroutine table_name_ch

!**********************************************************************
!*
!> set the directory for Green's function tables and machine params
!*
!**********************************************************************
   subroutine set_table_dir()
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'

      implicit none
#if defined(USEMPI)
      include 'mpif.h'
#endif
      integer*4 :: i,reason,nFiles,itmp,imin
      character(100) :: ftmp
      integer*4 :: shot_tables(100)
      character(100), dimension(:), allocatable :: filenames

      ! use directories specified in efit.setup if available
      if (link_efit(1:1).ne.'') then
        table_dir=trim(link_efit)//'green/'
        input_dir=trim(link_efit)
      endif
      if (link_store(1:1).ne.'')  store_dir=trim(link_store)

      ! recalculate length of default directories in case any change
      ltbdi2=0
      lindir=0
      lstdir=0
      table_di2 = table_dir
      do i=1,len(table_di2)
        if(table_di2(i:i).ne.' ') ltbdi2=ltbdi2+1
        if(input_dir(i:i).ne.' ') lindir=lindir+1
        if(store_dir(i:i).ne.' ') lstdir=lstdir+1
      enddo

      if (rank==0) then
        ! get the files
        call system('ls '//trim(table_di2)//' > shot_tables.txt')
        open(31,FILE='shot_tables.txt',action="read")

        !how many
        i = 0
        do
          read(31,FMT='(I10)',iostat=reason) shot_tables(i+1)
          if(reason/=0) EXIT
          i = i+1
        enddo
        nfiles = i 
        allocate(filenames(nfiles))
        rewind(31)
        do i = 1,nfiles
          read(31,'(a)') filenames(i)
        enddo
        close(31,status='delete')

        ! sort shot_tables just in case
        do i=1, nfiles
          imin = minloc(shot_tables(i:nfiles),dim=1)+i-1
          itmp = shot_tables(imin)
          ftmp = filenames(imin)
          shot_tables(imin) = shot_tables(i)
          filenames(imin) = filenames(i)
          shot_tables(i) = itmp
          filenames(i) = ftmp
        enddo 

        ! set the table directory
        do i=1, nfiles
          if(ishot.ge.shot_tables(i)) &
            table_di2 = table_di2(1:ltbdi2)//trim(filenames(i))//'/'
        enddo
        ltbdi2 = len(trim(table_di2))
        write(nttyo,10000) table_di2(1:ltbdi2)
      endif

#if defined(USEMPI)
      if (nproc > 1) then
        call MPI_BCAST(ltbdi2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(table_di2,ltbdi2,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      endif
#endif

10000 format(/,' table_di2 = <',a,'>')

   end subroutine set_table_dir

!*********************************************************************
!! 
!> reads Green's function tables and mhdin.dat data
!!
!*********************************************************************
   subroutine read_tables()
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'

      implicit none
      integer*4  :: i,j,istat,mcontr,mw,mh
      integer*4 vsid(nvesel) ! unused
      real*8 ecturn(necoil) ! unused
      character(1000) :: line
      parameter(mcontr=35)
 
      namelist/in3/mpnam2,xmp2,ymp2,amp2,smp2,patmp2, &
                   rsi,zsi,wsi,hsi,as,as2,rsisvs,lpname, &
                   rvs,zvs,wvs,hvs,avs,avs2,vsid,vsname, &
                   racoil,zacoil,wacoil,hacoil, &
                   rf,zf,wf,hf,af,af2,fcid,fcturn,turnfc, &
                   re,ze,we,he,ecid,ecturn
        
      ! initialize variables
      gridec = 0.0
      rmp2ec = 0.0
      rmp2vs = 0.0
      rsilvs = 0.0
      rsilec = 0.0
      fcid = 1
      ! for DIIID
      if(nfcoil.eq.18 .and. nfsum.eq.18) &
        fcid = (/1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 , &
                 10, 11, 12, 13, 14, 15, 16, 17, 18/)
      fcturn = 1.0
      rsi(1) = -1.
      re(1) = -1.
      rf(1) = -1.
      rvs(1) = -1.
!---------------------------------------------------------------------
!--   Read Green's tables from table_di2                            --
!---------------------------------------------------------------------
      open(unit=mcontr,status='old',form='unformatted', &
           file=table_di2(1:ltbdi2)//'ec'//trim(ch1)//trim(ch2)//'.ddd')
      read (mcontr) mw,mh

      if (.not.allocated(rgrid)) then
        allocate(rgrid(mw),stat=iallocate_stat)
        if(iallocate_stat/=0) stop "*** Not enough space for rgrid ***"
      endif
      if (.not.allocated(zgrid)) then
        allocate(zgrid(mh),stat=iallocate_stat)
        if(iallocate_stat/=0) stop "*** Not enough space for zgrid ***"
      endif
      read (mcontr) rgrid,zgrid

      if (.not.allocated(gridfc)) then
        allocate(gridfc(mw*mh,nfsum),stat=iallocate_stat)
        if(iallocate_stat/=0) stop "*** Not enough space for gridfc ***"
      endif
      read (mcontr) gridfc

      if (.not.allocated(gridpc)) then
        allocate(gridpc(mw*mh,mw),stat=iallocate_stat)
        if(iallocate_stat/=0) stop "*** Not enough space for gridpc ***"
      endif
      read (mcontr) gridpc
      close(unit=mcontr)
!----------------------------------------------------------------------
!--   read in the f coil response functions                          --
!----------------------------------------------------------------------
      open(unit=mcontr,form='unformatted', &
           status='old',file=table_di2(1:ltbdi2)//'rfcoil.ddd')
      if (.not.allocated(rsilfc)) then
        allocate(rsilfc(nsilop,nfsum),stat=iallocate_stat)
        if(iallocate_stat/=0) stop "*** Not enough space for rsilfc ***"
      endif

      read (mcontr) rsilfc
      if (.not.allocated(rmp2fc)) then
        allocate(rmp2fc(magpri,nfsum),stat=iallocate_stat)
        if(iallocate_stat/=0) stop "*** Not enough space for rmp2fc ***"
      endif
      read (mcontr) rmp2fc
      close(unit=mcontr)

      open(unit=mcontr,status='old',form='unformatted', &
           file=table_di2(1:ltbdi2)//'ep'//trim(ch1)//trim(ch2)//'.ddd')
      if (.not.allocated(gsilpc)) then
        allocate(gsilpc(nsilop,mw*mh),stat=iallocate_stat)
        if(iallocate_stat/=0) stop "*** Not enough space for gsilpc ***"
      endif
      read (mcontr) gsilpc

      if (.not.allocated(gmp2pc)) then
        allocate(gmp2pc(magpri,mw*mh),stat=iallocate_stat)
        if(iallocate_stat/=0) stop "*** Not enough space for gmp2pc ***"
      endif
      read (mcontr) gmp2pc
      close(unit=mcontr)
!----------------------------------------------------------------------
!--   read in the E-coil response function                           --
!----------------------------------------------------------------------
      if (iecurr.gt.0) then
        open(unit=mcontr,status='old',form='unformatted', &
          file=table_di2(1:ltbdi2)//'re'//trim(ch1)//trim(ch2)//'.ddd')
        read (mcontr) rsilec
        read (mcontr) rmp2ec
        read (mcontr) gridec
        close(unit=mcontr)
      endif
!----------------------------------------------------------------------
!--   read in the vessel response function                           --
!----------------------------------------------------------------------
      if (ifitvs.eq.1 .or. ivesel.gt.0) then
        open(unit=mcontr,status='old',form='unformatted', &
           file=table_di2(1:ltbdi2)//'rv'//trim(ch1)//trim(ch2)//'.ddd')
        read (mcontr) rsilvs
        read (mcontr) rmp2vs
        read (mcontr) gridvs
        close(unit=mcontr)
        ! try to protect against FP errors in read (answers can be sensitve...)
        if (ifitvs.eq.1 .or. ivesel.eq.3) then
          do i=1,nsilop
            do j=1,nvesel
              if(abs(rsilvs(i,j)).lt.1.e-5_dp) rsilvs(i,j)=0.0
            enddo
          enddo
          do i=1,magpri
            do j=1,nvesel
              if(abs(rmp2vs(i,j)).lt.1.e-5_dp) rmp2vs(i,j)=0.0
            enddo
          enddo
          do i=1,nwnh
            do j=1,nvesel
              if(abs(gridvs(i,j)).lt.1.e-5_dp) gridvs(i,j)=0.0
            enddo
          enddo
        endif
      endif

      open(unit=mcontr,status='old',file=table_di2(1:ltbdi2)//'mhdin.dat')

      read (mcontr,in3, iostat=istat)
      if (istat>0) then
        backspace(mcontr)
        read(mcontr,fmt='(A)') line
        write(*,'(A)') &
          'Invalid line in namelist in3: '//trim(line)
        stop
      endif
 
      if(rf(1).lt.0)&
        read (mcontr,10200) (rf(i),zf(i),wf(i),hf(i),af(i),af2(i), &
                             i=1,nfcoil)
      if(rsi(1).lt.0.) &
        read (mcontr,10200) (rsi(i),zsi(i),wsi(i),hsi(i),as(i),as2(i), &
                             i=1,nsilop)
      if(re(1).lt.0.) &
        read (mcontr,10220) (re(i),ze(i),we(i),he(i),ecid(i), &
                             i=1,necoil)
     
      if(rvs(1).lt.0 .and. (ifitvs.eq.1.or.ivesel.eq.3.or.icutfp.eq.2)) then
        read (mcontr,10200) (rvs(i),zvs(i),wvs(i),hvs(i), &
                             avs(i),avs2(i),i=1,nvesel)
      endif
      close(unit=mcontr)

10200 format (6e12.6)
10220 format (4e10.4,i5)
   end subroutine read_tables
