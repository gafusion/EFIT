   subroutine efit_read_tables
      use set_kinds
      include 'eparmdud129.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      integer*4  :: istat
      character(len=1000) :: line
 
      namelist/in3/mpnam2,xmp2,ymp2,amp2,smp2,rsi,zsi,wsi,hsi,as, &
        as2,lpname,rsisvs,vsname,turnfc,patmp2,racoil,zacoil, &
        hacoil,wacoil,rf,zf,fcid,wf,hf,wvs,hvs,avs,avs2,af,af2,fcturn, &
        re,ze,ecid,ecturn,vsid,rvs,zvs,we,he 
        
      mcontr = 35
!---------------------------------------------------------------------
!-- Read Green's tables from table_di2            --
!---------------------------------------------------------------------
      open(unit=mcontr,status='old',form='unformatted', &
           file=table_di2(1:ltbdi2)//'ec'//trim(ch1)//trim(ch2)//'.ddd')
      read (mcontr) mw,mh

      if(.not.allocated(rgrid)) then
        allocate(rgrid(mw),stat=iallocate_stat)
        if(iallocate_stat/=0) stop "*** Not enough space for rgrid ***"
      endif
      if(.not.allocated(zgrid)) then
        allocate(zgrid(mh),stat=iallocate_stat)
        if(iallocate_stat/=0) stop "*** Not enough space for zgrid ***"
      endif
      read (mcontr) rgrid,zgrid

      if(.not.allocated(gridfc)) then
        allocate(gridfc(mw*mh,nfcoil),stat=iallocate_stat)
        if(iallocate_stat/=0) stop "*** Not enough space for gridfc ***"
      endif
      read (mcontr) gridfc

      if(.not.allocated(gridpc)) then
        allocate(gridpc(mw*mh,mw),stat=iallocate_stat)
        if(iallocate_stat/=0) stop "*** Not enough space for gridpc ***"
      endif
      read (mcontr) gridpc
      close(unit=mcontr)

!----------------------------------------------------------------------
!-- read in the f coil response functions                            --
!----------------------------------------------------------------------
      open(unit=mcontr,form='unformatted', &
           status='old',file=table_di2(1:ltbdi2)//'rfcoil.ddd')
      if(.not.allocated(rsilfc)) then
        allocate(rsilfc(nsilop,nfcoil),stat=iallocate_stat)
        if(iallocate_stat/=0) stop "*** Not enough space for rsilfc ***"
      endif

      read (mcontr) rsilfc
      if(.not.allocated(rmp2fc)) then
        allocate(rmp2fc(magpri,nfcoil),stat=iallocate_stat)
        if(iallocate_stat/=0) stop "*** Not enough space for rmp2fc ***"
      endif
      read (mcontr) rmp2fc
      close(unit=mcontr)

      open(unit=mcontr,status='old',form='unformatted', &
           file=table_di2(1:ltbdi2)//'ep'//trim(ch1)//trim(ch2)//'.ddd')
      if(.not.allocated(gsilpc)) then
        allocate(gsilpc(nsilop,mw*mh),stat=iallocate_stat)
        if(iallocate_stat/=0) stop "*** Not enough space for gsilpc ***"
      endif
      read (mcontr) gsilpc

      if(.not.allocated(gmp2pc)) then
        allocate(gmp2pc(magpri,mw*mh),stat=iallocate_stat)
        if(iallocate_stat/=0) stop "*** Not enough space for gmp2pc ***"
      endif
      read (mcontr) gmp2pc
      close(unit=mcontr)

!----------------------------------------------------------------------
!-- read in the E-coil response function                             --
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
!-- read in the vessel response function                             --
!----------------------------------------------------------------------
      if (ivesel.gt.0) then
        open(unit=mcontr,status='old',form='unformatted', &
           file=table_di2(1:ltbdi2)//'rv'//trim(ch1)//trim(ch2)//'.ddd')
        read (mcontr) rsilvs
        read (mcontr) rmp2vs
        read (mcontr) gridvs
        close(unit=mcontr)
      endif
!
      open(unit=mcontr,status='old',file=table_di2(1:ltbdi2)//'dprobe.dat')
      rsi(1)=-1.
      re(1) = -1.
      rf(1) = -1.
      rvs(1) = -1.

      read (mcontr,in3, iostat=istat)

      if (istat/=0) then
        backspace(mcontr)
        read(mcontr,fmt='(A)') line
        write(*,'(A)') &
         'Invalid line in namelist: '//trim(line)
      end if
 
      if (rf(1).lt.0)&
      read (mcontr,10200) (rf(i),zf(i),wf(i),hf(i),af(i),af2(i), &
                i=1,mfcoil)
      if (rsi(1).lt.0.) &
        read (mcontr,10200) (rsi(i),zsi(i),wsi(i),hsi(i),as(i),as2(i), &
                i=1,nsilop)
      if (re(1).lt.0.) &
        read (mcontr,10220) (re(i),ze(i),we(i),he(i),ecid(i), &
                                        i=1,necoil)
     
      if (rvs(1).lt.0 .and. (ifitvs.gt.0.or.icutfp.eq.2)) then
        read (mcontr,10200) (rvs(i),zvs(i),wvs(i),hvs(i), &
                                        avs(i),avs2(i),i=1,nvesel)
      endif
      close(unit=mcontr)


10200 format (6e12.6)
10220 format (5e10.4)
   end subroutine efit_read_tables


   subroutine set_table_dir
      use set_kinds
      include 'eparmdud129.inc'
      include 'modules2.inc'
      include 'modules1.inc'

      implicit integer*4 (i-n), real*8 (a-h,o-z)
      integer :: i, reason, nFiles, itmp, imin 
      character(LEN=100), dimension(:), allocatable :: filenames
      integer:: shot_tables(100)

!----------------------------------------------------------------------
!--   recalculate length of default directories in case any change   --
!----------------------------------------------------------------------
      ltbdi2=0
      lindir=0
      lstdir=0
      table_di2 = table_dir
      do i=1,len(table_di2)
         if (table_di2(i:i).ne.' ') ltbdi2=ltbdi2+1
         if (input_dir(i:i).ne.' ') lindir=lindir+1
         if (store_dir(i:i).ne.' ') lstdir=lstdir+1
      enddo

      ! get the files
      call system('ls '//trim(table_di2)//' > shot_tables.txt')
      open(31,FILE='shot_tables.txt',action="read")

      !how many
      i = 0
      do
        read(31,FMT='(I10)',iostat=reason) shot_tables(i+1)
        if (reason/=0) EXIT
        i = i+1
      end do
      nfiles = i 
      allocate(fileNames(nfiles))
      rewind(31)
      do i = 1,nfiles
        read(31,'(a)') filenames(i)
      end do
      close(31)

      ! sort shot_tables just in case
      do i=1, nfiles
        imin = minloc(shot_tables(i:nfiles),dim=1)+i-1
        itmp = shot_tables(imin)
        shot_tables(imin) = shot_tables(i)
        shot_tables(i)  = itmp
      enddo 

      do i=1, nfiles
        if (ishot.ge. shot_tables(i)) then
          table_di2 = table_di2(1:ltbdi2)//trim(filenames(i))//'/'
        endif
      enddo
      
      ltbdi2 = len(trim(table_di2))

      if (rank == 0) then
        write(*,*) 'table_di2 = <',table_di2(1:ltbdi2),'>'
      endif

   end subroutine set_table_dir
