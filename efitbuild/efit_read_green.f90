   subroutine efit_read_green
      use set_kinds
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      namelist/in3/mpnam2,xmp2,ymp2,amp2,smp2,rsi,zsi,wsi,hsi,as, &
        as2,lpname,rsisvs,vsname,turnfc,patmp2,racoil,zacoil, &
        hacoil,wacoil
!---------------------------------------------------------------------
!-- Read Green's tables from table_dir            --
!---------------------------------------------------------------------
      open(unit=mcontr,status='old',form='unformatted', &
           file=table_dir(1:ltbdir)//'ec'//trim(ch1)//trim(ch2)//'.ddd')
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
           status='old',file=table_dir(1:ltbdir)//'rfcoil.ddd')
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
           file=table_dir(1:ltbdir)//'ep'//trim(ch1)//trim(ch2)//'.ddd')
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
           file=table_dir(1:ltbdir)//'re'//trim(ch1)//trim(ch2)//'.ddd')
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
           file=table_dir(1:ltbdir)//'rv'//trim(ch1)//trim(ch2)//'.ddd')
        read (mcontr) rsilvs
        read (mcontr) rmp2vs
        read (mcontr) gridvs
        close(unit=mcontr)
      endif
!
      open(unit=mcontr,status='old',file=table_di2(1:ltbdi2)//'dprobe.dat')
      rsi(1)=-1.

      read (mcontr,in3)
      read (mcontr,10200) (rf(i),zf(i),wf(i),hf(i),af(i),af2(i), &
                i=1,mfcoil)
      if (rsi(1).lt.0.) &
        read (mcontr,10200) (rsi(i),zsi(i),wsi(i),hsi(i),as(i),as2(i), &
                i=1,nsilop)
      read (mcontr,10220) (re(i),ze(i),we(i),he(i),ecid(i), &
                                        i=1,necoil)
      if (ifitvs.gt.0.or.icutfp.eq.2) then
        read (mcontr,10200) (rvs(i),zvs(i),wvs(i),hvs(i), &
                                        avs(i),avs2(i),i=1,nvesel)
      endif
      close(unit=mcontr)

    print *, 'end     '
10200 format (6e12.6)
10220 format (5e10.4)
   end subroutine efit_read_green
