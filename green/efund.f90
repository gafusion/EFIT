!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          efun generates the necessary response                   **
!**          functions used by efit for reconstruction of the        **
!**          magnetic surfaces and current density profile.          **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1) d.w. swain and g.h. neilson, nucl. fusion           **
!**              22 (1982) 1015.                                     **
!**                                                                  **
!**********************************************************************
      program efund
      use machine, only: nw,nh,nwnh
      use nio
      use grid, only: efund_grid
      implicit none
      character inp1*4,inp2*4
      integer*4 ioerr

!----------------------------------------------------------------------
!--   Read in grid size from command line and set global variables   --
!----------------------------------------------------------------------
      nw = 0
      nh = 0
      call getarg(1,inp1)
      call getarg(2,inp2)
      read (inp1,'(i4)',iostat=ioerr) nw
      if(ioerr.eq.0) read (inp2,'(i4)') nh
      ! Ensure grid size is defined
      if (nw == 0) then
        write(*,*) 'Must specify grid dimensions as arguments'
        stop
      endif
      if(nh == 0) nh = nw
      nwnh = nw * nh
      ! make the file names for green-table
      call inp_file_ch(nw,nh,ch1,ch2)

      write(*,*) 'Reading namelist'
      call efund_getsizes
      call efund_getset
      write(*,*) 'Calling matrix subroutine'
      call efund_matrix
      write(*,*) 'Calling grid subroutine'
      call efund_grid

      stop 'GREEN TABLE GENERATED!'
      end program efund
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          matrix calculates the appropriate response functions.   **
!**                                                                  **
!**********************************************************************
      subroutine efund_matrix
      use utils, only: pi,flux
      use nio
      use grid
      use vessel
      use acoil
      use siloop
      use mprobe
      use fcoil
      use ecoil
      use rogowl
      implicit none
      integer*4 i,ii,j,jj
      real*8 taz,taz2 
      real*8 rfcfc(nfcoil,nfcoil)
      real*8 rsilfc(nsilop,nfcoil),rmp2fc(magpri,nfcoil), &
             rgowfc(nrogow,nfcoil)
      real*8 gsilfc(nsilop,nfsum),gmp2fc(magpri,nfsum)
      real*8 rsilec(nsilop,nesum),rmp2ec(magpri,nesum), &
             rfcec(nfcoil,nesum), &
             recec(nesum,nesum),rsisec(nesum)
      real*8 rsilvs(nsilop,nvesel),rmp2vs(magpri,nvesel), &
             rfcvs(nfcoil,nvesel), &
             rvsec(nvesel,necoil),rvsfc(nvesel,nfcoil), &
             rvsvs(nvesel,nvesel),tav(nvesel),tav2(nvesel)
      real*8 gsilvs(nsilop,nvsum),gmp2vs(magpri,nvsum), &
             gfcvs(nfsum,nvsum),gecvs(nesum,nvsum),gvsvs(nvsum,nvsum)
      real*8 taf(nfcoil),taf2(nfcoil)
      real*8 rsilac(nsilop,nacoil),rmp2ac(magpri,nacoil)
      real*8 xdum(1),ydum(1)
      real*8, dimension(:,:), allocatable :: rfcpc,brgrfc,bzgrfc, &
                                             rsilpc,rmp2pc,rgowpc, &
                                             gridec,gridvs,ggridvs, &
                                             gridac
!
      if (.not. allocated(rfcpc)) then
        allocate(rfcpc(nfcoil,nwnh))
        rfcpc = 0.0
      endif
      if (.not. allocated(brgrfc)) then
        allocate(brgrfc(nwnh,nfsum))
        brgrfc = 0.0
      endif
      if (.not. allocated(bzgrfc)) then
        allocate(bzgrfc(nwnh,nfsum))
        bzgrfc = 0.0
      endif
      if (.not. allocated(rsilpc)) then
        allocate(rsilpc(nsilop,nwnh))
        rsilpc = 0.0
      endif
      if (.not. allocated(rmp2pc)) then
        allocate(rmp2pc(magpri,nwnh))
        rmp2pc = 0.0
      endif
      if (.not. allocated(rgowpc)) then
        allocate(rgowpc(nrogow,nwnh))
        rgowpc = 0.0
      endif
      if (.not. allocated(gridec)) then
        allocate(gridec(nwnh,nesum))
        gridec = 0.0
      endif
      if (.not. allocated(gridvs)) then
        allocate(gridvs(nwnh,nvesel))
        gridvs = 0.0
      endif
      if (.not. allocated(ggridvs)) then
        allocate(ggridvs(nwnh,nvsum))
        ggridvs = 0.0
      endif
      if (.not. allocated(gridac)) then
        allocate(gridac(nwnh,nacoil))
        gridac = 0.0
      endif
    
      if (ifcoil.eq.1) then
!----------------------------------------------------------------------
!--      calculate the response function of psi loops due to f coils --
!----------------------------------------------------------------------
         do i=1,nfcoil
            taf(i)=tan(af(i)*pi/180.)
            taf2(i)=tan(af2(i)*pi/180.)
         enddo 
         if (islpfc.eq.1) then
            do i=1,nfcoil
               do j=1,nfcoil
                  call flux(rf(i),zf(i),wf(i),hf(i),taf(i),taf2(i), &
                            rf(j),zf(j),wf(j),hf(j),taf(j),taf2(j), &
                            rfcfc(j,i))
                  rfcfc(j,i)=rfcfc(j,i)*0.5/pi

               enddo
               ii=i
               call gsilop(rgrid,nw,zgrid,nh,rfcpc,ii,rf,zf,wf,hf,af,af2, &
                           nfcoil)
            enddo

            print*,'file name : ','fc'//trim(ch1)//trim(ch2)//'.ddd' 
!vasorg      open(unit=nrspfc,status='unknown',file='fcfcpc.dat', &
            open(unit=nrspfc,status='unknown',file='fc'//trim(ch1)// & 
                 trim(ch2)//'.ddd' , &
                 form='unformatted')
            write (nrspfc) rfcfc
            write (nrspfc) rfcpc
            close(unit=nrspfc)
         endif
!---------------------------------------------------------------------
!--      flux loops                                                 --
!---------------------------------------------------------------------
         if (nsilop.gt.1) then
!---------------------------------------------------------------------
!           finite size flux loops                                  --
!---------------------------------------------------------------------
            if (isize.gt.0) then
               do i=1,nfcoil
                  do j=1,isize
                     taz=tan(as(j)*pi/180.)
                     taz2=tan(as2(j)*pi/180.)
                     call flux(rf(i),zf(i),wf(i),hf(i),taf(i),taf2(i), &
                               rsi(j),zsi(j),wsi(j),hsi(j),taz,taz2, &
                               rsilfc(j,i))
                     rsilfc(j,i)=rsilfc(j,i)*0.5/pi
                  enddo 
               enddo
            endif
!---------------------------------------------------------------------
!           thin flux loops                                         --
!---------------------------------------------------------------------
            if (isize.lt.nsilop) then
               do i=1,nfcoil
                  ii=i
                  do j=isize+1,nsilop
                     jj=j
                     call m1coef(xdum,xdum,nsilop,nfcoil,rsilfc,jj,ii)
                  enddo 
               enddo 
            endif 
         endif 

         !
         if (.not. allocated(brgridfc)) then
            allocate(brgridfc(nwnh,nfcoil))
            brgridfc(:,:) = 0.0
         endif
         if (.not. allocated(bzgridfc)) then
            allocate(bzgridfc(nwnh,nfcoil))
            bzgridfc(:,:) = 0.0
         endif
!
!-----------------------------------------------------------------------
!--      compute the response function of magnetic probes due to f coils
!-----------------------------------------------------------------------
         if(magpri.gt.1) &
            call m2coef(xdum,0,ydum,0,rmp2fc,magpri,nfcoil)
!----------------------------------------------------------------------
!--      compute the response function of partial rogowski loops due to
!--      f coils
!----------------------------------------------------------------------
         if(nrogow.gt.1) &
            call rogowc(xdum,0,ydum,0,rgowfc,nrogow,nfcoil)
!----------------------------------------------------------------------
!--      write f coil response functions                             --
!----------------------------------------------------------------------
         gsilfc=0.0
         gmp2fc=0.0
         do i=1,nfcoil
            do j=1,nsilop
               gsilfc(j,fcid(i))=gsilfc(j,fcid(i))+fcturn(i)*rsilfc(j,i)
            enddo
            do j=1,magpri
               gmp2fc(j,fcid(i))=gmp2fc(j,fcid(i))+fcturn(i)*rmp2fc(j,i)
            enddo
         enddo
!
         print*,'file name : ','rfcoil.ddd' 
!vasorg      open(unit=nrspfc,status='unknown',file='rfcoil.dat', &
         open(unit=nrspfc,status='unknown',file='rfcoil.ddd', &
              form='unformatted')
         write (nrspfc) gsilfc
         write (nrspfc) gmp2fc
         close(unit=nrspfc)
!
         brgrfc=0.0
         bzgrfc=0.0
         do i=1,nfcoil
            do j=1,nwnh
               brgrfc(j,fcid(i))=brgrfc(j,fcid(i))+fcturn(i)*brgridfc(j,i)
               bzgrfc(j,fcid(i))=bzgrfc(j,fcid(i))+fcturn(i)*bzgridfc(j,i)
            enddo
         enddo
!
         open(unit=nrspfc,status='unknown',file='brzgfc.dat', &
              form='unformatted')
         write (nrspfc) brgrfc
         write (nrspfc) bzgrfc
         close(unit=nrspfc)
      endif
!----------------------------------------------------------------------
!--   plasma response functions                                      --
!----------------------------------------------------------------------
      if (igrid.eq.1) then
         if (nsilop.gt.1) then
!----------------------------------------------------------------------
!--         filament plasma current model                            --
!----------------------------------------------------------------------
            if (isize.gt.0) then
               do j=1,isize
                  jj=j
                  call gsilop(rgrid,nw,zgrid,nh,rsilpc,jj, &
                              rsi,zsi,wsi,hsi,as,as2,nsilop)
               enddo
            endif
         endif
         if (isize.lt.nsilop) then
            do j=isize+1,nsilop
               jj=j
               call m1coef(rgrid,zgrid,nw,nh,rsilpc,jj,0)
            enddo 
         endif
         if(magpri.gt.1) &
            call m2coef(rgrid,nw,zgrid,nh,rmp2pc,magpri,nwnh)
         if(nrogow.gt.1) &
            call rogowc(rgrid,nw,zgrid,nh,rgowpc,nrogow,nwnh)
!----------------------------------------------------------------------
!--      write the plasma response function                          --
!----------------------------------------------------------------------
         print*,'file name : ','ep'//trim(ch1)//trim(ch2)//'.ddd' 
!vasorg      open(unit=nrsppc,status='unknown',file='eplasm.dat', &
         open(unit=nrsppc,status='unknown',file='ep'//trim(ch1)// & 
              trim(ch2)//'.ddd' , &
              form='unformatted')
         write (nrsppc) rsilpc
         write (nrsppc) rmp2pc
         close(unit=nrsppc)
!
      endif
      if (iecoil.eq.1) then
         call gecoil(rsilec,rmp2ec,gridec,rgrid,nw, &
                     zgrid,nh,rfcec,recec,rsisec)
         print*,'file name : ','re'//trim(ch1)//trim(ch2)//'.ddd' 
!vasorg      open(unit=nrsppc,status='unknown',file='recoil.dat', &
         open(unit=nrsppc,status='unknown',file='re'//trim(ch1)// & 
                     trim(ch2)//'.ddd', &
              form='unformatted')
         write (nrsppc) rsilec
         write (nrsppc) rmp2ec
         write (nrsppc) gridec
         close(unit=nrsppc)
      endif
!
      if (ivesel.eq.1) then
         call gvesel(rsilvs,rmp2vs,gridvs,rgrid,nw, &
                     zgrid,nh,rfcvs,rvsfc,rvsec)
         do i=1,nvesel
            tav(i)=tan(avs(i)*pi/180.)
            tav2(i)=tan(avs2(i)*pi/180.)
         enddo 
         do i=1,nvesel
            do j=1,nvesel
               call flux(rvs(i),zvs(i),wvs(i),hvs(i),tav(i),tav2(i), &
                         rvs(j),zvs(j),wvs(j),hvs(j),tav(j),tav2(j), &
                         rvsvs(j,i))
               rvsvs(j,i)=rvsvs(j,i)*0.5/pi
            enddo 
         enddo 
!
         gsilvs=0.0
         gmp2vs=0.0
         ggridvs=0.0
         gfcvs=0.0
         gecvs=0.0
         gvsvs=0.0
         do i=1,nvesel
            do j=1,nsilop
               gsilvs(j,vsid(i))=gsilvs(j,vsid(i))+rsilvs(j,i)
            enddo
            do j=1,magpri
               gmp2vs(j,vsid(i))=gmp2vs(j,vsid(i))+rmp2vs(j,i)
            enddo
            do j=1,nwnh
               ggridvs(j,vsid(i))=ggridvs(j,vsid(i))+gridvs(j,i)
            enddo
            do j=1,nfcoil
               gfcvs(fcid(j),vsid(i))=gfcvs(fcid(j),vsid(i))+rvsfc(i,j)
            enddo
            do j=1,necoil
               gecvs(ecid(j),vsid(i))=gecvs(ecid(j),vsid(i))+rvsec(i,j)
            enddo
            do j=1,nvesel
               gvsvs(vsid(j),vsid(i))=gvsvs(vsid(j),vsid(i))+rvsvs(i,j)
            enddo
         enddo
!
!vas
         print*,'file name : ','rv'//trim(ch1)//trim(ch2)//'.ddd' 
!vasorg         open(unit=nrsppc,status='unknown',file='rvesel.dat', &
         open(unit=nrsppc,status='unknown',file='rv'//trim(ch1)// & 
              trim(ch2)//'.ddd' , &
              form='unformatted')
         write (nrsppc) gsilvs
         write (nrsppc) gmp2vs
         write (nrsppc) ggridvs
         write (nrsppc) gfcvs
         write (nrsppc) gecvs
         write (nrsppc) gvsvs
         close(unit=nrsppc)
      endif
!---------------------------------------------------------------------
!--   advance divertor coil                                         --
!---------------------------------------------------------------------
      if (iacoil.eq.1) then
         call gacoil(rsilac,rmp2ac,gridac,rgrid,nw, &
                     zgrid,nh)
!vas
         print*,'file name : ','ra'//trim(ch1)//trim(ch2)//'.ddd' 
!vasorg         open(unit=nrsppc,status='unknown',file='racoil.dat', &
         open(unit=nrsppc,status='unknown',file='ra'//trim(ch1)// & 
              trim(ch2)//'.ddd' , &
              form='unformatted')
         write (nrsppc) gridac
         write (nrsppc) rsilac
         write (nrsppc) rmp2ac
         close(unit=nrsppc)
      endif
!
      if (allocated(rfcpc)) deallocate(rfcpc)
      if (allocated(brgrfc)) deallocate(brgrfc)
      if (allocated(bzgrfc)) deallocate(bzgrfc)
      if (allocated(rsilpc)) deallocate(rsilpc)
      if (allocated(rmp2pc)) deallocate(rmp2pc)
      if (allocated(rgowpc)) deallocate(rgowpc)
      if (allocated(gridec)) deallocate(gridec)
      if (allocated(gridvs)) deallocate(gridvs)
      if (allocated(ggridvs)) deallocate(ggridvs)
      if (allocated(gridac)) deallocate(gridac)
!
      return
      end subroutine efund_matrix
