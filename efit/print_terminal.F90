#include "config.f"
!**********************************************************************
!>
!!    prints stats to the terminal
!!
!!    @param it : time index
!!
!**********************************************************************
      subroutine print_stats(it)
      use commonblocks, only: worka2
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
#if defined(USEMPI)
      include 'mpif.h'
#endif
      integer*4, intent(in) :: it
      integer*4 i,im1,ip1,k,m,mb,mbb,ioerr,nves2
      real*8 cbetap,cbetat,cipmp3,dampz1,dampz2, &
             dang270,dang90,danm270,danm90,danp270,danp90, &
             dsi,nzz,saisq,sinow,sm1,sm2,ssiref, &
             sumif,sumifb,sumifs,sumift,tin,tout, &
             xdum,xnorm,xtest,ytest,xxm,yym,xxp,yyp
      character*30 sfname
      character(1000) :: line
      real*8 xrsp(npcurn)
      real*8 patmpz(magpri),xmpz(magpri),ympz(magpri),ampz(magpri), &
             rsisfc(nfcoil)
      integer*4, dimension(:), allocatable :: ishotall
      real*8, dimension(:), allocatable :: ch2all,timeall
      real*8, parameter :: zetafc=2.5e-08_dp
      ! DIIID specific parameters for writing mhdin.new file
      integer*4, parameter :: magpri67=29,magpri322=31

      namelist/in3/mpnam2,xmp2,ymp2,amp2,smp2,patmp2, &
                   rsi,zsi,wsi,hsi,as,as2,rsisvs,lpname, &
                   rvs,zvs,wvs,hvs,avs,avs2,vsname, &
                   racoil,zacoil,wacoil,hacoil, &
                   rf,zf,wf,hf,af,af2,fcid,fcturn,turnfc, &
                   re,ze,we,he,ecid!,ecturn,vsid ! not used or saved in efit

      ! avoid write garbage when gapin or gapout not update (e.g. 1.e10)
      tin=gapin(it)
      if(tin.gt.1.e9_dp) tin=0.0
      tout=gapout(it)
      if(tout.gt.1.e9_dp) tout=0.0
!
!#if defined(USEMPI)
!      ! TODO: Currently this ONLY works if nproc == total number of time slices.
!      ! If running in parallel, write out key info for all ranks
!      if (nproc > 1) then
!        if (rank==0) then
!          allocate(ishotall(nproc))
!          allocate(timeall(nproc))
!          allocate(ch2all(nproc))
!        end if
!        call mpi_gather(ishot, 1, MPI_INTEGER, ishotall, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!        call mpi_gather(time(it), 1, MPI_DOUBLE_PRECISION, timeall, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
!        call mpi_gather(chisq(it), 1, MPI_DOUBLE_PRECISION, ch2all, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
!        if (rank==0) then
!          write(nttyo,'(/,a)') '  Summary of all runs'
!          write(nttyo,'(a)') '   Shot    Time   chi^2'
!          do i = 1,nproc
!            write(nttyo,'(i8,1x,i6,1x,f10.7)') ishotall(i),int(timeall(i)),ch2all(i)
!            !write(nttyo,*) ishotall(i),int(timeall(i)),ch2all(i)
!          end do
!          if (allocated(ishotall)) deallocate(ishotall)
!          if (allocated(ch2all)) deallocate(ch2all)
!          if (allocated(timeall)) deallocate(timeall)
!        end if
!      else
!        !write(nttyo,*) ishot,int(time(it)),chisq(it) ! handy for debugging
!      end if
!#endif
!      if ((ivacum.eq.0).and.(itek.le.0)) then
      if (ivacum.eq.0) then
        mpi_rank: if (rank == 0) then
          !write (nttyo,10000) trim(ch1),trim(ch2)
          saisq=chisq(it)
          ! avoid roundoff
          if(abs(saisq).lt.1.e-20_dp) saisq=0_dp
          write (nttyo,9300)
          write (nttyo,9320) ipsi(it)
          write (nttyo,9340) imag2(it)
          write (nttyo,9380) iplasm(it)
          write (nttyo,9385) idlopc(it)
          write (nttyo,10480)
          write (nttyo,10500) ishot,int(time(it)),saisq
          write (nttyo,10520) betat(it),betap(it),li(it)
          write (nttyo,10540) volume(it),rout(it),zout(it)
          write (nttyo,10560) elong(it),utri(it),ltri(it)
          write (nttyo,10580) aminor(it),tin,tout
          write (nttyo,10600) gaptop(it),qstar(it),rcurrt(it)
          write (nttyo,10610) zcurrt(it),bcentr(it),qout(it)
        endif mpi_rank
      endif
!
! --- delete fitout.dat if IOUT does not contain 1
!
      if (iand(iout,1).eq.0) then ! if iout is even, including zero
        mpi_rank0: if (rank == 0) then
          close(unit=nout,status='delete')
        endif mpi_rank0
      else
        !write (nout,10000) trim(ch1),trim(ch2)
        saisq=chisq(it)
        write (nout,9300)
        write (nout,9320) ipsi(it)
        write (nout,9340) imag2(it)
        write (nout,9380) iplasm(it)
        write (nout,9385) idlopc(it)
        write (nout,10480)
        write (nout,10500) ishot,int(time(it)),saisq
        write (nout,10520) betat(it),betap(it),li(it)
        write (nout,10540) volume(it),rout(it),zout(it)
        write (nout,10560) elong(it),utri(it),ltri(it)
        write (nout,10580) aminor(it),tin,tout
        write (nout,10600) gaptop(it),qstar(it),rcurrt(it)
        write (nout,10610) zcurrt(it),bcentr(it),qout(it)
        write (nout,10620) sepin(it),sepout(it),septop(it)
        write (nout,10623) betat2

        if (icalbet.gt.0) &
          write (nout,10622) betat(it),vbtvac,vbtot2,vbtvac2,vbtor2,vbeta0
        if (icalbet.gt.0) &
          write (nout,10624) vbtmag,btvvac2,btvtor2,btvtot2

        if (scalea) then
          rowcnd=1./rowcnd
          colcnd=1./colcnd
          write (nout,11001)
          write (nout,11004) infosc,rowcnd,colcnd,arspmax
        endif

        write (nout,11000)
        write (nout,11020) (brsp(i)/turnfc(i),i=1,nfsum)
        write (nout,11030)
        write (nout,11020) (brsp(i),i=1,nfsum)

        sumif=0.0
        do i=1,5
          sumif=sumif+brsp(i)+brsp(i+9)
        end do
        sumif=sumif+brsp(8)+brsp(17)
        sumift=0.0
        sumifs=0.0
        do i=1,nfsum
          sumift=sumift+brsp(i)
          sumifs=sumifs+brsp(i)**2
        end do
        sumifs=sqrt(sumifs/(nfsum-1))
        write (nout,10020) sumif,sumift,sumifs

        rsisfc=-1.
        do i=1,nfcoil 
          if(rsisfc(i).le.-1.0) & 
            rsisfc(i)=(turnfc(fcid(i))*fcturn(i))**2 &
                      *twopi*rf(i)/wf(i)/hf(i)*zetafc 
        enddo
        write (nout,11010)
        write (nout,11020) (rsisfc(i),i=1,nfsum)

        if (ivacum.eq.0.and.(icurrt.eq.2.or.icurrt.eq.5)) then
          xnorm=brsp(nfsum+1)
          write (nout,11040) xnorm
          xrsp=0.0
          if (xnorm.ne.0.0) then
            do i=1,kwcurn
              xrsp(i)=brsp(nfsum+i)/xnorm
            end do
          endif
          write (nout,11020) (xrsp(i),i=1,kwcurn)
          xnorm=darea
          write (nout,11043) xnorm
          do i=1,kwcurn
            xrsp(i)=brsp(nfsum+i)/xnorm
          end do
          write (nout,11020) (xrsp(i),i=1,kwcurn)
          if (keecur.gt.0) then
            write (nout,11032)
            write (nout,11020) (cerer(i),i=1,keecur)
          endif
          if (kedgep.gt.0) then
            write (nout,11034) pedge
          endif
          if (kedgef.gt.0) then
            write (nout,11036) f2edge
          endif
        endif

        write (nout,11005)
        write (nout,11020) (cecurr(i),i=1,nesum)
        write (nout,11008)
        write (nout,11020) (rsisec(i),i=1,nesum)

        if (ivesel.gt.0) then
          sumif=sum(vcurrt)
          nves2=nvesel/2
          sumift=sum(vcurrt(1:nves2))
          sumifb=sum(vcurrt(nves2+1:nvesel))
          write (nout,11060) sumif,sumift,sumifb
          write (nout,11020) (vcurrt(i),i=1,nvesel)
          write (nout,11022) fzpol
          write (nout,11020) (vforcep(i),i=1,nvesel)
          write (nout,11024) fztor
          write (nout,11020) (vforcet(i),i=1,nvesel)
        endif
      endif

!-----------------------------------------------------------------
!--   Create a new machine file
!--   Warning: this is only intended for DIIID
!-----------------------------------------------------------------
      newmhdin: if (device == "DIII-D") then
      if (patmp2(1).le.0.0 .and. xmp2(1).le.0.0) &
        write(*,*) "bad code detected"
      if (xmp2(1).gt.0.0) then
       if (patmp2(1).le.0.0) then
        xmin=xmp2(1)
        xmax=xmin
        ymin=ymp2(1)
        ymax=ymin
        do i=1,magpri67
          xmpz(i)=xmp2(i)
          ympz(i)=ymp2(i)
          xmin=min(xmin,xmp2(i))
          xmax=max(xmax,xmp2(i))
          ymin=min(ymin,ymp2(i))
          ymax=max(ymax,ymp2(i))
        enddo
        xtest=(xmin+xmax)/2.
        ytest=(ymin+ymax)/2.
        nzz=0
        call packps(xmpz,ympz,magpri67,xtest,ytest,nzz)
        do i=1,magpri67
          do k=1,magpri67
            if ((xmp2(k).eq.xmpz(i)).and.(ymp2(k).eq.ympz(i))) &
             ampz(i)=amp2(k)
          enddo
        enddo
        do i=1,magpri67
          ip1=i+1
          im1=i-1
          if (i.eq.1) im1=magpri67
          if (i.eq.magpri67) ip1=1
          dang90=abs(abs(ampz(i))-90.)
          dang270=abs(abs(ampz(i))-270.)
          danm90=abs(abs(ampz(im1))-90.)
          danm270=abs(abs(ampz(im1))-270.)
          danp90=abs(abs(ampz(ip1))-90.)
          danp270=abs(abs(ampz(ip1))-270.)
          if (dang90.lt.1.0.or.dang270.lt.1.0) then
           xxm=xmpz(i)
           xxp=xmpz(i)
           if (danm90.lt.1.0.or.danm270.lt.1.0) then
            yym=(ympz(i)+ympz(im1))/2.
           else
            sm1=tan(ampz(im1))*180./pi
            yym=sm1*(xxm-xmpz(im1))+ympz(im1)
           endif
           if (danp90.lt.1.0.or.danp270.lt.1.0) then
            yyp=(ympz(i)+ympz(ip1))/2.
           else
            sm2=tan(ampz(ip1))*180./pi
            yym=sm2*(xxm-xmpz(ip1))+ympz(ip1)
           endif
          else
           if (danm90.lt.1.0.or.danm270.lt.1.0) then
            xxm=xmpz(im1)
            sm2=tan(ampz(i))*180./pi
            yym=sm2*(xxm-xmpz(i))+ympz(i)
           else
            dampz1=abs(ampz(im1)-ampz(i))
            dampz2=abs(dampz1-360.)
            if (dampz1.lt.1.0.or.dampz2.lt.1.0) then
             xxm=(xmpz(i)+xmpz(im1))/2.
             yym=(ympz(i)+ympz(im1))/2.
            else
             sm1=tan(ampz(im1))*180./pi
             sm2=tan(ampz(i))*180./pi
             xxm=(sm1*xmpz(im1)-sm2*xmpz(i)-ympz(im1)+ympz(i))/(sm1-sm2)
             yym=sm1*(xxm-xmpz(im1))+ympz(im1)
            endif
           endif
           if (danp90.lt.1.0.or.danp270.lt.1.0) then
            xxp=xmpz(ip1)
            sm1=tan(ampz(i))*180./pi
            yyp=sm1*(xxp-xmpz(i))+ympz(i)
           else
            dampz1=abs(ampz(ip1)-ampz(i))
            dampz2=abs(dampz1-360.)
            if (dampz1.lt.1.0.or.dampz2.lt.1.0) then
             xxp=(xmpz(i)+xmpz(ip1))/2.
             yyp=(ympz(i)+ympz(ip1))/2.
            else
             sm1=tan(ampz(i))*180./pi
             sm2=tan(ampz(ip1))*180./pi
             xxp=(sm1*xmpz(i)-sm2*xmpz(ip1)-ympz(i)+ympz(ip1))/(sm1-sm2)
             yyp=sm1*(xxp-xmpz(i))+ympz(i)
            endif
           endif
          endif
          patmpz(i)=sqrt((xxp-xxm)**2+(yyp-yym)**2)
        enddo
        do i=1,magpri67
          do k=1,magpri67
            if ((xmpz(k).eq.xmp2(i)).and.(ympz(k).eq.ymp2(i))) &
             patmp2(i)=patmpz(k)
          enddo
        enddo
       endif
       cipmp2=0.0
       do i=1,magpri67
         cipmp2=cipmp2+cmpr2(i,it)*patmp2(i)
       enddo
       cipmp2=cipmp2/tmu/twopi
      endif
!-----------------------------------------------------------------
!--   322 degree probes                                         --
!-----------------------------------------------------------------
      mb=magpri67+1
      mbb=magpri322
      if (xmp2(mb).gt.0.0) then
       if (patmp2(mb).le.0.0) then
        xmin=xmp2(mb)
        xmax=xmin
        ymin=ymp2(mb)
        ymax=ymin
        do i=mb,magpri67+magpri322
          xmpz(i)=xmp2(i)
          ympz(i)=ymp2(i)
          xmin=min(xmin,xmp2(i))
          xmax=max(xmax,xmp2(i))
          ymin=min(ymin,ymp2(i))
          ymax=max(ymax,ymp2(i))
        enddo
        xtest=(xmin+xmax)/2.
        ytest=(ymin+ymax)/2.
        nzz=0
        call packps(xmpz(mb),ympz(mb),mbb,xtest,ytest,nzz)
        do i=mb,magpri67+magpri322
          do k=mb,magpri67+magpri322
            if ((xmp2(k).eq.xmpz(i)).and.(ymp2(k).eq.ympz(i))) &
             ampz(i)=amp2(k)
          enddo
        enddo
        do i=mb,magpri67+magpri322
          ip1=i+1
          im1=i-1
          if (i.eq.mb) im1=magpri67+magpri322
          if (i.eq.magpri67+magpri322) ip1=mb
          dang90=abs(abs(ampz(i))-90.)
          dang270=abs(abs(ampz(i))-270.)
          danm90=abs(abs(ampz(im1))-90.)
          danm270=abs(abs(ampz(im1))-270.)
          danp90=abs(abs(ampz(ip1))-90.)
          danp270=abs(abs(ampz(ip1))-270.)
          if (dang90.lt.1.0.or.dang270.lt.1.0) then
           xxm=xmpz(i)
           xxp=xmpz(i)
           if (danm90.lt.1.0.or.danm270.lt.1.0) then
            yym=(ympz(i)+ympz(im1))/2.
           else
            sm1=tan(ampz(im1))*180./pi
            yym=sm1*(xxm-xmpz(im1))+ympz(im1)
           endif
           if (danp90.lt.1.0.or.danp270.lt.1.0) then
            yyp=(ympz(i)+ympz(ip1))/2.
           else
            sm2=tan(ampz(ip1))*180./pi
            yym=sm2*(xxm-xmpz(ip1))+ympz(ip1)
           endif
          else
           if (danm90.lt.1.0.or.danm270.lt.1.0) then
            xxm=xmpz(im1)
            sm2=tan(ampz(i))*180./pi
            yym=sm2*(xxm-xmpz(i))+ympz(i)
           else
            dampz1=abs(ampz(im1)-ampz(i))
            dampz2=abs(dampz1-360.)
            if (dampz1.lt.1.0.or.dampz2.lt.1.0) then
             xxm=(xmpz(i)+xmpz(im1))/2.
             yym=(ympz(i)+ympz(im1))/2.
            else
             sm1=tan(ampz(im1))*180./pi
             sm2=tan(ampz(i))*180./pi
             xxm=(sm1*xmpz(im1)-sm2*xmpz(i)-ympz(im1)+ympz(i))/(sm1-sm2)
             yym=sm1*(xxm-xmpz(im1))+ympz(im1)
            endif
           endif
           if (danp90.lt.1.0.or.danp270.lt.1.0) then
            xxp=xmpz(ip1)
            sm1=tan(ampz(i))*180./pi
            yyp=sm1*(xxp-xmpz(i))+ympz(i)
           else
            dampz1=abs(ampz(ip1)-ampz(i))
            dampz2=abs(dampz1-360.)
            if (dampz1.lt.1.0.or.dampz2.lt.1.0) then
             xxp=(xmpz(i)+xmpz(ip1))/2.
             yyp=(ympz(i)+ympz(ip1))/2.
            else
             sm1=tan(ampz(i))*180./pi
             sm2=tan(ampz(ip1))*180./pi
             xxp=(sm1*xmpz(i)-sm2*xmpz(ip1)-ympz(i)+ympz(ip1))/(sm1-sm2)
             yyp=sm1*(xxp-xmpz(i))+ympz(i)
            endif
           endif
          endif
          patmpz(i)=sqrt((xxp-xxm)**2+(yyp-yym)**2)
        enddo
        do i=mb,magpri67+magpri322
          do k=mb,magpri67+magpri322
            if ((xmpz(k).eq.xmp2(i)).and.(ympz(k).eq.ymp2(i))) &
             patmp2(i)=patmpz(k)
          enddo
        enddo
!
        if (rank == 0) then
          open(unit=80,status='old',file='mhdin.new',iostat=ioerr)
          if (ioerr.eq.0) close(unit=80,status='delete')
          open(unit=80,status='new',file='mhdin.new',delim='quote')
          write (80,in3)
          close(unit=80)
        endif
       endif
!
       cipmp3=0.0
       do i=magpri67+1,magpri67+magpri322
         cipmp3=cipmp3+cmpr2(i,it)*patmp2(i)
       enddo
       cipmp3=cipmp3/tmu/twopi
      endif
      endif newmhdin
!
      if (.not.fitsiref) then
       ssiref=csilop(iabs(nslref),it)
       do i=1,nsilop
         workb_jw4(i)=csilop(i,it)-ssiref
       enddo
      else
       ssiref=0.0
       do i=1,nsilop
         workb_jw4(i)=csilop(i,it)+csiref
       enddo
      endif
!
      if (iand(iout,1).ne.0) then
!
        write (nout,11100) ssiref
        write (nout,11020) (workb_jw4(i),i=1,nsilop)
        write (nout,11100) csiref
        write (nout,11020) (csilop(i,it),i=1,nsilop)
        workb_jw4(1:nsilop)=csilop(1:nsilop,it)*twopi
        write (nout,11101) csiref*twopi
        write (nout,11020) (workb_jw4(i),i=1,nsilop)
        write (nout,11102)
        write (nout,11020) (csilopv(i,it),i=1,nsilop)
        write (nout,11120) cipmp2,cipmp3
        write (nout,11020) (cmpr2(i,it),i=1,magpri)
        write (nout,11122)
        write (nout,11020) (cmpr2v(i,it),i=1,magpri)
        if (kstark.gt.0) then
          write (nout,11140)
          write (nout,11020) (cmgam(i,it),i=1,nstark)
        endif
        write (nout,11160) ipmhd(it)
        write (nout,11180) cdflux(it),sbpp,delbp(it),sbppa
        if (kecebz.gt.0) write(nout,11185) cmecebz(it)
        if (kece.gt.0)  then
              write(nout,11186)
              write(nout,11020)  (cmece(m,it),m=1,nece)
        endif
!
        write (nout,11200) psiref(it)
        write (nout,11020) (silopt(it,i),i=1,nsilop)
        write (nout,11220)
        write (nout,11020) (expmpi(it,i),i=1,magpri)
        if (kstark.gt.0) then
          write (nout,11240)
          write (nout,11020) (tangam(it,i),i=1,nstark)
        endif
        write (nout,11260) ipmeas(it)
        write (nout,11280) diamag(it)
        write (nout,11270) vloopt(it)
        write (nout,11292)
        write (nout,11020) (fccurt(it,i),i=1,nfsum)
        write (nout,11294)
        write (nout,11020) (eccurt(it,i),i=1,nesum)
!
      if (abs(sigdia(it)).gt.1.0e-08) &
        write (nout,11300) chidflux
      if (iconvr.eq.3) then
        write (nout,11320) emf,emp,enf,enp,betap0,rzero
        cbetap=0.0
        cli=0.0
        cqqxis=0.0
        cbetat=0.0
        ci0=0.0
        write (nout,11330) cbetap,cli,cqqxis,cbetat,ci0 ! all unset...
      endif
      if (nbdry.gt.0) then
        write (nout,11324) erbmax,erbave
        write (nout,11326) (erbloc(i),i=1,nbdry)
      endif
!
      if (kvtor.gt.0) then
        write (nout,13000)
        write (nout,13020) betatw(it),betapw(it),wplasw(it)
      endif
!
      write (nout,12000)
      do i=1,nitera
        write (nout,12020) i,cerror(i),csibry(i),csimag(i),cvolp(i), &
                  crmaxi(i),czmaxi(i),cemaxi(i),cqmaxi(i),cchisq(i)
      enddo
      if ((kwripre.gt.0).and.(kwripre.le.9)) then
          call setfnmd('n',ishot,itime,sfname)
          sfname=sfname(1:13)//'_error'
          open(unit=74,status='old',file=sfname,iostat=ioerr)
          if (ioerr.eq.0) close(unit=74,status='delete')
          open(unit=74,status='new',file=sfname)
          do i=1,nitera
           write (74,*) i,cerror(i),xdum,xdum
          enddo
          close(unit=74)
      endif
      write (nout,12010)
      do i=1,nitera
        write (nout,12025) i,aveerr(i),csumip(i),tratio(i),iermax(i), &
                  jermax(i)
      enddo
!
      if ((kecebz.gt.0).or.(kece.gt.0)) then
        write (nout,12015)
       do i=1,nitera
        write (nout,12020) i,receoi(i),(recemi(i,k),k=1,nece)
       enddo
        write (nout,12016)
       do i=1,nitera
        write (nout,12020) i,(recepi(i,k),k=1,nece)
       enddo
      endif
!
      dsi=1.0_dp/(nw-1)
      write (nout,12040)
      do i=1,nw
        sinow=dsi*(i-1)
        write (nout,12020) i,sinow,volp(i),pprime(i),curmid(i),ffprim(i) &
                ,pres(i),fpol(i),qpsi(i),rpres(i)
      enddo
!
      write (nout,12043)
      do i=1,nw
        sinow=dsi*(i-1)
        write (nout,12020) i,sinow,cjor(i)
      enddo
!
      endif
      return
 9300 format (/,4x,'   data used:   ')
 9320 format (1x,i2,' flux loops')
 9340 format (1x,i2,' magnetic probes(i)')
 9360 format (1x,i2,' partial rogowskis')
 9380 format (1x,i2,' full rogowski')
 9385 format (1x,i2,' diamagnetic loop')
 9390 format (1x,' bt0(t)   = ',f10.3)
!10000 format(/,6x,20('*'),'EFITAI',a3,' x ',a3,'  output ',20('*'))
10020 format (1x,'  sumif(amp) = ',e10.3,' sumift(amp) = ',e10.3, &
           ' sumifs(amp) = ',e10.3)
10480 format (1x,/)
10500 format(' shot #   = ',i10,' time(ms) = ',i10, &
             ' chi**2   = ',1pe10.3)
10520 format(' betat(%) = ',f10.3,' betap    = ',f10.3, &
             ' li       = ',f10.3)
10540 format(' vol(cm3) = ',1pe10.3,' rout(cm) = ',0pf10.3, &
             ' zout(cm) = ',0pf10.3)
10560 format(' elong    = ',f10.3,' utriang  = ',f10.3, &
             ' ltriang  = ',f10.3)
10580 format(' a(cm)    = ',f10.3,' lin(cm)  = ',f10.3, &
             ' lout(cm) = ',f10.3)
10600 format(' ltop(cm) = ',f10.3,' q*       = ',f10.3, &
             ' rc(cm)   = ',f10.3)
10610 format(' zc(cm)   = ',f10.3,' bt0(t)   = ',f10.3, &
             ' qout     = ',f10.3)
10620 format(' lins(cm) = ',f10.3,' louts(cm)= ',f10.3, &
             ' ltops(cm)= ',f10.3)
10623 format(' beta*(%) = ',f10.3)
10622 format(//, &
      ' betat, betat-btvac, beta-total, beta-btvac2, beta-btv :',/, &
      ' betat0 :',/, &
             5(2x,e12.4,2x),/,1(2x,e12.4,2x))
10624 format(//, &
      ' betatm, betat-vbtvac2, beta-vbtor2, beta-vbtot2:',/, &
             4(2x,e12.4,2x))
11000 format(//,22x,'F-coils currents (Amp)')
11001 format(//,22x,'A matrix condition    ')
11002 format(//,22x,'E-coils phases  ')
11004 format(1x,' info = ',i4,' row = ',1pe10.3,' col = ',1pe10.3, &
                ' max = ',1pe10.3)
11005 format(//,22x,'E-coils currents')
11008 format(//,22x,'E-coils resistance(Ohms)')
11010 format(//,22x,'F-coils resistance(Ohms)')
11017 format(//,2x,'power supply current (A) = ',e12.5, &
                5x,'phase (degree) = ',e12.5,/, &
                2x,'resistance (Ohm)         = ',e12.5, &
                5x,'inductance (H) = ',e12.5)
11020 format(4e15.6)
11022 format(/,12x,' vessel vertical forces (p, newton) ',e15.6)
11024 format(/,12x,' vessel vertical forces (t, newton) ',e15.6)
11030 format(//,22x,'F-coils currents (Amp-turn)')
11032 format(//,12x,'Electrostatic potential derivative PIEPRIM:')
11034 format(//,12x,'Hyperbolic P:',e15.6)
11036 format(//,12x,'Hyperbolic FF:',e15.6)
11040 format(//,22x,'plasma currents','  normalization = ',e15.6)
11043 format(//,22x,'plasma coeffics','  normalization = ',e15.6)
11060 format(//,22x,'vessel currents',' sum(amp) = ',e10.3, &
                '  top   =   ',e10.3,'  bot     = ',e10.3)
11100 format(//,16x,'calculated psi-loops signals', &
        ' ssiref(vs/rad) = ',e12.5)
11101 format(//,16x,'calculated psi-loops signals', &
        ' ssiref(vs) = ',e12.5)
11102 format(//,16x,'calculated vacuum psi-loops signals')
11120 format(//,4x,'calculated magnetic probes(i) signals', &
        ' ipmp2(A) = ',e12.5,1x,e12.5)
11122 format(//,4x,'calculated vacuum magnetic probes(i) signals')
11140 format(//,14x,'calculated polarimetry signals ')
11160 format(//,14x,'calculated total plasma current ',/,16x,e15.6)
11180 format(//,14x,'calculated diamagnetic flux(vs) ',/,16x,e15.6,/, &
           '     bpdia = ',e10.3,'     delbp = ',e10.3, &
           ' app bpdia = ',e10.3)
11185 format(//,14x,'calculated Bz(receo,zeceo) (T)  ',/,16x,e15.6)
11186 format(//,22x,'calculated psi(R-)-psi(R+) (VS/rad)')
11200 format(//,16x,'  measured psi-loops signals', &
        ' psiref(vs/rad) = ',e12.5)
11220 format(//,12x,'  measured magnetic probes(i) signals')
11240 format(//,14x,'  measured polarimetry signals ')
11260 format(//,14x,'  measured total plasma current ',/,16x,e15.6)
11270 format(//,14x,'  measured loop voltage (V)     ',/,16x,e15.6)
11280 format(//,14x,'  measured diamagnetic flux     ',/,16x,e15.6)
11292 format(//,14x,'  measured F-coil currents(A)   ')
11294 format(//,14x,'  measured E-coil currents(A)   ')
11300 format(//,14x,'    chisqr diamagnetic flux  =  ',e10.3)
11320 format (//,'   emf  =   ',e10.3,'   emp  =   ',e10.3, &
           '   enf  =   ',e10.3,'   enp  =   ',e10.3,/, &
           ' betap0 =   ',e10.3,' rzero  =   ',e10.3)
11324 format (//,'  erbmax =  ',e12.5,'   erbave = ',e12.5)
11326 format (8(1x,e12.5))
11330 format (//,'  cbetap =  ',e12.5,'   cli    = ',e12.5, &
                 '  cqqxis =  ',e12.5,'   cbetat = ',e12.5,/, &
                 '  ci0    =  ',e12.5)
12000 format (//,' iteration summary:',/,4h  i ,2x,'   error    ',2x, &
        '   psibry   ',2x,'   psimag   ',2x,'   volume   ',2x, &
        '   rmaxis   ',2x,'   zmaxis   ','   emaxis   ',2x, &
        '   qmaxis   ',2x,'   chisqr   ')
12010 format (//,'  i ',2x,'   errave   ',2x, &
        '   current  ',2x,'   cratio   ',2x,'   iermax   ',2x, &
        '   jermax   ',2x,'            ','            ',2x, &
        '            ',2x,'            ')
12015 format (//,' ECE iteration summary:',/,'  i ',2x,'   receo    ',2x, &
        '   recem    ')
12016 format (//,'  i ',2x,'   recep    ')
12020 format (i4,10(1x,e12.5))
12025 format (i4,3(2x,e12.5),2(6x,i4,4x),4(2x,e12.5))
12040 format (//,'    plasma summary:',/,'  i ',1x,'   pflux    ',1x, &
        '   vol(m3)  ',1x,'   pprime   ',1x,'  current   ',1x, &
        ' ffprime    ',1x,'  pressure  ','   fpol     ',1x, &
        '    q       ',1x,'     rm     ')
12043 format (//,'  i ',1x,'    pflux   ',1x, &
        '   <j>      ',1x,'            ',1x,'            ',1x, &
        '            ',1x,'            ','            ',1x, &
        '            ',1x,'            ')
13000 format (//,16x,'    toroidal rotation         ')
13020 format ('  betatw =  ',1pe12.5,' betapw   = ',1pe12.5, &
              '  W(J)   =  ',1pe12.5)
      end subroutine print_stats

!**********************************************************************
!>
!!    Prints header  
!!
!**********************************************************************
      subroutine print_header()
      include 'eparm.inc'
      include 'modules1.inc'

      if (itek.le.0) then
        if(rank == 0) write(nttyo,10001) nw,nh
      endif

      ! if iout is odd
      if(iand(iout,1).ne.0) write(nout,10000) trim(ch1),trim(ch2)

10000 format(/,6x,20('*'),'EFIT  ',a3,' x ',a3,'  output ',20('*'))
10001 format(1x,20('*'),' EFIT ',i0,' x ',i0,' grid ',20('*'))
      end subroutine print_header
