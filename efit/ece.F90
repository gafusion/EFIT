#include "config.f"
!**********************************************************************
!>
!!    setece obtains the response matrix
!!    for ECE measurement
!!    
!!
!!    @param jtime : time index
!!
!!    @param kerror : error flag
!!
!********************************************************************* 
      subroutine setece(jtime,kerror)
      use commonblocks,only: c,wk,bkx,bky
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: jtime
      integer*4, intent(out) :: kerror
      integer*4 i,idesto,ier,ii,iio,iname,ioerr,isplit,itot,iw, &
                j,jh,jj,k,kk,kkm,kgeteceb,km,kp,ksetece,l,m,n,nj,nk,nnnece
      real*8 a,bz,bzct,ddm,ddo,ddp,ddsiddro,dho,dist,distt,drzg,fitot, &
             gbzt,psical,rrreceo,r,rdif,rh,rsum,rselsum,rw,r1, &
             tdata,tdata1,tdata2,z,zdif,zzx,z1
      real*8 rrrecep(nnece),rrrecem(nnece)
      real*8 zece(nnece),pds(6),rsplt(2500),zsplt(2500),csplt(2500)
      real*8,dimension(:),allocatable :: dsidr,ddsiddr
      character*40 filenmme
      logical isopen
      data ksetece/0/,kgeteceb/0/
      integer*4, parameter :: nset=20
      real*8, parameter :: cdum=1.0

      kerror = 0
!
#ifdef DEBUG_LEVEL3
      write (6,*) 'Enter SETECE, nw = ',nw
#endif
      ALLOCATE(dsidr(nw),ddsiddr(nw))
!
#ifdef DEBUG_LEVEL3
      write (6,*) 'Enter SETECE, ksetece = ',ksetece
#endif
      kece=0
      kecebz=0
!
      ksetece0: if (ksetece.eq.0) then
!------------------------------------------------------------------
!--   zeceo from k-file, jo--zeceo on zgrid
!------------------------------------------------------------------
      zeceo=zteo
      receo=rteo
      dho=abs(zgrid(1)-zeceo)
      jo=1
      do jh=2,nh
        if (abs(zgrid(jh)-zeceo).lt.dho) then
          dho=abs(zgrid(jh)-zeceo)
          jo=jh
        endif
      enddo
#ifdef DEBUG_LEVEL3
      write (6,*) 'SETECE, receo/zeteo = ',receo,zeceo
#endif
!-------------------------------------------------------------------
!--   kfixro=1, receo  from k-file
!------------------------------------------------------------------- 
      if(kfixro.eq.1) receo=rteo
      if ((kfixro.ne.1).or.(kfitece.ne.1)) then
!------------------------------------------------------------------
!--     kfixrece=1, R+(recep) R-(recem) from k-file
!------------------------------------------------------------------
        if (kfixrece.eq.1) then
          do k=1,nnece
            if ((rtep(k).gt.1.E-5_dp).and.(rtem(k).gt.1.E-5_dp)) then
              recep(k)=rtep(k)
              recem(k)=rtem(k)
              nece=k
            endif
          enddo
        endif
!-----------------------------------------------------------------
!--     kfitece=1 or 2, receo or R+ R- from ECE data
!--     kfitece=3, receo or R+ R- from gettir
!-----------------------------------------------------------------
        if (kfitece.le.2 .and. nfit.gt.0) then
!-----------------------------------------------------------------
!--       kfixro=0 or kfixrece=0, receo or R+ R- from getecer
!-----------------------------------------------------------------
          if ((kfixro.eq.0).or.(kfixrece.eq.0)) then
            call getecer(jtime,kerror)
            if(kerror.gt.0) return
          endif
!-----------------------------------------------------------------
!--       kfixro=-1 or kfixrece=-1, receo or R+ R- from geteceb
!-----------------------------------------------------------------
          if ((kfixro.eq.-1).or.(kfixrece.eq.-1)) then
            call geteceb(jtime,kerror)
            if(kerror.gt.0) return
          endif
        elseif (kfitece.eq.3) then
          call gettir(jtime,kerror)
          if(kerror.gt.0) return
        endif
      endif
!----------------------------------------------------------------
!--   get iwo iwp iwm (receo, R+ R- on rgrid)
!----------------------------------------------------------------
      ddo=abs(rgrid(1)-receo)
      iwo=1
      do iw=2,nw
        if (abs(rgrid(iw)-receo).lt.ddo) then
          ddo=abs(rgrid(iw)-receo)
          iwo=iw
        endif
      enddo
      if (kfitece.ne.1) then
        zece(1:nece)=zeceo
        do k=1,nece
          ddp=abs(rgrid(1)-recep(k))
          ddm=abs(rgrid(1)-recem(k))
          iwp(k)=1
          iwm(k)=1
          do iw=2,nw
            if (abs(rgrid(iw)-recep(k)).lt.ddp) then
              ddp=abs(rgrid(iw)-recep(k))
              iwp(k)=iw
            endif
          enddo
          do iw=2,nw
            if (abs(rgrid(iw)-recem(k)).lt.ddm) then
              ddm=abs(rgrid(iw)-recem(k))
              iwm(k)=iw
            endif
          enddo
        enddo  
      endif
!-------------------------------------------------------------------
!--   try to read response function from recefile
!--------------------------------------------------------------------
      do iname=1,nset
        if (iname.le.9) then
          write(filenmme,6530) iname
        else
          write(filenmme,6540) iname
        endif
        isopen=.false.
        open(unit=nffile, status='old',form='unformatted', &
             file=filenmme,iostat=ioerr)
        if (ioerr.ne.0) then
          open(unit=nffile, status='old',form='unformatted', &
               file=table_dir(1:ltbdir)//filenmme,iostat=ioerr)
          if(ioerr.ne.0) exit
        endif
        read(nffile,iostat=ioerr) nnnece
        if(ioerr.ne.0) exit
        if(nnnece.ne.nece) then
          close(unit=nffile)
          cycle
        endif
        read(nffile,iostat=ioerr) rrrecem
        if(ioerr.ne.0) exit
        read(nffile,iostat=ioerr) rrrecep
        if(ioerr.ne.0) exit
        read(nffile,iostat=ioerr) rrreceo
        if(ioerr.ne.0) exit
        do ii=1,nece
          if(abs(rrrecem(ii)-recem(ii)).gt.1.e-4_dp) cycle
          if(abs(rrrecep(ii)-recep(ii)).gt.1.e-4_dp) cycle
        enddo
        if(abs(rrreceo-receo).gt.1.e-4_dp) cycle
        read(nffile,iostat=ioerr) recebzfc
        if(ioerr.ne.0) exit
        read(nffile,iostat=ioerr) gecebzpc
        if(ioerr.ne.0) exit
        read(nffile,iostat=ioerr) recebzec
        if(ioerr.ne.0) exit
        read(nffile,iostat=ioerr) recefc
        if(ioerr.ne.0) exit
        read(nffile,iostat=ioerr) gecepc
        if(ioerr.ne.0) exit
        read(nffile,iostat=ioerr) receec
        if(ioerr.ne.0) exit
        close(unit=nffile)
        return
      enddo
      if(isopen) close(unit=nffile)
!-------------------------------------------------------------------       
!--   compute the response function about Te peak point constraint    
!--    response due to F coils --- recebzfc(ncoil)
!-------------------------------------------------------------------
#ifdef DEBUG_LEVEL3
      write (6,*) 'SETECE, rf/zf = ',rf(1),zf(1)
#endif
      isplit=10
      itot=isplit*isplit
      fitot=itot 
      recebzfc = 0.0
      do k=1,nfcoil
        bzct=0 
        call splitc(isplit,rsplt,zsplt,csplt, &
                    rf(k),zf(k),wf(k),hf(k),af(k),af2(k),cdum)
        r1=receo
        do l=1,itot
          a=rsplt(l)  
          z1=zeceo-zsplt(l)
          bzct=bzct+bz(a,r1,z1)      
        enddo
        kkm=fcid(k)
        recebzfc(kkm)=recebzfc(kkm)+fcturn(k)*bzct/fitot*tmu
      enddo
!---------------------------------------------------------------------
!--   plasma response   gecebzpc(nwnh)                              --
!---------------------------------------------------------------------
      zzx=0.0
      drzg=(drgrid**2+dzgrid**2)*25.
      isplit=18
      if (receo.gt.1.e-8_dp) then
        r=receo
        do ii=1,nw
          a=rgrid(ii)
          do jj=1,nh
            z=zeceo-zgrid(jj)
            kk=(ii-1)*nh+jj
            dist=(a-r)**2+z**2
            if (dist.lt.drzg) then
   68         call splitc(isplit,rsplt,zsplt,csplt, &
                    rgrid(ii),zgrid(jj),drgrid,dzgrid,zzx,zzx,cdum)
              gbzt=0.0
              itot=isplit*isplit
              do k=1,itot
                a=rsplt(k)
                z=zeceo-zsplt(k)
                distt=(a-r)**2+z**2
                if (distt.lt.1.e-8_dp.and.isplit.lt.49) then
                  isplit=isplit+2
                  go to 68
                endif
                gbzt=gbzt+bz(a,r,z)
              enddo
              gecebzpc(kk)=gbzt*tmu/itot
            else
              gecebzpc(kk)=bz(a,r,z)*tmu
            endif
          enddo
        enddo
      endif
!---------------------------------------------------------------------
!--   E coils response   recebzec                                   --
!---------------------------------------------------------------------
      if (iecurr.gt.0) then
        isplit=10
        itot=isplit*isplit
        fitot=real(itot,dp)
        do i=1,nesum
          recebzec(i)=0.0
        enddo
        if (receo.gt.1.e-8_dp) then
          r1=receo
          do k=1,necoil
            bzct=0.0
            call splitc(isplit,rsplt,zsplt,csplt, &
                        re(k),ze(k),we(k),he(k),zzx,zzx,cdum)
            do l=1,itot
              a=rsplt(l)
              z1=zeceo-zsplt(l)      
              bzct=bzct+bz(a,r1,z1)
            enddo
            bzct=bzct*tmu/fitot
            kkm=ecid(k)  
            recebzec(kkm)=recebzec(kkm)+bzct
          enddo
        endif
      endif
!-------------------------------------------------------------------------
!--   compute the response function about R+ R- constraint     
!--   F coils response -----recefc(nece,nfsum)
!-------------------------------------------------------------------------
      ECE: if (kfitece.ne.1) then
      do n=1,nfsum
        call sets2d(gridfc(1,n),c,rgrid,nw,bkx,lkx,zgrid,nh,bky, &
                       lky,wk,ier)
        do i=1,nece
          if (eceiter.eq.'pair') then
            call seva2d(bkx,lkx,bky,lky,c,recem(i),zece(i),pds,ier,n111)
            recefc(i,n)=pds(1)
          endif
          if (eceiter.eq.'flux') recefc(i,n)=0.0
        enddo
        do i=1,nece
          call seva2d(bkx,lkx,bky,lky,c,recep(i),zece(i),pds,ier,n111)
          recefc(i,n)=pds(1)-recefc(i,n)
        enddo
      enddo
!------------------------------------------------------------------------------
!--   plasma response -----gecempc(nece,nwnh),geceppc(nece,nwnh)             --
!--                 R-, R+ flux from pc,gecepc=geceppc-gecempc               --
!------------------------------------------------------------------------------
      if (eceiter.eq.'pair') then
        do m=1,nece
          do i=1,nw
            do j=1,nh
              k=(i-1)*nh+j
              rdif=recem(m)-rgrid(i)
              zdif=zece(m)-zgrid(j)
              rsum=rdif**2+zdif**2
              if (rsum.le.dselsum) then
!                mk=(i-1)*nh+1
!                gecempc(m,k)=gridpc(mk,i)
                zdif=dselsum
                rselsum=rgrid(i)-dselsum
                gecempc(m,k)=psical(recem(m),rselsum,zdif)*tmu
              else
                gecempc(m,k)=psical(recem(m),rgrid(i),zdif)*tmu
              endif
            enddo
          enddo
        enddo
      endif
      if (eceiter.eq.'flux') gecempc = 0.0
      do m=1,nece
        do i=1,nw
          do j=1,nh
            k=(i-1)*nh+j
            rdif=recep(m)-rgrid(i)
            zdif=zece(m)-zgrid(j)
            rsum=rdif**2+zdif**2
            if (rsum.le.dselsum) then
!              mk=(i-1)*nh+1
!              geceppc(m,k)=gridpc(mk,i)
              zdif=dselsum
              rselsum=rgrid(i)-dselsum
              geceppc(m,k)=psical(recep(m),rselsum,zdif)*tmu
            else
              geceppc(m,k)=psical(recep(m),rgrid(i),zdif)*tmu
            endif
            gecepc(m,k)=geceppc(m,k)-gecempc(m,k)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
!--   Ohmic coils   receec(nece,nesum)                                --
!-----------------------------------------------------------------------
      ohmic: if (iecurr.gt.0) then
        do n=1,nesum
          call sets2d(gridec(1,n),c,rgrid,nw,bkx,lkx,zgrid,nh,bky, &
                      lky,wk,ier)
          if (eceiter.eq.'pair') then
            do i=1,nece
              call seva2d(bkx,lkx,bky,lky,c,recem(i),zece(i),pds,ier,n111)
              receec(i,n)=pds(1)
            enddo
          endif
          if (eceiter.eq.'flux') then
            do i=1,nece
              receec(i,n)=0.0
            enddo
          endif
          do i=1,nece
            call seva2d(bkx,lkx,bky,lky,c,recep(i),zece(i),pds,ier,n111)
            receec(i,n)=pds(1)-receec(i,n)
          enddo
        enddo
      endif ohmic
      endif ECE
!----------------------------------------------------------------------
      if (rank.eq.0) then
        open(unit=nffile,status='old',form='unformatted',iostat=ioerr, &
             file='recexx.dat')
        if(ioerr.eq.0) close(unit=nffile,status='delete')
        open(unit=nffile,status='new',form='unformatted', &
             file='recexx.dat')
        nnnece=nece
        write(nffile) nnnece 
        do ii=1,nece
          rrrecem(ii)=recem(ii)
          rrrecep(ii)=recep(ii)
        enddo
        rrreceo=receo
        write(nffile) rrrecem
        write(nffile) rrrecep
        write(nffile) rrreceo
        write(nffile) recebzfc
        write(nffile) gecebzpc
        write(nffile) recebzec
        write(nffile) recefc
        write(nffile) gecepc 
        write(nffile) receec
        close(unit=nffile)
      endif
      endif ksetece0
!-----------------------------------------------------------------------
!--   do every time from here
!-----------------------------------------------------------------------
      ksetece=ksetece+1
      if (ksetece.eq.mtxece) ksetece=0
!-------------------------------------------------------------------
!--   get d(psi)/dr and d2(psi)/dr2
!-------------------------------------------------------------------
      call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
      do iw=1,nw
        rw=rgrid(iw)
        rh=zgrid(jo)
        kk=(iw-1)*nh+jo
        call seva2d(bkx,lkx,bky,lky,c,rw,rh,pds,ier,n555)
        if (ier.ne.0) then
          write (nttyo,9910) ier,rw,rh
          return
        endif
        dsidr(iw)=pds(2)
        ddsiddr(iw)=pds(5)
      enddo
!
      if (eceiter.eq.'flux') then
        do i=1,nnece
          rw=recem(i)
          rh=zece(i)
          call seva2d(bkx,lkx,bky,lky,c,rw,rh,pds,ier,n555)
          if (ier.ne.0) then
            write (nttyo,9910) ier,rw,rh
            return
          endif
          brspece(jtime,i)=pds(1)
        enddo
#ifdef DEBUG_LEVEL3
        write (6,*) 'SETECE, recem/brspece = ',(recem(ii), &
            brspece(jtime,ii),ii=1,nece)
#endif
      endif
 9910 format('  error in spline =',i4,' (r,z)= ( ',f5.2,',',f5.2,')')
!----------------------------------------------------------------------
!--   give error (ecebit,ecebzbit)                                   --
!--   if kfixro kfixrece=0, robit and ecebit  from getecer           --
!--   if kfixro kfixrece=1, robit and rmbit, rpbit from k-file       --
!--     ecebit=sqrt((dpsi/dR|m*rmbit)**2+(dpsi/dR|p*rpbit)**2)       --
!--     ecebzbit=(d2psi/dR2)*robit/receo                             --
!----------------------------------------------------------------------
      if (kfitece.ne.1) then
        if (kfixrece.eq.1) then
          do k=1,nece
            kp=iwp(k)
            km=iwm(k)
            ecebit(k)=(dsidr(kp)*rpbit(k))**2
            ecebit(k)=ecebit(k)+(dsidr(km)*rmbit(k))**2
            ecebit(k)=sqrt(ecebit(k))
          enddo
        endif
      endif
      ddsiddro=ddsiddr(iwo)
      ecebzbit=ddsiddro*robit/receo
      if (kfitece.ne.1) then
        do m=1,nece
!          tdata1=serror*abs(brspece(jtime,m))
          tdata1=eceerror*abs(brspece(jtime,m))
          tdata2=abs(ecebit(m))
          tdata=max(tdata1,tdata2)
          if(tdata.gt.1.0e-10_dp) fwtece(m)=fwtece0(m)/tdata
          if(tdata.le.1.0e-10_dp) fwtece(m)=0.0
        enddo
        do i=1,nece
          if (fwtece(i).gt.0.0) kece=kece+1
        enddo
      endif
      tdata1=serror*abs(brspecebz(jtime))
      tdata2=abs(ecebzbit)
      tdata=max(tdata1,tdata2)
      if(tdata.gt.1.0e-10_dp) fwtecebz=fwtecebz0/tdata
      if(tdata.le.1.0e-10_dp) fwtecebz=0.0
      if(fwtecebz.gt.0.0) kecebz=kecebz+1
      receoi(nitera)=receo
      do k=1,nece
         if (fwtece(k).gt.1.e-10_dp) then
           recemi(nitera,k)=recem(k)
           recepi(nitera,k)=recep(k)
         else
           recemi(nitera,k)=1.e-10_dp
           recepi(nitera,k)=1.e-10_dp
         endif
      enddo
!
#ifdef DEBUG_LEVEL3
      write (6,*) 'SETECE, serror/eceerror = ',serror,eceerror
      write (6,*) 'SETECE, nece/fwtece = ',nece,(fwtece(i),i=1,nece)
#endif

      DEALLOCATE(dsidr,ddsiddr)
!
      return
 6530 format ('recefile_0',i1,'.ddd')
 6540 format ('recefile_',i2,'.ddd')
      end subroutine setece

!**********************************************************************
!>
!!    geteceb obtains the receo, R+ R-
!!    from ECE measurement data, (fitting T(B))
!!    if kfixro kfixrece = -1, called in setece
!!
!!    @param jtime : time index
!!    @param kerror : error flag
!**********************************************************************
      subroutine geteceb(jtime,kerror)
      use commonblocks,only: c,wk,bkx,bky
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: jtime
      integer*4, intent(out) :: kerror
      real*8 fpcurr,seval
      integer*4, parameter :: nnn1=1,kbre=5
      integer*4 i,idesto,ier,ii,iieerr,iio,imk,iout1,ioutk,ioutk1,ipk,iw, &
                j,k,kk,kgeteceb,kmin,kmax,m,mnow,n,nj,nk,nlowf
      real*8 baa,bbx,bbxm,bbxp,bbx1,bbx2,beceo,binmax,binmin,bitm,bitp, &
             bobit,b00,dbbf,dbdro,delsi,dest,desto,dsidrm,dsidrp,dtero, &
             fnow,fwtcm,ppteppbo,ptpr1,ptpr2,rh,rw,siii,sumf,t,teeceo,toler
      real*8 pds(6),nnout(kbre),bmink(kbre),bmaxk(kbre)
      real*8 arspfit(nnecein,nfit),brspfit(nnecein), &
             s(nfit),tte(nnecein),x(nfit),an(nnecein), &
             tebit(nnecein)
      real*8 telowf(nnnte),blowf(nnnte),bb(nnnte),cc(nnnte),dd(nnnte), &
             teece(nnece),pteprm(nnece),pteprp(nnece), &
             idestp(nnece),idestm(nnece),becem(nnece),becep(nnece), &
             dbdrp(nnece),dbdrm(nnece)
      real*8 rrgrid(kbre,nw),bfield(nw),rrout(kbre,nw), &
             bout(kbre,nw),babs(kbre,nw),bbb(nw),ccc(nw), &
             ddd(nw),btttt(nw),dsidr(nw),ddsiddr(nw),bx(nw), &
             ry(nw),bbk(nw),dbdr(nw)
      data kgeteceb/0/

      kerror = 0
!-------------------------------------------------------------------
      telowf=0
      blowf=0;bb=0;cc=0;dd=0
      teece=0;pteprm=0;pteprp=0
      idestp=0;idestm=0;becem=0;becep=0
      dbdrp=0;dbdrm=0
      fwtece0=swtece
      fwtecebz0=swtecebz
      babs=0.0
      bout=0.0
      rrout=0.0
      rrgrid=0.0
!-------------------------------------------------------------------
!--   input kgeteceb=0 from input file
!-------------------------------------------------------------------
      getECE: if (kgeteceb.le.0) then
      kgeteceb=kgeteceb+1
!---------------------------------------------------------------------
!--   Calculation of |B| array from fe array (harmonic nharm)       --
!--     becein(necein), fe(GHz), |B|(T)                             --
!--     !!! becein from low field to high field !!!                 --
!---------------------------------------------------------------------
      becein(1:necein)=0.001_dp*6.0*9.1095_dp*pi/4.8032_dp &
                      *feece0(1:necein)/nharm
!EALW      write(*,*)'becein'
!EALW      write(*,*)becein
!--------------------------------------------------------------------
!--   fitting data from teecein0,errorece0 and becein (nnecein)    --
!--     bbx=(B-b00)/baa                                            --
!--     Te=x(1)+x(2)*bbx+x(3)*bbx**2+...+x(nfit)*bbx**(nfit-1)     --
!--------------------------------------------------------------------
      binmin=becein(1)
      binmax=becein(necein)
      baa=0.5_dp*(binmax-binmin)
      b00=0.5_dp*(binmax+binmin)
      an(1:necein)=(becein(1:necein)-b00)/baa
      do nj=1,necein
        tebit(nj)=max(errorece0(nj),1.e-4_dp)
        do nk=1,nfit
          if (nk.eq.1) then
            arspfit(nj,nk)=1./tebit(nj)
          else
            arspfit(nj,nk)=an(nj)**(nk-1)/tebit(nj)
          endif
        enddo
      enddo
!---------------------------------------------------------------------
!--   teecein0,errorece0  from low field to high field
!---------------------------------------------------------------------
      brspfit(1:necein)=teecein0(1:necein)/tebit(1:necein)
!
      mnow=necein
      if (kcmin.gt.0) then
        fwtnow=0.001_dp ! input option, why hardcoded here?
        fwtcm=1.0
        do j=1,nfit
          mnow=mnow+1
          do k=1,nfit
            if (j.ne.k) then
              arspfit(necein+j,k)=0.0
            else
              arspfit(necein+j,k)=fwtcm/fwtnow
            endif
          enddo
          brspfit(necein+j)=0.0
        enddo
      endif
!
      iieerr=0
      call sdecm(arspfit,nnecein,mnow,nfit,brspfit,nnecein,nnn1,s,wk,iieerr)
      if (iieerr.eq.129) then
        kerror=1
        call errctrl_msg('geteceb','sdecm failed to converge')
        return
      endif
      toler=1.0e-06_dp*s(1)
      do i=1,nfit
        t=0.0
        if(s(i).gt.toler) t=brspfit(i)/s(i)
        brspfit(i)=t
      enddo
      do i=1,nfit
        x(i)=sum(arspfit(i,1:nfit)*brspfit(1:nfit))
      enddo
      xfit(1:nfit)=x(1:nfit)
!EALW      write(*,*)'x'
!EALW      write(*,*)x
      chisqfit=0.0
      do k=1,necein
        tte(k)=0.
        do nk=1,nfit
          tte(k)=tte(k)+x(nk)*an(k)**(nk-1)
        enddo
        chisqfit=chisqfit+(tte(k)-teecein0(k))**2/tebit(k)
      enddo
!EALW      write(*,*) 'chisqfit='
!EALW      write(*,*) chisqfit
!EALW      write(*,*)'tte'
!EALW      write(*,*)tte
!--------------------------------------------------------------------
!--   get Teeceb(bbf) in ECE data region, (Te(B)), bbf-B           --
!--------------------------------------------------------------------
      dbbf=(becein(necein)-becein(1))/(nnnte-1)
      do i=1,nnnte
        bbf(i)=becein(1)+dbbf*(i-1)
        bbx=(bbf(i)-b00)/baa
        teeceb(i)=0.
        do nk=1,nfit
          teeceb(i)=teeceb(i)+x(nk)*bbx**(nk-1)
        enddo
      enddo
!---------------------------------------------------------------------
!--   find beceo which is the B value of Te peak point             --
!---------------------------------------------------------------------
      fixro: if (kfixro.ne.1) then
        teeceo=teeceb(1)
        iio=1
        do i=2,nnnte
          if (teeceb(i).gt.teeceo) then
            iio=i
            teeceo=teeceb(i)
          endif
        enddo
        beceo=bbf(iio)
!EALW           write(*,*) 'find beceo, iio,bbf(iio),teeceo'
!EALW           write(*,*) iio,bbf(iio),teeceo
!EALW        write(*,*)'beceo'
!EALW        write(*,*)beceo
!--------------------------------------------------------------------
!--     find becein(idesto), it close to beceo                     --
!--       dTe on beceo from tebit(idesto)                          --
!--------------------------------------------------------------------
        desto=abs(beceo-becein(1))
        idesto=1
        do i=2,necein
          if (abs(beceo-becein(i)).lt.desto) then
            desto=abs(beceo-becein(i))
            idesto=i
          endif
        enddo
!EALW        write(*,*)'idesto'
!EALW        write(*,*)idesto
!--------------------------------------------------------------------
!--     get bobit=dB=sqrt(dTe/Te'')                                --
!--     Te''-- (d2Te/dB2) at beceo--ppteppbo, dTe=tebit(idesto)    --
!--------------------------------------------------------------------
        bbx1=(bbf(iio+1)-b00)/baa
        bbx2=(bbf(iio-1)-b00)/baa
        ptpr1=0.
        ptpr2=0.
        do nk=2,nfit
          ptpr1=ptpr1+x(nk)*bbx1**(nk-2)
          ptpr2=ptpr2+x(nk)*bbx2**(nk-2)
        enddo
        ptpr1=ptpr1/baa
        ptpr2=ptpr2/baa
        ppteppbo=abs(0.5_dp*(ptpr1-ptpr2)/dbbf)
        dtero=abs(tebit(idesto))
        bobit=sqrt(dtero/ppteppbo)
!EALW        write(*,*)'bobit'
!EALW        write(*,*)bobit
      endif fixro
!---------------------------------------------------------------------
!--   take B+ (becep) from becein>beceo and get B- (becem)          --
!--   find nece (the number of B+)                                  --
!--       B+, B- from centre to edge                                --
!---------------------------------------------------------------------
      ECE: if ((kfitece.ne.1).and.(kfixrece.ne.1)) then
      ii=0
      do k=1,necein
        if ((beceo-becein(k)).lt.0.) then
          ii=ii+1
          becep(ii)=becein(k)
        endif
      enddo
      nece=ii
      do k=1,nece
        bbx=(becep(k)-b00)/baa
        teece(k)=0.
        do nk=1,nfit
          teece(k)=teece(k)+x(nk)*bbx**(nk-1)
        enddo
      enddo
!
      ii=0
      do i=1,nnnte
        if (bbf(i).lt.beceo) then
          ii=ii+1
          blowf(ii)=bbf(i)
          telowf(ii)=teeceb(i)
        endif
      enddo
!
      nlowf=ii
      call zpline(nlowf,telowf,blowf,bb,cc,dd)
      do k=1,nece
        becem(k)=seval(nlowf,teece(k),telowf,blowf,bb,cc,dd)
        if ((becem(k).ge.becein(1)).and.(becem(k).lt.beceo)) &
          cycle
        fwtece0(k)=0.0
        becem(k)=1.E-6_dp
      enddo
!--------------------------------------------------------------------
!--   idestm(nece)- the point becein(idestm) close to B-(nece)
!--------------------------------------------------------------------
      do k=1,nece
        dest=abs(becem(k)-becein(1))
        idestm(k)=1
        do i=2,necein
          if (abs(becem(k)-becein(i)).lt.dest) then
            dest=abs(becem(k)-becein(i))
            idestm(k)=i
          endif
        enddo
      enddo
!---------------------------------------------------------------------
!--   Calculation of |B| on rgrid (z=zeceo)   bfield(nw)            --
!---------------------------------------------------------------------
      endif ECE
      select case (icurrt)
      case (1)
        ffprim(1)=cratio*srma*2.*salpha/darea*twopi*tmu
        ffprim(nw)=ffprim(1)
      case (2,5)
        ffprim(nw)=fpcurr(x111,kffcur)/darea*twopi*tmu
        ffprim(1)=fpcurr(x000,kffcur)/darea*twopi*tmu
      case (4)
        call currnt(n222,jtime,n222,kerror)
        if (kerror.gt.0) return
        ffprim(1)=rbetap*cratio*rzero*twopi*tmu/darea
        ffprim(nw)=ffprim(1)*gammaf
      end select
      do i=2,nw-1
        siii=sigrid(i)
        select case (icurrt)
        case (1)
          ffprim(i)=ffprim(1)
        case (2,5)
          ffprim(i)=fpcurr(siii,kffcur)/darea*twopi*tmu
        case (4)
          ffprim(i)=ffprim(1)*(1.-siii**enp)**emp*(1.-gammap)+gammap
        end select
      enddo
      fpol(nw)=fbrdy*tmu
!EALW      write(*,*)'fpol(nw)'
!EALW      write(*,*)fpol(nw)
      sumf=fpol(nw)**2/2.
!EALW      write(*,*)'psibry'
!EALW      write(*,*)psibry
!EALW      write(*,*)'simag'
!EALW      write(*,*)simag
      delsi=-(psibry-simag)/(nw-1)
      do i=1,nw-1
        sumf=sumf+0.5_dp*delsi*(ffprim(nw-i+1)+ffprim(nw-i))
        if (sumf .ge. 0.0) then
          fpol(nw-i)=sqrt(2.*sumf)*fpol(nw)/abs(fpol(nw))
        else
          fpol(nw-i)=fpol(nw)
        endif
      enddo
!EALW      write(*,*)'fpol'
!EALW      write(*,*)fpol
!EALW      write(*,*)'sigrid'
!EALW      write(*,*)sigrid
      call zpline(nw,sigrid,fpol,bbb,ccc,ddd)
      do iw=1,nw
        kk=(iw-1)*nh+jo
        if (xpsi(kk).gt.1.0.or.ivacum.eq.1) then
          fnow=fbrdy*tmu
!EALW            write(*,*)'iw, fnow'
!EALW            write(*,*)iw,fnow
        else
          fnow=seval(nw,xpsi(kk),sigrid,fpol,bbb,ccc,ddd)
!EALW            write(*,*)'iw, xpsi(kk),fnow'
!EALW            write(*,*)iw, xpsi(kk),fnow
        endif
        btttt(iw)=fnow/rgrid(iw)
      enddo
      ier = 0
      call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
      do iw=1,nw
        rw=rgrid(iw)
        rh=zgrid(jo)
        kk=(iw-1)*nh+jo
        call seva2d(bkx,lkx,bky,lky,c,rw,rh,pds,ier,n555)
        if (ier.ne.0) then
          write (nttyo,9910) ier,rw,rh
!EALW            write(*,*) 'ier,rw,rh'
!EALW            write(*,*) ier,rw,rh
          return
        endif
        bfield(iw)=sqrt(pds(2)**2+pds(3)**2)/rgrid(iw)
        dsidr(iw)=pds(2)
        ddsiddr(iw)=pds(5)
        bfield(iw)=sqrt(bfield(iw)**2+btttt(iw)**2)
      enddo
!EALW       write(*,*)'bfield'
!EALW       write(*,*)bfield
 9910 format('  error in spline =',i4,' (r,z)= ( ',f5.2,',',f5.2,')')
!--------------------------------------------------------------------
!--   get d|B|/dr on rgrid
!--------------------------------------------------------------------
      do i=2,nw-1
        dbdr(i)=(bfield(i+1)-bfield(i-1))/(rgrid(i+1)-rgrid(i-1))
      enddo
      dbdr(1)=(bfield(2)-bfield(1))/(rgrid(2)-rgrid(1))
      dbdr(nw)=(bfield(nw)-bfield(nw-1))/(rgrid(nw)-rgrid(nw-1))
!--------------------------------------------------------------------
!--   find the number of bmink and bmaxk (local min and max)
!--------------------------------------------------------------------
      kmin=0
      kmax=0
      do k=2,nw-1
        if ((bfield(k).lt.bfield(k+1)).and. &
            (bfield(k).lt.bfield(k-1))) then
          kmin=kmin+1
          bmink(kmin)=bfield(k)
        endif
        if ((bfield(k).gt.bfield(k+1)).and. &
            (bfield(k).gt.bfield(k-1))) then
          kmax=kmax+1
          bmaxk(kmax)=bfield(k)
        endif
      enddo
      if (kmin.eq.(kmax+1)) then
        kmax=kmax+1
        bmaxk(kmax)=bfield(nw)
      endif
      if (kmin.eq.(kmax-1)) then
        kmin=kmin+1
        bmink(1)=bfield(1)
        do k=2,kmin
          bbk(k)=bmink(k-1)
        enddo
        do k=2,kmin
          bmink(k)=bbk(k)
        enddo
      endif
!
      k_max: if (kmax.ne.0) then
!
!--------------------------------------------------------------------
!--   get babsk array (kmax+1) is strictly increasing order
!--------------------------------------------------------------------
      do k=1,kmax
        if (bmaxk(k).lt.bmink(k)) then
          kerror = 1
          call errctrl_msg('geteceb','bmax < bmin')
          return
        endif
      enddo
!
      iout1=0
      do i=1,nw
        if (bfield(i).gt.bmaxk(1)) then
          iout1=iout1+1
          bout(1,iout1)=bfield(i)
          rrout(1,iout1)=rgrid(i)
        endif
      enddo
      nnout(1)=iout1
      do k=2,kmax
        ioutk=0
        do i=1,nw
          if ((bfield(i).gt.bmaxk(k)).and. &
              (bfield(i).lt.bmink(k-1))) then
            ioutk=ioutk+1
            bout(k,ioutk)=bfield(i)
            rrout(k,ioutk)=rgrid(i)
          endif
        enddo
        nnout(k)=ioutk
      enddo
      ioutk1=0
      do i=1,nw
        if (bfield(i).lt.bmink(kmax)) then
          ioutk1=ioutk1+1
          bout(kmax+1,ioutk1)=bfield(i)
          rrout(kmax+1,ioutk1)=rgrid(i)
        endif    
      enddo
      nnout(kmax+1)=ioutk1
!
      do k=1,kmax+1
        if (nnout(k).gt.3) then
          n=nnout(k)
          do i=1,n
            babs(k,i)=bout(k,n-i+1)
            rrgrid(k,i)=rrout(k,n-i+1)
          enddo
        endif
      enddo
!EALW       write(*,*)'kmax,kmin'
!EALW       write(*,*)kmax,kmin
!EALW       write(*,*)'bmaxk,bmink'
!EALW       write(*,*)bmaxk,bmink
!EALW       write(*,*)'nnout'
!EALW       write(*,*)nnout
!EALW       write(*,*) 'babs, rrgrid'
      do k=1,kmax+1
        n=nnout(k)
!EALW       do i=1,n
!EALW         write(*,*) babs(k,i)
!EALW         write(*,*)rrgrid(k,i)
!EALW       enddo
      enddo
!-------------------------------------------------------------------
!--   get R-,R+,Ro  where |B| = B+,B-,Bo                          --
!-------------------------------------------------------------------
      do m=1,nece
        recem(m)=1.E-6_dp
        recep(m)=1.E-6_dp
      enddo
      receo=1.e-6_dp
!
      if (nnout(1).gt.3) then
        n=nnout(1)
        do i=1,n
          bx(i)=babs(1,i)
          ry(i)=rrgrid(1,i)
        enddo
        call zpline(n,bx,ry,bbb,ccc,ddd)
!
        if(beceo.ge.bmaxk(1)) &
          receo=seval(n,beceo,bx,ry,bbb,ccc,ddd)
        do m=1,nece
          if(becep(m).ge.bmaxk(1)) &
            recem(m)=seval(n,becep(m),bx,ry,bbb,ccc,ddd)
          if(becem(m).ge.bmaxk(1)) &
            recep(m)=seval(n,becem(m),bx,ry,bbb,ccc,ddd)
        enddo
      endif
!
      do k=2,kmax
        if (nnout(k).gt.3) then
          n=nnout(k)
          do i=1,n
            bx(i)=babs(k,i)
            ry(i)=rrgrid(k,i)
          enddo
          call zpline(n,bx,ry,bbb,ccc,ddd)
          if((beceo.ge.bmaxk(k)).and.(beceo.le.bmink(k-1))) &
            receo=seval(n,beceo,bx,ry,bbb,ccc,ddd)
          do m=1,nece
            if((becep(m).ge.bmaxk(k)).and.(becep(m).le.bmink(k-1))) &
              recem(m)=seval(n,becep(m),bx,ry,bbb,ccc,ddd)
            if((becem(m).ge.bmaxk(k)).and.(becem(m).le.bmink(k-1))) &
              recep(m)=seval(n,becem(m),bx,ry,bbb,ccc,ddd)
          enddo
        endif
      enddo
!
      if (nnout(kmax+1).gt.3) then
        n=nnout(kmax+1)
        do i=1,n
          bx(i)=babs(kmax+1,i)
          ry(i)=rrgrid(kmax+1,i)
        enddo
        call zpline(n,bx,ry,bbb,ccc,ddd)
        if(beceo.le.bmink(kmax)) &
           receo=seval(n,beceo,bx,ry,bbb,ccc,ddd)
        do m=1,nece
          if(becep(m).le.bmink(kmax)) &
            recem(m)=seval(n,becep(m),bx,ry,bbb,ccc,ddd)
          if(becem(m).le.bmink(kmax)) &
            recep(m)=seval(n,becem(m),bx,ry,bbb,ccc,ddd)
        enddo
      endif
!
      do m=1,nece
        if ((recep(m).lt.1.E-5_dp).or.(recem(m).lt.1.E-5_dp)) then
          fwtece0(m)=0.0
          recep(m)=1.0E-6_dp
          recem(m)=1.0E-6_dp
        endif
      enddo
      if (receo.lt.1.E-5_dp) fwtecebz0=0.0
!
      else k_max
      do i=1,nw
        bx(i)=bfield(nw-i+1)
        ry(i)=rgrid(nw-i+1)
      enddo
      call zpline(nw,bx,ry,bbb,ccc,ddd)
      do m=1,nece
        recep(m)=seval(nw,becem(m),bx,ry,bbb,ccc,ddd)
        recem(m)=seval(nw,becep(m),bx,ry,bbb,ccc,ddd)
      enddo
      receo=seval(nw,beceo,bx,ry,bbb,ccc,ddd)
      endif k_max
!EALW      write(*,*)'recem'
!EALW      write(*,*)recem
!EALW      write(*,*)'recep'
!EALW      write(*,*)recep
!EALW      write(*,*)'receo'
!EALW      write(*,*)receo
!EALW      write(*,*)'nece'
!EALW      write(*,*)nece
!------------------------------------------------------------------
!--   get dB/dr at receo (dbdro) and recep,recem (dbdrp,dbdrm)
!------------------------------------------------------------------
      call zpline(nw,rgrid,dbdr,bbb,ccc,ddd)
      if(fwtecebz0.gt.1.e-6_dp) &
        dbdro=seval(nw,receo,rgrid,dbdr,bbb,ccc,ddd)
      do k=1,nece
        if (fwtece0(k).gt.1.e-6_dp) then
          dbdrp(k)=seval(nw,recep(k),rgrid,dbdr,bbb,ccc,ddd)
          dbdrm(k)=seval(nw,recem(k),rgrid,dbdr,bbb,ccc,ddd)
        endif
      enddo  
!--------------------------------------------------------------------
!--   get robit-- from bobit and dB/dr, (robit=dB/(dB/dr),bobit--dB)
!--------------------------------------------------------------------
      robit=bobit/dbdro
!--------------------------------------------------------------------
!--   get ecebit from sqrt(bitm**2+bitp**2) ,
!--         bit(m,p)=dTe *(dpsi/dR)/(dTe/dB)/(dB/dR)
!--         dTe(m)=tebit(idestm),  dTe(p)=tebit(idestp)
!--         (dTe/dB)*(dB/dR)--pteprp,pteprm
!--         (dpsi/dR)--dsidrp,dsidrm
!---------------------------------------------------------------------
      do k=1,nece
        if (fwtece0(k).gt.1.e-6_dp) then
          bbxp=(becep(k)-b00)/baa
          bbxm=(becem(k)-b00)/baa
          pteprp(k)=0.
          pteprm(k)=0.
          do nk=2,nfit
            pteprp(k)=pteprp(k)+x(nk)*bbx**(nk-2)
            pteprm(k)=pteprm(k)+x(nk)*bbx**(nk-2)
          enddo
          pteprp(k)=pteprp(k)/baa
          pteprm(k)=pteprm(k)/baa
          pteprp(k)=pteprp(k)*dbdrp(k)
          pteprm(k)=pteprm(k)*dbdrm(k)
        endif
      enddo
!
      call zpline(nw,rgrid,dsidr,bbb,ccc,ddd)
      do k=1,nece
        if (fwtece0(k).gt.1.e-6_dp) then
          dsidrm=seval(nw,recem(k),rgrid,dsidr,bbb,ccc,ddd)
          dsidrp=seval(nw,recep(k),rgrid,dsidr,bbb,ccc,ddd)
          if ((abs(pteprm(k)).gt.1.E-10_dp).and.(abs(pteprp(k)).gt.1.E-10_dp)) &
            then
            imk=idestm(k)
            rmbit(k)=tebit(imk)/pteprm(k)
            bitm=rmbit(k)*dsidrm
            ipk=idestp(k)
            rpbit(k)=tebit(ipk)/pteprp(k)
            bitp=rpbit(k)*dsidrp
            ecebit(k)=sqrt(bitm**2+bitp**2)
          else
            fwtece0(k)=0.
          endif
        endif
      enddo
!EALW       write(*,*)'rmbit'
!EALW       write(*,*)rmbit
!EALW       write(*,*)'rpbit'
!EALW       write(*,*)rpbit
      endif getECE
!
      return
      end subroutine geteceb

!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          getecer obtains the receo, R+ R-                        **
!**          from ECE measurement data T(R)                          **
!**          if kfixro  kfixrece = 0 called in setece                **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          11/98..........first created, Cheng Zhang               **
!**     2013/08/07..........Update for real-space Ti                 **
!**                                                                  **
!**********************************************************************
      subroutine getecer(jtime,kerror)
      use commonblocks,only: c,wk,bkx,bky
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: jtime
      integer*4, intent(out) :: kerror
      real*8 fpcurr,seval
      integer*4, parameter :: nnn1=1,kbre=5
      integer*4 i,idesto,ier,ii,iieerr,iio,imk,ipk,iout1,ioutk,ioutk1,iw, &
                j,k,kk,kmin,kmax,kout,m,mnow,n,nj,nk,nlowf
      real*8 bitm,bitp,delsi,dest,desto,drrr,dsidrm,dsidrp,dtero, &
             fnow,fwtcm,ppteppro,ptpr1,ptpr2, &
             raa,rh,rmax,rmin,rw,rx,rx1,rx2,r00,siii,sumf,t,teeceo,toler
      real*8 pds(6),nnout(kbre),bmink(kbre),bmaxk(kbre)
      real*8 arspfit(nnecein,nfit),brspfit(nnecein), &
             s(nfit),tte(nnecein),x(nfit),an(nnecein), &
             tebit(nnecein)
      real*8 telowf(nnnte),rlowf(nnnte),bb(nnnte),cc(nnnte),dd(nnnte), &
             teece(nnece),pteprm(nnece),pteprp(nnece), &
             idestp(nnece),idestm(nnece)
      real*8 rrgrid(kbre,nw),bfield(nw),rrout(kbre,nw), &
             bout(kbre,nw),babs(kbre,nw),bbb(nw),ccc(nw), &
             ddd(nw),btttt(nw),dsidr(nw),ddsiddr(nw),bx(nw), &
             ry(nw),bbk(nw)
!
      kerror = 0
!
      fwtece0=swtece 
      fwtecebz0=swtecebz
      babs=0.0
      bout=0.0
      rrout=0.0
      rrgrid=0.0
!---------------------------------------------------------------------
!--   Calculation of |B| array from fe array (harmonic nharm)       --
!--     becein(necein), fe(GHz), |B|(T), becein from H.f to L.f     --
!---------------------------------------------------------------------
      becein(1:necein)=0.001_dp*6.0*9.1095_dp*pi/4.8032_dp &
                      *feece(1:necein)/nharm
!---------------------------------------------------------------------
!--   Calculation of |B| on rgrid (z=zeceo)   bfield(nw)            --
!---------------------------------------------------------------------
      select case (icurrt)
      case (1)
        ffprim(1)=cratio*srma*2.*salpha/darea*twopi*tmu
        ffprim(nw)=ffprim(1)
      case (2,5)
        ffprim(nw)=fpcurr(x111,kffcur)/darea*twopi*tmu
        ffprim(1)=fpcurr(x000,kffcur)/darea*twopi*tmu
      case (4)
        call currnt(n222,jtime,n222,kerror)
        if(kerror.gt.0) return
        ffprim(1)=rbetap*cratio*rzero*twopi*tmu/darea
        ffprim(nw)=ffprim(1)*gammaf
      end select
      do i=2,nw-1
        siii=sigrid(i)
        select case (icurrt)
        case (1)
          ffprim(i)=ffprim(1)
        case (2,5)
          ffprim(i)=fpcurr(siii,kffcur)/darea*twopi*tmu
        case (4)
          ffprim(i)=ffprim(1)*(1.-siii**enp)**emp*(1.-gammap)+gammap
        end select
      enddo
      fpol(nw)=fbrdy*tmu
      sumf=fpol(nw)**2/2.
      delsi=-(psibry-simag)/(nw-1)
      do i=1,nw-1
        sumf=sumf+0.5_dp*delsi*(ffprim(nw-i+1)+ffprim(nw-i))
        if (sumf .ge. 0.0) then
          fpol(nw-i)=sqrt(2.*sumf)*fpol(nw)/abs(fpol(nw))
        else
          fpol(nw-i)=fpol(nw)
        endif
      enddo
      call zpline(nw,sigrid,fpol,bbb,ccc,ddd)
      do iw=1,nw
        kk=(iw-1)*nh+jo
        if (xpsi(kk).gt.1.0.or.ivacum.eq.1) then
          fnow=fbrdy*tmu
        else
          fnow=seval(nw,xpsi(kk),sigrid,fpol,bbb,ccc,ddd)
        endif
        btttt(iw)=fnow/rgrid(iw)
      enddo
!
      ier = 0
      call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
      do iw=1,nw
        rw=rgrid(iw)
        rh=zgrid(jo)
        kk=(iw-1)*nh+jo
        call seva2d(bkx,lkx,bky,lky,c,rw,rh,pds,ier,n555)
        if (ier.ne.0) then
          write (nttyo,9910) ier,rw,rh
          return
        endif
        bfield(iw)=sqrt(pds(2)**2+pds(3)**2)/rgrid(iw)
        dsidr(iw)=pds(2)
        ddsiddr(iw)=pds(5)            
        bfield(iw)=sqrt(bfield(iw)**2+btttt(iw)**2)
      enddo
 9910 format('  error in spline =',i4,' (r,z)= ( ',f5.2,',',f5.2,')')
!--------------------------------------------------------------------
!--   find the number of bmink and bmaxk (local min and max)
!--------------------------------------------------------------------
      kmin=0
      kmax=0
      do k=2,nw-1
        if ((bfield(k).lt.bfield(k+1)).and. &
            (bfield(k).lt.bfield(k-1))) then
          kmin=kmin+1
          bmink(kmin)=bfield(k)
        endif
        if ((bfield(k).gt.bfield(k+1)).and. &
            (bfield(k).gt.bfield(k-1))) then
          kmax=kmax+1
          bmaxk(kmax)=bfield(k)
        endif
      enddo
      if (kmin.eq.(kmax+1)) then
        kmax=kmax+1  
        bmaxk(kmax)=bfield(nw)
      endif
      if (kmin.eq.(kmax-1)) then
        kmin=kmin+1
        bmink(1)=bfield(1)
        do k=2,kmin
          bbk(k)=bmink(k-1)
        enddo
        do k=2,kmin
          bmink(k)=bbk(k)
        enddo
      endif 
!
#ifdef DEBUG_LEVEL3
      write (6,*) 'GETECER, kmin/kmax = ', kmin, kmax
#endif
      kmax0: if (kmax.ne.0) then
!--------------------------------------------------------------------
!--   get babsk array (kmax+1) is strictly increasing order
!-------------------------------------------------------------------- 
      do k=1,kmax
        if(bmaxk(k).lt.bmink(k)) stop
      enddo
!
      iout1=0
      do i=1,nw
        if (bfield(i).gt.bmaxk(1)) then
          iout1=iout1+1
          bout(1,iout1)=bfield(i)
          rrout(1,iout1)=rgrid(i)
        endif
      enddo
      nnout(1)=iout1
      do k=2,kmax
        ioutk=0
        do i=1,nw
          if ((bfield(i).gt.bmaxk(k)).and. &
              (bfield(i).lt.bmink(k-1))) then
            ioutk=ioutk+1
            bout(k,ioutk)=bfield(i)
            rrout(k,ioutk)=rgrid(i)
          endif
        enddo
        nnout(k)=ioutk
      enddo
      ioutk1=0
      do i=1,nw
        if (bfield(i).lt.bmink(kmax)) then
          ioutk1=ioutk1+1
          bout(kmax+1,ioutk1)=bfield(i)
          rrout(kmax+1,ioutk1)=rgrid(i)
        endif       
      enddo
      nnout(kmax+1)=ioutk1
!
      do k=1,kmax+1
        if (nnout(k).gt.3) then
          n=nnout(k)
          do i=1,n
            babs(k,i)=bout(k,n-i+1)
            rrgrid(k,i)=rrout(k,n-i+1)
          enddo
        endif
      enddo 
!-------------------------------------------------------------------
!--   get recein, at which |B| = becein                           --
!--     recein(mecein)                                            --
!-------------------------------------------------------------------
      kout=0
!
      if (nnout(1).gt.3) then
      n=nnout(1)
      do i=1,n
        bx(i)=babs(1,i)
        ry(i)=rrgrid(1,i)
      enddo
      call zpline(n,bx,ry,bbb,ccc,ddd)
      do m=1,necein
        if (becein(m).ge.bmaxk(1)) then
          kout=kout+1
          recein(kout)=seval(n,becein(m),bx,ry,bbb,ccc,ddd) 
          teeceinr(kout)=teecein(m)
          tebit(kout)=errorece(m)         
        endif
      enddo
      endif
!
      do k=2,kmax
        if (nnout(k).gt.3) then
          n=nnout(k)
          do i=1,n
            bx(i)=babs(k,i)
            ry(i)=rrgrid(k,i)
          enddo
          call zpline(n,bx,ry,bbb,ccc,ddd)
          do m=1,necein
            if ((becein(m).ge.bmaxk(k)).and.(becein(m).le.bmink(k-1))) then
              kout=kout+1
              recein(kout)=seval(n,becein(m),bx,ry,bbb,ccc,ddd)
              teeceinr(kout)=teecein(m)
              tebit(kout)=errorece(m)
            endif
          enddo
        endif
      enddo
!
      if (nnout(kmax+1).gt.3) then
        n=nnout(kmax+1)
        do i=1,n
          bx(i)=babs(kmax+1,i)
          ry(i)=rrgrid(kmax+1,i)
        enddo
        call zpline(n,bx,ry,bbb,ccc,ddd)
        do m=1,necein
          if (becein(m).le.bmink(kmax)) then
            kout=kout+1
            recein(kout)=seval(n,becein(m),bx,ry,bbb,ccc,ddd)
            teeceinr(kout)=teecein(m)
            tebit(kout)=errorece(m)
          endif
        enddo
      endif
! 
      mecein=kout
!
!
      else kmax0
      do i=1,nw
        bx(i)=bfield(nw-i+1)
        ry(i)=rgrid(nw-i+1)
      enddo
      call zpline(nw,bx,ry,bbb,ccc,ddd)
      do m=1,necein
        recein(m)=seval(nw,becein(m),bx,ry,bbb,ccc,ddd)
      enddo
      teeceinr(1:necein)=teecein(1:necein)
      tebit(1:necein)=errorece(1:necein)
      mecein=necein
!
      endif kmax0
!--------------------------------------------------------------------
!--   fitting data from teeceinr,tebit and recein (nnecein)        --
!--     rx=(R-r00)/raa                                             --
!--     Te=x(1)+x(2)*rx+x(3)*rx**2+...+x(nfit)*rx**(nfit-1)        --
!--------------------------------------------------------------------
      rmin=recein(1)
      rmax=recein(mecein)
      raa=0.5_dp*(rmax-rmin)
      r00=0.5_dp*(rmax+rmin)
      do i=1,mecein
        an(i)=(recein(i)-r00)/raa
        tebit(i)=max(tebit(i),1.e-4_dp)
      enddo
      do nj=1,mecein
        do nk=1,nfit
          if (nk.eq.1) then
            arspfit(nj,nk)=1./tebit(nj)
          else
            arspfit(nj,nk)=an(nj)**(nk-1)/tebit(nj)
          endif
        enddo
      enddo
      brspfit(1:mecein)=teeceinr(1:mecein)/tebit(1:mecein)
!
      mnow=mecein
      if (kcmin.gt.0) then
        fwtnow=0.001_dp ! input option, why hardcoded here?
        fwtcm=1.0
        do j=1,nfit
          mnow=mnow+1
          do k=1,nfit
            if (j.ne.k) then
              arspfit(mecein+j,k)=0.0
            else
              arspfit(mecein+j,k)=fwtcm/fwtnow
            endif
          enddo
          brspfit(mecein+j)=0.0
        enddo
      endif
!
      iieerr=0
      call sdecm(arspfit,nnecein,mnow,nfit,brspfit,nnecein,nnn1,s,wk,iieerr)
      if (iieerr.eq.129) then
        kerror = 1
        call errctrl_msg('getecer','sdecm failed to converge')
        return
      end if
      toler=1.0e-06_dp*s(1)
      do i=1,nfit
        t=0.0
        if(s(i).gt.toler) t=brspfit(i)/s(i)
        brspfit(i)=t
      enddo
      do i=1,nfit
        x(i)=sum(arspfit(i,1:nfit)*brspfit(1:nfit))
      enddo
      xfit(1:nfit)=x(1:nfit)
      chisqfit=0.0
      do k=1,mecein
        tte(k)=0.
        do nk=1,nfit
          tte(k)=tte(k)+x(nk)*an(k)**(nk-1)
        enddo
        chisqfit=chisqfit+(tte(k)-teeceinr(k))**2/tebit(k)
      enddo
!--------------------------------------------------------------------
!--   get Teecer(rrr) in ECE data region                           --
!--------------------------------------------------------------------
      drrr=(recein(mecein)-recein(1))/(nnnte-1)
      do i=1,nnnte
        rrr(i)=recein(1)+drrr*(i-1)
        rx=(rrr(i)-r00)/raa
        teecer(i)=0.
        do nk=1,nfit
          teecer(i)=teecer(i)+x(nk)*rx**(nk-1)
        enddo
      enddo 
!---------------------------------------------------------------------
!--   find receo which is Te peak point                            --
!---------------------------------------------------------------------
      fixro: if (kfixro.ne.1) then      
        teeceo=teecer(1)
        iio=1
        do i=2,nnnte
          if (teecer(i).gt.teeceo) then
            iio=i
            teeceo=teecer(i)
          endif
        enddo
        receo=rrr(iio)
!--------------------------------------------------------------------
!--     find recein(idesto), it close to receo                     --
!--       dTe on receo from tebit(idesto)                          --
!--------------------------------------------------------------------
        desto=abs(receo-recein(1))
        idesto=1
        do i=2,mecein
          if (abs(receo-recein(i)).lt.desto) then
            desto=abs(receo-recein(i))
            idesto=i
          endif
        enddo
!--------------------------------------------------------------------
!--     get robit,  robit=sqrt(dTe/Te'')                           --
!--     Te''-- (d2Te/dR2) at receo--ppteppro, dTe=tebit(idesto)    --
!--------------------------------------------------------------------
        rx1=(rrr(iio+1)-r00)/raa
        rx2=(rrr(iio-1)-r00)/raa
        ptpr1=0.
        ptpr2=0.
        do nk=2,nfit
          ptpr1=ptpr1+x(nk)*rx1**(nk-2)
          ptpr2=ptpr2+x(nk)*rx2**(nk-2)
        enddo
        ptpr1=ptpr1/raa
        ptpr2=ptpr2/raa
        ppteppro=abs(0.5_dp*(ptpr1-ptpr2)/drrr)
        dtero=abs(tebit(idesto))
        robit=sqrt(dtero/ppteppro)
      endif fixro
!---------------------------------------------------------------------
!--   take R- and get R+                                            --
!--       nece=the number of R-,  recem(nece), recep(nece)          --
!---------------------------------------------------------------------
      ECE: if ((kfitece.ne.1).and.(kfixrece.ne.1)) then
        ii=0
        do k=mecein,1,-1
          if ((receo-recein(k)).gt.0.) then
            ii=ii+1
            recem(ii)=recein(k)
          endif
        enddo
        nece=ii
        do k=1,nece
          rx=(recem(k)-r00)/raa
          teece(k)=0.
          pteprm(k)=0.
          do nk=1,nfit
            teece(k)=teece(k)+x(nk)*rx**(nk-1)
          enddo
          do nk=2,nfit
            pteprm(k)=pteprm(k)+x(nk)*rx**(nk-2)
          enddo         
          pteprm(k)=pteprm(k)/raa
        enddo
!
        ii=0
        do i=nnnte,1,-1
          if (rrr(i).gt.receo) then
            ii=ii+1
            rlowf(ii)=rrr(i)
            telowf(ii)=teecer(i)
          endif
        enddo
        nlowf=ii
        call zpline(nlowf,telowf,rlowf,bb,cc,dd)
        do k=1,nece
          recep(k)=seval(nlowf,teece(k),telowf,rlowf,bb,cc,dd)
          if((recep(k).gt.receo).and.(recep(k).lt.recein(mecein))) &
            cycle
          fwtece0(k)=0.0
        enddo
!--------------------------------------------------------------------
!--     idestp(nece)- the point recein(idestp) close to R+(nece)
!--     idestm(nece)- the point recein(idestm) close to R-(nece)
!---------------------------------------------------------------------
        do k=1,nece
          dest=abs(recep(k)-recein(1))
          idestp(k)=1
          do i=2,mecein
            if (abs(recep(k)-recein(i)).lt.dest) then
              dest=abs(recep(k)-recein(i))
              idestp(k)=i
            endif
          enddo
        enddo 
! 
        do k=1,nece
          dest=abs(recem(k)-recein(1))
          idestm(k)=1
          do i=2,mecein
            if (abs(recem(k)-recein(i)).lt.dest) then
              dest=abs(recem(k)-recein(i))
              idestm(k)=i
            endif
          enddo
        enddo 
!--------------------------------------------------------------------
!--     get ecebit from sqrt(bitm**2+bitp**2) ,
!--         bit(m,p)=dTe *(dpsi/dR)/(dTe/dR)
!--         dTe(m)=tebit(idestm),  dTe(p)=tebit(idestp)
!---------------------------------------------------------------------
        do k=1,nece
          rx=(recep(k)-r00)/raa
          pteprp(k)=0.
          do nk=2,nfit
            pteprp(k)=pteprp(k)+x(nk)*rx**(nk-2)
          enddo
          pteprp(k)=pteprp(k)/raa
        enddo
!
        call zpline(nw,rgrid,dsidr,bbb,ccc,ddd)
        do k=1,nece
          dsidrm=seval(nw,recem(k),rgrid,dsidr,bbb,ccc,ddd)
          dsidrp=seval(nw,recep(k),rgrid,dsidr,bbb,ccc,ddd)
          if((abs(pteprm(k)).gt.1.E-10_dp).and.(abs(pteprp(k)).gt.1.E-10_dp)) &
            then
            imk=idestm(k)
            rmbit(k)=tebit(imk)/pteprm(k)
            bitm=rmbit(k)*dsidrm
            ipk=idestp(k)
            rpbit(k)=tebit(ipk)/pteprp(k)  
            bitp=rpbit(k)*dsidrp
            ecebit(k)=sqrt(bitm**2+bitp**2) 
          else
            fwtece0(k)=0.
          endif
        enddo
      endif ECE
!
      return
      end subroutine getecer

!**********************************************************************
!>
!!    gettir obtains the receo, R+ R-
!!    from Ti data
!!    kfitece = 3, called from setece
!!    CALLING ARGUMENTS:
!!    
!!    RECORD OF MODIFICATION:
!!    2013/08/07..........Update for real-space Ti based on ECE/Te
!!    
!!
!!    @param jtime : time index
!!
!!    @param kerror : error flag
!!
!**********************************************************************
      subroutine gettir(jtime,kerror)
      use commonblocks,only: c,wk,bkx,bky
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: jtime
      integer*4, intent(out) :: kerror
      real*8 seval
      integer*4, parameter :: nnn1=1,kbre=5
      integer*4 i,idesto,ier,ii,iieerr,iio,imk,ipk,iout1,ioutk,ioutk1,iw, &
                j,k,kk,kmin,kmax,kout,m,mmmte,mnow,n,nj,nk,nlowf
      real*8 bitm,bitp,delsi,dest,desto,drrr,dsidrm,dsidrp,dtero, &
             fnow,fwtcm,ppteppro,ptpr1,ptpr2, &
             raa,rh,rmax,rmin,rw,rx,rx1,rx2,r00,siii,sumf,t,teeceo,toler
      real*8 pds(6),nnout(kbre),bmink(kbre),bmaxk(kbre)
      real*8 arspfit(nnecein,nfit),brspfit(nnecein), &
             s(nfit),tte(nnecein),x(nfit),an(nnecein), &
             tebit(nnecein)
      real*8 telowf(nnnte),rlowf(nnnte),bb(nnnte),cc(nnnte),dd(nnnte), &
             teece(nnece),pteprm(nnece),pteprp(nnece), &
             idestp(nnece),idestm(nnece)
      real*8 rrgrid(kbre,nw),bfield(nw),rrout(kbre,nw), &
             bout(kbre,nw),babs(kbre,nw),bbb(nw),ccc(nw), &
             ddd(nw),btttt(nw),dsidr(nw),ddsiddr(nw),bx(nw), &
             ry(nw),bbk(nw)
!
#ifdef DEBUG_LEVEL3
      write (6,*) 'Enter GETTIR, kfitece/kfixrece = ',&
         kfitece, kfixrece
#endif
      kerror = 0
!
      fwtece0=swtece 
      fwtecebz0=swtecebz
      babs=0.0
      bout=0.0
      rrout=0.0
      rrgrid=0.0
!---------------------------------------------------------------------
!--   Copy Ti array                                                 --
!---------------------------------------------------------------------
      becein(1:necein)=feece(1:necein)
      recein(1:necein)=becein(1:necein)
      teeceinr(1:necein)=teecein(1:necein)
      tebit(1:necein)=errorece(1:necein)
!---------------------------------------------------------------------
!--   Calculation of |B| on rgrid (z=zeceo)   bfield(nw)            --
!---------------------------------------------------------------------
      do iw=1,nw
        rw=rgrid(iw)
        rh=zteo
        kk=(iw-1)*nh+jo
        call seva2d(bkx,lkx,bky,lky,c,rw,rh,pds,ier,n555)
        if (ier.ne.0) then
          write (nttyo,9910) ier,rw,rh
          return
        endif
        dsidr(iw)=pds(2)
        ddsiddr(iw)=pds(5)            
      enddo
 9910 format('  error in getecer/spline = ',i4,' (r,z)= ( ', &
                f5.2,',',f5.2,')')
!--------------------------------------------------------------------
!--   fitting data from teeceinr,tebit and recein (nnecein)        --
!--     rx=(R-r00)/raa                                             --
!--     Te=x(1)+x(2)*rx+x(3)*rx**2+...+x(nfit)*rx**(nfit-1)        --
!--------------------------------------------------------------------
      mecein=necein
      rmin=recein(1)
      rmax=recein(mecein)
      raa=0.5_dp*(rmax-rmin)
      r00=0.5_dp*(rmax+rmin)
      do i=1,mecein
        an(i)=(recein(i)-r00)/raa
        tebit(i)=max(tebit(i),1.e-4_dp)
      enddo
      do nj=1,mecein
        do nk=1,nfit
          if (nk.eq.1) then
            arspfit(nj,nk)=1./tebit(nj)
          else
            arspfit(nj,nk)=an(nj)**(nk-1)/tebit(nj)
          endif
        enddo
      enddo
      do nj=1,mecein
        brspfit(nj)=teeceinr(nj)/tebit(nj)
      enddo
!
#ifdef DEBUG_LEVEL3
      write (6,*) 'GETTIR/SDECM, mecein/raa/r00 = ',&
        mecein, raa, r00, nfit
      write (6,*) 'GETTIR teeceinr = ',(teeceinr(i),i=1,mecein)
      write (6,*) 'GETTIR tebit = ',(tebit(i),i=1,mecein)
#endif
      iieerr=0
      mnow=mecein
      call sdecm(arspfit,nnecein,mnow,nfit,brspfit,nnecein,nnn1,s,wk,iieerr)
      if (iieerr.eq.129) then
        kerror = 1
        call errctrl_msg('gettir','sdecm failed to converge')
        return
      end if
      toler=1.0e-06_dp*s(1)
      do i=1,nfit
        t=0.0
        if (s(i).gt.toler) t=brspfit(i)/s(i)
        brspfit(i)=t
      enddo
      do i=1,nfit
        x(i)=sum(arspfit(i,1:nfit)*brspfit(1:nfit))
      enddo
!
      xfit(1:nfit)=x(1:nfit)
      chisqfit=0.0
      do k=1,mecein
        tte(k)=0.
        do nk=1,nfit
          tte(k)=tte(k)+x(nk)*an(k)**(nk-1)
        enddo
        chisqfit=chisqfit+(tte(k)-teeceinr(k))**2/tebit(k)
      enddo
      mmmte = nnnte
#ifdef DEBUG_LEVEL3
      write (6,*) 'GETTIR chisqfit/kfixro/mnow = ', &
        chisqfit, kfixro, mnow
      write (6,*) 'GETTIR tte = ',(tte(i),i=1,mecein)
#endif
!--------------------------------------------------------------------
!--   get Teecer(rrr) in ECE data region                           --
!--------------------------------------------------------------------
      drrr=(recein(mecein)-recein(1))/(nnnte-1)
      do i=1,nnnte
        rrr(i)=recein(1)+drrr*(i-1)
        rx=(rrr(i)-r00)/raa
        teecer(i)=0.
        do nk=1,nfit
          teecer(i)=teecer(i)+x(nk)*rx**(nk-1)
        enddo
      enddo 
!---------------------------------------------------------------------
!--   find receo which is Te peak point                            --
!---------------------------------------------------------------------
      fixro: if (kfixro.ne.1) then
        teeceo=teecer(1)
        iio=1
        do i=2,nnnte
          if (teecer(i).gt.teeceo) then
            iio=i
            teeceo=teecer(i)
          endif
        enddo
        receo=rrr(iio)
#ifdef DEBUG_LEVEL3
        write (6,*) 'GETTIR teece, receo, iio = ',teeceo, receo, iio
#endif
!--------------------------------------------------------------------
!--     find recein(idesto), it close to receo                     --
!--       dTe on receo from tebit(idesto)                          --
!--------------------------------------------------------------------
        desto=abs(receo-recein(1))
        idesto=1
        do i=2,mecein
          if (abs(receo-recein(i)).lt.desto) then
            desto=abs(receo-recein(i))
            idesto=i
          endif
        enddo
!--------------------------------------------------------------------
!--     get robit,  robit=sqrt(2dTe/Te'')                          --
!--     Te''-- (d2Te/dR2) at receo--ppteppro, dTe=tebit(idesto)    --
!--------------------------------------------------------------------
        rx1=(rrr(iio+1)-r00)/raa
        rx2=(rrr(iio-1)-r00)/raa
        ptpr1=0.
        ptpr2=0.
        do nk=2,nfit
          ptpr1=ptpr1+x(nk)*rx1**(nk-2)*(nk-1)
          ptpr2=ptpr2+x(nk)*rx2**(nk-2)*(nk-1)
        enddo
        ptpr1=ptpr1/raa
        ptpr2=ptpr2/raa
        ppteppro=abs(0.5_dp*(ptpr1-ptpr2)/drrr)
        dtero=abs(tebit(idesto))
        robit=sqrt(2*dtero/ppteppro)
      endif fixro
!---------------------------------------------------------------------
!--   take R- and get R+                                            --
!--       nece=the number of R-,  recem(nece), recep(nece)          --
!---------------------------------------------------------------------
#ifdef DEBUG_LEVEL3
      write (6,*) 'GETTIR R-, kfitece/kfixrece  = ',kfitece, kfixrece
#endif
      ECE: if ((kfitece.ne.1).and.(kfixrece.ne.1)) then
        ii=0
!       do k=mecein,1,-1
        do k=1,mecein
          if ((receo-recein(k)).gt.0.) then
            ii=ii+1
            recem(ii)=recein(k)
          endif
        enddo
        nece=ii
        do k=1,nece
          rx=(recem(k)-r00)/raa
          teece(k)=0.
          pteprm(k)=0.
          do nk=1,nfit
            teece(k)=teece(k)+x(nk)*rx**(nk-1)
          enddo
          do nk=2,nfit
            pteprm(k)=pteprm(k)+x(nk)*rx**(nk-2)*(nk-1)
          enddo         
          pteprm(k)=pteprm(k)/raa
        enddo
#ifdef DEBUG_LEVEL3
        write (6,*) 'GETTIR R-, nece = ',nece
        write (6,*) 'GETTIR R-, recem = ',(recem(i),i=1,nece)
        write (6,*) 'GETTIR R-, teece = ',(teece(i),i=1,nece)
        write (6,*) 'GETTIR R-, pteprm = ',(pteprm(i),i=1,nece)
#endif
!
        ii=0
        do i=nnnte,1,-1
!       do i=1,nnnte
          if (rrr(i).gt.receo) then
            ii=ii+1
            rlowf(ii)=rrr(i)
            telowf(ii)=teecer(i)
          endif
        enddo
        nlowf=ii
!
        call zpline(nlowf,telowf,rlowf,bb,cc,dd)
        do k=1,nece
          recep(k)=seval(nlowf,teece(k),telowf,rlowf,bb,cc,dd)
          if ((recep(k).gt.receo).and.(recep(k).lt.recein(mecein))) &
            cycle
          fwtece0(k)=0.0
        enddo
#ifdef DEBUG_LEVEL3
        write (6,*) 'GETTIR R+, recep = ',(recep(i),i=1,nece)
#endif
!--------------------------------------------------------------------
!--     idestp(nece)- the point recein(idestp) close to R+(nece)
!--     idestm(nece)- the point recein(idestm) close to R-(nece)
!---------------------------------------------------------------------
        do k=1,nece
          dest=abs(recep(k)-recein(1))
          idestp(k)=1
          do i=2,mecein
            if (abs(recep(k)-recein(i)).lt.dest) then
              dest=abs(recep(k)-recein(i))
              idestp(k)=i
            endif
          enddo
        enddo 
! 
        do k=1,nece
          dest=abs(recem(k)-recein(1))
          idestm(k)=1
          do i=2,mecein
            if (abs(recem(k)-recein(i)).lt.dest) then
              dest=abs(recem(k)-recein(i))
              idestm(k)=i
            endif
          enddo
        enddo 
!--------------------------------------------------------------------
!--     get ecebit from sqrt(bitm**2+bitp**2) ,
!--         bit(m,p)=dTe *(dpsi/dR)/(dTe/dR)
!--         dTe(m)=tebit(idestm),  dTe(p)=tebit(idestp)
!---------------------------------------------------------------------
        do k=1,nece
          rx=(recep(k)-r00)/raa
          pteprp(k)=0.
          do nk=2,nfit
            pteprp(k)=pteprp(k)+x(nk)*rx**(nk-2)*(nk-1)
          enddo
          pteprp(k)=pteprp(k)/raa
        enddo
!
        call zpline(nw,rgrid,dsidr,bbb,ccc,ddd)
        do k=1,nece
          dsidrm=seval(nw,recem(k),rgrid,dsidr,bbb,ccc,ddd)
          dsidrp=seval(nw,recep(k),rgrid,dsidr,bbb,ccc,ddd)
          if ((abs(pteprm(k)).gt.1.E-10_dp).and.(abs(pteprp(k)).gt.1.E-10_dp)) &
            then
            imk=idestm(k)
            rmbit(k)=tebit(imk)/pteprm(k)
            bitm=rmbit(k)*dsidrm
            ipk=idestp(k)
            rpbit(k)=tebit(ipk)/pteprp(k)  
            bitp=rpbit(k)*dsidrp
            ecebit(k)=sqrt(bitm**2+bitp**2) 
          else
            fwtece0(k)=0.
          endif
        enddo
      endif ECE
#ifdef DEBUG_LEVEL3
      write (6,*) 'GETTIR, ecebit = ',(ecebit(i),i=1,nece)
#endif

!
      return
      end subroutine gettir
