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
!********************************************************************** 
      subroutine setece(jtime,kerror)
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      dimension rrrecep(nnece),rrrecem(nnece),iwp(nnece),iwm(nnece)
      dimension zece(nnece),pds(6),rsplt(2500),zsplt(2500),csplt(2500)
      real*8,dimension(:),allocatable :: dsidr,ddsiddr
      character*40 filenmme
      data nset/20/,cdum/1.0/
      save nset
      integer, intent(inout) :: kerror
      kerror = 0
!
      if (idebug.ge.3) write (6,*) 'Enter SETECE, nw = ',nw
      ALLOCATE(dsidr(nw),ddsiddr(nw))
!
      if (idebug.ge.3) write (6,*) 'Enter SETECE, ksetece = ',ksetece
      kece=0
      kecebz=0
!
      if (ksetece.ne.0) go to 60000
!------------------------------------------------------------------
!--     zeceo from k-file, jo--zeceo on zgrid
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
      if (idebug.ge.3) write (6,*) 'SETECE, receo/zeteo = ',receo,zeceo
!-------------------------------------------------------------------
!--     kfixro=1, receo  from k-file
!------------------------------------------------------------------- 
      if (kfixro.eq.1) then
         receo=rteo
         if (kfitece.eq.1) go to 20
      endif
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
      !--   kfixro=0 or kfixrece=0, receo or R+ R- from getecer
      !-----------------------------------------------------------------
      if (kfitece.le.2) then
        if (kfixrece.eq.0) then
          call getecer(jtime,kerror)
          if (kerror.gt.0) return
        end if
        if ((kfixro.eq.-1).or.(kfixrece.eq.-1)) then
          call geteceb(jtime,kerror)
          if (kerror.gt.0) return
        end if
      endif
      if (kfitece.eq.3) then
        call gettir(jtime,kerror)
        if (kerror.gt.0) return
      endif
20    continue
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
        do k=1,nece
          zece(k)=zeceo
        enddo
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
!-- try read responce function from recefile
!--------------------------------------------------------------------
      do iname=1,nset
        if (iname.le.9) then
          write(filenmme,6530) iname
        else
          write(filenmme,6540) iname
        endif
        open(unit=nffile, status='old',form='unformatted', &
             file=filenmme,err=10000)
        go to 15000
10000   continue
        open(unit=nffile, status='old',form='unformatted', &
             file=table_dir(1:ltbdir)//filenmme,err=30)
15000   continue
        read (nffile,err=30) nnnece
        if (nnnece.ne.nece) then
          close(unit=nffile)
          go to 20000
        endif
        read (nffile,err=30)  rrrecem
        read (nffile,err=30)  rrrecep
        read (nffile,err=30)  rrreceo
        do ii=1,nece
         if (abs(rrrecem(ii)-recem(ii)).gt.1.e-4_dp) go to 20000
         if (abs(rrrecep(ii)-recep(ii)).gt.1.e-4_dp) go to 20000
        enddo
        if (abs(rrreceo-receo).gt.1.e-4_dp) go to 20000
        read (nffile,err=30) recebzfc
        read (nffile,err=30) gecebzpc
        read (nffile,err=30) recebzec
        read (nffile,err=30) recefc
        read (nffile,err=30) gecepc
        read (nffile,err=30) receec
        close(unit=nffile)
        go to 25000
      enddo
20000 continue
      go to 30
25000 return
!-------------------------------------------------------------------       
!-- compute the response function about Te peak point constraint    
!--  response due to F coils --- recebzfc(ncoil)
!-------------------------------------------------------------------
30    continue
      if (idebug.ge.3) write (6,*) 'SETECE, rf/zf = ',rf(1),zf(1)
      isplit=10
      itot=isplit*isplit
      fitot=itot 
      recebzfc = 0.0
      do k=1,mfcoil
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
        recebzfc(kkm)=recebzfc(kkm)+bzct/fitot*tmu
      enddo
!---------------------------------------------------------------------
!--   plasma response   gecebzpc(nwnh)                              --
!---------------------------------------------------------------------
      zzx=0.0
      drzg=(drgrid**2+dzgrid**2)*25.
      isplit=18
      if (receo.le.1.e-8_dp)  go to 90
      r=receo
      do ii=1,nw
        a=rgrid(ii)
        do jj=1,nh
          z=zeceo-zgrid(jj)
          kk=(ii-1)*nh+jj
          dist=(a-r)**2+z**2
          if (dist.lt.drzg) then
   68       call splitc(isplit,rsplt,zsplt,csplt, &
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
   90 continue
!---------------------------------------------------------------------
!-- E coils response   recebzec                                     --
!---------------------------------------------------------------------
      if (iecurr.le.0) go to 201 
      isplit=10
      itot=isplit*isplit
      fitot=real(itot,dp)
      do 170 i=1,nesum
         recebzec(i)=0.0
  170 continue
      if (receo.le.1.e-8_dp) go to 201
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
  201 continue
!-------------------------------------------------------------------------
!--  compute the response function about R+ R- constraint     
!--  F coils response -----recefc(nece,nfcoil)
!-------------------------------------------------------------------------
      if (kfitece.eq.1) go to 1710
      do 1500 n=1,nfcoil
        call sets2d(gridfc(1,n),c,rgrid,nw,bkx,lkx,zgrid,nh,bky, &
                       lky,wk,ier)
        do 1450 i=1,nece
        if (eceiter.eq.'pair') then
        call seva2d(bkx,lkx,bky,lky,c,recem(i),zece(i),pds,ier,n111)
        recefc(i,n)=pds(1)
        endif
        if (eceiter.eq.'flux') recefc(i,n)=0.0
 1450   continue
        do 1460 i=1,nece
         call seva2d(bkx,lkx,bky,lky,c,recep(i),zece(i),pds,ier,n111)
         recefc(i,n)=pds(1)-recefc(i,n)
 1460   continue
 1500 continue 
!------------------------------------------------------------------------------
!--    plasma response -----gecempc(nece,nwnh),geceppc(nece,nwnh)            --
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
          if (rsum.gt.dselsum) go to 72054
!          mk=(i-1)*nh+1
!          gecempc(m,k)=gridpc(mk,i)
           zdif=dselsum
           rselsum=rgrid(i)-dselsum
           gecempc(m,k)=psical(recem(m),rselsum,zdif)*tmu
          go to 72056
72054     continue
          gecempc(m,k)=psical(recem(m),rgrid(i),zdif)*tmu
72056    continue
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
          if (rsum.gt.dselsum) go to 73054
!          mk=(i-1)*nh+1
!          geceppc(m,k)=gridpc(mk,i)
           zdif=dselsum
           rselsum=rgrid(i)-dselsum
           geceppc(m,k)=psical(recep(m),rselsum,zdif)*tmu
          go to 73056
73054     continue
          geceppc(m,k)=psical(recep(m),rgrid(i),zdif)*tmu
73056    continue
          gecepc(m,k)=geceppc(m,k)-gecempc(m,k)
        enddo
        enddo
      enddo
!-----------------------------------------------------------------------
!-- Ohmic coils   receec(nece,nesum)                                 --
!-----------------------------------------------------------------------
      if (iecurr.le.0) go to 1710
      do 1700 n=1,nesum
        call sets2d(gridec(1,n),c,rgrid,nw,bkx,lkx,zgrid,nh,bky, &
                       lky,wk,ier)
        if (eceiter.eq.'pair') then
        do 1650 i=1,nece
        call seva2d(bkx,lkx,bky,lky,c,recem(i),zece(i),pds,ier,n111)
        receec(i,n)=pds(1)
 1650   continue
        endif
        if (eceiter.eq.'flux') then
        do i=1,nece
          receec(i,n)=0.0
        enddo
        endif
        do 1660 i=1,nece
        call seva2d(bkx,lkx,bky,lky,c,recep(i),zece(i),pds,ier,n111)
        receec(i,n)=pds(1)-receec(i,n)
 1660   continue
 1700 continue
 1710 continue
!----------------------------------------------------------------------
      open(unit=nffile,status='old',form='unformatted',err=12917, &
           file='recexx.dat')
      close(unit=nffile,status='delete')
12917  continue
      open(unit=nffile,status='new',form='unformatted', &
           file='recexx.dat')
      nnnece=nece
      write (nffile) nnnece 
      do 2150 ii=1,nece
        rrrecem(ii)=recem(ii)
        rrrecep(ii)=recep(ii)
 2150 continue
      rrreceo=receo
      write (nffile) rrrecem
      write (nffile) rrrecep
      write (nffile) rrreceo
      write (nffile) recebzfc
      write (nffile) gecebzpc
      write (nffile) recebzec
      write (nffile) recefc
      write (nffile) gecepc 
      write (nffile) receec
      close(unit=nffile)
!-----------------------------------------------------------------------
!--    do every time from here
!-----------------------------------------------------------------------
60000  continue
      ksetece=ksetece+1
      if (ksetece.eq.mtxece) ksetece=0
!-------------------------------------------------------------------
!--    get d(psi)/dr and d2(psi)/dr2
!-------------------------------------------------------------------
      call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
      do 208 iw=1,nw
         rw=rgrid(iw)
         rh=zgrid(jo)
         kk=(iw-1)*nh+jo
         call seva2d(bkx,lkx,bky,lky,c,rw,rh,pds,ier,n555)
         if (ier.eq.0) go to 206
         write (nttyo,9910) ier,rw,rh
         return
 206     continue
         dsidr(iw)=pds(2)
         ddsiddr(iw)=pds(5)
 208   continue
!
      if (eceiter.eq.'flux') then
        do i=1,nnece
         rw=recem(i)
         rh=zece(i)
         call seva2d(bkx,lkx,bky,lky,c,rw,rh,pds,ier,n555)
         if (ier.eq.0) go to 211
         write (nttyo,9910) ier,rw,rh
         return
 211     continue
         brspece(jtime,i)=pds(1)
        enddo
        if (idebug.ge.3) write (6,*) 'SETECE, recem/brspece = ',(recem(ii), &
            brspece(jtime,ii),ii=1,nece)
      endif
 9910 format('  error in spline =',i4,' (r,z)= ( ',f5.2,',',f5.2,')')
!----------------------------------------------------------------------
!-- give error (ecebit,ecebzbit)                                     --
!-- if kfixro kfixrece=0, robit and ecebit  from getecer             --
!-- if kfixro kfixrece=1, robit and rmbit, rpbit from k-file         --
!--   ecebit=sqrt((dpsi/dR|m*rmbit)**2+(dpsi/dR|p*rpbit)**2)         --
!--   ecebzbit=(d2psi/dR2)*robit/receo                               --
!----------------------------------------------------------------------
      if (kfitece.eq.1) go to 350
      if (kfixrece.eq.1) then
        do k=1,nece
         kp=iwp(k)
         km=iwm(k)
         ecebit(k)=(dsidr(kp)*rpbit(k))**2
         ecebit(k)=ecebit(k)+(dsidr(km)*rmbit(k))**2
         ecebit(k)=sqrt(ecebit(k))
        enddo
      endif
350   continue
      ddsiddro=ddsiddr(iwo)
      ecebzbit=ddsiddro*robit/receo
      if (kfitece.eq.1) go to 468
      do 360 m=1,nece
!       tdata1=serror*abs(brspece(jtime,m))
        tdata1=eceerror*abs(brspece(jtime,m))
        tdata2=abs(ecebit(m))
        tdata=max(tdata1,tdata2)
        if (tdata.gt.1.0e-10_dp) fwtece(m)=fwtece0(m)/tdata
        if (tdata.le.1.0e-10_dp) fwtece(m)=0.0
  360 continue
      do 466 i=1,nece
        if (fwtece(i).gt.0.0) kece=kece+1
  466 continue
  468 continue
      tdata1=serror*abs(brspecebz(jtime))
      tdata2=abs(ecebzbit)
      tdata=max(tdata1,tdata2)
      if (tdata.gt.1.0e-10_dp) fwtecebz=fwtecebz0/tdata
      if (tdata.le.1.0e-10_dp) fwtecebz=0.0
      if (fwtecebz.gt.0.0) kecebz=kecebz+1
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
      if (idebug.ge.3) then 
         write (6,*) 'SETECE, serror/eceerror = ',serror,eceerror
         write (6,*) 'SETECE, nece/fwtece = ',nece,(fwtece(i),i=1,nece)
      endif

      DEALLOCATE(dsidr,ddsiddr)
!
      return
 6530 format ('recefile_0',i1,'.ddd')
 6540 format ('recefile_',i2,'.ddd')
      end subroutine setece

!**********************************************************************
!>
!!    setff sets up the basis functions for Er.
!!    
!!
!!    @param ypsi :
!!
!!    @param xpsii :
!!
!**********************************************************************
      subroutine seter(ypsi,xpsii)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8, intent(inout) :: ypsi
      real*8, intent(inout) :: xpsii(keecur)

      do i=1,keecur
         xpsii(i) = bserel(keefnc,i,ypsi)
      enddo
      return
      end subroutine seter


!**********************************************************************
!>
!!    seterp computes derivative of er.
!!    
!!
!!    @param ypsi :
!!
!!    @param xpsii :
!!
!**********************************************************************
      subroutine seterp(ypsi,xpsii)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8, intent(inout) :: ypsi
      real*8, intent(inout) :: xpsii(kffcur)

      do i=1,kffcur
        xpsii(i) = bserpel(keefnc,i,ypsi)
      enddo
      return
      end subroutine seterp

!**********************************************************************
!>
!!    This subroutine sets the P' and FF' basis funciton parameters
!!    
!!
!**********************************************************************
      subroutine set_basis_params()

      use set_kinds
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

!---------------------------------------------------------------------
!--  specific choice of current profile                             --
!--       ICPROF=1  no edge current density allowed                 --
!--       ICPROF=2  free edge current density                       --
!--       ICPROF=3  weak edge current density constraint            --
!---------------------------------------------------------------------
      if (icprof.eq.1) then
        kffcur=2
        kppcur=2
        fcurbd=1.
        pcurbd=1.
        fwtbp=1.
        fwtqa=0.
        qvfit=0.
      elseif (icprof.eq.2) then
        kffcur=2
        kppcur=2
        fcurbd=0.
        pcurbd=0.
        fwtbp=0.
        fwtqa=0.
        qvfit=0.
      elseif (icprof.eq.3) then
        kffcur=3
        kppcur=2
        fcurbd=0.
        pcurbd=0.
        fwtbp=0.
        fwtqa=0.
        qvfit=0.
        kcalpa=1
        calpa(1,1)=0.1_dp
        calpa(2,1)=0.1_dp
        calpa(3,1)=0.1_dp
        xalpa(1)=0.0
        kcgama=1
        cgama(1,1)=0.1_dp
        cgama(2,1)=0.1_dp
        cgama(3,1)=0.1_dp
        xgama(1)=0.0
      endif
      if(mse_usecer .eq. 1)keecur = 0
      if(mse_usecer .eq. 2 .and. keecur .eq. 0) then
           keecur = 2
           keefnc = 0
           itek = 5
      endif
      if (imagsigma.gt.0) then
         do_spline_fit=.false.
         saimin=300.
      endif
!---------------------------------------------------------------------
!-- adjust fit parameters based on basis function selected          --
!---------------------------------------------------------------------
       if (kppfnc .eq. 3) then
          kppcur = 4 * (kppknt - 1)
       endif
       if (kppfnc .eq. 4) then
          kppcur = 4 * (kppknt - 1)
       endif
       if (kppfnc .eq. 5) then
          kppcur = kppcur * (kppknt - 1)
       endif
       if (kppfnc .eq. 6) then
          kppcur = kppknt * 2
       endif
       if (kfffnc .eq. 3) then
          kffcur = 4 * (kffknt - 1)
       endif
       if (kfffnc .eq. 4) then
          kffcur = 4 * (kffknt - 1)
       endif
       if (kfffnc .eq. 5) then
          kffcur = kffcur * (kffknt - 1)
       endif
       if (kfffnc .eq. 6) then
          kffcur = kffknt * 2
       endif
       if (kwwfnc .eq. 3) then
          kwwcur = 4 * (kwwknt - 1)
       endif
       if (kwwfnc .eq. 4) then
          kwwcur = 4 * (kwwknt - 1)
       endif
       if (kwwfnc .eq. 5) then
          kwwcur = kwwcur * (kwwknt - 1)
       endif
       if (kwwfnc .eq. 6) then
          kwwcur = kwwknt * 2
       endif
       if (keecur.gt.0) then
       if (keefnc .eq. 3) then
          keecur = 4 * (keeknt - 1)
       endif
       if (keefnc .eq. 4) then
          keecur = 4 * (keeknt - 1)
       endif
       if (keefnc .eq. 5) then
          keecur = keecur * (keeknt - 1)
       endif
       if (keefnc .eq. 6) then
          keecur = keeknt * 2
       endif
       endif
!
      if (kzeroj.eq.1.and.sizeroj(1).lt.0.0) sizeroj(1)=psiwant

      end subroutine set_basis_params


!**********************************************************************
!>
!!    setff sets up the basis functions for ff-ff(1).
!!    
!!
!!    @param ypsi :
!!
!!    @param xpsii :
!!
!**********************************************************************
      subroutine setff(ypsi,xpsii)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8, intent(inout) :: ypsi
      real*8, intent(inout) :: xpsii(kffcur)

      do i=1,kffcur
        xpsii(i) = bsffin(kfffnc,i,ypsi)
      enddo
      return
      end subroutine setff


!**********************************************************************
!>
!!    setfp sets up the basis functions for ffp.
!!    
!!
!!    @param ypsi :
!!
!!    @param xpsii :
!!
!**********************************************************************
      subroutine setfp(ypsi,xpsii)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8, intent(inout) :: ypsi
      real*8, intent(inout) :: xpsii(kffcur)

      do i=1,kffcur
        xpsii(i) = bsffel(kfffnc,i,ypsi)
      enddo
      return
      end subroutine setfp


!**********************************************************************
!>
!!    setfpp computes derivative of ffp.
!!    
!!
!!    @param ypsi :
!!
!!    @param xpsii :
!!
!**********************************************************************
      subroutine setfpp(ypsi,xpsii)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8, intent(inout) :: ypsi
      real*8, intent(inout) :: xpsii(kffcur)

      do i=1,kffcur
        xpsii(i) = bsffpel(kfffnc,i,ypsi)
      enddo
      return
      end subroutine setfpp


!**********************************************************************
!>
!!    setpp sets up the basis functions for pp.
!!    
!!
!!    @param ypsi :
!!
!!    @param xpsii :
!!
!**********************************************************************
      subroutine setpp(ypsi,xpsii)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8, intent(inout) :: ypsi
      real*8, intent(inout) :: xpsii(kppcur)

      do i=1,kppcur
        xpsii(i) = bsppel(kppfnc,i,ypsi)
      enddo
      return
      end subroutine setpp


!**********************************************************************
!>
!!    setppp computes derivative of pp.
!!    
!!
!!    @param ypsi :
!!
!!    @param xpsii :
!!
!**********************************************************************
      subroutine setppp(ypsi,xpsii)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8, intent(inout) :: ypsi
      real*8, intent(inout) :: xpsii(kppcur)

      do i=1,kppcur
        xpsii(i) = bspppel(kppfnc,i,ypsi)
      enddo
      return
      end subroutine setppp


!**********************************************************************
!>
!!    setpr sets up the basis functions for p-p(1).
!!    
!!
!!    @param ypsi :
!!
!!    @param xpsii :
!!
!**********************************************************************
      subroutine setpr(ypsi,xpsii)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8, intent(inout) :: ypsi
      real*8, intent(inout) :: xpsii(kppcur)

      do i=1,kppcur
        xpsii(i) = bsppin(kppfnc,i,ypsi)
      enddo
      return
      end subroutine setpr


!**********************************************************************
!>
!!    setpw sets up the basis functions for pw-pw(1).
!!    
!!
!!    @param ypsi :
!!
!!    @param xpsii :
!!
!**********************************************************************
      subroutine setpw(ypsi,xpsii)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8, intent(inout) :: ypsi
      real*8, intent(inout) :: xpsii(kwwcur)

      do i=1,kwwcur
        xpsii(i) = bswwin(kwwfnc,i,ypsi)
      enddo
      end subroutine setpw

!**********************************************************************
!>
!!    setpwp sets up the basis functions for pwp.
!!    
!!
!!    @param ypsi :
!!
!!    @param xpsii :
!!
!**********************************************************************
      subroutine setpwp(ypsi,xpsii)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8, intent(inout) :: ypsi
      real*8, intent(inout) :: xpsii(kwwcur)

      do i=1,kwwcur
        xpsii(i) = bswwel(kwwfnc,i,ypsi)
      enddo
      return
      end subroutine setpwp


!**********************************************************************
!>
!!    setpwpp computes derivative of pwp.
!!    
!!
!!    @param ypsi :
!!
!!    @param xpsii :
!!
!**********************************************************************
      subroutine setpwpp(ypsi,xpsii)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8, intent(inout) :: ypsi
      real*8, intent(inout) :: xpsii(kwwcur)

      do i=1,kwwcur
        xpsii(i) = bswwpel(kwwfnc,i,ypsi)
      enddo
      return
      end subroutine setpwpp


!**********************************************************************
!>
!!    setstark obtains the response matrix
!!    for polarimetry measurement
!!    
!!
!!    @param jtime : time index
!!
!**********************************************************************
      subroutine setstark(jtime)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      dimension rsplt(2500),zsplt(2500),csplt(2500) &
                ,rrgamin(nstark),zzgamin(nstark)
      character*40 filenmme
      data nset/20/,cdum/1.0/
      save nset
!---------------------------------------------------------------------
!--  try read in, first locally, then efit area                     --
!---------------------------------------------------------------------
      do 20 iname=1,nset
        if (iname.le.9) then
          write(filenmme,6530) iname
        else
          write(filenmme,6540) iname
        endif
        open(unit=nffile, &
             status='old',form='unformatted', &
             file=filenmme,err=10)
        go to 15
   10   continue
        open(unit=nffile, &
             status='old',form='unformatted', &
             file=table_dir(1:ltbdir)//filenmme,err=30)
   15   continue
        read (nffile,err=30) kkstark
        if (kkstark.ne.nstark) then
          close(unit=nffile)
          go to 20
        endif
        read (nffile,err=30)  rrgamin
        read (nffile,err=30)  zzgamin
        do ii=1,nstark
          if (abs(rrgamin(ii)-rrgam(jtime,ii)).gt.1.e-4_dp) go to 20
          if (abs(zzgamin(ii)-zzgam(jtime,ii)).gt.1.e-4_dp) go to 20
        enddo
        read (nffile,err=30) rbrfc
        read (nffile,err=30) rbzfc
        read (nffile,err=30) gbrpc
        read (nffile,err=30) gbzpc
        read (nffile,err=30) rbrec
        read (nffile,err=30) rbzec
        close(unit=nffile)
        go to 25
   20 continue
      go to 30
   25 return
!---------------------------------------------------------------------
!-- response due to F coils                                         --
!---------------------------------------------------------------------
   30 isplit=10
      itot=isplit*isplit
      fitot=itot
      rbrfc = 0.0
      rbzfc = 0.0
      do k=1,mfcoil
      call splitc(isplit,rsplt,zsplt,csplt, &
                  rf(k),zf(k),wf(k),hf(k),af(k),af2(k),cdum)
      do mmm=1,nstark
        if (rrgam(jtime,mmm).le.1.e-8_dp)  cycle
        brct=0.0
        bzct=0.0
        r1=rrgam(jtime,mmm)
        do l=1,itot
            a=rsplt(l)
            z1=zzgam(jtime,mmm)-zsplt(l)
            brct=brct+br(a,r1,z1)
            bzct=bzct+bz(a,r1,z1)
        enddo
        kkm=fcid(k)
        rbrfc(mmm,kkm)=rbrfc(mmm,kkm)+brct/fitot*tmu
        rbzfc(mmm,kkm)=rbzfc(mmm,kkm)+bzct/fitot*tmu
      enddo
      enddo
!---------------------------------------------------------------------
!--   plasma response                                               --
!---------------------------------------------------------------------
      zzx=0.0
      drzg=(drgrid**2+dzgrid**2)*25.
      isplit=18
      kk = 1
      do mmm=1,nstark
        gbrpc(mmm,kk)=0.0
        gbzpc(mmm,kk)=0.0
        if (rrgam(jtime,mmm).le.1.e-8_dp)  go to 90
        r=rrgam(jtime,mmm)
        do ii=1,nw
          a=rgrid(ii)
          do jj=1,nh
            z=zzgam(jtime,mmm)-zgrid(jj)
            kk=(ii-1)*nh+jj
            dist=(a-r)**2+z**2
            if (dist.lt.drzg) then
   68        call splitc(isplit,rsplt,zsplt,csplt, &
                  rgrid(ii),zgrid(jj),drgrid,dzgrid,zzx,zzx,cdum)
             gbrt=0.0
             gbzt=0.0
             itot=isplit*isplit
             do k=1,itot
               a=rsplt(k)
               z=zzgam(jtime,mmm)-zsplt(k)
               distt=(a-r)**2+z**2
               if (distt.lt.1.e-8_dp.and.isplit.lt.49) then
                 isplit=isplit+2
                 go to 68
               endif
               gbrt=gbrt+br(a,r,z)
               gbzt=gbzt+bz(a,r,z)
             enddo
             gbrpc(mmm,kk)=gbrt*tmu/itot
             gbzpc(mmm,kk)=gbzt*tmu/itot
            else
             gbrpc(mmm,kk)=br(a,r,z)*tmu
             gbzpc(mmm,kk)=bz(a,r,z)*tmu
            endif
          enddo
        enddo
      enddo
   90 continue
!---------------------------------------------------------------------
!-- E coils response                                                --
!---------------------------------------------------------------------
      isplit=10
      itot=isplit*isplit
      fitot=real(itot,dp)
      do 201 m=1,nstark
        do i=1,nesum
          rbrec(m,i)=0.0
          rbzec(m,i)=0.0
        enddo
      if (rrgam(jtime,m).le.1.e-8_dp) go to 201
      r1=rrgam(jtime,m)
      do 200 k=1,necoil
        brct=0.0
        bzct=0.0
        call splitc(isplit,rsplt,zsplt,csplt, &
                    re(k),ze(k),we(k),he(k),zzx,zzx,cdum)
        do l=1,itot
          a=rsplt(l)
          z1=zzgam(jtime,m)-zsplt(l)
          brct=brct+br(a,r1,z1)
          bzct=bzct+bz(a,r1,z1)
        enddo
      brct=brct*tmu/fitot
      bzct=bzct*tmu/fitot
      kkm=ecid(k)
      rbrec(m,kkm)=rbrec(m,kkm)+brct
      rbzec(m,kkm)=rbzec(m,kkm)+bzct
  200 continue
  201 continue
!
! --- write out rstarkxx.dat if flag IOUT contains 8.
!
      if (iand(iout,8).ne.0) then
      open(unit=nffile,status='old',form='unformatted',err=12917, &
           file='rstarkxx.dat')
      close(unit=nffile,status='delete')
12917  continue
      open(unit=nffile,status='new',form='unformatted', &
           file='rstarkxx.dat')
      kkstark=nstark
      write (nffile) kkstark
      do 215 ii=1,nstark
        rrgamin(ii)=rrgam(jtime,ii)
        zzgamin(ii)=zzgam(jtime,ii)
  215 continue
      write (nffile) rrgamin
      write (nffile) zzgamin
      write (nffile) rbrfc
      write (nffile) rbzfc
      write (nffile) gbrpc
      write (nffile) gbzpc
      write (nffile) rbrec
      write (nffile) rbzec
      close(unit=nffile)
      endif
      return
 6530 format ('rs129129_0',i1,'.ddd')
 6540 format ('rs129129_',i2,'.ddd')
      end subroutine setstark

!END_SET_BASIS


