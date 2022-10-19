#include "config.f"
!**********************************************************************
!>
!!    This subroutine gets the beam pressure.
!!
!**********************************************************************
      subroutine getbeam()
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      dimension bwork(ndata),cwork(ndata),dwork(ndata)
!
      if (nbeam.lt.0) return
!-----------------------------------------------------------------
!--   interpolate beam pressure into Thomson grid               --
!-----------------------------------------------------------------
      call zpline(nbeam,sibeam,pbeam,bwork,cwork,dwork)
      do i=1,npress
        xn=-rpress(i)
        pbimth(i)=seval(nbeam,xn,sibeam,pbeam,bwork,cwork,dwork)
      enddo
      pbeamb=seval(nbeam,x111,sibeam,pbeam,bwork,cwork,dwork)
      pbimpb=speval(nbeam,x111,sibeam,pbeam,bwork,cwork,dwork)
!-----------------------------------------------------------------
!--   interpolate beam ion density into Thomson grid            --
!-----------------------------------------------------------------
      call zpline(nbeam,sibeam,dnbeam,bwork,cwork,dwork)
      do i=1,npress
        xn=-rpress(i)
        dnbthom(i)=seval(nbeam,xn,sibeam,dnbeam,bwork,cwork,dwork)
      enddo
      return
!------------------------------------------------------------------
!--   compute beam pressure analytically                         --
!------------------------------------------------------------------
      return
      end subroutine getbeam

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
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      parameter (nn=30)
      parameter (kbre=5)
      dimension pds(6),nnout(kbre),bmink(kbre),bmaxk(kbre) 
      real*8,allocatable :: rrgrid(:,:),bfield(:),rrout(:,:), &
          bout(:,:),babs(:,:),bbb(:),ccc(:),ddd(:),btttt(:), &
          dsidr(:),ddsiddr(:),bx(:),ry(:),bbk(:),dbdr(:)
      dimension arspfit(nnecein,nn),brspfit(nnecein) &
           ,s(nn),tte(nnecein),x(nn),an(nnecein) &
           ,tebit(nnecein)
      dimension telowf(nnnte) &
         ,blowf(nnnte),bb(nnnte),cc(nnnte),dd(nnnte) &
         ,teece(nnece),pteprm(nnece),pteprp(nnece) &
         ,idestp(nnece),idestm(nnece),becem(nnece),becep(nnece) &
         ,dbdrp(nnece),dbdrm(nnece)
      integer*4, intent(inout) :: kerror

      kerror = 0
!-------------------------------------------------------------------
      allocate(rrgrid(kbre,nw),bfield(nw),rrout(kbre,nw), &
         bout(kbre,nw),babs(kbre,nw),bbb(nw),ccc(nw), &
         ddd(nw),btttt(nw),dsidr(nw),ddsiddr(nw),bx(nw), &
         ry(nw),bbk(nw),dbdr(nw))
!
      telowf=0 &
         ;blowf=0;bb=0;cc=0;dd=0 &
         ;teece=0;pteprm=0;pteprp=0 &
         ;idestp=0;idestm=0;becem=0;becep=0 &
         ;dbdrp=0;dbdrm=0

!-------------------------------------------------------------------
      do k=1,nnece
         fwtece0(k)=swtece(k)
      enddo
      fwtecebz0=swtecebz
      do k=1,kbre
       do i=1,nw
         babs(k,i)=0.0
         bout(k,i)=0.0
         rrout(k,i)=0.0
         rrgrid(k,i)=0.0
       enddo
      enddo
!-------------------------------------------------------------------
!--   input kgeteceb=0 from input file
!-------------------------------------------------------------------
      getECE: if (kgeteceb.le.0) then
      kgeteceb=kgeteceb+1
!---------------------------------------------------------------------
!--   Calculation of |B| array from fe array ( harmonic nharm)      --
!--     becein(necein),   fe(GHz),|B|(T)                            --
!--     !!! becein  from Low field to high field !!!                --
!---------------------------------------------------------------------
      do k=1,necein
        becein(k)=0.001_dp*6.0*9.1095_dp*pi/4.8032_dp*feece0(k)/nharm
      enddo
!EALW      write(*,*)'becein'
!EALW      write(*,*)becein
!--------------------------------------------------------------------
!--   fitting data from teecein0,errorece0 and becein (nnecein)    --
!--     bbx=(B-b00)/baa                                            --
!--     Te=x(1)+x(2)*bbx+x(3)*bbx**2+...+x(nfit)*bbx**(nfit-1)     --
!--------------------------------------------------------------------
!heng          mm--nnecein   m---necein  nn--parameter, n--nfit input
      binmin=becein(1)
      binmax=becein(necein)
      baa=0.5_dp*(binmax-binmin)
      b00=0.5_dp*(binmax+binmin)
      do i=1,necein
        an(i)=(becein(i)-b00)/baa
        tebit(i)=max(errorece0(i),1.e-4_dp)
      enddo
      do nj=1,necein
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
      do nj=1,necein
        brspfit(nj)=teecein0(nj)/tebit(nj)
      enddo
!
      mnow=necein
      if (kcmin.gt.0) then
        fwtnow=0.001_dp
        fwtcm =1.0
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
      nnn1=1
      iieerr=0
      call sdecm(arspfit,nnecein,mnow,nfit,brspfit,nnecein,nnn1,s,wk,iieerr)
      if (iieerr.eq.129) then
        kerror = 1
        call errctrl_msg('geteceb','sdecm failed to converge')
        return
      endif
      toler=1.0e-06_dp*s(1)
      do I = 1,nfit
        T = 0.0
        if (S(I).gt.toler) T = Brspfit(I)/S(I)
        Brspfit(I) = T
      enddo
      do I = 1, Nfit
        X(I) = 0.0
        do J = 1,nfit
          X(I) = X(I) + Arspfit(I,J)*Brspfit(J)
        enddo
      enddo
      do k=1,nfit
        xfit(k)=x(k)
      enddo
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
!
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
        call currnt(n222,jtime,n222,n222,kerror)
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
        if(sumf .ge. 0.0) then
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
        if (xpsi(kk).gt.1.0.or.ivacum.gt.0) then
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
        if (beceo.ge.bmaxk(1)) then
          receo=seval(n,beceo,bx,ry,bbb,ccc,ddd)
        endif
        do m=1,nece
          if (becep(m).ge.bmaxk(1)) then
            recem(m)=seval(n,becep(m),bx,ry,bbb,ccc,ddd)
          endif
          if (becem(m).ge.bmaxk(1)) then
            recep(m)=seval(n,becem(m),bx,ry,bbb,ccc,ddd)
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
          if ((beceo.ge.bmaxk(k)).and.(beceo.le.bmink(k-1))) then
            receo=seval(n,beceo,bx,ry,bbb,ccc,ddd)
          endif
          do m=1,nece
            if ((becep(m).ge.bmaxk(k)).and.(becep(m).le.bmink(k-1))) then
              recem(m)=seval(n,becep(m),bx,ry,bbb,ccc,ddd)
            endif
            if ((becem(m).ge.bmaxk(k)).and.(becem(m).le.bmink(k-1))) then
              recep(m)=seval(n,becem(m),bx,ry,bbb,ccc,ddd)
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
        if (beceo.le.bmink(kmax)) then
           receo=seval(n,beceo,bx,ry,bbb,ccc,ddd)
        endif
        do m=1,nece
          if (becep(m).le.bmink(kmax)) then
            recem(m)=seval(n,becep(m),bx,ry,bbb,ccc,ddd)
          endif
          if (becem(m).le.bmink(kmax)) then
            recep(m)=seval(n,becem(m),bx,ry,bbb,ccc,ddd)
          endif
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
      if (fwtecebz0.gt.1.e-6_dp) then
        dbdro=seval(nw,receo,rgrid,dbdr,bbb,ccc,ddd)
      endif
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
      deallocate(rrgrid,bfield,rrout,bout,babs,bbb,ccc,ddd,btttt, &
                 dsidr,ddsiddr,bx,ry,bbk,dbdr)
!
      return
      end subroutine geteceb

!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          getecer obtains the receo, R+ R-                        **
!**          from ECE measurement data                               **
!**          if kfixro  kfixrece = 0 called in setece                **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          11/98..........first created, Cheng Zhang               **
!**     2013/08/07..........Update for real-space Ti                 **
!**                                                                  **
!**********************************************************************
      subroutine getecer(jtime,kerror)
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      parameter (nn=30)
      parameter (kbre=5)
      dimension pds(6),nnout(kbre),bmink(kbre),bmaxk(kbre)
      real*8,allocatable :: rrgrid(:,:),bfield(:),rrout(:,:), &
          bout(:,:),babs(:,:),bbb(:),ccc(:),ddd(:),btttt(:), &
          dsidr(:),ddsiddr(:),bx(:),ry(:),bbk(:)
      dimension arspfit(nnecein,nn),brspfit(nnecein) &
           ,s(nn),tte(nnecein),x(nn),an(nnecein) &
           ,tebit(nnecein)
      dimension telowf(nnnte) &
         ,rlowf(nnnte),bb(nnnte),cc(nnnte),dd(nnnte) &
         ,teece(nnece),pteprm(nnece),pteprp(nnece) &
         ,idestp(nnece),idestm(nnece)
      integer*4, intent(inout) :: kerror
!
      kerror = 0
      allocate(rrgrid(kbre,nw),bfield(nw),rrout(kbre,nw), &
         bout(kbre,nw),babs(kbre,nw),bbb(nw),ccc(nw), &
         ddd(nw),btttt(nw),dsidr(nw),ddsiddr(nw),bx(nw), &
         ry(nw),bbk(nw))
!
      do k=1,nnece
        fwtece0(k)=swtece(k) 
      enddo
      fwtecebz0=swtecebz
      do k=1,kbre
        do i=1,nw
          babs(k,i)=0.0
          bout(k,i)=0.0
          rrout(k,i)=0.0
          rrgrid(k,i)=0.0
        enddo
      enddo
!---------------------------------------------------------------------
!--   Calculation of |B| array from fe array ( harmonic nharm)      --
!--     becein(necein),   fe(GHz),|B|(T) becein form H.f to L.f     --
!---------------------------------------------------------------------
      do k=1,necein
        becein(k)=0.001_dp*6.0*9.1095_dp*pi/4.8032_dp*feece(k)/nharm
      enddo
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
        call currnt(n222,jtime,n222,n222,kerror)
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
      sumf=fpol(nw)**2/2.
      delsi=-(psibry-simag)/(nw-1)
      do i=1,nw-1
        sumf=sumf+0.5_dp*delsi*(ffprim(nw-i+1)+ffprim(nw-i))
        if(sumf .ge. 0.0) then
          fpol(nw-i)=sqrt(2.*sumf)*fpol(nw)/abs(fpol(nw))
        else
          fpol(nw-i)=fpol(nw)
        endif
      enddo
      call zpline(nw,sigrid,fpol,bbb,ccc,ddd)
      do iw=1,nw
        kk=(iw-1)*nh+jo
        if (xpsi(kk).gt.1.0.or.ivacum.gt.0) then
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
      if (kmax.ne.0) then
!--------------------------------------------------------------------
!--   get babsk array (kmax+1) is strictly increasing order
!-------------------------------------------------------------------- 
      do k=1,kmax
        if (bmaxk(k).lt.bmink(k)) then
          stop
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
      else ! kmax.eq.0
      do i=1,nw
        bx(i)=bfield(nw-i+1)
        ry(i)=rgrid(nw-i+1)
      enddo
      call zpline(nw,bx,ry,bbb,ccc,ddd)
      do m=1,necein
        recein(m)=seval(nw,becein(m),bx,ry,bbb,ccc,ddd)
        teeceinr(m)=teecein(m)
        tebit(m)=errorece(m)
      enddo
      mecein=necein
!
      endif ! kmax.eq.0
!--------------------------------------------------------------------
!--   fitting data from teeceinr,tebit and recein (nnecein)        --
!--     rx=(R-r00)/raa                                             --
!--     Te=x(1)+x(2)*rx+x(3)*rx**2+...+x(nfit)*rx**(nfit-1)        --
!--------------------------------------------------------------------
!heng          mm--nnecein   m---mecein  nn--parameter, n--nfit input
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
      mnow=mecein
      if (kcmin.gt.0) then
        fwtnow=0.001_dp
        fwtcm =1.0
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
      nnn1=1
      iieerr=0
      call sdecm(arspfit,nnecein,mnow,nfit,brspfit,nnecein,nnn1,s,wk,iieerr)
      if (iieerr.eq.129) then
        kerror = 1
        call errctrl_msg('getecer','sdecm failed to converge')
        return
      end if
      toler=1.0e-06_dp*s(1)
      do I = 1,nfit
        T = 0.0
        if (S(I).gt.toler) T = Brspfit(I)/S(I)
        Brspfit(I) = T
      enddo
      do I = 1, Nfit
        X(I) = 0.0
        do J = 1,nfit
          X(I) = X(I) + Arspfit(I,J)*Brspfit(J)
        enddo
      enddo
      do k=1,nfit
        xfit(k)=x(k)
      enddo
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
          if ((recep(k).gt.receo).and.(recep(k).lt.recein(mecein))) &
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
      deallocate(rrgrid,bfield,rrout,bout,babs,bbb,ccc,ddd,btttt, &
                 dsidr,ddsiddr,bx,ry,bbk)
!
      return
      end subroutine getecer

!**********************************************************************
!>
!!    gettir obtains the receo, R+ R-
!!    from Ti data
!!    kfixro = 0, kfixrece = 3, called from setece
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
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      parameter (nn=30)
      parameter (kbre=5)
      dimension pds(6),nnout(kbre),bmink(kbre),bmaxk(kbre)
      real*8,allocatable :: rrgrid(:,:),bfield(:),rrout(:,:), &
          bout(:,:),babs(:,:),bbb(:),ccc(:),ddd(:),btttt(:), &
          dsidr(:),ddsiddr(:),bx(:),ry(:),bbk(:)
      dimension arspfit(nnecein,nn),brspfit(nnecein) &
           ,s(nn),tte(nnecein),x(nn),an(nnecein) &
           ,tebit(nnecein)
      dimension telowf(nnnte) &
         ,rlowf(nnnte),bb(nnnte),cc(nnnte),dd(nnnte) &
         ,teece(nnece),pteprm(nnece),pteprp(nnece) &
         ,idestp(nnece),idestm(nnece)
      integer*4, intent(inout) :: kerror
!
#ifdef DEBUG_LEVEL3
      write (6,*) 'Enter GETTIR, kfitece/kfixrece = ',&
         kfitece, kfixrece
#endif
      kerror = 0
      allocate(rrgrid(kbre,nw),bfield(nw),rrout(kbre,nw), &
         bout(kbre,nw),babs(kbre,nw),bbb(nw),ccc(nw), &
         ddd(nw),btttt(nw),dsidr(nw),ddsiddr(nw),bx(nw), &
         ry(nw),bbk(nw))
!
      do k=1,nnece
        fwtece0(k)=swtece(k) 
      enddo
      fwtecebz0=swtecebz
      do k=1,kbre
        do i=1,nw
          babs(k,i)=0.0
          bout(k,i)=0.0
          rrout(k,i)=0.0
          rrgrid(k,i)=0.0
        enddo
      enddo
!---------------------------------------------------------------------
!--   Copy Ti  array                                                --
!---------------------------------------------------------------------
      do k=1,necein
        becein(k)=feece(k)
        recein(k)=becein(k)
        teeceinr(k)=teecein(k)
        tebit(k)=errorece(k)
      enddo
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
!--     mm--nnecein   m---mecein  nn--parameter, n--nfit input     --
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
      nnn1=1
      iieerr=0
      mnow=mecein
      call sdecm(arspfit,nnecein,mnow,nfit,brspfit,nnecein,nnn1,s,wk,iieerr)
      if (iieerr.eq.129) then
        kerror = 1
        call errctrl_msg('gettir','sdecm failed to converge')
        return
      end if
      toler=1.0e-06_dp*s(1)
      do I = 1,nfit
        T = 0.0
        if (S(I).gt.toler) T = Brspfit(I)/S(I)
        Brspfit(I) = T
      enddo
      do I = 1, Nfit
        X(I) = 0.0
        do J = 1,nfit
          X(I) = X(I) + Arspfit(I,J)*Brspfit(J)
        enddo
      enddo
!
      do k=1,nfit
        xfit(k)=x(k)
      enddo
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
      deallocate(rrgrid,bfield,rrout,bout,babs,bbb,ccc,ddd,btttt, &
                 dsidr,ddsiddr,bx,ry,bbk)
!
      return
      end subroutine gettir

!**********************************************************************
!>
!!    fixstark adjusts the internal pitch angles
!!    based on spatial averaging data
!!    
!!
!!    @param jtime :
!!
!!    @param kerror :
!!
!**********************************************************************
      subroutine fixstark(jtime,kerror)
      use commonblocks,only: ct,wkt,bkrt,bkzt
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      dimension pds(6)
      real*8,dimension(:),allocatable :: bwork,cwork,dwork
      integer*4, intent(inout) :: kerror

      kerror = 0
!
      allocate(bwork(nw),cwork(nw),dwork(nw))
!
      if (keecur .gt. 0 .and. ifirst .eq. 0) then
        write(6,*) "Spatial averaging correction of MSE data"
        write(6,*) "not supported with Er fit at this time."
        write(6,*) "Going ahead with correction anyway,at this"
        write(6,*) "time, on channels requested."
!        write(6,*) "Calculating but not applying correction."
!        where(mse_spave_on .ne. 0) mse_spave_on = -1
      endif

      ltime = jtime
      if (jtime .lt. 0) then
        ltime = -jtime
!       do ichan=1,nstark
!         tangam(ltime,ichan) = save_tangam(ltime,ichan)
!       enddo
      endif

      if ( ifirst .eq. 0) then
        ifirst = 1
        write(6,*)"Calculate pitch angle corrections", &
          " using spatial averaging"
        do ichan=1,nstark
          save_gam(ltime,ichan) = atan(tangam(ltime,ichan))
          save_tangam(ltime,ichan) = tangam(ltime,ichan)
        enddo
        return
      endif

      ip_sign = - cpasma(ltime) / abs(cpasma(ltime))

      call sets2d(psi,ct,rgrid,nw,bkrt,lkrt,zgrid,nh,bkzt,lkzt,wkt,ier)
  
      if (pasmat(ltime).gt.0.0) then
        ssimag=-simag
        ssibry=-psibry
      else
        ssimag=simag
        ssibry=psibry
      endif

      if (jtime .gt. 0.0) then
!---------------------------------------------------------------------
!--   set up P' and FF', then integration                             --
!--   ffprim = (RBt) * d/dpsi(RBt)                                    --
!---------------------------------------------------------------------
      select case (icurrt)
      case (1)
        pprime(1)=cratio*sbeta/darea/srma
        ffprim(1)=cratio*srma*2.*salpha/darea*twopi*tmu
        pprime(nw)=pprime(1)
        ffprim(nw)=ffprim(1)
      case (2,5)
        pprime(nw)=ppcurr(x111,kppcur)/darea
        ffprim(nw)=fpcurr(x111,kffcur)/darea*twopi*tmu
        pprime(1)=ppcurr(x000,kppcur)/darea
        ffprim(1)=fpcurr(x000,kffcur)/darea*twopi*tmu
        if (kfffnc.eq.8) then
          ffprec(nw)=fpecrr(x111,kffcur)/darea*twopi*tmu
          ffprec(1)=fpecrr(x000,kffcur)/darea*twopi*tmu
        else
          ffprec(nw)=0.0
          ffprec(1)=0.0
        endif
      case (4)
        call currnt(n222,iges,n222,n222,kerror)
        if (kerror.gt.0) return
        pprime(1)=cratio/darea/rzero
        ffprim(1)=rbetap*cratio*rzero*twopi*tmu/darea
        ffprim(nw)=ffprim(1)*gammaf
        pprime(nw)=pprime(1)*gammap
      end select

      do i=2,nw-1
        ii=nw-i+1
        siii=1.0_dp-1.0_dp/(nw-1)*(i-1)
        sigrid(ii)=siii
        select case (icurrt)
        case(1)
          pprime(ii)=pprime(1)
          ffprim(ii)=ffprim(1)
        case (2,5)
          pprime(ii)=ppcurr(siii,kppcur)/darea
          ffprim(ii)=fpcurr(siii,kffcur)/darea*twopi*tmu
          if (kfffnc.eq.8) then
            ffprec(ii)=fpecrr(siii,kffcur)/darea*twopi*tmu
          else
            ffprec(ii)=0.0
          endif
        case(4)
          pprime(ii)=(1.-siii**enp)**emp*(1.-gammap)+gammap
          ffprim(ii)=ffprim(1)*pprime(ii)
          pprime(ii)=pprime(1)*pprime(ii)
        end select
      enddo
      endif ! jtime .gt. 0.0

      fpol(nw)=fbrdy*tmu
      sumf=fpol(nw)**2/2.
      ! TODO: psimag is undefined... was ssimag intended?
!      delsi=-(psibry+psimag)/(nw-1)
      delsi=-(psibry+ssimag)/(nw-1)
      do i=1,nw-1
        sumf=sumf+0.5_dp*delsi*(ffprim(nw-i+1)+ffprim(nw-i))
        if (sumf .ge. 0.0) then
          fpol(nw-i)=sqrt(2.*sumf)*fpol(nw)/abs(fpol(nw))
        else
          fpol(nw-i)=fpol(nw)
        endif
      enddo

      call zpline(nw,sigrid,fpol,bwork,cwork,dwork)

      do ichan = 1,nmtark
        if (mse_spave_on(ichan) .ne. 0) then

          ttl = 0.0
          rl = rrgam(ltime,ichan)
          zl = zzgam(ltime,ichan)
          call seva2d(bkrt,lkrt,bkzt,lkzt,ct,rl,zl,pds,ier,n333)
          brl = -pds(3) / rl
          bzl = pds(2) / rl
          psi_norm = (ssimag -pds(1)/ip_sign)/(ssimag-ssibry)
          btl = seval(nw,abs(psi_norm),sigrid,fpol,bwork, &
            cwork,dwork) / rl
          tglocal = (bzl * a1gam(ltime,ichan)) /  &
            (btl * a2gam(ltime,ichan) + brl * a3gam(ltime,ichan) &
            + bzl * a4gam(ltime,ichan))


          do i = 1,ngam_u
            do j = 1,ngam_w
              rl = spatial_avg_gam(ichan,1,i,j)
              zl = spatial_avg_gam(ichan,2,i,j)
              call seva2d(bkrt,lkrt,bkzt,lkzt,ct,rl,zl,pds,ier,n333)
              brl = -pds(3) / rl
              bzl = pds(2) / rl
              psi_norm = (ssimag -pds(1)/ip_sign)/(ssimag-ssibry)
              btl = seval(nw,abs(psi_norm),sigrid,fpol,bwork, &
                cwork,dwork) / rl
              tl = 0.0
              tl = tl + spatial_avg_gam(ichan,4,i,j) * btl
              tl = tl + spatial_avg_gam(ichan,5,i,j) * brl
              tl = tl + spatial_avg_gam(ichan,6,i,j) * bzl
              tl = spatial_avg_gam(ichan,3,i,j) * bzl / tl
              ttl = ttl + tl
              !if(jtime .lt. 0) write(7,'(I2,6F13.8)') ichan,rl,zl,btl,brl,bzl,tl
            enddo
          enddo
          !spatial_fix(ichan,ltime) =
          ! atan(cmgam(ichan,ltime)) - atan(ttl)
          spatial_fix(ichan,ltime) =  &
            atan(tglocal) - atan(ttl)

          if (jtime.gt.0.and.mse_spave_on(ichan) .eq. 1) then
            tangam(ltime,ichan) = tan(save_gam(ltime,ichan) &
              - spatial_fix(ichan,ltime))
          endif
        endif
      enddo

      if (jtime .lt. 0) then
        ifirst = 0
      endif
      
      deallocate(bwork,cwork,dwork)

      return
      end subroutine fixstark

!**********************************************************************
!>
!!    GETSIGMA is the control for getting the uncertainty
!!    in Magnetic Data
!!    
!!
!!    @param jtime :  time slice number
!!
!!    @param nitera :  iteration number
!!
!**********************************************************************
      subroutine getsigma(jtimex,niterax)
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      integer*4 jtimex,niterax
      real*8 :: gradsdr,gradsdz,brdr,brdz,bzdr,bzdz,cost,sint,oldfit
      dimension pds(6)
!----------------------------------------------------------------------
!--   BR=-1/R dpsi/dZ           BZ=1/R dpsi/dR                       --
!--            Bprobe= BR cost + BZ sint                             --
!--            gradsdr= dBRdR cost + dBZdR sint         dBprobe/dR   --
!--            gradsdz= dBRdZ cost + dBZdZ sint         dBrpobe/dZ   --
!--            gradsmp= sqrt (gradsdr**2 + gradsdz**2)               --
!--            gradsfl= |flux loop gradient|                         --
!--            bpermp= -BR sint + BZ cost                            --
!----------------------------------------------------------------------
      if (jtimex*niterax.eq.1) then
        open(unit=33,file='getsigma.log',form='formatted',status='unknown')
        open(unit=34,file='mprobe.err',form='formatted',status='unknown')
        open(unit=35,file='siloop.err',form='formatted',status='unknown')
        open(unit=36,file='fcoil.err',form='formatted',status='unknown')
        open(unit=37,file='ecoil.err',form='formatted',status='unknown')
        open(unit=38,file='bandcur.err',form='formatted',status='unknown')
      endif
      do i=1,magpri
        call seva2d(bkx,lkx,bky,lky,c,xmp2(i),ymp2(i),pds,ier,n666)
!----------------------------------------------------------------------
!--     Calculate dBR/dr, dBR/dz, dBZ/dr, dBZ/dz                     --
!----------------------------------------------------------------------
        brdr=(-pds(4)+pds(3)/xmp2(i))/xmp2(i)
        brdz=(-pds(6))/xmp2(i)
        bzdr=(pds(5)-(pds(2)/xmp2(i)))/xmp2(i)
        bzdz=pds(4)/xmp2(i)
!----------------------------------------------------------------------
!--     Form dBprobe/dR and dBprobe/dZ, then gradient                --
!----------------------------------------------------------------------
        sinm=sin(radeg*amp2(i))
        cosm=cos(radeg*amp2(i))
        gradsdr=brdr*cosm + bzdr*sinm
        gradsdz=brdz*cosm + bzdz*sinm
        gradsmp(jtimex,i)=sqrt(gradsdr**2+gradsdz**2)
!----------------------------------------------------------------------
!--     Calculate B perpendicular to magnetic probe                  --
!----------------------------------------------------------------------
        bpermp(jtimex,i)=(pds(2)*cosm + pds(3)*sinm)/xmp2(i)
      enddo
!----------------------------------------------------------------------
!--   calc gradient of flux loops                                    --
!----------------------------------------------------------------------
      do i=1,nsilop
        call seva2d(bkx,lkx,bky,lky,c,rsi(i),zsi(i),pds,ier,n333)
        gradsfl(jtimex,i)=sqrt(pds(2)**2+pds(3)**2)
      enddo
      write(33,*) '#', ishot,time(jtimex),jtimex
      write(33,*)'#magprobe no.,gradsmp,bpermp'
      do i=1,magpri
        write(33,'(i5,2e13.5)')i,gradsmp(jtimex,i),bpermp(jtimex,i)
      enddo
      write(33,*)'#fluxloop no., gradsfl'
      do i=1,nsilop
        write(33,'(i5,2e13.5)') i,gradsfl(jtimex,i)
      enddo
      call magsigma(ishot,time(jtimex),jtimex,gradsmp,gradsfl, &
                    bpermp,sigmaf,sigmab,sigmae,sigmafl,sigmamp)
!----------------------------------------------------------------------
!--   Set fitting weights                                            --
!----------------------------------------------------------------------
      do i=33,38
        write(i,*) '#', ishot,time(jtimex),jtimex
        write(i,*) '#errorm=', errorm
        write(i,*) '#errmag=', errmag
      enddo
      write(38,*) '#sigmab'
      write(38,'(e13.5)') sigmab(jtimex)
      write(38,*) '#sigmaip0,sigmaip'
      write(38, '(2e13.5)') sigmaip0,sigmaip(jtimex)
      write(34,*) '#sigmamp0,sigmamp'
      do i=1,magpri
        write(34,'(i5,2e13.5)') i,sigmamp0(i),sigmamp(jtimex,i)
      end do
      write(34,*) ' '
      write(34,*) '#t0mp,mp_k, mprcg, vresmp, devmp'
      do i=1,magpri
        write(34,'(6e13.5)') t0xmp(i),xmp_k(i),xmprcg(i),vresxmp(i), &
                 devxmp(jtimex,i),rnavxmp(jtimex,i)
      enddo
      write(34,*) ' '
      write(36,*) '#sigmaf0,sigmaf'
      do i=1,nfcoil
        write(36,'(i5,2e13.5)') i,sigmaf0(i),sigmaf(jtimex,i)
      end do
      write(36,*) ' '
      write(37,*) '#sigmae0,sigmae'
      do i=1,nesum
        write(37,'(i5, 2e13.5)') i, sigmae0(i),sigmae(jtimex,i)
      end do
      write(37,*) ' '
      write(35,*) '#sigmafl0,sigmafl'
      do i=1,nsilop
        write(35,'(i5,2e13.5)') i, sigmafl0(i),sigmafl(jtimex,i)
      enddo
      write(35,*) ' '
      write(35,*) '#t0psi,psi_k, psircg, vrespsi, devpsi'
      do i=1,nsilop
        write(35,'(6e13.5)') t0psi(i),psi_k(i),psircg(i),vrespsi(i), &
               devpsi(jtimex,i),rnavpsi(jtimex,i)
      enddo
      write(35,*) ' '
      do i=1,nsilop
        oldfit = fwtsi(i)
        if (sigmafl(jtimex,i)/=0.0) then
          fwtsi(i)=swtsi(i)/sigmafl(jtimex,i)
        else
          fwtsi(i)=0.0
        endif
        write (99,*) i, swtsi(i), oldfit, fwtsi(i)
      enddo
      do i=1,magpri
        oldfit = fwtmp2(i)
        if (sigmamp(jtimex,i)/=0.0) then
          fwtmp2(i)=swtmp2(i)/sigmamp(jtimex,i)
        else
          fwtmp2(i)=0.0
         endif
         write (99,*) i, swtmp2(i), oldfit, fwtmp2(i)
      enddo
      do i=1,nesum
        oldfit = fwtec(i)
        if (sigmae(jtimex,i)/=0.0) then
          fwtec(i)=swtec(i)/sigmae(jtimex,i)
        else
          fwtec(i)=0.0
        endif
        write (99,*) i, swtec(i), oldfit, fwtec(i)
      enddo
      do i=1,nfcoil
        oldfit = fwtfc(i)
        if (sigmaf(jtimex,i)/=0.0) then
          fwtfc(i)=swtfc(i)/sigmaf(jtimex,i)
        else
          fwtfc(i)=0.0
        endif
        write (99,*) i, swtfc(i), oldfit, fwtfc(i)
      enddo
      oldfit = fwtcur
      if (sigmaip(jtimex)/=0.0) then
        fwtcur=swtcur/sigmaip(jtimex)
      else
        fwtcur=0.0
      endif
      write (99,*) swtcur, oldfit, fwtcur
!
      return
      end subroutine getsigma

!*************************************************************************
!>
!!    This subroutine calculates the uncertainties for the magnetic
!!    diagnostics.  It is based on estimates described in
!!    DIII-D Physics Memo D3DPM 0202, "Estimating the Uncertainty of DIII-D
!!    Magnetic Data," by E.J. Strait (Aug. 30, 2002).
!!    
!!    The following sources of uncertainty are included. (Further explanation
!!    is given under the corresponding item numbers in the Physics Memo.)
!!    The individual terms are combined in quadrature to give
!!    the total uncertainty for each signal.
!!    
!!    1) Loop calibration   dS = a1 S
!!    2) Loop cal. - Long-term change  dS = a2 S
!!    3) Integrator calibration  dS = a3 S
!!    4) Int. cal. - Long-term change  dS = a4 S
!!    5) Integrator drift    dS = a5 K(RC/G) T
!!    6) Loop position    dS = a6 grad(S)
!!    7) Loop tilt angle   dS = a7 Bperp
!!    8) Bt pickup     dS = a8 dBt
!!    9) C-coil pickup   dS = a9 Cn
!!    10) Bp pickup in leads   dS = a10 K integral(Bperp^2)ds
!!    11) Bp pickup in ports   dS = a11 K Bport
!!    12) Detector nonlinearity  dS = a12 S/Bt
!!    13) Noise     dS = a13 (/<S^2> - <S>^2/)^0.5 / N^0.5
!!    14) Digitizer resolution  dS = a14 K(RC/G)(Vres/2) / N^0.5
!!    where
!!    S  = measured signal: flux, field, or current (in physics units)
!!    dS  = estimated uncertainty in the measured signal
!!    grad(S)  = gradient of measured flux or field (physics units/meter)
!!    N   = number of time samples averaged
!!    Vres  = one-bit resolution of the digitizer (volts)
!!    K   = inherent number (physics units/volt-second)
!!    RC  = integrator time constant (seconds)
!!    G   = integrator gain
!!    T   = time elapsed since the start of integration (seconds)
!!    dBt  = change in toroidal field since start of integration (seconds)
!!    Bperp  = poloidal field normal to axis of mag. probe or leads (Tesla)
!!    Cn  = current in a C-coil pair (Amps)
!!    integral(f)ds = integral along leads from probe to connector (meters)
!!    Bport  = poloidal field in port, perpendicular to leads (Tesla)
!!    an  = numerical coefficient: units vary with n,
!!    and values vary between types of signals
!!    
!!    Note that items 10 and 11 will not be implemented immediately due to
!!    the additional difficulty in calculating the path integral (10) and in
!!    estimating Bport outside the efit grid (11).
!!
!!
!!    @param ishotx :  shot number
!!
!!    @param timexy : time slice
!!
!!    @param jtimex : time slice index
!!
!!    @param gradsmpx : Grad(S) of magnetic probe \n 
!!    s   = BR cost + BZ sint
!!
!!    @param gradsflx : Grad(S) of flux loop
!!
!!    @param bpermpx : B perpendicular to the magnetic probe
!!
!!    @param sigmafx : uncertainty of f coil
!!
!!    @param sigmabx : uncertainty of b coil
!!
!!    @param sigmaex : uncertainty of e coil
!!
!!    @param sigmaflx : uncert. of flux loop
!!
!!    @param sigmampx : uncert. of magnetic probes
!!
!**********************************************************************
      subroutine magsigma(ishotx,timexy,jtimex,gradsmpx,gradsflx, &
                          bpermpx,sigmafx,sigmabx,sigmaex, &
                          sigmaflx,sigmampx)

      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!
      dimension gradsmpx(ntime,magpri),gradsflx(ntime,nsilop)
      dimension bpermpx(ntime,magpri)
      dimension sigmafx(ntime,nfcoil),sigmabx(ntime)
      dimension sigmaex(ntime,nesum),sigmaflx(ntime,nsilop)
      dimension sigmampx(ntime,magpri)

      !-----
      ! LOCAL VARIABLES
      !-----
      dimension dd(8)
      equivalence (dd(1), dpol), &
                  (dd(2), drdp) , &
                  (dd(3), dfl) , &
                  (dd(4), df67) , &
                  (dd(5), dvv) , &
                  (dd(6), dfc) , &
                  (dd(7), dbe) , &
                  (dd(8), dip)

      ! rindex diagnostic  units

      ! 1 Pol. Mag. Probes   (Tesla)
      ! 2 RDP Mag. Probes   (Tesla)
      ! 3 F-coil Flux Loops (V-s/radian)
      ! 4 F6, F7 Flux Loops (V-s/radian)
      ! 5 Vessel Flux Loops (V-s/radian)
      ! 6 F-coil Rogowskis  (Amps)
      ! 7 E and B Rogowskis (Amps)
      ! 8 Ip Rogowskis   (Amps)

      !-----
      ! ARRAYS TO ACCUMULATE THE UNCERTAINTIES
      !-----
      real*8 :: sigmp(magpri),sigfl(nsilop),sigfc(nfcoil),sige(nesum)
      real*8 :: maskpol(magpri),maskrdp(magpri),maskff(nsilop), &
                maskinf(nsilop),maskoutf(nsilop),maskvv(nsilop)
      !-----
      ! MASKS FOR SUBSETS WITHIN ARRAYS
      !-----
!-- poloidal array probes
      maskpol = (/ 60*1., 16*0./)
!-- RDP probes
      maskrdp = (/ 60*0., 16*1./)
!-- F-coil flux loops
      maskff = (/18*1.,7*0.,1., 0., 1.,7*0.,1.,0.,1.,6*0./)
!-- inner F-coil flux loops
      maskinf = (/ 5*1., 2*0., 7*1., 2*0., 2*1., 26*0. /)
!-- outer F-coil  loops
      maskoutf = (/5*0., 2*1., 7*0., 2*1., 2*0., &
                   7*0., 1., 0., 1., 7*0., 1., 0., 1., 6*0./)
!-- vacuum vessel flux loops
      maskvv=(/18*0.,7*1.,0., 1., 0., 7*1., 0., 1.,0.,6*1./)
      !-----
      ! INITIALIZE ARRAYS
      !-----
      sigmp = (/ (0., i=1, magpri) /)
      sigfl = (/ (0., i=1, nsilop) /)
      sigfc = (/ (0., i=1, nfcoil) /)
      sige =  (/ (0., i=1, nesum)  /)
      sigb =  0.
      sigip = 0.

      !-----
      !  TESTING DIFFERENT CONTRIBUTIONS TO SIGMA
      !-----
      timex = timexy / 1000.
      if (ksigma .eq. 0 .or. ksigma .eq. 1) then
        !-------------------------------------------------------------------------
        !***** (1) LOOP CALIBRATION
        !-------------------------------------------------------------------------
        dd = (/.0019, .0019, 0., 0., 0., .0017, .0017, .003/)
        !-------------------------------------------------------------------------
        sigmp = sigmp + ( expmpi(jtimex,:) * dpol )**2
        sigfl = sigfl + 0
        sigfc = sigfc + ( fccurt(jtimex,:) * dfc  )**2
        sige  = sige  + ( eccurt(jtimex,:) * dbe  )**2
        sigb  = sigb  + ( bcentr(jtimex)   * dbe  )**2
        sigip = sigip + ( pasmat(jtimex)   * dip  )**2
      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 2) then
        !-------------------------------------------------------------------------
        !***** (2) LOOP CALIBRATION: LONG-TERM
        !-------------------------------------------------------------------------
        dd = (/.002, .002, 0., 0., 0., .003, .003, .006/)
        !-------------------------------------------------------------------------
        sigmp = sigmp + ( expmpi(jtimex,:) * dpol )**2
        sigfl = sigfl + 0
        sigfc = sigfc + ( fccurt(jtimex,:) * dfc  )**2
        sige  = sige  + ( eccurt(jtimex,:) * dbe  )**2
        sigb  = sigb  + ( bcentr(jtimex)   * dbe  )**2
        sigip = sigip + ( pasmat(jtimex)   * dip  )**2
      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 3) then
        !-------------------------------------------------------------------------
        !***** (3) INTEGRATOR CALIBRATION
        !-------------------------------------------------------------------------
        dd = (/.0013, .0013, .0013, .0013, .0013, .0013, .0013, .0013/)
        !-------------------------------------------------------------------------
        sigmp = sigmp + ( expmpi(jtimex,:) * dpol )**2
        sigfl = sigfl + ( silopt(jtimex,:) * dfl  )**2
        sigfc = sigfc + ( fccurt(jtimex,:) * dfc  )**2
        sige  = sige  + ( eccurt(jtimex,:) * dbe  )**2
        sigb  = sigb  + ( bcentr(jtimex)   * dbe  )**2
        sigip = sigip + ( pasmat(jtimex)   * dip  )**2
      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 4) then
        !-------------------------------------------------------------------------
        !***** (4) INTEGRATOR CAL.: LONG-TERM
        !-------------------------------------------------------------------------
        dd = (/.0017, .0017, .0017, .0017, .0017, .0017, .0017, .0017/)
        !-------------------------------------------------------------------------
        sigmp = sigmp + ( expmpi(jtimex,:) * dpol )**2
        sigfl = sigfl + ( silopt(jtimex,:) * dfl  )**2
        sigfc = sigfc + ( fccurt(jtimex,:) * dfc  )**2
        sige  = sige  + ( eccurt(jtimex,:) * dbe  )**2
        sigb  = sigb  + ( bcentr(jtimex)   * dbe  )**2
        sigip = sigip + ( pasmat(jtimex)   * dip  )**2
      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 5) then
        !-------------------------------------------------------------------------
        !***** (5) INTEGRATOR DRIFT
        !-------------------------------------------------------------------------
        dd = (/.0007, .0007, .0007, .0007, .0007, .0022, .0007, .0007/)
        !-------------------------------------------------------------------------
        sigmp = sigmp + (xmp_k  *xmprcg /vresxmp * (timex - t0xmp) &
          * dpol )**2
        sigfl = sigfl + ( psi_k * psircg/vrespsi * (timex - t0psi) &
          * dfl  )**2
        sigfc = sigfc + ( fc_k  * fcrcg /vresfc  * (timex - t0fc) &
          * dfc  )**2
        sige  = sige  + ( e_k   * ercg  /vrese   * (timex - t0e) &
          * dbe  )**2
        sigb  = sigb  + ( bc_k  * bcrcg /vresbc  * (timex - t0bc) &
          * dbe  )**2
        sigip = sigip + ( p_k   * prcg  /vresp   * (timex - t0p) &
          * dip  )**2
      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 6) then
        !-------------------------------------------------------------------------
        !***** (6) LOOP POSITION
        !-------------------------------------------------------------------------
        dd = (/.0020, .0020, .0030, .0045, .0020, 0., 0., 0./)
        !-------------------------------------------------------------------------
        !vas f90 modifi
        !vas sigmp = sigmp + ( gradsmp(jtimex,:) *
        !vas 1  (maskpol * dpol + maskrdp * drdp ) )**2
        !vas sigfl = sigfl + ( gradsfl(jtimex,:) *
        !vas 1  (maskinf * dfl + maskoutf * df67 + maskvv * dvv ) )**2
        sigmp = sigmp + ( gradsmp(jtimex,:) * &
          (maskpol * dpol + maskrdp * drdp ) )**2
        sigfl = sigfl + ( gradsfl(jtimex,:) * &
          (maskinf * dfl + maskoutf * df67 + maskvv * dvv ) )**2
        sigfc = sigfc + 0
        sige  = sige  + 0
        sigb  = sigb  + 0
        sigip = sigip + 0
      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 7) then
        !-------------------------------------------------------------------------
        !***** (7) LOOP TILT ANGLE - revised
        !-------------------------------------------------------------------------
        dd = (/.017, .017, 0., 0., 0., 0., 0., 0./)
        !-------------------------------------------------------------------------
        sigmp = sigmp + ( bpermp(jtimex,:) * dpol )**2
        sigfl = sigfl + 0
        sigfc = sigfc + 0
        sige  = sige  + 0
        sigb  = sigb  + 0
        sigip = sigip + 0
      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 8) then
        !-------------------------------------------------------------------------
        !***** (8) BT PICKUP
        !-------------------------------------------------------------------------
        dd = (/.003, .003, .00044, .00044, .00044, 0., 0., 1.3e4/)
        !-------------------------------------------------------------------------
        sigmp = sigmp + ( bti322(jtimex) * dpol )**2
        sigfl = sigfl + ( bti322(jtimex) * dfl  )**2
        sigfc = sigfc + 0
        sige  = sige  + 0
        sigb  = sigb  + 0
        sigip = sigip + ( bti322(jtimex) * dip  )**2
      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 9) then
        !-------------------------------------------------------------------------
        !**** (9a) C79 PICKUP
        !-------------------------------------------------------------------------
        dd  = (/2.3e-08,  1.4e-08,  5.1e-08,  5.1e-08, &
          3.6e-08, 0., 0., 0./)
        !-------------------------------------------------------------------------
        !vas f90 modifi
        sigmp = sigmp + ( curc79(jtimex) * &
          (maskpol * dpol + maskrdp * drdp ) )**2
        sigfl = sigfl + ( curc79(jtimex) * &
          (maskff * dfl   + maskvv * dvv   ) )**2

        !-------------------------------------------------------------------------
        !***** (9b) C139 PICKUP
        !-------------------------------------------------------------------------
        dd = (/8.1e-08,  2.3e-07,  4.1e-08,   4.1e-08, &
          3.2e-08, 0., 0., 0./)
        !-------------------------------------------------------------------------
        !vas f90 modifi
        sigmp = sigmp + ( curc139(jtimex) * &
          (maskpol * dpol + maskrdp * drdp ) )**2
        sigfl = sigfl + ( curc139(jtimex) * &
          (maskff * dfl   + maskvv * dvv   ) )**2

        !-------------------------------------------------------------------------
        !***** (9c) C199 PICKUP
        !-------------------------------------------------------------------------
        dd = (/3.1e-08,  1.1e-07,  5.2e-08,   5.2e-08, &
          3.9e-08, 0., 0., 0./)
        !-------------------------------------------------------------------------
        !vas f90 modifi
        sigmp = sigmp + ( curc199(jtimex) * &
          (maskpol * dpol + maskrdp * drdp ) )**2
        sigfl = sigfl + ( curc199(jtimex) * &
          (maskff * dfl   + maskvv * dvv   ) )**2
        sigfc = sigfc + 0
        sige  = sige  + 0
        sigb  = sigb  + 0
        sigip = sigip + 0
      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 10) then
        !-------------------------------------------------------------------------
        !***** (10) BP PICKUP IN LEADS
        !-------------------------------------------------------------------------
        dd = (/.00016,  .00016,  0.,  0., .00016, 0., 0., 0./)
        !-------------------------------------------------------------------------
        Brleads = 0.
        sigmp = sigmp + ( Brleads *xmp_k  * dpol )**2
        sigfl = sigfl + ( Brleads * psi_k * &
          (maskff * dfl + maskvv * dvv ) )**2
        sigfc = sigfc + 0
        sige  = sige  + 0
        sigb  = sigb  + 0
        sigip = sigip + 0
      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 11) then
        !-------------------------------------------------------------------------
        !***** (11) BP PICKUP IN PORT
        !-------------------------------------------------------------------------
        dd = (/.0002,  .0002,  0., 0., .0002, 0., 0., 0./)
        !-------------------------------------------------------------------------
        Bport = 0.
        sigmp = sigmp + ( Bport *xmp_k  * dpol )**2
        sigfl = sigfl + ( Bport * psi_k * &
          (maskff * dfl + maskvv * dvv ) )**2
        sigfc = sigfc + 0
        sige  = sige  + 0
        sigb  = sigb  + 0
        sigip = sigip + 0
      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 12) then
        !-------------------------------------------------------------------------
        !***** (12) DETECTOR NONLINEARITY
        !-------------------------------------------------------------------------
        dd = (/.0025, .0025, 0., 0., 0., 0., 0., 0./)
        !-------------------------------------------------------------------------
        if (abs(bcentr(jtimex)) .gt. 0.1) then
          sigmp = sigmp + ( expmpi(jtimex,:) &
            / bcentr(jtimex) * dpol )**2
        endif
        sigfl = sigfl + 0
        sigfc = sigfc + 0
        sige  = sige  + 0
        sigb  = sigb  + 0
        sigip = sigip + 0
      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 13) then
        !-------------------------------------------------------------------------
        !***** (13) NOISE
        !-------------------------------------------------------------------------
        dd = (/1., 1., 1., 1., 1., 1., 1., 1./)
        !-------------------------------------------------------------------------
        do i=1,magpri
          if (rnavxmp(jtimex,i).ne.0.0) &
            sigmp(i) = sigmp(i) + ( devxmp(jtimex,i) / &
            sqrt(rnavxmp(jtimex,i)) * dpol )**2
        enddo
        do i=1,nsilop
          if (rnavpsi(jtimex,i).ne.0.0) &
            sigfl(i) = sigfl(i) + ( devpsi(jtimex,i)/ &
            sqrt(rnavpsi(jtimex,i))* dfl  )**2
        enddo
        do i=1,nfcoil
          if (rnavfc(jtimex,i).ne.0.0) &
            sigfc(i) = sigfc(i) + ( devfc(jtimex,i) / &
            sqrt(rnavfc(jtimex,i)) * dfc  )**2
        enddo
        do i=1,nesum
          if (rnavec(jtimex,i).ne.0.0) &
            sige(i)  = sige(i)  + ( deve(jtimex,i)  / &
            sqrt(rnavec(jtimex,i))  * dbe  )**2
        enddo
        if (rnavbc(jtimex).ne.0.0) &
          sigb  = sigb  + ( devbc(jtimex)   / &
          sqrt(rnavbc(jtimex))   * dbe  )**2
        if (rnavp(jtimex).ne.0.0) &
          sigip = sigip + ( devp(jtimex)    / &
          sqrt(rnavp(jtimex))    * dip  )**2
      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 14) then
        !-------------------------------------------------------------------------
        !***** (14) DIGITIZER RESOLUTION
        !-------------------------------------------------------------------------
        dd = (/1., 1., 1., 1., 1., 1., 1., 1./)
        !-------------------------------------------------------------------------
        do i=1,magpri
          if (rnavxmp(jtimex,i).ne.0.0) &
            sigmp(i) = sigmp(i) + (xmp_k(i)  *xmprcg(i)  /2/ &
            sqrt(rnavxmp(jtimex,i)) * dpol )**2
        enddo
        do i=1,nsilop
          if (rnavpsi(jtimex,i).ne.0.0) &
            sigfl(i) = sigfl(i) + ( psi_k(i) * psircg(i) /2/ &
            sqrt(rnavpsi(jtimex,i))* dfl  )**2
        enddo
        do i=1,nfcoil
          if (rnavfc(jtimex,i).ne.0.0) &
            sigfc(i) = sigfc(i) + ( fc_k(i)  * fcrcg(i)  /2/ &
            sqrt(rnavfc(jtimex,i)) * dfc  )**2
        enddo
        do i=1,nesum
          if (rnavec(jtimex,i).ne.0.0) &
            sige(i)  = sige(i)  + ( e_k(i)   * ercg(i)   /2/ &
            sqrt(rnavec(jtimex,i))  * dbe  )**2
        enddo
        if (rnavbc(jtimex).ne.0.0) &
          sigb  = sigb  + ( bc_k  * bcrcg  /2/ &
          sqrt(rnavbc(jtimex))   * dbe  )**2
        if (rnavp(jtimex).ne.0.0) &
          sigip = sigip + ( p_k   * prcg   /2/ &
          sqrt(rnavp(jtimex))    * dip  )**2
      !------------------------------------------------------------
      endif

      sigmamp(jtimex,:) = sqrt(sigmp)
      sigmafl(jtimex,:) = sqrt(sigfl)
      sigmaf(jtimex,:)  = sqrt(sigfc)
      sigmae(jtimex,:) = sqrt(sige)
      sigmab(jtimex)  = sqrt(sigb)
      sigmaip(jtimex) = sqrt(sigip)

      ! print *,'sigmamp'
      ! print 99,sigmamp
      ! print *,'sigmafl'
      ! print 99,sigmafl
      ! print *,'sigmaf'
      ! print 99,sigmaf
      ! print *,'sigmae'
      ! print 99,sigmae
      ! print *,'sigmab'
      ! print 99,sigmab
      ! print *,'sigmaip'
      ! print 99,sigmaip

!99    format(5e12.4)

      return
      end subroutine magsigma

!**********************************************************************
!>
!!    GETTE gets the electron temperature
!!    profiles.
!!    
!!
!!    @param kerror : error flag
!!
!**********************************************************************
      subroutine gette(kerror)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      integer*4, intent(inout) :: kerror

      kerror = 0
!----------------------------------------------------------------------
!--   singular decomposition                                         --
!----------------------------------------------------------------------
      do nj=1,npress
        do nk=1,nptef
          xn=-rpress(nj)
          arsp_cw2(nj,nk)=xn**(nk-1)/sgteth(nj)
        enddo
        bdata_cw2(nj)=tethom(nj)/sgteth(nj)
      enddo
      ntedat=npress
      if (cstabte.gt.0.0) then
        nj=npress
        do jj=ncstte,nptef
          nj=nj+1
          do nk=1,nptef
            arsp_cw2(nj,nk)=0.0
            if (jj.eq.nk) arsp_cw2(nj,nk)=cstabte
          enddo
          bdata_cw2(nj)=0.0
        enddo
        ntedat=ntedat+nptef-ncstte+1
      endif
!---------------------------------------------------------------------
!--   form error matrix                                             --
!---------------------------------------------------------------------
      do i=1,nptef
        do j=1,nptef
          ematrix_cw2(i,j)=0.0
          do k=1,npress
            ematrix_cw2(i,j)=ematrix_cw2(i,j)+arsp_cw2(k,i)*arsp_cw2(k,j)
          enddo
        enddo
      enddo

      call sdecm(arsp_cw2,ndata,ntedat,nptef,bdata_cw2,ntedat,n111,wrsp_cw2,work_cw2,ier)
      if (ier.eq.129) then
        kerror = 1
        call errctrl_msg('gette','sdecm failed to converge')
        return
      end if

      cond=ier
      toler=1.0e-06_dp*wrsp_cw2(1)
      do i=1,nptef
        t=0.0
        if (wrsp_cw2(i).gt.toler) t=bdata_cw2(i)/wrsp_cw2(i)
        work_cw2(i)=t
      enddo
      do i=1,nptef
        tefit(i)=0.0
        do j=1,nptef
          tefit(i)=tefit(i)+arsp_cw2(i,j)*work_cw2(j)
        enddo
      enddo
!------------------------------------------------------------------
!--   compute chi square                                         --
!------------------------------------------------------------------
      chisqte=0.0
      do i=1,npress
        tenow=0.0
        xn=-rpress(i)
        do j=1,nptef
          tenow=tenow+tefit(j)*xn**(j-1)
        enddo
        chisqte=chisqte+((tenow-tethom(i))/sgteth(i))**2
      enddo
!---------------------------------------------------------------------
!--   get inverse of error matrix                                   --
!---------------------------------------------------------------------
      call linv1f(ematrix_cw2,nptef,nppcur,einv_cw2,n444,work_cw2,ier)
!----------------------------------------------------------------------
!--   boundary values                                                --
!----------------------------------------------------------------------
      tebdry=0.0
      stebdry=0.0
      tepbry=0.0
      sigtepb=0.0
      do j=1,nptef
        do i=1,nptef
          stebdry=stebdry+einv_cw2(i,j)
          sigtepb=sigtepb+(i-1)*(j-1)*einv_cw2(i,j)
        enddo
        tepbry=tepbry+(j-1)*tefit(j)
        tebdry=tebdry+tefit(j)
      enddo
      if (stebdry.gt.0.0) then
        stebdry=sqrt(stebdry)
      else
        stebdry=tebdry
      endif
      if (sigtepb.gt.0.0) then
        sigtepb=sqrt(sigtepb)
      else
        sigtepb=tepbry
      endif
!
      return
      end subroutine gette

!**********************************************************************
!>
!!    GETTION gets the ion temperature profile.
!!    
!!
!!    @param kerror :
!!
!**********************************************************************
      subroutine gettion(kerror)
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      dimension pds(6),bwork(ndata),cwork(ndata),dwork(ndata)
      integer*4, intent(inout) :: kerror

      kerror = 0

      if (rion(2).lt.0.0) then
        do i=1,nption
          xsiion(i)=-rion(i)
          sigti(i)=sgtimin*tionex(i)
        enddo
        call zpline(nption,xsiion,tionex,bwork,cwork,dwork)
        do i=1,npress
          xn=-rpress(i)
          tithom(i)=seval(nption,xn,xsiion,tionex,bwork,cwork,dwork)
          stitho(i)=sgtimin*tithom(i)
        enddo
        tibdry=seval(nption,x111,xsiion,tionex,bwork,cwork,dwork)
        tipbry=speval(nption,x111,xsiion,tionex,bwork,cwork,dwork)
        stibdry=sgtimin*tibdry
        sigtipb=sgtimin*tipbry
        return
      endif
!----------------------------------------------------------------------
!--   singular decomposition                                         --
!----------------------------------------------------------------------
      do i=1,nption
        call seva2d(bkx,lkx,bky,lky,c,rion(i),zion(i),pds,ier,n111)
        xsiion(i)=(simag-pds(1))/sidif
      enddo
      need=nption+1
      xsiion(need)=1.0
      tionex(need)=tebdry
      sigti(need)=stebdry
!
      do nj=1,need
        do nk=1,nptionf
          xn=xsiion(nj)
          arsp_cw2(nj,nk)=xn**(nk-1)/sigti(nj)
        enddo
        bdata_cw2(nj)=tionex(nj)/sigti(nj)
      enddo
!---------------------------------------------------------------------
!--   form error matrix                                             --
!---------------------------------------------------------------------
      do i=1,nptionf
        do j=1,nptionf
          ematrix_cw2(i,j)=0.0
          do k=1,need
            ematrix_cw2(i,j)=ematrix_cw2(i,j)+arsp_cw2(k,i)*arsp_cw2(k,j)
          enddo
        enddo
      enddo
!
      nnn=1
      call sdecm(arsp_cw2,ndata,need,nptionf,bdata_cw2,need,nnn,wrsp_cw2,work_cw2,ier)
      if (ier.eq.129) then
        kerror = 1
        call errctrl_msg('gettion','sdecm failed to converge')
        return
      end if

      cond=ier
      toler=1.0e-06_dp*wrsp_cw2(1)
      do i=1,nptionf
        t=0.0
        if (wrsp_cw2(i).gt.toler) t=bdata_cw2(i)/wrsp_cw2(i)
        work_cw2(i)=t
      enddo
      do i=1,nptionf
        tifit(i)=0.0
        do j=1,nptionf
          tifit(i)=tifit(i)+arsp_cw2(i,j)*work_cw2(j)
        enddo
      enddo
!------------------------------------------------------------------
!--   compute chi square                                         --
!------------------------------------------------------------------
      chisqti=0.0
      do i=1,need
        tinow=0.0
        do j=1,nptionf
          tinow=tinow+tifit(j)*xsiion(i)**(j-1)
        enddo
        chisqti=chisqti+((tinow-tionex(i))/sigti(i))**2
      enddo
!---------------------------------------------------------------------
!--   get inverse of error matrix                                   --
!---------------------------------------------------------------------
      call linv1f(ematrix_cw2,nptionf,nppcur,einv_cw2,n444,work_cw2,ier)
!---------------------------------------------------------------------
!--   project ion temperature into Thompson flux grid               --
!---------------------------------------------------------------------
      do i=1,npress
        tithom(i)=0.0
        stitho(i)=0.0
        xn=-rpress(i)
        do j=1,nptionf
          do k=1,nptionf
            stitho(i)=stitho(i)+einv_cw2(k,j)*xn**(j-1)*xn**(k-1)
          enddo
          tithom(i)=tithom(i)+tifit(j)*xn**(j-1)
        enddo
        if (stitho(i).gt.0.0) then
          stitho(i)=sqrt(stitho(i))
        else
          stitho(i)=0.5_dp*tithom(i)
        endif
      enddo
!----------------------------------------------------------------------
!--   boundary values                                                 --
!----------------------------------------------------------------------
      tibdry=0.0
      stibdry=0.0
      tipbry=0.0
      sigtipb=0.0
      do j=1,nptionf
        do i=1,nptionf
          stibdry=stibdry+einv_cw2(i,j)
          sigtipb=sigtipb+(i-1)*(j-1)*einv_cw2(i,j)
        enddo
        tipbry=tipbry+(j-1)*tifit(j)
        tibdry=tibdry+tifit(j)
      enddo
      if (stibdry.gt.0.0) then
        stibdry=sqrt(stibdry)
      else
        stibdry=tibdry
      endif
      if (sigtipb.gt.0.0) then
        sigtipb=sqrt(sigtipb)
      else
        sigtipb=tipbry
      endif
!
      return
      end subroutine gettion

!**********************************************************************
!>
!!    erpote computes the stream function for the
!!    radial electric field. eradial computes the
!!    radial electric field.
!!    
!!
!!    @param ypsi :
!!
!!    @param nnn :
!!
!**********************************************************************
      function erpote(ypsi,nnn)
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      dimension xpsii(nercur)
      data init/0/
!
      if (abs(ypsi).gt.1.0) then
        erpote=0.0
        return
      endif
      call seter(ypsi,xpsii)
      erpote=sum(cerer(1:nnn)*xpsii(1:nnn))
      return
!
      entry erppote(ypsi,nnn)
      if (abs(ypsi).gt.1.0) then
        erppote=0.0
        return
      endif
      call seterp(ypsi,xpsii)
      erppote=-sum(cerer(1:nnn)*xpsii(1:nnn))/sidif
      return
      end function erpote

!**********************************************************************
!>
!!    eradial computes the radial electric field.
!!    
!!
!!    @param ypsi :
!!
!!    @param nnn :
!!
!!    @param reee :
!!
!!    @param zeee :
!!
!**********************************************************************
      function eradial(ypsi,nnn,reee,zeee)
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      dimension pds(6)
      data init/0/
!
      if (abs(ypsi).ge.1.0) then
        eradial=0.0
        return
      endif
      call seva2d(bkx,lkx,bky,lky,c,reee,zeee,pds,ier,n333)
      eradial=erpote(ypsi,nnn)
      eradial=-pds(2)*eradial
      return
!
      entry esradial(ypsi,nnn,reee,zeee)
      if (abs(ypsi).ge.1.0) then
        esradial=0.0
        return
      endif
      call seva2d(bkx,lkx,bky,lky,c,reee,zeee,pds,ier,n333)
      esradial=erppote(ypsi,nnn)
      esradial=-(reee*pds(2))**2*esradial
      fnow=seval(nw,ypsi,sigrid,fpol,bbfpol,ccfpol,ddfpol)
      bbnow=sqrt(fnow**2+pds(2)**2+pds(3)**2)/reee
      esradial=esradial/bbnow
      return
      end function eradial

!**********************************************************************
!>
!!    fpcurr computes the radial derivative
!!    of the poloidal current ff. ffcurr computes
!!    the poloidal current F=twopi RBt/mu0
!!    
!!
!!    @param upsi :
!!
!!    @param nnn :
!!
!**********************************************************************
      function fpcurr(upsi,nnn)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      dimension xpsii(nffcur)
      !$omp declare target 

      if (icutfp.gt.0) then
        ypsi=upsi*xpsimin
      else
        ypsi=upsi
      endif
      if (abs(ypsi).gt.1.0) then
        fpcurr=0.0
        return
      endif
      if (npsi_ext > 0) then
        !fpcurr = seval(npsi_ext,ypsi,psin_ext,ffprim_ext,bfp_ext,cfp_ext,dfp_ext)
        fpcurr = fpcurr * cratiof_ext
        return
      endif
      fpcurr=0.0
      !call setfp(ypsi,xpsii)
      do iiij=nbase+1,nbase+nnn
        iijj=iiij-nbase
        fpcurr=fpcurr+brsp(iiij)*bsffel(kfffnc,iijj,ypsi)
      enddo
!----------------------------------------------------------------------
!--   edge hyperbolic tangent component                              --
!----------------------------------------------------------------------
      if (kedgef.eq.0) return
      siedge=(ypsi-fe_psin)/fe_width
      f0edge=f2edge
      f0back=f0edge/fe_width/sidif
      fpcurr=fpcurr+f0back/cosh(siedge)**2
      return
      end function fpcurr
!
      function fpecrr(upsi,nnn)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      dimension xpsii(nffcur)

      !entry fpecrr(upsi,nnn)
      if (icutfp.gt.0) then
        ypsi=upsi*xpsimin
      else
        ypsi=upsi
      endif
      if (abs(ypsi).gt.1.0) then
        fpecrr=0.0
        return
      endif
      fpecrr=0.0
      call setfp(ypsi,xpsii)
      do iiij=nbase+nnn,nbase+nnn
        iijj=iiij-nbase
        fpecrr=fpecrr+brsp(iiij)*xpsii(iijj)
      enddo
      return
      end function fpecrr
!
      function ffcurr(upsi,nnn)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      dimension xpsii(nffcur)

      !entry ffcurr(upsi,nnn)
      if (icutfp.gt.0) then
        ypsi=upsi*xpsimin
        xsidif=-sidif/xpsimin
      else
        ypsi=upsi
        xsidif=-sidif
      endif
      if (abs(ypsi).ge.1.0) then
        ffcurr=fbrdy
        return
      endif
      ffcurr=0.0
      call setff(ypsi,xpsii)
      do i=nbase+1,nbase+nnn
        nn=i-nbase
        ffcurr=ffcurr+brsp(i)*xpsii(nn)
      enddo
      fb22=fbrdy**2
      ff22=fb22+xsidif*ffcurr/constf2
      if (ff22.lt.0.0) ff22=fb22
      ffcurr=sqrt(ff22)*fbrdy/abs(fbrdy)
!----------------------------------------------------------------------
!--   edge hyperbolic tangent component                              --
!----------------------------------------------------------------------
      if (kedgef.eq.0) return
      siedge=(ypsi-fe_psin)/fe_width
      f0edge=f2edge/constf2
      ff22=ff22+f0edge*(tfedge-tanh(siedge))
      if (ffcurr.lt.0.0) ffcurr=fb22
      ffcurr=sqrt(ff22)*fbrdy/abs(fbrdy)
      return
      end function ffcurr
