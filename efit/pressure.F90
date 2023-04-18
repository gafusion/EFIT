!**********************************************************************
!>
!!    presur computes the relevant parameters for pressure
!!    profile fitting.
!!    
!!    @param jtime : time index
!!    @param niter : iteration number
!!    @param kerror : error flag
!!
!**********************************************************************
      subroutine presur(jtime,niter,kerror)
      use commonblocks,only: c,bkx,bky
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      real*8 ppcurr
      integer*4, intent(in) :: jtime,niter
      integer*4, intent(out) :: kerror
      integer*4 i,ier,m
      real*8 xn,siedge,sgggmin
      real*8 pds(6),xnsi(nppcur),xnsp(nppcur)
      character*50 edatname
      namelist/edat/npress,rpress,zpress,pressr,sigpre

      kerror=0
      kdofit=1
      input_pressure: if (kprfit.ne.2) then
!---------------------------------------------------------------------
!--   input pressure profile                                        --
!---------------------------------------------------------------------
      do i=1,npress
        xn=-rpress(i)
        if (rpress(i).gt.0.0) then
          call seva2d(bkx,lkx,bky,lky,c,rpress(i),zpress(i),pds,ier,n333)
          xn=(simag-pds(1))/sidif
        endif
        call setpr(xn,xnsi)
        call setpp(xn,xnsp)
        do m=1,kppcur
          rprepc(i,m)=-sidif*xnsi(m)/darea
        enddo
!----------------------------------------------------------------------
!--     response for hyperbolic tangent component                    --
!----------------------------------------------------------------------
        if (kedgep.gt.0) then
          siedge=(xn-pe_psin)/pe_width
          rprepe(i)=(tpedge-tanh(siedge))/darea
        endif
!----------------------------------------------------------------------
!--     response for DELTAZ                                          --
!----------------------------------------------------------------------
        if (fitdelz.and.niter.ge.ndelzon) then
          if (rpress(i).le.0.0) then
          rpredz(i)=0.0
          else
          rpredz(i)=pds(3)*ppcurr(xn,kppcur)/darea
          endif
        endif
      enddo
      return
      endif input_pressure
!----------------------------------------------------------------
!--   construct pressure from kinetic data                     --
!--   use psi grid defined by Thompson data                    --
!----------------------------------------------------------------
      npress=npteth
      do i=1,npteth
        call seva2d(bkx,lkx,bky,lky,c,rteth(i),zteth(i),pds,ier,n111)
        xn=(simag-pds(1))/sidif
        if (xn.ge.1.0) then
          npress=i-1
          exit
        endif
        rpress(i)=-xn
        call setpr(xn,xnsi)
        call setpp(xn,xnsp)
        do m=1,kppcur
          rprepc(i,m)=-sidif*xnsi(m)/darea
        enddo
!----------------------------------------------------------------------
!--     response for hyperbolic tangent component                    --
!----------------------------------------------------------------------
        if (kedgep.gt.0) then
          siedge=(xn-pe_psin)/pe_width
          rprepe(i)=(tpedge-tanh(siedge))/darea
        endif
      enddo
!---------------------------------------------------------------
!--   npress=0 can cause problems                             --
!---------------------------------------------------------------
      if (npress.le.0) then
        kdofit=0
        return
      endif
!----------------------------------------------------------------
!--   get ion temperature and density profile                  --
!----------------------------------------------------------------
      if (nptef.ne.0) then
        call gette(kerror)
        if (kerror.gt.0) return
      endif

      call getne(jtime,kerror)
      if (kerror.gt.0) return

      call gettion(kerror)
      if (kerror.gt.0) return

      if (nbeam.ne.0) call getbeam
!----------------------------------------------------------------
!--   construct pressure profile                               --
!----------------------------------------------------------------
      pressb=1.602e+03_dp*(dibdry*tibdry+debdry*tebdry)+pbeamb
      prespb=1.602e+03_dp*(dipbry*tibdry+dibdry*tipbry &
                          +depbry*tebdry+debdry*tepbry)+pbimpb
      sigpreb=(dibdry**2*stibdry**2+sdibdry**2*tibdry**2 &
              +debdry**2*stebdry**2+sdebdry**2*tebdry**2)
      sigpreb=1.602e+03_dp*sqrt(sigpreb)
      sigppb=dibdry**2*sigtipb**2+sdibdry**2*tipbry**2 &
            +dipbry**2*stibdry**2+sigdipb**2*tibdry**2 &
            +debdry**2*sigtepb**2+sdebdry**2*tepbry**2 &
            +depbry**2*stebdry**2+sigdepb**2*tebdry**2
      sigppb=1.602e+03_dp*sqrt(sigppb)
      do i=1,npress
        pressr(i)=1.602e+03_dp*(dnitho(i)*tithom(i)+dnethom(i)*tethom(i)) &
                  +pbimth(i)
        sigpre(i)=(snitho(i)**2*tithom(i)**2+dnitho (i)**2*stitho(i)**2 &
                  +sgneth(i)**2*tethom(i)**2+dnethom(i)**2*sgteth(i)**2)
        sigpre(i)=1.602e+03_dp*sqrt(sigpre(i))
      enddo
      sgggmin=sgprmin
      if(sgprmin.lt.0.0) sgggmin=abs(sgprmin)*pressr(1)
      do i=1,npress
        sigpre(i)=max(sigpre(i),sgggmin)
      enddo
      if (kpressb.eq.1) then
        npress=npress+1
        pressr(npress)=pressbi
        sigpre(npress)=sigprebi
        rpress(npress)=-1.0
      endif
      if (ndokin.ge.100) then
        call setfnmeq(itimeu,'k',ishot,itime,edatname)
        edatname='edat_'//edatname(2:7)// &
                       '_'//edatname(9:13)//'.pressure'
        open(unit=nin,status='old',file=edatname,iostat=ier)
        if (ier.eq.0) close(unit=nin,status='delete')
        open(unit=nin,status='new',file=edatname,delim='quote')
        write (nin,edat)
        close(unit=nin)
      endif
      return
      end subroutine presur

!**********************************************************************
!>
!!    ppcurr computes the radial derivative of the
!!    pressure.
!!
!!    @param ypsi : flux location to produce result
!!    @param nnn : number of coefficients used in representation
!!
!**********************************************************************
      real*8 function ppcurr(ypsi,nnn)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit none
      real*8 prcurr,seval
      integer*4, intent(in) :: nnn
      real*8, intent(in) :: ypsi
      integer*4 i,iijj,iiij,nn
      real*8 p0back,siedge
      real*8 xpsii(nffcur)
!
      if (abs(ypsi).gt.1.0) then
        ppcurr=0.0
        return
      endif
      if (npsi_ext > 0) then
        ppcurr = seval(npsi_ext,ypsi,psin_ext,pprime_ext,bpp_ext,cpp_ext,dpp_ext)
        ppcurr = ppcurr * cratiop_ext
        return
      endif
      ppcurr=0.0
      call setpp(ypsi,xpsii)
      do iiij=nfsum+1,nnn+nfsum
        iijj=iiij-nfsum
        ppcurr=ppcurr+brsp(iiij)*xpsii(iijj)
      enddo
!----------------------------------------------------------------------
!--   edge hyperbolic tangent component                              --
!----------------------------------------------------------------------
      if (kedgep.eq.0) return
      siedge=(ypsi-pe_psin)/pe_width
      p0back=pedge/pe_width/sidif
      ppcurr=ppcurr+p0back/cosh(siedge)**2
      return
!
      entry prcurr(ypsi,nnn)
      if (abs(ypsi).gt.1.0) then
        prcurr=0.0
        return
      endif
      prcurr=0.0
      call setpr(ypsi,xpsii)
      do i=nfsum+1,nfsum+nnn
        nn=i-nfsum
        prcurr=prcurr+brsp(i)*xpsii(nn)
      enddo
      prcurr=-sidif*prcurr/darea+prbdry
!----------------------------------------------------------------------
!--   edge hyperbolic tangent component                              --
!----------------------------------------------------------------------
      if (kedgep.eq.0) return
      siedge=(ypsi-pe_psin)/pe_width
      prcurr=prcurr+pedge/darea*(tpedge-tanh(siedge))
      return
      end function ppcurr

!**********************************************************************
!>
!!    prcur4 computes the plasma pressure by integration.
!!
!!    @param n1set :
!!    @param ypsi : flux location to produce result
!!
!**********************************************************************
      real*8 function prcur4(n1set,ypsi)
      use commonblocks,only: sifpre,bwpre,cwpre,dwpre,sfpre,sprep
      include 'eparm.inc'
      include 'modules1.inc'
      implicit none
      real*8 ppcur4,seval
      integer*4, intent(in) :: n1set
      real*8, intent(in) :: ypsi
      integer*4 i,mw
      real*8 delsi,siii

      if (abs(ypsi).gt.1.0) then
        prcur4=0.0
        return
      endif
!
      mw=nw
      if (n1set.gt.0) then
        do i=2,nw-1
         siii=real(i-1,dp)/(nw-1)
         sifpre(i)=siii
        enddo
        sifpre(1)=0.0
        sifpre(nw)=1.0
        do i=1,nw
          sprep(i)=ppcur4(sifpre(i),kppcur)/darea
        enddo
!
        sfpre(nw)=prbdry
        delsi=sidif/(nw-1)
        do i=1,nw-1
          sfpre(nw-i)=sfpre(nw-i+1)+0.5_dp*(sprep(nw-i+1)+sprep(nw-i))*delsi
        enddo
!
        call zpline(mw,sifpre,sfpre,bwpre,cwpre,dwpre)
      endif
      prcur4=seval(mw,ypsi,sifpre,sfpre,bwpre,cwpre,dwpre)
      return
      end function prcur4

!**********************************************************************
!>
!!    ppcur4 computes the radial derivative of the
!!    pressure based on icurrt.
!!    
!!    @param ypsi : flux location to produce result
!!    @param nnn : number of coefficients used in representation
!!
!**********************************************************************
      real*8 function ppcur4(ypsi,nnn)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit none
      real*8 ppcurr,seval
      integer*4, intent(in) :: nnn
      real*8, intent(in) :: ypsi

      if (abs(ypsi).gt.1.0) then
        ppcur4=0.0
        return
      endif
      select case (icurrt)
      case (1)
        ppcur4=cratio*sbeta/srma
      case (4)
        ppcur4=((1.-ypsi**enp)**emp*(1.-gammap)+gammap)*cratio/rzero
      case default
        ppcur4=ppcurr(ypsi,nnn)
      end select
      return
      end function ppcur4

!**********************************************************************
!>
!!    presurw computes the relevant parameters for rotational
!!    pressure profile fitting.
!!    
!!    @param jtime : time index
!!    @param niter : iteration number
!!
!**********************************************************************
      subroutine presurw(jtime,niter)
      use commonblocks,only: c,bkx,bky
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      real*8 seval,pwpcur
      integer*4, intent(in) :: jtime,niter
      integer*4 init,i,ier,m
      real*8 pds(6),xn,xbaseb,rnow,znow,dmnow
      real*8,dimension(:),allocatable ::ybase,ybaseb
      real*8,dimension(:),allocatable :: bwork,cwork,dwork
      data init/0/

      allocate(ybase(nwwcur),ybaseb(nwwcur),bwork(nmass),cwork(nmass),dwork(nmass))
!
      kdofit=1
      
!--   construct rotational pressure from kinetic data has not been setup
      if (kprfit.eq.4) return
      if (nmass*nomegat.le.0) then
!------------------------------------------------------------------------
!--     specify the rotational pressure profile directly               --
!------------------------------------------------------------------------
        do i=1,npresw
          xn=-rpresw(i)
          if (rpresw(i).gt.0.0) then
            call seva2d(bkx,lkx,bky,lky,c,rpresw(i),zpresw(i),pds,ier,n333)
            xn=(simag-pds(1))/sidif
          endif
          rpresws(i)=xn
          call setpw(xn,ybase)
          call setpwp(xn,ybaseb)
          xbaseb=ybaseb(kwwcur)*xn**2
          do m=1,kwwcur
            rprwpc(i,m)=-sidif*ybase(m)/darea
          enddo
!----------------------------------------------------------------------
!--       response for DELTAZ                                        --
!----------------------------------------------------------------------
          if (fitdelz.and.niter.ge.ndelzon) then
            if (rpresw(i).le.0.0) then
              rprwdz(i)=0.0
            else
              rprwdz(i)=pds(3)*pwpcur(xn,kwwcur)/darea
            endif
          endif
        enddo
        return
      endif
!------------------------------------------------------------------------
!--   form rotational pressure from mass density and rotaional frequency
!------------------------------------------------------------------------
      if (init.eq.0) then
!------------------------------------------------------------------------
!--     set up interpolation                                           --
!------------------------------------------------------------------------
        call zpline(nmass,sibeam,dmass,bwork,cwork,dwork)
        init=1
      endif
      do i=1,nomegat
        xn=-romegat(i)
        if (romegat(i).gt.0.0) then
          rnow=romegat(i)
          znow=zomegat(i)
          call seva2d(bkx,lkx,bky,lky,c,rnow,znow,pds,ier,n333)
          xn=(simag-pds(1))/sidif
        endif
        rpresw(i)=-xn
        rpresws(i)=xn
        dmnow=seval(nmass,xn,sibeam,dmass,bwork,cwork,dwork)
        
        presw(i)=dmnow*omegat(i)*rvtor**2
        sigprw(i)=abs(presw(i))*sigome(i)
        presw(i)=0.5_dp*presw(i)*omegat(i)
        call setpw(xn,ybase)
        call setpwp(xn,ybaseb)
        xbaseb=ybaseb(kwwcur)*xn**2
        do m=1,kwwcur
          rprwpc(i,m)=-sidif*ybase(m)/darea
        enddo
!----------------------------------------------------------------------
!--     response for DELTAZ                                          --
!----------------------------------------------------------------------
        if (fitdelz.and.niter.ge.ndelzon) then
          if (romegat(i).le.0.0) then
            rprwdz(i)=0.0
          else
            rprwdz(i)=pds(3)*pwpcur(xn,kwwcur)/darea
          endif
        endif
      enddo
      npresw=nomegat
!
      return
      end subroutine presurw

!**********************************************************************
!>
!!    pwpcur computes the radial derivative of the
!!    rotational pressure
!!
!!    @param ypsi : flux location to produce result
!!    @param nnn : number of coefficients used in representation
!!
!**********************************************************************
      real*8 function pwpcur(ypsi,nnn)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: nnn
      real*8, intent(in) :: ypsi
      integer*4 iiij,iijj
      real*8 xpsii(nwwcur)

      if (abs(ypsi).gt.1.0) then
        pwpcur=0.0
        return
      endif
      pwpcur=0.0
      call setpwp(ypsi,xpsii)
      do iiij=nfnpcr+1,nnn+nfnpcr
        iijj=iiij-nfnpcr
        pwpcur=pwpcur+brsp(iiij)*xpsii(iijj)
      enddo
      return
      end function pwpcur

!**********************************************************************
!>
!!    pwcur4 computes the rotational pressure by integration.
!!    
!!    @param n1set :
!!    @param ypsi : flux location to produce result
!!
!**********************************************************************
      real*8 function pwcur4(n1set,ypsi)
      use commonblocks, only: sifprw,bwprw,cwprw,dwprw,sfprw,sprwp
      include 'eparm.inc'
      include 'modules1.inc'
      implicit none
      real*8 pwpcu4,seval
      integer*4, intent(in) :: n1set
      real*8, intent(in) :: ypsi
      integer*4 i,mw
      real*8 delsi,siii

      if (abs(ypsi).gt.1.0) then
        pwcur4=0.0
        return
      endif
!
      mw=nw
      if (n1set.gt.0) then
        do i=2,nw-1
          siii=real(i-1,dp)/(nw-1)
          sifprw(i)=siii
        enddo
        sifprw(1)=0.0
        sifprw(nw)=1.0
        do i=1,nw
          sprwp(i)=pwpcu4(sifprw(i),kwwcur)/darea
        enddo
!
        sfprw(nw)=preswb
        delsi=sidif/(nw-1)
        do i=1,nw-1
          sfprw(nw-i)=sfprw(nw-i+1)+0.5_dp*(sprwp(nw-i+1)+sprwp(nw-i))*delsi
        enddo
!
        call zpline(mw,sifprw,sfprw,bwprw,cwprw,dwprw)
      endif
      pwcur4=seval(mw,ypsi,sifprw,sfprw,bwprw,cwprw,dwprw)
      return
      end function pwcur4

!**********************************************************************
!>
!!    pwpcu4 computes the radial derivative of the
!!    rotational pressure based on icurrt.
!!
!!    @param ypsi : flux location to produce result
!!    @param nnn : number of coefficients used in representation
!!
!**********************************************************************
      real*8 function pwpcu4(ypsi,nnn)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit none
      real*8 pwpcur
      integer*4, intent(in) :: nnn
      real*8, intent(in) :: ypsi

      if (abs(ypsi).gt.1.0) then
        pwpcu4=0.0
        return
      endif

      if (icurrt.eq.4) then
        pwpcu4=((1.-ypsi**enw)**emw*(1.-gammaw)+gammaw)*rbetaw &
          *cratio/rzero
      elseif (icurrt.eq.1) then
        pwpcu4=sbetaw*cratio/srma
      else
        pwpcu4=pwpcur(ypsi,nnn)
      endif

      return
      end function pwpcu4

!**********************************************************************
!>
!!    pwcurr computes the ???
!!    
!!    @param ypsi : flux location to produce result
!!    @param nnn : number of coefficients used in representation
!!
!**********************************************************************
      real*8 function pwcurr(ypsi,nnn)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: nnn
      real*8, intent(in) :: ypsi
      integer*4 i,nn
      real*8 xpsii(nwwcur)

      if (abs(ypsi).gt.1.0) then
        pwcurr=0.0
        return
      endif
      pwcurr=0.0
      call setpw(ypsi,xpsii)
      do i=nfnpcr+1,nfnpcr+nnn
        nn=i-nfnpcr
        pwcurr=pwcurr+brsp(i)*xpsii(nn)
      enddo
      pwcurr=-sidif*pwcurr/darea+preswb
      return
      end function pwcurr
