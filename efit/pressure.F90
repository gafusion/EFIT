!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          prcur4 computes the plasma pressure by integration.     **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          94/03/01..........first created                         **
!**                                                                  **
!**********************************************************************
      function prcur4(n1set,ypsi,nnn)
      use commonblocks,only: sifpre,bwpre,cwpre,dwpre,sfpre,sprep
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!      include 'ecomdu1.f90'

      if (abs(ypsi).gt.1.0) then
        prcur4=0.0
        return
      endif
!
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
       do 1000 i=1,nw-1
       sfpre(nw-i)=sfpre(nw-i+1)+0.5_dp*(sprep(nw-i+1)+sprep(nw-i))*delsi
 1000  continue
!
       mw=nw
       call zpline(mw,sifpre,sfpre,bwpre,cwpre,dwpre)
      endif
      prcur4=seval(mw,ypsi,sifpre,sfpre,bwpre,cwpre,dwpre)
      return
      end function prcur4
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          ppcur4 computes the radial derivative of the            **
!**          pressure based on icurrt.                               **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          25/07/85..........revised                               **
!**                                                                  **
!**********************************************************************
      function ppcur4(ypsi,nnn)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!      include 'ecomdu1.f90'

      if (abs(ypsi).gt.1.0) then
        ppcur4=0.0
        return
      endif
      if (icurrt.eq.4) go to 2000
      if (icurrt.eq.1) go to 3000
      ppcur4=ppcurr(ypsi,nnn)
      return
!
 2000 continue
      ppcur4= ((1.-ypsi**enp)**emp*(1.-gammap)+gammap)*cratio &
                /rzero
      return
 3000 continue
      ppcur4= cratio*sbeta/srma
      return
      end function ppcur4
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          ppcurr computes the radial derivative of the            **
!**          pressure.                                               **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          25/07/85..........revised                               **
!**          94/03/11..........revised                               **
!**                                                                  **
!**********************************************************************
      function ppcurr(ypsi,nnn)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!      include 'ecomdu1.f90'
      dimension xpsii(nppcur)
!
! jm.s
      !real*8 :: linearinterp
! jm.e
      if (abs(ypsi).gt.1.0) then
        ppcurr=0.0
        return
      endif
! jm.s
      if (npsi_ext > 0) then
!        ppcurr = linearinterp(ypsi,psin_ext,pprime_ext,npsi_ext)
        ppcurr = seval(npsi_ext,ypsi,psin_ext,pprime_ext,bpp_ext,cpp_ext,dpp_ext)
        ppcurr = ppcurr * cratiop_ext
        return
      endif
! jm.e
      ppcurr=0.0
      call setpp(ypsi,xpsii)
      do 1400 iiij=nfcoil+1,nnn+nfcoil
        iijj=iiij-nfcoil
        ppcurr=ppcurr+brsp(iiij)*xpsii(iijj)
 1400 continue
!----------------------------------------------------------------------
!-- edge hyperbolic tangent component                                --
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
      brspp=0.0
      prcurr=0.0
      call setpr(ypsi,xpsii)
      do 1600 i=nfcoil+1,nfcoil+nnn
        nn=i-nfcoil
        prcurr=prcurr+brsp(i)*xpsii(nn)
 1600 continue
      prcurr=-sidif*prcurr/darea+prbdry
!----------------------------------------------------------------------
!-- edge hyperbolic tangent component                                --
!----------------------------------------------------------------------
      if (kedgep.eq.0) return
      siedge=(ypsi-pe_psin)/pe_width
      prcurr=prcurr+pedge/darea*(tpedge-tanh(siedge))
      return
      end function ppcurr
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          presurw computes the relevant parameters for rotational **
!**          pressure profile fitting.                               **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          94/03/11..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine presurw(jtime,niter)
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      dimension pds(6),ybase(nwwcur),ybaseb(nwwcur)
      dimension bwork(ndata),cwork(ndata),dwork(ndata)
      data init/0/
      save init
!
      kdofit=1
        select case (kprfit-2)
        case (1)
          go to 1000
        case (2)
          go to 2000
        end select
 1000 continue
      if (nmass*nomegat.gt.0) go to 1300
!------------------------------------------------------------------------
!--  specify the rotational pressure profile directly                  --
!------------------------------------------------------------------------
      do 1200 i=1,npresw
        xn=-rpresw(i)
        if (rpresw(i).le.0.0) go to 1030
        call seva2d(bkx,lkx,bky,lky,c,rpresw(i),zpresw(i),pds,ier,n333)
        xn=(simag-pds(1))/sidif
 1030   continue
        rpresws(i)=xn
        call setpw(xn,ybase)
        call setpwp(xn,ybaseb)
        xbaseb=ybaseb(kwwcur)*xn**2
        do 1150 m=1,kwwcur
          rprwpc(i,m)=-sidif*ybase(m)/darea
 1150   continue
!----------------------------------------------------------------------
!-- response for DELTAZ                                              --
!----------------------------------------------------------------------
        if (fitdelz.and.niter.ge.ndelzon) then
          if (rpresw(i).le.0.0) then
          rprwdz(i)=0.0
          else
          rprwdz(i)=pds(3)*pwpcur(xn,kwwcur)/darea
          endif
        endif
 1200 continue
      return
 1300 continue
!------------------------------------------------------------------------
!-- form rotational pressure from mass density and rotaional frequency --
!------------------------------------------------------------------------
      if (init.eq.0) then
!------------------------------------------------------------------------
!--  set up interpolation                                              --
!------------------------------------------------------------------------
        call zpline(nmass,sibeam,dmass,bwork,cwork,dwork)
        init=1
      endif
      do i=1,nomegat
        xn=-romegat(i)
        if (romegat(i).le.0.0) go to 1330
        rnow=romegat(i)
        znow=zomegat(i)
        call seva2d(bkx,lkx,bky,lky,c,rnow,znow,pds,ier,n333)
        xn=(simag-pds(1))/sidif
 1330   continue
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
!-- response for DELTAZ                                              --
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
      return
!------------------------------------------------------------------------
!--  construct rotational pressure from kinetic data                   --
!------------------------------------------------------------------------
 2000 continue
      return
      end subroutine presurw
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          presur computes the relevant parameters for pressure    **
!**          profile fitting.                                        **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          09/09/85..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine presur(jtime,niter,kerror)
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      dimension pds(6),xnsi(nppcur),xnsp(nppcur)
      character*50 edatname
      namelist/edat/npress,rpress,zpress,pressr,sigpre
      integer, intent(inout) :: kerror

      kerror = 0
      kdofit=1
      select case (kprfit)
      case (1)
        go to 1000
      case (2)
        go to 2000
      case (3)
        go to 1000
      end select
 1000 continue
!---------------------------------------------------------------------
!--  input pressure profile                                         --
!---------------------------------------------------------------------
      do 1200 i=1,npress
        xn=-rpress(i)
        if (rpress(i).le.0.0) go to 1030
        call seva2d(bkx,lkx,bky,lky,c,rpress(i),zpress(i),pds,ier,n333)
        xn=(simag-pds(1))/sidif
 1030   continue
        call setpr(xn,xnsi)
        call setpp(xn,xnsp)
        do 1150 m=1,kppcur
          rprepc(i,m)=-sidif*xnsi(m)/darea
 1150   continue
!----------------------------------------------------------------------
!-- response for hyperbolic tangent component                        --
!----------------------------------------------------------------------
        if (kedgep.gt.0) then
          siedge=(xn-pe_psin)/pe_width
          rprepe(i)=(tpedge-tanh(siedge))/darea
        endif
!----------------------------------------------------------------------
!-- response for DELTAZ                                              --
!----------------------------------------------------------------------
        if (fitdelz.and.niter.ge.ndelzon) then
          if (rpress(i).le.0.0) then
          rpredz(i)=0.0
          else
          rpredz(i)=pds(3)*ppcurr(xn,kppcur)/darea
          endif
        endif
 1200 continue
      return
!----------------------------------------------------------------
!--  construct pressure from kinetic data                      --
!----------------------------------------------------------------
 2000 continue
!----------------------------------------------------------------
!--  use psi grid defined by Thompson data                     --
!----------------------------------------------------------------
      npress=npteth
      do 2200 i=1,npteth
        call seva2d(bkx,lkx,bky,lky,c,rteth(i),zteth(i),pds,ier,n111)
        xn=(simag-pds(1))/sidif
        if (xn.ge.1.0) then
          npress=i     -1
          go to 2220
        endif
        rpress(i)=-xn
        call setpr(xn,xnsi)
        call setpp(xn,xnsp)
        do 2150 m=1,kppcur
          rprepc(i,m)=-sidif*xnsi(m)/darea
 2150   continue
!----------------------------------------------------------------------
!-- response for hyperbolic tangent component                        --
!----------------------------------------------------------------------
        if (kedgep.gt.0) then
          siedge=(xn-pe_psin)/pe_width
          rprepe(i)=(tpedge-tanh(siedge))/darea
        endif
 2200 continue
 2220 continue
!---------------------------------------------------------------
!--  npress=0 can cause problems                              --
!---------------------------------------------------------------
      if (npress.le.0) then
        kdofit=0
        return
      endif
!----------------------------------------------------------------
!--  get ion temperature and density profile                   --
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
!--  construct pressure profile                                --
!----------------------------------------------------------------
      pressb=1.602e+03_dp*(dibdry*tibdry+debdry*tebdry)+pbeamb
      prespb=1.602e+03_dp*(dipbry*tibdry+dibdry*tipbry &
                       +depbry*tebdry+debdry*tepbry) +pbimpb
      sigpreb=(dibdry**2*stibdry**2+sdibdry**2*tibdry**2 &
              +debdry**2*stebdry**2+sdebdry**2*tebdry**2)
      sigpreb=1.602e+03_dp*sqrt(sigpreb)
      sigppb = dibdry**2*sigtipb**2+sdibdry**2*tipbry**2 &
              +dipbry**2*stibdry**2+sigdipb**2*tibdry**2 &
              +debdry**2*sigtepb**2+sdebdry**2*tepbry**2 &
              +depbry**2*stebdry**2+sigdepb**2*tebdry**2
      sigppb =1.602e+03_dp*sqrt(sigppb)
      do 2500 i=1,npress
        pressr(i)=1.602e+03_dp*(dnitho(i)*tithom(i)+dnethom(i)*tethom(i)) &
                  +pbimth(i)
        sigpre(i)=(snitho(i)**2*tithom(i)**2+dnitho(i)**2*stitho(i)**2 &
                +sgneth(i)**2*tethom(i)**2+dnethom(i)**2*sgteth(i)**2)
        sigpre(i)=1.602e+03_dp*sqrt(sigpre(i))
 2500 continue
      sgggmin=sgprmin
      if (sgprmin.lt.0.0) sgggmin=abs(sgprmin)*pressr(1)
      do 2600 i=1,npress
        sigpre(i)=max(sigpre(i),sgggmin)
 2600 continue
      if (kpressb.eq.1) then
      npress=npress+1
      pressr(npress)=pressbi
      sigpre(npress)=sigprebi
      rpress(npress)=-1.0
      endif
      if (ndokin.ge.100) then
        call getfnmu(itimeu,'k',ishot,itime,edatname)
        edatname='edat_'//edatname(2:7)// &
                       '_'//edatname(9:13)//'.pressure'
        open(unit=nin,status='old',file=edatname,err=12916)
        close(unit=nin,status='delete')
12916 continue
        open(unit=nin,status='new',file=edatname,delim='quote')
        write (nin,edat)
        close(unit=nin)
      endif
      return
      end subroutine presur

!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          pwcur4 computes the rotational pressure by integration. **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          94/02/28..........first created                         **
!**                                                                  **
!**********************************************************************
      function pwcur4(n1set,ypsi,nnn)
      use commonblocks,only: sifprw,bwprw,cwprw,dwprw,sfprw,sprwp
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!      include 'ecomdu1.f90'

      if (abs(ypsi).gt.1.0) then
        pwcur4=0.0
        return
      endif
!
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
       mw=nw
       call zpline(mw,sifprw,sfprw,bwprw,cwprw,dwprw)
      endif
      pwcur4=seval(mw,ypsi,sifprw,sfprw,bwprw,cwprw,dwprw)
      return
      end function pwcur4
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          pwpcu4 computes the radial derivative of the            **
!**          rotational pressure based on icurrt.                    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          94/03/08..........first created                         **
!**                                                                  **
!**********************************************************************
      function pwpcu4(ypsi,nnn)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!      include 'ecomdu1.f90'

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
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          ppcurr computes the radial derivative of the            **
!**          rotational pressure.                                    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          94/03/08..........first created                         **
!**                                                                  **
!**********************************************************************
      function pwpcur(ypsi,nnn)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!      include 'ecomdu1.f90'
      dimension xpsii(nwwcur)

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
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          pwcurr 
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**********************************************************************
      function pwcurr(ypsi,nnn)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      dimension xpsii(nwwcur)
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
