!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          pflux computes the poloidal fluxes on the r-z grid.     **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**          09/03/90..........Lazarus vertical feedback             **
!**          28/03/93..........fix relax                             **
!**                                                                  **
!**********************************************************************
      subroutine pflux(niter,nnin,ntotal,jtime,kerror)
!vas  f90 modifi
      use var_bunemn
      use commonblocks,only: c,wk,copy,bkx,bky,psiold,psipold, &
                             psipp,work
      include 'eparmdud129.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      integer initresult
      dimension pds(6)
      real*8,dimension(:),allocatable :: psikkk,gfbsum
      data initfb/0/,init/0/

      kerror = 0
      ALLOCATE(psikkk(nwnh),gfbsum(nwnh))

      vfeed=(isetfb.ne.0).and.(init.ne.0).and.(niter.gt.2.or.nnin.gt.2)
      if (ivesel.gt.10) return
!----------------------------------------------------------------------------
!--  save flux from current iterations before update                       --
!----------------------------------------------------------------------------
      do 2100 kk=1,nwnh
        psiold(kk)=psi(kk)
        psipold(kk)=psipla(kk)
 2100 continue
!
      if (ibunmn.eq.1) go to 2000
      if ((ibunmn.eq.2).and.(errorm.gt.errcut)) go to 2000
      if (ibunmn.eq.3) go to 2000
      if ((ibunmn.eq.4).and.(errorm.gt.errcut)) go to 2000
!-----------------------------------------------------------------------------
!--  ibunmn=0, and 2,4 when errorm less than errcut                         --
!--  Green's integral method of obtaining flux, can be computationally      --
!--  intensive                                                              --
!-----------------------------------------------------------------------------
      do 1000 i=1,nw
      do 1000 j=1,nh
        kk=(i-1)*nh+j
        psi(kk)=0.0
        do 300 m=1,nfcoil
          psi(kk)=psi(kk)+gridfc(kk,m)*brsp(m)
  300   continue
        if (ivesel.le.0) go to 340
        do 335 m=1,nvesel
          psi(kk)=psi(kk)+gridvs(kk,m)*vcurrt(m)
  335   continue
  340   continue
        if (iecurr.ne.1) go to 400
        do 350 m=1,nesum
          psi(kk)=psi(kk)+gridec(kk,m)*ecurrt(m)
  350   continue
  400   continue
        if (iecurr.ne.2) go to 405
        do  m=1,nesum
          psi(kk)=psi(kk)+gridec(kk,m)*cecurr(m)
        enddo
  405   continue
        if (iacoil.gt.0) then
         do 421 m=1,nacoil
          psi(kk)=psi(kk)+gridac(kk,m)*caccurt(jtime,m)
  421    continue
        endif
        psipla(kk)=psi(kk)
        if (ivacum.gt.0) go to 1000
        do 500 ii=1,nw
        do 500 jj=1,nh
          kkkk=(ii-1)*nh+jj
          mj=abs(j-jj)+1
          mk=(i-1)*nh+mj
          psi(kk)=psi(kk)+gridpc(mk,ii)*pcurrt(kkkk)
  500   continue
        psipla(kk)=psi(kk)-psipla(kk)
 1000 continue
!
      ibunmn=0
!------------------------------------------------------------------------------
!--  relaxation before return                                                --
!------------------------------------------------------------------------------
      go to 6000
 2000 continue

!-----------------------------------------------------------------------------
!-- Buneman's method of obtaining psi at the inner grid points              --
!-- only plasma contribution                                                --
!-----------------------------------------------------------------------------
      if (init.gt.0) go to 2020
      nww=nw-1
      nhh=nh-1
      nbmdim=max(nw,nh)+1
!-----These must be brought into the integrals 
      rgrid1=rgrid(1)
      delrgrid=rgrid(2)-rgrid(1)
      delz=zgrid(2)-zgrid(1)
      drdz2=(delrgrid/delz)**2
!---------------------------------------------------------------------
!--  optional vertical feedback control                             --
!---------------------------------------------------------------------
      if (isetfb.ne.0) then
        if(nw.gt.30) ioffr=ioffr*((nw+1)/33)
        if(nh.gt.30) ioffz=ioffz*((nh+1)/33)
        imd=(nw*nh)/2+1+nh*ioffr+ishiftz
        kct1=imd-ioffz
        kct2=imd+ioffz
        deltaz=zgrid(ioffz+nh/2+1)-zgrid(-ioffz+nh/2+1)
      endif
      init=1
 2020 continue
!----------------------------------------------------------------------------
!--  obtain fluxes by inverting del*                                       --
!----------------------------------------------------------------------------
      if (vfeed) then
        dpsip_last=dpsip
        psic1=psi(kct1)
        psic2=psi(kct2)
        dpsip=psic2-psic1
        psibar=(psic2+psic1)/2
      endif
      if (abs(vcurfb(1)).gt.1.e-6_dp) then
        iinow=vcurfb(3)
        if (ntotal.lt.iinow) then
          if (ntotal.eq.1) then
            do 2027 kk=1,nwnh
              psikkk(kk)=0.0
 2027       continue
          endif
        else
          icurfb=vcurfb(2)
          if (mod(ntotal-iinow,icurfb).eq.0) then
            do 2033 kk=1,nwnh
              psikkk(kk)=psiold(kk)
 2033       continue
          endif
        endif
      endif

!-----------------------------------------------------------------------
!-- boundary terms                                                    --
!-----------------------------------------------------------------------
 2050 continue 
      do 2200 j=1,nh
        kk=(nw-1)*nh+j
        psipp(j)=0.
        psipp(kk)=0.
        do 2170 ii=1,nw
        do 2170 jj=1,nh
          kkkk=(ii-1)*nh+jj
          mj=abs(j-jj)+1
          mk=(nw-1)*nh+mj
          psipp(j)=psipp(j)-gridpc(mj,ii)*pcurrt(kkkk)
          psipp(kk)=psipp(kk)-gridpc(mk,ii)*pcurrt(kkkk)
          psi(j)=psipp(j)
          psi(kk)=psipp(kk)
 2170   continue
 2200 continue
      do 2400 i=2,nw-1
        kk1=(i-1)*nh
        kknh=kk1+nh
        kk1=kk1+1
        psipp(kk1)=0.
        psipp(kknh)=0.
        do 2370 ii=1,nw
        do 2370 jj=1,nh
          kkkk=(ii-1)*nh+jj
          mj1=abs(jj-1)+1
          mjnh=abs(nh-jj)+1
          mk1=(i-1)*nh+mj1
          mknh=(i-1)*nh+mjnh
          psipp(kk1)=psipp(kk1)-gridpc(mk1,ii)*pcurrt(kkkk)
          psipp(kknh)=psipp(kknh)-gridpc(mknh,ii)*pcurrt(kkkk)
          psi(kk1)=psipp(kk1)
          psi(kknh)=psipp(kknh)
 2370   continue
 2400 continue
!-------------------------------------------------------------------------
!-- get flux at inner points by inverting del*, only plasma flux        --
!-- first set up plasma currents, single cyclic method gets factor of 2 --
!-------------------------------------------------------------------------
      if (isolve.eq.0) then
!       original buneman solver method
        do 2600 i=2,nw-1
        do 2600 j=2,nh-1
          kk=(i-1)*nh+j
          psi(kk)=tmu2*pcurrt(kk)*rgrid(i)
 2600   continue
        call buneto(psi,nw,nh,work)
      else 
!       new faster single cyclic reduction method
        do 2700 i=2,nw-1
        do 2700 j=2,nh-1
          kk=(i-1)*nh+j
          psi(kk)=tmu2*pcurrt(kk)*rgrid(i)*2.0
 2700   continue
        call pflux_cycred(psi,work,kerror)
        if (kerror.gt.0) return
      endif
      do 3000 i=1,nwnh
        psi(i)=-psi(i)
 3000 continue

!----------------------------------------------------------------------------
!--  optional symmetrized solution                                         --
!----------------------------------------------------------------------------
      if (symmetrize) then
        difpsi=0.0
        do 3100 i=1,nw
        do 3100 j=1,nh/2
          k1=(i-1)*nh+j
          k2=i*nh-j+1
          val=(psi(k1)+psi(k2))/2
          if (sidif.ne.0.0) then
            difnow=(psi(k1)-psi(k2))/sidif
            difnow=abs(difnow)
            difpsi=max(difnow,difpsi)
          endif
          psi(k1)=val
          psi(k2)=val
3100    continue
      endif ! end symmetrize loop
!------------------------------------------------------------------------
!--  optional damping out the m=1 vertical eigen mode                  --
!------------------------------------------------------------------------
      if (abs(vcurfb(1)).lt.1.e-6_dp) go to 43000
      if (ntotal.lt.5) go to 43000
!-----------------------------------------------------------------------
!-- sum Green's functions for m=1 eigen mode                          --
!-----------------------------------------------------------------------
      if (initfb.eq.0) then
      idelr=nw/8
      idelz=nh/12
      jtop=(nh+1)/2+nh/4-1
      do 42200 i=1,nw
      do 42200 j=1,nh
        kk=(i-1)*nh+j
        gfbsum(kk)=0.0
        do 42170 ii=1+idelr,nw-idelr
        do 42170 jj=jtop-idelz,jtop+idelz
          mj=abs(j-jj)+1
          mk=(nw-1)*nh+mj
          gfbsum(kk)=gfbsum(kk)+gridpc(mj,ii)
42170   continue
42200 continue
      jlow=(nh+1)/2-nh/4+1
      do 42500 i=1,nw
      do 42500 j=1,nh
        kk=(i-1)*nh+j
        do 42400 ii=1+idelr,nw-idelr
        do 42400 jj=jlow-idelz,jlow+idelz
          mj=abs(j-jj)+1
          mk=(nw-1)*nh+mj
          gfbsum(kk)=gfbsum(kk)-gridpc(mj,ii)
42400   continue
42500 continue
      initfb=-1
      endif
!---------------------------------------------------------------------
!--  get damping currents                                           --
!---------------------------------------------------------------------
      tvfbrt(ntotal)=0.0
      do 42600  i=1+idelr,nw-idelr
      do 42600  j=jtop-idelz,jtop+idelz
          kk=(i-1)*nh+j
          tvfbrt(ntotal)=tvfbrt(ntotal)+(psikkk(kk)-psiold(kk))
42600 continue
      do 42620  i=1+idelr,nw-idelr
      do 42620  j=jlow-idelz,jlow+idelz
          kk=(i-1)*nh+j
          tvfbrt(ntotal)=tvfbrt(ntotal)-(psikkk(kk)-psiold(kk))
42620 continue
      if (initfb.eq.-1) then
        vcurfi=vcurfb(1)*cpasma(jtime)/abs(tvfbrt(ntotal))
        initfb=1
      endif
      tvfbrt(ntotal)=tvfbrt(ntotal)*vcurfi
      do 42800 kk=1,nwnh
        psi(kk)=psi(kk)+gfbsum(kk)*tvfbrt(ntotal)
42800 continue
43000 continue
!--------------------------------------------------------------------
!--  optional vertical feedback control                            --
!--    psi(kct1)  !psi at lower point                              --
!--    psi(kct2) ! at upper point                                  --
!--    psilu     !psi (per amp) at lower point due to upper coils  --
!--    brfb(1) is lower current                                    --
!--------------------------------------------------------------------
      if (.not.vfeed) go to 3401
      psiul_psiuu=-grdfdb(kct2,1)
      psilu_psill=grdfdb(kct1,1)
      zdwn_now=0
      zup_now=0
      zcontr=0
      zcurnow=0
      sumc=0
      rcurnow=0
      do 12300 icontr=1,nfound
        zcontr=zcontr+yout(icontr)/nfound
12300 continue
      do 12320 i=1,nw
      do 12320 j=1,nh
        kk=(i-1)*nh+j
        zcurnow=zcurnow+pcurrt(kk)*zgrid(j)
        rcurnow=rcurnow+pcurrt(kk)*rgrid(j)
        sumc=sumc+pcurrt(kk)
12320 continue
      zcurnow=zcurnow/sumc
      rcurnow=rcurnow/sumc
      if (zelip.eq.0.) then
        brfb(1)=-gainp*dpsip/(psiul_psiuu+psilu_psill) &
        -gain*(dpsip-dpsip_last)/(psiul_psiuu+psilu_psill)
      else
        brfb(1)=-gainp* &
        (dpsip-psibar*((zelip-ishiftz*delz)/deltaz)) &
        /(psiul_psiuu+psilu_psill) &
        -gain*(dpsip-dpsip_last)/(psiul_psiuu+psilu_psill)
      endif
      brfb(2)=-brfb(1)
      brfbc(ntotal)=brfb(1)
      if (.not.vfeed) then
          brfb(1)=0.
          brfb(2)=0.
      endif
 3401 continue
!----------------------------------------------------------------------------
!--  add flux from external coils                                          --
!----------------------------------------------------------------------------
      do 3600 kk=1,nwnh
        psipla(kk)=psi(kk)
        if (ivesel.le.0) go to 3140
        do 3135 m=1,nvesel
          psi(kk)=psi(kk)+gridvs(kk,m)*vcurrt(m)
 3135   continue
 3140   continue
        if (iecurr.ne.1) go to 3160
        do 3150 m=1,nesum
          psi(kk)=psi(kk)+gridec(kk,m)*ecurrt(m)
 3150   continue
 3160   continue
        if (iecurr.ne.2) go to 3165
        do  m=1,nesum
          psi(kk)=psi(kk)+gridec(kk,m)*cecurr(m)
        enddo
 3165   continue
        if (vfeed) then
          psi(kk)=psi(kk)+grdfdb(kk,1)*brfb(2)
        endif
        do 3200 m=1,nfcoil
          psi(kk)=psi(kk)+gridfc(kk,m)*brsp(m)
 3200   continue
        if (iacoil.gt.0) then
         do 3535 m=1,nacoil
          psi(kk)=psi(kk)+gridac(kk,m)*caccurt(jtime,m)
 3535    continue
        endif
 3600   continue
 6000 continue
!----------------------------------------------------------------------------
!-- rigid vertical shift correction ?                                      --
!----------------------------------------------------------------------------
      if (ifitdelz.eq.1) then
      if (fitdelz.and.ntotal.ge.ndelzon) then
      cdelznow=cdelz(ntotal-1)/100.
      do i=1,nw
      do j=1,nh
        kk=(i-1)*nh+j
        call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n333)
        psi(kk)=psi(kk)+cdelznow*pds(3)
      enddo
      enddo
      endif
      endif
!----------------------------------------------------------------------------
!--  relax flux if needed                                                  --
!----------------------------------------------------------------------------
      if (ntotal.le.1) go to 7000
      if (abs(relax-1.0).lt.1.0e-03_dp) go to 7000
      do 6500 kk=1,nwnh
        psi(kk)=relax*psi(kk)+(1.-relax)*psiold(kk)
        psipla(kk)=relax*psipla(kk)+(1.-relax)*psipold(kk)
 6500 continue
 7000 continue
!
      DEALLOCATE(psikkk,gfbsum)
!
      return
      end subroutine pflux
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**         residu computes the flux variations on the r-z grid.     **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**          28/03/93..........fixe relax                            **
!**                                                                  **
!**********************************************************************
      subroutine residu(nx,jtime)
      use commonblocks,only: psiold,psipold,psipp
      include 'eparmdud129.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!      include 'ecomdu1.f90'

      if (ivacum.gt.0) return
      errold=errorm
      errave=0.0
      errorm=0.0
      do i=1,nw
        do j=1,nh
          kk=(i-1)*nh+j
          change=abs(psi(kk)-psiold(kk))
          errorm=max(errorm,change)
          errave=errave+change
          if (errorm.le.change) then
            iermax(nx)=i
            jermax(nx)=j
          end if
        end do
      end do

      errorm=errorm/abs(sidif)

      aveerr(nx)=errave/abs(sidif)/nwnh
      cerror(nx)=errorm
      idone=0
      if (errorm.le.error) idone=1
      !----------------------------------------------------------------------
      !--  Turn on vertical stabilization if error small                   --
      !----------------------------------------------------------------------
      if ((errorm.le.errdelz).and.fitdelz) then
        ndelzon = 3
      else
        ndelzon = 999
      endif
      !----------------------------------------------------------------------
      !--  vertical stabilization and iteration information                --
      !----------------------------------------------------------------------
      if (itell.gt.0.and.isetfb.ge.0) then
        if (itell.eq.1) then
          !if (nx.eq.1) write (nttyo,10017) itime, rank
          if (nx.eq.1) write (nttyo,'(x)')
          if (nsol.eq.0) then
            if (mmbmsels.eq.0) then
              write (nttyo,10019) rank,itime,nx,tsaisq(jtime),zmaxis,errorm,delzmm, &
                sum(chigam)
            else
              write (nttyo,90019) rank,itime,nx,tsaisq(jtime),zmaxis,errorm,delzmm, &
                sum(chigam),tchimls
            endif
          else
            write (nttyo,10020) rank,itime,nx,tsaisq(jtime),zmaxis,errorm,delzmm, &
              sum(chigam),erbmax,erbsmax
          endif
        elseif (itell.eq.2) then
          write (nttyo,10021) rank,itime,nx,ali(jtime),abs(betatn),errorm,qsiw(1)
        elseif (itell.eq.3) then
          write (nttyo,10023) rank,itime,nx,difpsi,zmaxis,errorm,delzmm
        elseif (itell.eq.4) then
          if (nx.eq.1) then
            !write (nttyo,10017) itime, rank
            write (nttyo,'(x)')
            if (kstark.gt.0) then
              write (nttyo,80019) rank,itime,nx,tsaisq(jtime),zmaxis,errorm,delzmm, &
                cdelz(1),cdeljsum,sum(chigam)
            else
              write (nttyo,10025) rank,itime,nx,tsaisq(jtime),zmaxis,errorm &
                ,delzmm,cdelz(1),cdeljsum
            endif
          else
            if (kstark.gt.0) then
              write (nttyo,80019) rank,itime,nx,tsaisq(jtime),zmaxis,errorm,delzmm, &
                cdelz(nx-1),cdeljsum,sum(chigam)
            else
              write (nttyo,10025)  rank,itime,nx,tsaisq(jtime),zmaxis,errorm &
                ,delzmm,cdelz(nx-1),cdeljsum
            endif
          endif
        endif
      endif
      if (isetfb.ne.0) then
        write (4,10009) rank,itime,nx,tsaisq(jtime),zmaxis,errorm,delzmm &
          ,brfb(1)
        if (isetfb.lt.0) &
          write (6,10009) rank,itime,nx,tsaisq(jtime),zmaxis,errorm,delzmm &
          ,brfb(1)
      elseif (eelip.gt.2.25_dp .and. itell.eq.0) then
        write (6,10009) rank,itime,nx,tsaisq(jtime),zmaxis,errorm,delzmm &
          ,brfb(1)
      endif
      if (idebug /= 0) then
        write (nttyo,*) 'cratio,cratio_ext,cratiop_ext,cratiof_ext= ', &
          cratio,cratio_ext,cratiop_ext,cratiof_ext
        write (nttyo,*) 'scalepp_ext,scaleffp_ext= ', &
          scalepp_ext,scaleffp_ext
      endif
      return
10009 format (x,'r=',i3,1x,'t=',i6,1x,'iter',i3.3, &
      ' chsq=',1pe8.2,' zmag=',1pe9.2,' err=',1pe8.2,' dz=',1pe10.3, &
      ' Ifb=',1pe9.2)
!10017 format (/,x,' ----- time =',i6,' ms ----- (',i2,')')
10019 format (x,'r=',i3,1x,'t=',i6,1x,'it=',i3, &
      ' chi2=',1pe8.2,' zm=',1pe9.2,' err=',1pe9.3, &
      ' dz=',1pe10.3,' chigam=',1pe9.2)
80019 format (x,'r=',i3,1x,'t=',i6,1x,'it=',i3, &
      ' chi2=',1pe8.2,' zm=',1pe9.2,' err=',1pe9.3, &
      ' dz=',1pe10.3,' delz=',1pe10.3,' dj=',1pe9.3,' chigam=',1pe9.2)
90019 format (x,'r=',i3,1x,'t=',i6,1x,'it=',i3, &
      ' chi2=',1pe8.2,' zm=',1pe9.2,' err=',1pe9.3, &
      ' dz=',1pe10.3,' chigam=',1pe9.2,' chimls=',1pe9.2)
10020 format (x,'r=',i3,1x,'t=',i6,1x,'it=',i3, &
      ' chi2=',1pe8.2,' zm=',1pe9.2,' err=',1pe9.3, &
      ' dz=',1pe10.3,' chigam=',1pe9.2,' errb=',1pe9.2,' errbs=',1pe9.2)
10021 format (x,'r=',i3,1x,'t=',i6,1x,'it=',i3, &
      ' li=',1pe9.3,' betan=',1pe9.3,' err=',1pe9.3, &
      ' qs=',1pe9.3)
10023 format (x,'r=',i3,1x,'t=',i6,1x,'it=',i3, &
      ' dpsi=',1pe10.3,' zm=',1pe9.2,' err=',1pe9.3, &
      ' dz=',1pe10.3)
10025 format (x,'r=',i3,1x,'t=',i6,1x,'it=',i3, &
      ' chi2=',1pe8.2,' zm=',1pe9.2,' err=',1pe9.3, &
      ' dz=',1pe10.3,' delz=',1pe10.3,' dj=',1pe9.3)
      end subroutine residu
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          step computes the dimensionless poloidal fluxes for     **
!**          the r-z grid.                                           **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**                                                                  **
!**********************************************************************
      subroutine steps(ix,ixt,ixout,jtime,kerror)
      use commonblocks,only: c,wk,copy,bkx,bky,zeros,xouts,youts, &
           rsplt,zsplt,csplt
      include 'eparmdud129.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      common/wwork1/xlims(5),ylims(5),limtrs,xlmins
      dimension pds(6)
      integer iii
      real :: zmaxis_last = 0.0
      data isplit/8/,psitol/1.0e-04_dp/
      save xguess, yguess, xltrac, radbou
!
      if (ivacum.gt.0) return
      if (ixt.gt.1) go to 100
      xguess=(rgrid(1)+rgrid(nw))/2.
      yguess=(zgrid(1)+zgrid(nh))/2.
      if (zbound.ne.0.0) yguess=zbound
      if (rbound.ne.0.0) xguess=rbound
      xltrac=xlmin
      if (ibound.eq.-1) xltrac=xlmax
      radbou=(xguess+xltrac)/2.
!
  100 continue
!----------------------------------------------------------------------
!-- first set up bi-cubic spline interpolation in findax             --
!----------------------------------------------------------------------
      m10=10

      if (idebug >= 2) write (6,*) 'Entering findax'
      !print *, 'nw,nh,rgrid,zgrid',nw,nh,rgrid,zgrid
      !print *, 'rmaxis,zmaxis,simag', rmaxis,zmaxis,simag
      !print *, 'psibry,rseps(1,jtime),zseps(1,jtime),m10', psibry,rseps(1,jtime),zseps(1,jtime),m10
      !print *, 'xout,yout,nfound,psi', xout,yout,nfound,psi
      !print *, 'xmin,xmax,ymin,ymax',xmin,xmax,ymin,ymax
      !print *,  'zxmin,zxmax,rymin,rymax' , zxmin,zxmax,rymin,rymax
      !print *,  'dpsi,bpol,bpolz' , dpsi,bpol,bpolz
      !print *, 'limitr,xlim,ylim,limfag', limitr,xlim,ylim,limfag
      !print *, 'ixt,jtime,kerror', ixt,jtime,kerror
      call findax(nw,nh,rgrid,zgrid,rmaxis,zmaxis,simag, &
                  psibry,rseps(1,jtime),zseps(1,jtime),m10, &
                  xout,yout,nfound,psi,xmin,xmax,ymin,ymax, &
                  zxmin,zxmax,rymin,rymax,dpsi,bpol,bpolz, &
                  limitr,xlim,ylim,limfag,ixt,jtime,kerror)

      if (kerror.gt.0) return
      if (nsol.gt.0) then
   
          write (6,*) 'STEPS R,Z,Si,Err = ', rsol(1),zsol(1),wsisol,ier
          call seva2d(bkx,lkx,bky,lky,c,rbdry(1),zbdry(1),pds,ier,n111)
          write (6,*) 'STEPS R,Z,Si,Err = ', rbdry(1),zbdry(1),pds(1),ier
          call seva2d(bkx,lkx,bky,lky,c,rbdry(nbdry),zbdry(nbdry),pds,ier &
             ,n111)
          write (6,*) 'STEPS R,Z,Si,Err = ', rbdry(nbdry),zbdry(nbdry) &
            ,pds(1),ier

            write(6,*) 'rsplt(kk),zsplt(kk)',rbdry(nbdry),zbdry(nbdry)
            write(6,*) 'lkx, lky',lkx,lky
            write(6,*) 'pds,ier,n111', pds,ier,n111

        if (idebug >= 2) then

          call seva2d(bkx,lkx,bky,lky,c,xout(1),yout(1),pds,ier &
             ,n111)
          write (6,*) 'STEPS simag,psibry,n1,si = ',simag,psibry,n111,pds(1)
          write (6,*) 'STEPS Z,R = ', zgrid(33),(rgrid(i),i=45,45)
          write (6,*) 'STEPS si = ',(psi((i-1)*65+33),i=45,45)
          call seva2d(bkx,lkx,bky,lky,c,rgrid(45),zgrid(33),pds,ier &
             ,n111)
          write (6,*) 'STEPS R,Z,si = ', rgrid(45),zgrid(33),pds(1)
        endif
      endif
!-----------------------------------------------------------------------
!--  Trace boundary, first check for counter beam injection           --
!-----------------------------------------------------------------------
      if (pasmat(jtime).lt.-1.e3_dp) then
        nnerr=10000
      else
        nnerr=0
      endif
      call bound(psi,nw,nh,nwnh,psibry,xmin,xmax,ymin,ymax, &
                 zero,rgrid,zgrid,xguess,yguess,ixt,limitr,xlim,ylim, &
                 xout,yout,nfound,xltrac,npoint,rymin,rymax,dpsi, &
                 zxmin,zxmax,nnerr,ishot,itime, &
                 limfag,radbou,kbound,tolbndpsi)
      if (nnerr.gt.0) then
        kerror=1
        return
      endif
!----------------------------------------------------------------------
!--  find magnetic axis and poloidal flux at axis simag              --
!----------------------------------------------------------------------
      m20=20
      call findax(nw,nh,rgrid,zgrid,rmaxis,zmaxis,simag, &
                  psibry,rseps(1,jtime),zseps(1,jtime),m20, &
                  xout,yout,nfound,psi,xmin,xmax,ymin,ymax, &
                  zxmin,zxmax,rymin,rymax,dpsi,bpol,bpolz, &
                  limitr,xlim,ylim,limfag,ixt,jtime,kerror)
      if (kerror.gt.0) return
      sidif=simag-psibry
      eouter=(ymax-ymin)/(xmax-xmin)
      zplasm=(ymin+ymax)/2.
      aouter=(xmax-xmin)/2.

!-----------------------------------------------------------------------
!--   force free current in the scrape-off layer                      --
!-----------------------------------------------------------------------
      if (icutfp.eq.2) then
        xvsmaxo=xvsmax
        xvsmin=1.e10_dp
        xvsmax=-1.e10_dp
        if (limvs.eq.0) then
        itot=isplit*isplit
        do 51000 k=1,nvesel
          call splitc(isplit,rsplt,zsplt,csplt, &
                  rvs(k),zvs(k),wvs(k),hvs(k),avs(k),avs2(k),cdum)
          do 50900 kk=2,itot
            call seva2d(bkx,lkx,bky,lky,c,rsplt(kk),zsplt(kk), &
                        pds,ier,n111)
            write(6,*) 'rgrid(i),zgrid(j)',rsplt(kk),zsplt(kk)
            write(6,*) 'lkx, lky',lkx,lky
            write(6,*) 'pds,ier,n111', pds,ier,n111
            xvsmin=min(xvsmin,pds(1))
            xvsmax=max(xvsmax,pds(1))
50900     continue
51000   continue
        else
        do 51009 k=1,limitr-1
           delx=xlim(k+1)-xlim(k)
           dely=ylim(k+1)-ylim(k)
           dels=sqrt(delx**2+dely**2)
           nn=dels/0.002_dp
           nn=max(5,nn)
           delx=delx/(nn-1)
           dely=dely/(nn-1)
           do 51007 kk=2,nn
            xww=xlim(k)+delx *(kk-1)
            yww=ylim(k)+dely *(kk-1)
            call seva2d(bkx,lkx,bky,lky,c,xww,yww,pds,ier,n111)
            xvsmin=min(xvsmin,pds(1))
            xvsmax=max(xvsmax,pds(1))
51007      continue
51009   continue
        endif
!--------------------------------------------------------------
!--  exclude private region flux                             --
!--------------------------------------------------------------
        xvsmax=psibry
!--------------------------------------------------------------
!--  possible second separatrix                              --
!--------------------------------------------------------------
        rsepex=-999.
        yvs2=1000.
        if (kskipvs.eq.0) go to 10000
        avebp=cpasma(jtime)*tmu/aouter
        bpmin=avebp
        sibpold=sibpmin
        do 51100 i=1,nw
        do 51090 j=1,nh
          kk=(i-1)*nh+j
          if (zero(kk).gt.0.0005_dp.and.www(kk).lt.0.1_dp) then
           call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n333)
           write(6,*) 'rgrid(i),zgrid(j)',rgrid(i),zgrid(j)
           write(6,*) 'lkx, lky',lkx,lky
           write(6,*) 'pds,ier,n111', pds,ier,n333
           bpnow=sqrt(pds(2)**2+pds(3)**2)/rgrid(i)
           if (bpnow.le.bpmin) then
           if ((abs(dpsi).le.psitol).or.((abs(dpsi).gt.psitol).and. &
               (zgrid(j)*zseps(1,jtime).lt.0.0))) then
              bpmin=bpnow
              xs=rgrid(i)
              ys=zgrid(j)
              sibpmin=pds(1)
           endif
           endif
!          endif
          endif
51090   continue
51100   continue

        if (bpmin.eq.avebp) go to 9320
        relsi=abs((sibpmin-psibry)/sidif)
        if (bpmin.le.0.10_dp*avebp.and.relsi.gt.0.005_dp) then
!------------------------------------------------------------------
!-- find second separatrix                                       --
!------------------------------------------------------------------
          do 9300 j=1,40
          call seva2d(bkx,lkx,bky,lky,c,xs,ys,pds,ier,n666)
          write(6,*) 'xs,ys',xs, ys
          write(6,*) 'lkx, lky',lkx,lky
          write(6,*) 'pds,ier,n111', pds,ier,n111

          det=pds(5)*pds(6)-pds(4)*pds(4)
          if (abs(det).lt.1.0e-15_dp) go to 9305
          xerr=(-pds(2)*pds(6)+pds(4)*pds(3))/det
          yerr=(-pds(5)*pds(3)+pds(2)*pds(4))/det
          xs=xs+xerr
          ys=ys+yerr
          if (xerr*xerr+yerr*yerr.lt.1.0e-12_dp) go to 9310
 9300     continue
 9305     continue
          epssep=xerr*xerr+yerr*yerr
          write (nttyo,11001) epssep,ixt
          if (iand(iout,1).ne.0) write (nout,11001) epssep,ixt
          if (epssep.lt.1.0e-10_dp) go to 9310
          go to 9320
 9310     continue
         sibpmin=pds(1)
         yvs2=ys
         rsepex=xs
         relsi=abs((sibpmin-psibry)/sidif)
         if (relsi.gt.0.005_dp) then
         if (ixt.gt.1) sibpmin=sibpmin*(1.-vsdamp)+sibpold*vsdamp
         xvsmin=max(xvsmin,sibpmin)
         endif
        endif
10000   continue
 9320   continue
        if (alphafp.ge.0.0) then
           xpsimin=xvsmin+alphafp*(xvsmax-xvsmin)
        else
           xpsimino=xpsimins
           xpsimin=abs(alphafp)*xvsmax
           xpsimins=xpsimin
           if (ixt.gt.1) xpsimin=xpsimin*(1.-vsdamp)+xpsimino*vsdamp
           alphamu=(xpsimin-xvsmin)/(xvsmax-xvsmin)
        endif
        xpsialp=xpsimin
        xpsimin=sidif/(simag-xpsimin)
      endif 
!-----------------------------------------------------------------------
!-- get normalized flux function XPSI                                 --
!-----------------------------------------------------------------------
      do 1000 i=1,nw
      do 1000 j=1,nh
        kk=(i-1)*nh+j
        if (icutfp.eq.0) then
          xpsi(kk)=1.1_dp
          if ((rgrid(i).lt.xmin).or.(rgrid(i).gt.xmax)) go to 1000
          if ((zgrid(j).lt.ymin).or.(zgrid(j).gt.ymax)) go to 1000
          xpsi(kk)=(simag-psi(kk))/sidif
        else
          if (zero(kk).gt.0.0005_dp) then
            xpsi(kk)=(simag-psi(kk))/sidif
            if (xpsi(kk)*xpsimin.le.1.0.and.xpsi(kk)*xpsimin.ge.0.0) then
              if ((rgrid(i).lt.rminvs).or.(rgrid(i).gt.rmaxvs)) &
                   xpsi(kk)=1000.
              if ((zgrid(j).lt.zminvs).or.(zgrid(j).gt.zmaxvs)) then
                   xpsi(kk)=1000.
              endif
              if (xpsi(kk).lt.1.0.and.zgrid(j).lt.ymin) xpsi(kk)=1000.
              if (xpsi(kk).lt.1.0.and.zgrid(j).gt.ymax) xpsi(kk)=1000.
              if (abs(zgrid(j)).gt.abs(yvs2).and. &
                zgrid(j)*yvs2.gt.0.0) xpsi(kk)=1000.
            endif
          else
            xpsi(kk)=1000.
          endif
        endif
 1000 continue
!-----------------------------------------------------------------------
!-- get SOL flux if needed                                            --
!-----------------------------------------------------------------------
      if (nsol.gt.0) then
        call seva2d(bkx,lkx,bky,lky,c,rsol(1),zsol(1), &
                     pds,ier,n111)
        wsisol=pds(1)
        if (idebug >= 2) then
          write (6,*) 'STEPS R,Z,Si,Err = ', rsol(1),zsol(1),wsisol,ier 
          call seva2d(bkx,lkx,bky,lky,c,rbdry(1),zbdry(1),pds,ier,n111)
          write (6,*) 'STEPS R,Z,Si,Err = ', rbdry(1),zbdry(1),pds(1),ier
          call seva2d(bkx,lkx,bky,lky,c,rbdry(nbdry),zbdry(nbdry),pds,ier &
             ,n111)
          write (6,*) 'STEPS R,Z,Si,Err = ', rbdry(nbdry),zbdry(nbdry) &
            ,pds(1),ier
        endif

      endif
      if (idebug >= 2) then
        call seva2d(bkx,lkx,bky,lky,c,xout(1),yout(1),pds,ier &
             ,n111)
        write (6,*) 'STEPS simag,psibry,n1,si = ',simag,psibry,n111,pds(1)
        write (6,*) 'STEPS Z,R = ', zgrid(33),(rgrid(i),i=45,45)
        write (6,*) 'STEPS si = ',(psi((i-1)*nw+33),i=45,45)
        call seva2d(bkx,lkx,bky,lky,c,rgrid(45),zgrid(33),pds,ier &
             ,n111)
        write (6,*) 'STEPS R,Z,si = ', rgrid(45),zgrid(33),pds(1)
        write (6,*) 'STEPS lkx,lky = ',lkx,lky

      endif

!-----------------------------------------------------------------------
!-- get weighting function                                            --
!-----------------------------------------------------------------------

      call weight(rgrid,zgrid)
!-----------------------------------------------------------------------
!--  get response functions for MSE                                   --
!-----------------------------------------------------------------------
      if (kstark.gt.0.or.kdomse.gt.0) then
      do 50299 k=1,nstark
        if (rrgam(jtime,k).le.0.0) go to 50299
        call seva2d(bkx,lkx,bky,lky,c,rrgam(jtime,k) &
                    ,zzgam(jtime,k),pds,ier,n111)
        sistark(k)=pds(1)
        sisinow=(simag-pds(1))/sidif
        sigam(k)=sisinow
        fpnow=ffcurr(sisinow,kffcur)
        btgam(k)=fpnow*tmu/rrgam(jtime,k)
50299 continue
      endif

!-----------------------------------------------------------------------
!--  get response functions for MSE-LS                                --
!-----------------------------------------------------------------------
      if (mmbmsels.gt.0.or.kdomsels.gt.0) then
      do 60299 k=1,nmsels
        if (idebug >= 2) then
          write (6,*) 'STEPS MSE-LS k,rrmselt= ', k,rrmselt(jtime,k)
        endif
        if (rrmselt(jtime,k).le.0.0) go to 60299
        call seva2d(bkx,lkx,bky,lky,c,rrmselt(jtime,k) &
                    ,zzmselt(jtime,k),pds,ier,n333)
        simls(k)=pds(1)
        sisinow=(simag-pds(1))/sidif
        sinmls(k)=sisinow
        fpnow=ffcurr(sisinow,kffcur)
        btmls(k)=fpnow*tmu/rrmselt(jtime,k)
        brmls(k)=-pds(3)/rrmselt(jtime,k)
        bzmls(k)=pds(2)/rrmselt(jtime,k)
        if (idebug >= 2) then
          write (6,*) 'STEPS MSE-LS k,rrmselt,br,bz,bt= ', k,rrmselt(jtime,k) &
                       ,brmls(k),bzmls(k),btmls(k)
        endif
60299 continue

      endif
!
      do 51200 k=1,nfound
        xouts(k)=1./xout(k)**2
51200 continue
      nzz=0
      call fluxav(xouts,xout,yout,nfound,psi,rgrid,nw,zgrid,nh, &
                  r2bdry,nzz,sdlobp,sdlbp)
      do 51210 k=1,nfound
        xouts(k)=1./xout(k)
51210 continue
      call fluxav(xouts,xout,yout,nfound,psi,rgrid,nw,zgrid,nh, &
                  r1bdry,nzz,sdlobp,sdlbp)
!-----------------------------------------------------------------------
!--  get metric elements at PSIWANT for edge constraint if needed     --
!----------------------------------------------------------------------
      r1sdry(1)=r1bdry
      r2sdry(1)=r2bdry
      nnn=1
      if (abs(sizeroj(1)-1.0).le.1.e-05_dp.and.kzeroj.eq.1) go to 51977
      if (kzeroj.gt.0) then
       do i=1,kzeroj
       if (sizeroj(i).ge.1.0) sizeroj(i)=0.99999_dp
       siwant=simag+sizeroj(i)*(psibry-simag)
       call surfac(siwant,psi,nw,nh,rgrid,zgrid,rsplt,zsplt,nfounc, &
                    npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nnn, &
                    rmaxis,zmaxis,negcur,kerror)
       if (kerror.gt.0) return
       do 51900 k=1,nfounc
        csplt(k)=1./rsplt(k)**2
51900  continue
       call fluxav(csplt,rsplt,zsplt,nfounc,psi,rgrid,nw,zgrid,nh, &
                  r2sdry(i),nzz,sdlobp,sdlbp)
       do 51920 k=1,nfounc
        csplt(k)=1./rsplt(k)
51920  continue
       call fluxav(csplt,rsplt,zsplt,nfounc,psi,rgrid,nw,zgrid,nh, &
                  r1sdry(i),nzz,sdlobp,sdlbp)
       r2surs = r2sdry(i)*sdlobp
       fpnow = ffcurr(psiwant,kffcur)
       fpnow = fpnow*tmu
       enddo
      endif
51977 continue
!-----------------------------------------------------------------------
!--  get metric elements at PSIWANT for q constraint if needed        --
!-----------------------------------------------------------------------
      if (nqwant.gt.0) then
      do 53999 i=1,nqwant
       siwant=simag+siwantq(i)*(psibry-simag)
       call surfac(siwant,psi,nw,nh,rgrid,zgrid,rsplt,zsplt,nfounc, &
                    npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nnn, &
                    rmaxis,zmaxis,negcur,kerror)
       if (kerror.gt.0) return
       do 53810 k=1,nfounc
        csplt(k)=1./rsplt(k)**2
53810  continue
       call fluxav(csplt,rsplt,zsplt,nfounc,psi,rgrid,nw,zgrid,nh, &
                  r2qdry,nzz,sdlobp,sdlbp)
       do 53820 k=1,nfounc
        csplt(k)=1./rsplt(k)
53820  continue
       call fluxav(csplt,rsplt,zsplt,nfounc,psi,rgrid,nw,zgrid,nh, &
                  r1qdry,nzz,sdlobp,sdlbp)
       r2surq = r2qdry*sdlobp
       fpnow = ffcurr(siwantq(i),kffcur)
       fpnow = fpnow*tmu
       qsiw(i)= abs(fpnow)/twopi*r2surq
       pasmsw(i)=sdlbp/tmu/twopi
53999 continue
      endif
!
      cvolp(ixt)=0.0
      carea=0.0
      xym=xout(1)*yout(1)
      xyma=yout(1)
      do 1450 i=2,nfound
        xyp=xout(i)*yout(i)
        xypa=yout(i)
        dx=xout(i)-xout(i-1)
        cvolp(ixt)=cvolp(ixt)+twopi*(xyp+xym)/2.0*dx
        carea=carea+(xyma+xypa)/2.0*dx
        xym=xyp
        xyma=xypa
 1450 continue
      carea=abs(carea)
      if (iconvr.eq.3) go to 1470
      if ((ix.gt.1).or.(ixout.le.1)) call chisqr(jtime)
 1470 continue
      cvolp(ixt)=abs(cvolp(ixt))*1.0e+06_dp
      vout(jtime)=cvolp(ixt)
      rout(jtime)=(xmax+xmin)/2.*100.
      aout(jtime)=100.*(xmax-xmin)/2.0
      csimag(ixt)=simag
      csibry(ixt)=psibry
      crmaxi(ixt)=rmaxis*100.
      czmaxi(ixt)=zmaxis*100.
      cemaxi(ixt)=emaxis
      cchisq(ixt)=tsaisq(jtime)
      csumip(ixt)=sumip
      tratio(ixt)=cratio
!---------------------------------------------------------------------
!--  get beta and li for constraints                                --
!---------------------------------------------------------------------
      if (fli.gt.0.0.or.fbetan.gt.0.0) then
           call betsli(jtime,rgrid,zgrid,kerror)
           if (kerror.gt.0) return
      endif
!----------------------------------------------------------------------
!--  vertical stabilization information                              --
!----------------------------------------------------------------------
      if (itell.gt.0.and.isetfb.ge.0) then
        delzmm = zmaxis - zmaxis_last
        zmaxis_last=zmaxis
      endif
      if (isetfb.ne.0) then
        fb_plasma(jtime)=abs(brfb(1)*nfbcoil/cpasma(jtime))
        delzmm = zmaxis - zmaxis_last
        zmaxis_last=zmaxis
      elseif (eelip.gt.2.25_dp .and. itell.eq.0) then
        delzmm = zmaxis - zmaxis_last
        zmaxis_last=zmaxis
      endif
!----------------------------------------------------------------------
!--  magnetic axis parameters if needed                              --
!----------------------------------------------------------------------
      if (icinit.gt.0) then
        if ((iconvr.ne.3).and.(ixout.le.1)) go to 1580
      endif
        select case (icurrt)
        case (1)
          go to 1570
        case (2)
          go to 1590
        case (3)
          go to 1595
        case (4)
          go to 1580
        case (5)
          go to 1590
        end select
!
 1570 continue
      rdiml=rmaxis/srma
      cjmaxi=cratio*(sbeta*rdiml+2.*salpha/rdiml)/darea
      if (kvtor.gt.0) then
          cjmaxi=cjmaxi+cratio/darea*sbetaw*rdiml*(rdiml**2-1.)
      endif
      go to 1600
 1580 continue
      nqend=1
      n22=2
      if (errorm.lt.0.1_dp.and.icurrt.eq.4) nqend=nqiter
      do 1585 i=1,nqend
      if (i.gt.1) then
        call currnt(n22,jtime,n22,n22,kerror)
        if (kerror.gt.0) return
      end if

      fcentr=fbrdy**2+sidif*dfsqe
      if (fcentr.lt.0.0) fcentr=fbrdy**2
      fcentr=sqrt(fcentr)*fbrdy/abs(fbrdy)
      rdiml=rmaxis/rzero
      cjmaxi=cratio/darea*(rdiml+rbetap/rdiml)
      if (kvtor.eq.1) then
        rgmvt=(rmaxis/rvtor)**2-1.
        cjmaxi=cjmaxi+cratio/darea*rdiml*rbetaw*rgmvt
      elseif (kvtor.eq.11) then
         ypsm=0.0
         n1set=1
         pres0=prcur4(n1set,ypsm,kppcur)
         prew0=pwcur4(n1set,ypsm,kwwcur)
         rgmvt=(rmaxis/rvtor)**2-1.
         pwop0=prew0/pres0
         ptop0=exp(pwop0*rgmvt)
         pp0= 1.-pwop0*rgmvt
         ppw=rbetaw*rgmvt
         cjmaxi=cjmaxi+(pp0+ppw)*rdiml*ptop0
      endif
      cqmaxi(ixt)=(emaxis**2+1.)*abs(fcentr)/twopi/emaxis &
                  /rmaxis**2/abs(cjmaxi)
      qmaxis=cqmaxi(ixt)
      if (icurrt.eq.4) then
      if (qenp.gt.0.0) enp=enp*qenp/qmaxis
      if (qemp.gt.0.0) emp=emp*qmaxis/qemp
      enf=enp
      emf=emp
      endif
 1585 continue

      return
 1590 continue
      if ((ixt.le.1).and.(icinit.gt.0)) go to 1580
      fcentr=ffcurr(x000,kffcur)
      if (kvtor.eq.0) then
        cjmaxi=(rmaxis*ppcurr(x000,kppcur) &
             +fpcurr(x000,kffcur)/rmaxis)*cratio/darea
      else
        cjmaxi=0.0
        do j=1,kppcur
          cjmaxi=rjjjx(j)*brsp(nfcoil+j)+cjmaxi
        enddo
        do j=1,kffcur
          cjmaxi=rjjfx(j)*brsp(nfcoil+kppcur+j)+cjmaxi
        enddo
        do j=1,kwwcur
          cjmaxi=rjjwx(j)*brsp(nfnpcr+j)+cjmaxi
        enddo
        cjmaxi=cjmaxi/darea
      endif
      go to 1600
 1595 continue
 1600 continue
      cqmaxi(ixt)=(emaxis**2+1.)*abs(fcentr)/twopi/emaxis &
                   /rmaxis**2/abs(cjmaxi)
      qmaxis=cqmaxi(ixt)
      return
 2100 format(/,1x,'shot',i6,' at ',i6,' ms ','** Problem in BOUND **')
11001 format(/,1x,'** 2nd seperatrix **',2x,e10.3,2x,i4)
      end subroutine steps
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          inicur initializes the current density distribution.    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**           ks..............time slice number                      **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine inicur(ks)
      include 'eparmdud129.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      save isicinit,zelips
!
      if (ivacum.gt.0) return
!----------------------------------------------------------------------
!--  icinit=1 uniform and elliptical flux surfaces                   --
!--         2 parabolic and elliptical                               --
!--        -2 use ESAVE.DAT at first slice and previous slice later  --
!--       -12 parabolic and elliptical first slice and previous      --
!--           slice subsequently                                     --
!----------------------------------------------------------------------
      if (ks.eq.1) isicinit=icinit
      if (isicinit.lt.0)   then
        if (ks.gt.1) then
          icinit=isicinit
          return
        else
          if (isicinit.eq.-12) then
            icinit=2
            go to 1200
          endif
          if (icinit.eq.-2)  go to 1100
        endif
      endif
        select case (abs(icinit))
        case (1)
          go to 100
        case (2)
          go to 1100
        end select
      return
  100 continue
      if (aelip.gt.0.0) go to 150
      aelip=0.50_dp
      do 140 i=1,limitr
        erho=sqrt((xlim(i)-relip)**2+((ylim(i)-zelip)/eelip)**2)
        aelip=min(aelip,erho)
  140 continue
  150 continue
      delcur=1.
      sumi=0.0
      do 300 i=1,nw
      do 300 j=1,nh
        kk=(i-1)*nh+j
        erho=sqrt((rgrid(i)-relip)**2+((zgrid(j)-zelip)/eelip)**2)
        if (erho.gt.aelip) go to 300
        pcurrt(kk)=delcur*zero(kk)/rgrid(i)
        sumi=sumi+pcurrt(kk)
  300 continue
      cratio=pasmat(ks)/sumi
      do 320 i=1,nwnh
        pcurrt(i)=pcurrt(i)*cratio*zero(i)
  320 continue
      return
!
 1100 continue
      if (icinit.gt.0) go to 1200
      open(unit=nsave,form='unformatted',file='esave.dat', &
           status='old',err=1200)
      read (nsave) mw,mh
      read (nsave) xpsi
      read (nsave) brsp
      read (nsave) www
      read (nsave) emaxis, rmaxis, fcentr
      close(unit=nsave)
      return
 1200 continue
      if (ks.eq.1)  zelips=zelip
      if (zelip.gt.1.e5_dp .or. zelips.gt.1.e5_dp) then
        zelip=1.447310_dp*(expmpi(ks,37)-expmpi(ks,43)) &
             +0.692055_dp*(expmpi(ks,57)-expmpi(ks,53)) &
             +0.728045_dp*(silopt(ks,27)-silopt(ks,37)) &
             +2.047150_dp*(silopt(ks,2) -silopt(ks,11))
        zelip=zelip*1.e6_dp/pasmat(ks)
        zbound=zelip
        eelip=1.5_dp
      endif
!----------------------------------------------------------------
!-- set zelip=0.0 if bad signals              96/06/24         --
!----------------------------------------------------------------
      if (abs(fwtmp2(37)).le.1.0e-30_dp)  zelip=0.0
      if (abs(fwtmp2(43)).le.1.0e-30_dp)  zelip=0.0
      if (abs(fwtmp2(57)).le.1.0e-30_dp)  zelip=0.0
      if (abs(fwtmp2(53)).le.1.0e-30_dp)  zelip=0.0
      if (abs(fwtsi(27)).le.1.0e-30_dp)  zelip=0.0
      if (abs(fwtsi(37)).le.1.0e-30_dp)  zelip=0.0
      if (abs(fwtsi(2)).le.1.0e-30_dp)  zelip=0.0
      if (abs(fwtsi(11)).le.1.0e-30_dp)  zelip=0.0
!
      do 1300 i=1,nw
      do 1300 j=1,nh
        kk=(i-1)*nh+j
        erho=sqrt((rgrid(i)-relip)**2+((zgrid(j)-zelip)/eelip)**2)
        xpsi(kk)=(erho/aelip)**2
 1300 continue
      return
      end subroutine inicur

!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          vescur computes the currents induced in the vessel      **
!**          segments due to E coils and F coils.                    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/11/85..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine vescur(jtime)
      include 'eparmdud129.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
!
      if (ivesel.eq.5) return
!
 2000 continue
      if (ivesel.eq.1) return
      do 2100 i=1,nvesel
 2100 vcurrt(i)=vloopt(jtime)/rsisvs(i)
      return
      end subroutine vescur


