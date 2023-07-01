!**********************************************************************
!>
!!    pflux computes the poloidal fluxes on the r-z grid.
!!
!!
!!    @param niter : current profile (outer) loop iteration index
!!
!!    @param nnin : equilibrium (inner) loop iteration index
!!
!!    @param ntotal : total iteration index (current+equilibirum loops)
!!
!!    @param jtime : time index
!!
!!    @param kerror : error flag
!!
!**********************************************************************
      subroutine pflux(niter,nnin,ntotal,jtime,kerror)
      use var_buneman, only: rgrid1,delrgrid,delz,drdz2
      use commonblocks, only: c,bkx,bky,psiold,psipold
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      integer*4 ef_init_cycred_data
      integer*4, intent(in) :: niter,nnin,ntotal,jtime
      integer*4, intent(out) :: kerror
      integer*4 i,j,ii,jj,kk,kkkk,kknh,k1,kk1,k2,mj,mk,mjnh,mknh,mj1,mk1, &
                iinow,imd,idelr,idelz,jtop,jlow,ier,icurfb
      real*8 tempsum1,tempsum2,val,cdelznow,deltaz,difnow,psibar, &
             psic1,psic2,psilu_psill,psiul_psiuu,rcurnow,sumc,vcurfi, &
             zdwn_now,zup_now
      real*8 pds(6)
      real*8,dimension(nwnh) :: psikkk,gfbsum
      logical vfeed

      kerror = 0
      vfeed=(isetfb.ne.0).and.(niter.gt.2.or.nnin.gt.2)
!----------------------------------------------------------------------------
!--   save flux from current iterations before update                      --
!----------------------------------------------------------------------------
      psiold=psi
      psipold=psipla
!
      buneman_green: if ((ibunmn.eq.1).or. &
                        ((ibunmn.eq.2).and.(errorm.gt.errcut))) then
!-----------------------------------------------------------------------------
!--   Buneman's method of obtaining psi at the inner grid points            --
!--   only plasma contribution                                              --
!-----------------------------------------------------------------------------
!---  These must be brought into the integrals 
      rgrid1=rgrid(1)
      delrgrid=rgrid(2)-rgrid(1)
      delz=zgrid(2)-zgrid(1)
      drdz2=(delrgrid/delz)**2
!---------------------------------------------------------------------
!--   optional vertical feedback control                          --
!---------------------------------------------------------------------
      if (isetfb.ne.0) then
        if(nw.gt.30) ioffr=ioffr*((nw+1)/33)
        if(nh.gt.30) ioffz=ioffz*((nh+1)/33)
        imd=(nw*nh)/2+1+nh*ioffr+ishiftz
        kct1=imd-ioffz
        kct2=imd+ioffz
        deltaz=zgrid(ioffz+nh/2+1)-zgrid(-ioffz+nh/2+1)
      endif
!----------------------------------------------------------------------------
!--   obtain fluxes by inverting del*                                      --
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
          if(ntotal.eq.1) psikkk=0.0
        else
          icurfb=vcurfb(2)
          if(mod(ntotal-iinow,icurfb).eq.0) psikkk=psiold
        endif
      endif
!-----------------------------------------------------------------------
!--   boundary terms                                                  --
!-----------------------------------------------------------------------
      do j=1,nh
        kk=(nw-1)*nh+j
        tempsum1=0.
        tempsum2=0.
        do ii=1,nw
         do jj=1,nh
          kkkk=(ii-1)*nh+jj
          mj=abs(j-jj)+1
          mk=(nw-1)*nh+mj
          tempsum1=tempsum1-gridpc(mj,ii)*pcurrt(kkkk)
          tempsum2=tempsum2-gridpc(mk,ii)*pcurrt(kkkk)
         enddo
        enddo
        psi(j) =tempsum1
        psi(kk)=tempsum2
      enddo
      do i=2,nw-1
        kk1=(i-1)*nh
        kknh=kk1+nh
        kk1=kk1+1
        tempsum1=0.
        tempsum2=0.
        do ii=1,nw
         do jj=1,nh
          kkkk=(ii-1)*nh+jj
          mj1=abs(jj-1)+1
          mjnh=abs(nh-jj)+1
          mk1=(i-1)*nh+mj1
          mknh=(i-1)*nh+mjnh
          tempsum1=tempsum1-gridpc(mk1 ,ii)*pcurrt(kkkk)
          tempsum2=tempsum2-gridpc(mknh,ii)*pcurrt(kkkk)
         enddo
        enddo
        psi(kk1 )=tempsum1
        psi(kknh)=tempsum2
      enddo
!-------------------------------------------------------------------------
!--   get flux at inner points by inverting del*, only plasma flux
!--   first set up plasma currents, single cyclic method gets factor of 2
!-------------------------------------------------------------------------
      if (isolve.eq.0) then
!       original buneman solver method
        do i=2,nw-1
          do j=2,nh-1
            kk=(i-1)*nh+j
            psi(kk)=tmu2*pcurrt(kk)*rgrid(i)
          enddo   
        enddo   
        call buneto(psi,nw,nh)
      else 
!       new faster single cyclic reduction method
        do i=2,nw-1
          do j=2,nh-1
            kk=(i-1)*nh+j
            psi(kk)=tmu2*pcurrt(kk)*rgrid(i)*2.0
          enddo
        enddo
        if (ntotal.le.1) then
!         call the initialization and check the result
          kerror=ef_init_cycred_data()
          if(kerror.eq.1) return
        endif
        call pflux_cycred(psi,kerror)
        if(kerror.gt.0) return
      endif
      psi=-psi
!----------------------------------------------------------------------------
!--   optional symmetrized solution                                        --
!----------------------------------------------------------------------------
      if (symmetrize) then
        difpsi=0.0
        do i=1,nw
          do j=1,nh/2
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
          enddo
        enddo
      endif
!------------------------------------------------------------------------
!--   optional damping out the m=1 vertical eigen mode                 --
!------------------------------------------------------------------------
      if ((abs(vcurfb(1)).ge.1.e-6_dp) .and. (ntotal.ge.5)) then
!-----------------------------------------------------------------------
!--     sum Green's functions for m=1 eigen mode                      --
!-----------------------------------------------------------------------
        idelr=nw/8
        idelz=nh/12
        jtop=(nh+1)/2+nh/4-1
        do i=1,nw
         do j=1,nh
          kk=(i-1)*nh+j
          gfbsum(kk)=0.0
          do ii=1+idelr,nw-idelr
           do jj=jtop-idelz,jtop+idelz
            mj=abs(j-jj)+1
            mk=(nw-1)*nh+mj
            gfbsum(kk)=gfbsum(kk)+gridpc(mj,ii)
           enddo
          enddo
         enddo
        enddo
        jlow=(nh+1)/2-nh/4+1
        do i=1,nw
         do j=1,nh
          kk=(i-1)*nh+j
          do ii=1+idelr,nw-idelr
           do jj=jlow-idelz,jlow+idelz
            mj=abs(j-jj)+1
            mk=(nw-1)*nh+mj
            gfbsum(kk)=gfbsum(kk)-gridpc(mj,ii)
           enddo
          enddo
         enddo
        enddo
!---------------------------------------------------------------------
!--     get damping currents                                        --
!---------------------------------------------------------------------
        tvfbrt(ntotal)=0.0
        do i=1+idelr,nw-idelr
          do j=jtop-idelz,jtop+idelz
            kk=(i-1)*nh+j
            tvfbrt(ntotal)=tvfbrt(ntotal)+(psikkk(kk)-psiold(kk))
          enddo
        enddo
        do i=1+idelr,nw-idelr
          do j=jlow-idelz,jlow+idelz
            kk=(i-1)*nh+j
            tvfbrt(ntotal)=tvfbrt(ntotal)-(psikkk(kk)-psiold(kk))
          enddo
        enddo
        vcurfi=vcurfb(1)*ipmhd(jtime)/abs(tvfbrt(ntotal))
        tvfbrt(ntotal)=tvfbrt(ntotal)*vcurfi
        psi(1:nwnh)=psi(1:nwnh)+gfbsum(1:nwnh)*tvfbrt(ntotal)
      endif
!--------------------------------------------------------------------
!--   optional vertical feedback control                           --
!--    psi(kct1)  !psi at lower point                              --
!--    psi(kct2) ! at upper point                                  --
!--    psilu     !psi (per amp) at lower point due to upper coils  --
!--    brfb(1) is lower current                                    --
!--------------------------------------------------------------------
      if (vfeed) then
        psiul_psiuu=-grdfdb(kct2,1)
        psilu_psill=grdfdb(kct1,1)
        zdwn_now=0
        zup_now=0
        zcurnow=0
        rcurnow=0
        zcontr=sum(yout(1:nfound))/nfound
        sumc=sum(pcurrt(1:nwnh))
        do i=1,nw
         do j=1,nh
          kk=(i-1)*nh+j
          zcurnow=zcurnow+pcurrt(kk)*zgrid(j)
          rcurnow=rcurnow+pcurrt(kk)*rgrid(j)
         enddo
        enddo
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
      endif
!----------------------------------------------------------------------------
!--   add flux from external coils                                         --
!----------------------------------------------------------------------------
      do kk=1,nwnh
        psipla(kk)=psi(kk)
        if(ivesel.gt.0) &
          psi(kk)=psi(kk)+sum(gridvs(kk,:)*vcurrt)
        if (iecurr.eq.1) then
          psi(kk)=psi(kk)+sum(gridec(kk,:)*ecurrt)
        elseif (iecurr.eq.2) then
          psi(kk)=psi(kk)+sum(gridec(kk,:)*cecurr)
        endif
        if (vfeed) psi(kk)=psi(kk)+grdfdb(kk,1)*brfb(2)
        psi(kk)=psi(kk)+sum(gridfc(kk,:)*brsp(1:nfsum))
        if (iacoil.gt.0) &
          psi(kk)=psi(kk)+sum(gridac(kk,:)*caccurt(jtime,:))
      enddo

      else buneman_green
!-----------------------------------------------------------------------------
!--   ibunmn=0 and 2 when errorm less than errcut                        --
!--   Green's integral method of obtaining flux, can be computationally     --
!--   intensive                                                             --
!-----------------------------------------------------------------------------
      psi=0.0
      do i=1,nw
       do j=1,nh
        kk=(i-1)*nh+j
        psi(kk)=psi(kk)+sum(gridfc(kk,:)*brsp(1:nfsum))
        if(ivesel.gt.0) &
          psi(kk)=psi(kk)+sum(gridvs(kk,:)*vcurrt)
        if (iecurr.eq.1) then
          psi(kk)=psi(kk)+sum(gridec(kk,:)*ecurrt)
        elseif (iecurr.eq.2) then
          psi(kk)=psi(kk)+sum(gridec(kk,:)*cecurr)
        endif
        if (iacoil.gt.0) &
          psi(kk)=psi(kk)+sum(gridac(kk,:)*caccurt(jtime,:))
        psipla(kk)=psi(kk)
        if (ivacum.eq.0) then
          do ii=1,nw
            do jj=1,nh
              kkkk=(ii-1)*nh+jj
              mj=abs(j-jj)+1
              mk=(i-1)*nh+mj
              psi(kk)=psi(kk)+gridpc(mk,ii)*pcurrt(kkkk)
            enddo
          enddo
          psipla(kk)=psi(kk)-psipla(kk)
        endif 
       enddo
      enddo
!
      ibunmn=0
!------------------------------------------------------------------------------
!--   relaxation before return                                               --
!------------------------------------------------------------------------------
      endif buneman_green
!----------------------------------------------------------------------------
!--   rigid vertical shift correction ?                                    --
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
!--   relax flux if needed                                                 --
!----------------------------------------------------------------------------
      if ((ntotal.gt.1) .and. (abs(relax-1.0).ge.1.0e-03_dp)) then
        psi(1:nwnh)=relax*psi(1:nwnh)+(1.-relax)*psiold(1:nwnh)
        psipla(1:nwnh)=relax*psipla(1:nwnh)+(1.-relax)*psipold(1:nwnh)
      endif
!
      return
      end subroutine pflux
