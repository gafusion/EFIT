#include "config.f"
!**********************************************************************
!>
!!    pflux computes the poloidal fluxes on the r-z grid.
!!
!!
!!    @param niter inner equilibrium loop iteration index
!!
!!    @param nnin current profile loop iteration index?
!!
!!    @param ntotal current profile loop iteration index (how is this 
!!                  different than nnin?)
!!
!!    @param jtime  time index
!!
!!    @param kerror error flag
!!
!**********************************************************************
      subroutine pflux(niter,nnin,ntotal,jtime,kerror)
      use set_kinds 
      use var_buneman, only: rgrid1,delrgrid,delz,drdz2
      use commonblocks, only: c,wk,copy,bkx,bky,psiold,psipold,psipp
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      real(dp), dimension(:), allocatable ::   work
      integer*4 initresult
      dimension pds(6)
      real*8,dimension(:),allocatable :: psikkk,gfbsum
      data initfb/0/,init/0/

      kerror = 0
      ALLOCATE(psikkk(nwnh),gfbsum(nwnh))

      vfeed=(isetfb.ne.0).and.(init.ne.0).and.(niter.gt.2.or.nnin.gt.2)
      if(ivesel.gt.10) return
!----------------------------------------------------------------------------
!--   save flux from current iterations before update                      --
!----------------------------------------------------------------------------
      !$omp target teams distribute parallel do collapse(2)
      do i=1,nw
        do j=1,nh
          kk=(i-1)*nh+j
          psiold(kk)=psi(kk)
          psipold(kk)=psipla(kk)
        enddo
      enddo
!
      buneman_green: if ((ibunmn.eq.1).or.(ibunmn.eq.3).or. &
         ((ibunmn.eq.2).and.(errorm.gt.errcut)).or. &
         ((ibunmn.eq.4).and.(errorm.gt.errcut))) then
!-----------------------------------------------------------------------------
!--   Buneman's method of obtaining psi at the inner grid points            --
!--   only plasma contribution                                              --
!-----------------------------------------------------------------------------
      if(init.le.0) then
!-----  These must be brought into the integrals 
        rgrid1=rgrid(1)
        delrgrid=rgrid(2)-rgrid(1)
        delz=zgrid(2)-zgrid(1)
        drdz2=(delrgrid/delz)**2
!---------------------------------------------------------------------
!--     optional vertical feedback control                          --
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
          if(ntotal.eq.1) then
            !$omp target teams distribute parallel do collapse(2)
            do i=1,nw
              do j=1,nh
                kk=(i-1)*nh+j
                psikkk(kk)=0.0
              enddo
            enddo 
          endif
        else
          icurfb=vcurfb(2)
          if(mod(ntotal-iinow,icurfb).eq.0) then
            !$omp target teams distribute parallel do collapse(2)
            do i=1,nw
              do j=1,nh
                kk=(i-1)*nh+j
                psikkk(kk)=psiold(kk)
              enddo
            enddo
          endif
        endif
      endif
!-----------------------------------------------------------------------
!--   boundary terms                                                  --
!-----------------------------------------------------------------------
      !$omp target teams distribute reduction(+:tempsum1,tempsum2)
      do j=1,nh
        kk=(nw-1)*nh+j
        tempsum1=0.
        tempsum2=0.
#ifdef USE_OPENMP_NV
        !$omp parallel do reduction(+:tempsum1,tempsum2) collapse(2)
#elif defined (USE_OPENMP_AMD)
        !$omp loop bind(teams) reduction(+:tempsum1,tempsum2) collapse(2)
#endif
        do ii=1,nw
         do jj=1,nh
          kkkk=(ii-1)*nh+jj
          mj=abs(j-jj)+1
          mk=(nw-1)*nh+mj
          tempsum1=tempsum1-gridpc(mj,ii)*pcurrt(kkkk)
          tempsum2=tempsum2-gridpc(mk,ii)*pcurrt(kkkk)
         enddo
        enddo
        psi(j)=tempsum1
        psi(kk)=tempsum2
      enddo

      !$omp target teams distribute reduction(+:tempsum1,tempsum2)
      do i=2,nw-1
        kk1=(i-1)*nh
        kknh=kk1+nh
        kk1=kk1+1
        tempsum1=0.
        tempsum2=0.
#ifdef USE_OPENMP_NV
        !$omp parallel do reduction(+:tempsum1,tempsum2) collapse(2)
#elif defined (USE_OPENMP_AMD)
        !$omp loop bind(thread) reduction(+:tempsum1,tempsum2) collapse(2)
#endif
        do ii=1,nw
         do jj=1,nh
          kkkk=(ii-1)*nh+jj
          mj1=abs(jj-1)+1
          mjnh=abs(nh-jj)+1
          mk1=(i-1)*nh+mj1
          mknh=(i-1)*nh+mjnh
          tempsum1=tempsum1-gridpc(mk1,ii)*pcurrt(kkkk)
          tempsum2=tempsum2-gridpc(mknh,ii)*pcurrt(kkkk)
         enddo
        enddo
        psi(kk1)=tempsum1
        psi(kknh)=tempsum2
      enddo
!-------------------------------------------------------------------------
!--   get flux at inner points by inverting del*, only plasma flux
!--   first set up plasma currents, single cyclic method gets factor of 2
!-------------------------------------------------------------------------
      allocate(work(SIZE(psi)))
      if (isolve.eq.0) then
!       original buneman solver method
        !$omp target teams distribute parallel do collapse(2)
        do i=2,nw-1
          do j=2,nh-1
            kk=(i-1)*nh+j
            psi(kk)=tmu2*pcurrt(kk)*rgrid(i)
          enddo   
        enddo   
        call buneto(psi,INT(nw,i4),INT(nh,i4),work)
      else 
!       new faster single cyclic reduction method
        !$omp target teams distribute parallel do collapse(2)
        do i=2,nw-1
          do j=2,nh-1
            kk=(i-1)*nh+j
            psi(kk)=tmu2*pcurrt(kk)*rgrid(i)*2.0
          enddo
        enddo
        call pflux_cycred(psi,work,kerror)
        if (kerror.gt.0) return
      endif
      deallocate(work)
      !$omp target teams distribute parallel do collapse(2)
      do i=1,nw
        do j=1,nh
          kk=(i-1)*nh+j
          psi(kk)=-psi(kk)
        enddo
      enddo
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
        if (initfb.eq.0) then
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
          initfb=-1
        endif
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
        if (initfb.eq.-1) then
          vcurfi=vcurfb(1)*cpasma(jtime)/abs(tvfbrt(ntotal))
          initfb=1
        endif
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
        ! TODO: unreachable code?
!        if (.not.vfeed) then
!          brfb(1)=0.
!          brfb(2)=0.
!        endif
      endif
!----------------------------------------------------------------------------
!--   add flux from external coils                                         --
!----------------------------------------------------------------------------
      !$omp target teams distribute parallel do
      do kk=1,nwnh
        psipla(kk)=psi(kk)
        if(ivesel.gt.0) &
          psi(kk)=psi(kk)+sum(gridvs(kk,1:nvesel)*vcurrt(1:nvesel))
        if (iecurr.eq.1) &
          psi(kk)=psi(kk)+sum(gridec(kk,1:nesum)*ecurrt(1:nesum))
        if (iecurr.eq.2) &
          psi(kk)=psi(kk)+sum(gridec(kk,1:nesum)*cecurr(1:nesum))
        if (vfeed) psi(kk)=psi(kk)+grdfdb(kk,1)*brfb(2)
        psi(kk)=psi(kk)+sum(gridfc(kk,1:nfcoil)*brsp(1:nfcoil))
        if (iacoil.gt.0) &
          psi(kk)=psi(kk)+sum(gridac(kk,1:nacoil)*caccurt(jtime,1:nacoil))
      enddo

      else buneman_green
!-----------------------------------------------------------------------------
!--   ibunmn=0, and 2,4 when errorm less than errcut                        --
!--   Green's integral method of obtaining flux, can be computationally     --
!--   intensive                                                             --
!-----------------------------------------------------------------------------
      psi(:) = 0.0
      !$omp target teams distribute parallel do simd collapse(2) 
      do i=1,nw
       do j=1,nh
        tempsum1 = 0.
        kk=(i-1)*nh+j
        psi(kk)=psi(kk)+sum(gridfc(kk,1:nfcoil)*brsp(1:nfcoil))
        if(ivesel.gt.0) &
          psi(kk)=psi(kk)+sum(gridvs(kk,1:nvesel)*vcurrt(1:nvesel))
        if (iecurr.eq.1) &
          psi(kk)=psi(kk)+sum(gridec(kk,1:nesum)*ecurrt(1:nesum))
        if (iecurr.eq.2) &
          psi(kk)=psi(kk)+sum(gridec(kk,1:nesum)*cecurr(1:nesum))
        if (iacoil.gt.0) &
          psi(kk)=psi(kk)+sum(gridac(kk,1:nacoil)*caccurt(jtime,1:nacoil))
        psipla(kk)=psi(kk)
        if (ivacum.le.0) then
          !$omp simd reduction(+:tempsum1)
          do ii=1,nw
            do jj=1,nh
              kkkk=(ii-1)*nh+jj
              mj=abs(j-jj)+1
              mk=(i-1)*nh+mj
              tempsum1=tempsum1+gridpc(mk,ii)*pcurrt(kkkk)
            enddo
          enddo
          psi(kk)=psi(kk)+tempsum1
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
      DEALLOCATE(psikkk,gfbsum)
!
      return
      end subroutine pflux

!**********************************************************************
!>
!!    residu computes the flux variations on the r-z grid.
!!
!!
!!    @param nx :
!!
!!    @param jtime : time index
!!
!**********************************************************************
      subroutine residu(nx,jtime)
      use commonblocks,only: psiold,psipold,psipp
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      if(ivacum.gt.0) return
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
      if(errorm.le.error) idone=1
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
          if(nx.eq.1) write (nttyo,'(x)')
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
      call flush(6)
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
      call flush(6)
#ifdef DEBUG_LEVEL1
      write (nttyo,*) 'cratio,cratio_ext,cratiop_ext,cratiof_ext= ', &
        cratio,cratio_ext,cratiop_ext,cratiof_ext
      write (nttyo,*) 'scalepp_ext,scaleffp_ext= ', &
        scalepp_ext,scaleffp_ext
#endif
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
!>
!!    step computes the dimensionless poloidal fluxes for
!!    the r-z grid
!!
!!
!!    @param ix :
!!
!!    @param ixt :
!!
!!    @param ixout :
!!
!!    @param jtime : time index
!!
!!    @param kerror : error flag
!!
!**********************************************************************
      subroutine steps(ix,ixt,ixout,jtime,kerror)
      use commonblocks,only: c,wk,copy,bkx,bky,zeros,xouts,youts, &
           rsplt,zsplt,csplt
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      dimension pds(6)
      integer*4 iii
      real*8 :: zmaxis_last = 0.0
      data isplit/8/,psitol/1.0e-04_dp/,cdum/1.0/
      save xguess, yguess, xltrac, radbou
!
      ! initialize variables
      zmaxis=0.0
      simag=0.0
      vout=0.0
      rout=0.0
      aout=0.0
!
      if(ivacum.gt.0) return
      if (ixt.le.1) then
        xguess=(rgrid(1)+rgrid(nw))/2.
        yguess=(zgrid(1)+zgrid(nh))/2.
        if(zbound.ne.0.0) yguess=zbound
        if(rbound.ne.0.0) xguess=rbound
        xltrac=xlmin
        if(ibound.eq.-1) xltrac=xlmax ! hardcoded option only...
        radbou=(xguess+xltrac)/2.
      endif
!----------------------------------------------------------------------
!--   first set up bi-cubic spline interpolation in findax           --
!----------------------------------------------------------------------
      m10=10

#ifdef DEBUG_LEVEL2
      write (6,*) 'Entering findax'
#endif
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
                  psibry,rseps(:,jtime),zseps(:,jtime),m10, &
                  xout,yout,nfound,psi,xmin,xmax,ymin,ymax, &
                  zxmin,zxmax,rymin,rymax,dpsi,bpol,bpolz, &
                  limitr,xlim,ylim,limfag,ixt,jtime,kerror)

      if(kerror.gt.0) return
      if (nsol.gt.0) then

        ier=0 ! set only for consistent output? (not useful...)   
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

#ifdef DEBUG_LEVEL2
        call seva2d(bkx,lkx,bky,lky,c,xout(1),yout(1),pds,ier,n111)
        write (6,*) 'STEPS simag,psibry,n1,si = ',simag,psibry,n111,pds(1)
        write (6,*) 'STEPS Z,R = ', zgrid(33),(rgrid(i),i=45,45)
        write (6,*) 'STEPS si = ',(psi((i-1)*65+33),i=45,45)
        call seva2d(bkx,lkx,bky,lky,c,rgrid(45),zgrid(33),pds,ier,n111)
        write (6,*) 'STEPS R,Z,si = ', rgrid(45),zgrid(33),pds(1)
#endif
      endif
!-----------------------------------------------------------------------
!--   Trace boundary, first check for counter beam injection          --
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
!--   find magnetic axis and poloidal flux at axis simag             --
!----------------------------------------------------------------------
      m20=20
#ifdef DEBUG_LEVEL2
      write (6,*) 'Entering findax after m20 set'
#endif
      call findax(nw,nh,rgrid,zgrid,rmaxis,zmaxis,simag, &
                  psibry,rseps(1,jtime),zseps(1,jtime),m20, &
                  xout,yout,nfound,psi,xmin,xmax,ymin,ymax, &
                  zxmin,zxmax,rymin,rymax,dpsi,bpol,bpolz, &
                  limitr,xlim,ylim,limfag,ixt,jtime,kerror)
      if(kerror.gt.0) return
      sidif=simag-psibry
      eouter=(ymax-ymin)/(xmax-xmin)
      zplasm=(ymin+ymax)/2.
      aouter=(xmax-xmin)/2.
!-----------------------------------------------------------------------
!--   force free current in the scrape-off layer                      --
!-----------------------------------------------------------------------
      SOL_curr: if (icutfp.eq.2) then
        xvsmin=1.e10_dp
        xvsmax=-1.e10_dp
        if (limvs.eq.0) then
          itot=isplit*isplit
          do k=1,nvesel
            call splitc(isplit,rsplt,zsplt,csplt, &
                  rvs(k),zvs(k),wvs(k),hvs(k),avs(k),avs2(k),cdum)
            do kk=2,itot
              call seva2d(bkx,lkx,bky,lky,c,rsplt(kk),zsplt(kk), &
                          pds,ier,n111)
              write(6,*) 'rgrid(i),zgrid(j)',rsplt(kk),zsplt(kk)
              write(6,*) 'lkx, lky',lkx,lky
              write(6,*) 'pds,ier,n111', pds,ier,n111
              xvsmin=min(xvsmin,pds(1))
              xvsmax=max(xvsmax,pds(1))
            enddo
          enddo
        else
          do k=1,limitr-1
            delx=xlim(k+1)-xlim(k)
            dely=ylim(k+1)-ylim(k)
            dels=sqrt(delx**2+dely**2)
            nn=dels/0.002_dp
            nn=max(5,nn)
            delx=delx/(nn-1)
            dely=dely/(nn-1)
            do kk=2,nn
              xww=xlim(k)+delx *(kk-1)
              yww=ylim(k)+dely *(kk-1)
              call seva2d(bkx,lkx,bky,lky,c,xww,yww,pds,ier,n111)
              xvsmin=min(xvsmin,pds(1))
              xvsmax=max(xvsmax,pds(1))
            enddo
          enddo
        endif
!--------------------------------------------------------------
!--     exclude private region flux                          --
!--------------------------------------------------------------
        xvsmax=psibry
!--------------------------------------------------------------
!--     possible second separatrix                           --
!--------------------------------------------------------------
        rsepex=-999.
        yvs2=1000.
        skipvs: if (kskipvs.ne.0) then
        avebp=cpasma(jtime)*tmu/aouter
        bpmin=avebp
        ! TODO: sibpmin has not yet been defined...
!        sibpold=sibpmin
        do i=1,nw
          do j=1,nh
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
            endif
          enddo
        enddo

        bmin_ave: if (bpmin.ne.avebp) then
          relsi=abs((sibpmin-psibry)/sidif)
          if (bpmin.le.0.10_dp*avebp.and.relsi.gt.0.005_dp) then
!------------------------------------------------------------------
!--         find second separatrix                               --
!------------------------------------------------------------------
            do j=1,40
              call seva2d(bkx,lkx,bky,lky,c,xs,ys,pds,ier,n666)
              write(6,*) 'xs,ys',xs, ys
              write(6,*) 'lkx, lky',lkx,lky
              write(6,*) 'pds,ier,n111', pds,ier,n111

              det=pds(5)*pds(6)-pds(4)*pds(4)
              if (abs(det).lt.1.0e-15_dp) exit
              xerr=(-pds(2)*pds(6)+pds(4)*pds(3))/det
              yerr=(-pds(5)*pds(3)+pds(2)*pds(4))/det
              xs=xs+xerr
              ys=ys+yerr
              if (xerr*xerr+yerr*yerr.lt.1.0e-12_dp) exit
            enddo
            if (xerr*xerr+yerr*yerr.ge.1.0e-12_dp) then
              epssep=xerr*xerr+yerr*yerr
              write (nttyo,11001) epssep,ixt
              if (iand(iout,1).ne.0) write (nout,11001) epssep,ixt
            else
              epssep=-1.0
            endif
            if (epssep.lt.1.0e-10_dp) then
              sibpmin=pds(1)
              yvs2=ys
              rsepex=xs
              relsi=abs((sibpmin-psibry)/sidif)
              if (relsi.gt.0.005_dp) then
                ! TODO: what is the intention here? sibpold is unset...
!                if (ixt.gt.1) sibpmin=sibpmin*(1.-vsdamp)+sibpold*vsdamp
                xvsmin=max(xvsmin,sibpmin)
              endif
            endif
          endif
        endif bmin_ave
        endif skipvs
        if (alphafp.ge.0.0) then
          xpsimin=xvsmin+alphafp*(xvsmax-xvsmin)
        else
          ! TODO: xpsipmins has not yet been defined...
!          xpsimino=xpsimins
          xpsimin=abs(alphafp)*xvsmax
          xpsimins=xpsimin
          ! TODO: what is the intention here? xpsimino is unset...
!          if (ixt.gt.1) xpsimin=xpsimin*(1.-vsdamp)+xpsimino*vsdamp
        endif
        xpsialp=xpsimin
        xpsimin=sidif/(simag-xpsimin)
      endif SOL_curr 
!-----------------------------------------------------------------------
!--   get normalized flux function XPSI                               --
!-----------------------------------------------------------------------
      do i=1,nw
       do j=1,nh
        kk=(i-1)*nh+j
        if (icutfp.eq.0) then
          xpsi(kk)=1.1_dp
          if((rgrid(i).lt.xmin).or.(rgrid(i).gt.xmax)) cycle
          if((zgrid(j).lt.ymin).or.(zgrid(j).gt.ymax)) cycle
          xpsi(kk)=(simag-psi(kk))/sidif
        else
          if (zero(kk).gt.0.0005_dp) then
            xpsi(kk)=(simag-psi(kk))/sidif
            if (xpsi(kk)*xpsimin.le.1.0.and.xpsi(kk)*xpsimin.ge.0.0) then
              if((rgrid(i).lt.rminvs).or.(rgrid(i).gt.rmaxvs)) &
                   xpsi(kk)=1000.
              if ((zgrid(j).lt.zminvs).or.(zgrid(j).gt.zmaxvs)) then
                   xpsi(kk)=1000.
              endif
              if(xpsi(kk).lt.1.0.and.zgrid(j).lt.ymin) xpsi(kk)=1000.
              if(xpsi(kk).lt.1.0.and.zgrid(j).gt.ymax) xpsi(kk)=1000.
              if(abs(zgrid(j)).gt.abs(yvs2).and. &
                zgrid(j)*yvs2.gt.0.0) xpsi(kk)=1000.
            endif
          else
            xpsi(kk)=1000.
          endif
        endif
       enddo
      enddo
!-----------------------------------------------------------------------
!--   get SOL flux if needed                                          --
!-----------------------------------------------------------------------
      if (nsol.gt.0) then
        call seva2d(bkx,lkx,bky,lky,c,rsol(1),zsol(1), &
                     pds,ier,n111)
        wsisol=pds(1)
#ifdef DEBUG_LEVEL2
        write (6,*) 'STEPS R,Z,Si,Err = ', rsol(1),zsol(1),wsisol,ier 
        call seva2d(bkx,lkx,bky,lky,c,rbdry(1),zbdry(1),pds,ier,n111)
        write (6,*) 'STEPS R,Z,Si,Err = ', rbdry(1),zbdry(1),pds(1),ier
        call seva2d(bkx,lkx,bky,lky,c,rbdry(nbdry),zbdry(nbdry),pds,ier &
             ,n111)
        write (6,*) 'STEPS R,Z,Si,Err = ', rbdry(nbdry),zbdry(nbdry) &
            ,pds(1),ier
#endif
      endif

#ifdef DEBUG_LEVEL2
      call seva2d(bkx,lkx,bky,lky,c,xout(1),yout(1),pds,ier,n111)
      write (6,*) 'STEPS simag,psibry,n1,si = ',simag,psibry,n111,pds(1)
      write (6,*) 'STEPS Z,R = ', zgrid(33),(rgrid(i),i=45,45)
      write (6,*) 'STEPS si = ',(psi((i-1)*nw+33),i=45,45)
      call seva2d(bkx,lkx,bky,lky,c,rgrid(45),zgrid(33),pds,ier,n111)
      write (6,*) 'STEPS R,Z,si = ', rgrid(45),zgrid(33),pds(1)
      write (6,*) 'STEPS lkx,lky = ',lkx,lky
#endif
!-----------------------------------------------------------------------
!--   get weighting function                                          --
!-----------------------------------------------------------------------
      call weight(rgrid,zgrid)
!-----------------------------------------------------------------------
!--   get response functions for MSE                                  --
!-----------------------------------------------------------------------
      if (kstark.gt.0.or.kdomse.gt.0) then
        do k=1,nstark
          if (rrgam(jtime,k).le.0.0) cycle
          call seva2d(bkx,lkx,bky,lky,c,rrgam(jtime,k) &
                      ,zzgam(jtime,k),pds,ier,n111)
          sistark(k)=pds(1)
          sisinow=(simag-pds(1))/sidif
          sigam(k)=sisinow
          fpnow=ffcurr(sisinow,kffcur)
          btgam(k)=fpnow*tmu/rrgam(jtime,k)
        enddo
      endif
!-----------------------------------------------------------------------
!--   get response functions for MSE-LS                               --
!-----------------------------------------------------------------------
      if (mmbmsels.gt.0.or.kdomsels.gt.0) then
        do k=1,nmsels
#ifdef DEBUG_LEVEL2
          write (6,*) 'STEPS MSE-LS k,rrmselt= ', k,rrmselt(jtime,k)
#endif
          if (rrmselt(jtime,k).le.0.0) cycle
          call seva2d(bkx,lkx,bky,lky,c,rrmselt(jtime,k) &
                      ,zzmselt(jtime,k),pds,ier,n333)
          simls(k)=pds(1)
          sisinow=(simag-pds(1))/sidif
          sinmls(k)=sisinow
          fpnow=ffcurr(sisinow,kffcur)
          btmls(k)=fpnow*tmu/rrmselt(jtime,k)
          brmls(k)=-pds(3)/rrmselt(jtime,k)
          bzmls(k)=pds(2)/rrmselt(jtime,k)
#ifdef DEBUG_LEVEL2
          write (6,*) 'STEPS MSE-LS k,rrmselt,br,bz,bt= ',  &
                      k,rrmselt(jtime,k),brmls(k),bzmls(k),btmls(k)
#endif
        enddo
      endif
!
      xouts(1:nfound)=1./xout(1:nfound)**2
      nzz=0
      call fluxav(xouts,xout,yout,nfound,psi,rgrid,nw,zgrid,nh, &
                  r2bdry,nzz,sdlobp,sdlbp)
      xouts(1:nfound)=1./xout(1:nfound)
      call fluxav(xouts,xout,yout,nfound,psi,rgrid,nw,zgrid,nh, &
                  r1bdry,nzz,sdlobp,sdlbp)
!-----------------------------------------------------------------------
!--   get metric elements at PSIWANT for edge constraint if needed    --
!-----------------------------------------------------------------------
      r1sdry(1)=r1bdry
      r2sdry(1)=r2bdry
      nnn=1
      if (abs(sizeroj(1)-1.0).gt.1.e-05_dp.or.kzeroj.ne.1) then
        if (kzeroj.gt.0) then
          do i=1,kzeroj
            if(sizeroj(i).ge.1.0) sizeroj(i)=0.99999_dp
            siwant=simag+sizeroj(i)*(psibry-simag)
            call surfac(siwant,psi,nw,nh,rgrid,zgrid,rsplt,zsplt,nfounc, &
                        npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nnn, &
                        rmaxis,zmaxis,negcur,kerror)
            if(kerror.gt.0) return
            csplt(1:nfounc)=1./rsplt(1:nfounc)**2
            call fluxav(csplt,rsplt,zsplt,nfounc,psi,rgrid,nw,zgrid,nh, &
                        r2sdry(i),nzz,sdlobp,sdlbp)
            csplt(1:nfounc)=1./rsplt(1:nfounc)
            call fluxav(csplt,rsplt,zsplt,nfounc,psi,rgrid,nw,zgrid,nh, &
                        r1sdry(i),nzz,sdlobp,sdlbp)
            r2surs = r2sdry(i)*sdlobp
            fpnow = ffcurr(psiwant,kffcur)
            fpnow = fpnow*tmu
          enddo
        endif
      endif
!-----------------------------------------------------------------------
!--   get metric elements at PSIWANT for q constraint if needed       --
!-----------------------------------------------------------------------
      if (nqwant.gt.0) then
        do i=1,nqwant
          siwant=simag+siwantq(i)*(psibry-simag)
          call surfac(siwant,psi,nw,nh,rgrid,zgrid,rsplt,zsplt,nfounc, &
                      npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nnn, &
                      rmaxis,zmaxis,negcur,kerror)
          if(kerror.gt.0) return
          csplt(1:nfounc)=1./rsplt(1:nfounc)**2
          call fluxav(csplt,rsplt,zsplt,nfounc,psi,rgrid,nw,zgrid,nh, &
                      r2qdry,nzz,sdlobp,sdlbp)
          csplt(1:nfounc)=1./rsplt(1:nfounc)
          call fluxav(csplt,rsplt,zsplt,nfounc,psi,rgrid,nw,zgrid,nh, &
                      r1qdry,nzz,sdlobp,sdlbp)
          r2surq = r2qdry*sdlobp
          fpnow = ffcurr(siwantq(i),kffcur)
          fpnow = fpnow*tmu
          qsiw(i)= abs(fpnow)/twopi*r2surq
          pasmsw(i)=sdlbp/tmu/twopi
        enddo
      endif
!
      cvolp(ixt)=0.0
      carea=0.0
      xym=xout(1)*yout(1)
      xyma=yout(1)
      do i=2,nfound
        xyp=xout(i)*yout(i)
        xypa=yout(i)
        dx=xout(i)-xout(i-1)
        cvolp(ixt)=cvolp(ixt)+twopi*(xyp+xym)/2.0*dx
        carea=carea+(xyma+xypa)/2.0*dx
        xym=xyp
        xyma=xypa
      enddo
      carea=abs(carea)
      if (iconvr.ne.3) then
        if((ix.gt.1).or.(ixout.le.1)) call chisqr(jtime)
      endif
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
!--   get beta and li for constraints                               --
!---------------------------------------------------------------------
      if (fli.gt.0.0.or.fbetan.gt.0.0) then
        call betsli(jtime,rgrid,zgrid,kerror)
        if(kerror.gt.0) return
      endif
!----------------------------------------------------------------------
!--   vertical stabilization information                             --
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
!--   magnetic axis parameters if needed                             --
!----------------------------------------------------------------------
      if (((icinit.gt.0).and.(iconvr.ne.3).and.(ixout.le.1)).or.(icurrt.eq.4).or. &
          (((icurrt.eq.2).or.(icurrt.eq.5)).and.(ixt.le.1).and.(icinit.gt.0))) then
        nqend=1
        n22=2
        if(errorm.lt.0.1_dp.and.icurrt.eq.4) nqend=nqiter
        do i=1,nqend
          if (i.gt.1) then
            call currnt(n22,jtime,n22,n22,kerror)
            if(kerror.gt.0) return
          endif

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
            if(qenp.gt.0.0) enp=enp*qenp/qmaxis
            if(qemp.gt.0.0) emp=emp*qmaxis/qemp
            enf=enp
            emf=emp
          endif
        enddo
        return
      endif
      select case (icurrt)
      case default ! 1
        rdiml=rmaxis/srma
        cjmaxi=cratio*(sbeta*rdiml+2.*salpha/rdiml)/darea
        if(kvtor.gt.0) &
          cjmaxi=cjmaxi+cratio/darea*sbetaw*rdiml*(rdiml**2-1.)
      case (2,5)
        fcentr=ffcurr(x000,kffcur)
        if (kvtor.eq.0) then
          cjmaxi=(rmaxis*ppcurr(x000,kppcur) &
                +fpcurr(x000,kffcur)/rmaxis)*cratio/darea
        else
          cjmaxi=(sum(rjjjx(1:kppcur)*brsp(nfcoil+1:nfcoil+kppcur)) &
                 +sum(rjjfx(1:kffcur)*brsp(nfcoil+kppcur+1:nfcoil+kppcur+kffcur)) &
                 +sum(rjjwx(1:kwwcur)*brsp(nfnpcr+1:nfnpcr+kwwcur))) &
                 /darea
        endif
      case (3)
        ! continue
      ! case (4) handled by preceeding if statement
      end select
!
      cqmaxi(ixt)=(emaxis**2+1.)*abs(fcentr)/twopi/emaxis &
                   /rmaxis**2/abs(cjmaxi)
      qmaxis=cqmaxi(ixt)
      return
 2100 format(/,1x,'shot',i6,' at ',i6,' ms ','** Problem in BOUND **')
11001 format(/,1x,'** 2nd seperatrix **',2x,e10.3,2x,i4)
      end subroutine steps

!**********************************************************************
!>
!!    inicur initializes the current density distribution
!!
!!
!!    @param ks : time index
!!
!**********************************************************************
      subroutine inicur(ks)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      save isicinit,zelips
      character(14) :: sfile
!
      pcurrt=0.0
      if(ivacum.gt.0) return
!----------------------------------------------------------------------
!--   icinit=1 uniform and elliptical flux surfaces                  --
!--          2 parabolic and elliptical                              --
!--         -2 use ESAVE.DAT at first slice and previous slice later --
!--        -12 parabolic and elliptical first slice and previous     --
!--            slice subsequently                                    --
!----------------------------------------------------------------------
      if(ks.eq.1) isicinit=icinit
      if (isicinit.lt.0) then
        if (ks.gt.1) then
          icinit=isicinit
          return
        endif
        if(isicinit.eq.-12) icinit=2
      endif
      select case (abs(icinit))
      case (1)
        if (aelip.le.0.0) &
          aelip=min(0.50_dp, &
                minval(sqrt((xlim(1:limitr)-relip)**2 &
                          +((ylim(1:limitr)-zelip)/eelip)**2)))
        delcur=1.
        sumi=0.0
        do i=1,nw
          do j=1,nh
            kk=(i-1)*nh+j
            pcurrt(kk)=delcur*zero(kk)/rgrid(i)
            sumi=sumi+pcurrt(kk)
          enddo
        enddo
        cratio=pasmat(ks)/sumi
        pcurrt(1:nwnh)=pcurrt(1:nwnh)*cratio*zero(1:nwnh)
!
      case (2)
        if (icinit.le.0) then
          if (nproc.gt.1) then
            WRITE(sfile,fmt='(i5.5)') rank
            sfile='esave'//TRIM(sfile)//'.dat'
          else
            sfile='esave.dat'
          endif
          open(unit=nsave,form='unformatted',file=sfile, &
               status='old',iostat=ioerr)
          if (ioerr.eq.0) then 
            read (nsave) mw,mh
            read (nsave) xpsi
            read (nsave) brsp
            read (nsave) www
            read (nsave) emaxis, rmaxis, fcentr
            close(unit=nsave)
          endif
        else
          if(ks.eq.1) zelips=zelip
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
!--       set zelip=0.0 if bad signals              96/06/24   --
!----------------------------------------------------------------
          if (size(fwtmp2).ge.37) then
            if(abs(fwtmp2(37)).le.1.0e-30_dp) zelip=0.0
          endif
          if (size(fwtmp2).ge.43) then
            if(abs(fwtmp2(43)).le.1.0e-30_dp) zelip=0.0
          endif
          if (size(fwtmp2).ge.57) then
            if(abs(fwtmp2(57)).le.1.0e-30_dp) zelip=0.0
          endif
          if (size(fwtmp2).ge.53) then
            if(abs(fwtmp2(53)).le.1.0e-30_dp) zelip=0.0
          endif
          if (size(fwtsi).ge.27) then
            if(abs(fwtsi(27)).le.1.0e-30_dp)  zelip=0.0
          endif
          if (size(fwtsi).ge.37) then
            if(abs(fwtsi(37)).le.1.0e-30_dp)  zelip=0.0
          endif
          if (size(fwtsi).ge.2) then
            if(abs(fwtsi(2)).le.1.0e-30_dp)   zelip=0.0
          endif
          if (size(fwtsi).ge.11) then
            if(abs(fwtsi(11)).le.1.0e-30_dp)  zelip=0.0
          endif
!
          do i=1,nw
            do j=1,nh
              kk=(i-1)*nh+j
              erho=sqrt((rgrid(i)-relip)**2+((zgrid(j)-zelip)/eelip)**2)
              xpsi(kk)=(erho/aelip)**2
            enddo
          enddo
        endif
      end select
      return
      end subroutine inicur

!**********************************************************************
!>
!!    vescur computes the currents induced in the vessel
!!    segments due to E coils and F coils
!!
!!
!!    @param jtime : time index
!!
!**********************************************************************
      subroutine vescur(jtime)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!
      if(ivesel.eq.5) return
!
      if(ivesel.eq.1) return
!
      vcurrt(1:nvesel)=vloopt(jtime)/rsisvs(1:nvesel)
      return
      end subroutine vescur
