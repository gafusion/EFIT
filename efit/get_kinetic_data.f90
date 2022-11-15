!**********************************************************************
!>
!!    GETNE gets the electron density
!!    profile
!!
!!    @param jtime : time index
!! 
!!    @param kerror : error flag
!!
!**********************************************************************
      subroutine getne(jtime,kerror)
      use commonblocks,only: c,bkx,bky
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      integer*4, intent(in) :: jtime
      integer*4, intent(out) :: kerror
      dimension pds(6)

      kerror = 0

      ne_constraint: if (npnef.ne.0) then
!----------------------------------------------------------------------
!--   singular decomposition                                         --
!----------------------------------------------------------------------
      do nj=1,npress
        do nk=1,npnef
          xn=-rpress(nj)
          arsp_cw2(nj,nk)=xn**(nk-1)/sgneth(nj)
        enddo
        bdata_cw2(nj)=dnethom(nj)/sgneth(nj)
      enddo
      nnedat=npress
      if (cstabne.gt.0.0) then
        nj=npress
        do jj=ncstne,npnef
          nj=nj+1
          do nk=1,npnef
            arsp_cw2(nj,nk)=0.0
            if(jj.eq.nk) arsp_cw2(nj,nk)=cstabne
          enddo
          bdata_cw2(nj)=0.0
        enddo
        nnedat=nnedat+npnef-ncstne+1
      endif
!---------------------------------------------------------------------
!--   form error matrix                                             --
!---------------------------------------------------------------------
      do i=1,npnef
        do j=1,npnef
          ematrix_cw2(i,j)=0.0
          do k=1,npress
            ematrix_cw2(i,j)=ematrix_cw2(i,j)+arsp_cw2(k,i)*arsp_cw2(k,j)
          enddo
        enddo
      enddo

      nnn=1
      call sdecm(arsp_cw2,ndata,nnedat,npnef,bdata_cw2,nnedat,nnn,wrsp_cw2,work_cw2,ier)
      if (ier.eq.129) then
        kerror = 1
        call errctrl_msg('getne','sdecm failed to converge')
        return
      end if

      cond=ier
      toler=1.0e-06*wrsp_cw2(1)
      do i=1,npnef
        t=0.0
        if (wrsp_cw2(i).gt.toler) t=bdata_cw2(i)/wrsp_cw2(i)
        work_cw2(i)=t
      enddo
      do i=1,npnef
        defit(i)=0.0
        do j=1,npnef
          defit(i)=defit(i)+arsp_cw2(i,j)*work_cw2(j)
        enddo
      enddo
!------------------------------------------------------------------
!--   compute chi square                                         --
!------------------------------------------------------------------
      chisqne=0.0
      do i=1,npress
        denow=0.0
        xn=-rpress(i)
        do j=1,npnef
          denow=denow+defit(j)*xn**(j-1)
        enddo
        chisqne=chisqne+((denow-dnethom(i))/sgneth(i))**2
      enddo
!---------------------------------------------------------------------
!--   get inverse of error matrix                                   --
!---------------------------------------------------------------------
      n44=4
      call linv1f(ematrix_cw2,npnef,nppcur,einv_cw2,n44,work_cw2,ier)
!----------------------------------------------------------------------
!--   boundary values                                                --
!----------------------------------------------------------------------
      debdry=0.0
      sdebdry=0.0
      depbry=0.0
      sigdepb=0.0
      do j=1,npnef
        do i=1,npnef
          sdebdry=sdebdry+einv_cw2(i,j)
          sigdepb=sigdepb+(i-1)*(j-1)*einv_cw2(i,j)
        enddo
        depbry=depbry+(j-1)*defit(j)
        debdry=debdry+defit(j)
      enddo
      if (sdebdry.gt.0.0) then
        sdebdry=sqrt(sdebdry)
      else
        sdebdry=debdry
      endif
      if (sigdepb.gt.0.0) then
        sigdepb=sqrt(sigdepb)
      else
        sigdepb=depbry
      endif
!----------------------------------------------------------------------
!--   get co2 v2 chord normalization factor                          --
!----------------------------------------------------------------------
      call lenco2(xout,yout,nfound,jtime)
      delz=(zuperts(jtime)-zlowerts)/(nh-1)*0.01_dp
      dneco2=debdry
      do i=2,nh-1
        znow=zlowerts*0.01_dp+delz*(i-1)
         call seva2d(bkx,lkx,bky,lky,c,rmajts,znow,pds,ier,n111)
        denow=0.0
        xn=(simag-pds(1))/sidif
        do j=1,npnef
          denow=denow+defit(j)*xn**(j-1)
        enddo
        dneco2=dneco2+denow
      enddo
      dneco2=dneco2/(nh-1)
      fco2now=dco2v(jtime,2)*1.e-13/dneco2
      fco2ne=fco2now*fco2ne
      do i=1,npnef
        defit(i)=defit(i)*fco2now
      enddo
      do i=1,npneth
        dnethom(i)=dnethom(i)*fco2now
        sgneth(i)=sgneth(i)*fco2now
      enddo
      debdry=debdry*fco2now
      sdebdry=sdebdry*fco2now
      dibdry=debdry
      dipbry=depbry
      sigdipb=sigdepb
      sdibdry=sdebdry
!
      endif ne_constraint
!-----------------------------------------------------------------------
!--   get ion density profile from zeff                               --
!-----------------------------------------------------------------------
      factor=(1.+zlowimp-zeffvs)/zlowimp
      do i=1,npress
        dnitho(i)=factor*dnethom(i)
!------------------------------------------------------------------------
!--     correct for dilution factor due to beam                        --
!------------------------------------------------------------------------
        dnitho(i)=dnitho(i)-dnbthom(i)
        snitho(i)=factor*sgneth(i)
      enddo
!
      return
      end subroutine getne

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
            if(jj.eq.nk) arsp_cw2(nj,nk)=cstabte
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
      use commonblocks,only: c,bkx,bky
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
