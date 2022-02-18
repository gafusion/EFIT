#include "config.f"
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
      use commonblocks,only: cjrf,c,wk,copy,bkx,bky
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
#if defined(USEMPI)
      include 'mpif.h'
#endif
      integer*4 :: ier_all
      integer*4, intent(inout) :: kerror
      dimension pds(6)

      kerror = 0

      if (npnef.ne.0) then
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
          if (jj.eq.nk) arsp_cw2(nj,nk)=cstabne
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
!
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
      factor=1.0 ! TODO: this appears to be unset (should it be removed?)
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
      dibdry=factor*debdry
      dipbry=factor*depbry
      sigdipb=factor*sigdepb
      sdibdry=factor*sdebdry
!
      endif ! npnef.ne.0
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
      end
