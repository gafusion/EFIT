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
      integer :: ier_all
      integer, intent(inout) :: kerror
      dimension pds(6)

      kerror = 0

      if (npnef.eq.0) go to 1950
!----------------------------------------------------------------------
!--  singular decomposition                                          --
!----------------------------------------------------------------------
      do 1100 nj=1,npress
        do 1000 nk=1,npnef
          xn=-rpress(nj)
          arsp_cw2(nj,nk)=xn**(nk-1)/sgneth(nj)
 1000   continue
        bdata_cw2(nj)=dnethom(nj)/sgneth(nj)
 1100 continue
      nnedat=npress
      if (cstabne.gt.0.0) then
        nj=npress
        do 22200 jj=ncstne,npnef
         nj=nj+1
         do 22190 nk=1,npnef
          arsp_cw2(nj,nk)=0.0
          if (jj.eq.nk) arsp_cw2(nj,nk)=cstabne
22190    continue
         bdata_cw2(nj)=0.0
22200   continue
        nnedat=nnedat+npnef-ncstne+1
      endif
!---------------------------------------------------------------------
!-- form error matrix                                               --
!---------------------------------------------------------------------
      do 1120 i=1,npnef
      do 1120 j=1,npnef
        ematrix_cw2(i,j)=0.0
        do 1120 k=1,npress
          ematrix_cw2(i,j)=ematrix_cw2(i,j)+arsp_cw2(k,i)*arsp_cw2(k,j)
 1120   continue
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
      do 1600 i=1,npnef
        t=0.0
        if (wrsp_cw2(i).gt.toler) t=bdata_cw2(i)/wrsp_cw2(i)
        work_cw2(i)=t
 1600 continue
      do 1650 i=1,npnef
        defit(i)=0.0
        do 1650 j=1,npnef
          defit(i)=defit(i)+arsp_cw2(i,j)*work_cw2(j)
 1650   continue
!------------------------------------------------------------------
!-- compute chi square                                           --
!------------------------------------------------------------------
      chisqne=0.0
      do 1700 i=1,npress
        denow=0.0
        xn=-rpress(i)
        do 1670 j=1,npnef
 1670     denow=denow+defit(j)*xn**(j-1)
        chisqne=chisqne+((denow-dnethom(i))/sgneth(i))**2
 1700 continue
!---------------------------------------------------------------------
!-- get inverse of error matrix                                     --
!---------------------------------------------------------------------
      n44=4
      call linv1f(ematrix_cw2,npnef,nppcur,einv_cw2,n44,work_cw2,ier)
!----------------------------------------------------------------------
!--  boundary values                                                 --
!----------------------------------------------------------------------
        debdry=0.0
        sdebdry=0.0
        depbry=0.0
        sigdepb=0.0
        do 1850 j=1,npnef
          do 22300 i=1,npnef
            sdebdry=sdebdry+einv_cw2(i,j)
            sigdepb=sigdepb+(i-1)*(j-1)*einv_cw2(i,j)
22300     continue
          depbry=depbry+(j-1)*defit(j)
 1850     debdry=debdry+defit(j)
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
!-- get co2 v2 chord normalization factor                             -
!----------------------------------------------------------------------
      call lenco2(xout,yout,nfound,jtime)
      delz=(zuperts(jtime)-zlowerts)/(nh-1)*0.01_dp
      dneco2=debdry
      do 1900 i=2,nh-1
        znow=zlowerts*0.01_dp+delz*(i-1)
         call seva2d(bkx,lkx,bky,lky,c,rmajts,znow,pds,ier,n111)
        denow=0.0
        xn=(simag-pds(1))/sidif
        do 1870 j=1,npnef
 1870     denow=denow+defit(j)*xn**(j-1)
        dneco2=dneco2+denow
 1900 continue
      dneco2=dneco2/(nh-1)
      fco2now=dco2v(jtime,2)*1.e-13/dneco2
      fco2ne=fco2now*fco2ne
      do 22340 i=1,npnef
        defit(i)=defit(i)*fco2now
22340 continue
      do 22350 i=1,npneth
        dnethom(i)=dnethom(i)*fco2now
        sgneth(i)=sgneth(i)*fco2now
22350 continue
      debdry=debdry*fco2now
      sdebdry=sdebdry*fco2now
      dibdry=factor*debdry
      dipbry=factor*depbry
      sigdipb=factor*sigdepb
      sdibdry=factor*sdebdry
!
 1950 continue
!-----------------------------------------------------------------------
!-- get ion density profile from zeff                                 --
!-----------------------------------------------------------------------
      factor=(1.+zlowimp-zeffvs)/zlowimp
      do 2000 i=1,npress
        dnitho(i)=factor*dnethom(i)
!------------------------------------------------------------------------
!--  correct for dilution factor due to beam                           --
!------------------------------------------------------------------------
        dnitho(i)=dnitho(i)-dnbthom(i)
        snitho(i)=factor*sgneth(i)
 2000 continue
!
      return
      end
