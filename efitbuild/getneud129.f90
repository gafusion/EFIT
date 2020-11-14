  subroutine getne(jtime,kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          GETNE gets the electron density                         **
!**          profile.                                                **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          15/09/87..........first created                         **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: cjrf,c,wk,copy,bkx,bky
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
! MPI >>>
#if defined(USEMPI)
      include 'mpif.h'
#endif
      integer :: ier_all
      integer, intent(inout) :: kerror
! MPI <<<
      common/cwork2/arsp(ndata,nppcur),wrsp(nppcur),work(ndata), &
                    bdata(ndata),ematrix(nppcur,nppcur), &
                    einv(nppcur,nppcur)
      common/cwork3/lkx,lky
      dimension pds(6)

      kerror = 0

      if (npnef.eq.0) go to 1950
!----------------------------------------------------------------------
!--  singular decomposition                                          --
!----------------------------------------------------------------------
      do 1100 nj=1,npress
        do 1000 nk=1,npnef
          xn=-rpress(nj)
          arsp(nj,nk)=xn**(nk-1)/sgneth(nj)
 1000   continue
        bdata(nj)=dnethom(nj)/sgneth(nj)
 1100 continue
      nnedat=npress
      if (cstabne.gt.0.0) then
        nj=npress
        do 22200 jj=ncstne,npnef
         nj=nj+1
         do 22190 nk=1,npnef
          arsp(nj,nk)=0.0
          if (jj.eq.nk) arsp(nj,nk)=cstabne
22190    continue
         bdata(nj)=0.0
22200   continue
        nnedat=nnedat+npnef-ncstne+1
      endif
!---------------------------------------------------------------------
!-- form error matrix                                               --
!---------------------------------------------------------------------
      do 1120 i=1,npnef
      do 1120 j=1,npnef
        ematrix(i,j)=0.0
        do 1120 k=1,npress
          ematrix(i,j)=ematrix(i,j)+arsp(k,i)*arsp(k,j)
 1120   continue
!
      nnn=1
      call sdecm(arsp,ndata,nnedat,npnef,bdata,nnedat,nnn,wrsp,work,ier)
      if (ier.eq.129) then
        kerror = 1
        call errctrl_msg('getne','sdecm failed to converge')
        return
      end if
      cond=ier
      toler=1.0e-06*wrsp(1)
      do 1600 i=1,npnef
        t=0.0
        if (wrsp(i).gt.toler) t=bdata(i)/wrsp(i)
        work(i)=t
 1600 continue
      do 1650 i=1,npnef
        defit(i)=0.0
        do 1650 j=1,npnef
          defit(i)=defit(i)+arsp(i,j)*work(j)
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
      call linv1f(ematrix,npnef,nppcur,einv,n44,work,ier)
!----------------------------------------------------------------------
!--  boundary values                                                 --
!----------------------------------------------------------------------
        debdry=0.0
        sdebdry=0.0
        depbry=0.0
        sigdepb=0.0
        do 1850 j=1,npnef
          do 22300 i=1,npnef
            sdebdry=sdebdry+einv(i,j)
            sigdepb=sigdepb+(i-1)*(j-1)*einv(i,j)
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
