      subroutine matrix(jtime,iter,ichisq,nniter,kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          matrix calculates the appropriate response matrix and   **
!**          invert it to get the plasma current strengths.  note    **
!**          that npcurn=nffcur+nppcur, nrsmat=nfcoil+npcurn+        **
!**          number of constraints.                                  **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**          24/07/85..........revised                               **
!**          94/03/08..........revised                               **
!**        2015/03/27..........revised MSE-LS                        **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      include 'basiscomdu.f90'
      dimension arsp(nrsmat,mfnpcr),wrsp(mfnpcr)
      dimension brsold(nrsmat),work(nrsma2),vcurrto(nvesel)
      dimension xpsfp(nffcur),xpspp(nppcur),xpspwp(nwwcur)
      dimension crsp(4*(npcurn-2)+6+npcurn*npcurn,nrsmat)
      dimension b(nrsmat),z(4*(npcurn-2)+6+npcurn*npcurn)
      common/cwork3/lkx,lky
      dimension pds(6)
      dimension rxxx(ndata),rxxxf(ndata),rxx2(ndata),rxxw(ndata)
!---------------------------------------------------------------------
!--   relax saimin=50 from 30               04/27/90                --
!--                60 from 50               03/31/93                --
!---------------------------------------------------------------------
      data iupdat/0/,minite/8/,ten24/1.e4/,z04/1.0e-04/
      save z04
!
      integer, intent(inout) :: kerror
      kerror = 0
      if (iconvr.eq.3) go to 6000
!----------------------------------------------------------------------
!-- Variable fitdelz                                                 --
!----------------------------------------------------------------------
      if (fitdelz) scadelz=scaledz
!----------------------------------------------------------------------
!--  set up fitting weight for boundary constraints                  --
!----------------------------------------------------------------------
      if (nbdry.gt.0) then
        fwtbdr=abs(errbry)*max(abs(sidif),z04)
        fwtbdr=1.0/fwtbdr
        do i=1,nbdry
          if (sigrbd(i).lt.1.e10.and.sigzbd(i).lt.1.e10) then
          call seva2d(bkx,lkx,bky,lky,c,rbdry(i),zbdry(i),pds,ier,n666)
          fwtbdr=sqrt((sigrbd(i)*pds(2))**2+(sigzbd(i)*pds(3))**2)
          fwtbdr=1.0/fwtbdr
          endif
          fwtbry(i)=fwtbdr*fwtbdry(i)
        enddo
      endif
!----------------------------------------------------------------------
!-- set up the response matrix arsp                                  --
!----------------------------------------------------------------------
      ichisq=0
      nsq=2
      saiold=saisq
      do 1000 i=1,nrsmat
        brsold(i)=brsp(i)
 1000 continue
      if (ifitvs.gt.0) then
        do 1010 i=1,nvesel
          vcurrto(i)=vcurrt(i)
 1010   continue
      endif
!----------------------------------------------------------------------
!--  singular decomposition, first F-coil currents, set up arsp      --
!----------------------------------------------------------------------
      do 2100 nk=1,nfcoil
        nj=0
        do 2020 m=1,nsilop
          if (fwtsi(m).le.0.0) go to 2020
          nj=nj+1
          arsp(nj,nk)=fwtsi(m)*rsilfc(m,nk)
 2020   continue
        do 2040 m=1,magpri
          if (fwtmp2(m).le.0.0) go to 2040
          nj=nj+1
          arsp(nj,nk)=fwtmp2(m)*rmp2fc(m,nk)
 2040   continue
        do 2060 m=1,nstark
          if (fwtgam(m).le.0.0) go to 2060
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgamfc(m,nk)
 2060   continue
        do 92060 m=1,nmsels
          if (fwtbmselt(jtime,m).le.0.0) go to 92060
          nj=nj+1
          arsp(nj,nk)=0.0
92060   continue
        do 92070 m=1,nmsels
          if (fwtemselt(jtime,m).le.0.0) go to 92070
          nj=nj+1
          arsp(nj,nk)=0.0
92070   continue
        do 2070 m=1,nece
          if (fwtece(m).le.0.0) go to 2070
          nj=nj+1
          arsp(nj,nk)=fwtece(m)*recefc(m,nk)
 2070   continue
          if (fwtecebz.le.0.0) go to 2075
          nj=nj+1
          arsp(nj,nk)=fwtecebz*recebzfc(nk)
 2075   continue
        if (fwtcur.le.0.0) go to 2080
        nj=nj+1
        arsp(nj,nk)=0.0
 2080   continue
        if (fwtqa.le.0.0) go to 2085
        nj=nj+1
        arsp(nj,nk)=0.0
 2085 continue
      if (fwtbp.le.0.0) go to 2090
      do 2087 m=2,kffcur
        nj=nj+1
        arsp(nj,nk)=0.0
 2087 continue
 2090 continue
      if (fwtdlc.le.0.0) go to 2092
      nj=nj+1
      arsp(nj,nk)=0.0
 2092 continue
 2096 continue
      do 2097 m=1,nfcoil
        if (fwtfc(m).le.0.0) go to 2097
        nj=nj+1
        arsp(nj,nk)=0.
        if (m.eq.nk) arsp(nj,nk)=fwtfc(m)
 2097 continue
!--------------------------------------------------------------------------
!--  pressure data                                                       --
!--------------------------------------------------------------------------
      if (kprfit.le.0.or.kdofit.eq.0) go to 2099
      if (npress.le.0) goto 2099
      do 2098 m=1,npress
        if (fwtpre(m).le.0.0) goto 2098
        nj=nj+1
        arsp(nj,nk)=0.0
 2098 continue
 2099 continue
      if (kpressb.eq.2) then
        nj=nj+1
        arsp(nj,nk)=0.0
      endif
!---------------------------------------------------------------------------
!-- rotational pressure data                                              --
!---------------------------------------------------------------------------
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          if (fwtprw(m).gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
          endif
        enddo
      endif
!---------------------------------------------------------------------------
!--  J(PSIWANT) constraint                                                --
!---------------------------------------------------------------------------
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
        arsp(nj,nk)=0.0
       enddo
      endif
!----------------------------------------------------------------------------
!-- P', FF', and rotational constraints                                    --
!----------------------------------------------------------------------------
      if (kcalpa.gt.0) then
        do 30050 j=1,kcalpa
          nj=nj+1
          arsp(nj,nk)=0.0
30050   continue
      endif
      if (kcgama.gt.0) then
        do 30060 j=1,kcgama
          nj=nj+1
          arsp(nj,nk)=0.0
30060   continue
      endif
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
!------------------------------------------------------------------------------
!-- Boundary constraints, data                                               --
!------------------------------------------------------------------------------
      if (nbdry.gt.0) then
        do j=1,nbdry
          if (fwtbdry(j).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=fwtbry(j)*rbdrfc(j,nk)
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!--  E coil currents, data                                                   --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
       do m=1,nesum
        if (fwtec(m).gt.0.0) then
        nj=nj+1
        arsp(nj,nk)=0.
        endif
       enddo
      endif
!------------------------------------------------------------------------------
!--  Reference flux, data                                                    --
!------------------------------------------------------------------------------
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
      endif
!------------------------------------------------------------------------
!-- Summation of F-coils currents
!------------------------------------------------------------------------
      if (fitfcsum) then
         nj = nj + 1
         arsp(nj,nk) = fwtfcsum(nk)
      endif
 2100 continue
!-----------------------------------------------------------------------
!--  plasma current P', FF', and Pw', set up response matrix arsp     --
!-----------------------------------------------------------------------
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        ysiwant=sizeroj(i)
        call setfp(ysiwant,xpsfp)
        call setpp(ysiwant,xpspp)
!----------------------------------------------------------------------
!--  local or flux surface averaged J constraint ?                   --
!----------------------------------------------------------------------
        if (rzeroj(i).gt.0.0) rxxx(i)=rzeroj(i)
        if (rzeroj(i).lt.0.0) rxxx(i)=rseps(1,jtime)/100.
        rxxxf(i)=rxxx(i)
        if (rzeroj(i).eq.0.0) then
             rxxx(i)=1./r1sdry(i)
             rxxxf(i)=r1sdry(i)/r2sdry(i)
        endif
        if (rxxx(i).le.0.0) then
            rxxx(i)=rcentr
            rxxxf(i)=rxxx(i)
        endif
        brspmin=max(ten24,abs(brsp(nfcoil+1)))
        fwtxxzj=fwtxxj*1000./brspmin
        if (kvtor.gt.0) then
          call setpwp(ysiwant,xpspwp)
          if (nniter.gt.0) then
          if (kvtor.eq.2.or.kvtor.eq.3) then
            prew0=pwcurr(ysiwant,kwwcur)
            pres0=prcurr(ysiwant,kppcur)
            if (abs(pres0).gt.1.e-10) then
               pwop0=prew0/pres0
            else
               pwop0=0.0
            endif
          endif
          endif
          if (rzeroj(i).gt.0.0) then
              rxx2(i)=(rzeroj(i)/rvtor)**2-1.
              rxxw(i)=rxx2(i)*rzeroj(i)
              if (nniter.gt.0) then
              if (kvtor.eq.2) then
                rxxw(i)=rxxw(i)*(1.+pwop0*rxx2(i))
                rxxx(i)=rxxx(i)*(1.-0.5*(pwop0*rxx2(i))**2)
              endif
              if (kvtor.eq.3) then
                pwp0r2=pwop0*rxx2(i)
                ptop0=exp(pwp0r2)
                rxxw(i)=rxxw(i)*ptop0
                rxxx(i)=rxxx(i)*ptop0*(1.-pwp0r2)
              endif
              endif
          endif
          if (rzeroj(i).lt.0.0) then
              rxxw(i)=rseps(1,jtime)/100.
              rxx2(i)=(rxxw(i)/rvtor)**2-1.
              rxxw(i)=rxx2(i)*rxxw(i)
              if (nniter.gt.0) then
              if (kvtor.eq.2) then
                rxxw(i)=rxxw(i)*(1.+pwop0*rxx2(i))
                rxxx(i)=rxxx(i)*(1.-0.5*(pwop0*rxx2(i))**2)
              endif
              if (kvtor.eq.3) then
                pwp0r2=pwop0*rxx2(i)
                ptop0=exp(pwp0r2)
                rxxw(i)=rxxw(i)*ptop0
                rxxx(i)=rxxx(i)*ptop0*(1.-pwp0r2)
              endif
              endif
          endif
          if (rzeroj(i).eq.0.0) then
              rxx2(i)=  r2wdry/rvtor**2-1.
              rxxw(i)=  rxx2(i)/r1sdry(i)
              if (nniter.gt.0) then
              if (kvtor.eq.2) then
                rxxw(i)=rxxw(i)+pwop0*r4wdry/r1sdry(i)
                rxxx(i)=rxxx(i)-0.5*pwop0**2*r4wdry/r1sdry(i)
              endif
              if (kvtor.eq.3) then
                rxxx(i)=(rpwdry-pwop0*rp2wdry)/r1sdry(i)
                rxxw(i)=rp2wdry/r1sdry(i)
              endif
              endif
          endif
        endif
       enddo
      endif
!----------------------------------------------------------------------
!--  start loop for plasma fitting parameters: P', FF', Pw'          --
!----------------------------------------------------------------------
      do 2210 nk=nfcoil+1,nfnwcr
        n=nk-nfcoil
        nj=0
        do 2120 m=1,nsilop
          if (fwtsi(m).le.0.0) go to 2120
          nj=nj+1
          arsp(nj,nk)=fwtsi(m)*rsilpc(m,n)
 2120   continue
        do 2140 m=1,magpri
          if (fwtmp2(m).le.0.0) go to 2140
          nj=nj+1
          arsp(nj,nk)=fwtmp2(m)*rmp2pc(m,n)
 2140   continue
        do 2160 m=1,nstark
          if (fwtgam(m).le.0.0) go to 2160
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgampc(m,n)
 2160   continue
        do 92160 m=1,nmsels
          if (fwtbmselt(jtime,m).le.0.0) go to 92160
          nj=nj+1
          arsp(nj,nk)=fwtbmselt(jtime,m)*rmlspc(m,n)
92160   continue
        if (jdebug.eq.'MSEL') then
             m=3
             write (6,*) 'MATRIX jtime, m,n,fwtbmselt,rmlspc = ',jtime, &
                  m,n,fwtbmselt(jtime,m),rmlspc(m,n)
             write (6,*) 'nj,nk,arsp= ',nj,nk,arsp(nj,nk)
        endif
        do 92170 m=1,nmsels
          if (fwtemselt(jtime,m).le.0.0) go to 92170
          nj=nj+1
          arsp(nj,nk)=0.0
92170   continue
        do 2170 m=1,nece
          if (fwtece(m).le.0.0) go to 2170
          nj=nj+1
          arsp(nj,nk)=fwtece(m)*recepc(m,n)
 2170   continue
          if (fwtecebz.le.0.0) go to 2175
          nj=nj+1
          arsp(nj,nk)=fwtecebz*recebzpc(n)
 2175   continue
        if (fwtcur.le.0.0) go to 2180
        nj=nj+1
        arsp(nj,nk)=fwtcur*fgowpc(n)
 2180   continue
        if (fwtqa.le.0.0) go to 2190
        nj=nj+1
        arsp(nj,nk)=0.0
 2190 continue
      if (fwtbp.le.0.0) go to 2194
      do 2192 m=2,kffcur
        nj=nj+1
        arsp(nj,nk)=0.0
 2192 continue
 2194 continue
      if (fwtdlc.le.0.0) go to 2196
      nj=nj+1
      arsp(nj,nk)=0.0
 2196 continue
 2202 continue
      do 2203 m=1,nfcoil
        if (fwtfc(m).le.0.0) go to 2203
        nj=nj+1
        arsp(nj,nk)=0.
 2203 continue
!--------------------------------------------------------------------
!--  pressure                                                      --
!--------------------------------------------------------------------
      if (kprfit.le.0.or.kdofit.eq.0) go to 2206
      if (npress.le.0) goto 2206
      do 2204 m=1,npress
        if (fwtpre(m).le.0.0) goto 2204
        nj=nj+1
        arsp(nj,nk)=0.0
        if (n.le.kppcur) arsp(nj,nk)=rprepc(m,n)/sigpre(m)*fwtpre(m)
 2204 continue
 2206 continue
!--------------------------------------------------------------------
!-- P'(1)                                                          --
!--------------------------------------------------------------------
      if (kpressb.eq.2) then
        nj=nj+1
        if (n.le.kppcur) arsp(nj,nk)=1./sigppb/darea
      endif
!---------------------------------------------------------------------
!-- rotational pressure                                             --
!---------------------------------------------------------------------
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          if (fwtprw(m).gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
          if (n.gt.kpcurn) then
            nmw=n-kpcurn
            arsp(nj,nk)=rprwpc(m,nmw)/sigprw(m)*fwtprw(m)
          endif
          endif
        enddo
      endif
!----------------------------------------------------------------------
!--  J(PSIWANT) constraint                                           --
!----------------------------------------------------------------------
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
        ysiwant=sizeroj(i)
        if (n.le.kppcur) then
            call setpp(ysiwant,xpspp)
            xjj=xpspp(n)
            arsp(nj,nk)=rxxx(i)*fwtxxzj*xjj
        elseif (n.le.kpcurn) then
            call setfp(ysiwant,xpsfp)
            xjj=xpsfp(n-kppcur)
            arsp(nj,nk)=fwtxxzj/rxxxf(i)*xjj
        elseif (kvtor.gt.0) then
            xjj=xpspwp(n-kpcurn)
            arsp(nj,nk)=rxxw(i)*fwtxxzj*xjj
        endif
       enddo
      endif
!-------------------------------------------------------------------------
!--  p' and ff' constraints                                             --
!-------------------------------------------------------------------------
      if (kcalpa.gt.0) then
        brspmin=max(ten24,abs(brsp(nfcoil+1)))
        fwtxxa=fwtxx*1000./brspmin
        do 30140 j=1,kcalpa
          nj=nj+1
          if (n.le.kppcur) then
            arsp(nj,nk)=calpa(n,j)*fwtxxa
          else
            arsp(nj,nk)=0.0
          endif
30140   continue
      endif
      if (kcgama.gt.0) then
        brspmin=max(ten24,abs(brsp(nfcoil+1)))
        fwtxxf=fwtxx*1000./brspmin
        do 30150 j=1,kcgama
          nj=nj+1
          if (n.le.kppcur) then
            arsp(nj,nk)=0.0
          elseif (n.le.kpcurn) then
            arsp(nj,nk)=cgama(n-kppcur,j)*fwtxxf
          else
            arsp(nj,nk)=0.0
          endif
30150   continue
      endif
!---------------------------------------------------------------------
!-- rotational constraints                                          --
!---------------------------------------------------------------------
      if (kcomega.gt.0) then
        brspmin=max(ten24,abs(brsp(nfcoil+1)))
        fwtxxo=fwtxx*1000./brspmin
        do j=1,kcomega
          nj=nj+1
          if (n.le.kpcurn) then
            arsp(nj,nk)=0.0
          else
            arsp(nj,nk)=comega(n,j)*fwtxxo
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!-- Boundary constraints due to plasma contributions                         --
!------------------------------------------------------------------------------
      if (nbdry.gt.0) then
        do j=1,nbdry
          if (fwtbdry(j).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=fwtbry(j)*gbdrpc(j,n)
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!--  E coil currents                                                         --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
       do m=1,nesum
        if (fwtec(m).gt.0.0) then
        nj=nj+1
        arsp(nj,nk)=0.
        endif
       enddo
      endif
!------------------------------------------------------------------------------
!--  Reference flux                                                          --
!------------------------------------------------------------------------------
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
      endif
!-------------------------------------------------------------------------
!-- Summation of F-coils currents
!-------------------------------------------------------------------------
      if (fitfcsum) then
         nj = nj + 1
         arsp(nj,nk) = 0.0
      endif
 2210 continue
      need=nfnwcr
!----------------------------------------------------------------------
!-- fit vessel currents                                              --
!----------------------------------------------------------------------
      if (ifitvs.le.0) go to 2310
      if (nfourier.gt.1) then
       need=need+nfourier*2+1
      else
       need=need+nvesel
      endif
      do 2300 nk=nfnwcr+1,need
        mk=nk-nfnwcr
        nj=0
        do 2220 m=1,nsilop
          if (fwtsi(m).le.0.0) go to 2220
          nj=nj+1
          if(nfourier.gt.1)then
          temp=0.
          do i=1,nvesel
          temp=temp+rsilvs(m,i)*vecta(mk,i)
          enddo
           arsp(nj,nk)=fwtsi(m)*temp
          else
          arsp(nj,nk)=fwtsi(m)*rsilvs(m,mk)
          endif
 2220   continue
        do 2240 m=1,magpri
          if (fwtmp2(m).le.0.0) go to 2240
          nj=nj+1
          if(nfourier.gt.1)then
          temp=0.
          do i=1,nvesel
          temp=temp+rmp2vs(m,i)*vecta(mk,i)
          enddo
           arsp(nj,nk)=fwtmp2(m)*temp
          else
          arsp(nj,nk)=fwtmp2(m)*rmp2vs(m,mk)
          endif
 2240   continue
        do 2260 m=1,nstark
          if (fwtgam(m).le.0.0) go to 2260
          nj=nj+1
          if(nfourier.gt.1)then
          temp=0.
          do i=1,nvesel
          temp=temp+rgamvs(m,i)*vecta(mk,i)
          enddo
           arsp(nj,nk)=fwtgam(m)*temp
          else
          arsp(nj,nk)=fwtgam(m)*rgamvs(m,mk)
          endif
 2260   continue
        do 92260 m=1,nmsels
          if (fwtbmselt(jtime,m).le.0.0) go to 92260
          nj=nj+1
          arsp(nj,nk)=0.0
92260   continue
        do 92270 m=1,nmsels
          if (fwtemselt(jtime,m).le.0.0) go to 92270
          nj=nj+1
          arsp(nj,nk)=0.0
92270   continue
        do 2270 m=1,nece
          if (fwtece(m).le.0.0) go to 2270  ! add by qian for ece
          nj=nj+1
          if(nfourier.gt.1)then
          temp=0.
          do i=1,nvesel
          temp=temp+recevs(m,i)*vecta(mk,i)
          enddo
           arsp(nj,nk)=fwtece(m)*temp
          else
          arsp(nj,nk)=fwtece(m)*recevs(m,mk)
          endif
 2270   continue
          if (fwtecebz.le.0.0) go to 2275
          nj=nj+1
          arsp(nj,nk)=fwtecebz*recebzvs(mk)
 2275   continue
        if (fwtcur.le.0.0) go to 2280
        nj=nj+1
          if(nfourier.gt.1)then
          temp=0.
          do i=1,nvesel
          temp=temp+vecta(mk,i)
          enddo
           arsp(nj,nk)=fwtcur*temp
          else
          arsp(nj,nk)=fwtcur
          endif
 2280   continue
        if (fwtqa.le.0.0) go to 2285
        nj=nj+1
        arsp(nj,nk)=0.0
 2285 continue
      if (fwtbp.le.0.0) go to 2290
      do 2287 m=2,kffcur
        nj=nj+1
        arsp(nj,nk)=0.0
 2287 continue
 2290 continue
      if (fwtdlc.le.0.0) go to 2292
      nj=nj+1
      arsp(nj,nk)=0.0
 2292 continue
 2296 continue
      do 2297 m=1,nfcoil
        if (fwtfc(m).le.0.0) go to 2297
        nj=nj+1
        arsp(nj,nk)=0.
 2297 continue
      if (kprfit.le.0.or.kdofit.eq.0) go to 2299
      if (npress.le.0) goto 2299
      do 2298 m=1,npress
        if (fwtpre(m).le.0.0) goto 2298
        nj=nj+1
        arsp(nj,nk)=0.0
 2298 continue
 2299 continue
      if (kpressb.eq.2) then
        nj=nj+1
        arsp(nj,nk)=0.0
      endif
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          if (fwtprw(m).gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
          endif
        enddo
      endif
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
        arsp(nj,nk)=0.0
       enddo
      endif
      if (kcalpa.gt.0) then
        do 30190 j=1,kcalpa
          nj=nj+1
          arsp(nj,nk)=0.0
30190   continue
      endif
      if (kcgama.gt.0) then
        do 30200 j=1,kcgama
          nj=nj+1
          arsp(nj,nk)=0.0
30200   continue
      endif
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
!------------------------------------------------------------------------------
!-- Boundary constraints                                                     --
!------------------------------------------------------------------------------
      if (nbdry.gt.0) then
        do j=1,nbdry
          if (fwtbdry(j).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=fwtbry(j)*rbdrvs(j,mk)
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!--  E coil currents                                                         --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
       do m=1,nesum
        if (fwtec(m).gt.0.0) then
        nj=nj+1
        arsp(nj,nk)=0.
        endif
       enddo
      endif
!------------------------------------------------------------------------------
!--  Reference flux                                                          --
!------------------------------------------------------------------------------
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
      endif
!------------------------------------------------------------------------------
!--  add by qian  sum of vessel sement currents=IPV1-IP1                    --
!------------------------------------------------------------------------------
      if (isumvesel.gt.0) then
          nj=nj+1
          if(nfourier.gt.0)then
           temp=0.
            do i=1,nvesel
             temp=temp+vecta(mk,i)
            enddo
             arsp(nj,nk)=fwtcur*temp/100.
           else
             arsp(nj,nk)=1.0
          endif
      endif
 2300 continue
 2310 continue
!-----------------------------------------------------------------------
!-- boundary pressure term for kinetic fitting P(1)                   --
!-- P(1) is an additional fitting parameter in kinetic fitting        --
!-----------------------------------------------------------------------
        if (kprfit.le.0.or.kdofit.eq.0) go to 9210
        if (npress.le.0) goto 9210
        need=need+1
        nk=need
        nj=0
        do 9120 m=1,nsilop
          if (fwtsi(m).le.0.0) go to 9120
          nj=nj+1
          arsp(nj,nk)=0.0
 9120   continue
        do 9140 m=1,magpri
          if (fwtmp2(m).le.0.0) go to 9140
          nj=nj+1
          arsp(nj,nk)=0.0
 9140   continue
        do 9160 m=1,nstark
          if (fwtgam(m).le.0.0) go to 9160
          nj=nj+1
          arsp(nj,nk)=0.0
 9160   continue
        do 99160 m=1,nmsels
          if (fwtbmselt(jtime,m).le.0.0) go to 99160
          nj=nj+1
          arsp(nj,nk)=0.0
99160   continue
        do 99170 m=1,nmsels
          if (fwtemselt(jtime,m).le.0.0) go to 99170
          nj=nj+1
          arsp(nj,nk)=0.0
99170   continue
        do 9170 m=1,nece
          if (fwtece(m).le.0.0) go to 9170
          nj=nj+1
          arsp(nj,nk)=0.0
 9170   continue
          if (fwtecebz.le.0.0) go to 9175
          nj=nj+1
          arsp(nj,nk)=0.0
 9175   continue
        if (fwtcur.le.0.0) go to 9180
        nj=nj+1
        arsp(nj,nk)=0.0
 9180   continue
        if (fwtqa.le.0.0) go to 9190
        nj=nj+1
        arsp(nj,nk)=0.0
 9190 continue
      if (fwtbp.le.0.0) go to 9194
      do 9192 m=2,kffcur
        nj=nj+1
        arsp(nj,nk)=0.0
 9192 continue
 9194 continue
      if (fwtdlc.le.0.0) go to 9196
      nj=nj+1
      arsp(nj,nk)=0.0
 9196 continue
 9202 continue
      do 9203 m=1,nfcoil
        if (fwtfc(m).le.0.0) go to 9203
        nj=nj+1
        arsp(nj,nk)=0.
 9203 continue
      do 9204 m=1,npress
        if (fwtpre(m).le.0.0) goto 9204
        nj=nj+1
        arsp(nj,nk)=1./sigpre(m)*fwtpre(m)
 9204 continue
!-----------------------------------------------------------------------
!--  boundary pressure constraint on P'(1), no coupling to P(1)       --
!-----------------------------------------------------------------------
      if (kpressb.eq.2) then
        nj=nj+1
        arsp(nj,nk)=0.0
      endif
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          if (fwtprw(m).gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
          endif
        enddo
      endif
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
        arsp(nj,nk)=0.0
       enddo
      endif
      if (kcalpa.gt.0) then
        do 30250 j=1,kcalpa
          nj=nj+1
          arsp(nj,nk)=0.0
30250   continue
      endif
      if (kcgama.gt.0) then
        do 30260 j=1,kcgama
          nj=nj+1
          arsp(nj,nk)=0.0
30260   continue
      endif
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
      if (nbdry.gt.0) then
        do j=1,nbdry
          if (fwtbdry(j).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.0
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!--  E coil currents                                                         --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
       do m=1,nesum
        if (fwtec(m).gt.0.0) then
        nj=nj+1
        arsp(nj,nk)=0.
        endif
       enddo
      endif
!------------------------------------------------------------------------------
!--  Reference flux                                                          --
!------------------------------------------------------------------------------
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
      endif
!-----------------------------------------------------------------------------
!-- Summation of Fcoils current
!-----------------------------------------------------------------------------
      if (fitfcsum) then
         nj = nj + 1
         arsp(nj,nk) = 0.0
      endif
 9210 continue
!-----------------------------------------------------------------------
!-- DELZ rigid vertical shift   96/01                                 --
!-----------------------------------------------------------------------
      if (fitdelz.and.nniter.ge.ndelzon) then
        need=need+1
        nsavdz=need
        nk=need
        nj=0
        do 39120 m=1,nsilop
          if (fwtsi(m).le.0.0) go to 39120
          nj=nj+1
          arsp(nj,nk)=fwtsi(m)*gsildz(m)*scadelz
39120   continue
        do 39140 m=1,magpri
          if (fwtmp2(m).le.0.0) go to 39140
          nj=nj+1
          arsp(nj,nk)=fwtmp2(m)*gmp2dz(m)*scadelz
39140   continue
        do 39160 m=1,nstark
          if (fwtgam(m).le.0.0) go to 39160
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgamdz(m)*scadelz
39160   continue
        do 99261 m=1,nmsels
          if (fwtbmselt(jtime,m).le.0.0) go to 99261
          nj=nj+1
          arsp(nj,nk)=fwtbmselt(jtime,m)*rmlsdz(m)*scadelz
99261   continue
        do 99271 m=1,nmsels
          if (fwtemselt(jtime,m).le.0.0) go to 99271
          nj=nj+1
          arsp(nj,nk)=fwtemselt(jtime,m)*relsdz(m)*scadelz
99271   continue
        do 39170 m=1,nece
          if (fwtece(m).le.0.0) go to 39170
          nj=nj+1
          arsp(nj,nk)=fwtece(m)*recedz(m)*scaledz
39170   continue
          if (fwtecebz.le.0.0) go to 39175
          nj=nj+1
          arsp(nj,nk)=fwtecebz*recebzdz*scaledz
39175   continue
        if (fwtcur.le.0.0) go to 39180
          nj=nj+1
          arsp(nj,nk)=fwtcur*fgowdz*scadelz
39180   continue
        if (fwtqa.le.0.0) go to 39190
          nj=nj+1
          arsp(nj,nk)=0.0
39190   continue
        if (fwtbp.le.0.0) go to 39194
          do 39192 m=2,kffcur
            nj=nj+1
            arsp(nj,nk)=0.0
39192     continue
39194   continue
        if (fwtdlc.le.0.0) go to 39196
          nj=nj+1
          arsp(nj,nk)=0.0
39196   continue
        do 39203 m=1,nfcoil
          if (fwtfc(m).le.0.0) go to 39203
          nj=nj+1
          arsp(nj,nk)=0.
39203   continue
        if (kprfit.le.0.or.kdofit.eq.0) go to 39206
        if (npress.le.0) goto 39206
        do 39204 m=1,npress
          if (fwtpre(m).le.0.0) goto 39204
          nj=nj+1
          arsp(nj,nk)=rpredz(m)/sigpre(m)*fwtpre(m)*scadelz
39204   continue
39206   continue
        if (kpressb.eq.2) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        if (kprfit.ge.3.and.npresw.gt.0) then
          do m=1,npresw
            if (fwtprw(m).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=rprwdz(m)/sigprw(m)*fwtprw(m)*scadelz
            endif
          enddo
        endif
        if (kzeroj.gt.0) then
          do i=1,kzeroj
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (kcalpa.gt.0) then
          do 39250 j=1,kcalpa
            nj=nj+1
            arsp(nj,nk)=0.0
39250     continue
        endif
        if (kcgama.gt.0) then
          do 39260 j=1,kcgama
            nj=nj+1
            arsp(nj,nk)=0.0
39260     continue
        endif
        if (kcomega.gt.0) then
          do j=1,kcomega
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (nbdry.gt.0) then
          do j=1,nbdry
            if (fwtbdry(j).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=fwtbry(j)*gbdrdz(j)*scadelz
            endif
          enddo
        endif
        if (iecurr.eq.2) then
         do m=1,nesum
          if (fwtec(m).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.
          endif
         enddo
        endif
        if (fitsiref) then
          if (fwtref.gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.0
          endif
        endif
!-----------------------------------------------------------------------------
!-- Summation of Fcoils current
!-----------------------------------------------------------------------------
        if (fitfcsum) then
           nj = nj + 1
           arsp(nj,nk) = 0.0
        endif
      endif
!-----------------------------------------------------------------------
!--  set up response matrix for advanced divertor coil                --
!--  lacking boundary constraints                                     --
!-----------------------------------------------------------------------
      if (iacoil.le.0) go to 9410
      nkb=need+1
      need=need+nacoil
      do 9400 nk=nkb,need
        mk=nk-nkb+1
        nj=0
        do 9320 m=1,nsilop
          if (fwtsi(m).le.0.0) go to 9320
          nj=nj+1
          arsp(nj,nk)=fwtsi(m)*rsilac(m,mk)
 9320   continue
        do 9340 m=1,magpri
          if (fwtmp2(m).le.0.0) go to 9340
          nj=nj+1
          arsp(nj,nk)=fwtmp2(m)*rmp2ac(m,mk)
 9340   continue
        do 9360 m=1,nstark
          if (fwtgam(m).le.0.0) go to 9360
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgamac(m,mk)
 9360   continue
        do 9370 m=1,nece
          if (fwtece(m).le.0.0) go to 9370
          nj=nj+1
          arsp(nj,nk)=fwtece(m)*receac(m,mk)
 9370   continue
          if (fwtecebz.le.0.0) go to 9375
          nj=nj+1
          arsp(nj,nk)=fwtecebz*recebzac(mk)
 9375   continue
        if (fwtcur.le.0.0) go to 9380
        nj=nj+1
        arsp(nj,nk)=0.0
 9380   continue
        if (fwtqa.le.0.0) go to 9385
        nj=nj+1
        arsp(nj,nk)=0.0
 9385 continue
      if (fwtbp.le.0.0) go to 9390
      do 9387 m=2,kffcur
        nj=nj+1
        arsp(nj,nk)=0.0
 9387 continue
 9390 continue
      if (fwtdlc.le.0.0) go to 9392
      nj=nj+1
      arsp(nj,nk)=0.0
 9392 continue
      do 9397 m=1,nfcoil
        if (fwtfc(m).le.0.0) go to 9397
        nj=nj+1
        arsp(nj,nk)=0.
 9397 continue
      if (kprfit.le.0.or.kdofit.eq.0) go to 9399
      if (npress.le.0) goto 9399
      do 9398 m=1,npress
        if (fwtpre(m).le.0.0) goto 9398
        nj=nj+1
        arsp(nj,nk)=0.0
 9398 continue
 9399 continue
      if (kpressb.eq.2) then
        nj=nj+1
        arsp(nj,nk)=0.0
      endif
!---------------------------------------------------------------------------
!-- rotational pressure                                                   --
!---------------------------------------------------------------------------
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          if (fwtprw(m).gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
          endif
        enddo
      endif
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
        arsp(nj,nk)=0.0
       enddo
      endif
      if (kcalpa.gt.0) then
        do 39390 j=1,kcalpa
          nj=nj+1
          arsp(nj,nk)=0.0
39390   continue
      endif
      if (kcgama.gt.0) then
        do 39400 j=1,kcgama
          nj=nj+1
          arsp(nj,nk)=0.0
39400   continue
      endif
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
!------------------------------------------------------------------------------
!--  E coil currents                                                         --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
       do m=1,nesum
        if (fwtec(m).gt.0.0) then
        nj=nj+1
        arsp(nj,nk)=0.
        endif
       enddo
      endif
!------------------------------------------------------------------------------
!--  Reference flux                                                          --
!------------------------------------------------------------------------------
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
      endif
!-----------------------------------------------------------------------------
!-- Summation of Fcoils current
!-----------------------------------------------------------------------------
      if (fitfcsum) then
         nj = nj + 1
         arsp(nj,nk) = 0.0
      endif
 9400 continue
 9410 continue
!-----------------------------------------------------------------------
!--  set up response matrix for E coils                               --
!-----------------------------------------------------------------------
      if (iecurr.ne.2) go to 48410
      nkb=need+1
      need=need+nesum
      do 48400 nk=nkb,need
        mk=nk-nkb+1
        nj=0
        do 48320 m=1,nsilop
          if (fwtsi(m).le.0.0) go to 48320
          nj=nj+1
          arsp(nj,nk)=fwtsi(m)*rsilec(m,mk)
48320   continue
        do 48340 m=1,magpri
          if (fwtmp2(m).le.0.0) go to 48340
          nj=nj+1
          arsp(nj,nk)=fwtmp2(m)*rmp2ec(m,mk)
48340   continue
        do 48360 m=1,nstark
          if (fwtgam(m).le.0.0) go to 48360
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgamec(m,mk)
48360   continue
        do 98360 m=1,nmsels
          if (fwtbmselt(jtime,m).le.0.0) go to 98360
          nj=nj+1
          arsp(nj,nk)=0.0
98360   continue
        do 98370 m=1,nmsels
          if (fwtemselt(jtime,m).le.0.0) go to 98370
          nj=nj+1
          arsp(nj,nk)=0.0
98370   continue
        do 48370 m=1,nece
          if (fwtece(m).le.0.0) go to 48370
          nj=nj+1
          arsp(nj,nk)=fwtece(m)*receec(m,mk)
48370   continue
          if (fwtecebz.le.0.0) go to 48375
          nj=nj+1
          arsp(nj,nk)=fwtecebz*recebzec(mk)
48375   continue
        if (fwtcur.le.0.0) go to 48380
        nj=nj+1
        arsp(nj,nk)=0.0
48380   continue
        if (fwtqa.le.0.0) go to 48385
        nj=nj+1
        arsp(nj,nk)=0.0
48385 continue
      if (fwtbp.le.0.0) go to 48390
      do 48387 m=2,kffcur
        nj=nj+1
        arsp(nj,nk)=0.0
48387 continue
48390 continue
      if (fwtdlc.le.0.0) go to 48392
      nj=nj+1
      arsp(nj,nk)=0.0
48392 continue
      do 48397 m=1,nfcoil
        if (fwtfc(m).le.0.0) go to 48397
        nj=nj+1
        arsp(nj,nk)=0.
48397 continue
      if (kprfit.le.0.or.kdofit.eq.0) go to 48399
      if (npress.le.0) goto 48399
      do 48398 m=1,npress
        if (fwtpre(m).le.0.0) goto 48398
        nj=nj+1
        arsp(nj,nk)=0.0
48398 continue
48399 continue
      if (kpressb.eq.2) then
        nj=nj+1
        arsp(nj,nk)=0.0
      endif
!---------------------------------------------------------------------------
!-- rotational pressure                                                   --
!---------------------------------------------------------------------------
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          if (fwtprw(m).gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
          endif
        enddo
      endif
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
        arsp(nj,nk)=0.0
       enddo
      endif
      if (kcalpa.gt.0) then
        do 48490 j=1,kcalpa
          nj=nj+1
          arsp(nj,nk)=0.0
48490   continue
      endif
      if (kcgama.gt.0) then
        do 48495 j=1,kcgama
          nj=nj+1
          arsp(nj,nk)=0.0
48495   continue
      endif
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
!------------------------------------------------------------------------------
!-- Boundary constraints                                                     --
!------------------------------------------------------------------------------
      if (nbdry.gt.0) then
        do j=1,nbdry
          if (fwtbdry(j).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=fwtbry(j)*rbdrec(j,mk)
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!--  E coil currents                                                         --
!------------------------------------------------------------------------------
      do 48497 m=1,nesum
        if (fwtec(m).le.0.0) go to 48497
        nj=nj+1
        arsp(nj,nk)=0.
        if (m.eq.mk) arsp(nj,nk)=fwtec(m)
48497 continue
!------------------------------------------------------------------------------
!--  Reference flux                                                          --
!------------------------------------------------------------------------------
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
      endif
!-----------------------------------------------------------------------------
!-- Summation of Fcoils current
!-----------------------------------------------------------------------------
      if (fitfcsum) then
         nj = nj + 1
         arsp(nj,nk) = 0.0
      endif
48400 continue
48410 continue
!------------------------------------------------------------------------------
!--  fitting relative flux, set up response for fitted reference flux        --
!------------------------------------------------------------------------------
        if (.not.fitsiref) goto 52501
        need=need+1
        nj=0
        nk=need
        do 52020 m=1,nsilop
          if (fwtsi(m).le.0.0) go to 52020
          nj=nj+1
          arsp(nj,nk)=-fwtsi(m)*scalesir
52020   continue
        do 52040 m=1,magpri
          if (fwtmp2(m).le.0.0) go to 52040
          nj=nj+1
          arsp(nj,nk)=0.0
52040   continue
        do 52060 m=1,nstark
          if (fwtgam(m).le.0.0) go to 52060
          nj=nj+1
          arsp(nj,nk)=0.0
52060   continue
        do 92161 m=1,nmsels
          if (fwtbmselt(jtime,m).le.0.0) go to 92161
          nj=nj+1
          arsp(nj,nk)=0.0
92161   continue
        do 92171 m=1,nmsels
          if (fwtemselt(jtime,m).le.0.0) go to 92171
          nj=nj+1
          arsp(nj,nk)=0.0
92171   continue
        do 52070 m=1,nece
          if (fwtece(m).le.0.0) go to 52070
          nj=nj+1
          arsp(nj,nk)=0.0
52070   continue
          if (fwtecebz.le.0.0) go to 52075
          nj=nj+1
          arsp(nj,nk)=0.0
52075   continue
        if (fwtcur.le.0.0) go to 52080
        nj=nj+1
        arsp(nj,nk)=0.0
52080   continue
        if (fwtqa.le.0.0) go to 52085
        nj=nj+1
        arsp(nj,nk)=0.0
52085 continue
      if (fwtbp.le.0.0) go to 52090
      do 52087 m=2,kffcur
        nj=nj+1
        arsp(nj,nk)=0.0
52087 continue
52090 continue
      if (fwtdlc.le.0.0) go to 52092
      nj=nj+1
      arsp(nj,nk)=0.0
52092 continue
52096 continue
      do 52097 m=1,nfcoil
        if (fwtfc(m).le.0.0) go to 52097
        nj=nj+1
        arsp(nj,nk)=0.
52097 continue
!--------------------------------------------------------------------------
!--  pressure                                                            --
!--------------------------------------------------------------------------
      if (kprfit.le.0.or.kdofit.eq.0) go to 52099
      if (npress.le.0) goto 52099
      do 52098 m=1,npress
        if (fwtpre(m).le.0.0) goto 52098
        nj=nj+1
        arsp(nj,nk)=0.0
52098 continue
52099 continue
      if (kpressb.eq.2) then
        nj=nj+1
        arsp(nj,nk)=0.0
      endif
!---------------------------------------------------------------------------
!-- rotational pressure                                                   --
!---------------------------------------------------------------------------
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          if (fwtprw(m).gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
          endif
        enddo
      endif
!---------------------------------------------------------------------------
!--  J(PSIWANT)                                                           --
!---------------------------------------------------------------------------
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
        arsp(nj,nk)=0.0
       enddo
      endif
!----------------------------------------------------------------------------
!-- P', FF', and rotational                                                --
!----------------------------------------------------------------------------
      if (kcalpa.gt.0) then
        do 53050 j=1,kcalpa
          nj=nj+1
          arsp(nj,nk)=0.0
53050   continue
      endif
      if (kcgama.gt.0) then
        do 53060 j=1,kcgama
          nj=nj+1
          arsp(nj,nk)=0.0
53060   continue
      endif
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
!------------------------------------------------------------------------------
!-- Boundary constraints                                                     --
!------------------------------------------------------------------------------
      if (nbdry.gt.0) then
        do j=1,nbdry
          if (fwtbdry(j).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.0
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!--  E coil currents                                                         --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
       do m=1,nesum
        if (fwtec(m).gt.0.0) then
        nj=nj+1
        arsp(nj,nk)=0.
        endif
       enddo
      endif
!------------------------------------------------------------------------------
!--  Reference flux                                                          --
!------------------------------------------------------------------------------
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=fwtref*scalesir
        endif
      endif
!-----------------------------------------------------------------------------
!-- Summation of Fcoils current                                             --
!-----------------------------------------------------------------------------
      if (fitfcsum) then
         nj = nj + 1
         arsp(nj,nk) = 0.0
      endif
52501 continue
!------------------------------------------------------------------------------
!-- set up response matrix for ER, fitting ER                                --
!------------------------------------------------------------------------------
      if (keecur.le.0.or.kdomse.gt.0) goto 72111
        needs=need
        need=need+keecur
      do 72100 nk=needs+1,need
        nkk=nk-needs
        nj=0
        do 72020 m=1,nsilop
          if (fwtsi(m).le.0.0) go to 72020
          nj=nj+1
          arsp(nj,nk)=0.0
72020   continue
        do 72040 m=1,magpri
          if (fwtmp2(m).le.0.0) go to 72040
          nj=nj+1
          arsp(nj,nk)=0.0
72040   continue
        do 72060 m=1,nstark
          if (fwtgam(m).le.0.0) go to 72060
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgamer(m,nkk)
72060   continue
        do 92760 m=1,nmsels
          if (fwtbmselt(jtime,m).le.0.0) go to 92760
          nj=nj+1
          arsp(nj,nk)=0.0
92760   continue
        do 92770 m=1,nmsels
          if (fwtemselt(jtime,m).le.0.0) go to 92770
          nj=nj+1
          arsp(nj,nk)=fwtemselt(jtime,m)*relser(m,nkk)
92770   continue
        do 72070 m=1,nece
          if (fwtece(m).le.0.0) go to 72070
          nj=nj+1
          arsp(nj,nk)=0.0
72070   continue
          if (fwtecebz.le.0.0) go to 72075
          nj=nj+1
          arsp(nj,nk)=0.0
72075   continue
        if (fwtcur.le.0.0) go to 72080
        nj=nj+1
        arsp(nj,nk)=0.0
72080   continue
        if (fwtqa.le.0.0) go to 72085
        nj=nj+1
        arsp(nj,nk)=0.0
72085 continue
      if (fwtbp.le.0.0) go to 72090
      do 72087 m=2,kffcur
        nj=nj+1
        arsp(nj,nk)=0.0
72087 continue
72090 continue
      if (fwtdlc.le.0.0) go to 72092
      nj=nj+1
      arsp(nj,nk)=0.0
72092 continue
      do 72097 m=1,nfcoil
        if (fwtfc(m).le.0.0) go to 72097
        nj=nj+1
        arsp(nj,nk)=0.
72097 continue
      if (kprfit.le.0.or.kdofit.eq.0) go to 72099
      if (npress.le.0) goto 72099
      do 72098 m=1,npress
        if (fwtpre(m).le.0.0) goto 72098
        nj=nj+1
        arsp(nj,nk)=0.0
72098 continue
72099 continue
      if (kpressb.eq.2) then
        nj=nj+1
        arsp(nj,nk)=0.0
      endif
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          if (fwtprw(m).gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
          endif
        enddo
      endif
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
        arsp(nj,nk)=0.0
       enddo
      endif
      if (kcalpa.gt.0) then
        do j=1,kcalpa
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
      if (kcgama.gt.0) then
        do j=1,kcgama
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
      if (nbdry.gt.0) then
        do j=1,nbdry
          if (fwtbdry(j).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.0
          endif
        enddo
      endif
      if (iecurr.eq.2) then
       do m=1,nesum
        if (fwtec(m).gt.0.0) then
        nj=nj+1
        arsp(nj,nk)=0.
        endif
       enddo
      endif
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
      endif
!-----------------------------------------------------------------------------
!-- Summation of Fcoils current
!-----------------------------------------------------------------------------
      if (fitfcsum) then
         nj = nj + 1
         arsp(nj,nk) = 0.0
      endif
72100 continue
72111 continue
!-----------------------------------------------------------------------
!-- Response for Pedge hyperbolic tangent                             --
!-----------------------------------------------------------------------
      if (kedgep.eq.0) go to 73311
        need=need+1
        nk=need
        npedge=need
        nj=0
        do 73120 m=1,nsilop
          if (fwtsi(m).le.0.0) go to 73120
          nj=nj+1
          arsp(nj,nk)=fwtsi(m)*rsilpe(m)
73120   continue
        do 73140 m=1,magpri
          if (fwtmp2(m).le.0.0) go to 73140
          nj=nj+1
          arsp(nj,nk)=fwtmp2(m)*rmp2pe(m)
73140   continue
        do 73160 m=1,nstark
          if (fwtgam(m).le.0.0) go to 73160
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgampe(m)
73160   continue
        do 93160 m=1,nmsels
          if (fwtbmselt(jtime,m).le.0.0) go to 93160
          nj=nj+1
          arsp(nj,nk)=0.0
93160   continue
        do 93170 m=1,nmsels
          if (fwtemselt(jtime,m).le.0.0) go to 93170
          nj=nj+1
          arsp(nj,nk)=0.0
93170   continue
        if (fwtcur.le.0.0) go to 73180
        nj=nj+1
        arsp(nj,nk)=fwtcur*fgowpe
73180   continue
        if (fwtqa.le.0.0) go to 73190
        nj=nj+1
        arsp(nj,nk)=0.0
73190   continue
        if (fwtbp.le.0.0) go to 73194
        do 73192 m=2,kffcur
          nj=nj+1
          arsp(nj,nk)=0.0
73192   continue
73194   continue
        if (fwtdlc.le.0.0) go to 73196
        nj=nj+1
        arsp(nj,nk)=0.0
73196   continue
        do 73203 m=1,nfcoil
          if (fwtfc(m).le.0.0) go to 73203
          nj=nj+1
          arsp(nj,nk)=0.
73203   continue
!--------------------------------------------------------------------
!--  pressure                                                      --
!--------------------------------------------------------------------
        if (kprfit.le.0.or.kdofit.eq.0) go to 73206
        if (npress.le.0) goto 73206
        do 73204 m=1,npress
          if (fwtpre(m).le.0.0) goto 73204
          nj=nj+1
          arsp(nj,nk)=rprepe(m)/sigpre(m)*fwtpre(m)
73204 continue
73206 continue
!--------------------------------------------------------------------
!-- P'(1)                                                          --
!--------------------------------------------------------------------
      if (kpressb.eq.2) then
        nj=nj+1
        arsp(nj,nk)=1./sigppb/darea/cosh(s1edge)**2/pe_width/sidif
      endif
!---------------------------------------------------------------------
!-- rotational pressure                                             --
!---------------------------------------------------------------------
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          if (fwtprw(m).gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
          endif
        enddo
      endif
!----------------------------------------------------------------------
!--  J(PSIWANT) constraint                                           --
!----------------------------------------------------------------------
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
        siedge=(sizeroj(i)-pe_psin)/pe_width
        arsp(nj,nk)=rxxx(i)*fwtxxzj/cosh(siedge)**2/pe_width/sidif
       enddo
      endif
!-------------------------------------------------------------------------
!--  p' and ff' constraints                                             --
!-------------------------------------------------------------------------
      if (kcalpa.gt.0) then
        do 73240 j=1,kcalpa
          nj=nj+1
          arsp(nj,nk)=calpa(kppcur+1,j)*fwtxxa
73240   continue
      endif
      if (kcgama.gt.0) then
        do 73250 j=1,kcgama
          nj=nj+1
          arsp(nj,nk)=0.0
73250   continue
      endif
!---------------------------------------------------------------------
!-- rotational constraints                                          --
!---------------------------------------------------------------------
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
!------------------------------------------------------------------------------
!-- Boundary constraints due to plasma contributions                         --
!------------------------------------------------------------------------------
      if (nbdry.gt.0) then
        do j=1,nbdry
          if (fwtbdry(j).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=fwtbry(j)*gbdrpe(j)
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!--  E coil currents                                                         --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
       do m=1,nesum
        if (fwtec(m).gt.0.0) then
        nj=nj+1
        arsp(nj,nk)=0.
        endif
       enddo
      endif
!------------------------------------------------------------------------------
!--  Reference flux                                                          --
!------------------------------------------------------------------------------
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
      endif
!-----------------------------------------------------------------------------
!-- Summation of Fcoils current
!-----------------------------------------------------------------------------
      if (fitfcsum) then
         nj = nj + 1
         arsp(nj,nk) = 0.0
      endif
73311 continue
!-----------------------------------------------------------------------
!-- Response for f2edge hyperbolic tangent                            --
!-----------------------------------------------------------------------
      if (kedgef.eq.0) go to 74311
        need=need+1
        nk=need
        nfedge=need
        nj=0
        do 74120 m=1,nsilop
          if (fwtsi(m).le.0.0) go to 74120
          nj=nj+1
          arsp(nj,nk)=fwtsi(m)*rsilfe(m)
74120   continue
        do 74140 m=1,magpri
          if (fwtmp2(m).le.0.0) go to 74140
          nj=nj+1
          arsp(nj,nk)=fwtmp2(m)*rmp2fe(m)
74140   continue
        do 74160 m=1,nstark
          if (fwtgam(m).le.0.0) go to 74160
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgamfe(m)
74160   continue
        do 94160 m=1,nmsels
          if (fwtbmselt(jtime,m).le.0.0) go to 94160
          nj=nj+1
          arsp(nj,nk)=fwtbmselt(jtime,m)*rmlsfe(m)
94160   continue
        do 94170 m=1,nmsels
          if (fwtemselt(jtime,m).le.0.0) go to 94170
          nj=nj+1
          arsp(nj,nk)=0.0
94170   continue
        if (fwtcur.le.0.0) go to 74180
        nj=nj+1
        arsp(nj,nk)=fwtcur*fgowfe
74180   continue
        if (fwtqa.le.0.0) go to 74190
        nj=nj+1
        arsp(nj,nk)=0.0
74190   continue
        if (fwtbp.le.0.0) go to 74194
        do 74192 m=2,kffcur
          nj=nj+1
          arsp(nj,nk)=0.0
74192   continue
74194   continue
        if (fwtdlc.le.0.0) go to 74196
        nj=nj+1
        arsp(nj,nk)=0.0
74196   continue
        do 74203 m=1,nfcoil
          if (fwtfc(m).le.0.0) go to 74203
          nj=nj+1
          arsp(nj,nk)=0.
74203   continue
!--------------------------------------------------------------------
!--  pressure                                                      --
!--------------------------------------------------------------------
        if (kprfit.le.0.or.kdofit.eq.0) go to 74206
        if (npress.le.0) goto 74206
        do 74204 m=1,npress
          if (fwtpre(m).le.0.0) goto 74204
          nj=nj+1
          arsp(nj,nk)=0.0
74204 continue
74206 continue
!--------------------------------------------------------------------
!-- P'(1)                                                          --
!--------------------------------------------------------------------
      if (kpressb.eq.2) then
        nj=nj+1
        arsp(nj,nk)=0.0
      endif
!---------------------------------------------------------------------
!-- rotational pressure                                             --
!---------------------------------------------------------------------
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          if (fwtprw(m).gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
          endif
        enddo
      endif
!----------------------------------------------------------------------
!--  J(PSIWANT) constraint                                           --
!----------------------------------------------------------------------
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
        siedge=(sizeroj(i)-fe_psin)/fe_width
        arsp(nj,nk)=fwtxxzj/rxxxf(i)/cosh(siedge)**2/fe_width/sidif
       enddo
      endif
!-------------------------------------------------------------------------
!--  p' and ff' constraints                                             --
!-------------------------------------------------------------------------
      if (kcalpa.gt.0) then
        do 74240 j=1,kcalpa
          nj=nj+1
          arsp(nj,nk)=0.0
74240   continue
      endif
      if (kcgama.gt.0) then
        do 74250 j=1,kcgama
          nj=nj+1
          arsp(nj,nk)=cgama(kffcur+1,j)*fwtxxf
74250   continue
      endif
!---------------------------------------------------------------------
!-- rotational constraints                                          --
!---------------------------------------------------------------------
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
!------------------------------------------------------------------------------
!-- Boundary constraints due to plasma contributions                         --
!------------------------------------------------------------------------------
      if (nbdry.gt.0) then
        do j=1,nbdry
          if (fwtbdry(j).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=fwtbry(j)*gbdrfe(j)
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!--  E coil currents                                                         --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
       do m=1,nesum
        if (fwtec(m).gt.0.0) then
        nj=nj+1
        arsp(nj,nk)=0.
        endif
       enddo
      endif
!------------------------------------------------------------------------------
!--  Reference flux                                                          --
!------------------------------------------------------------------------------
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
      endif
!-----------------------------------------------------------------------------
!-- Summation of Fcoils current
!-----------------------------------------------------------------------------
      if (fitfcsum) then
         nj = nj + 1
         arsp(nj,nk) = 0.0
      endif
74311 continue
!-----------------------------------------------------------------------
!-- stabilization constraint for dz, setup here for all fitting paras --
!-----------------------------------------------------------------------
      if (fitdelz.and.nniter.ge.ndelzon) then
        if (stabdz.gt.0.0) then
        nj=nj+1
        do i=1,need
          arsp(nj,i )=0.0
        enddo
        arsp(nj,nsavdz)=stabdz
        endif
      endif
!-----------------------------------------------------------------------
!-- user defined constraint equations
!-----------------------------------------------------------------------
      if (kccoils.gt.0) then
         do m=1,kccoils
           nj=nj+1
           do i=1,nfcoil
             arsp(nj,i)=ccoils(i,m)
           enddo
         enddo
      endif
!-----------------------------------------------------------------------
!-- set up the corresponding data, RHS                                --
!-----------------------------------------------------------------------
      nj=0
      do 2332 m=1,nsilop
        if (fwtsi(m).le.0.0) go to 2332
        nj=nj+1
        brsp(nj)=fwtsi(m)*silopt(jtime,m)
        if (iecurr.ne.1) go to 2320
        ework=0.0
        do 2315 ne=1,nesum
          ework=ework+rsilec(m,ne)*ecurrt(ne)
 2315   continue
        brsp(nj)=brsp(nj)-fwtsi(m)*ework
 2320   if (ivesel.le.0.or.ifitvs.gt.0) go to 2332
        ework=0.0
        do 2330 ne=1,nvesel
          ework=ework+rsilvs(m,ne)*vcurrt(ne)
 2330   continue
        brsp(nj)=brsp(nj)-fwtsi(m)*ework
 2332 continue
      do 2352 m=1,magpri
        if (fwtmp2(m).le.0.0) go to 2352
        nj=nj+1
        brsp(nj)=fwtmp2(m)*expmpi(jtime,m)
        if (iecurr.ne.1) go to 2340
        ework=0.0
        do 2335 ne=1,nesum
          ework=ework+rmp2ec(m,ne)*ecurrt(ne)
 2335   continue
        brsp(nj)=brsp(nj)-fwtmp2(m)*ework
 2340 continue
        if (ivesel.le.0.or.ifitvs.gt.0) go to 2352
        ework=0.0
        do 2345 ne=1,nvesel
          ework=ework+rmp2vs(m,ne)*vcurrt(ne)
 2345   continue
        brsp(nj)=brsp(nj)-fwtmp2(m)*ework
 2352 continue
      do 2360 m=1,nstark
        if (fwtgam(m).le.0.0) go to 2360
        nj=nj+1
        brsp(nj)=fwtgam(m)*rhsgam(jtime,m)
        if (iecurr.ne.1) go to 2356
        ework=0.0
        do 2354 ne=1,nesum
          ework=ework+rgamec(m,ne)*ecurrt(ne)
 2354   continue
        brsp(nj)=brsp(nj)-fwtgam(m)*ework
 2356 continue
        if (ivesel.le.0.or.ifitvs.gt.0) go to 2360
        ework=0.0
        do 2358 ne=1,nvesel
          ework=ework+rgamvs(m,ne)*vcurrt(ne)
 2358   continue
        brsp(nj)=brsp(nj)-fwtgam(m)*ework
 2360 continue
!-----------------------------------------------------------------------
!-- set up the corresponding MSE-LS data + MSE-LE                     --
!-----------------------------------------------------------------------
      do 98461 m=1,nmsels
        if (fwtbmselt(jtime,m).le.0.0) go to 98461
        nj=nj+1
        brsp(nj)=fwtbmselt(jtime,m)*(rhsmls(jtime,m)+bmselt(jtime,m)**2)
98461 continue
      do 98471 m=1,nmsels
        if (fwtemselt(jtime,m).le.0.0) go to 98471
        nj=nj+1
        brsp(nj)=fwtemselt(jtime,m)*emselt(jtime,m)
98471 continue
      do 2370 m=1,nece
        if (fwtece(m).le.0.0) go to 2370
        nj=nj+1
        brsp(nj)=fwtece(m)*brspece(jtime,m)
        if (iecurr.ne.1) go to 2366
        ework=0.0
        do 2364 ne=1,nesum
          ework=ework+receec(m,ne)*ecurrt(ne)
 2364   continue
        brsp(nj)=brsp(nj)-fwtece(m)*ework
 2366 continue
        if (ivesel.le.0.or.ifitvs.gt.0) go to 2370
        ework=0.0
        do 2368 ne=1,nvesel
          ework=ework+recevs(m,ne)*vcurrt(ne)
 2368   continue
        brsp(nj)=brsp(nj)-fwtece(m)*ework
 2370 continue
        if (fwtecebz.le.0.0) go to 2379
        nj=nj+1
        brsp(nj)=fwtecebz*brspecebz(jtime)
        if (iecurr.ne.1) go to 2376
        ework=0.0
        do 2374 ne=1,nesum
          ework=ework+recebzec(ne)*ecurrt(ne)
 2374   continue
        brsp(nj)=brsp(nj)-fwtecebz*ework
 2376 continue
        if (ivesel.le.0.or.ifitvs.gt.0) go to 2379
        ework=0.0
        do 2378 ne=1,nvesel
          ework=ework+recebzvs(ne)*vcurrt(ne)
 2378   continue
        brsp(nj)=brsp(nj)-fwtecebz*ework
 2379 continue
      if (fwtcur.le.0.0) go to 2380
      nj=nj+1
      brsp(nj)=fwtcur*pasmat(jtime)
 2380 continue
      if (fwtqa.le.0.0) go to 2400
      nj=nj+1
      do i=1,kppcur
        arsp(nj,nfcoil+i)=fwtqa*rqajx(i)
      enddo
      do i=1,kffcur
        arsp(nj,nbase+i)=fwtqa*rqafx(i)
      enddo
      if (kedgep.gt.0) then
        arsp(nj,npedge)=fwtqa*rqape
      endif
      if (kedgef.gt.0) then
        arsp(nj,nfedge)=fwtqa*rqafe
      endif
      if (kvtor.gt.0) then
        do i=1,kwwcur
          arsp(nj,nfcoil+kpcurn+i)=fwtqa*rqawx(i)
        enddo
      endif
      brsp(nj)=fwtqa/qvfit*pasmat(jtime)/abs(pasmat(jtime))
 2400 continue
      if (fwtbp.le.0.0) go to 2450
      fwtbpp=fwtbp/abs(brsold(nfcoil+1)*brsold(nbase+2))
      do 2420 i=2,kffcur
        nj=nj+1
        arsp(nj,nfcoil+1)=fwtbpp*brsold(nbase+i)
        arsp(nj,nfcoil+i)=-fwtbpp*brsold(nbase+1)
        arsp(nj,nbase+1)=-fwtbpp*brsold(nfcoil+i)
        arsp(nj,nbase+i)=fwtbpp*brsold(nfcoil+1)
        brsp(nj)=fwtbpp*(brsold(nfcoil+1)*brsold(nbase+i)- &
                        brsold(nfcoil+i)*brsold(nbase+1))
 2420 continue
 2450 continue
      if (fwtdlc.le.0.0) go to 2470
      nj=nj+1
      do 2460 i=1,kffcur
        arsp(nj,nbase+i)=fwtdlc*rspdlc(i)
 2460 continue
      if (kedgef.gt.0) then
        arsp(nj,nfedge)=fwtdlc*rdlcfe
      endif
      brsp(nj)=-fwtdlc*diamag(jtime)
 2470 continue
 2475 continue
      do 2476 i=1,nfcoil
        if (fwtfc(i).le.0.0) go to 2476
        nj=nj+1
        brsp(nj)=fccurt(jtime,i)*fwtfc(i)
 2476 continue
!------------------------------------------------------------------
!--  pressure                                                    --
!------------------------------------------------------------------
      if (kprfit.le.0.or.kdofit.eq.0) go to 2480
      if (npress.le.0) goto 2480
      do 2477 i=1,npress
        if (fwtpre(i).le.0.0) goto 2477
        nj=nj+1
        brsp(nj)=pressr(i)/sigpre(i)*fwtpre(i)
 2477 continue
 2480 continue
      if (kpressb.eq.2) then
        nj=nj+1
        brsp(nj)=prespb/sigppb
      endif
!-------------------------------------------------------------------
!--  rotational pressure                                          --
!-------------------------------------------------------------------
      if (kprfit.ge.3.and.npresw.gt.0) then
        do i=1,npresw
          if (fwtprw(i).gt.0.0) then
          nj=nj+1
          brsp(nj)=(presw(i)-preswb)/sigprw(i)*fwtprw(i)
          endif
        enddo
      endif
!--------------------------------------------------------------------
!-- J(PSIWANT)                                                     --
!--------------------------------------------------------------------
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
        brsp(nj)=fwtxxzj*vzeroj(i)*darea*pasmat(jtime)/carea
       enddo
      endif
!--------------------------------------------------------------------
!--  P' and FF'                                                    --
!--------------------------------------------------------------------
      if (kcalpa.gt.0) then
        do 30430 j=1,kcalpa
          nj=nj+1
          brsp(nj)=fwtxxa*xalpa(j)*darea
30430   continue
      endif
      if (kcgama.gt.0) then
        do 30440 j=1,kcgama
          nj=nj+1
          brsp(nj)=fwtxxf*xgama(j)*darea
30440   continue
      endif
!---------------------------------------------------------------------
!--  rotational                                                     --
!---------------------------------------------------------------------
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          brsp(nj)=xomega(j)*darea*fwtxxo
        enddo
      endif
!------------------------------------------------------------------------------
!-- Boundary constraints                                                     --
!------------------------------------------------------------------------------
      if (nbdry.gt.0) then
        do j=1,nbdry
          if (fwtbdry(j).gt.0.0) then
            nj=nj+1
            brsp(nj)=0.0
            if (iecurr.eq.1) then
              do k=1,nesum
                brsp(nj)=brsp(nj)+rbdrec(j,k)*ecurrt(k)
              enddo
            endif
            brsp(nj)=fwtbry(j)*(psibry-brsp(nj))
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!--  E coil currents                                                         --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
      do i=1,nesum
        if (fwtec(i).gt.0.0) then
         nj=nj+1
         brsp(nj)=ecurrt(i)*fwtec(i)
        endif
      enddo
      endif
!------------------------------------------------------------------------------
!--  fitting relative flux                                                   --
!------------------------------------------------------------------------------
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          brsp(nj)=fwtref*psiref(jtime)
        endif
      endif
!-----------------------------------------------------------------------
!-- Summation of F-coils current
!----------------------------------------------------------------------
      if (fitfcsum) then
         nj = nj + 1
         brsp(nj) = 0.0
      endif
!-----------------------------------------------------------------------
!-- stabilization constraint for dz                                   --
!-----------------------------------------------------------------------
      if (fitdelz.and.nniter.ge.ndelzon) then
        if (stabdz.gt.0.0) then
        nj=nj+1
        brsp(nj)=0.0
        endif
      endif
!-----------------------------------------------------------------------
!-- user defined constraint equations
!-----------------------------------------------------------------------
      if (kccoils.gt.0) then
         do j=1,kccoils
           nj=nj+1
           brsp(nj)=xcoils(j)
         enddo
      endif
!
      nnn=1
      ncrsp = 0
      nfffff = nfcoil
      call ppcnst(ncrsp,crsp,z,nfffff)
      call ffcnst(ncrsp,crsp,z,nfffff)
      call wwcnst(ncrsp,crsp,z,nfffff)
      if (keecur.gt.0.and.kdomse.eq.0) then
      nfffff = needs-kppcur-kffcur
      needer = needs
      call eecnst(ncrsp,crsp,z,nfffff)
      endif
!---------------------------------------------------------------------
!-- preconditioning A matrix if need                                --
!---------------------------------------------------------------------
      if (scalea) then
         call dgeequ(nj,need,arsp,nrsmat,rowscale,colscale, &
                     rowcnd,colcnd,arspmax,infosc)
         do j = 1,nj
         do i = 1,need
            arsp(j,i) = arsp(j,i) * colscale(i)
         enddo
         enddo
      endif
      if (ncrsp .le. 0) then
         call sdecm(arsp,nrsmat,nj,need,brsp,nrsmat, &
                      nnn,wrsp,work,ier)
         if (ier.ne.129) go to 2500
         write (nttyo,8000) ier
! MPI >>>
#if defined(USEMPI)
         call mpi_stop
#else
         stop
#endif
! MPI <<<
 2500    continue
!-----------------------------------------------------------------------
!--  unfold fitting parameters                                        --
!-----------------------------------------------------------------------
         if ( wrsp(need).eq.0 ) then
           goto 2656
         end if
         condno=wrsp(1)/wrsp(need)
         toler=condin*wrsp(1)
         do 2600 i=1,need
           t=0.0
           if (wrsp(i).gt.toler) t=brsp(i)/wrsp(i)
           work(i)=t
 2600    continue
         do 2650 i=1,need
           brsp(i)=0.0
           do 2650 j=1,need
             brsp(i)=brsp(i)+arsp(i,j)*work(j)
 2650    continue
      else
         do 2655 j=1,nrsmat
           b(j) = brsp(j)
         2655 continue

         call dgglse(nj,need,ncrsp,arsp,nrsmat,crsp,4*(npcurn-2)+6+ &
                   npcurn*npcurn,b,z,brsp,work,nrsma2,info,condno)
         if (info.eq.0) goto 2657
         2656   continue
         write (nttyo,8000) info
! MPI >>>
#if defined(USEMPI)
         call mpi_stop
         stop
#else
         stop
#endif
! MPI <<<
 2657    continue
      endif
!----------------------------------------------------------------------
!--  rescale results if A is preconditioned                          --
!----------------------------------------------------------------------
      if (scalea) then
        do i=1, need
          brsp(i)=brsp(i)*colscale(i)
        enddo
      endif
      nload=nfnwcr
      if (ifitvs.gt.0) then
       if(nfourier.gt.1) then
        do j=1,nvesel
          temp=0.
         do 2671 i=1,(nfourier*2+1)
          temp=temp+brsp(i+nfnwcr)*vecta(i,j)
 2671    continue
         vcurrt(j)=temp
        enddo
       nload=nload+(nfourier*2+1)
       else
         do 2670 i=1,nvesel
          vcurrt(i)=brsp(i+nfnwcr)
 2670    continue
       nload=nload+nvesel
       endif
      endif
!
      if (kprfit.gt.0.and.kdofit.gt.0) then
        nload=nload+1
        prbdry=brsp(nload)
      endif
      if (fitdelz.and.nniter.ge.ndelzon) then
        nload=nload+1
        cdelz(nniter)=brsp(nload)*scaledz*100.*relaxdz
      else
        cdelz(nniter)=0.0
      endif
      if (iacoil.eq.0) then
          do 2711 i=1,nacoil
            caccurt(jtime,i)=0.0
 2711     continue
      else
          do 2717 i=1,nacoil
            caccurt(jtime,i)=brsp(nload+i)
 2717     continue
          nload=nload+nacoil
      endif
!------------------------------------------------------------------------------
!--  E coil currents                                                         --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
       do m=1,nesum
          nload=nload+1
          cecurr(m)=brsp(nload)
       enddo
      endif
!------------------------------------------------------------------------------
!--  Reference flux                                                          --
!------------------------------------------------------------------------------
      if (fitsiref) then
          nload=nload+1
          csiref=brsp(nload)*scalesir
      endif
!------------------------------------------------------------------------------
!--  Hyperbolic tangents                                                     --
!------------------------------------------------------------------------------
      if (kedgep.gt.0) then
        pedge=brsp(npedge)
      endif
      if (kedgef.gt.0) then
        f2edge=brsp(nfedge)
      endif
!
 3000 continue
!----------------------------------------------------------------------------
!-- Update poloidal flux due to vertical shift if necessary                --
!----------------------------------------------------------------------------
      if (ifitdelz.eq.2) then
      if (fitdelz.and.nniter.ge.ndelzon) then
      cdelznow=cdelz(nniter)/100.
      call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
      do i=1,nw
      do j=1,nh
        kk=(i-1)*nh+j
        call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n333)
        psi(kk)=psi(kk)+cdelznow*pds(3)
      enddo
      enddo
      endif
      endif
!----------------------------------------------------------------------
!-- calculate the fitting figure of merit saisq                      --
!----------------------------------------------------------------------
      saisq=0.0
      do 4600 m=1,nsilop
        cm=0.0
        do 4520 n=1,nfcoil
          cm=cm+rsilfc(m,n)*brsp(n)
 4520   continue
        if (ivesel.le.0) go to 4550
        do 4545 n=1,nvesel
          cm=cm+rsilvs(m,n)*vcurrt(n)
 4545   continue
 4550   continue
        if (iecurr.ne.1) go to 4565
        do 4560 n=1,nesum
          cm=cm+rsilec(m,n)*ecurrt(n)
 4560   continue
 4565   continue
        if (iacoil.le.0) go to 4571
        do 4569 n=1,nacoil
          cm=cm+rsilac(m,n)*caccurt(jtime,n)
 4569   continue
 4571   continue
        if (iecurr.eq.2) then
          do n=1,nesum
            cm=cm+rsilec(m,n)*cecurr(n)
          enddo
        endif
        cmv=cm
        do 4540 n=1,kwcurn
          cm=cm+rsilpc(m,n)*brsp(nfcoil+n)
 4540   continue
        if (fitsiref) then
            cm=cm-csiref
        endif
        if (kedgep.gt.0) then
          cm=cm+rsilpe(m)*pedge
        endif
        if (kedgef.gt.0) then
          cm=cm+rsilfe(m)*f2edge
        endif
        if (swtsi(m).ne.0.0) then
        saisil(m)=fwtsi(m)**nsq*(silopt(jtime,m)-cm)**2
        saisil(m)=saisil(m)/swtsi(m)**nsq
        else
        saisil(m)=0.0
        endif
        saisq=saisq+saisil(m)
        csilop(m,jtime)=cm
        csilopv(m,jtime)=cmv
 4600 continue
!
      do 4700 m=1,magpri
        cm=0.0
        do 4620 n=1,nfcoil
          cm=cm+rmp2fc(m,n)*brsp(n)
 4620   continue
        if (ivesel.le.0) go to 4648
        do 4644 n=1,nvesel
          cm=cm+rmp2vs(m,n)*vcurrt(n)
 4644   continue
 4648   continue
        if (iecurr.ne.1) go to 4655
        do 4650 n=1,nesum
          cm=cm+rmp2ec(m,n)*ecurrt(n)
 4650   continue
 4655   continue
        if (iacoil.le.0) go to 4671
        do 4669 n=1,nacoil
          cm=cm+rmp2ac(m,n)*caccurt(jtime,n)
 4669   continue
 4671   continue
        if (iecurr.eq.2) then
          do n=1,nesum
            cm=cm+rmp2ec(m,n)*cecurr(n)
          enddo
        endif
        cmv=cm
        do 4640 n=1,kwcurn
          cm=cm+rmp2pc(m,n)*brsp(nfcoil+n)
 4640   continue
        if (kedgep.gt.0) then
          cm=cm+rmp2pe(m)*pedge
        endif
        if (kedgef.gt.0) then
          cm=cm+rmp2fe(m)*f2edge
        endif
        if (swtmp2(m).ne.0.0) then
        saimpi(m)=fwtmp2(m)**nsq*(expmpi(jtime,m)-cm)**2
        saimpi(m)=saimpi(m)/swtmp2(m)**nsq
        else
        saimpi(m)=0.0
        endif
        saisq=saisq+saimpi(m)
        cmpr2(m,jtime)=cm
        cmpr2v(m,jtime)=cmv
 4700 continue
!
      if (kstark.gt.0) then
      chigamt=0.0
      do 4800 m=1,nstark
        chigam(m)=0.0
        cmgam(m,jtime)=0.0
        if (rrgam(jtime,m).le.0.0) goto 4800
        cmbr=0.0
        cmbz=0.0
        do 4720 n=1,nfcoil
          cmbr=cmbr+rbrfc(m,n)*brsp(n)
          cmbz=cmbz+rbzfc(m,n)*brsp(n)
 4720   continue
        do 4740 n=1,kwcurn
          cmbr=cmbr+rbrpc(m,n)*brsp(nfcoil+n)
          cmbz=cmbz+rbzpc(m,n)*brsp(nfcoil+n)
 4740   continue
        if (ivesel.le.0) go to 4748
        do 4744 n=1,nvesel
          cmbr=cmbr+rbrvs(m,n)*vcurrt(n)
          cmbz=cmbz+rbzvs(m,n)*vcurrt(n)
 4744   continue
 4748   continue
        if (iecurr.ne.1) go to 4755
        do 4750 n=1,nesum
          cmbr=cmbr+rbrec(m,n)*ecurrt(n)
          cmbz=cmbz+rbzec(m,n)*ecurrt(n)
 4750   continue
 4755   continue
        if (iacoil.le.0) go to 4761
        do 4759 n=1,nacoil
          cmbr=cmbr+rbrac(m,n)*caccurt(jtime,n)
          cmbz=cmbz+rbzac(m,n)*caccurt(jtime,n)
 4759   continue
 4761   continue
        if (iecurr.eq.2) then
          do n=1,nesum
            cmbr=cmbr+rbrec(m,n)*cecurr(n)
            cmbz=cmbz+rbzec(m,n)*cecurr(n)
          enddo
        endif
        if (kedgep.gt.0) then
          cmbr=cmbr+rbrpe(m)*pedge
          cmbz=cmbz+rbzpe(m)*pedge
        endif
        if (kedgef.gt.0) then
          cmbr=cmbr+rbrfe(m)*f2edge
          cmbz=cmbz+rbzfe(m)*f2edge
        endif
        cm=a2gam(jtime,m)*btgam(m)+a3gam(jtime,m)*cmbr+a4gam(jtime,m) &
             *cmbz
        bzmsec(m)=cmbz
        if (keecur.le.0) then
          bzmse(m)=(tangam(jtime,m)*cm-a8gam(jtime,m)*cmbr) &
                       /a1gam(jtime,m)
          cm=(a1gam(jtime,m)*cmbz+a8gam(jtime,m)*cmbr)/cm
        else
          erup=0.0
          erbot=0.0
          do i=1,keecur
            cerer(i)=brsp(needer+i)
            erup=erup+e1rbz(m,i)*cerer(i)
            erbot=erbot+(e2rbz(m,i)+e3rbr(m,i))*cerer(i)
          enddo
          cm= cm-erbot
          bzmse(m)=(tangam(jtime,m)*cm+erup-a8gam(jtime,m)*cmbr) &
                      /a1gam(jtime,m)
          cm=(a1gam(jtime,m)*cmbz+a8gam(jtime,m)*cmbr-erup)/cm
          ermse(m)=-erup/a5gam(jtime,m)
        endif
        if (swtgam(m).ne.0.0) then
        chigam(m)=fwtgam(m)**nsq*(tangam(jtime,m)-cm)**2
        chigam(m)=chigam(m)/swtgam(m)**nsq
        else
        chigam(m)=0.0
        endif
        cmgam(m,jtime)=cm
 4800 continue
      do m=1,nmtark
        chigamt=chigamt+chigam(m)
      enddo
      chi2gamt(jtime)=chigamt
      chilibt=0.0
      do m=nmtark+1,nstark
        chilibt=chilibt+chigam(m)
      enddo
      if (ishot.le.97400) then
        mcentral=15
      else
        mcentral=10
      endif
      do m=1,mcentral
        mp1=m+1
        drgam=rrgam(jtime,mp1)-rrgam(jtime,m)
        cjmse(m)=-(bzmse(mp1)-bzmse(m))/drgam/twopi/tmu
        cjmsec(m)=-(bzmsec(mp1)-bzmsec(m))/drgam/twopi/tmu
      enddo
      endif
!
       tchimls=0.0
       do 94880 m=1,nmsels
        cm=0.0
        do 94840 n=1,kwcurn
          cm=cm+rmlspc(m,n)*brsp(nfcoil+n)
94840   continue
        cm=cm-rhsmls(jtime,m)
        cm=sqrt(cm)
        if (swtbmsels(m).ne.0.0) then
        chimls(m)=fwtbmselt(jtime,m)**nsq*(bmselt(jtime,m)-cm)**2
        chimls(m)=chimls(m)/swtbmsels(m)**nsq
        else
        chimls(m)=0.0
        endif
        cmmls(jtime,m)=cm
        tchimls=tchimls+chimls(m)
94880 continue
      chi2mls(jtime)=tchimls
!
       tchiels=0.0
       do 95880 m=1,nmsels
        cm=0.0
        do 95840 n=1,keecur
          cm=cm+relser(m,n)*cecurr(n)
95840   continue
        if (swtemsels(m).ne.0.0) then
        chiels(m)=fwtemselt(jtime,m)**nsq*(emselt(jtime,m)-cm)**2
        chiels(m)=chiels(m)/swtemsels(m)**nsq
        else
        chiels(m)=0.0
        endif
        cmels(jtime,m)=cm
        tchiels=tchiels+chiels(m)
95880 continue
!
       tchiece=0.0
       do 4880 m=1,nece
        cm=0.0
        do 4820 n=1,nfcoil
          cm=cm+recefc(m,n)*brsp(n)
 4820   continue
        do 4840 n=1,kwcurn
          cm=cm+recepc(m,n)*brsp(nfcoil+n)
 4840   continue
        if (ivesel.le.0) go to 4848
        do 4844 n=1,nvesel
          cm=cm+recevs(m,n)*vcurrt(n)
 4844   continue
 4848   continue
        if (iecurr.ne.1) go to 4855
        do 4850 n=1,nesum
          cm=cm+receec(m,n)*ecurrt(n)
 4850   continue
 4855   continue
        if (iacoil.le.0) go to 4871
        do 4869 n=1,nacoil
          cm=cm+receac(m,n)*caccurt(jtime,n)
 4869   continue
 4871   continue
        if (iecurr.eq.2) then
          do n=1,nesum
            cm=cm+receec(m,n)*cecurr(n)
          enddo
        endif
        if (swtece(m).ne.0.0) then
        chiece(m)=fwtece(m)**nsq*(brspece(jtime,m)-cm)**2
        chiece(m)=chiece(m)/swtece(m)**nsq
        else
        chiece(m)=0.0
        endif
        cmece(m,jtime)=cm
        tchiece=tchiece+chiece(m)
 4880 continue
!
        cm=0.0
        do 84820 n=1,nfcoil
          cm=cm+recebzfc(n)*brsp(n)
84820   continue
        do 84840 n=1,kwcurn
          cm=cm+recebzpc(n)*brsp(nfcoil+n)
84840   continue
        if (ivesel.le.0) go to 84848
        do 84844 n=1,nvesel
          cm=cm+recevs(m,n)*vcurrt(n)
84844   continue
84848   continue
        if (iecurr.ne.1) go to 84855
        do 84850 n=1,nesum
          cm=cm+recebzec(n)*ecurrt(n)
84850   continue
84855   continue
        if (iacoil.le.0) go to 84871
        do 84869 n=1,nacoil
          cm=cm+receac(m,n)*caccurt(jtime,n)
84869   continue
84871   continue
        if (iecurr.eq.2) then
          do n=1,nesum
            cm=cm+recebzec(n)*cecurr(n)
          enddo
        endif
        if (swtecebz.ne.0.0) then
        chiecebz=fwtecebz**nsq*(brspecebz(jtime)-cm)**2
        chiecebz=chiecebz/swtecebz**nsq
        else
        chiecebz=0.0
        endif
        cmecebz(jtime)=cm
!
      cm=0.0
      do 4900 n=nfcoil+1,nfnwcr
        cm=cm+brsp(n)*fgowpc(n-nfcoil)
 4900 continue
        if (kedgep.gt.0) then
          cm=cm+fgowpe*pedge
        endif
        if (kedgef.gt.0) then
          cm=cm+fgowfe*f2edge
        endif
      cpasma(jtime)=cm
      if (kfffnc.eq.8) then
        cjeccd=brsp(nfnwcr)*fgowpc(nfnwcr-nfcoil)/1000.
      else
        cjeccd=0.0
      endif
      if (ifitvs.gt.0) then
        do 30500 i=1,nvesel
          cm=cm+vcurrt(i)
30500   continue
      endif
      if (swtcur.ne.0.0) then
      saiip=(fwtcur/swtcur)**nsq*(pasmat(jtime)-cm)**2
      else
      saiip=0.0
      endif
      saisq=saisq+saiip
!
      tsaifc=0.0
      do 30510 i=1,nfcoil
        saifc(i)=0.0
        if (fwtfc(i).gt.0.0) then
          saifc(i)=fwtfc(i)**nsq*(brsp(i)-fccurt(jtime,i))**2
          saifc(i)=saifc(i)/swtfc(i)**nsq
        endif
        saisq=saisq+saifc(i)
        tsaifc=tsaifc+saifc(i)
30510 continue
!
      if (iecurr.eq.2) then
      do i=1,nesum
        saiec(i)=0.0
        if (fwtec(i).gt.0.0) then
          saiec(i)=fwtec(i)**nsq*(cecurr(i)-ecurrt(i))**2
          saiec(i)=saiec(i)/swtec(i)**nsq
        endif
        saisq=saisq+saiec(i)
      enddo
      endif
!
      if (fitsiref) then
        saisref=0.0
        if (fwtref.gt.0.0) then
          saisref=fwtref**nsq*(psiref(jtime)-csiref)**2
          saisref=saisref/swtsi(nslref)**nsq
        endif
        saisq=saisq+saisref
      endif
      tsaisq(jtime)=saisq
!
      chipre=0.0
      if (kprfit.le.0.or.kdofit.eq.0) go to 4910
      if (npress.le.0) goto 4910
      do 4908 m=1,npress
        cm=0.0
        do 4906 n=1,kppcur
          cm=cm+rprepc(m,n)*brsp(nfcoil+n)
 4906   continue
        if (kedgep.gt.0) then
          cm=cm+rprepe(m)*pedge
        endif
        cm=cm+prbdry
        if (fwtpre(m).gt.0.0) then
        saipre(m)=((cm-pressr(m))/sigpre(m))**2
        else
        saipre(m)=0.0
        endif
        saipre2(m)=saipre(m)  ! preserve saipre - changed later by pltout
        chipre=chipre+saipre(m)
        precal(m)=cm
 4908 continue
 4910 continue
!
      chiprw=0.0
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          cm=0.0
          do n=1,kwwcur
            cm=cm+rprwpc(m,n)*brsp(nfnpcr+n)
          enddo
        cm=cm+preswb
        if (fwtprw(m).gt.0.0) then
        saiprw(m)=((cm-presw(m))/sigprw(m))**2
        else
        saiprw(m)=0.0
        endif
        chiprw=chiprw+saiprw(m)
        prwcal(m)=cm
        enddo
      endif
!
      if ((nniter.lt.minite).and.((eouter.gt.elomin).or.(fwtdlc.gt.0.0) &
      ))  go to 4950
      if ((nniter.lt.kcallece).and.(kfitece.gt.0.0))  go to 4950
      if ((errorm.gt.errmin).and.((eouter.gt.elomin).or.(fwtdlc.gt.0.0) &
      ))  go to 4950
      if (saisq.gt.saicon) go to 4950
      if (iconvr.ne.2) go to 4950
      if (abs(saisq-saiold).le.0.10) go to 4918
      if (saisq.lt.saiold) go to 4950
 4918 continue
      ichisq=1
      do 4920 i=1,nrsmat
        brsp(i)=brsold(i)
 4920 continue
      if (ifitvs.gt.0) then
        do 30520 i=1,nvesel
          vcurrt(i)=vcurrto(i)
30520   continue
      endif
      saisq=saiold
 4950 continue
!
      if (iand(iout,1).ne.0) then
      write (nout,7400) time(jtime),chipre,cpasma(jtime), &
                        nniter,condno,saisq,chigamt
      write (nout,97400) time(jtime),chimlst,chielst
      write (nout,7445) need
      write (nout,7450) (brsp(i),i=1,need)
      write (nout,7450) (wrsp(i),i=1,need)
      endif
!
      if (iupdat.gt.0) return
      if (saisq.gt.saimin) return
      tcrrnt=cpasma(jtime)
      iupdat=1
      return
!
 6000 continue
      return
 7400 format (/,2x,7htime = ,e12.5,2x,8hchipr = ,e12.5, &
              2x,10hcurrent = ,e12.5,/,2x,5hit = ,i5, &
              1x,8hcondn = ,1pe11.4, &
              1x,8hchisq = ,1pe11.4,1x,9hchigam = ,1pe11.4, &
              /,2x,11hchiecebz = ,1pe11.4,1x,14htchieceR+R- = , &
              1pe11.4)
97400 format (/,2x,7htime = ,e12.5,2x,8hchimls= ,e12.5, &
              2x,10hchiels = ,e12.5)
 7445 format (10x,22hfitting parameters:   ,i5)
 7450 format (8(1x,e12.5,1x))
 7460 format (10x,7hchi ip:,/,15x,e12.5)
 8000 format (/,'  ** Problem in Decomposition **',i10)
      end
