!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          fit carries out the fitting and equilibrium iterations. **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**          94/03/11..........revised                               **
!**                                                                  **
!**********************************************************************
      subroutine fit(jtime,kerror)
      include 'eparmdud129.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      data nzero/0/
      save nzero
!-----------------------------------------------------------------------
!--  inner equilibrium do loop and outer current profile parameters   --
!--  do loop                                                          --
!-----------------------------------------------------------------------
      if (idebug /= 0) write (6,*) 'Enter FIT'
      kerror=0
      jerror(jtime)=0
      nitera=0
      ndelzon=999
      iwantk=0
      iend=mxiter+1
      if (iconvr.eq.3) iend=1
      idosigma=1
      do i=1,iend
        ix=i
        if (i.gt.1) then
          iwantk=iwantk+1
          !------------------------------------------------------------------
          !--  nitera.ge.kcallece, then call setece                        --
          !--  mtxece=1 call setece every iteration                        --
          !--  mtxece>1 call setece per mtece time iteration               --
          !------------------------------------------------------------------
          if ((nitera.ge.kcallece).and.(kfitece.gt.0)) then
            nleft=abs(mxiter)-nitera
            if(nleft .ge. mtxece*nconstr) then
              if (idebug.ge.2) then
                write (6,*) 'Call FIT/SETECE',ix
                write (6,*) '  nitera/kcallece/kfitece/nleft/mtxece/nconstr = ', &
                  nitera,kcallece,kfitece,nleft,mtxece,nconstr
              endif
              call setece(jtime,kerror)
              if (kerror /= 0) then
                jerror(jtime) = 1
                return
              endif
            endif
          endif
          call green(nzero,jtime,nitera)
          if (kprfit.gt.0.and.iwantk.eq.ndokin) then
            call presur(jtime,nitera,kerror)
            if (kerror /= 0) then
              jerror(jtime) = 1
              return
            endif
            iwantk=0
          endif
          if (kprfit.ge.3) call presurw(jtime,nitera)
          if (errorm.lt.errmagb) then
            if ((imagsigma.eq.1) .AND. (errorm > errmag) ) &
              call getsigma(jtime,nitera)
            if ((imagsigma.eq.2).and.(idosigma.eq.1)) then
              call getsigma(jtime,nitera)
              idosigma=2
            endif
          endif
          if (idebug.ge.2) write (6,*) 'Call FIT/MATRIX',ix
          call matrix(jtime,ix,ichisq,nitera,kerror)
           
          if (kerror /= 0) then
            jerror(jtime) = 1
            return
          endif
          if ((iconvr.eq.2).and.(ichisq.gt.0)) then
            call errctrl_msg('fit','not converged properly',2)
            go to 2020
          end if
        end if

        do in=1,nxiter
          ixnn=in
          nitera=nitera+1

          if (idebug.ge.2) write(6,*) 'Entering currnt' 
          call currnt(ix,jtime,ixnn,nitera,kerror)
          if (kerror /= 0) then
            jerror(jtime) = 1
            return
          endif
           
          if (ivesel.ge.2) then
             if (idebug.ge.2) write(6,*) 'Entering vescur'
             call vescur(jtime)
          endif
          if ((i.le.1).or.(in.gt.1)) then
              if (idebug.ge.2) write(6,*) 'Entering fcurrt'
             call fcurrt(jtime,ix,nitera,kerror)
          endif 
          if (kerror /= 0) then
            jerror(jtime) = 1
            return
          endif
          if (idebug.ge.2) write(6,*) 'Entering pflux' 
          call pflux(ix,ixnn,nitera,jtime,kerror)
          if (kerror.gt.0) then
            jerror(jtime) = 1
            return
          endif

          if (idebug.ge.2) write(6,*) 'Entering steps'
          call steps(ixnn,nitera,ix,jtime,kerror)
          if (kerror.gt.0) then
            jerror(jtime) = 1
            return
          endif
          if (kmtark.gt.0) then
            if (kwaitmse.ne.0 .and. i.ge.kwaitmse) then
              if (idebug.ge.2) write(6,*) 'Entering fixstark'
              call fixstark(jtime,kerror)
              if (kerror.gt.0) then
                jerror(jtime) = 1
                return
              endif
            end if
          endif

          if (idebug.ge.2) write(6,*) 'Entering residu'
          call residu(nitera,jtime)
          if ((nitera.lt.kcallece).and.(kfitece.gt.0.0)) exit
          if ((in.eq.1).and.(idone.gt.0).and.(tsaisq(jtime).le.saimin)) then
            go to 2020
          end if
          if (idone.gt.0) exit
          if (i.eq.mxiter+1) exit
        end do ! in
      end do ! i
      call errctrl_msg('fit','not converged, reached max iterations',2)
2020  continue
      !---------------------------------------------------------------------
      !--  update pressure if needed                                      --
      !---------------------------------------------------------------------
      if (kprfit.gt.1) then

        if (idebug.ge.2) write(6,*) 'Entering presur'
        call presur(jtime,nitera,kerror)
        if (kerror /= 0) then
          jerror(jtime) = 1
          return
        endif
      end if
      return
      end subroutine fit

!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          chisqr computes the figure of merit for fitting         **
!**          chisq.                                                  **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      subroutine chisqr(jtime)
      include 'eparmdud129.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
!
      if (idebug>=2) then
        write (6,*) 'CHISQR, jtime= ',jtime
      endif
!
      if (ivesel.gt.10) return
      if (nbdry.le.0) go to 200
!----------------------------------------------------------------------
!-- read in the plasma response function                             --
!----------------------------------------------------------------------

      open(unit=nrsppc,status='old',form='unformatted', &
           file=table_di2(1:ltbdi2)//'ep'//trim(ch1)//trim(ch2)//'.ddd')
      read (nrsppc) gsilpc
      read (nrsppc) gmp2pc
      close(unit=nrsppc)
!
  200 continue
      nsq=2
      saisq=0.0
      chi2rm(jtime)=0.0
      do 4600 m=1,nsilop
        cm=0.0
        do 4520 n=1,nfcoil
          cm=cm+rsilfc(m,n)*brsp(n)
 4520   continue
        do 4540 n=1,nwnh
          cm=cm+gsilpc(m,n)*pcurrt(n)
 4540   continue
        if (ivesel.le.0) go to 4548
        do 4544 n=1,nvesel
          cm=cm+rsilvs(m,n)*vcurrt(n)
 4544   continue
 4548   continue
        if (iecurr.ne.1) go to 4560
        do 4550 n=1,nesum
          cm=cm+rsilec(m,n)*ecurrt(n)
 4550   continue
 4560   continue
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
        if (fitsiref) then
            cm=cm-csiref
        endif
        if (swtsi(m).ne.0.0) then
        saisil(m)=(fwtsi(m)/swtsi(m))**nsq*(silopt(jtime,m)-cm)**2
        else
        saisil(m)=0.0
        endif
        saisq=saisq+saisil(m)
        csilop(m,jtime)=cm
        chi2rm(jtime)=max(chi2rm(jtime),saisil(m))
 4600 continue
!
      do 4700 m=1,magpri
        cm=0.0
        do 4620 n=1,nfcoil
          cm=cm+rmp2fc(m,n)*brsp(n)
 4620   continue
        do 4640 n=1,nwnh
          cm=cm+gmp2pc(m,n)*pcurrt(n)
 4640   continue
        if (ivesel.le.0) go to 4648
        do 4644 n=1,nvesel
          cm=cm+rmp2vs(m,n)*vcurrt(n)
 4644   continue
 4648   continue
        if (iecurr.ne.1) go to 4660
        do 4650 n=1,nesum
          cm=cm+rmp2ec(m,n)*ecurrt(n)
 4650   continue
 4660   continue
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
        if (swtmp2(m).ne.0.0) then
        saimpi(m)=(fwtmp2(m)/swtmp2(m))**nsq*(expmpi(jtime,m)-cm)**2
        else
        saimpi(m)=0.0
        endif
        saisq=saisq+saimpi(m)
        cmpr2(m,jtime)=cm
        chi2rm(jtime)=max(chi2rm(jtime),saimpi(m))
 4700 continue
!-------------------------------------------------------------------
!--    calculate ECE chisqr (chiece, tchiece for psi(R+)=psi(R-))
!-------------------------------------------------------------------
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
!-------------------------------------------------------------------
!--    calculate ECE chisqr (chiecebz for Bz(receo)=0)
!-------------------------------------------------------------------
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
      do 4900 n=1,nwnh
        cm=cm+pcurrt(n)
 4900 continue
      cpasma(jtime)=cm
      if (ifitvs.gt.0.or.ivesel.gt.0) then
        do 4910 i=1,nvesel
          cm=cm+vcurrt(i)
 4910   continue
      endif
      if (swtcur.ne.0.0) then
      saiip=(fwtcur/swtcur)**nsq*(pasmat(jtime)-cm)**2
      else
      saiip=0.0
      endif
      saisq=saisq+saiip
      chi2rm(jtime)=max(chi2rm(jtime),saiip)
!
      tsaifc=0.0
      do 4920 i=1,nfcoil
        saifc(i)=0.0
        if (fwtfc(i).gt.0.0) then
          saifc(i)=fwtfc(i)**nsq*(brsp(i)-fccurt(jtime,i))**2
          saifc(i)=saifc(i)/swtfc(i)**nsq
        endif
        saisq=saisq+saifc(i)
        tsaifc=tsaifc+saifc(i)
        chi2rm(jtime)=max(chi2rm(jtime),saifc(i))
 4920 continue
      if (iecurr.eq.2) then
      do i=1,nesum
        saiec(i)=0.0
        if (fwtec(i).gt.0.0) then
          saiec(i)=fwtec(i)**nsq*(cecurr(i)-ecurrt(i))**2
          saiec(i)=saiec(i)/swtec(i)**nsq
        endif
        saisq=saisq+saiec(i)
        chi2rm(jtime)=max(chi2rm(jtime),saiec(i))
      enddo
      endif
      if (fitsiref) then
        saisref=0.0
        if (fwtref.gt.0.0) then
          saisref=fwtref**nsq*(psiref(jtime)-csiref)**2
          saisref=saisref/swtsi(nslref)**nsq
        endif
        saisq=saisq+saisref
        chi2rm(jtime)=max(chi2rm(jtime),saisref)
      endif
!
      do 4950 n=1,nfcoil
        ccbrsp(n,jtime)=brsp(n)
 4950 continue
!
      tsaisq(jtime)=saisq
      if (iand(iout,1).ne.0) then
      write (nout,7400) time(jtime),tsaisq(jtime),cpasma(jtime)
      write (nout,7420)
      write (nout,7450) (saisil(m),m=1,nsilop)
      write (nout,7430)
      write (nout,7450) (saimpi(m),m=1,magpri)
      write (nout,7460) saiip
      write (nout,7480) tsaifc
      write (nout,7450) (saifc(m),m=1,nfcoil)
      write (nout,7482) saisref
      write (nout,7485)
      write (nout,7450) (saiec(m),m=1,nesum)
      if (kprfit.gt.0) write (nout,7470) chipre
      if (kprfit.gt.0) write (nout,7450) (saipre(m),m=1,npress)
      if (kecebz.gt.0) write (nout,7486) chiecebz
      if (kece.gt.0) write (nout,7487) tchiece
      if (kece.gt.0) then
       write(nout,7488)
       write (nout,7450) (chiece(m),m=1,nece)
      endif
      endif
!-------------------------------------------------------------------------
!-- compute signals at MSE locations if requested                       --
!-------------------------------------------------------------------------
      if (kdomse.gt.0) then
      nnn=0
      call green(nnn,jtime,n222)
      do 4800 m=1,nstark
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
        if (kedgep.gt.0) then
          cmbr=cmbr+rbrpe(m)*pedge
          cmbz=cmbz+rbzpe(m)*pedge
        endif
        if (kedgef.gt.0) then
          cmbr=cmbr+rbrfe(m)*f2edge
          cmbz=cmbz+rbzfe(m)*f2edge
        endif
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
        if (iecurr.eq.2) then
          do n=1,nesum
            cmbr=cmbr+rbrec(m,n)*cecurr(n)
            cmbz=cmbz+rbzec(m,n)*cecurr(n)
          enddo
        endif
        if (iacoil.le.0) go to 4761
        do 4759 n=1,nacoil
          cmbr=cmbr+rbrac(m,n)*caccurt(jtime,n)
          cmbz=cmbz+rbzac(m,n)*caccurt(jtime,n)
 4759   continue
 4761   continue
        cm=a2gam(jtime,m)*btgam(m)+a3gam(jtime,m)*cmbr+a4gam(jtime,m) &
             *cmbz
        bzmsec(m)=cmbz
        if (keecur.le.0) then
        cm=a1gam(jtime,m)*cmbz/cm
        else
        ce1rbz=0.0
        ce2rbz=0.0
        ce3rbr=0.0
        do n=1,keecur
         ce1rbz=ce1rbz+e1rbz(m,n)*cerer(n)
         ce2rbz=ce2rbz+e2rbz(m,n)*cerer(n)
         ce3rbr=ce3rbr+e3rbr(m,n)*cerer(n)
        enddo
        cm=cm-ce2rbz-ce3rbr
        cm=(a1gam(jtime,m)*cmbz-ce1rbz)/cm
        endif
        cmgam(m,jtime)=cm
 4800 continue
      endif
!-------------------------------------------------------------------------
!-- compute signals at MSE-LS locations if requested                    --
!-------------------------------------------------------------------------
      if (kdomsels.gt.0) then
      call green(nnn,jtime,n222)
      do 94800 m=1,nmsels
        cmbr=0.0
        do 94720 n=1,kpcurn
          cmbr=cmbr+rmlspc(m,n)*brsp(nfcoil+n)
94720   continue
        cmbr=cmbr-rhsmls(jtime,m)
        cmmls(jtime,m)=sqrt(cmbr)
        cmmls2(jtime,m)=l1mselt(jtime,m)*btmls(m)**2+l2mselt(jtime,m)* &
           brmls(m)*btmls(m)+l3mselt(jtime,m)*brmls(m)**2+             &
          (l1mselt(jtime,m)+l3mselt(jtime,m))*bzmls(m)**2
        cmmls2(jtime,m)=sqrt(cmmls2(jtime,m))
94800 continue
      write (nout,7585)
      write (nout,7450) (cmmls(jtime,m),m=1,nmsels)
      write (nout,7588)
      write (nout,7450) (cmmls2(jtime,m),m=1,nmsels)
      endif
!
      if (mmbmsels.gt.0.or.kdomsels.gt.0) then
        do m=1,nmsels
          cmmlsv(jtime,m)=0.0
          if (rrmselt(jtime,m).gt.0.0) then
            cmmlsv(jtime,m)=abs(bcentr(jtime)*rcentr/rrmselt(jtime,m))
            cmmlsv(jtime,m)=cmmlsv(jtime,m)*sqrt(l1mselt(jtime,m))
          endif
        enddo 
        write (nout,7590)
        write (nout,7450) (cmmlsv(jtime,m),m=1,nmsels)
      endif
!
      return
 7400 format (/,2x,'time = ',e12.5,2x,'chisq = ',e12.5,2x,'current = ',e12.5)
 7420 format (10x,'chi psi loops:')
 7430 format (10x,'chi inner magnetic probes:')
 7450 format (8(1x,e12.5:,1x))
 7460 format (10x,'chi ip:',/,15x,e12.5)
 7470 format (10x,'chi pressure:         ',/,1x,e12.5)
 7480 format (10x,'chi F-coils:          ',/,10x,e12.5)
 7482 format (10x,'chi psiref:',/,15x,e12.5)
 7485 format (10x,'chi E-coils:          ',/,1x,e12.5)
 7486 format (10x,'chi ecebz:            ',/,1x,e12.5)
 7487 format (10x,'chi total eceR+R-:    ',/,1x,e12.5)
 7488 format (10x,'chi eceR+R-:          ')
 7585 format (10x,'Simulated MSE-LS (T): ',/)
 7588 format (10x,'Simulated MSE-LS2 (T): ',/)
 7590 format (10x,'Simulated MSE-LSV (T): ',/)
      end subroutine chisqr



