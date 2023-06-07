#include "config.f"
!**********************************************************************
!>
!!    fit carries out the fitting and equilibrium iterations.
!!
!!    @param jtime : time index
!!
!!    @param kerror : error flag
!!
!**********************************************************************
      subroutine fit(jtime,kerror)
      include 'eparm.inc'
      include 'modules1.inc'
      use var_cvesel, only: vcurrt
      implicit none
      integer*4, intent(in) :: jtime
      integer*4, intent(out) :: kerror
      integer*4 i,ii,ix,idone,iend,ichisq,iwantk,nleft
!-----------------------------------------------------------------------
!--   inner equilibrium do loop and outer current profile parameters
!--   do loop
!-----------------------------------------------------------------------
#ifdef DEBUG_LEVEL1
      write (6,*) 'Enter FIT'
#endif
      kerror=0
      jerror(jtime)=0
      nitera=0
      ndelzon=999
      iwantk=0
      cjeccd=0.0
      lflag=0
      iend=mxiter+1
      if (iconvr.eq.3) iend=1
      current_profile: do i=1,iend
        ix=i
        if (i.gt.1) then
          iwantk=iwantk+1
!------------------------------------------------------------------
!--       nitera.ge.kcallece, then call setece
!--       mtxece<=1 call setece every nconstr iterations
!--       mtxece>1 call setece every mtxece*nconstr time iterations
!------------------------------------------------------------------
          if ((nitera.ge.kcallece).and.(kfitece.gt.0)) then
            nleft=abs(mxiter)-nitera
            if(nleft .ge. mtxece*nconstr) then
#ifdef DEBUG_LEVEL2
              write (6,*) 'Call FIT/SETECE',ix
              write (6,*) &
                '  nitera/kcallece/kfitece/nleft/mtxece/nconstr = ', &
                nitera,kcallece,kfitece,nleft,mtxece,nconstr
#endif
              call setece(jtime,kerror)
              if (kerror /= 0) then
                jerror(jtime) = 1
                return
              endif
            endif
          endif
          call green(0,jtime,nitera)
          if (kprfit.gt.0.and.iwantk.eq.ndokin) then
            call presur(jtime,nitera,kerror)
            if (kerror /= 0) then
              jerror(jtime) = 1
              return
            endif
            iwantk=0
          endif
          if (kprfit.ge.3) call presurw(jtime,nitera)
#ifdef DEBUG_LEVEL2
          write (6,*) 'Call FIT/RESPONSE_MATRIX',ix
#endif
          call response_matrix(jtime,ix,ichisq,nitera,kerror)
           
          if (kerror /= 0) then
            jerror(jtime) = 1
            return
          endif
          if ((iconvr.eq.2).and.(ichisq.gt.0)) then
            call errctrl_msg('fit','iconvr=2 satisfied, exiting',2)
            go to 2020
          end if
        end if

        equilibrium: do ii=1,nxiter
          ixnn=ii
          nitera=nitera+1

#ifdef DEBUG_LEVEL2
          write(6,*) 'Entering currnt'
#endif
          call currnt(ix,jtime,nitera,kerror)
          if (kerror /= 0) then
            jerror(jtime) = 1
            return
          endif
           
          if(ivesel.eq.2) vcurrt=vloopt(jtime)/rsisvs
          if ((i.le.1).or.(ii.gt.1)) then
#ifdef DEBUG_LEVEL2
            write(6,*) 'Entering external_current'
#endif
            call external_current(jtime,ix,nitera,kerror)
          endif 
          if (kerror /= 0) then
            jerror(jtime) = 1
            return
          endif
#ifdef DEBUG_LEVEL2
          write(6,*) 'Entering pflux'
#endif
          call pflux(ix,ixnn,nitera,jtime,kerror)
          if (kerror.gt.0) then
            jerror(jtime) = 1
            return
          endif

#ifdef DEBUG_LEVEL2
          write(6,*) 'Entering update_params'
#endif
          call update_params(ixnn,nitera,ix,jtime,kerror)
          if (kerror.gt.0) then
            jerror(jtime) = 1
            return
          endif

          if (kmtark.gt.0) then
            if (kwaitmse.ne.0 .and. i.ge.kwaitmse) then
#ifdef DEBUG_LEVEL2
              write(6,*) 'Entering fixstark'
#endif
              call fixstark(jtime,kerror)
              if (kerror.gt.0) then
                jerror(jtime) = 1
                return
              endif
            end if
          endif

#ifdef DEBUG_LEVEL2
          write(6,*) 'Entering residu'
#endif
          call residu(nitera,jtime,idone)
          if((nitera.lt.kcallece).and.(kfitece.gt.0.0)) exit
          if((ii.eq.1).and.(idone.gt.0).and.(chisq(jtime).le.saimin)) &
            go to 2020
          if(idone.gt.0) exit
          if(i.eq.mxiter+1) exit
        end do equilibrium
      end do current_profile
      if ((nbdry.le.0).and.(ivacum.eq.0).and.(iconvr.ne.4)) then
        call errctrl_msg('fit','not converged, reached max iterations',2)
        lflag=1
      endif
2020  continue
      terror(jtime)=errorm
!---------------------------------------------------------------------
!--    update pressure if needed                                    --
!---------------------------------------------------------------------
      if (kprfit.gt.1) then
#ifdef DEBUG_LEVEL2
        write(6,*) 'Entering presur'
#endif
        call presur(jtime,nitera,kerror)
        if (kerror /= 0) then
          jerror(jtime) = 1
          return
        endif
      end if
      return
      end subroutine fit

!**********************************************************************
!>
!!    residu computes the flux variations on the r-z grid.
!!
!!
!!    @param nx : iteration number
!!
!!    @param jtime : time index
!!
!!    @param idone : flag for whether convergence criteria is met
!!
!**********************************************************************
      subroutine residu(nx,jtime,idone)
      use commonblocks,only: psiold,psipold
      include 'eparm.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: nx,jtime
      integer*4, intent(out) :: idone
      integer*4 i,j,kk
      real*8 errold,errave,change

      idone=0
      if(ivacum.eq.1) return
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

      errorm=errorm/abs(sidif)/relax

      aveerr(nx)=errave/abs(sidif)/nwnh
      cerror(nx)=errorm
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
              write (nttyo,10019) rank,itime,nx,chisq(jtime),zmaxis, &
                errorm,delzmm,sum(chigam)
            else
              write (nttyo,90019) rank,itime,nx,chisq(jtime),zmaxis, &
                errorm,delzmm,sum(chigam),tchimls
            endif
          else
            write (nttyo,10020) rank,itime,nx,chisq(jtime),zmaxis, &
              errorm,delzmm,sum(chigam),erbmax,erbsmax
          endif
        elseif (itell.eq.2) then
          write (nttyo,10021) rank,itime,nx,li(jtime),abs(betatn), &
            errorm,qsiw(1)
        elseif (itell.eq.3) then
          write (nttyo,10023) rank,itime,nx,difpsi,zmaxis,errorm,delzmm
        elseif (itell.eq.4) then
          if (nx.eq.1) then
            !write (nttyo,10017) itime, rank
            write (nttyo,'(x)')
            if (kstark.gt.0) then
              write (nttyo,80019) rank,itime,nx,chisq(jtime),zmaxis, &
                errorm,delzmm,cdelz(1),cdeljsum,sum(chigam)
            else
              write (nttyo,10025) rank,itime,nx,chisq(jtime),zmaxis, &
                errorm,delzmm,cdelz(1),cdeljsum
            endif
          else
            if (kstark.gt.0) then
              write (nttyo,80019) rank,itime,nx,chisq(jtime),zmaxis, &
                errorm,delzmm,cdelz(nx-1),cdeljsum,sum(chigam)
            else
              write (nttyo,10025)  rank,itime,nx,chisq(jtime),zmaxis, &
                errorm,delzmm,cdelz(nx-1),cdeljsum
            endif
          endif
        endif
      endif
      call flush(6)
      if (isetfb.ne.0) then
        write (4,10009) rank,itime,nx,chisq(jtime),zmaxis,errorm, &
          delzmm,brfb(1)
        if (isetfb.lt.0) &
          write (6,10009) rank,itime,nx,chisq(jtime),zmaxis,errorm, &
            delzmm,brfb(1)
      elseif (eelip.gt.2.25_dp .and. itell.eq.0) then
        write (6,10009) rank,itime,nx,chisq(jtime),zmaxis,errorm, &
          delzmm,brfb(1)
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
