#include "config.f"
!**********************************************************************
!*
!>          FLUXAV does the flux surface average
!!   
!!
!!     @param f : function array 
!!
!!     @param si : flux array
!!
!!     @param n : length of contour array
!!
!!     @param x : R coordinates of contour
!!
!!     @param y : Z coordinates of contour     
!!
!!     @param fave  : int(dl f/Bp) / sdlobp
!!
!!     @param sdlbp : int(dl Bp)
!!
!!     @param sdlobp: int(dl/Bp)
!!
!!     @param rx: 
!!
!!     @param ry: 
!!
!!     @param msx: 
!!
!!     @param msy: 
!!
!!     @param ns: 
!!  
!**********************************************************************
      subroutine fluxav(f,x,y,n,si,rx,msx,ry,msy,fave,ns,sdlobp,sdlbp)
      use set_kinds, only: dp
      use commonblocks,only: c,wk,bkx,bky
      include 'eparm.inc'
      use var_cwork3, only:lkx,lky
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      dimension f(n),x(n),y(n),si(n),pds(6),rx(n),ry(n)
      if (ns.ne.0) then
        call sets2d(si,c,rx,msx,bkx,lkx,ry,msy,bky,lky,wk,ier)
      endif
!------------------------------------------------------------------
!--   flux surface average of f                                   --
!------------------------------------------------------------------
      fave=0.0
      fnorm=0.0
      sdlbp=0.0
      do i=2,n
        xnow=0.5_dp*(x(i-1)+x(i))
        ynow=0.5_dp*(y(i-1)+y(i))
        fnow=0.5_dp*(f(i-1)+f(i))
        dxnow=x(i)-x(i-1)
        dynow=y(i)-y(i-1)
        dl=sqrt(dxnow**2+dynow**2)
        call seva2d(bkx,lkx,bky,lky,c,xnow,ynow,pds,ier,n333)
        bpol = sqrt(pds(2)**2+pds(3)**2)/xnow
        dlbpol = dl/bpol
        fnorm = fnorm + dlbpol
        fave = fave + dlbpol*fnow
        sdlbp = sdlbp + dl*bpol
      enddo
      fave = fave/fnorm
      sdlobp = fnorm
      return
      end subroutine fluxav


!**********************************************************************
!>
!!    this subroutine does ...
!!    
!!
!!    @param is :
!!
!!    @param rs :
!!
!!    @param zs :
!!
!!    @param cs :
!!
!!    @param rc :
!!
!!    @param zc :
!!
!!    @param wc :
!!
!!    @param hc :
!!
!!    @param ac :
!!
!!    @param ac2 :
!!
!!    @param cc :
!!
!**********************************************************************
      subroutine splitc(is,rs,zs,cs,rc,zc,wc,hc,ac,ac2,cc)
      use global_constants
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension rs(is*is),zs(is*is),cs(is*is)

!----------------------------------------------------------------------
!--   rectangle                                                        --
!----------------------------------------------------------------------
      if (ac+ac2.eq.0.) then
        wdelt=wc/is
        hdelt=hc/is
        rstrt=rc-wc/2.+wdelt/2.
        zstrt=zc-hc/2.+hdelt/2.
        zz=zstrt
        ic=0
        c=cc/(is*is)
        do ii=1,is
          rr=rstrt
          do jj=1,is
            ic=ic+1
            zs(ic)=zz
            rs(ic)=rr
            cs(ic)=c
            rr=rr+wdelt
          enddo
          zz=zz+hdelt
        enddo
!----------------------------------------------------------------------
      elseif (ac.ne.0.) then
        side=tan(radeg*ac)*wc
        hdelt=hc/is
        wdelt=wc/is
        zdelt=tan(radeg*ac)*wdelt
        rstrt=rc-wc/2.+wdelt/2.
        tsid=hc+side
        zstrt =zc-tsid/2.+tsid/2.*1./is
        rr=rstrt
        ic=0
        c=cc/(is*is)
        do ii=1,is
          zz=zstrt+(ii-1)*zdelt
          do jj=1,is
            ic=ic+1
            zs(ic)=zz
            rs(ic)=rr
            cs(ic)=c
            zz=zz+hdelt
          enddo
          rr=rr+wdelt
        enddo
!----------------------------------------------------------------------
      elseif (ac2.ne.0.) then
        side=hc/tan(radeg*ac2)
        hdelt=hc/is
        wdelt=wc/is
        zstrt=zc-hc/2.+hdelt/2.
        rdelt=hdelt/tan(radeg*ac2)
        rstrt=rc-side/2.-wc/2.+rdelt/2.+wdelt/2.
        side=hc/tan(radeg*ac2)
        wtot=side+wc
        whaf=(side+wc)/2.
        rcorn=rc-whaf
        rcornr=rc+whaf
        rcorn2=rcorn+wtot/is
        rstrt=(rcorn+rcorn2)/2.
        zz=zstrt
        ic=0
        c=cc/(is*is)
        do ii=1,is
          rr=rstrt+(ii-1)*rdelt
          do jj=1,is
            ic=ic+1
            zs(ic)=zz
            rs(ic)=rr
            cs(ic)=c
            rr=rr+wdelt
          enddo
          zz=zz+hdelt
        enddo
      endif
      return
      end subroutine splitc

      
!**********************************************************************
!>
!!    this subroutine reorders the z profile data to be
!!    in ascending order and sets the ne and te data to
!!    correspond to the new order.
!!    
!!
!!    @param mbox :
!!
!!    @param zprof :
!!
!!    @param nemprof :
!!
!!    @param temprof :
!!
!!    @param nerprof :
!!
!!    @param terprof :
!!
!**********************************************************************
      subroutine tsorder(mbox,zprof,nemprof,temprof,nerprof,terprof)
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension ztemp(40),temtemp(40),tertemp(40),zinc(40), &
                  zprof(1), temprof(1), terprof(1)
      real*8 nemtemp(40),nertemp(40),nemprof(1),nerprof(1)
!---------------------------------------------------------------------
!--   copy zprof to ztemp (temporary work space)                    --
!---------------------------------------------------------------------
      do j=1,mbox
        ztemp(j)=zprof(j)
      enddo
!---------------------------------------------------------------------
!--   find min z in ztemp                                           --
!---------------------------------------------------------------------
      do j=1,mbox
        zmin=999.
        do i=1,mbox-j+1
          if (ztemp(i).lt.zmin) then
            zmin=ztemp(i)
            kmin=i
          endif
        enddo
!---------------------------------------------------------------------
!--     put zmin into new vectors                                   --
!---------------------------------------------------------------------
        zinc(j)=zmin
        nemtemp(j)=nemprof(kmin)
        temtemp(j)=temprof(kmin)
        nertemp(j)=nerprof(kmin)
        tertemp(j)=terprof(kmin)
!---------------------------------------------------------------------
!--     create new ztemp with remaining data                        --
!---------------------------------------------------------------------
        k=0
        do i=1,mbox-j+1
          if (zmin.ne.ztemp(i)) then
            k=k+1
            ztemp(k)=ztemp(i)
          endif
        enddo
      enddo
!---------------------------------------------------------------------
!--   rewrite new vectors                                           --
!---------------------------------------------------------------------
      do j=1,mbox
        zprof(j)=zinc(j)
        nemprof(j)=nemtemp(j)
        temprof(j)=temtemp(j)
        nerprof(j)=nertemp(j)
        terprof(j)=tertemp(j)
      enddo
      !
      return
      end subroutine tsorder

!**********************************************************************
!!
!!    FITPP does polynomial fitting of p' by using a LAPACK singular
!!      value decomposition routine
!!
!!    @param y : input array being fitted, supposing X is evenly
!!               distributed between [0,1]
!!
!!    @param ny : number of points
!!
!!    @param alpa : return fitting coefficients in its first nalpa
!!                  elements.  Y(x)=a0+a1*x+a2*x^2+...
!!
!!    @param nalpa : number of coefficients
!!
!**********************************************************************
      subroutine fitpp(y,ny,alpa,nalpa)
      use set_kinds, only: dp
      use var_nio
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      parameter (malpa=30)
      dimension x(ny),y(ny),xpsii(nalpa)
      dimension arsp(ny,malpa),work(ny*2),wrsp(malpa),alpa(ny)

      y0=y(1)
      if(abs(y0).le.1.e-5) y0=y(ny)
      do i=1,ny
        x(i)=real(i-1,dp)/real(ny-1,dp)
        y(i)=y(i)/y0
      enddo
!--------------------------------------------------------------
!--  fitting                                                 --
!--------------------------------------------------------------
      arsp=0.0
      if (nalpa.gt.0) then
        do i=1,ny
          call setpp(x(i),arsp(i,1:nalpa))
          alpa(i)=y(i)
        enddo
        nfit=ny
!       write (6,*) nfit,nalpa,(x(i),i=1,nalpa)
!       write (6,*) y(1),y0,y(nfit)
!----------------------------------------------------------------
! ---   LAPACK wrapper routine to do singular value decomposition
!----------------------------------------------------------------
        call sdecm(arsp,ny,nfit,nalpa,alpa,ny,1,wrsp,work,ier)
        if (ier.eq.129) then
          write (nttyo,8000) ier
          stop
        endif
        cond=ier
        toler=1.0e-06*wrsp(1)
        do i=1,nalpa
          t=0.0
          if (wrsp(i).gt.toler) t=alpa(i)/wrsp(i)
          work(i)=t
        enddo
        do i=1,nalpa
          alpa(i)=0.0
          do j=1,nalpa
            alpa(i)=alpa(i)+arsp(i,j)*work(j)
          enddo
        enddo
        alpa0=1.
      endif
!     write (nttyo,3187)
! 3187 format (/,1x,' alphap = ')
!     write(nttyo,*) (alpa(i),i=1,nalpa)

      do i=1,nalpa
        alpa(i)=alpa(i)*y0
      enddo
!      write (nttyo,*) (alpa(i),i=1,nalpa)
      return
 8000 format (/,1x,'Problem in fitting pp')
      end subroutine fitpp

!**********************************************************************
!!
!!    FITFP does polynomial fitting of FF' by using a LAPACK singular
!!      value decomposition routine
!!
!!    @param y : input array being fitted, supposing X is evenly
!!               distributed between [0,1]
!!
!!    @param ny : number of points
!!
!!    @param alpa : return fitting coefficients in its first nalpa
!!                  elements.  Y(x)=a0+a1*x+a2*x^2+...
!!
!!    @param nalpa : number of coefficients
!!
!**********************************************************************
      subroutine fitfp(y,ny,alpa,nalpa)
      use set_kinds, only: dp
      use var_nio
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      parameter (malpa=30)
      dimension x(ny),y(ny),xpsii(nalpa)
      dimension arsp(ny,malpa),work(ny*2),wrsp(malpa),alpa(ny)

      y0=y(1)
      if(abs(y0).le.1.e-5) y0=y(ny)
      do i=1,ny
        x(i)=real(i-1,dp)/real(ny-1,dp)
        y(i)=y(i)/y0
      enddo
!--------------------------------------------------------------
!--  fitting                                                 --
!--------------------------------------------------------------
      arsp=0.0
      if (nalpa.gt.0) then
        do i=1,ny
          call setfp(x(i),arsp(i,1:nalpa))
          alpa(i)=y(i)
        enddo
        nfit=ny
!       write (6,*) nfit,nalpa,(x(i),i=1,nalpa)
!       write (6,*) y(1),y0,y(nfit)
!----------------------------------------------------------------
! ---   LAPACK wrapper routine to do singular value decomposition
!----------------------------------------------------------------
        call sdecm(arsp,ny,nfit,nalpa,alpa,ny,1,wrsp,work,ier)
        if (ier.eq.129) then
          write (nttyo,8000) ier
          stop
        endif
        cond=ier
        toler=1.0e-06*wrsp(1)
        do i=1,nalpa
          t=0.0
          if (wrsp(i).gt.toler) t=alpa(i)/wrsp(i)
          work(i)=t
        enddo
        do i=1,nalpa
          alpa(i)=0.0
          do j=1,nalpa
            alpa(i)=alpa(i)+arsp(i,j)*work(j)
          enddo
        enddo
        alpa0=1.
      endif
!     write (nttyo,3187)
! 3187 format (/,1x,' alphap = ')
!     write(nttyo,*) (alpa(i),i=1,nalpa)

      do i=1,nalpa
        alpa(i)=alpa(i)*y0
      enddo
!      write (nttyo,*) (alpa(i),i=1,nalpa)
      return
 8000 format (/,1x,'Problem in fitting Fp')
      end subroutine fitfp

!**********************************************************************
!!
!>    lenco2 calculates the co2 path lengths.
!!
!!    @param xplt :
!!    @param yplt :
!!    @param nplt :
!!    @param jges : time index
!!
!*********************************************************************
      subroutine lenco2(xplt,yplt,nplt,jges)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: jges,nplt
      real*8, intent(in) :: xplt(npoint),yplt(npoint)
      integer*4 i,k
      real*8 rpath,yr1,yr2,zoutjme,zpath
      real*8 zuper(nco2v),zlower(nco2v),rco2(nco2r),rco2in(nco2r), &
                izout(nco2v)
      real*8, parameter :: zcentr=0.,zlibim=-0.127_dp
!
      rco2=100.
      rco2in=0.
      zuper=100.
      zlower=0.
      izout=0.
      zuperts(jges)=100.
      zlowerts=0.
      zoutjme=zout(jges)*0.01
#ifdef DEBUG_LEVEL1
      write(6,*) ' ZOUT = ',zoutjme
#endif
      do i=1,nplt-1
        do k=1,nco2r
          yr1=chordr(k)-yplt(i)
          yr2=chordr(k)-yplt(i+1)
          if(yr1*yr2.gt.0.0) cycle
          rpath=xplt(i)+yr1*(xplt(i+1)-xplt(i))/(yplt(i+1)-yplt(i))
          if(rpath.ge.rcentr) rco2(k)=rpath
          if(rpath.le.rcentr) rco2in(k)=rpath
        enddo
        yr1=zlibim-yplt(i)
        yr2=zlibim-yplt(i+1)
        if (yr1*yr2.le.0.0) then
          rpath=xplt(i)+yr1*(xplt(i+1)-xplt(i))/(yplt(i+1)-yplt(i))
          if(rpath.ge.rcentr) rlibim(jges)=rpath*100.0
        endif
        do k=1,nco2v
          yr1=chordv(k)-xplt(i)
          yr2=chordv(k)-xplt(i+1)
          if(yr1*yr2.gt.0.0) cycle
          zpath=yplt(i)+yr1*(yplt(i+1)-yplt(i))/(xplt(i+1)-xplt(i))
!          if(zpath.ge.zcentr) zuper(k)=zpath
!          if(zpath.le.zcentr) zlower(k)=zpath
!          if(zpath.ge.zoutjme) zuper(k)=zpath
!          if(zpath.le.zoutjme) zlower(k)=zpath
          izout(k)=izout(k)+1
          if(izout(k).eq.1) zuper(k)=zpath
          if (izout(k).eq.2) then
            zlower(k)=zpath
            if (zuper(k).lt.zpath) then
              zlower(k)=zuper(k)
              zuper(k)=zpath
            endif
          endif
        enddo
        yr1=rmajts-xplt(i)
        yr2=rmajts-xplt(i+1)
        if (yr1*yr2.le.0.0) then
          zpath=yplt(i)+yr1*(yplt(i+1)-yplt(i))/(xplt(i+1)-xplt(i))
          if(zpath.ge.zcentr) zuperts(jges)=zpath*100.
          if(zpath.le.zcentr) zlowerts=zpath*100.
        endif
      enddo
      rco2v(:,jges)=100.0*(zuper-zlower)
      dco2v(jges,:)=denvt(jges,:)/rco2v(:,jges)
#ifdef DEBUG_LEVEL1
      do k=1,nco2v
        write(6,*) ' V, RCO2V, DCO2V, ZU, ZL = ',k, rco2v(k,jges), &
                   dco2v(jges,k), zuper(k), zlower(k)
      enddo
#endif

      rco2r(:,jges)=100.0*(rco2-rco2in)
      dco2r(jges,:)=denrt(jges,:)/rco2r(:,jges)
!
      return
      end subroutine lenco2
