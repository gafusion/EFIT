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
!**********************************************************************
      subroutine fluxav(f,x,y,n,si,rx,msx,ry,msy,fave,ns,sdlobp,sdlbp)
      use set_kinds
      use commonblocks,only: c,wk,copy,bkx,bky
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
      use set_kinds
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


