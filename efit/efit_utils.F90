!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          FLUXAV does the flux surface average.                   **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**          f     : function array                                  **
!**          si    : flux array                                      **
!**          x,y,n : n (R,Z) coordinates of contour                  **
!**                                                                  **
!**     RETURN ARGUMENTS:                                            **
!**          fave  : int(dl f/Bp) / sdlobp                           **
!**          sdlbp : int(dl Bp)                                      **
!**          sdlobp: int(dl/Bp)                                      **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          14/10/87..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine fluxav(f,x,y,n,si,rx,msx,ry,msy,fave,ns,sdlobp,sdlbp)
      use set_kinds
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparm.inc'
      use var_cwork3, only:lkx,lky
!vas
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      dimension f(n),x(n),y(n),si(n),pds(6),rx(n),ry(n)
      if (ns.ne.0) then
        call sets2d(si,c,rx,msx,bkx,lkx,ry,msy,bky,lky,wk,ier)
      endif
!------------------------------------------------------------------
!--  flux surface average of f                                   --
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
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**                                                                  **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine splitc(is,rs,zs,cs,rc,zc,wc,hc,ac,ac2,cc)
      use global_constants
      use set_kinds
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension rs(is*is),zs(is*is),cs(is*is)

      if(ac+ac2.eq.0.) go to 100
      if(ac.ne.0.) go to 200
      if(ac2.ne.0.) go to 300
!----------------------------------------------------------------------
!-- rectangle                                                        --
!----------------------------------------------------------------------
  100 continue
      wdelt=wc/is
      hdelt=hc/is
      rstrt=rc-wc/2.+wdelt/2.
      zstrt=zc-hc/2.+hdelt/2.
      zz=zstrt
      ic=0
      c=cc/(is*is)
      do 120 ii=1,is
      rr=rstrt
      do 110 jj=1,is
      ic=ic+1
      zs(ic)=zz
      rs(ic)=rr
      cs(ic)=c
      rr=rr+wdelt
  110 continue
      zz=zz+hdelt
  120 continue
      go to 900
!----------------------------------------------------------------------
!-- ac .ne. 0                                                        --
!----------------------------------------------------------------------
  200 continue
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
      do 220 ii=1,is
      zz=zstrt+(ii-1)*zdelt
      do 210 jj=1,is
      ic=ic+1
      zs(ic)=zz
      rs(ic)=rr
      cs(ic)=c
      zz=zz+hdelt
  210 continue
      rr=rr+wdelt
  220 continue
      go to 900
!----------------------------------------------------------------------
!-- ac2 .ne. 0                                                       --
!----------------------------------------------------------------------
  300 continue
!
  340 continue
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
      do 320 ii=1,is
      rr=rstrt+(ii-1)*rdelt
      do 310 jj=1,is
      ic=ic+1
      zs(ic)=zz
      rs(ic)=rr
      cs(ic)=c
      rr=rr+wdelt
  310 continue
      zz=zz+hdelt
  320 continue
      go to 900
!
  900 continue
      return
      end subroutine splitc
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          this subroutine reorders the z profile data to be       **
!**          in ascending order and sets the ne and te data to       **
!**          correspond to the new order.                            **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          27/11/90..........first created, T. Carlstrom           **
!**          08/07/91..........revised for EFIT, L. Lao              **
!**                                                                  **
!**********************************************************************
      subroutine tsorder(mbox,zprof,nemprof,temprof,nerprof,terprof)
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension ztemp(40),temtemp(40),tertemp(40),zinc(40), &
                  zprof(1), temprof(1), terprof(1)
      real*8 nemtemp(40),nertemp(40),nemprof(1),nerprof(1)
      !---------------------------------------------------------------------
      !-- copy zprof to ztemp (temporary work space)                  --
      !---------------------------------------------------------------------
      do 100 j=1,mbox
        ztemp(j)=zprof(j)
100   continue
      !---------------------------------------------------------------------
      !-- find min z in ztemp                                         --
      !---------------------------------------------------------------------
      do 1000 j=1,mbox
        zmin=999.
        do 800 i=1,mbox-j+1
          if(ztemp(i).lt.zmin) then
            zmin=ztemp(i)
            kmin=i
          end if
800     continue
        !---------------------------------------------------------------------
        !-- put zmin into new vectors                                   --
        !---------------------------------------------------------------------
        zinc(j)=zmin
        nemtemp(j)=nemprof(kmin)
        temtemp(j)=temprof(kmin)
        nertemp(j)=nerprof(kmin)
        tertemp(j)=terprof(kmin)
        !---------------------------------------------------------------------
        !-- create new ztemp with remaining data                        --
        !---------------------------------------------------------------------
        k=0
        do 900 i=1,mbox-j+1
          if(zmin.ne.ztemp(i))then
            k=k+1
            ztemp(k)=ztemp(i)
          end if
900     continue
1000  continue
      !---------------------------------------------------------------------
      !-- rewrite new vectors                                         --
      !---------------------------------------------------------------------
      do 1200 j=1,mbox
        zprof(j)=zinc(j)
        nemprof(j)=nemtemp(j)
        temprof(j)=temtemp(j)
        nerprof(j)=nertemp(j)
        terprof(j)=tertemp(j)
1200  continue
      !
      return
      end subroutine tsorder


