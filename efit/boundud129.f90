      subroutine bound(psi,nw,nh,nwh,psivl,xmin,xmax,ymin,ymax, &
           zero,x,y,xctr,yctr,ix,limitr,xlim,ylim,xcontr,ycontr, &
           ncontr,xlmin,npoint,rymin,rymax,dpsi,zxmin,zxmax,nerr, &
           ishot,itime,limfag,radold,kbound,tolbndpsi)
!**********************************************************************
!**                                                                  **
!**     main program:  mhd fitting code                              **
!**                                                                  **
!**                                                                  **
!**     subprogram description:                                      **
!**          bound finds the outermost contour of a function         **
!**          given in the array psi. note that the contoured         **
!**          function psi is assumed to be a decreasing function     **
!**          of distance away from the axis.                         **
!**          if ix<0, then BOUND just traces the field line          **
!**          starting at (xctr,yctr) either to closure or to the     **
!**          wall. ix=-1, trace clockwise. -2, counter clockwise.    **
!**                                                                  **
!**     calling arguments:                                           **
!**       psi.............function to be contoured                   **
!**       nw..............row dimension of psi                       **
!**       nh..............column dimension of psi                    **
!**       nwh.............nw x nh                                    **
!**       psivl...........psi value on the outermost contour         **
!**       xmin............minimum r value on contour                 **
!**       xmax............maximum r value                            **
!**       ymin............minimum z value                            **
!**       ymax............maximum z value                            **
!**       zero............weighting array                            **
!**       x...............r grid                                     **
!**       y...............z grid                                     **
!**       xctr............r guess of contour center                  **
!**       yctr............z guess of contour center                  **
!**       ix..............flag                                       **
!**       limitr..........number of limiter points                   **
!**       xlim............r coordinates of limiter points            **
!**       ylim............z coordinates of limiter points            **
!**       xcontr..........output r coordinates of contour            **
!**       ycontr..........output z coordinates of contour            **
!**       ncontr..........number of contour points found             **
!**       xlmin...........minimum or maximum limiter coordinates     **
!**       npoint..........maximum number of contour points           **
!**       dpsi............difference of psi values on contour        **
!**                       and on limiter.  dpsi is used to           **
!**                       distinguish between a limited and a        **
!**                       diverted plasma                            **
!**       rymin...........r at ymin                                  **
!**       rymax...........r at ymax                                  **
!**       zxmin...........z at xmin                                  **
!**       zxmax...........z at xmax                                  **
!**       tolbndpsi (in)..tolerance on psi                           **
!**                                                                  **
!**     references:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     record of modification:                                      **
!**          26/04/83..........first created                         **
!**          01/09/83..........replaced contouring routines          **
!**                            with new dec10 version                **
!**          12/03/84..........modifications to allow function as a  **
!**                            one field line tracer routine.        **
!**                            R Stambaugh.                          **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use set_kinds
      use error_control
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension psi(*),zero(*),x(*),y(*),xcontr(*),ycontr(*)
      dimension dist(5),xlim(*),ylim(*)
      data etolc,etol,nloop/1.e-06_dp,1.e-04_dp,60/
      data nttyo/6/,psitol/1.0e-04_dp/,mecopy/0/,n111/1/
      save dx,dy,area,rmid,mecopy

      save n111

      !! Debug tool: Write out the surface being contoured to
      !! verify it looks reasonable, sometimes it's not.
      !open(unit=93,file='junk-psi.txt',status='unknown')
      !do i = 1,nw
      !  do j = 1,nh
      !    k = (i-1)*nh+j
      !    write(93,*) x(i),y(j),psi(k)
      !  end do
      !end do
      !close(unit=93)

!----------------------------------------------------------------------
!--           nerr=10000, negative plasma current                    --
!----------------------------------------------------------------------
      nosign=0
      if (nerr .eq. 10000) then
        nosign=1
        do i=1,nwh
          psi(i)=-psi(i)
        end do
      endif
!------------------------------------------------------------------------
!--          BBrown's version of BOUND                                 --
!------------------------------------------------------------------------
      if (ix.gt.0 .and. kbound.ne.0) then
        !call old_new    (psi,nw,nh,nwh,psivl,xmin,xmax,ymin,ymax, &
        !     zero,x,y,xctr,yctr,ix,limitr,xlim,ylim,xcontr,ycontr, &
        !     ncontr,xlmin,npoint,rymin,rymax,dpsi,zxmin,zxmax,nerr, &
        !     ishot,itime,limfag,radold,kbound)
        go to 2000
      endif

      nerr=0
      psib0=-1.e+10_dp
      rad=radold
      rin=xctr
      rout=xlmin
!---------------------------------------------------------------------
!--  field line tracing ix < 0                                      --
!---------------------------------------------------------------------
      if (ix.lt.0) rad=xctr
      if (mecopy.le.0) then
        dx=x(2)-x(1)
        dy=y(2)-y(1)
        area=dx*dy
        rmid=1.02_dp*(x(1)+x(nw))/2.0
        mecopy=1
      end if

      !----------------------------------------------------------------------
      !--   find starting value of psi                                     --
      !----------------------------------------------------------------------
      do loop = 1,nloop
        i=1+(rad-x(1))/dx
        if(rad-x(i).lt.0.0)i=i-1
        j=1+(yctr-y(1))/(dy-0.000001_dp)
        jjj=j

        if ((ix.eq.-2).or.(ix.lt.-2.and.rad.gt.rmid)) then
          j=j-1
          jjj=j+1
        endif

        if ((ix.gt.0).and.(rad.gt.rmid)) then
          j=j-1
          jjj=j+1
        endif
        kstrt=(i-1)*nh+j
        kk=kstrt
        kold=kk

        if (ix.ge.0.or.ix.lt.-2) then
          dpsi=1.e10_dp
          xmin=1.e10_dp
          xmax=-1.e10_dp
          ymin=1.e10_dp
          ymax=-1.e10_dp
        end if

        xt=rad
        yt=y(jjj)
        ncontr=1
        xcontr(ncontr)=xt
        ycontr(ncontr)=yt
        a3=(xt-x(i))*dy
        a4=area-a3
        psivl=(psi(kk+nh+1)*a3+psi(kk+1)*a4)/area
        if (yt.eq.y(j)) psivl=psi(kk)+(psi(kk+nh)-psi(kk))*(xt-x(i))/dx

        do while (.true.) ! contr
          f1=psi(kk)
          f2=psi(kk+nh)
          f3=psi(kk+nh+1)
          f4=psi(kk+1)
          x1=x(i)
          x2=x(i+1)
          y1=y(j)
          y2=y(j+1)
          if(ncontr.eq.1) go to 100
          !----------------------------------------------------------------------
          !--   check for proximity to corner                                  --
          !----------------------------------------------------------------------
          dist(1)=(xt-x1)**2+(yt-y1)**2
          dist(2)=(xt-x2)**2+(yt-y1)**2
          dist(3)=(xt-x2)**2+(yt-y2)**2
          dist(4)=(xt-x1)**2+(yt-y2)**2
          dist(5)=min(dist(1),dist(2),dist(3),dist(4))
          if (dist(5).gt.etolc) go to 100
          do l=1,4
            kj=l
            if (dist(l).eq.dist(5)) exit
          end do
          !----------------------------------------------------------------------
          !--   kj points to appropriate corner                                --
          !----------------------------------------------------------------------
          call chkcrn(psi,nwh,psivl,kold,knew,kj,kk,nh,i,i1)
          kk=knew
          i=i1
          j=kk-(i-1)*nh
          f1=psi(kk)
          f2=psi(kk+nh)
          f3=psi(kk+nh+1)
          f4=psi(kk+1)
          x1=x(i)
          x2=x(i+1)
          y1=y(j)
          y2=y(j+1)
          !----------------------------------------------------------------------
          !--   check for limiter in cell                                      --
          !----------------------------------------------------------------------
100       zsum=zero(kk)+zero(kk+1)+zero(kk+nh)+zero(kk+nh+1)
          if(zsum.eq.0.0) exit ! contr

          if (abs(zsum-4.0).ge.1.e-03_dp) then
            !----------------------------------------------------------------------
            !--   from one to three corners of cell are inside limiter.  get max --
            !--   psi on line segment of limiter in cell and compare this max    --
            !--   with current value of psilim (or psisep)                       --
            !--   note: do loop index assumes point 'limitr+1' is the same as    --
            !--   point 'limitr'                                                 --
            !----------------------------------------------------------------------
            psilx=-1.e10_dp
            do k=1,limitr-1
              xc1=xlim(k)
              yc1=ylim(k)
              xc2=xlim(k+1)
              yc2=ylim(k+1)
              ik1=1+(xlim(k)-x(1))/dx
              jk1=1+(ylim(k)-y(1))/dy
              ik2=1+(xlim(k+1)-x(1))/dx
              jk2=1+(ylim(k+1)-y(1))/dy
              kij1=(ik1-1)*nh+jk1
              kij2=(ik2-1)*nh+jk2
              !----------------------------------------------------------------------
              !--  at least one limiter point is not in cell.  subroutine cellb    --
              !--  returns intersections of cell boundaries and line defined by    --
              !--  points k and k+1 or one cell boundary and one interior point.   --
              !--  ifail=1 if points k and k+1 do not intersect the current cell   --
              !--  of interest.                                                    --
              !----------------------------------------------------------------------
              if ((kij1.ne.kk).or.(kij2.ne.kk)) then
                ifail=0
                call cellb(xc1,yc1,xc2,yc2,x1,y1,x2,y2,ifail)
                if (ifail.eq.1) cycle ! line segment does not intersect cell
              end if
              ! psilm is largest psi value along line segment betw pts
              call maxpsi(xc1,yc1,xc2,yc2,x1,y1,x2,y2,f1,f2,f3,f4,area,psilm,xtry1,ytry1,nerr)
              if (nerr.gt.0) go to 2000
              psilx=max(psilm,psilx)
              if (psilx.le.psilm) then
                xtry=xtry1
                ytry=ytry1
              end if
            end do ! limitr

            if (psilx.eq.-1.e10_dp) then
              nerr=3
              call errctrl_msg('bound','Limiter points do not intersect cell')
              go to 2000
            end if

            dpsi=min(dpsi,abs(psivl-psilx))
            if (psilx-psivl.ge.tolbndpsi) then
              call zlim(zerol,n111,n111,limitr,xlim,ylim,xt,yt,limfag)
              if (zerol.le.0.01_dp) then
                exit ! contr
              end if
            end if
          end if

          call extrap(f1,f2,f3,f4,x1,y1,x2,y2,xt,yt,xt1,yt1,xt2,yt2, &
            psivl,area,dx,dy)
          !----------------------------------------------------------------------
          !--   decide which intersection (xt1,yt1) or (xt2,yt2) is required   --
          !----------------------------------------------------------------------
          dist1=(yt1-yt)**2+(xt1-xt)**2
          dist2=(yt2-yt)**2+(xt2-xt)**2
          if (dist1.lt.dist2) then
            yt=yt2
            xt=xt2
          else
            yt=yt1
            xt=xt1
          end if
          ncontr=ncontr+1
          if (ncontr.gt.npoint) then
            nerr=3
            call errctrl_msg('bound','Number of contour points greater than max allowed')
            go to 2000
          end if
          xcontr(ncontr)=xt
          ycontr(ncontr)=yt

          ! Debug tool: Write out the contour coordinates for each loop (iteration)
          !write(*,*) loop,ncontr,xcontr(ncontr),ycontr(ncontr)

          !----------------------------------------------------------------------
          !--   find next cell                                                 --
          !----------------------------------------------------------------------
          if (xt.eq.x2) i=i+1
          if (xt.eq.x1) i=i-1
          if (yt.eq.y2) j=j+1
          if (yt.eq.y1) j=j-1

          if (ix.ge.0.or.ix.lt.-2) then
            if(yt.lt.ymin)rymin=xt
            if(yt.gt.ymax)rymax=xt
            if(xt.lt.xmin)zxmin=yt
            if(xt.gt.xmax)zxmax=yt
            xmin=min(xmin,xt)
            xmax=max(xmax,xt)
            ymin=min(yt,ymin)
            ymax=max(yt,ymax)
          end if

          kold=kk
          !----------------------------------------------------------------------
          !     find new cell index
          !----------------------------------------------------------------------
          kk=(i-1)*nh+j
          if (kk.eq.kstrt) go to 1040
          dis2p=sqrt((xcontr(1)-xt)**2+(ycontr(1)-yt)**2)
          if((dis2p.lt.0.1_dp*dx).and.(ncontr.gt.5))go to 1040
        end do ! contr

        !----------------------------------------------------------------------
        !--  psi on boundary smaller than psi on limiter, decrease rad and   --
        !--  try again.                                                      --
        !----------------------------------------------------------------------
        psib0=psivl
        !
        if(ix.lt.0) go to 2000 ! ix, -1=trace clockwise, -2=counter clockwise
        !
        if(loop.ge.nloop) exit ! loop
        rout=rad
        rad=(rin+rout)*0.5_dp
        cycle ! loop

        !----------------------------------------------------------------------
        !--   check for convergence of boundary                              --
        !----------------------------------------------------------------------
1040    err=abs((psivl-psib0)/psivl)
        if(ix.lt.0) then
          if (ix.lt.-2) dpsi=1.e-06_dp
          go to 2000
        end if

        if (err.le.etol) exit ! loop
        if (loop.ge.nloop) exit ! loop
        !----------------------------------------------------------------------
        !--   new rad,psi and try again                                      --
        !----------------------------------------------------------------------
        psib0=psivl
        call zlim(zerol,n111,n111,limitr,xlim,ylim,rad,yctr,limfag)
        if (zerol.le.0.01_dp) then
          rout=rad
        else
          rin=rad
        endif
        rad=(rin+rout)*0.5_dp
      end do ! loop

      radold=rad
      psib0=psivl
      if ((abs(ycontr(1)-ycontr(ncontr)).gt.0.5_dp*dy) .or. &
          (abs(xcontr(1)-xcontr(ncontr)).gt.0.5_dp*dx)) then
         nerr=3
         call errctrl_msg('bound','First and last contour points are too far apart')
      end if

 2000 continue
      if (ncontr.lt.3) then
        nerr=3
        call errctrl_msg('bound','Less than 3 contour points found')
      end if
      if (nosign.eq.1) then
        do i=1,nwh
          psi(i)=-psi(i)
        end do
        psivl=-psivl
      endif
      return
      end

      subroutine cellb(xc1,yc1,xc2,yc2,x1,y1,x2,y2,ifail)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          cellb redefines (xc1,yc1) and/or (xc2,yc2) so that      **
!**          they are intersections of cell boundaries unless        **
!**          one of the points is an interior point in which         **
!**          case it is not disturbed.                               **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      if((xc1.lt.x1).and.(xc2.lt.x1))go to 1000
      if((xc1.gt.x2).and.(xc2.gt.x2))go to 1000
      if((yc1.lt.y1).and.(yc2.lt.y1))go to 1000
      if((yc1.gt.y2).and.(yc2.gt.y2))go to 1000
      dx=xc1-xc2
      dy=yc1-yc2
      if((dx.eq.0.0).and.(dy.eq.0.0))ifail=1
      if(ifail.eq.1)return
      if(dx.eq.0.0)go to 200
      if(dy.eq.0.0)go to 500
!----------------------------------------------------------------------
!--   line is inclined. get equation y=alpha*x+beta.                 --
!----------------------------------------------------------------------
      alpha=dy/dx
      beta=yc1-alpha*xc1
      y1b=alpha*x1+beta
      if(y1b.lt.y1)go to 20
      if(y1b.gt.y2)go to 30
      xs1=x1
      ys1=y1b
      go to 40
   20 x1b=(y1-beta)/alpha
      if ((x1b.lt.x1).or.(x1b.gt.x2))go to 1000
      xs1=x1b
      ys1=y1
      go to 40
   30 x1b=(y2-beta)/alpha
      if((x1b.lt.x1).or.(x1b.gt.x2))go to 1000
      xs1=x1b
      ys1=y2
   40 y2b=alpha*x2+beta
      if(y2b.lt.y1)go to 50
      if(y2b.gt.y2)go to 60
      xs2=x2
      ys2=y2b
      go to 70
   50 x2b=(y1-beta)/alpha
      if((x2b.gt.x2).or.(x2b.lt.x1))go to 1000
      xs2=x2b
      ys2=y1
      go to 70
   60 x2b=(y2-beta)/alpha
      if((x2b.gt.x2).or.(x2b.lt.x1))go to 1000
      xs2=x2b
      ys2=y2
!----------------------------------------------------------------------
!--   at this point we have intersection (xs1,ys1) and (xs2,ys2)     --
!--   check for interior point                                       --
!----------------------------------------------------------------------
   70 if((x1.lt.xc1).and.(xc1.lt.x2))go to 80
      if((x1.lt.xc2).and.(xc2.lt.x2))go to 90
      go to 100
   80 if((y1.lt.yc1).and.(yc1.lt.y2))go to 160
      go to 100
   90 if((y1.lt.yc2).and.(yc2.lt.y2))go to 170
  100 xc1=xs1
      yc1=ys1
      xc2=xs2
      yc2=ys2
      return
!----------------------------------------------------------------------
!--   point (xc1,yc1) is interior point                              --
!----------------------------------------------------------------------
  160 if(yc2.gt.yc1)go to 165
      xc2=xs1
      yc2=ys1
      if(ys1.gt.ys2)yc2=ys2
      if(ys1.gt.ys2)xc2=xs2
      return
  165 xc2=xs2
      yc2=ys2
      if(ys1.gt.ys2)yc2=ys1
      if(ys1.gt.ys2)xc2=xs1
      return
!----------------------------------------------------------------------
!--   point (xc2,yc2) is interior point                              --
!----------------------------------------------------------------------
  170 if(yc1.gt.yc2)go to 190
      xc1=xs1
      yc1=ys1
      if(ys1.gt.ys2)yc1=ys2
      if(ys1.gt.ys2)xc1=xs1
      return
  190 xc1=xs2
      yc1=ys2
      if(ys1.gt.ys2)yc1=ys1
      if(ys1.gt.ys2)xc1=xs1
      return
!----------------------------------------------------------------------
!--   check if intersection exists for vertical line                 --
!----------------------------------------------------------------------
  200 if((xc1.lt.x1).or.(xc1.gt.x2))go to 1000
!----------------------------------------------------------------------
!--   is there an interior point ?                                   --
!----------------------------------------------------------------------
      if((y1.le.yc1).and.(yc1.le.y2))go to 300
      if((y1.le.yc2).and.(yc2.le.y2))go to 400
!----------------------------------------------------------------------
!--   no interior points                                             --
!----------------------------------------------------------------------
      yc1=y1
      yc2=y2
      return
!----------------------------------------------------------------------
!--   point (xc1,yc1) is interior point                              --
!----------------------------------------------------------------------
  300 if(yc2.gt.yc1)yc2=y2
      if(yc2.lt.yc1)yc2=y1
      return
!----------------------------------------------------------------------
!--   point  (xc2,yc2) is interior point                             --
!----------------------------------------------------------------------
  400 if(yc2.gt.yc1)yc1=y1
      if(yc2.lt.yc1)yc1=y2
      return
!----------------------------------------------------------------------
!--   check if intersection exists for horizontal line               --
!----------------------------------------------------------------------
  500 if((yc1.lt.y1).or.(yc1.gt.y2))go to 1000
!----------------------------------------------------------------------
!--   is there an interior point ?                                   --
!----------------------------------------------------------------------
      if((x1.lt.xc1).and.(xc1.le.x2))go to 550
      if((x1.le.xc2).and.(xc2.le.x2))go to 600
!----------------------------------------------------------------------
!--   no interior points                                             --
!----------------------------------------------------------------------
      xc1=x1
      xc2=x2
      return
!----------------------------------------------------------------------
!--   point (xc1,yc1) is interior point                              --
!----------------------------------------------------------------------
  550 if(xc1.lt.xc2)xc2=x2
      if(xc1.gt.xc2)xc2=x1
      return
!----------------------------------------------------------------------
!--   point (xc2,yc2) is interior point                              --
!----------------------------------------------------------------------
  600 if(xc2.gt.xc1)xc1=x1
      if(xc2.lt.xc1)xc1=x2
      return
 1000 ifail=1
      return
      end
      subroutine chkcrn(psi,nwh,psivl,kold,knew,icrnr,kk,nh,i,i1)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          chkcrn decides which cell is next when an error in      **
!**          the cell step is likely.                                **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          07/09/83..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension psi(nwh)
      knew=0
      i1=0
!----------------------------------------------------------------------
!--   cell #1                                                        --
!----------------------------------------------------------------------
      kn=kk
      if(kn.eq.kold)go to 10
      call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
      if(iflag.eq.1)go to 10
      i1=i
      return
   10 continue
      select case (icrnr)
      case (1)
        go to 100
      case (2)
        go to 200
      case (3)
        go to 300
      case (4)
        go to 400
      end select
!----------------------------------------------------------------------
!--   corner #1                                                      --
!--   cell #2                                                        --
!----------------------------------------------------------------------
  100 kn=kk-nh
      if(kn.eq.kold)go to 110
      call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
      if(iflag.eq.1)go to 110
      i1=i-1
      return
!----------------------------------------------------------------------
!--   cell #3                                                        --
!----------------------------------------------------------------------
  110 kn=kk-nh-1
      if(kold.eq.kn)go to 120
      call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
      if(iflag.eq.1) go to 120
      i1=i-1
      return
!----------------------------------------------------------------------
!--   cell #4                                                        --
!----------------------------------------------------------------------
  120 kn=kk-1
      call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
      i1=i-1
      return
!----------------------------------------------------------------------
!--   corner #2                                                      --
!--   cell #2                                                        --
!----------------------------------------------------------------------
  200 kn=kk-1
      if(kn.eq.kold)go to 210
      call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
      if(iflag.eq.1)go to 210
      i1=i
      return
!----------------------------------------------------------------------
!--   cell #3                                                        --
!----------------------------------------------------------------------
  210 kn=kk+nh-1
      if(kn.eq.kold)go to 220
      call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
      if(iflag.eq.1)go to 220
      i1=i+1
      return
!----------------------------------------------------------------------
!--   cell #4                                                        --
!----------------------------------------------------------------------
  220 kn=kk+nh
      call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
      i1=i+1
      return
!----------------------------------------------------------------------
!--   corner #3                                                      --
!--   cell #2                                                        --
!----------------------------------------------------------------------
  300 kn=kk+nh
      if(kn.eq.kold)go to 310
      call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
      if(iflag.eq.1)go to 310
      i1=i+1
      return
!----------------------------------------------------------------------
!--   cell #3                                                        --
!----------------------------------------------------------------------
  310 kn=kk+nh+1
      if(kn.eq.kold)go to 320
      call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
      if(iflag.eq.1)go to 320
      i1=i+1
      return
!----------------------------------------------------------------------
!--   cell #4                                                        --
!----------------------------------------------------------------------
  320 kn=kk+1
      call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
      i1=i
      return
!----------------------------------------------------------------------
!--   corner #4                                                      --
!--   cell #2                                                        --
!----------------------------------------------------------------------
  400 kn=kk+1
      if(kn.eq.kold)go to 410
      call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
      if(iflag.eq.1)go to 410
      i1=i
      return
!----------------------------------------------------------------------
!--   cell #3                                                        --
!----------------------------------------------------------------------
  410 kn=kk-nh+1
      if(kn.eq.kold)go to 420
      call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
      if(iflag.eq.1)go to 420
      i1=i-1
      return
!----------------------------------------------------------------------
!--   cell #4                                                        --
!----------------------------------------------------------------------
  420 kn=kk-nh
      call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
      i1=i-1
      return
      end

      subroutine cntour(xaxd,yaxd,psivl,xemin,xemax,yemin,yemax, &
      yxmin,yxmax,xymin,xymax,dang,arcl,bperr,dx,dy,xmin,xmax,ymin,ymax, &
      iauto,iautoc,xc,yc,ipts,x,nw,y,nh,cspln,n2cspln,nh2,itty,iptsm, &
      negcur,bkx,lkx,bky,lky,kerror)
!------------------------------------------------------------------------------
!---                                                                         --
!---given (xaxd,yaxd) and a psi value, psivl, generate a contour of          --
!---ordered points,(xc(i),yc(i)),i=1,ipts,   which has (xaxd,yaxd)           --
!---as an interior point. the contour must fully encircle (xaxd,yaxd).       --
!---the search is limited to a rectangle specified by xmin,xmax,ymin,ymax.   --
!---dx and dy determine the basic increment in x and y for the coarse        --
!---grid search.  they should be picked so that over the distance dx or dy   --
!---psi is single valued.a sliding interval search rather than a binary      --
!---search is used for the coarse grid so that non-monotonic psi can be      --
!---handled.dx and dy are in meters. dx=dy=0.01 is suggested.                --
!---dang is angular step for the rays wich emmanate from (xaxd,yaxd)         --
!---in degrees.  dang should be set according to shape of equilibrium.       --
!---if dang is too small many unnecessary points will be generated           --
!---in (xc,yc). if dang is too large for bperr(see below) the routine        --
!---wastes computation time in cutting dang down to size.                    --
!---dang=10 deg. near psilim and dang=30 deg. near psimax is suggested.      --
!---for highly elongated plasmas control with dang is difficult to set for all-
!---contours. therefore arcl may be used to set dang internally. arcl is the  -
!---arc length in meters taken from the current point to get to the next point-
!---set arcl to a large number (i.e. 10 meters) to overide this option and guse
!---dang only. the angle increment will be the minimum of (arcl/rad,dang).   --
!---bperr is relative change in poloidal b field between (xc(i),yc(i))       --
!---and (xc(i+1),yc(i+1)). bperr=0.03 is suggested.  note: if dang yields    --
!---bperr<bperr requested then dang is used. otherwise dang is succesively   --
!---reduced by 0.5 until this condition is meet.                             --
!---xemin,xemax,yemin,yemax are min and max values of contour                --
!---yxmin,yxmax,are y values at x=xemin and x=xemax respectively.            --
!---xymin,xymax are x values at y=yemin and y=yemax respectively             --
!---iauto is used as an "error recovery switch". that is,if a flux surface   --
!---is found to pass outside the search box an error exit is taken by way    --
!---of label 1000 if iauto=0. if iauto=1 the input value of psi,psivl,       --
!---is increased toward the value of psi at (xaxd,yaxd) until a closed       --
!---surface inside the search box is found. obviously this option makes      --
!---sense only for the limiter flux surface. it could be used to find        --
!---the plasma boundary if an appopriate search box were defined.            --
!---however its primary use here is to slightly modify psilim so that        --
!---a closed limiter flux surface is found. the option is necessary          --
!---due to the fact some codes generate eqdsks using linear interpolation,   --
!---which differs from the bicubic interpolation used here.                  --
!---the user is informed that psivl was changed by returning iautoc=1.       --
!---if no change occured iautoc=0.                                           --
!------------------------------------------------------------------------------
      use global_constants
      use set_kinds
      use error_control
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension pds(6),xc(*),yc(*)
      dimension cspln(kubicx,lubicx,kubicy,lubicy)
      real*8 piov2,piov4,fpiov4,spiov4,tpiov4,tpiov2
      data n111/1/,n333/3/
      save n111,n333
      integer, intent(inout) :: kerror

      kerror = 0
      ier = 0
      piov2 = pi/2.0
      piov4 = pi/4.0
      fpiov4 = 1.25_dp*pi
      spiov4 = 1.75_dp*pi
      tpiov4 = 0.75_dp*pi
      tpiov2 = 1.5_dp*pi

      if (negcur.gt.0) psivl=-psivl
      iautoc=0
 1040 xemin=1.e10_dp
      xemax=-1.e10_dp
      yemin=1.e10_dp
      yemax=-1.e10_dp
      xymin=0.0
      xymax=0.0
      yxmin=0.0
      yxmax=0.0
      ipts=0
      thet=0.0
      dthet0=twopi*dang/360.
      dthet=0.0
      serrt=3.5e-06_dp
      derrt=0.5e-07_dp
!-----------------------------------------------------------------------------
!---serrt is absolute error convergence criteria for newtons method below.  --
!---get psi at (xaxd,yaxd)                                                  --
!-----------------------------------------------------------------------------
      call seva2d(bkx,lkx,bky,lky,cspln,xaxd,yaxd,pds,ier,n111)
      if (negcur.eq.0) then
        psiaxd=pds(1)
      else
        psiaxd=-pds(1)
      endif
      if (psiaxd.lt.psivl) then
        kerror = 1
        call errctrl_msg('cntour','psiaxd.lt.psivl')
        return
      endif
!----------------------------------------------------------------------
!---loop over theta from 0 to twopi                                  --
!----------------------------------------------------------------------
   10 thet=thet+dthet
      thet=min(thet,twopi)
      if(thet.eq.twopi)go to 200
!----------------------------------------------------------------------
!---get equation of ray emanating from (xaxd,yaxd)                   --
!----------------------------------------------------------------------
      if((piov4.le.thet).and.(thet.le.tpiov4))go to 20
      if((fpiov4.le.thet).and.(thet.le.spiov4))go to 20
!----------------------------------------------------------------------
!---y as a function of x            y=a*x+bincp                      --
!----------------------------------------------------------------------
      isgn=-1
      if((thet.lt.piov4).or.(thet.gt.spiov4))isgn=1
      a=tan(thet)
      iflg=0
      bincp=yaxd-a*xaxd
      go to 30
!----------------------------------------------------------------------
!---x as a function of y            x=a*y+bincp                      --
!----------------------------------------------------------------------
   20 isgn=1
      if(thet.gt.pi)isgn=-1
      if(isgn.eq.-1)go to 22
      thet1=piov2-thet
      if(thet.gt.piov2)thet1=twopi-abs(thet1)
      go to 25
   22 thet1=tpiov2-thet
      if(thet.gt.tpiov2)thet1=pi-abs(thet1)
   25 a=tan(thet1)
      iflg=1
      bincp=xaxd-a*yaxd
   30 continue
!-----------------------------------------------------------------------
!---now have y=a*x+bincp    (iflg=0)   or                             --
!---x=a*y+bincp             (iflg=1)                                  --
!-----------------------------------------------------------------------
      x1=xaxd
      y1=yaxd
      cost=cos(thet)
      sint=sin(thet)
      psi1=psiaxd
!------------------------------------------------------------------------
!---sliding interval search. max width of interval ~1.41*(dx or dy)    --
!------------------------------------------------------------------------
   40 if(iflg.eq.1)go to 50
!---search in x
      x2=x1+isgn*dx
      y2=a*x2+bincp
      go to 60
!---search in y
   50 y2=y1+isgn*dy
      x2=a*y2+bincp
   60 if((x2.lt.xmin).or.(x2.gt.xmax))go to 1000
      if((y2.lt.ymin).or.(y2.gt.ymax))go to 1000
      call seva2d(bkx,lkx,bky,lky,cspln,x2,y2,pds,ier,n111)
      if (negcur.eq.0) then
        psi2=pds(1)
      else
        psi2=-pds(1)
      endif
      dpsi=(psivl-psi1)*(psivl-psi2)
      if(dpsi.le.0.0)go to 70
      x1=x2
      y1=y2
      psi1=psi2
      go to 40
!--------------------------------------------------------------------------
!---now have psivl between psi1 and psi2,converge using newton-raphson   --
!--------------------------------------------------------------------------
   70 newti=0
      if(iflg.eq.1)go to 75
      xn=x1+isgn*dx*0.5_dp
      yn=a*xn+bincp
      go to 80
   75 yn=y1+isgn*dy*0.5_dp
      xn=a*yn+bincp
   80 call seva2d(bkx,lkx,bky,lky,cspln,xn,yn,pds,ier,n333)
      if (negcur.eq.0) then
       dpsids=pds(2)*cost+pds(3)*sint
       dpsi=pds(1)-psivl
      else
       dpsids=-(pds(2)*cost+pds(3)*sint)
       dpsi=-pds(1)-psivl
      endif
      serr=-dpsi/dpsids
      if(abs(serr).lt.serrt)go to 90
      if(abs(dpsi).lt.derrt)go to 90
      delx=serr*cost
      dely=serr*sint
      xn=xn+delx
      yn=yn+dely
      newti=newti+1
      if (newti.ge.20) then
        kerror = 1
        call errctrl_msg('cntour','newti.ge.20')
        return
      endif
      go to 80
!-----------------------------------------------------------------------------
!---end of newton iteration                                                 --
!---check for sufficient accuracy in point spacing as determined by thet    --
!---accuracy test is based on a relative error in poloidal b field of bperr --
!-----------------------------------------------------------------------------
   90 bp2=(1./xn)*sqrt(pds(2)**2+pds(3)**2)
      if(thet.eq.0.0)go to 100
      if(abs((bp2-bp1)/max(bp2,bp1)).lt.bperr)go to 100
!---avoid accumulation at x points
      ihalf=ihalf+1
      if(ihalf.gt.4)go to 100
!---spacing too large for grad psi. decrease theta and try again
      thet=thet-dthet
      dthet=dthet*0.5_dp
      go to 10
  100 bp1=bp2
      ipts=ipts+1
      if (ipts.gt.iptsm-1) then
        kerror = 1
        call errctrl_msg('cntour','ipts.gt.iptsm-1')
        return
      endif
      xc(ipts)=xn
      yc(ipts)=yn
      xemin=min(xemin,xn)
      xemax=max(xemax,xn)
      if(xemax.eq.xn)yxmax=yn
      if(xemin.eq.xn)yxmin=yn
      yemin=min(yemin,yn)
      yemax=max(yemax,yn)
      if(yemax.eq.yn)xymax=xn
      if(yemin.eq.yn)xymin=xn
!---limit angle increment to approximately arcl meters per step
      rad1=sqrt((xn-xaxd)**2+(yn-yaxd)**2)
      dthet=min(dthet0,arcl/rad1)
      ihalf=0
      go to 10
!--------------------------------------------------------------------
!---close contour                                                  --
!--------------------------------------------------------------------
  200 ipts=ipts+1
      xc(ipts)=xc(1)
      yc(ipts)=yc(1)
      return
!---------------------------------------------------------
!---errors                                              --
!---------------------------------------------------------
 1000 continue
      if (iauto.ne.1) then
        kerror = 1
        call errctrl_msg('cntour','flux surface is outside search area, iauto.ne.1')
        return
      end if

      psivl0=psivl
      dapsi=psiaxd-psivl0
      psivl=psivl0+dapsi*0.0005_dp
      iautoc=1
      write (itty,1020) psivl0,psivl
 1020 format(2x,'boundary search, will change psilim from', &
             /,e16.8,'  to  ',e16.8,'  and try again')
      go to 1040
      end
      subroutine extrap(f1,f2,f3,f4,x1,y1,x2,y2,xt,yt,xt1,yt1, &
                        xt2,yt2,psivl,area,dx,dy)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          extrap extrapolates across a (x,y) cell.                **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**          07/09/83..........replaced                              **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension xp(2),yp(2)
      ip=0
      fmin1=min(f1,f2)
      fmax1=max(f1,f2)
      fmin2=min(f2,f3)
      fmax2=max(f2,f3)
      fmin3=min(f3,f4)
      fmax3=max(f3,f4)
      fmin4=min(f1,f4)
      fmax4=max(f1,f4)
      a=y2*(f2-f1)+y1*(f4-f3)
      b=x1*(f2-f3)+x2*(f4-f1)
      c=f1-f2+f3-f4
!----------------------------------------------------------------------
!--   vertical asymptote is x=-b/c                                   --
!--   horizontal asymptote is y=-a/c                                 --
!--   there can only be vertical and horizontal asymptotes           --
!----------------------------------------------------------------------
      if(c.eq.0.0)go to 950
      xasym=-b/c
      yasym=-a/c
      if((x1.le.xasym).and.(xasym.le.x2))go to 840
      if((yasym.lt.y1).or.(yasym.gt.y2))go to 950
!----------------------------------------------------------------------
!--   there is a horizontal asymptote                                --
!--   find psi value on asymptote at x1 and x2                       --
!----------------------------------------------------------------------
      psias4=((f4-f1)*(yasym-y1)/dy)+f1
      psias2=((f3-f2)*(yasym-y1)/dy)+f2
      if(yt.lt.yasym)go to 810
!----------------------------------------------------------------------
!--   upper section of cell                                          --
!----------------------------------------------------------------------
      fpmin4=min(psias4,f4)
      fpmax4=max(f4,psias4)
      if((psivl.lt.fpmin4).or.(psivl.gt.fpmax4))go to 790
      ip=ip+1
      xp(ip)=x1
      yp(ip)=((psivl-f1)/(f4-f1))*dy+y1
      if(psivl.eq.f4)yp(ip)=y2
  790 fpmin2=min(f3,psias2)
      fpmax2=max(f3,psias2)
      if((psivl.lt.fpmin2).or.(psivl.gt.fpmax2))go to 800
      ip=ip+1
      xp(ip)=x2
      yp(ip)=((psivl-f2)/(f3-f2))*dy+y1
      if(psivl.eq.f3)yp(ip)=y2
      if(ip.eq.2)go to 990
  800 ip=ip+1
      yp(ip)=y2
      xp(ip)=((psivl-f4)/(f3-f4))*dx+x1
      go to 990
!----------------------------------------------------------------------
!--   lower section of cell                                          --
!----------------------------------------------------------------------
  810 fpmin4=min(f1,psias4)
      fpmax4=max(f1,psias4)
      if((psivl.lt.fpmin4).or.(psivl.gt.fpmax4))go to 820
      ip=ip+1
      xp(ip)=x1
      yp(ip)=((psivl-f1)/(f4-f1))*dy+y1
  820 fpmin2=min(f2,psias2)
      fpmax2=max(f2,psias2)
      if((psivl.lt.fpmin2).or.(psivl.gt.fpmax2))go to 830
      ip=ip+1
      xp(ip)=x2
      yp(ip)=((psivl-f2)/(f3-f2))*dy+y1
      if(ip.eq.2)go to 990
  830 ip=ip+1
      yp(ip)=y1
      xp(ip)=((psivl-f1)/(f2-f1))*dx+x1
      go to 990
  840 if((y1.lt.yasym).and.(yasym.lt.y2))go to 900
!----------------------------------------------------------------------
!--   vertical asymptote                                             --
!     find psi value on asymptote at y1 and y2
!----------------------------------------------------------------------
      psias1=((f2-f1)*(xasym-x1)/dx)+f1
      psias3=((f3-f4)*(xasym-x1)/dx)+f4
      if(xt.lt.xasym)go to 870
!----------------------------------------------------------------------
!--   right side of cell                                             --
!----------------------------------------------------------------------
      fpmin1=min(psias1,f2)
      fpmax1=max(psias1,f2)
      if((psivl.lt.fpmin1).or.(psivl.gt.fpmax1))go to 850
      ip=ip+1
      yp(ip)=y1
      xp(ip)=((psivl-f1)/(f2-f1))*dx+x1
      if(psivl.eq.f2)xp(ip)=x2
  850 fpmin3=min(f3,psias3)
      fpmax3=max(f3,psias3)
      if((psivl.lt.fpmin3).or.(psivl.gt.fpmax3))go to 860
      ip=ip+1
      yp(ip)=y2
      xp(ip)=((psivl-f4)/(f3-f4))*dx+x1
      if(psivl.eq.f3)xp(ip)=x2
      if(ip.eq.2)go to 990
  860 ip=ip+1
      xp(ip)=x2
      yp(ip)=((psivl-f2)/(f3-f2))*dy+y1
      go to 990
!----------------------------------------------------------------------
!--   left side of cell                                              --
!----------------------------------------------------------------------
  870 fpmin1=min(f1,psias1)
      fpmax1=max(f1,psias1)
      if((psivl.lt.fpmin1).or.(psivl.gt.fpmax1))go to 880
      ip=ip+1
      yp(ip)=y1
      xp(ip)=((psivl-f1)/(f2-f1))*dx+x1
  880 fpmin3=min(f4,psias3)
      fpmax3=max(f4,psias3)
      if((psivl.lt.fpmin3).or.(psivl.gt.fpmax3))go to 890
      ip=ip+1
      yp(ip)=y2
      xp(ip)=((psivl-f4)/(f3-f4))*dx+x1
      if(ip.eq.2)go to 990
  890 ip=ip+1
      xp(ip)=x1
      yp(ip)=((psivl-f1)/(f4-f1))*dy+y1
      go to 990
!----------------------------------------------------------------------
!--   both horizontal and vertical asymptotes are present            --
!--   find psi on asymptotes                                         --
!----------------------------------------------------------------------
  900 psias1=((f2-f1)*(xasym-x1)/dx)+f1
      psias2=((f3-f2)*(yasym-y1)/dy)+f2
      psias4=((f4-f1)*(yasym-y1)/dy)+f1
      psias3=((f3-f4)*(xasym-x1)/dx)+f4
      if(xt.gt.xasym)go to 920
      if (yt.gt.yasym) go to 910
!----------------------------------------------------------------------
!--   xt.lt.xasym and yt.lt.yasym                                    --
!----------------------------------------------------------------------
      fpmin1=min(f1,psias1)
      fpmax1=max(f1,psias1)
      yp(1)=y1
      xp(1)=((psivl-f1)/(f2-f1))*dx+x1
      fpmin4=min(f1,psias4)
      fpmax4=max(f1,psias4)
      xp(2)=x1
      yp(2)=((psivl-f1)/(f4-f1))*dy+y1
      go to 990
!----------------------------------------------------------------------
!--   xt.lt.xasym and yt.gt.yasym                                    --
!----------------------------------------------------------------------
  910 fpmin4=min(f4,psias4)
      fpmax4=max(f4,psias4)
      xp(1)=x1
      yp(1)=((psivl-f1)/(f4-f1))*dy+y1
      fpmin3=min(psias3,f4)
      fpmax3=max(psias3,f4)
      yp(2)=y2
      xp(2)=((psivl-f4)/(f3-f4))*dx+x1
      go to 990
  920 if(yt.gt.yasym)go to 930
!----------------------------------------------------------------------
!--   xt.gt.xasym and yt.lt.yasym                                    --
!----------------------------------------------------------------------
      fpmin1=min(f2,psias1)
      fpmax1=max(f2,psias1)
      yp(1)=y1
      xp(1)=((psivl-f1)/(f2-f1))*dx+x1
      fpmin2=min(f2,psias2)
      fpmax2=max(f2,psias2)
      xp(2)=x2
      yp(2)=((psivl-f2)/(f3-f2))*dy+y1
      go to 990
!----------------------------------------------------------------------
!--   xt.gt.xasym and yt.ft.yasym                                    --
!----------------------------------------------------------------------
  930 fpmin2=min(f3,psias2)
      fpmax2=max(f3,psias2)
      xp(1)=x2
      yp(1)=((psivl-f2)/(f3-f2))*dy+y1
      fpmin3=min(f3,psias3)
      fpmax3=max(f3,psias3)
      yp(2)=y2
      xp(2)=((psivl-f4)/(f3-f4))*dx+x1
      go to 990
!----------------------------------------------------------------------
!--   no asymptotes                                                  --
!----------------------------------------------------------------------
  950 if((psivl.lt.fmin4).or.(psivl.gt.fmax4))go to 960
      ip=ip+1
      xp(ip)=x1
      yp(ip)=((psivl-f1)/(f4-f1))*dy+y1
      if(psivl.eq.f4)yp(ip)=y2
  960 if((psivl.lt.fmin2).or.(psivl.gt.fmax2))go to 970
      ip=ip+1
      xp(ip)=x2
      yp(ip)=((psivl-f2)/(f3-f2))*dy+y1
      if(psivl.eq.f3)yp(ip)=y2
      if(ip.eq.2)go to 990
  970 if((psivl.le.fmin1).or.(psivl.gt.fmax1))go to 980
      ip=ip+1
      yp(ip)=y1
      xp(ip)=((psivl-f1)/(f2-f1))*dx+x1
      if(ip.eq.2)go to 990
  980 ip=ip+1
      yp(ip)=y2
      xp(ip)=((psivl-f4)/(f3-f4))*dx+x1
  990 xt1=xp(1)
      xt2=xp(2)
      yt1=yp(1)
      yt2=yp(2)
      return
      end
      subroutine findax(nx,nz,x,y,xax,yax,psimx,psiout,xseps,yseps, &
        kaxis,xxout,yyout,kfound,psipsi,rmin,rmax, &
        zmin,zmax,zrmin,zrmax,rzmin,rzmax,dpsipsi, &
        bpoo,bpooz,limtrv,xlimv,ylimv,limfagv,ifit,jtime,kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          findax finds magnetic axis for arbitrary position       **
!**          greater than 3 cells from boundary of grid.             **
!**          note that nh2=2*nh, nwrk=2*(nw+1)*nh.                   **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**          kaxis = 10, bicubic spline only                         **
!**                  20, find axis and separatrix if any             **
!**                <  0, seperatrix only                             **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**          05/10/20..........Add Gradient Ascent Method as option  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky
      use set_kinds
      include 'eparmdud129.f90'
      include 'modules1.f90'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      dimension x(nx),y(nz),pds(6),xxout(*),yyout(*),psipsi(*)
      dimension xseps(1),yseps(1) ! this is an address of a location inside a 2-d array
      dimension bpoo(*),bpooz(*),pdss(6),xlimv(*),ylimv(*)
      dimension pdsold(6)
      data psitol/1.0e-04_dp/
      character(len=80) :: strtmp
      logical :: dodebugplts = .false. ! write surface files for debugging/plotting. Serial only, not parallel
      kerror = 0
      orelax = 1.0 ! Newton's Method relaxation constant (0.0-1.0)
      niter = 20   ! Number of iterations
      n111=1
      xseps(1)=-999.
      yseps(1)=-999.
      xseps(2)=-999.
      yseps(2)=-999.
      if (iabs(kaxis).lt.20) then
        !----------------------------------------------------------------------
        !--   fit 2-d spline to psi                                          --
        !----------------------------------------------------------------------
        !     psipsi - (in) psi function to spline, 1-d array (nx by nz)
        !     x,y - (in) 1-d array of coordinate values for function
        !     c - (out) 4-d array of spline coefficients
        !     bkx, bky - (out) interval coefficients w/ lkx,lky terms

        if (idebug >= 2) write (6,*)  'Entering sets2d'
        call sets2d(psipsi,c,x,nx,bkx,lkx,y,nz,bky,lky,wk,ier)
        if (idebug >= 2) then
          write (6,*) 'FINDAX Z,R = ', y(33),(x(i),i=45,45)
          write (6,*) 'FINDAX si = ',(psipsi((i-1)*nx+33),i=45,45)
          call seva2d(bkx,lkx,bky,lky,c,x(45),y(33),pds,ier,n111)
          write (6,*) 'FINDAX R,Z,si = ', x(45),y(33),pds(1)
          write (6,*) 'FINDAX lkx,lky = ',lkx,lky
          write (6,*) 'FINDAX lkx,lky,c = ',bkx(33),bky(33),c(1,33,1,33)
        endif
        if (kaxis.eq.10) return ! if only the bicubic spline is wanted return
      end if

      do n=1,kfound
        ! xxout,yyout - (in) interp points outlining (psipsi=0) the raised mag flux region
        ! pds - (out) interpolation value
        ! pds(1)=f, pds(2)=fx, pds(3)=fy, pds(4)=fxy, pds(5)=fxx, pds(6)=fyy
        call seva2d(bkx,lkx,bky,lky,c,xxout(n),yyout(n),pds,ier,n333)
        bpooz(n)=pds(2)/xxout(n)
        bpoo(n)=sqrt(bpooz(n)**2+(pds(3)/xxout(n))**2)
      end do

      sumip=0.
      do i=2,kfound
        delx=xxout(i)-xxout(i-1)
        dely=yyout(i)-yyout(i-1)
        dell=sqrt(delx**2+dely**2)
        abpol=(bpoo(i-1)+bpoo(i))/2.0
        sumip=sumip+abpol*dell
      end do
      sumip=sumip/tmu/twopi
      if (kaxis.le.0) go to 1000
!----------------------------------------------------------------------
!--   find magnetic axis, its elongation, and flux value             --
!----------------------------------------------------------------------
      if (negcur.eq.0) then
        psimx=-1.0e+10_dp
      else
        psimx=1.0e+10_dp
      endif
      ! Find psi max/min w/in r,z limits depending on current sign
      do i=1,nx
        do j=1,nz
          kk=(i-1)*nz+j
          if ((x(i).lt.rmin).or.(x(i).gt.rmax)) cycle
          if ((y(j).lt.zmin).or.(y(j).gt.zmax)) cycle
          if (psipsi(kk).le.psimx.and.negcur.eq.0) cycle
          if (psipsi(kk).ge.psimx.and.negcur.eq.1) cycle
          psimx=psipsi(kk)
          xax=x(i)
          yax=y(j)
        end do
      end do

      xs=xax
      ys=yax
      ps=psimx

      if (dodebugplts) then ! for debugging
        write(strtmp,'(a,i0.2,a,i0.2,a)') 'debug-surf',jtime,'-',ifit,'.txt'
        open(unit=99,file=trim(strtmp),status='replace')
        do iyplt = 1,nz
          do ixplt = 1,nx
            call seva2d(bkx,lkx,bky,lky,c,x(ixplt),y(iyplt),pds,ier,n666)
            write(99,'(3(1x,1pe12.5))') x(ixplt),y(iyplt),pds(1)
          end do
        end do
        close(unit=99)
        write(strtmp,'(a,i0.2,a,i0.2,a)') 'debug-conv',jtime,'-',ifit,'.txt'
        open(unit=99,file=trim(strtmp),status='replace')
      end if

      if (ifindopt==2) then
        xaxold = xax
        yaxold = yax
        pdsold = pds
        signcur = -1.0
        if (negcur.eq.0) signcur = 1.0
      end if

      errtmp = 0.0
      do j=1,niter
        ! pds(1)=f, pds(2)=fx, pds(3)=fy, pds(4)=fxy, pds(5)=fxx, pds(6)=fyy
        call seva2d(bkx,lkx,bky,lky,c,xax,yax,pds,ier,n666)
        if (dodebugplts) write(99,'(3(1x,1pe12.5))') xax,yax,pds(1) ! for debugging

        ! Gradient Ascent Method - better for sharp peaks
        if (ifindopt==2) then
          xerr=signcur*pds(2) ! find max or min depending on current direction
          yerr=signcur*pds(3)
          ! Adapt step size using Barzilai and Borwein approach
          dfx = pds(2) - pdsold(2)
          dfy = pds(3) - pdsold(3)
          if (j==1 .or. dfx**2+dfy**2<1.0e-15_dp) then
            gamman = 0.001_dp
          else
            gamman = abs((xax-xaxold)*dfx + (yax-yaxold)*dfy)/(dfx**2+dfy**2)
          endif
          xaxold = xax
          yaxold = yax
          pdsold = pds
          xax = xax + gamman*xerr
          yax = yax + gamman*yerr
          errtmp = gamman**2*(xerr**2+yerr**2)
          if (errtmp.lt.1.0e-12_dp) go to 310

        ! Original Newton's Method for optimization, xn+1 = xn - f'/f''
        else ! ifindopt==1
          det=pds(5)*pds(6)-pds(4)*pds(4)
          if (abs(det).lt.1.0e-15_dp) then
            kerror = 1 ! rls - previously continued to 305 with no error
            call errctrl_msg('findax','Newtons method to find magnetic axis has det=0')
            if (iand(iout,1).ne.0) write (nout,5000) xax,yax
            return
            !go to 305
          end if
          xerr=(-pds(2)*pds(6)+pds(4)*pds(3))/det
          yerr=(-pds(5)*pds(3)+pds(2)*pds(4))/det
          xax=xax+orelax*xerr
          yax=yax+orelax*yerr
          errtmp = xerr*xerr+yerr*yerr
          if (xax<x(1) .or. xax>x(nx) .or. yax<y(1) .or. yax>y(nz)) then
            kerror = 1 ! rls - previously continued to 305 with no error
            call errctrl_msg('findax','Newtons method to find separatrix point is off grid')
            if (iand(iout,1).ne.0) write (nout,5000) xax,yax
            return
            !go to 305
          end if
          if ((abs(pds(2)).lt.1.0e-06_dp).and.(abs(pds(3)).lt.1.0e-06_dp)) go to 310
          if (errtmp.lt.1.0e-12_dp) go to 310
        end if
      enddo
      if (errtmp.gt.1.0e-6_dp) then
        call errctrl_msg('findax','Iterative method to find magnetic axis reached max iterations',2)
      end if
  305 continue
      !if (iand(iout,1).ne.0) write (nout,5000) xax,yax
      xax=xs
      yax=ys
      psimx=pds(1)
      emaxis=1.3_dp
      go to 1000
  310 continue
      psimx=pds(1)
!----------------------------------------------------------------------
!--   compute elongation on axis                                     --
!----------------------------------------------------------------------
      thet=2.0*pds(4)/(pds(5)-pds(6))
      thet=0.5_dp*atan(thet)
      sint=sin(thet)
      cost=cos(thet)
      sint2=sint**2
      cost2=cost**2
      scost=sint*cost*2.0
      ar=pds(5)*cost2+pds(4)*scost+pds(6)*sint2
      az=pds(5)*sint2-pds(4)*scost+pds(6)*cost2
      siar=-0.5_dp*ar
      siaz=-0.5_dp*az
      emaxis=ar/az
      if (emaxis.gt.0.0) emaxis=sqrt(emaxis)
      if (emaxis.le.0.0) emaxis=1.3_dp
 1000 continue
      if (dodebugplts) then ! for debugging
        close(unit=99)
      end if
      delrmax1=0.40_dp
      delrmax2=0.40_dp
      sifsep=-1.e10_dp
      sissep=-1.e10_dp
      rfsep=-89.0
      zfsep=-89.0
      rssep=-89.0
      zssep=-89.0
      if (abs(dpsipsi).le.0.5_dp*psitol) return
!----------------------------------------------------------------------
!--   find the separatrix                                            --
!--   relaxd criteria for searching, 02/23/90                        --
!----------------------------------------------------------------------
      bpols=bpoo(1)
      ns=1
      xs=xxout(1)
      ys=yyout(1)
      do n=2,kfound
        if (bpoo(n).ge.bpols) cycle
        bpols=bpoo(n)
        xs=xxout(n)
        ys=yyout(n)
        ns=n
      end do

      errtmp = 0.0
      do j=1,niter
        if (xs.le.x(2) .or. xs.ge.x(nx-1) .or. ys.le.y(2) .or. ys.ge.y(nz-1)) then
          kerror = 1 ! rls - previously returned with no error
          call errctrl_msg('findax','1st separatrix point is off grid')
          if (iand(iout,1).ne.0) write (nout,5020) xs,ys
          return
        end if
        call seva2d(bkx,lkx,bky,lky,c,xs,ys,pds,ier,n666)
        det=pds(5)*pds(6)-pds(4)*pds(4)
        if (abs(det).lt.1.0e-15_dp) then
          kerror = 1 ! rls - previously returned with no error
          call errctrl_msg('findax','1st separatrix has det=0')
          if (iand(iout,1).ne.0) write (nout,5020) xs,ys
          return
        end if
        xerr=(-pds(2)*pds(6)+pds(4)*pds(3))/det
        yerr=(-pds(5)*pds(3)+pds(2)*pds(4))/det
        xs=xs+orelax*xerr
        ys=ys+orelax*yerr
        errtmp = xerr*xerr+yerr*yerr
        if (errtmp.lt.1.0e-12_dp*100.0) go to 1310
      end do ! j
      if (errtmp.gt.1.0e-6_dp*100.0) then
        call errctrl_msg('findax','Iterative method to find 1st separatrix reached max iterations',2)
      end if
      ! rls - previously returned with no error here, just because max iter=20 was reached
 1310 continue
!-----------------------------------------------------------------------
!--  found x separatrix point, check to see if inside vessel                     --
!-----------------------------------------------------------------------
      call zlim(zeross,n111,n111,limtrv,xlimv,ylimv,xs,ys,limfagv)
      if (zeross.le.0.1_dp) then
        kerror = 1 ! rls - previously returned with no error
        call errctrl_msg('findax','Separatrix point is not inside vessel, zeross.le.0.1')
        return
      end if
      xseps(1)=xs*100.
      yseps(1)=ys*100.
!-----------------------------------------------------------------------
!--  consider x separatrix point on surface if psi/dpsi/dR < 0.004 a             --
!-----------------------------------------------------------------------
      anow=(rmax-rmin)*0.5_dp
      znow=0.5_dp*(zmin+zmax)
      relpsi=abs((pds(1)-psiout))
      call seva2d(bkx,lkx,bky,lky,c,rmax,znow,pdss,ier,n333)
      delrmax1=relpsi/abs(pdss(2))
      relpsi=relpsi/abs((psimx-psiout))
      if (delrmax1.gt.0.004_dp*anow) then
        !kerror = 1 ! rls - previously returned with no error
        call errctrl_msg('findax','Separatrix point is not on surface, delrmax1.gt.0.004_dp*anow',2)
        return
      end if
      sifsep=pds(1)
      rfsep=xs
      zfsep=ys
      psiout=pds(1)
      xxout(ns)=xs
      yyout(ns)=ys
      xxout(ns-1)=0.5_dp*(xxout(ns)+xxout(ns-2))
      yyout(ns-1)=0.5_dp*(yyout(ns)+yyout(ns-2))
      xxout(ns+1)=0.5_dp*(xxout(ns)+xxout(ns+2))
      yyout(ns+1)=0.5_dp*(yyout(ns)+yyout(ns+2))
      rmin=xxout(1)
      rmax=xxout(1)
      zmin=yyout(1)
      zmax=yyout(1)
      rzmin=xxout(1)
      rzmax=xxout(1)
      zrmin=yyout(1)
      zrmax=yyout(1)
      bpave=0.
      do i=2,kfound
        if (xxout(i).lt.rmin) zrmin=yyout(i)
        if (xxout(i).gt.rmax) zrmax=yyout(i)
        if (yyout(i).lt.zmin) rzmin=xxout(i)
        if (yyout(i).gt.zmax) rzmax=xxout(i)
        bpave=bpave+bpoo(i)
        rmin=min(rmin,xxout(i))
        rmax=max(rmax,xxout(i))
        zmin=min(zmin,yyout(i))
        zmax=max(zmax,yyout(i))
      end do
!----------------------------------------------------------------
!-- find tracing points                                        --
!----------------------------------------------------------------
      jwant=(zrmax-y(1)+1.e-6_dp)/(y(2)-y(1))+1
      rminmax=0.5_dp*(rmin+rmax)
      zmaxfs=y(jwant)
      zminfs=zmaxfs
      do i=1,kfound-1
        zminus=zmaxfs-yyout(i)
        zplus=zmaxfs-yyout(i+1)
        if (xxout(i).gt.rminmax.and.zminus*zplus.le.0.0) then
         rmaxfs=xxout(i)+(xxout(i+1)-xxout(i))/(yyout(i+1)- &
                  yyout(i))*zminus
        elseif (xxout(i).lt.rminmax.and.zminus*zplus.le.0.0) then
         rminfs=xxout(i)+(xxout(i+1)-xxout(i))/(yyout(i+1)- &
                  yyout(i))*zminus
        endif
      enddo
!
      znow=(zmax+zmin)/2.
      anow=(rmax-rmin)/2.
      bpave=bpave/real(kfound-1,dp)
!-----------------------------------------------------------------------
!-- find possible second separatrix                                   --
!-----------------------------------------------------------------------
      bpmins=10.
      ns=-1
      do i=2,kfound
        if ((ys-znow)*(yyout(i)-znow).lt.0.0.and.bpoo(i).lt.bpmins) then
          bpmins=bpoo(i)
          ns=i
        endif
      end do
      if (ns.eq.-1) then
        kerror = 1 ! rls - previously returned with no error
        call errctrl_msg('findax','2nd separatrix not found, ns.eq.-1')
        return
      end if
      xs=xxout(ns)
      ys=yyout(ns)
!
      errtmp = 0.0
      do j=1,niter
        if (xs.le.x(2) .or. xs.ge.x(nx-1) .or. ys.le.y(2) .or. ys.ge.y(nz-1)) then
          !kerror = 1 ! rls - previously returned with no error
          call errctrl_msg('findax','2nd separatrix point is off grid',2)
          if (iand(iout,1).ne.0) write (nout,5025) xs,ys
          return
        end if
        call seva2d(bkx,lkx,bky,lky,c,xs,ys,pds,ier,n666)
        det=pds(5)*pds(6)-pds(4)*pds(4)
        if (abs(det).lt.1.0e-15_dp) then
          kerror = 1 ! rls - previously returned with no error
          call errctrl_msg('findax','2nd separatrix has det=0')
          if (iand(iout,1).ne.0) write (nout,5025) xs,ys
          return
        end if
        xerr=(-pds(2)*pds(6)+pds(4)*pds(3))/det
        yerr=(-pds(5)*pds(3)+pds(2)*pds(4))/det
        xs=xs+orelax*xerr
        ys=ys+orelax*yerr
        errtmp = xerr*xerr+yerr*yerr
        if (errtmp.lt.1.0e-12_dp*100.0) go to 9310
      end do
      if (errtmp.gt.1.0e-6_dp*100.0) then
        call errctrl_msg('findax','Iterative method to find 2nd separatrix reached max iterations',2)
      end if
      ! rls - previously returned with no error here, just because max iter=20 was reached
 9310 continue
!-----------------------------------------------------------------------
!--  make sure 2nd seperatrix inside vessel                               --
!-----------------------------------------------------------------------
      call zlim(zeross,n111,n111,limtrv,xlimv,ylimv,xs,ys,limfagv)
      ! TODO: Previously this check was allowed to return without error and no message.
      ! It occurs a lot, so no error and supress message for now.
      if (zeross.le.0.1_dp) then
        !kerror = 1 ! rls - previously returned with no error
        call errctrl_msg('findax','2nd seperatrix point is not inside vessel, zeross.le.0.1',2)
        return
      end if
      if (abs(ys*100.-yseps(1)).lt.2.0*anow) then
        !kerror = 1 ! rls - previously returned with no error
        call errctrl_msg('findax','2nd seperatrix too far away',2)
        return
      end if
      ! TODO: If 2nd separatrix errors out (returns) above, the following variables
      ! (xseps=-999, rssep=-89, sissep=-1.0e-10, delrmax2=0.4) have default values that
      ! are used elsewhere. Check that they do not cause problems.
      xseps(2)=xs*100.
      yseps(2)=ys*100.
      rssep=xs
      zssep=ys
      sissep=pds(1)
!-----------------------------------------------------------------------
!--  consider x point on surface if psi/dpsi/dR < 0.004 a             --
!-----------------------------------------------------------------------
      relpsi=abs((pds(1)-psiout))
      delrmax2=relpsi/abs(pdss(2))
      relpsi=relpsi/abs((psimx-psiout))
      if (delrmax2.gt.0.004_dp*anow) then
        !kerror = 1 ! rls - previously returned with no error. Occurs a lot, so no error, supress message
        !call errctrl_msg('findax','2nd separatrix point is not on surface, delrmax2.gt.0.004_dp*anow',2)
        return
      end if

      xxout(ns)=xs
      yyout(ns)=ys
      xxout(ns-1)=0.5_dp*(xxout(ns)+xxout(ns-2))
      yyout(ns-1)=0.5_dp*(yyout(ns)+yyout(ns-2))
      xxout(ns+1)=0.5_dp*(xxout(ns)+xxout(ns+2))
      yyout(ns+1)=0.5_dp*(yyout(ns)+yyout(ns+2))
      rmin=xxout(1)
      rmax=xxout(1)
      zmin=yyout(1)
      zmax=yyout(1)
      rzmin=xxout(1)
      rzmax=xxout(1)
      zrmin=yyout(1)
      zrmax=yyout(1)
      do i=2,kfound
        if (xxout(i).lt.rmin) zrmin=yyout(i)
        if (xxout(i).gt.rmax) zrmax=yyout(i)
        if (yyout(i).lt.zmin) rzmin=xxout(i)
        if (yyout(i).gt.zmax) rzmax=xxout(i)
        rmin=min(rmin,xxout(i))
        rmax=max(rmax,xxout(i))
        zmin=min(zmin,yyout(i))
        zmax=max(zmax,yyout(i))
      end do

      return

 5000 format (/,1x,'no convergence to magnetic axis, rax, yax = ', &
              2(1x,e10.3))
 5020 format (/,1x,'no convergence to separatrix, rs, ys = ', &
              2(1x,e10.3))
 5025 format (/,1x,'no convergence to 2nd septrx, rs, ys = ', &
              2(1x,e10.3))
      end

      subroutine fqlin(x1,y1,x2,y2,f1,f2,f3,f4,x,y,area,psivl)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      a1=(x2-x)*(y2-y)
      a2=(x-x1)*(y2-y)
      a3=(x-x1)*(y-y1)
      a4=(x2-x)*(y-y1)
      psivl=(a1*f1+a2*f2+a3*f3+a4*f4)/area
      return
      end

      subroutine maxpsi(xl1,yl1,xl2,yl2,x1,y1,x2,y2,f1,f2,f3,f4, &
                        area,psimax,xtry,ytry,nerr)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          maxpsi finds the largest psi value along the line       **
!**          segment y=alpha*x+beta joining the two points           **
!**          (xl1,yl1) and (xl2,yl2).                                **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**          24/07/85..........revised                               **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use set_kinds
      use error_control
      implicit integer*4 (i-n), real*8 (a-h, o-z)
!
      dx=xl2-xl1
      if(dx.eq.0.0)go to 100
      dy=yl2-yl1
      if(dy.eq.0.0)go to 100
      c=f1+f3-(f2+f4)
      if(c.eq.0.0)go to 100
      a=y2*(f2-f1)+y1*(f4-f3)
      b=x1*(f2-f3)+x2*(f4-f1)
      alpha=dy/dx
      secder=2.*alpha*c
      if(secder.gt.0.0)go to 100
      beta=yl1-alpha*xl1
      xcrit=-(b*alpha+c*beta+a)/secder
      if((xcrit.gt.xl2).or.(xcrit.lt.xl1))go to 100
      ycrit=alpha*xcrit+beta
      if((ycrit.lt.y1).or.(ycrit.gt.y2)) then
        nerr = 3
        call errctrl_msg('maxpsi','ycrit is out of bounds')
        return
      end if
      xl2=xcrit
      yl2=ycrit
      psip1=-1.0e+35_dp
      go to 110
  100 call fqlin(x1,y1,x2,y2,f1,f2,f3,f4,xl1,yl1,area,psip1)
  110 call fqlin(x1,y1,x2,y2,f1,f2,f3,f4,xl2,yl2,area,psip2)
      psimax=max(psip1,psip2)
      if (psimax.ne.psip1) go to 120
      xtry=xl1
      ytry=yl1
      return
  120 xtch=xl2
      ytch=yl2
      return
      end

      subroutine minmax(psi,nwh,nh,kn,psivl,iflag,knew)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          minmax finds minimum and maximum value of psi in a      **
!**          cell.                                                   **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          07/09/83..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension psi(nwh)
      iflag=0
      fmin=min(psi(kn),psi(kn+1),psi(kn+nh),psi(kn+nh+1))
      fmax=max(psi(kn),psi(kn+1),psi(kn+nh),psi(kn+nh+1))
      if ((psivl.lt.fmin).or.(psivl.gt.fmax)) iflag=1
      knew=kn
      return
      end
      subroutine order(xp,yp,np)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          order puts the points (xp,yp) in increasing order       **
!**          of yp.                                                  **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          22/03/84..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension xp(*),yp(*)
      nptr=np
   80 is=0
      nptr=nptr-1
      do 90 k=1,nptr
      if (yp(k+1).ge.yp(k)) go to 90
      is=1
      xs=xp(k+1)
      ys=yp(k+1)
      xp(k+1)=xp(k)
      yp(k+1)=yp(k)
      xp(k)=xs
      yp(k)=ys
   90 continue
      if (is.eq.1) go to 80
      return
      end
      subroutine packps(xp,yp,np,rm,zm,kadd)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          pack orders the points (xp,yp) in sequencial order.     **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          22/03/84..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: cjrf,wxin,wyin,wxout,wyout
      include 'eparmdud129.f90'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      dimension xp(*),yp(*)
      data iflag/2/
!
      if (iflag.eq.2) go to 2000
      ymin=yp(1)
      ymax=yp(1)
      rymin=xp(1)
      rymax=xp(1)
      do 1350 i=2,np
        if (yp(i).lt.ymin) rymin=xp(i)
        if (yp(i).gt.ymax) rymax=xp(i)
        ymin=min(ymin,yp(i))
        ymax=max(ymax,yp(i))
 1350 continue
      slope=(rymax-rymin)/(ymax-ymin)
      kin=0
      kout=0
      do 1500 i=1,np
        rcut=rymin+(yp(i)-ymin)*slope
        if (xp(i).ge.rcut) go to 1400
        kin=kin+1
        wxin(kin)=xp(i)
        wyin(kin)=yp(i)
        go to 1500
 1400   kout=kout+1
        wxout(kout)=xp(i)
        wyout(kout)=yp(i)
 1500 continue
      call order(wxin,wyin,kin)
      call order(wxout,wyout,kout)
      do 1600 k=1,kin
        xp(k)=wxin(k)
        yp(k)=wyin(k)
 1600 continue
      do 1700 k=1,kout
        kk=k+kin
        mm=kout-k+1
        xp(kk)=wxout(mm)
        yp(kk)=wyout(mm)
 1700 continue
      if (kadd.eq.0) return
      np=kk+1
      xp(np)=xp(1)
      yp(np)=yp(1)
      return
!
 2000 continue
      do 2200 i=1,np
        wxin(i)=atan2d(yp(i)-zm,rm-xp(i))
        if (wxin(i).lt.0.0) wxin(i)=wxin(i)+360.0
 2200 continue
      nptr=np
 2800 is=0
      nptr=nptr-1
      do 2900 k=1,nptr
      if (wxin(k+1).ge.wxin(k)) go to 2900
      is=1
      xs=xp(k+1)
      ys=yp(k+1)
      ws=wxin(k+1)
      xp(k+1)=xp(k)
      yp(k+1)=yp(k)
      wxin(k+1)=wxin(k)
      xp(k)=xs
      yp(k)=ys
      wxin(k)=ws
 2900 continue
      if (is.eq.1) go to 2800
      if (kadd.eq.0) return
      np=np+1
      xp(np)=xp(1)
      yp(np)=yp(1)
      return
      end
      subroutine qfit(k,x1,x2,x3,y1,y2,y3,x,y,yp,ierr)
!**********************************************************************
!**                                                                  **
!**     main program:  mhd fitting code                              **
!**                                                                  **
!**     subprogram description:                                      **
!**          QFIT is a quadratic fitter.                             **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**         1. (x1,y1); (x2,y2); (x3,y3)                             **
!**             three points used to determine parabola.             **
!**              note: x1.lt.x2.lt.x3.                               **
!**         2. function switch                                       **
!**            k=1  find y,yp given x                                **
!**            k=2  find x,yp given y                                **
!**            k=3  find a,b,c (returned as x,y,yp)                  **
!**            where y(x)=a*x**2+b*x+c                               **
!**                  yp  = derivative of parabola at (x,y)           **
!**         3. ierr=1  error return                                  **
!**                                                                  **
!**     REFERENCES:                                                  **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          04/08/86..........first created                         **
!**                                                                  **
!**********************************************************************
      use set_kinds
      use error_control
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      ierr=0
      if (x1.eq.x2.or.x2.eq.x3.or.x1.eq.x3) then
        ierr=1
        call errctrl_msg('qfit','at least one pair of x-coordinates are the same')
        return
      endif
      alpha=y1/((x1-x2)*(x1-x3))
      beta=y2/((x2-x1)*(x2-x3))
      gamma=y3/((x3-x1)*(x3-x2))
!
      a=alpha+beta+gamma
      b=-(gamma+beta)*x1-(gamma+alpha)*x2-(alpha+beta)*x3
      c=x2*x3*alpha+x1*x3*beta+x1*x2*gamma
!
      select case (k)
      case (1)
        go to 10
      case (2)
        go to 20
      case (3)
        go to 30
      end select

10    y=a*x*x+b*x+c
      go to 4
20    rad=sqrt(b*b-4.0*a*(c-y))
      root1=(-b+rad)*.5_dp/a
      root2=(-b-rad)*.5_dp/a
      t1=(root1-x1)*(x3-root1)
      t2=(root2-x1)*(x3-root2)
      zero=-x1*1.0e-7_dp
      if (t1.ge.zero) go to 1
      if (t2.ge.zero) go to 2
      ierr=1
      call errctrl_msg('qfit','t1<0 or t2<0')
      return
1     if (t2.ge.zero) go to 3
      x=root1
      go to 4
2     x=root2
      go to 4
3     x=min(root1,root2)
4     yp=2.0*a*x+b
      return
30    x=a
      y=b
      yp=c
      return
      end

      subroutine surfac(siwant,psi,nw,nh,rgrid,zgrid,xout,yout, &
                    nfound,npoint,drgrid,dzgrid,xmin, &
                    xmax,ymin,ymax,ipack,rmaxis,zmaxis,negcur,kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          surfac generates a contour of constant psi of value     **
!**          siwant.                                                 **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use set_kinds
      use error_control
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension xout(*),yout(*),psi(*),rgrid(*),zgrid(*)

      kerror = 0
      n111=1
      if (negcur.eq.0) then
        curneg=1.
      else
        curneg=-1.
      endif
      nfound=0
      do 200 i=2,nw-1
      do 200 j=1,nh
        if (rgrid(i).lt.xmin) go to 200
        if (rgrid(i).gt.xmax) go to 200
        if (zgrid(j).lt.ymin) go to 200
        if (zgrid(j).gt.ymax) go to 200
        kk1=(i-1)*nh+j
        df1=siwant-psi(kk1)
        if (df1*curneg.lt.0.0.and.rgrid(i-1).lt.xmin) then
          kk2x=kk1
          df2x=df1
          kk1x=kk1-nh
          df1x=siwant-psi(kk1x)
          if (df1x*df2x.le.0.0) then
            if (nfound+1.gt.npoint-1) go to 200
            nfound=nfound+1
            xout(nfound)=rgrid(i-1)+df1x*drgrid/(psi(kk2x)-psi(kk1x))
            yout(nfound)=zgrid(j)
          endif
        endif
        kk2=i*nh+j
        df2=siwant-psi(kk2)
        if (df1*df2.gt.0.0) go to 200
        if (nfound+1.gt.npoint-1) go to 200
        nfound=nfound+1
        xout(nfound)=rgrid(i)+df1*drgrid/(psi(kk2)-psi(kk1))
        yout(nfound)=zgrid(j)
  200 continue
      do 300 i=1,nw
      do 300 j=2,nh-1
        if (rgrid(i).lt.xmin) go to 300
        if (rgrid(i).gt.xmax) go to 300
        if (zgrid(j).lt.ymin) go to 300
        if (zgrid(j).gt.ymax) go to 300
        kk1=(i-1)*nh+j
        df1=siwant-psi(kk1)
        if (df1*curneg.lt.0.0.and.zgrid(j-1).lt.ymin) then
          kk2x=kk1
          df2x=df1
          kk1x=kk1-1
          df1x=siwant-psi(kk1x)
          if (df1x*df2x.le.0.0) then
            if (nfound+1.gt.npoint-1) go to 300
            nfound=nfound+1
            xout(nfound)=rgrid(i)
            yout(nfound)=zgrid(j-1)+df1x*dzgrid/(psi(kk2x)-psi(kk1x))
          endif
        endif
        kk2=(i-1)*nh+j+1
        df2=siwant-psi(kk2)
        if (df1*df2.gt.0.0) go to 300
        if (nfound+1.gt.npoint-1) go to 300
        nfound=nfound+1
        xout(nfound)=rgrid(i)
        yout(nfound)=zgrid(j)+df1*dzgrid/(psi(kk2)-psi(kk1))
  300 continue
      if (ipack.gt.0) call packps(xout,yout,nfound,rmaxis,zmaxis,n111)
      if (nfound.lt.3) then
        kerror = 1
        call errctrl_msg('surfac','Less than 3 contour points found')
        return
      end if
      return
      end

      subroutine zlim(zero,nw,nh,limitr,xlim,ylim,x,y,iflag)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          zlim determines whether points on the (x,y) grid are    **
!**          inside or outside of the boundary set by the limiters.  **
!**                                                                  **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       zero (out)......1 if inside and 0 otherwise                **
!**       nw..............dimension of x                             **
!**       nh..............dimension of y                             **
!**       limitr..........number of limiter points                   **
!**       xlim............r coordinates of limiter                   **
!**       ylim............z coordinates of limiter                   **
!**       x...............r grid                                     **
!**       y...............z grid                                     **
!**       iflag...........1 convex geometry                          **
!**                       2 general geometry                         **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**          13/08/85..........iflag added                           **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use set_kinds
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      real*8 :: zero(nw*nh),x(nw),y(nh) ! sometimes an array, other times a constant
      dimension xlim(*),ylim(*)
      logical b,c,d,inside,bold

      select case (iflag)
        case (1)
          go to 10
        case (2)
          go to 200
      end select

10    continue
      kk = 0
      do i = 1,nw
      do j = 1,nh
        kk = kk + 1
        zero(kk) = 1.
        ncross = 0
        do 20 k = 1,limitr-1
          if ((ylim(k).lt.y(j)) .and. (ylim(k+1).lt.y(j))) go to 20
          if (x(i) .eq. xlim(k))  go to 20
          t = x(i) - xlim(k)
          s = xlim(k+1) - x(i)
          if ((t*s) .lt. 0.) go to 20
          di = (ylim(k+1)-ylim(k)) / (xlim(k+1)-xlim(k))
          f = ylim(k) + di*(x(i)-xlim(k))
          if (f .lt. y(j)) go to 20
          ncross = ncross + 1
20      continue
        mcross = .5_dp*ncross ! truncates to integer
        mcross = 2*mcross
        if (ncross .eq. mcross) zero(kk) = 0.
      end do
      end do
      return
!
  200 continue
      kk=0
      do i=1,nw
        do j=1,nh
          kk=kk+1
          d=.false.
          b=.true.
          n=0
          inside=.false.
          bold=b
          do k=1,limitr-1
            c=.false.
            !---------------------------------------------------------------------------
            !--  fixed if test logic, for ge and le per Wolfe of MIT, 93/09/02        --
            !--       if (y(j).le.ylim(k).and.y(j).ge.ylim(k+1)                       --
            !--  .        .or.y(j).ge.ylim(k).and.y(j).le.ylim(k+1)) then             --
            !---------------------------------------------------------------------------
            if (y(j).le.ylim(k).and.y(j).gt.ylim(k+1) &
              .or.y(j).ge.ylim(k).and.y(j).lt.ylim(k+1)) then
              c=.true.
              d=.true.
              n=n+1
            endif
            if(c.and. &
              (y(j)-ylim(k))*(xlim(k+1)-xlim(k))- &
              (ylim(k+1)-ylim(k))*(x(i)-xlim(k)).gt.0.) &
              b=.not.b
            if (n.eq.2) then
              n=0
              if (bold.eqv.b) inside=.true.
              bold=b
            endif
          end do
          zero(kk)=0.0
          if (inside.and.d) zero(kk)=1.0
        end do
      end do
      return
      end
