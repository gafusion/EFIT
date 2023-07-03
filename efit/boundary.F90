#include "config.f"
!*********************************************************************
!>     Bound finds the outermost contour of a function        
!!     given in the array psi. note that the contoured     
!!     function psi is assumed to be a decreasing function   
!!     of distance away from the axis.                 
!!     if ix<0, then BOUND just traces the field line  
!!     starting at (xctr,yctr) either to closure or to the 
!!     wall. ix=-1, trace clockwise. -2, counter clockwise. 
!!                                                                  
!!     calling arguments:
!!       @param psi : function to be contoured 
!!       @param nw : row dimension of psi
!!       @param nh : column dimension of psi  
!!       @param nwh : nw x nh
!!       @param psivl : psi value on the outermost contour
!!       @param xmin : minimum r value on contour
!!       @param xmax : maximum r value
!!       @param ymin : minimum z value
!!       @param ymax : maximum z value
!!       @param zero : weighting array
!!       @param x : r grid
!!       @param y : z grid
!!       @param xctr : r guess of contour center
!!       @param yctr : z guess of contour center
!!       @param ix : flag
!!       @param limitr : number of limiter points
!!       @param xlim : r coordinates of limiter points
!!       @param ylim : z coordinates of limiter points
!!       @param xcontr: output r coordinates of contour
!!       @param ycontr : output z coordinates of contour
!!       @param ncontr : number of contour points found
!!       @param xlmin : minimum or maximum limiter coordinates
!!       @param npoint : maximum number of contour points
!!       @param dpsi : difference of psi values on contour
!!                       and on limiter.  dpsi is used to
!!                       distinguish between a limited and a
!!                       diverted plasma 
!!       @param rymin : r at ymin 
!!       @param rymax : r at ymax
!!       @param zxmin : z at xmin
!!       @param zxmax : z at xmax
!!       @param tolbndpsi : tolerance on psi
!!
!*********************************************************************
      subroutine bound(psi,nw,nh,nwh,psivl,xmin,xmax,ymin,ymax, &
           zero,x,y,xctr,yctr,ix,limitr,xlim,ylim,xcontr,ycontr, &
           ncontr,xlmin,npoint,rymin,rymax,dpsi,zxmin,zxmax,nerr, &
           ishot,itime,limfag,radold,kbound,tolbndpsi)
      use error_control
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      double precision, dimension(nwh) :: psi, zero
      double precision, dimension(nwh) :: x,y
      double precision, dimension(npoint) :: xcontr,ycontr
      double precision, dimension(5) :: dist(5)
      double precision, dimension(limitr) :: xlim,ylim
      double precision, dimension(1) :: zerol(1),xt(1),yt(1),rad(1),yctr(1)
      parameter (n111=1,nloop=60,etolc=1.e-06_dp,etol=1.e-04_dp)
      data mecopy/0/
      save dx,dy,carea,rmid

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
!--   nerr=10000, negative plasma current                            --
!----------------------------------------------------------------------
      nosign=0
      if (nerr .eq. 10000) then
        nosign=1
        do i=1,nwh
          psi(i)=-psi(i)
        end do
      end if
!------------------------------------------------------------------------
!--   BBrown's version of BOUND                                        --
!------------------------------------------------------------------------
      if (ix.ge.0 .and. kbound.ne.0) then
        !call old_new    (psi,nw,nh,nwh,psivl,xmin,xmax,ymin,ymax, &
        !     zero,x,y,xctr,yctr,ix,limitr,xlim,ylim,xcontr,ycontr, &
        !     ncontr,xlmin,npoint,rymin,rymax,dpsi,zxmin,zxmax,nerr, &
        !     ishot,itime,limfag,radold,kbound)
        if (ncontr.lt.3) then
          nerr=3
          call errctrl_msg('bound','Less than 3 contour points found')
        end if
        if (nosign.eq.1) then
          do i=1,nwh
            psi(i)=-psi(i)
          end do
          !psivl=-psivl ! not yet set
        end if
        return
      end if

      nerr=0
      psib0=-1.e+10_dp
      rad(1)=radold
      rin=xctr
      rout=xlmin
!---------------------------------------------------------------------
!--   field line tracing ix < 0                                     --
!---------------------------------------------------------------------
      if(ix.lt.0) rad(1)=xctr
      if (mecopy.le.0) then
        dx=x(2)-x(1)
        dy=y(2)-y(1)
        carea=dx*dy
        rmid=1.02_dp*(x(1)+x(nw))/2.0
        mecopy=1
      end if
!----------------------------------------------------------------------
!--   find starting value of psi                                     --
!----------------------------------------------------------------------
      psiloop: do loop=1,nloop
        i=1+(rad(1)-x(1))/dx
        if(rad(1)-x(i).lt.0.0) i=i-1
        j=1+(yctr(1)-y(1))/(dy-0.000001_dp)
        jjj=j

        if ((ix.eq.-2).or.(ix.lt.-2.and.rad(1).gt.rmid)) then
          j=j-1
          jjj=j+1
        end if

        if ((ix.ge.0).and.(rad(1).gt.rmid)) then
          j=j-1
          jjj=j+1
        end if
        j=min(max(j,1),nh)
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

        xt(1)=rad(1)
        jjj=min(max(jjj,1),nh)
        yt(1)=y(jjj)
        ncontr=1
        xcontr(ncontr)=xt(1)
        ycontr(ncontr)=yt(1)
        a3=(xt(1)-x(i))*dy
        a4=carea-a3
        psivl=(psi(kk+nh+1)*a3+psi(kk+1)*a4)/carea
        if(yt(1).eq.y(j)) psivl=psi(kk)+(psi(kk+nh)-psi(kk))*(xt(1)-x(i))/dx
        zsum=1.e10_dp
        zerol=1.e10_dp

        contr: do while (.true.)
          f1=psi(kk)
          f2=psi(kk+nh)
          f3=psi(kk+nh+1)
          f4=psi(kk+1)
          x1=x(i)
          x2=x(i+1)
          y1=y(j)
          y2=y(j+1)
          if (ncontr.ne.1) then
!----------------------------------------------------------------------
!--         check for proximity to corner                            --
!----------------------------------------------------------------------
            dist(1)=(xt(1)-x1)**2+(yt(1)-y1)**2
            dist(2)=(xt(1)-x2)**2+(yt(1)-y1)**2
            dist(3)=(xt(1)-x2)**2+(yt(1)-y2)**2
            dist(4)=(xt(1)-x1)**2+(yt(1)-y2)**2
            dist(5)=min(dist(1),dist(2),dist(3),dist(4))
            if (dist(5).le.etolc) then
              do l=1,4
                kj=l
                if(dist(l).eq.dist(5)) exit
              end do
!----------------------------------------------------------------------
!--           kj points to appropriate corner                        --
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
            end if
          end if
!----------------------------------------------------------------------
!--       check for limiter in cell (fieldline hit the wall)
!----------------------------------------------------------------------
          zsum=zero(kk)+zero(kk+1)+zero(kk+nh)+zero(kk+nh+1)
          if(zsum.eq.0.0) exit contr

          if (abs(zsum-4.0).ge.1.e-03_dp) then
!-------------------------------------------------------------------------
!--         from one to three corners of cell are inside limiter.  get max
!--         psi on line segment of limiter in cell and compare this max
!--         with current value of psilim (or psisep)
!--         note: do loop index assumes point 'limitr' is the same as
!--         point 'limitr-1'
!-------------------------------------------------------------------------
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
!--------------------------------------------------------------------------
!--           at least one limiter point is not in cell.  subroutine cellb
!--           returns intersections of cell boundaries and line defined by
!--           points k and k+1 or one cell boundary and one interior point.
!--           ifail=1 if points k and k+1 do not intersect the current cell
!--           of interest.
!--------------------------------------------------------------------------
              if ((kij1.ne.kk).or.(kij2.ne.kk)) then
                ifail=0
                call cellb(xc1,yc1,xc2,yc2,x1,y1,x2,y2,ifail)
                if(ifail.eq.1) cycle ! line segment does not intersect cell
              end if
              ! psilm is largest psi value along line segment betw pts
              call maxpsi(xc1,yc1,xc2,yc2,x1,y1,x2,y2,f1,f2,f3,f4,carea,psilm,xtry1,ytry1,nerr)
              if (nerr.gt.0) then
                if (ncontr.lt.3) then
                  nerr=3
                  call errctrl_msg('bound', &
                    'Less than 3 contour points found')
                end if
                if (nosign.eq.1) then
                  do i=1,nwh
                    psi(i)=-psi(i)
                  end do
                  psivl=-psivl
                end if
                return
              end if
              psilx=max(psilm,psilx)
              if (psilx.le.psilm) then
                xtry=xtry1
                ytry=ytry1
              end if
            end do ! limitr

            if (psilx.eq.-1.e10_dp) then
              nerr=3
              call errctrl_msg('bound', &
                'Limiter points do not intersect cell')
              if (ncontr.lt.3) then
                nerr=3
                call errctrl_msg('bound', &
                  'Less than 3 contour points found')
              end if
              if (nosign.eq.1) then
                do i=1,nwh
                  psi(i)=-psi(i)
                end do
                psivl=-psivl
              end if
              return
            end if

            ! check if the current point is outside the limiter
            dpsi=min(dpsi,abs(psivl-psilx))
            if (psilx-psivl.ge.tolbndpsi) then
              call zlim(zerol,n111,n111,limitr,xlim,ylim,xt,yt,limfag)
              if(zerol(1).le.0.01_dp) exit contr
            end if
          end if
          call extrap(f1,f2,f3,f4,x1,y1,x2,y2,xt(1),yt(1),xt1,yt1,xt2,yt2, &
                      psivl,dx,dy)
!----------------------------------------------------------------------
!--       decide which intersection (xt1,yt1) or (xt2,yt2) is required
!----------------------------------------------------------------------
          dist1=(yt1-yt(1))**2+(xt1-xt(1))**2
          dist2=(yt2-yt(1))**2+(xt2-xt(1))**2
          if (dist1.lt.dist2) then
            yt(1)=yt2
            xt(1)=xt2
          else
            yt(1)=yt1
            xt(1)=xt1
          end if
          ncontr=ncontr+1
          if (ncontr.gt.npoint) then
            nerr=3
            call errctrl_msg('bound', &
              'Number of contour points greater than max allowed')
            if (ncontr.lt.3) then
              nerr=3
              call errctrl_msg('bound', &
                'Less than 3 contour points found')
            end if
            if (nosign.eq.1) then
              do i=1,nwh
                psi(i)=-psi(i)
              end do
              psivl=-psivl
            end if
            return
          end if
          xcontr(ncontr)=xt(1)
          ycontr(ncontr)=yt(1)

          ! Debug tool: Write out the contour coordinates for each loop (iteration)
          !write(*,*) loop,ncontr,xcontr(ncontr),ycontr(ncontr)

!----------------------------------------------------------------------
!--       find next cell                                             --
!----------------------------------------------------------------------
          if(xt(1).eq.x2) i=i+1
          if(xt(1).eq.x1) i=i-1
          if(yt(1).eq.y2) j=j+1
          if(yt(1).eq.y1) j=j-1

          if (ix.ge.0.or.ix.lt.-2) then
            if(yt(1).lt.ymin) rymin=xt(1)
            if(yt(1).gt.ymax) rymax=xt(1)
            if(xt(1).lt.xmin) zxmin=yt(1)
            if(xt(1).gt.xmax) zxmax=yt(1)
            xmin=min(xmin,xt(1))
            xmax=max(xmax,xt(1))
            ymin=min(yt(1),ymin)
            ymax=max(yt(1),ymax)
          end if

          kold=kk
!----------------------------------------------------------------------
!         find new cell index                                        --
!----------------------------------------------------------------------
          kk=(i-1)*nh+j
          if(kk.eq.kstrt) exit contr
          dis2p=sqrt((xcontr(1)-xt(1))**2+(ycontr(1)-yt(1))**2)
          if((dis2p.lt.0.1_dp*dx).and.(ncontr.gt.5)) exit contr
        end do contr

!----------------------------------------------------------------------
!--     psi on boundary smaller than psi on limiter, decrease rad and--
!--     try again.                                                   --
!----------------------------------------------------------------------
        if ((kk.ne.kstrt).and.((dis2p.ge.0.1_dp*dx).or.(ncontr.le.5))) then
          psib0=psivl
          !
          if (ix.lt.0) then ! ix, -1=trace clockwise, -2=counter clockwise
            if (ncontr.lt.3) then
              nerr=3
              call errctrl_msg('bound', &
                'Less than 3 contour points found')
            end if
            if (nosign.eq.1) then
              do i=1,nwh
                psi(i)=-psi(i)
              end do
              psivl=-psivl
            end if
            return
          end if
          !
          if(loop.ge.nloop) exit psiloop
          rout=rad(1)
          rad(1)=(rin+rout)*0.5_dp
          cycle psiloop
        end if

        if(loop.ge.nloop) exit psiloop
!----------------------------------------------------------------------
!--     check for convergence of boundary                            --
!----------------------------------------------------------------------
        err=abs((psivl-psib0)/psivl)
        if (ix.lt.0) then
          if(ix.lt.-2) dpsi=1.e-06_dp
          if (ncontr.lt.3) then
            nerr=3
            call errctrl_msg('bound', &
              'Less than 3 contour points found')
          end if
          if (nosign.eq.1) then
            do i=1,nwh
              psi(i)=-psi(i)
            end do
            psivl=-psivl
          end if
          return
        end if
        if (err.le.etol) then
          ! if the fieldline hit the wall try moving inside separatrix
          if (zsum.eq.0.0 .or. zerol(1).le.0.01_dp) then
            if (rad(1).gt.xctr) then
              rad(1)=rad(1)-1.5*etol*rad(1)
            else
              rad(1)=rad(1)+1.5*etol*rad(1)
            endif
          else
            exit psiloop
          endif
        endif
!----------------------------------------------------------------------
!--     new rad,psi and try again                                    --
!----------------------------------------------------------------------
        psib0=psivl
        call zlim(zerol,n111,n111,limitr,xlim,ylim,rad,yctr,limfag)
        if (zerol(1).le.0.01_dp) then
          rout=rad(1)
        else
          rin=rad(1)
        end if
        rad(1)=(rin+rout)*0.5_dp
      end do psiloop

      radold=rad(1)
      psib0=psivl
      if ((abs(ycontr(1)-ycontr(ncontr)).gt.0.5_dp*dy) .or. &
          (abs(xcontr(1)-xcontr(ncontr)).gt.0.5_dp*dx)) then
         nerr=3
         call errctrl_msg('bound', &
           'First and last contour points are too far apart')
      end if

      if (ncontr.lt.3) then
        nerr=3
        call errctrl_msg('bound','Less than 3 contour points found')
      end if
      if (nosign.eq.1) then
        do i=1,nwh
          psi(i)=-psi(i)
        end do
        psivl=-psivl
      end if
      return
      end subroutine bound

!*********************************************************************
!>    cellb redefines (xc1,yc1) and/or (xc2,yc2) so that
!!    they are intersections of cell boundaries unless
!!    one of the points is an interior point in which
!!    case it is not disturbed.                                    
!!                                                                  
!!    @param xc1 : R of point on limiter 
!!    @param yc1 : Z of point on limiter
!!    @param xc2 : R of adjacent point on limiter
!!    @param yc2 : Z of adjacent point on limiter
!!    @param x1 : R of left side of cell
!!    @param y1 : Z of bottom of cell
!!    @param x2 : R of right side of cell
!!    @param y2 : Z of top of cell
!!    @param ifail : 0 if limiter points do not intersect cell
!!                   1 if limiter points do not intersect cell
!!                                                                 
!*********************************************************************
      subroutine cellb(xc1,yc1,xc2,yc2,x1,y1,x2,y2,ifail)
      implicit none
      integer*4, intent(out) :: ifail
      real*8, intent(in) :: x1,x2,y1,y2
      real*8, intent(inout) :: xc1,xc2,yc1,yc2
      real*8 alpha,beta,dx,dy,x1b,x2b,xs1,xs2,y1b,y2b,ys1,ys2

      ifail=0
      if((xc1.lt.x1).and.(xc2.lt.x1)) ifail=1 ! left of cell
      if((xc1.gt.x2).and.(xc2.gt.x2)) ifail=1 ! right of cell
      if((yc1.lt.y1).and.(yc2.lt.y1)) ifail=1 ! below cell
      if((yc1.gt.y2).and.(yc2.gt.y2)) ifail=1 ! above cell
      if(ifail.eq.1) return
      dx=xc1-xc2
      dy=yc1-yc2
      if((dx.eq.0.0).and.(dy.eq.0.0)) ifail=1 ! duplicate limiter point
      if(ifail.eq.1) return
      if (dx.ne.0.0) then
      if (dy.eq.0.0) then
!----------------------------------------------------------------------
!--     check if intersection exists for horizontal line            --
!----------------------------------------------------------------------
        ! this is a duplicate check that should never be
        ! encountered...
        !if ((yc1.lt.y1).or.(yc1.gt.y2)) then
        !  ifail=1
        !  return
        !end if
!----------------------------------------------------------------------
!--     is there an interior point ?                                 --
!----------------------------------------------------------------------
        if ((x1.lt.xc1).and.(xc1.le.x2)) then
!----------------------------------------------------------------------
!--       point (xc1,yc1) is interior point                          --
!----------------------------------------------------------------------
          if(xc1.lt.xc2) xc2=x2
          if(xc1.gt.xc2) xc2=x1
        else if ((x1.le.xc2).and.(xc2.le.x2)) then
!----------------------------------------------------------------------
!--       point (xc2,yc2) is interior point                          --
!----------------------------------------------------------------------
          if(xc2.gt.xc1) xc1=x1
          if(xc2.lt.xc1) xc1=x2
        else
!----------------------------------------------------------------------
!--       no interior points                                         --
!----------------------------------------------------------------------
          xc1=x1
          xc2=x2
        end if
        return
      end if ! dy.eq.0.0
!----------------------------------------------------------------------
!--   line is inclined. get equation y=alpha*x+beta.                 --
!----------------------------------------------------------------------
      alpha=dy/dx
      beta=yc1-alpha*xc1
      y1b=alpha*x1+beta ! Z of line at left side of cell
      if (y1b.lt.y1) then
        x1b=(y1-beta)/alpha ! R of line at bottom of cell
        if ((x1b.lt.x1).or.(x1b.gt.x2)) then
          ! line goes below cell
          ifail=1
          return
        end if
        ! line intersects bottom of cell
        xs1=x1b
        ys1=y1
      else if (y1b.gt.y2) then
        x1b=(y2-beta)/alpha ! R of line at top of cell
        if ((x1b.lt.x1).or.(x1b.gt.x2)) then
          ! line goes above cell
          ifail=1
          return
        end if
        ! line intersects top of cell
        xs1=x1b
        ys1=y2
      else
        ! line intersects left side of cell
        xs1=x1
        ys1=y1b
      end if
      y2b=alpha*x2+beta ! Z of line at right side of cell
      if (y2b.lt.y1) then
        x2b=(y1-beta)/alpha ! R of line at bottom of cell
        if ((x2b.gt.x2).or.(x2b.lt.x1)) then
          ! line goes below cell
          ifail=1
          return
        end if
        ! line intersects bottom of cell
        xs2=x2b
        ys2=y1
      else if (y2b.gt.y2) then
        x2b=(y2-beta)/alpha ! R of line at top of cell
        if ((x2b.gt.x2).or.(x2b.lt.x1)) then
          ! line goes above cell
          ifail=1
          return
        end if
        ! line intersects top of cell
        xs2=x2b
        ys2=y2
      else
        ! line intersects right side of cell
        xs2=x2
        ys2=y2b
      end if
!----------------------------------------------------------------------
!--   at this point we have intersection (xs1,ys1) and (xs2,ys2)     --
!--   check for interior point                                       --
!----------------------------------------------------------------------
      if ((x1.lt.xc1).and.(xc1.lt.x2)) then
        if ((y1.lt.yc1).and.(yc1.lt.y2)) then
!----------------------------------------------------------------------
!--       point (xc1,yc1) is interior point                          --
!----------------------------------------------------------------------
          if (yc2.gt.yc1) then
            xc2=xs2
            yc2=ys2
            if(ys1.gt.ys2) yc2=ys1
            if(ys1.gt.ys2) xc2=xs1
          else
            xc2=xs1
            yc2=ys1
            if(ys1.gt.ys2) yc2=ys2
            if(ys1.gt.ys2) xc2=xs2
          end if
          return
        end if
      else if ((x1.lt.xc2).and.(xc2.lt.x2)) then
        if ((y1.lt.yc2).and.(yc2.lt.y2)) then
!----------------------------------------------------------------------
!--       point (xc2,yc2) is interior point                          --
!----------------------------------------------------------------------
          if (yc1.gt.yc2) then
            xc1=xs2
            yc1=ys2
            if(ys1.gt.ys2) yc1=ys1
            if(ys1.gt.ys2) xc1=xs1
          else
            xc1=xs1
            yc1=ys1
            if(ys1.gt.ys2) yc1=ys2
            if(ys1.gt.ys2) xc1=xs1
          end if
          return
        end if
      end if
      xc1=xs1
      yc1=ys1
      xc2=xs2
      yc2=ys2
      return
      end if ! dx.ne.0.0
!----------------------------------------------------------------------
!--   check if intersection exists for vertical line                 --
!----------------------------------------------------------------------
      ! this is a duplicate check that should never be
      ! encountered...
      !if ((xc1.lt.x1).or.(xc1.gt.x2)) then
      !  ifail=1
      !  return
      !end if
!----------------------------------------------------------------------
!--   is there an interior point ?                                   --
!----------------------------------------------------------------------
      if ((y1.le.yc1).and.(yc1.le.y2)) then
!----------------------------------------------------------------------
!--     point (xc1,yc1) is interior point                            --
!----------------------------------------------------------------------
        if(yc2.gt.yc1) yc2=y2
        if(yc2.lt.yc1) yc2=y1
      else if ((y1.le.yc2).and.(yc2.le.y2)) then
!----------------------------------------------------------------------
!--     point  (xc2,yc2) is interior point                           --
!----------------------------------------------------------------------
        if(yc2.gt.yc1) yc1=y1
        if(yc2.lt.yc1) yc1=y2
      else
!----------------------------------------------------------------------
!--     no interior points                                           --
!----------------------------------------------------------------------
        yc1=y1
        yc2=y2
      end if
      return
      end subroutine cellb


!*********************************************************************
!>    chkcrn decides which cell is next when an error in
!!    the cell step is likely.
!!     
!!    @param  psi : nw x nh psi grid 
!!
!!    @param  nwh : total number of grid point nw * nh
!!
!!    @param  psivl : value of psi at boundary
!!
!!    @param  kold : 
!!
!!    @param  knew : 
!!
!!    @param  icrnr : 
!!
!!    @param  kk : 
!!
!!    @param  nh : Number of horizontal grid points
!!
!!    @param  i : 
!!
!!    @param  il :  
!! 
!*********************************************************************
      subroutine chkcrn(psi,nwh,psivl,kold,knew,icrnr,kk,nh,i,i1)
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      integer :: nwh
      double precision, dimension(nwh) :: psi
      knew=0
      i1=0
!----------------------------------------------------------------------
!--   cell #1                                                        --
!----------------------------------------------------------------------
      kn=kk
      if (kn.ne.kold) then
        call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
        if (iflag.ne.1) then
          i1=i
          return
        end if
      end if
      select case (icrnr)
      case default ! 1
!----------------------------------------------------------------------
!--     corner #1                                                    --
!--     cell #2                                                      --
!----------------------------------------------------------------------
        kn=kk-nh
        if (kn.ne.kold) then
          call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
          if (iflag.ne.1) then
            i1=i-1
            return
          end if
        end if
!----------------------------------------------------------------------
!--     cell #3                                                      --
!----------------------------------------------------------------------
        kn=kk-nh-1
        if (kold.ne.kn) then
          call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
          if (iflag.ne.1) then
            i1=i-1
            return
          end if
        end if
!----------------------------------------------------------------------
!--     cell #4                                                      --
!----------------------------------------------------------------------
        kn=kk-1
        call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
        i1=i-1
      case (2)
!----------------------------------------------------------------------
!--     corner #2                                                    --
!--     cell #2                                                      --
!----------------------------------------------------------------------
        kn=kk-1
        if (kn.ne.kold) then
          call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
          if (iflag.ne.1) then
            i1=i
            return
          end if
        end if
!----------------------------------------------------------------------
!--     cell #3                                                      --
!----------------------------------------------------------------------
        kn=kk+nh-1
        if (kn.ne.kold) then
          call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
          if (iflag.ne.1) then
            i1=i+1
            return
          end if
        end if
!----------------------------------------------------------------------
!--     cell #4                                                      --
!----------------------------------------------------------------------
        kn=kk+nh
        call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
        i1=i+1
      case (3)
!----------------------------------------------------------------------
!--     corner #3                                                    --
!--     cell #2                                                      --
!----------------------------------------------------------------------
        kn=kk+nh
        if(kn.ne.kold) then
          call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
          if (iflag.ne.1) then
            i1=i+1
            return
          end if
        end if
!----------------------------------------------------------------------
!--     cell #3                                                      --
!----------------------------------------------------------------------
        kn=kk+nh+1
        if (kn.ne.kold) then
          call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
          if (iflag.ne.1) then
            i1=i+1
            return
          end if
        end if
!----------------------------------------------------------------------
!--     cell #4                                                      --
!----------------------------------------------------------------------
        kn=kk+1
        call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
        i1=i
      case (4)
!----------------------------------------------------------------------
!--     corner #4                                                    --
!--     cell #2                                                      --
!----------------------------------------------------------------------
        kn=kk+1
        if (kn.ne.kold) then
          call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
          if (iflag.ne.1) then
            i1=i
            return
          end if
        end if
!----------------------------------------------------------------------
!--     cell #3                                                      --
!----------------------------------------------------------------------
        kn=kk-nh+1
        if (kn.ne.kold) then
          call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
          if (iflag.ne.1) then
            i1=i-1
            return
          end if
        end if
!----------------------------------------------------------------------
!--     cell #4                                                      --
!----------------------------------------------------------------------
        kn=kk-nh
        call minmax(psi,nwh,nh,kn,psivl,iflag,knew)
        i1=i-1
      end select
      return
      end subroutine chkcrn

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
      subroutine cntour(xaxd,yaxd,psivl,xemin,xemax,yemin,yemax, &
      yxmin,yxmax,xymin,xymax,dang,arcl,bperr,dx,dy,xmin,xmax, &
      ymin,ymax,iauto,iautoc,xc,yc,ipts,x,nw,y,nh,cspln,n2cspln,&
      nh2,itty,iptsm,negcur,bkx,lkx,bky,lky,kerror)
      use global_constants
      use error_control
      use eparm, only: kubicx,lubicx,kubicy,lubicy
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      integer*4, intent(inout) :: kerror
      real*8 :: xc(*),yc(*)
      real*8, dimension(6) :: pds(6)
      real*8, dimension(lubicx+1) :: bkx(lubicx+1),bky(lubicy+1)
      dimension cspln(kubicx,lubicx,kubicy,lubicy)
      real*8 piov2,piov4,fpiov4,spiov4,tpiov4,tpiov2
      parameter (n111=1,n333=3)

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
!---  serrt is absolute error convergence criteria for newtons method below.--
!---  get psi at (xaxd,yaxd)                                                --
!-----------------------------------------------------------------------------
      call seva2d(bkx,lkx,bky,lky,cspln,xaxd,yaxd,pds,ier,n111)
      if (negcur.eq.0) then
        psiaxd=pds(1)
      else
        psiaxd=-pds(1)
      end if
      if (psiaxd.lt.psivl) then
        kerror = 1
        call errctrl_msg('cntour','psiaxd.lt.psivl')
        return
      end if
!----------------------------------------------------------------------
!---  loop over theta from 0 to twopi                                --
!----------------------------------------------------------------------
   10 thet=thet+dthet
      thet=min(thet,twopi)
      if (thet.eq.twopi) then
!--------------------------------------------------------------------
!---    close contour                                              --
!--------------------------------------------------------------------
        ipts=ipts+1
        xc(ipts)=xc(1)
        yc(ipts)=yc(1)
        return
      end if
!----------------------------------------------------------------------
!---  get equation of ray emanating from (xaxd,yaxd)                 --
!----------------------------------------------------------------------
      if (((piov4.le.thet).and.(thet.le.tpiov4)).or. &
          ((fpiov4.le.thet).and.(thet.le.spiov4))) then
!----------------------------------------------------------------------
!---    x as a function of y            x=a*y+bincp                  --
!----------------------------------------------------------------------
        isgn=1
        if (thet.gt.pi) isgn=-1
        if (isgn.eq.-1) then
          thet1=tpiov2-thet
          if (thet.gt.tpiov2) thet1=pi-abs(thet1)
        else
          thet1=piov2-thet
          if (thet.gt.piov2) thet1=twopi-abs(thet1)
        end if
        a=tan(thet1)
        iflg=1
        bincp=xaxd-a*yaxd
      else
!----------------------------------------------------------------------
!---    y as a function of x            y=a*x+bincp                  --
!----------------------------------------------------------------------
        isgn=-1
        if ((thet.lt.piov4).or.(thet.gt.spiov4)) isgn=1
        a=tan(thet)
        iflg=0
        bincp=yaxd-a*xaxd
      end if
!-----------------------------------------------------------------------
!---  now have y=a*x+bincp    (iflg=0)   or                           --
!---  x=a*y+bincp             (iflg=1)                                --
!-----------------------------------------------------------------------
      x1=xaxd
      y1=yaxd
      cost=cos(thet)
      sint=sin(thet)
      psi1=psiaxd
!------------------------------------------------------------------------
!---  sliding interval search. max width of interval ~1.41*(dx or dy)  --
!------------------------------------------------------------------------
   40 if (iflg.eq.1) then
!---    search in y
        y2=y1+isgn*dy
        x2=a*y2+bincp
      else
!---    search in x
        x2=x1+isgn*dx
        y2=a*x2+bincp
      end if
      if ((x2.lt.xmin).or.(x2.gt.xmax).or. &
          (y2.lt.ymin).or.(y2.gt.ymax)) then
!---------------------------------------------------------
!---    errors                                          --
!---------------------------------------------------------
        if (iauto.ne.1) then
          kerror = 1
          call errctrl_msg('cntour', &
            'flux surface is outside search area, iauto.ne.1')
          return
        end if

        psivl0=psivl
        dapsi=psiaxd-psivl0
        psivl=psivl0+dapsi*0.0005_dp
        iautoc=1
        write (itty,1020) psivl0,psivl
 1020   format(2x,'boundary search, will change psilim from', &
               /,e16.8,'  to  ',e16.8,'  and try again')
        go to 1040
      end if
      call seva2d(bkx,lkx,bky,lky,cspln,x2,y2,pds,ier,n111)
      if (negcur.eq.0) then
        psi2=pds(1)
      else
        psi2=-pds(1)
      end if
      dpsi=(psivl-psi1)*(psivl-psi2)
      if (dpsi.gt.0.0) then
        x1=x2
        y1=y2
        psi1=psi2
        go to 40
      end if
!--------------------------------------------------------------------------
!---  now have psivl between psi1 and psi2,converge using newton-raphson --
!--------------------------------------------------------------------------
      newti=0
      if (iflg.eq.1) then
        yn=y1+isgn*dy*0.5_dp
        xn=a*yn+bincp
      else
        xn=x1+isgn*dx*0.5_dp
        yn=a*xn+bincp
      end if
   80 call seva2d(bkx,lkx,bky,lky,cspln,xn,yn,pds,ier,n333)
      if (negcur.eq.0) then
       dpsids=pds(2)*cost+pds(3)*sint
       dpsi=pds(1)-psivl
      else
       dpsids=-(pds(2)*cost+pds(3)*sint)
       dpsi=-pds(1)-psivl
      end if
      serr=-dpsi/dpsids
      if ((abs(serr).ge.serrt).and.(abs(dpsi).ge.derrt)) then
        delx=serr*cost
        dely=serr*sint
        xn=xn+delx
        yn=yn+dely
        newti=newti+1
        if (newti.ge.20) then
          kerror = 1
          call errctrl_msg('cntour','newti.ge.20')
          return
        end if
        go to 80
      end if
!-----------------------------------------------------------------------------
!---  end of newton iteration
!---  check for sufficient accuracy in point spacing as determined by thet
!---  accuracy test is based on a relative error in poloidal b field of bperr
!-----------------------------------------------------------------------------
      bp2=(1./xn)*sqrt(pds(2)**2+pds(3)**2)
      if (thet.ne.0.0) then
        if (abs((bp2-bp1)/max(bp2,bp1)).ge.bperr) then
!---      avoid accumulation at x points
          ihalf=ihalf+1
          if (ihalf.le.4) then
!---        spacing too large for grad psi. decrease theta and try again
            thet=thet-dthet
            dthet=dthet*0.5_dp
            go to 10
          end if
        end if
      end if
      bp1=bp2
      ipts=ipts+1
      if (ipts.gt.iptsm-1) then
        kerror = 1
        call errctrl_msg('cntour','ipts.gt.iptsm-1')
        return
      end if
      xc(ipts)=xn
      yc(ipts)=yn
      xemin=min(xemin,xn)
      xemax=max(xemax,xn)
      if (xemax.eq.xn) yxmax=yn
      if (xemin.eq.xn) yxmin=yn
      yemin=min(yemin,yn)
      yemax=max(yemax,yn)
      if (yemax.eq.yn) xymax=xn
      if (yemin.eq.yn) xymin=xn
!---  limit angle increment to approximately arcl meters per step
      rad1=sqrt((xn-xaxd)**2+(yn-yaxd)**2)
      dthet=min(dthet0,arcl/rad1)
      ihalf=0
      go to 10
      end subroutine cntour

!**********************************************************************
!>
!!    extrap extrapolates across a (x,y) cell.
!!    
!!    @param f1 :
!!    @param f2 :
!!    @param f3 :
!!    @param f4 :
!!    @param x1 :
!!    @param y1 :
!!    @param x2 :
!!    @param y2 :
!!    @param xt :
!!    @param yt :
!!    @param xt1 :
!!    @param yt1 :
!!    @param xt2 :
!!    @param yt2 :
!!    @param psivl :
!!    @param dx :
!!    @param dy :
!**********************************************************************
      subroutine extrap(f1,f2,f3,f4,x1,y1,x2,y2,xt,yt,xt1,yt1, &
                        xt2,yt2,psivl,dx,dy)
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
      if (c.ne.0.0) then
      xasym=-b/c
      yasym=-a/c
      if ((x1.le.xasym).and.(xasym.le.x2)) then
        if ((y1.lt.yasym).and.(yasym.lt.y2)) then
!----------------------------------------------------------------------
!--       both horizontal and vertical asymptotes are present        --
!--       find psi on asymptotes                                     --
!----------------------------------------------------------------------
          psias1=((f2-f1)*(xasym-x1)/dx)+f1
          psias2=((f3-f2)*(yasym-y1)/dy)+f2
          psias4=((f4-f1)*(yasym-y1)/dy)+f1
          psias3=((f3-f4)*(xasym-x1)/dx)+f4
          if (xt.gt.xasym) then
            if (yt.gt.yasym) then
!----------------------------------------------------------------------
!--           xt.gt.xasym and yt.gt.yasym                            --
!----------------------------------------------------------------------
              fpmin2=min(f3,psias2)
              fpmax2=max(f3,psias2)
              xp(1)=x2
              yp(1)=((psivl-f2)/(f3-f2))*dy+y1
              fpmin3=min(f3,psias3)
              fpmax3=max(f3,psias3)
              yp(2)=y2
              xp(2)=((psivl-f4)/(f3-f4))*dx+x1
            else
!----------------------------------------------------------------------
!--           xt.gt.xasym and yt.lt.yasym                            --
!----------------------------------------------------------------------
              fpmin1=min(f2,psias1)
              fpmax1=max(f2,psias1)
              yp(1)=y1
              xp(1)=((psivl-f1)/(f2-f1))*dx+x1
              fpmin2=min(f2,psias2)
              fpmax2=max(f2,psias2)
              xp(2)=x2
              yp(2)=((psivl-f2)/(f3-f2))*dy+y1
            end if
          else if (yt.gt.yasym) then
!----------------------------------------------------------------------
!--         xt.lt.xasym and yt.gt.yasym                              --
!----------------------------------------------------------------------
            fpmin4=min(f4,psias4)
            fpmax4=max(f4,psias4)
            xp(1)=x1
            yp(1)=((psivl-f1)/(f4-f1))*dy+y1
            fpmin3=min(psias3,f4)
            fpmax3=max(psias3,f4)
            yp(2)=y2
            xp(2)=((psivl-f4)/(f3-f4))*dx+x1
          else
!----------------------------------------------------------------------
!--         xt.lt.xasym and yt.lt.yasym                              --
!----------------------------------------------------------------------
            fpmin1=min(f1,psias1)
            fpmax1=max(f1,psias1)
            yp(1)=y1
            xp(1)=((psivl-f1)/(f2-f1))*dx+x1
            fpmin4=min(f1,psias4)
            fpmax4=max(f1,psias4)
            xp(2)=x1
            yp(2)=((psivl-f1)/(f4-f1))*dy+y1
          end if
        else
!----------------------------------------------------------------------
!--       vertical asymptote                                         --
!         find psi value on asymptote at y1 and y2
!----------------------------------------------------------------------
          psias1=((f2-f1)*(xasym-x1)/dx)+f1
          psias3=((f3-f4)*(xasym-x1)/dx)+f4
          if (xt.lt.xasym) then
!----------------------------------------------------------------------
!--         left side of cell                                        --
!----------------------------------------------------------------------
            fpmin1=min(f1,psias1)
            fpmax1=max(f1,psias1)
            if ((psivl.ge.fpmin1).and.(psivl.le.fpmax1)) then
              ip=ip+1
              yp(ip)=y1
              xp(ip)=((psivl-f1)/(f2-f1))*dx+x1
            end if
            fpmin3=min(f4,psias3)
            fpmax3=max(f4,psias3)
            if ((psivl.ge.fpmin3).and.(psivl.le.fpmax3)) then
              ip=ip+1
              yp(ip)=y2
              xp(ip)=((psivl-f4)/(f3-f4))*dx+x1
              if (ip.eq.2) then
                xt1=xp(1)
                xt2=xp(2)
                yt1=yp(1)
                yt2=yp(2)
                return
              end if
            end if
            ip=ip+1
            xp(ip)=x1
            yp(ip)=((psivl-f1)/(f4-f1))*dy+y1
          else
!----------------------------------------------------------------------
!--         right side of cell                                       --
!----------------------------------------------------------------------
            fpmin1=min(psias1,f2)
            fpmax1=max(psias1,f2)
            if ((psivl.ge.fpmin1).and.(psivl.le.fpmax1)) then
              ip=ip+1
              yp(ip)=y1
              xp(ip)=((psivl-f1)/(f2-f1))*dx+x1
              if (psivl.eq.f2) xp(ip)=x2
            end if
            fpmin3=min(f3,psias3)
            fpmax3=max(f3,psias3)
            if ((psivl.ge.fpmin3).and.(psivl.le.fpmax3)) then
              ip=ip+1
              yp(ip)=y2
              xp(ip)=((psivl-f4)/(f3-f4))*dx+x1
              if (psivl.eq.f3) xp(ip)=x2
              if (ip.eq.2) then
                xt1=xp(1)
                xt2=xp(2)
                yt1=yp(1)
                yt2=yp(2)
                return
              end if
            end if
            ip=ip+1
            xp(ip)=x2
            yp(ip)=((psivl-f2)/(f3-f2))*dy+y1
          end if
        end if
        xt1=xp(1)
        xt2=xp(2)
        yt1=yp(1)
        yt2=yp(2)
        return
      end if
      end if ! c.ne.0.0
      if ((yasym.lt.y1).or.(yasym.gt.y2)) then
!----------------------------------------------------------------------
!--     no asymptotes                                                --
!----------------------------------------------------------------------
        if ((psivl.ge.fmin4).and.(psivl.le.fmax4)) then
          ip=ip+1
          xp(ip)=x1
          yp(ip)=((psivl-f1)/(f4-f1))*dy+y1
          if (psivl.eq.f4) yp(ip)=y2
        end if
        if ((psivl.ge.fmin2).and.(psivl.le.fmax2)) then
          ip=ip+1
          xp(ip)=x2
          yp(ip)=((psivl-f2)/(f3-f2))*dy+y1
          if (psivl.eq.f3)yp(ip)=y2
          if (ip.eq.2) then
            xt1=xp(1)
            xt2=xp(2)
            yt1=yp(1)
            yt2=yp(2)
            return
          end if
        end if
        if ((psivl.gt.fmin1).and.(psivl.le.fmax1)) then
          ip=ip+1
          yp(ip)=y1
          xp(ip)=((psivl-f1)/(f2-f1))*dx+x1
          if (ip.eq.2) then
            xt1=xp(1)
            xt2=xp(2)
            yt1=yp(1)
            yt2=yp(2)
            return
          end if
        end if
        ip=ip+1
        yp(ip)=y2
        xp(ip)=((psivl-f4)/(f3-f4))*dx+x1
      else
!----------------------------------------------------------------------
!--     there is a horizontal asymptote                              --
!--     find psi value on asymptote at x1 and x2                     --
!----------------------------------------------------------------------
        psias4=((f4-f1)*(yasym-y1)/dy)+f1
        psias2=((f3-f2)*(yasym-y1)/dy)+f2
        if (yt.ge.yasym) then
!----------------------------------------------------------------------
!--       upper section of cell                                      --
!----------------------------------------------------------------------
          fpmin4=min(psias4,f4)
          fpmax4=max(f4,psias4)
          if ((psivl.ge.fpmin4).and.(psivl.le.fpmax4)) then
            ip=ip+1
            xp(ip)=x1
            yp(ip)=((psivl-f1)/(f4-f1))*dy+y1
            if (psivl.eq.f4) yp(ip)=y2
          end if
          fpmin2=min(f3,psias2)
          fpmax2=max(f3,psias2)
          if ((psivl.ge.fpmin2).and.(psivl.le.fpmax2)) then
            ip=ip+1
            xp(ip)=x2
            yp(ip)=((psivl-f2)/(f3-f2))*dy+y1
            if (psivl.eq.f3) yp(ip)=y2
            if (ip.eq.2) then
              xt1=xp(1)
              xt2=xp(2)
              yt1=yp(1)
              yt2=yp(2)
              return
            end if
          end if
          ip=ip+1
          yp(ip)=y2
          xp(ip)=((psivl-f4)/(f3-f4))*dx+x1
        else
!----------------------------------------------------------------------
!--       lower section of cell                                      --
!----------------------------------------------------------------------
          fpmin4=min(f1,psias4)
          fpmax4=max(f1,psias4)
          if ((psivl.ge.fpmin4).and.(psivl.le.fpmax4)) then
            ip=ip+1
            xp(ip)=x1
            yp(ip)=((psivl-f1)/(f4-f1))*dy+y1
          end if
          fpmin2=min(f2,psias2)
          fpmax2=max(f2,psias2)
          if ((psivl.ge.fpmin2).and.(psivl.le.fpmax2)) then
            ip=ip+1
            xp(ip)=x2
            yp(ip)=((psivl-f2)/(f3-f2))*dy+y1
            if (ip.eq.2) then
              xt1=xp(1)
              xt2=xp(2)
              yt1=yp(1)
              yt2=yp(2)
              return
            end if
          end if
          ip=ip+1
          yp(ip)=y1
          xp(ip)=((psivl-f1)/(f2-f1))*dx+x1
        end if
      end if
      xt1=xp(1)
      xt2=xp(2)
      yt1=yp(1)
      yt2=yp(2)
      return
      end subroutine extrap

!>*********************************************************************
!!
!!    findax finds magnetic axis for arbitrary position 
!!    greater than 3 cells from boundary of grid.
!!    note that nh2=2*nh, nwrk=2*(nw+1)*nh.
!! 
!!    @param nx :
!!    @param nz :
!!    @param x : 1-d array of coordinate values for psipsi
!!    @param y : 1-d array of coordinate values for psipsi
!!    @param xax :
!!    @param yax :
!!    @param psimx :
!!    @param psiout : ??
!!    @param xseps : coordinates of separtrices (if requested/found)
!!    @param yseps : coordinates of separtrices (if requested/found)
!!    @param kaxis : kaxis = 10, bicubic spline only 
!!                           20, find axis and separatrix if any 
!!                           <0, seperatrix only 
!!    @param xxout : coordinates of raised mag flux region outline (psipsi=0)
!!    @param yyout : coordinates of raised mag flux region outline (psipsi=0)
!!    @param kfound : number of points defining the raised mag flux region outline (psipsi=0)
!!    @param psipsi : psi function, 1-d array (nx by nz)
!!    @param rmin :
!!    @param rmax :
!!    @param zmin :
!!    @param zmax :
!!    @param zrmin :
!!    @param zrmax :
!!    @param rzmin :
!!    @param rzmax :
!!    @param dpsipsi : difference of psi values on contour
!!                       and on limiter.  dpsipsi is used to
!!                       distinguish between a limited and a
!!                       diverted plasma 
!!    @param bpoo :
!!    @param bpooz :
!!    @param limtrv :
!!    @param xlimv :
!!    @param ylimv :
!!    @param limfagv :
!!    @param ifit :
!!    @param jtime : time index
!!    @param kerror : error flag
!*********************************************************************
      subroutine findax(nx,nz,x,y,xax,yax,psimx,psiout,xseps,yseps, &
        kaxis,xxout,yyout,kfound,psipsi,rmin,rmax, &
        zmin,zmax,zrmin,zrmax,rzmin,rzmax,dpsipsi, &
        bpoo,bpooz,limtrv,xlimv,ylimv,limfagv,ifit,jtime,kerror)
      use commonblocks,only: c,wk,bkx,bky
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      real*8, intent(out) :: xseps(2),yseps(2)
      dimension x(nx),y(nz),pds(6),xxout(kfound),yyout(kfound),psipsi(nx*nz)
      dimension bpoo(kfound),bpooz(kfound),pdss(6),xlimv(limtrv),ylimv(limtrv)
      dimension pdsold(6),zeross(1),xs(1),ys(1)
      character(len=80) :: strtmp
      parameter (psitol=1.0e-04_dp)
      kerror=0
      orelax=1.0 ! Newton's Method relaxation constant (0.0-1.0)
      niter=20   ! Number of iterations
      xseps(1)=-999.
      yseps(1)=-999.
      xseps(2)=-999.
      yseps(2)=-999.
      if (iabs(kaxis).lt.20) then
!----------------------------------------------------------------------
!--     fit 2-d spline to psi                                        --
!----------------------------------------------------------------------
!       psipsi - (in) psi function to spline, 1-d array (nx by nz)
!       x,y - (in) 1-d array of coordinate values for function
!       c - (out) 4-d array of spline coefficients
!       bkx, bky - (out) interval coefficients w/ lkx,lky terms

#ifdef DEBUG_LEVEL2
        write (6,*)  'Entering sets2d'
#endif
        call sets2d(psipsi,c,x,nx,bkx,lkx,y,nz,bky,lky,wk,ier)
#ifdef DEBUG_LEVEL2
        write (6,*) 'FINDAX Z,R = ', y(33),(x(i),i=45,45)
        write (6,*) 'FINDAX si = ',(psipsi((i-1)*nx+33),i=45,45)
        call seva2d(bkx,lkx,bky,lky,c,x(45),y(33),pds,ier,n111)
        write (6,*) 'FINDAX R,Z,si = ', x(45),y(33),pds(1)
        write (6,*) 'FINDAX lkx,lky = ',lkx,lky
        write (6,*) 'FINDAX lkx,lky,c = ',bkx(33),bky(33),c(1,33,1,33)
#endif
        if(kaxis.eq.10) return ! if only the bicubic spline is wanted return
      endif

      do n=1,kfound
        ! xxout,yyout - (in) interp points outlining (psipsi=0) the raised mag flux region
        ! pds - (out) interpolation value
        ! pds(1)=f, pds(2)=fx, pds(3)=fy, pds(4)=fxy, pds(5)=fxx, pds(6)=fyy
        call seva2d(bkx,lkx,bky,lky,c,xxout(n),yyout(n),pds,ier,n333)
        bpooz(n)=pds(2)/xxout(n)
        bpoo(n)=sqrt(bpooz(n)**2+(pds(3)/xxout(n))**2)
      enddo

      sumip=0.
      do i=2,kfound
        delx=xxout(i)-xxout(i-1)
        dely=yyout(i)-yyout(i-1)
        dell=sqrt(delx**2+dely**2)
        abpol=(bpoo(i-1)+bpoo(i))/2.0
        sumip=sumip+abpol*dell
      enddo
      sumip=sumip/tmu/twopi
      if(kaxis.le.0) go to 1000
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
          if((x(i).lt.rmin).or.(x(i).gt.rmax)) cycle
          if((y(j).lt.zmin).or.(y(j).gt.zmax)) cycle
          if(psipsi(kk).le.psimx.and.negcur.eq.0) cycle
          if(psipsi(kk).ge.psimx.and.negcur.eq.1) cycle
          psimx=psipsi(kk)
          xax=x(i)
          yax=y(j)
        enddo
      enddo

      xs(1)=xax
      ys(1)=yax
      ps=psimx

#ifdef DEBUG_PLTS 
      ! for debugging
      write(strtmp,'(a,i0.2,a,i0.2,a)') 'debug-surf',jtime,'-', &
                                        ifit,'.txt'
      open(unit=99,file=trim(strtmp),status='replace')
      do iyplt = 1,nz
        do ixplt = 1,nx
          call seva2d(bkx,lkx,bky,lky,c,x(ixplt),y(iyplt),pds,ier,n666)
          write(99,'(3(1x,1pe12.5))') x(ixplt),y(iyplt),pds(1)
        end do
      end do
      close(unit=99)
      write(strtmp,'(a,i0.2,a,i0.2,a)') 'debug-conv',jtime,'-', &
                                        ifit,'.txt'
      open(unit=99,file=trim(strtmp),status='replace')
#endif

      if (ifindopt==2) then
        xaxold = xax
        yaxold = yax
        pdsold = pds
        signcur = -1.0
        if(negcur.eq.0) signcur = 1.0
      end if

      errtmp = 0.0
      do j=1,niter
        ! pds(1)=f, pds(2)=fx, pds(3)=fy, pds(4)=fxy, pds(5)=fxx, pds(6)=fyy
        call seva2d(bkx,lkx,bky,lky,c,xax,yax,pds,ier,n666)
#ifdef DEBUG_PLTS 
        write(99,'(3(1x,1pe12.5))') xax,yax,pds(1) ! for debugging
#endif

        search_method: if (ifindopt==2) then
          ! Gradient Ascent Method - better for sharp peaks
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
          if(errtmp.lt.1.0e-12_dp) go to 310

        else search_method
          ! Original Newton's Method for optimization, xn+1 = xn - f'/f''
          det=pds(5)*pds(6)-pds(4)*pds(4)
          if (abs(det).lt.1.0e-15_dp) then
            kerror = 1
            call errctrl_msg('findax', &
              'Newtons method to find magnetic axis has det=0')
            if(iand(iout,1).ne.0) write (nout,5000) xax,yax
            return
          endif
          xerr=(-pds(2)*pds(6)+pds(4)*pds(3))/det
          yerr=(-pds(5)*pds(3)+pds(2)*pds(4))/det
          xax=xax+orelax*xerr
          yax=yax+orelax*yerr
          errtmp = xerr*xerr+yerr*yerr
          if (xax<x(1) .or. xax>x(nx) .or. yax<y(1) .or. yax>y(nz)) then
            kerror = 1
            call errctrl_msg('findax', &
              'Newtons method to find separatrix point is off grid')
            if(iand(iout,1).ne.0) write (nout,5000) xax,yax
            return
          endif
          if ((abs(pds(2)).lt.1.0e-06_dp).and.(abs(pds(3)).lt.1.0e-06_dp)) go to 310
          if (errtmp.lt.1.0e-12_dp) go to 310
        endif search_method
      enddo
      if (errtmp.gt.1.0e-6_dp) then
        call errctrl_msg('findax', &
         'Iterative method to find magnetic axis reached max iterations',2)
      endif
      !if (iand(iout,1).ne.0) write (nout,5000) xax,yax
      xax=xs(1)
      yax=ys(1)
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
      if(emaxis.gt.0.0) emaxis=sqrt(emaxis)
      if(emaxis.le.0.0) emaxis=1.3_dp
 1000 continue
#ifdef DEBUG_PLTS 
      close(unit=99) ! for debugging
#endif
      delrmax1=0.40_dp
      delrmax2=0.40_dp
      sifsep=-1.e10_dp
      sissep=-1.e10_dp
      rfsep=-89.0
      zfsep=-89.0
      rssep=-89.0
      zssep=-89.0
      if(abs(dpsipsi).le.0.5_dp*psitol) return
!----------------------------------------------------------------------
!--   find the separatrix                                            --
!--   relaxed criteria for searching, 02/23/90                        --
!----------------------------------------------------------------------
      bpols=bpoo(1)
      ns=1
      do n=2,kfound
        if(bpoo(n).ge.bpols) cycle
        bpols=bpoo(n)
        ns=n
      enddo
      xs(1)=xxout(ns)
      ys(1)=yyout(ns)

      errtmp = 0.0
      do j=1,niter
        if (xs(1).le.x(2) .or. xs(1).ge.x(nx-1) .or. &
            ys(1).le.y(2) .or. ys(1).ge.y(nz-1)) then
          kerror = 1
          call errctrl_msg('findax','1st separatrix point is off grid')
          if(iand(iout,1).ne.0) write (nout,5020) xs,ys
          return
        endif
        call seva2d(bkx,lkx,bky,lky,c,xs(1),ys(1),pds,ier,n666)
        det=pds(5)*pds(6)-pds(4)*pds(4)
        if (abs(det).lt.1.0e-15_dp) then
          kerror = 1
          call errctrl_msg('findax','1st separatrix has det=0')
          if(iand(iout,1).ne.0) write (nout,5020) xs,ys
          return
        endif
        xerr=(-pds(2)*pds(6)+pds(4)*pds(3))/det
        yerr=(-pds(5)*pds(3)+pds(2)*pds(4))/det
        xs(1)=xs(1)+orelax*xerr
        ys(1)=ys(1)+orelax*yerr
        errtmp = xerr*xerr+yerr*yerr
        if(errtmp.lt.1.0e-12_dp*100.0) exit
      enddo
      if (errtmp.gt.1.0e-6_dp*100.0) then
       call errctrl_msg('findax', &
        'Iterative method to find 1st separatrix reached max iterations',2)
      endif
!-----------------------------------------------------------------------
!--   found x separatrix point, check to see if inside vessel         --
!-----------------------------------------------------------------------
      call zlim(zeross,n111,n111,limtrv,xlimv,ylimv,xs,ys,limfagv)
      if (zeross(1).le.0.1_dp) then
        kerror = 1
        call errctrl_msg('findax', &
          'Separatrix point is not inside vessel, zeross.le.0.1')
        return
      endif
      xseps(1)=xs(1)*100.
      yseps(1)=ys(1)*100.
!-----------------------------------------------------------------------
!--   consider x separatrix point on surface if psi/dpsi/dR < 0.004 a --
!-----------------------------------------------------------------------
      anow=(rmax-rmin)*0.5_dp
      znow=0.5_dp*(zmin+zmax)
      relpsi=abs((pds(1)-psiout))
      call seva2d(bkx,lkx,bky,lky,c,rmax,znow,pdss,ier,n333)
      delrmax1=relpsi/abs(pdss(2))
      relpsi=relpsi/abs((psimx-psiout))
      if (delrmax1.gt.0.004_dp*anow) then
        !kerror = 1
        call errctrl_msg('findax', &
         'Separatrix point is not on surface, delrmax1.gt.0.004_dp*anow',2)
        return
      endif
      sifsep=pds(1)
      rfsep=xs(1)
      zfsep=ys(1)
      psiout=pds(1)
      xxout(ns)=xs(1)
      yyout(ns)=ys(1)
      xxout(ns-1)=0.5_dp*(xxout(ns)+xxout(ns-2))
      yyout(ns-1)=0.5_dp*(yyout(ns)+yyout(ns-2))
      xxout(ns+1)=0.5_dp*(xxout(ns)+xxout(ns+2))
      yyout(ns+1)=0.5_dp*(yyout(ns)+yyout(ns+2))
      loc=minloc(xxout(1:kfound),1)
      rmin=xxout(loc)
      zrmin=yyout(loc)
      loc=maxloc(xxout(1:kfound),1)
      rmax=xxout(loc)
      zrmax=yyout(loc)
      loc=minloc(yyout(1:kfound),1)
      zmin=yyout(loc)
      rzmin=xxout(loc)
      loc=maxloc(yyout(1:kfound),1)
      zmax=yyout(loc)
      rzmax=xxout(loc)
!----------------------------------------------------------------
!--   find tracing points                                      --
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
!-----------------------------------------------------------------------
!--   find possible second separatrix                                 --
!-----------------------------------------------------------------------
      bpmins=10.
      ns=-1
      do i=2,kfound
        if ((ys(1)-znow)*(yyout(i)-znow).lt.0.0.and.bpoo(i).lt.bpmins) then
          bpmins=bpoo(i)
          ns=i
        endif
      enddo
      if (ns.eq.-1) then
        kerror = 1
        call errctrl_msg('findax','2nd separatrix not found, ns.eq.-1')
        return
      endif
      xs(1)=xxout(ns)
      ys(1)=yyout(ns)
!
      errtmp = 0.0
      do j=1,niter
        if (xs(1).le.x(2) .or. xs(1).ge.x(nx-1) .or. ys(1).le.y(2) .or. ys(1).ge.y(nz-1)) then
          !kerror = 1
          call errctrl_msg('findax','2nd separatrix point is off grid',2)
          if(iand(iout,1).ne.0) write (nout,5025) xs,ys
          return
        endif
        call seva2d(bkx,lkx,bky,lky,c,xs(1),ys(1),pds,ier,n666)
        det=pds(5)*pds(6)-pds(4)*pds(4)
        if (abs(det).lt.1.0e-15_dp) then
          kerror = 1
          call errctrl_msg('findax','2nd separatrix has det=0')
          if(iand(iout,1).ne.0) write (nout,5025) xs,ys
          return
        endif
        xerr=(-pds(2)*pds(6)+pds(4)*pds(3))/det
        yerr=(-pds(5)*pds(3)+pds(2)*pds(4))/det
        xs(1)=xs(1)+orelax*xerr
        ys(1)=ys(1)+orelax*yerr
        errtmp = xerr*xerr+yerr*yerr
        if(errtmp.lt.1.0e-12_dp*100.0) exit
      enddo
      if (errtmp.gt.1.0e-6_dp*100.0) then
       call errctrl_msg('findax', &
        'Iterative method to find 2nd separatrix reached max iterations',2)
      endif
!-----------------------------------------------------------------------
!--   make sure 2nd seperatrix inside vessel                          --
!-----------------------------------------------------------------------
      call zlim(zeross,n111,n111,limtrv,xlimv,ylimv,xs,ys,limfagv)
      ! Previously this check was allowed to return without error and no message.
      ! It occurs a lot, so no error and supress message for now.
      if (zeross(1).le.0.1_dp) then
        !kerror = 1
        call errctrl_msg('findax', &
          '2nd separatrix point is not inside vessel, zeross.le.0.1',2)
        return
      endif
      if (abs(ys(1)*100.-yseps(1)).lt.2.0*anow) then
        !kerror = 1
        call errctrl_msg('findax','2nd seperatrix too far away',2)
        return
      endif
      ! If 2nd separatrix errors out (returns) above, the following variables
      ! (xseps=-999, rssep=-89, sissep=-1.0e-10, delrmax2=0.4) have default values that
      ! are used elsewhere.
      xseps(2)=xs(1)*100.
      yseps(2)=ys(1)*100.
      rssep=xs(1)
      zssep=ys(1)
      sissep=pds(1)
!-----------------------------------------------------------------------
!--   consider x point on surface if psi/dpsi/dR < 0.004 a            --
!-----------------------------------------------------------------------
      relpsi=abs((pds(1)-psiout))
      delrmax2=relpsi/abs(pdss(2))
      relpsi=relpsi/abs((psimx-psiout))
      if (delrmax2.gt.0.004_dp*anow) then
        !kerror = 1
        !call errctrl_msg('findax','2nd separatrix point is not on surface, delrmax2.gt.0.004_dp*anow',2)
        return
      endif

      xxout(ns)=xs(1)
      yyout(ns)=ys(1)
      xxout(ns-1)=0.5_dp*(xxout(ns)+xxout(ns-2))
      yyout(ns-1)=0.5_dp*(yyout(ns)+yyout(ns-2))
      xxout(ns+1)=0.5_dp*(xxout(ns)+xxout(ns+2))
      yyout(ns+1)=0.5_dp*(yyout(ns)+yyout(ns+2))
      loc=minloc(xxout(1:kfound),1)
      rmin=xxout(loc)
      zrmin=yyout(loc)
      loc=maxloc(xxout(1:kfound),1)
      rmax=xxout(loc)
      zrmax=yyout(loc)
      loc=minloc(yyout(1:kfound),1)
      zmin=yyout(loc)
      rzmin=xxout(loc)
      loc=maxloc(yyout(1:kfound),1)
      zmax=yyout(loc)
      rzmax=xxout(loc)

      return

 5000 format (/,1x,'no convergence to magnetic axis, rax, yax = ', &
              2(1x,e10.3))
 5020 format (/,1x,'no convergence to separatrix, rs, ys = ', &
              2(1x,e10.3))
 5025 format (/,1x,'no convergence to 2nd septrx, rs, ys = ', &
              2(1x,e10.3))
      end subroutine findax

!*********************************************************************
!>    This subroutine does ...
!!       
!!    @param x1 :
!!    @param y1 :
!!    @param x2 :
!!    @param y2 :
!!    @param f1 :
!!    @param f2 :
!!    @param f3 :
!!    @param f4 :
!!    @param x :
!!    @param y :
!!    @param carea : cell area
!!    @param psivl : value of psi at boundary
!*********************************************************************
      subroutine fqlin(x1,y1,x2,y2,f1,f2,f3,f4,x,y,carea,psivl)
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      a1=(x2-x)*(y2-y)
      a2=(x-x1)*(y2-y)
      a3=(x-x1)*(y-y1)
      a4=(x2-x)*(y-y1)
      psivl=(a1*f1+a2*f2+a3*f3+a4*f4)/carea
      return
      end subroutine fqlin

!*********************************************************************
!>    maxpsi finds the largest psi value along the line 
!!    segment y=alpha*x+beta joining the two points
!!    (xl1,yl1) and (xl2,yl2).
!!
!!    @param xl1 : R of first point in cell 
!!    @param yl1 : Z of first point in cell
!!    @param xl2 : R of second point in cell
!!    @param yl2 : Z of second point in cell
!!    @param x1 : R of left side of cell
!!    @param y1 : Z of bottom of cell
!!    @param x2 : R of right side of cell
!!    @param y2 : Z of top of cell
!!    @param f1 :
!!    @param f2 :
!!    @param f3 :
!!    @param f4 :
!!    @param carea : cell area
!!    @param psimax : largest value of psi on line segment
!!    @param xtry :
!!    @param ytry :
!!    @param nerr : 0 if the max psi is successfully computed
!!                  3 if the line segment is not within the cell
!!********************************************************************
      subroutine maxpsi(xl1,yl1,xl2,yl2,x1,y1,x2,y2,f1,f2,f3,f4, &
                        carea,psimax,xtry,ytry,nerr)
      use error_control
      implicit integer*4 (i-n), real*8 (a-h, o-z)
!
      nerr=0
      dx=xl2-xl1
      if (dx.ne.0.0) then
        dy=yl2-yl1
        if (dy.ne.0.0) then
          c=f1+f3-(f2+f4)
          if (c.ne.0.0) then
            a=y2*(f2-f1)+y1*(f4-f3)
            b=x1*(f2-f3)+x2*(f4-f1)
            alpha=dy/dx
            secder=2.*alpha*c
            if (secder.le.0.0) then
              beta=yl1-alpha*xl1
              xcrit=-(b*alpha+c*beta+a)/secder
              if ((xcrit.le.xl2).and.(xcrit.ge.xl1)) then
                ycrit=alpha*xcrit+beta
                if ((ycrit.lt.y1).or.(ycrit.gt.y2)) then
                  nerr = 3
                  call errctrl_msg('maxpsi','ycrit is out of bounds')
                  return
                end if
                xl2=xcrit
                yl2=ycrit
                psip1=-1.0e+35_dp
                go to 110
              end if
            end if
          end if
        end if
      end if
      call fqlin(x1,y1,x2,y2,f1,f2,f3,f4,xl1,yl1,carea,psip1)
  110 call fqlin(x1,y1,x2,y2,f1,f2,f3,f4,xl2,yl2,carea,psip2)
      psimax=max(psip1,psip2)
      if (psimax.eq.psip1) then
        xtry=xl1
        ytry=yl1
      else
        xtch=xl2
        ytch=yl2
      end if
      return
      end subroutine maxpsi

!*********************************************************************
!!
!>    minmax finds minimum and maximum value of psi in a cell.
!! 
!!    @param psi : poloidal flux on nw x nh grid 
!!    @param nwh :
!!    @param nh  :
!!    @param kn  :
!!    @param psivl  :
!!    @param iflag :
!!    @param knew :
!!                                                                  
!*******************************************************************
      subroutine minmax(psi,nwh,nh,kn,psivl,iflag,knew)
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension psi(nwh)
      iflag=0
      fmin=min(psi(kn),psi(kn+1),psi(kn+nh),psi(kn+nh+1))
      fmax=max(psi(kn),psi(kn+1),psi(kn+nh),psi(kn+nh+1))
      if ((psivl.lt.fmin).or.(psivl.gt.fmax)) iflag=1
      knew=kn
      return
      end subroutine minmax

!*********************************************************************
!!                                                                  
!>    order puts the points (xp,yp) in increasing order of yp. 
!!                                                                 
!!    @param xp :
!!    @param yp :
!!    @param np :                                                              
!*********************************************************************
      subroutine order(xp,yp,np)
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension xp(*),yp(*)
      nptr=np
   80 is=0
      nptr=nptr-1
      do k=1,nptr
      if (yp(k+1).ge.yp(k)) cycle
      is=1
      xs=xp(k+1)
      ys=yp(k+1)
      xp(k+1)=xp(k)
      yp(k+1)=yp(k)
      xp(k)=xs
      yp(k)=ys
      end do
      if (is.eq.1) go to 80
      return
      end subroutine order

!*********************************************************************
!!  
!>    packps orders the points (xp,yp) in sequencial order.
!!
!!    @param xp :
!!    @param yp :
!!    @param np :
!!    @param rm :
!!    @param zm :
!!    @param kadd :                                                   
!*********************************************************************
      subroutine packps(xp,yp,np,rm,zm,kadd)
      use global_constants,only: pi
      use commonblocks,only: cjrf,wxin,wyin,wxout,wyout
      include 'eparm.inc'
      implicit none
      integer*4, intent(inout) :: np
      integer*4, intent(in) :: kadd
      real*8, dimension(npoint), intent(inout) :: xp,yp
      real*8, intent(in) :: rm,zm
      integer*4 i,is,k,kin,kk,kout,mm,nptr
      real*8 rcut,rymax,rymin,slope,ws,xs,ys,ymax,ymin
      integer*4, parameter :: iflag=2
!
      if (iflag.ne.2) then
        ymin=yp(1)
        ymax=yp(1)
        rymin=xp(1)
        rymax=xp(1)
        do i=2,np
          if(yp(i).lt.ymin) rymin=xp(i)
          if(yp(i).gt.ymax) rymax=xp(i)
          ymin=min(ymin,yp(i))
          ymax=max(ymax,yp(i))
        end do
        slope=(rymax-rymin)/(ymax-ymin)
        kin=0
        kout=0
        do i=1,np
          rcut=rymin+(yp(i)-ymin)*slope
          if (xp(i).lt.rcut) then
            kin=kin+1
            wxin(kin)=xp(i)
            wyin(kin)=yp(i)
          else
            kout=kout+1
            wxout(kout)=xp(i)
            wyout(kout)=yp(i)
          end if
        end do
        call order(wxin,wyin,kin)
        call order(wxout,wyout,kout)
        do k=1,kin
          xp(k)=wxin(k)
          yp(k)=wyin(k)
        end do
        do k=1,kout
          kk=k+kin
          mm=kout-k+1
          xp(kk)=wxout(mm)
          yp(kk)=wyout(mm)
        end do
        if(kadd.eq.0) return
        np=kk+1
        xp(np)=xp(1)
        yp(np)=yp(1)
        return
      end if
!
      do i=1,np
        wxin(i)=atan2((yp(i)-zm)/180.0*pi,(rm-xp(i))/180.0*pi)
        if(wxin(i).lt.0.0) wxin(i)=wxin(i)+360.0
      end do
      nptr=np
      is=1
      do while (is.eq.1)
        is=0
        nptr=nptr-1
        do k=1,nptr
          if(wxin(k+1).ge.wxin(k)) cycle
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
        end do
      end do
      if(kadd.eq.0) return
      np=np+1
      xp(np)=xp(1)
      yp(np)=yp(1)
      return
      end subroutine packps

!*********************************************************************
!!                                          
!>    QFIT is a quadratic fitter from three points
!!
!!    @param k  : flag for method type
!!    @param x1 : x value of point 1
!!    @param x2 : x value of point 2
!!    @param x3 : x value of point 3
!!    @param y1 : y value of point 1
!!    @param y2 : y value of point 2
!!    @param y3 : y value of point 3
!!    @param x :
!!    @param y :
!!    @param yp :
!!    @param ierr : error flag 
!*********************************************************************
      subroutine qfit(k,x1,x2,x3,y1,y2,y3,x,y,yp,ierr)
      use error_control
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      ierr=0
      if (x1.eq.x2.or.x2.eq.x3.or.x1.eq.x3) then
        ierr=1
        call errctrl_msg('qfit', &
          'at least one pair of x-coordinates are the same')
        return
      end if
      alpha=y1/((x1-x2)*(x1-x3))
      beta=y2/((x2-x1)*(x2-x3))
      gamma=y3/((x3-x1)*(x3-x2))
!
      a=alpha+beta+gamma
      b=-(gamma+beta)*x1-(gamma+alpha)*x2-(alpha+beta)*x3
      c=x2*x3*alpha+x1*x3*beta+x1*x2*gamma
!
      select case (k)
      case default ! 1
        y=a*x*x+b*x+c
      case (2)
        rad=sqrt(b*b-4.0*a*(c-y))
        root1=(-b+rad)*.5_dp/a
        root2=(-b-rad)*.5_dp/a
        t1=(root1-x1)*(x3-root1)
        t2=(root2-x1)*(x3-root2)
        zero=-x1*1.0e-7_dp
        if (t1.ge.zero) then
          if (t2.ge.zero) then
            x=min(root1,root2)
          else
            x=root1
          end if
        else if (t2.ge.zero) then
          x=root2
        else
          ierr=1
          call errctrl_msg('qfit','t1<0 or t2<0')
          return
        end if
      case (3)
        x=a
        y=b
        yp=c
        return
      end select
      yp=2.0*a*x+b
      return
      end subroutine qfit

!**********************************************************************
!!
!>    surfac generates a contour of constant psi of value
!!    siwant. 
!!
!!    @param siwant : value of psi to contour
!!    @param psi : 2d grid of psi
!!    @param nw : row dimension of psi 
!!    @param nh : column dimension of psi
!!    @param rgrid : R values on grid points
!!    @param zgrid : Z values on grid points
!!    @param xout : output R contour values
!!    @param yout : output Z contour values
!!    @param nfound : number of contour points
!!    @param npoint : maximum number of contour points
!!    @param drgrid : R grid spacing
!!    @param dzgrid : Z grid spacing 
!!    @param xmin : minimum R to contour
!!    @param xmax : maximum R to contour
!!    @param ymin : minimum Z to contour
!!    @param ymax : minimum Z to contour
!!    @param ipack : if >0, put contour in sequential order
!!    @param rmaxis : R value of magnetic axis
!!    @param zmaxis : Z value of magnetic axis
!!    @param negcur : negative current flag
!!    @param kerror : error flag
!!    @param err_type : type of error to be sent to errctrl_msg
!*********************************************************************
      subroutine surfac(siwant,psi,nw,nh,rgrid,zgrid,xout,yout, &
                        nfound,npoint,drgrid,dzgrid,xmin, &
                        xmax,ymin,ymax,ipack,rmaxis,zmaxis,negcur, &
                        kerror,err_type)
      use error_control
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      integer*4 :: err_type
      dimension xout(npoint),yout(npoint),psi(nw*nh),rgrid(nw),zgrid(nh)
      parameter (n111=1)

      kerror = 0
      if (negcur.eq.0) then
        curneg=1.
      else
        curneg=-1.
      end if
      nfound=0
      do i=2,nw-1
        do j=1,nh
          if (rgrid(i).lt.xmin) cycle
          if (rgrid(i).gt.xmax) cycle
          if (zgrid(j).lt.ymin) cycle
          if (zgrid(j).gt.ymax) cycle
          kk1=(i-1)*nh+j
          df1=siwant-psi(kk1)
          if (df1*curneg.lt.0.0.and.rgrid(i-1).lt.xmin) then
            kk2x=kk1
            df2x=df1
            kk1x=kk1-nh
            df1x=siwant-psi(kk1x)
            if (df1x*df2x.le.0.0) then
              if (nfound+1.gt.npoint-1) cycle
              nfound=nfound+1
              xout(nfound)=rgrid(i-1)+df1x*drgrid/(psi(kk2x)-psi(kk1x))
              yout(nfound)=zgrid(j)
            end if
          end if
          kk2=i*nh+j
          df2=siwant-psi(kk2)
          if (df1*df2.gt.0.0) cycle
          if (nfound+1.gt.npoint-1) cycle
          nfound=nfound+1
          xout(nfound)=rgrid(i)+df1*drgrid/(psi(kk2)-psi(kk1))
          yout(nfound)=zgrid(j)
        end do
      end do
      do i=1,nw
        do j=2,nh-1
          if (rgrid(i).lt.xmin) cycle
          if (rgrid(i).gt.xmax) cycle
          if (zgrid(j).lt.ymin) cycle
          if (zgrid(j).gt.ymax) cycle
          kk1=(i-1)*nh+j
          df1=siwant-psi(kk1)
          if (df1*curneg.lt.0.0.and.zgrid(j-1).lt.ymin) then
            kk2x=kk1
            df2x=df1
            kk1x=kk1-1
            df1x=siwant-psi(kk1x)
            if (df1x*df2x.le.0.0) then
              if (nfound+1.gt.npoint-1) cycle
              nfound=nfound+1
              xout(nfound)=rgrid(i)
              yout(nfound)=zgrid(j-1)+df1x*dzgrid/(psi(kk2x)-psi(kk1x))
            end if
          end if
          kk2=(i-1)*nh+j+1
          df2=siwant-psi(kk2)
          if (df1*df2.gt.0.0) cycle
          if (nfound+1.gt.npoint-1) cycle
          nfound=nfound+1
          xout(nfound)=rgrid(i)
          yout(nfound)=zgrid(j)+df1*dzgrid/(psi(kk2)-psi(kk1))
        end do
      end do
      if (ipack.gt.0) call packps(xout,yout,nfound,rmaxis,zmaxis,n111)
      if (nfound.lt.3) then
        kerror = 1
        call errctrl_msg('surfac','Less than 3 contour points found',err_type)
        return
      end if
      return
      end subroutine surfac

!**********************************************************************
!>    zlim determines whether points on the (x,y) grid are
!!    inside or outside of the boundary set by the limiters.\n
!!
!!    13/07/21..........WARNING added, needs fixing!
!!
!!    @param zero : 1 if inside and 0 otherwise 
!!    @param nw : dimension of x 
!!    @param nh : dimension of y
!!    @param limitr : number of limiter points
!!    @param xlim : r coordinates of limiter
!!    @param ylim : z coordinates of limiter
!!    @param x : r grid 
!!    @param y : z grid 
!!    @param iflag :  1 convex geometry\n 
!!                    2 general geometry
!*********************************************************************
      subroutine zlim(zero,nw,nh,limitr,xlim,ylim,x,y,iflag)
      use set_kinds, only: dp
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      real*8 :: zero(nw*nh),x(nw),y(nh) 
      dimension xlim(limitr),ylim(limitr)
      logical b,c,d,inside,bold

      select case (iflag)
        case (1)
          kk = 0
          do i = 1,nw
          do j = 1,nh
            kk = kk + 1
            zero(kk) = 1.
            ncross = 0
            do k = 1,limitr-1
              if ((ylim(k).lt.y(j)) .and. (ylim(k+1).lt.y(j))) cycle
              if (x(i) .eq. xlim(k))  cycle
              t = x(i) - xlim(k)
              s = xlim(k+1) - x(i)
              if ((t*s) .lt. 0.) cycle
              di = (ylim(k+1)-ylim(k)) / (xlim(k+1)-xlim(k))
              f = ylim(k) + di*(x(i)-xlim(k))
              if (f .lt. y(j)) cycle
              ncross = ncross + 1
            end do
            mcross = .5_dp*ncross ! truncates to integer
            mcross = 2*mcross
            if (ncross .eq. mcross) zero(kk) = 0.
          end do
          end do
        case (2)
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
!--             fixed if test logic, for ge and le per Wolfe of MIT       --
!--               if (y(j).le.ylim(k).and.y(j).ge.ylim(k+1)               --
!--                   .or.y(j).ge.ylim(k).and.y(j).le.ylim(k+1)) then     --
!---------------------------------------------------------------------------
!--**WARNING:   if(abs(y(j)-ylim(k)).lt.1.e-10_dp .or. 
!                  abs(y(j)-ylim(k+1)).lt.1.e-10_dp)
!               then roundoff errors can affect whether points are inside or
!               out with this algorithm... more robust method needed
                if (y(j).le.ylim(k).and.y(j).gt.ylim(k+1) &
                  .or.y(j).ge.ylim(k).and.y(j).lt.ylim(k+1)) then
                  c=.true.
                  d=.true.
                  n=n+1
                end if
                if (c) then
                  if((y(j)-ylim(k))*(xlim(k+1)-xlim(k))- &
                    (ylim(k+1)-ylim(k))*(x(i)-xlim(k)).gt.0.) &
                    b=.not.b
                end if
                if (n.eq.2) then
                  n=0
                  if (bold.eqv.b) inside=.true.
                  bold=b
                end if
              end do
              zero(kk)=0.0
              if (inside.and.d) zero(kk)=1.0
            end do
          end do
      end select
      return
      end subroutine zlim
