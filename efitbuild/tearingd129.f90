      subroutine tearing(jtime,ktear)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**         routine computes the value of delta prime                **
!**             and stability parameter lambda for a                 **
!**             select set of modes.                                 **
!**              and eventually the neoclassical terms               **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       input:                                                     **
!**               jtime                                              **
!**               ktear=1, calcuate for the default modes            **
!**                 output data to: the tearing page                 **
!**                 the a-file, and fitout.dat.                      **
!**               ktear=2, specify a set of modes in a namelist file **
!**                 output data to: the tearing page                 **
!**                 the a-file, and fitout.dat.                      **
!**               ktear=3, specify a set of modes in a namelist file **
!**                 output data to the t-file.                       **
!**       output:                                                    **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)  C.C. Hegna and J.D. Callen,                        **
!**                    Phys. Plasmas 1 (1994) 2308.                  **
!**          (2)   neoclassical reference                            **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          2/24/96 ..........first created by T.A. Gianakon        **
!**          9/22/98 ..........T.A.Gianakon, updated for consistency **
!**                            with the standalone version.          **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: psiold,psipold,psipp,cw,wkw,copyw,bwx, &
            bwy,sifprw,bwprw,cwprw,dwprw,sfprw,sprwp,xsisii,f_a, &
            f_b,f_c,f_d,pp_a,pp_b,pp_c,pp_d,c,wk,copy,bkx,bky,chi_c, &
            chi_wk,chi_copy,chi_bkx,chi_bky
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      common/cww/lwx,lwy
      real*8,allocatable :: gpsipsi_bar(:),bdotj_bar(:),gchichi_bar(:),&
         eps_avg_a(:),eps_avg_b(:),eps_avg_c(:),eps_avg_d(:), &
         gpsipsi_b(:),gpsipsi_c(:),gpsipsi_d(:),gchichi_b(:), &
         gchichi_c(:),gchichi_d(:),chi_a0(:,:),chi_b0(:,:),chi_c0(:,:), &
         chi_d0(:,:),bdotj_b(:),bdotj_c(:),bdotj_d(:),qpsi_a(:), &
         qpsi_b(:),qpsi_c(:),qpsi_d(:),chi(:),chi_psi(:),chi_phi(:), &
         gpsipsi(:),bdotj(:),gchichi(:),chi_geom(:),eps_avg(:), &
         gchi_b(:),gchi_c(:),gchi_d(:),r_cnt(:),z_cnt(:),rho_a(:), &
         rho_b(:),rho_c(:),rho_d(:)
!
      common/cwork3/lkx,lky
      parameter(num_rs=100,num_tear=33,ntdsk=17)
      common/tear_out/psi_rsg(num_rs),psi_norm(num_rs), &
                      val_lambda(num_rs), &
                      delta_prime(num_rs),w_nc(num_rs), &
                      w_star1(num_rs),w_star2(num_rs), &
                      r_rsg(num_rs),w_sat(num_rs), &
                      m_dp(num_rs),n_dp(num_rs),num_rsd
      dimension psi_root(3),mdp(num_tear),ndp(num_tear)
      dimension chi_pds(6),pds(6)
      dimension chi_a1(5),chi_b1(5),chi_c1(5),chi_d1(5)
      namelist/tear_anal/num_dp,mdp,ndp
      namelist/tear_dump/psi_rsg,psi_norm,val_lambda,delta_prime,w_nc, &
                      w_star1,w_star2,r_rsg,w_sat, &
                      m_dp,n_dp,num_rsd

!
      ALLOCATE( gpsipsi_bar(nw),bdotj_bar(nw),gchichi_bar(nw),&
         eps_avg_a(nw),eps_avg_b(nw),eps_avg_c(nw), &
         eps_avg_d(nw),gpsipsi_b(nw),gpsipsi_c(nw), &
         gpsipsi_d(nw),gchichi_b(nw),gchichi_c(nw), &
         gchichi_d(nw),chi_a0(nh,5),chi_b0(nh,5), &
         chi_c0(nh,5),chi_d0(nh,5),bdotj_b(nw),bdotj_c(nw), &
         bdotj_d(nw),qpsi_a(nw),qpsi_b(nw),qpsi_c(nw), &
         qpsi_d(nw),chi(5*nh),chi_psi(nw),chi_phi(nh),rho_a(nw), &
         rho_b(nw),rho_c(nw),rho_d(nw),gpsipsi(npoint),bdotj(npoint), &
         gchichi(npoint),chi_geom(npoint),eps_avg(npoint), &
         gchi_b(npoint),gchi_c(npoint),gchi_d(npoint), r_cnt(npoint), &
         z_cnt(npoint))

      mw=nw
      mh=nh
      pi=-4.0*atan(-1.0)
      tmu=4.0e-07*pi
      nnn=1
      nin=25
!
!          This is a hack for the
!     default mode list
      do i=1,100
         psi_rsg(i)=0.0
         psi_norm(i)=0.0
         r_rsg(i)=0.0
         val_lambda(i)=0.0
         delta_prime(i)=0.0
         w_nc(i)=0.0
         w_star1(i)=0.0
         w_star2(i)=0.0
         m_dp(i)=0
         n_dp(i)=0
      enddo
      num_dp=3
      mdp(1)=2
      mdp(2)=3
      mdp(3)=5
      ndp(1)=1
      ndp(2)=2
      ndp(3)=4
      !
      !      Read and/or create the file
      !      mode_list which contains the
      !                namelist tear_anal
      if((ktear.eq.2).or.(ktear.eq.3))then
        !
        !           Check current directory
        !      for modelist
        open(unit=nin,status='old',file='mode_list',err=12933)
        goto 12934
      !
      !     mode_list doesn't exist in cwd
12933   continue
        !     Check one directory above
        !     cwd for modelist
        open(unit=nin,status='old',file='../mode_list',err=12935)
        goto 12934
12935   open(unit=nin,status='new',file='mode_list')
        write(nin,tear_anal)
        goto 106
12934   read (nin,tear_anal,err=11219,end=106)
106     continue
11219   close(unit=nin)
      endif
!
!
!           Determine the largest value
!      of q requested
      q_big=float(mdp(1))/float(ndp(1))
      do i=2,num_dp
         q_test=float(mdp(i))/float(ndp(i))
         if(q_test.gt.q_big)q_big=q_test
      enddo
!
!           Spline the flux surface data
      call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
!
!           formulate grids:
!      (spline requires monotonicaly
!      increasing grid)
!                   psi grid (xsisii, rgrid_psi)
!      in terms of the rgrid
!                   geometric angle grid (chi_phi)
      dsi=(psibry-simag)/float(mw-1)
      if(dsi.gt.0)then
        cur_neg=-1.
      else
        cur_neg=1.
      endif
      do i=1,mw
          xsisii(i)=cur_neg*(psibry-(i-1)*dsi)
          f_a(i)=fpol(mw+1-i)
!       (pprime negative since
!       order is reversed.)
          pp_a(i)=-cur_neg*pprime(mw+1-i)
          qpsi_a(i)=qpsi(mw+1-i)
      enddo
      d_phi=2.0*pi/float(mh-1)
      do j=1,mh
        chi_phi(j)=-pi+(j-1)*d_phi
      enddo
!
!           perform additional splines
!      required for jdotb.f
      call zpline(mw,xsisii,f_a,f_b,f_c,f_d)
      call zpline(mw,xsisii,pp_a,pp_b,pp_c,pp_d)
      call zpline(mw,xsisii,qpsi_a,qpsi_b,qpsi_c,qpsi_d)
!
!           Determine domain of q_big,
!               since  q_psi needs to be
!      inside the separatrix
      i=1
      do while(q_big.lt.qpsi_a(i))
         i=i+1
      enddo
      n_qbig=i-1
      if(n_qbig.lt.(mw-10))n_qbig=mw-10
!
!      Separatrix finder
      ilast=1
      iflg=1
      j=mh/3
      z_min=0.0
      r_min=1000.0
      do while((j.gt.2)) 
        istart=2
        kk0=(istart-1)*mh+j
        kkm=kk0-mh
        do while((psi(kkm).lt.psi(kk0)).and.(istart.lt.mw))
          kk0=(istart-1)*mh+j
          kkm=kk0-mh
          istart=istart+1
        enddo
        if(r_min.gt.rgrid(istart))then
            z_min=zgrid(j)
            r_min=rgrid(istart)
        endif
        j=j-1
      enddo
      z_min=z_min-2.0*drgrid
      if(z_min.ge.zgrid(mh/3))z_min=zgrid(1)
!
!           Determine bounds for the
!      contour searches
!      (inside the separatrix)
      call seva2d(bkx,lkx,bky,lky,c,r_min,z_min,pds,ier,1)
      call surfac(pds(1),psi,mw,mh,rgrid,zgrid,r_cnt,z_cnt,nfounc &
                    ,npoint,drgrid,dzgrid,xmin,xmax,z_min,ymax,nnn, &
                    rmaxis,zmaxis,negcur)
      xmin=r_cnt(1)
      xmax=r_cnt(1)
      ymin=z_cnt(1)
      ymax=z_cnt(1)
      do iiii=2,nfounc
        if(r_cnt(iiii).lt.xmin)xmin=r_cnt(iiii)
        if(r_cnt(iiii).gt.xmax)xmax=r_cnt(iiii)
        if(z_cnt(iiii).lt.ymin)ymin=z_cnt(iiii)
        if(z_cnt(iiii).gt.ymax)ymax=z_cnt(iiii)
      enddo
      xmin=xmin-drgrid
      xmax=xmax+drgrid
      ymin=ymin-dzgrid
      ymax=ymax+dzgrid
!
!      Find the X-point
      test_length=rgrid(mw)-rgrid(1)
      do iiii=mh/3,2,-1
        psi_test=cur_neg*xsisii(1)
        call surfac(psi_test,psi,mw,mh,rgrid,zgrid &
                    ,r_cnt,z_cnt,nfounc &
                    ,npoint,drgrid,dzgrid &
                    ,xmin,xmax,zgrid(iiii-1),zgrid(iiii),nnn &
                    ,rmaxis,zmaxis,negcur)
        if(nfounc.gt.3)then
          xmin1=r_cnt(1)
          xmax1=r_cnt(1)
          jjjjmax=1
          jjjjmin=1
          do jjjj=2,nfounc
            if(r_cnt(jjjj).gt.xmax1)then
              jjjjmax=jjjj
              xmax1=r_cnt(jjjj)
            endif
            if(r_cnt(jjjj).lt.xmin1)then
              jjjjmin=jjjj
              xmin1=r_cnt(jjjj)
            endif
          enddo
          alpha1=r_cnt(jjjjmax)-r_cnt(jjjjmin)
          if(alpha1.lt.test_length)then
            ymin=zgrid(iiii)
            test_length=alpha1
          endif
        endif
      enddo
!
!
!      Calculate toroidal flux
      do iiii=1,mw
        rhox=0.0
        do i=1,mw
        do j=1,mh
          if(zgrid(j).ge.ymin)then
            kk=(i-1)*mh+j
            psi_test=cur_neg*psi(kk)
            if(psi_test.gt.xsisii(iiii))then
              fpx=seval(mw,psi_test,xsisii,f_a,f_b,f_c,f_d)
              rhox=rhox+fpx/rgrid(i)
            endif
           endif
        enddo
        enddo
        rhox=rhox*drgrid*dzgrid
        rho_a(iiii)=(abs(rhox))**0.5
      enddo
      rho_norm=rho_a(1)
      do iiii=1,mw
        rho_a(iiii)=rho_a(iiii)/rho_norm
      enddo
      call zpline(mw,xsisii,rho_a,rho_b,rho_c,rho_d)
!      Final bounds check
      psi_test=xsisii(1)*cur_neg
      call surfac(psi_test,psi,mw,mh,rgrid,zgrid &
                    ,r_cnt,z_cnt,nfounc &
                    ,npoint,drgrid,dzgrid &
                    ,xmin,xmax,ymin,ymax,nnn &
                    ,rmaxis,zmaxis,negcur)
      xmin=r_cnt(1)
      xmax=r_cnt(1)
      ymin=z_cnt(1)
      ymax=z_cnt(1)
      write(7,*)r_cnt(1),z_cnt(1)
      do iiii=2,nfounc
        write(7,*)r_cnt(iiii),z_cnt(iiii)
        if(r_cnt(iiii).lt.xmin)xmin=r_cnt(iiii)
        if(r_cnt(iiii).gt.xmax)xmax=r_cnt(iiii)
        if(z_cnt(iiii).lt.ymin)ymin=z_cnt(iiii)
        if(z_cnt(iiii).gt.ymax)ymax=z_cnt(iiii)
      enddo
! ------------------------------------------------------------
!
!           Loop through the modes.
      num_rsd=0
      do k=1,num_dp
       qwant=float(mdp(k))/float(ndp(k))
       do i=1,mw-1
!
!      check each zone for multiple qwant
!      by solving the cubic equation for
!      multiple roots in the zone over which
!      the cubic spline is valid
        call q_search(num_root,psi_root,qwant, &
                     qpsi_a(i),qpsi_b(i),qpsi_c(i),qpsi_d(i))
!
!             Loop through the cubic roots
        do j=1,num_root
          psi_rs=xsisii(i)+psi_root(j)
!
!      Insure that the rational surface
!      is bounded by the region of
!      validity for the cubic spline.
         if((psi_rs.ge.xsisii(i)).and.(psi_rs.le.xsisii(i+1)))then
          num_rsd=num_rsd+1
          dpsi_work=psi_rs*0.001
!***********************************************************************
          do jjjj=1,5
!
!                                       Determine psi contour
           psi_work=psi_rs+dpsi_work*float(jjjj-3)
           chi_psi(jjjj)=psi_work
           psi_test=psi_work*cur_neg
           call surfac(psi_test,psi,mw,mh,rgrid,zgrid, &
                    r_cnt,z_cnt,nfounc, &
                    npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nnn, &
                    rmaxis,zmaxis,negcur)
!
!                                       Print all the contours
!          write(6,*)jjjj,psi_work
!          write(7,*)
!          do iiii=1,nfounc
!            write(7,*)r_cnt(iiii),z_cnt(iiii)
!          enddo
!
!      Functional values on the contour
!      and flux averaged quantities.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
           qwant=seval(mw,psi_work,xsisii, &
                      qpsi_a,qpsi_b,qpsi_c,qpsi_d)
           fwant=seval(mw,psi_work,xsisii, &
                      f_a,f_b,f_c,f_d)
          if(nfounc.eq.1)then
             eps_avg_a(jjjj)=0.0
             bdotj_bar(jjjj)=0.0
             gpsipsi_bar(jjjj)=0.0
             gchichi_bar(jjjj)=0.0
          else
!
!                                       Evaluate values for contour integral
           do kk=1,nfounc
             call jdotb(r_cnt(kk),z_cnt(kk),bdotj(kk), &
                        gpsipsi(kk),mw,mh,cur_neg)
           enddo
!
!                                       Evaluate flux-surface averages
           call fluxav(bdotj,r_cnt,z_cnt,nfounc,psi,rgrid,mw,zgrid,mh, &
                  bdotj_bar(jjjj),nnn ,sdlobp,sdlbp)
           call fluxav(gpsipsi,r_cnt,z_cnt,nfounc,psi,rgrid,mw,zgrid,mh, &
                  gpsipsi_bar(jjjj),nnn ,sdlobp,sdlbp)
          endif
!
!                                       Since r_cnt and z_cnt are ordered 
!      calculate the min&max magnetic fields.
          b_small=0.5*(1.0/(r_cnt(1)*r_cnt(1)) &
                    +1.0/(r_cnt(nfounc-1)*r_cnt(nfounc-1)))
          kk=1
          do while((z_cnt(kk)*z_cnt(kk+1)).gt.0.0)
            kk=kk+1
          enddo
          b_large=0.5*(1.0/(r_cnt(kk)*r_cnt(kk)) &
                    +1.0/(r_cnt(kk+1)*r_cnt(kk+1)))
          eps_avg_a(jjjj)=0.5*abs(b_small-b_large)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!           Straight field line calculator:

!                                       Determine integrand on
!                                       the flux contour
           do kk=1,nfounc
             gchichi(kk)=fwant/(qwant*r_cnt(kk)*r_cnt(kk))
             ang_test=atan((z_cnt(kk)-zmaxis)/(r_cnt(kK)-rmaxis))
             if((r_cnt(kk)-rmaxis).lt.0.0)then
               if((z_cnt(kk)-zmaxis).lt.0.0)then
                 ang_test=-pi+ang_test
               else
                 ang_test=pi+ang_test
               endif
             endif
             chi_geom(kk)=ang_test
           enddo
           chi_geom(nfounc)=chi_geom(nfounc)-2.0*pi
 
           call zpline(nfounc,chi_geom, &
                    gchichi,gchi_b,gchi_c,gchi_d)

!
!                                       Loop through the geometric
!                                       angles and integrate up
!                                       the straight field line angle.
           kk=(jjjj-1)*mh+1
           chi(kk)=-pi
           do jj=2,mh-1
              kk=(jjjj-1)*mh+jj
              call arcav(chi_phi(jj-1),chi_phi(jj), &
                  chi_geom,gchichi,gchi_b,gchi_c,gchi_d, &
                  r_cnt,z_cnt,nfounc,psi,rgrid,mw,zgrid,mh, &
                  chi(kk),nnn ,sdlobp,sdlbp)
              chi(kk)=chi(kk-1)+chi(kk)
           enddo
           kk=jjjj*mh
           chi(kk)=pi
!
!                                       Cleanup chi(kk) so that it normalizes 
!      to zero at chi_phi=0.0 by assuming a
!                                       linear accumulation of errors.
           mh_2=mh/2+1
           delta_cor=chi((jjjj-1)*mh+mh_2)/dble(mh_2-1)
           do jj=2,mh-1
             kk=(jjjj-1)*mh+jj
             chi(kk)=chi(kk)-delta_cor*dble(jj-1)
           enddo
!                                       Print stright field angles
!                                       versus geometric angle.
!          write(7,*)
!          do jj=1,mh
!            kk=(jjjj-1)*mh+jj
!            write(7,*)chi_phi(jj),chi(kk)
!          enddo
!
!      Transfer to alternate storage
           do jj=1,mh
             kk=(jjjj-1)*mh+jj
             chi_a0(jj,jjjj)=chi(kk)
           enddo
          enddo
!***********************************************************************
!     ================================================================
!
!                                       Compute required splines
      call zpline(5,chi_psi,bdotj_bar,bdotj_b,bdotj_c,bdotj_d)
      call zpline(5,chi_psi,gpsipsi_bar,gpsipsi_b,gpsipsi_c,gpsipsi_d)
      call zpline(5,chi_psi,eps_avg_a,eps_avg_b,eps_avg_c,eps_avg_d)
      do jjjj=1,5
        call zpline(mh,chi_phi,chi_a0(1,jjjj),chi_b0(1,jjjj), &
                               chi_c0(1,jjjj),chi_d0(1,jjjj))
      enddo
!     ================================================================
!
!                 Evaluate parameters at psi_rs
            psi_rsg(num_rsd)=psi_rs
            m_dp(num_rsd)=mdp(k)
            n_dp(num_rsd)=ndp(k)
            psi_test=psi_rs*cur_neg
!
!      gttwant is a messy one:
            call surfac(psi_test,psi,mw,mh, &
                    rgrid,zgrid,r_cnt,z_cnt,nfounc, &
                    npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nnn, &
                    rmaxis,zmaxis,negcur)

            do kk=1,nfounc
              chi_ang=atan((z_cnt(kk)-zmaxis)/(r_cnt(kk)-rmaxis))
              if((r_cnt(kk)-rmaxis).lt.0.0)then
                if((z_cnt(kk)-zmaxis).lt.0.0)then
                  chi_ang=-pi+chi_ang
                else
                  chi_ang=pi+chi_ang
                endif
              endif
            

              do jjjj=1,5
                chi_a1(jjjj)=seval(mh,chi_ang,chi_phi, &
                              chi_a0(1,jjjj),chi_b0(1,jjjj), &
                              chi_c0(1,jjjj),chi_d0(1,jjjj))
              enddo
              call zpline(5,chi_psi,chi_a1,chi_b1,chi_c1,chi_d1)
              chi_pds(1)=chi_a1(3)
              chi_pds(2)=speval(5,psi_rs,chi_psi, &
                                chi_a1,chi_b1,chi_c1,chi_d1)
              chi_pds(3)=speval(mh,chi_ang,chi_phi, &
                      chi_a0(1,3),chi_b0(1,3),chi_c0(1,3),chi_d0(1,3))
              call gchichi_calc(r_cnt(kk),z_cnt(kk), &
                                gchichi(kk),chi_pds,cur_neg)
            enddo

            call fluxav(gchichi,r_cnt,z_cnt,nfounc, &
                  psi,rgrid,mw,zgrid,mh, &
                  gttwant,nnn ,sdlobp,sdlbp)

!
!      Evaluate the other terms
!      from their splines.
            qpwant=speval(mw,psi_rs,xsisii, &
                      qpsi_a,qpsi_b,qpsi_c,qpsi_d)
            fwant=seval(mw,psi_rs,xsisii, &
                      f_a,f_b,f_c,f_d)
            ppwant=seval(mw,psi_rs,xsisii, &
                      pp_a,pp_b,pp_c,pp_d)
            sigwant=speval(5,psi_rs,chi_psi, &
                      bdotj_bar,bdotj_b,bdotj_c,bdotj_d)
            grrwant=seval(5,psi_rs,chi_psi, &
                      gpsipsi_bar,gpsipsi_b,gpsipsi_c,gpsipsi_d)
            eps_want=seval(5,psi_rs,chi_psi, &
                      eps_avg_a,eps_avg_b,eps_avg_c,eps_avg_d)
            if(eps_want.lt.0.0)eps_want=0.0
            psi_norm(num_rsd)=seval(mw,psi_rs,xsisii, &
                      rho_a,rho_b,rho_c,rho_d)
            r_rsg(num_rsd)=1.0/(gttwant**0.5)
!
!                 Evaluate lambda
!           write(6,*)fwant,qwant,sigwant,qpwant,grrwant,gttwant
            val_lambda(num_rsd)= fwant*qwant*sigwant &
                      /(2.0*float(mdp(k))*qpwant) &
                      /(grrwant*gttwant)**0.5
            val_lambda(num_rsd)=abs(val_lambda(num_rsd))*tmu
            if(val_lambda(num_rsd).gt.0.9999) &
               val_lambda(num_rsd)=0.9999

!
!                 Evaluate saturated island
!      half-width (m) if appropriate.
            if(val_lambda(num_rsd).gt.0.5)then
              w_sat(num_rsd)=2.04*(val_lambda(num_rsd)-0.5) &
                             /(m_dp(num_rsd)*gttwant**0.5)
            else
              w_sat(num_rsd)=0.0
            endif
!
!      Evaluate delta prime (m)
            delta_prime(num_rsd)=-2*float(mdp(k))*gttwant**0.5 &
           *pi*val_lambda(num_rsd)/tan(pi*val_lambda(num_rsd))
!
!      Evaluate neoclassical
!      width w_nc
            w_nc(num_rsd)=pi*17.2e-11*eps_want**0.5 &
                 *qwant*ppwant*(rout(jtime)**2)/(qpwant*grrwant)
!
!                     Evaluate w_star factor in:
!                 Collisional limit
            den_want=1.0e+19
            w_star1(num_rsd)=3.275e+05 &
                     *(qwant*ppwant/qpwant)**2 &
                     /(grrwant/rout(jtime)**2)**2 &
                     /den_want
!
!                  Collisionless limit
            w_star2(num_rsd)=w_star1(num_rsd)*eps_want**1.5
          endif
        enddo
       enddo
      enddo
!
      DEALLOCATE( gpsipsi_bar,bdotj_bar,gchichi_bar,eps_avg_a, &
         eps_avg_b,eps_avg_c,eps_avg_d,gpsipsi_b,gpsipsi_c, &
         gpsipsi_d,gchichi_b,gchichi_c,gchichi_d,chi_a0,chi_b0, &
         chi_c0,chi_d0,bdotj_b,bdotj_c,bdotj_d,qpsi_a,qpsi_b,qpsi_c, &
         qpsi_d,chi,chi_psi,chi_phi,rho_a,rho_b,rho_c,rho_d,gpsipsi, &
         bdotj,gchichi,chi_geom,eps_avg,gchi_b,gchi_c,gchi_d,r_cnt, &
         z_cnt)
!
      return
      end

      subroutine jdotb(r_loci,z_loci,bdotj_val,gpsipsi,mw,mh,cur_neg)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**         routine computes the local value of J dot B / B^2        **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       input:                                                     **
!**         r_loci    -->  radial grid location                      **
!**         z_loci    -->  z grid location                           **
!**         cur_neg   -->  The current is negative
!**   Spline fits of F,p' and  bicubic spline fit of psi      **
!**       output:                                                    **
!**         bdotj_val -->  Computed value of J dot B / B^2           **
!**         gpsipsi   -->  Local value of the g^psi^psi metric elmt  **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)  T.A. Gianakon                                      **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          2/23/96/..........first created                         **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: cw,wkw,copyw,bwx, &
            bwy,sifprw,bwprw,cwprw,dwprw,sfprw,sprwp,xsisii,f_a, &
            f_b,f_c,f_d,pp_a,pp_b,pp_c,pp_d,c,wk,copy,bkx,bky
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      common/cww/lwx,lwy
      dimension pds(6)
      
      common/cowrk3/lkx,lky
!           Calculate B^2
        call seva2d(bkx,lkx,bky,lky,c,r_loci,z_loci,pds,ier,n333)
      psi_test=cur_neg*pds(1)
      fwant=seval(mw,psi_test,xsisii,f_a,f_b,f_c,f_d)
      gpsipsi=pds(2)*pds(2)+pds(3)*pds(3)
      bdotb=(fwant*fwant+gpsipsi)/(r_loci*r_loci)
!
!            Calculate B dot J
      fpwant=speval(mw,psi_test,xsisii,f_a,f_b,f_c,f_d)
      ppwant=seval(mw,psi_test,xsisii,pp_a,pp_b,pp_c,pp_d)
      bdotj_val=-(fpwant/tmu+fwant*ppwant/bdotb)

      return
      end


      subroutine gchichi_calc(r_loci,z_loci,gchichi,chi_pds,cur_neg)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**         routine computes the local value of gchichi              **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       input:                                                     **
!**         r_loci    -->  radial grid location   (m)                **
!**         z_loci    -->  z grid location        (m)                **
!**         chi_pds   --> 
!**         cur_neg   -->
!**   Spline fits of F,p' and  bicubic spline fit of psi             **
!**       output:                                                    **
!**         gchichi   -->  Local value of the g^chi^chi metric elmt  **
!**                             units of (m^{-2})                    **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          2/23/96/..........first created  by T.A. Gianakon       **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky,chi_c,chi_wk,chi_copy, &
                             chi_bkx,chi_bky
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      dimension pds(6),chi_pds(6)
      
!
!           Various pre-calculated
!      common blocks
      common/cwork3/lkx,lky
!
!           First convert the geometric
!      locations to data which can be
!      used to evaluate the chi-spline
      chi_ang=atan((z_loci-zmaxis)/(r_loci-rmaxis))
      if((r_loci-rmaxis).lt.0.0)then
        if((z_loci-zmaxis).lt.0.0)then
          chi_ang=-pi+chi_ang
        else
          chi_ang=pi+chi_ang
        endif
      endif
!
!      Calculate psi, d psi /dR
!      and d psi /dZ
      call seva2d(bkx,lkx,bky,lky,c, &
                  r_loci,z_loci,pds,ier,n333)
!
!      Calculate d chi /d psi
!      and d chi /d Theta
!     call seva2d(chi_bkx,lkx,chi_bky,lky,chi_c,
!    &            psi_test,chi_ang,chi_pds,ier,n333)
          pds(1)=cur_neg*pds(1)
          pds(2)=cur_neg*pds(2)
          pds(3)=cur_neg*pds(3)

!
!           Calculate gchichi
!      chi_term1: grad R component
!      chi_term2: grad Z component
      ang_mag=(r_loci-rmaxis)**2+(z_loci-zmaxis)**2
      chi_term1=chi_pds(2)*pds(2)- &
                chi_pds(3)*(z_loci-zmaxis)/ang_mag
      chi_term2=chi_pds(2)*pds(3)+ &
                chi_pds(3)*(r_loci-rmaxis)/ang_mag
      gchichi=chi_term1*chi_term1+chi_term2*chi_term2
      return
      end


      subroutine arcav(ang_start,ang_end, &
            chi_geom,gchi_a,gchi_b,gchi_c,gchi_d, &
            x,y,n,si,rx,msx,ry,msy, &
            fave,ns,sdlobp,sdlbp)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          arcav computes an integral on a flux surface.           **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**          ang_start: starting angle  (geometric)                  **
!**          ang_end: ending angle      (geometric)                  **
!**          chi_geom: geometric angle associated with spline of     **
!**          gchi_? : function and spline to be integrated           **
!**          si    : flux array                                      **
!**          x,y,n : n (R,Z) coordinates of contour                  **
!**                                                                  **
!**     RETURN ARGUMENTS:                                            **
!**          fave  : int(dl f/Bp) / sdlobp                           **
!**          sdlbp : int(dl Bp)                                      **
!**          sdlobp: int(dl/Bp)                                      **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          02/29/96..........first modified from fluxav            **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparmdud129.f90'
      common/cowrk3/lkx,lky
       dimension pds(6)
      real*8,dimension(:),allocatable :: x,y,si,rx,ry,chi_geom, &
              gchi_a,gchi_b,gchi_c,gchi_d
!
      ALLOCATE(x(npoint),y(npoint),si(npoint),rx(npoint), &
         ry(npoint),chi_geom(npoint),gchi_a(npoint),gchi_b(npoint), &
         gchi_c(npoint),gchi_d(npoint))
!
      if (ns.ne.0) then
      call sets2d(si,c,rx,msx,bkx,lkx,ry,msy,bky,lky,wk,ier)
      endif
!
!
!          Determine the range of flux
!      surface points which span
!           the interval of the flux
!      surface arc length integral.
      k=1
      pi=-4.0*atan(-1.0)
      ifnd=1
      n_start=n-1
      n_end=n-1
      do while(k.lt.n)
        ang_test=chi_geom(k)
        if(ifnd.gt.0)then
          if(ang_test.lt.ang_end)then
            ifnd=-1
            n_end=k-1
!
!                  Is ang_start in this interval?
            if(ang_test.lt.ang_start)then
              n_start=k-1
              k=n
            endif
          endif
        else
          if(ang_test.lt.ang_start)then
            n_start=k-1
            k=n
          endif
        endif
      k=k+1
      enddo
!     write(6,*)chi_geom(n_end),ang_end,chi_geom(n_end+1)
!     write(6,*)chi_geom(n_start),ang_start,chi_geom(n_start+1)
!     write(6,*)n_start,n_end,n
!      Implies that
!       chi_geom(n_end)   >   ang_end   >   chi_geom(n_end+1)
!       chi_geom(n_start) >   ang_start >   chi_geom(n_start+1)
!       ang_end > ang_start
!------------------------------------------------------------------
!          Integral of f on a portion
!      of a flux surface
!------------------------------------------------------------------
!
!     
      if(n_end.gt.0)then
       if(n_start.lt.n)then
        if(n_start.eq.n_end)then
          xnow=0.5*(x(n_start)+x(n_start+1))
          ynow=0.5*(y(n_start)+y(n_start+1))
          fnow=0.5*(gchi_a(n_start)+gchi_a(n_start+1))
          dxnow=x(n_start+1)-x(n_start)
          dynow=y(n_start+1)-y(n_start)
          dl_2=chi_geom(n_start)-chi_geom(n_start+1)
          dl_1=(ang_end-ang_start)/dl_2
          dl=sqrt(dxnow**2+dynow**2)*dl_1
          call seva2d(bkx,lkx,bky,lky,c,xnow,ynow,pds,ier,n333)
          bpol = sqrt(pds(2)**2+pds(3)**2)/xnow
          dlbpol = dl/bpol
          fnorm = dlbpol
          fave = fnow*dlbpol
          sdlbp = dl*bpol
        else
          fave=0.0
          fnorm=0.0
          sdlbp=0.0
!
!               The first region integral
!      (ang_end to chi_geom(n_end+1))
          xnow=0.5*(x(n_end)+x(n_end+1))
          ynow=0.5*(y(n_end)+y(n_end+1))
          fnow=0.5*(gchi_a(n_end)+gchi_a(n_end+1))
          dxnow=x(n_end+1)-x(n_end)
          dynow=y(n_end+1)-y(n_end)
          dl_2=chi_geom(n_end)-chi_geom(n_end+1)
          dl_1=(ang_end-chi_geom(n_end+1))/dl_2
          dl=sqrt(dxnow**2+dynow**2)*dl_1
          call seva2d(bkx,lkx,bky,lky,c,xnow,ynow,pds,ier,n333)
          bpol = sqrt(pds(2)**2+pds(3)**2)/xnow
          dlbpol = dl/bpol
          fnorm = fnorm + dlbpol
          fave = fave + dlbpol*fnow
          sdlbp = sdlbp + dl*bpol
!
!               Integrate between
!      chi_geom(nstart+1) and
!      chi_geom(n_end-1)
          i=n_end+1
          do while(i.lt.n_start)
            xnow=0.5*(x(i)+x(i+1))
            ynow=0.5*(y(i)+y(i+1))
            fnow=0.5*(gchi_a(i)+gchi_a(i+1))
            dxnow=x(i+1)-x(i)
            dynow=y(i+1)-y(i)
            dl=sqrt(dxnow**2+dynow**2)
            call seva2d(bkx,lkx,bky,lky,c,xnow,ynow,pds,ier,n333)
            bpol = sqrt(pds(2)**2+pds(3)**2)/xnow
            dlbpol = dl/bpol
            fnorm = fnorm + dlbpol
            fave = fave + dlbpol*fnow
            sdlbp = sdlbp + dl*bpol
            i=i+1
          enddo
!
!               Integrate the last region
          xnow=0.5*(x(n_start)+x(n_start+1))
          ynow=0.5*(y(n_start)+y(n_start+1))
          fnow=0.5*(gchi_a(n_start)+gchi_a(n_start+1))
          dxnow=x(n_start+1)-x(n_start)
          dynow=y(n_start+1)-y(n_start)
          dl_2=chi_geom(n_start)-chi_geom(n_start+1)
          dl_1=(chi_geom(n_start)-ang_start)/dl_2
          dl=sqrt(dxnow**2+dynow**2)*dl_1
          call seva2d(bkx,lkx,bky,lky,c,xnow,ynow,pds,ier,n333)
          bpol = sqrt(pds(2)**2+pds(3)**2)/xnow
          dlbpol = dl/bpol
          fnorm = fnorm + dlbpol
          fave = fave + dlbpol*fnow
          sdlbp = sdlbp + dl*bpol
          sdlobp = fnorm
        endif
       endif
!     -------------
      else
!     -------------
        if(n_start.eq.n_end)then
          n_temp=n-1
          xnow=0.5*(x(n_temp)+x(n_temp+1))
          ynow=0.5*(y(n_temp)+y(n_temp+1))
          fnow=0.5*(gchi_a(n_temp)+gchi_a(n_temp+1))
          dxnow=x(n_temp+1)-x(n_temp)
          dynow=y(n_temp+1)-y(n_temp)
          dl_2=chi_geom(n_temp)-chi_geom(n_temp+1)
          dl_1=(ang_end-ang_start)/dl_2
          dl=sqrt(dxnow**2+dynow**2)*dl_1
          call seva2d(bkx,lkx,bky,lky,c,xnow,ynow,pds,ier,n333)
          bpol = sqrt(pds(2)**2+pds(3)**2)/xnow
          dlbpol = dl/bpol
          fnorm = dlbpol
          fave = fnow*dlbpol
          sdlbp = dl*bpol
        else
          fave=0.0
          fnorm=0.0
          sdlbp=0.0
!
!               The first region integral
!      (ang_end to chi_geom(n))
          n_temp=n-1
          xnow=0.5*(x(n_temp)+x(n_temp+1))
          ynow=0.5*(y(n_temp)+y(n_temp+1))
          fnow=0.5*(gchi_a(n_temp)+gchi_a(n_temp+1))
          dxnow=x(n_temp+1)-x(n_temp)
          dynow=y(n_temp+1)-y(n_temp)
          dl_2=chi_geom(n_temp)-chi_geom(n_temp+1)
          dl_1=(ang_end-chi_geom(1))/dl_2
!         dl_1=(ang_end-chi_geom(n_temp+1))/dl_2
          dl=sqrt(dxnow**2+dynow**2)*dl_1
          call seva2d(bkx,lkx,bky,lky,c,xnow,ynow,pds,ier,n333)
          bpol = sqrt(pds(2)**2+pds(3)**2)/xnow
          dlbpol = dl/bpol
          fnorm = fnorm + dlbpol
          fave = fave + dlbpol*fnow
          sdlbp = sdlbp + dl*bpol
!
!               Integrate between
!      chi_geom(nstart+1) and
!      chi_geom(n_end-1)
          i=n_end+1
          do while(i.lt.n_start)
            xnow=0.5*(x(i)+x(i+1))
            ynow=0.5*(y(i)+y(i+1))
            fnow=0.5*(gchi_a(i)+gchi_a(i+1))
            dxnow=x(i+1)-x(i)
            dynow=y(i+1)-y(i)
            dl=sqrt(dxnow**2+dynow**2)
            call seva2d(bkx,lkx,bky,lky,c,xnow,ynow,pds,ier,n333)
            bpol = sqrt(pds(2)**2+pds(3)**2)/xnow
            dlbpol = dl/bpol
            fnorm = fnorm + dlbpol
            fave = fave + dlbpol*fnow
            sdlbp = sdlbp + dl*bpol
            i=i+1
          enddo
!
!          Integrate the last region
          xnow=0.5*(x(n_start)+x(n_start+1))
          ynow=0.5*(y(n_start)+y(n_start+1))
          fnow=0.5*(gchi_a(n_start)+gchi_a(n_start+1))
          dxnow=x(n_start+1)-x(n_start)
          dynow=y(n_start+1)-y(n_start)
          dl_2=chi_geom(n_start)-chi_geom(n_start+1)
          dl_1=(chi_geom(n_start)-ang_start)/dl_2
          dl=sqrt(dxnow**2+dynow**2)*dl_1
          call seva2d(bkx,lkx,bky,lky,c,xnow,ynow,pds,ier,n333)
          bpol = sqrt(pds(2)**2+pds(3)**2)/xnow
          dlbpol = dl/bpol
          fnorm = fnorm + dlbpol
          fave = fave + dlbpol*fnow
          sdlbp = sdlbp + dl*bpol
          sdlobp = fnorm
        endif
      endif
!
      DEALLOCATE(x,y,si,rx,ry,chi_geom,gchi_a,gchi_b, &
         gchi_c,gchi_d)
!
      return
      end


      subroutine q_search(num_root,psi_root,qwant,qa,qb,qc,qd)
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
      dimension psi_root(3)
        a_zero=(qa-qwant)/qd
        a_one=qb/qd
        a_two=qc/qd

        q_cubic=a_one/3.0-a_two*a_two/9.0
        r_cubic=(a_one*a_two-3.0*a_zero)/6.0-a_two**3/27.0
        dis_cubic=q_cubic**3+r_cubic**2

        if(dis_cubic.gt.0.0)then
          num_root=1
!               Implies one real root and
!      two complex conjugates roots
          dis_cubic=(dis_cubic)**0.5
          s_one_real=r_cubic+dis_cubic
          if(s_one_real.lt.0.0)then
            s_one_real=-(-s_one_real)**(1.0/3.00)
          else
            s_one_real=(s_one_real)**(1.0/3.00)
          endif
          s_two_real=r_cubic-dis_cubic
          if(s_two_real.lt.0.0)then
            s_two_real=-(-s_two_real)**(1.0/3.00)
          else
            s_two_real=(s_two_real)**(1.0/3.00)
          endif
          psi_root(1)=-a_two/3.0+(s_one_real+s_two_real)
        else
          num_root=3
!               Implies three real roots
!
          cubic_mag=(r_cubic*r_cubic+abs(dis_cubic))**(1.0/6.0)
          dis_cubic=(abs(dis_cubic))**0.5
          cubic_angle=atan(dis_cubic/abs(r_cubic))/3.0
          if(r_cubic.gt.0)then
           real_mag=cubic_mag*cos(cubic_angle)
           unreal_mag=cubic_mag*sin(cubic_angle)
           s_one_real=real_mag
           s_one_imag=unreal_mag
           s_two_real=real_mag
           s_two_imag=-unreal_mag
          else
           cubic_angle=pi/3.0-cubic_angle
           real_mag=cubic_mag*cos(cubic_angle)
           unreal_mag=cubic_mag*sin(cubic_angle)
           s_one_real=real_mag
           s_one_imag=-unreal_mag
           s_two_real=real_mag
           s_two_imag=unreal_mag
          endif
          psi_root(1)=-a_two/3.0+2.0*s_one_real
          psi_root(2)=-a_two/3.0-s_one_real+unreal_mag*(3.0)**0.5
          psi_root(3)=-a_two/3.0-s_one_real-unreal_mag*(3.0)**0.5
        endif
      return
      end

      subroutine wtear(ktear,jtime)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**         routine prints out data from the tearing module to       **
!**           to the desired set of files.      **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       input:                                                     **
!**               jtime            **
!**               ktear=1, output data to: the tearing page      **
!**                                the a-file, and fitout.dat.      **
!**               ktear=2, output data to: the tearing page      **
!**                               the a-file, and fitout.dat.      **
!**               ktear=3, output data to the t-file.      **
!**       output:                                                    **
!**                                                                  **
!**     REFERENCES:                                                  **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          5/31/96 ..........first created by T.A. Gianakon        **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
      include 'curve2d129.inc'
      include 'env2d.inc'
      parameter(num_rs=100,num_tear=33,ntdsk=17)
      character*10 case(6)
      character let
      character*72  text
      character eqdsk*72,wform*20
      character*50 plotname
      character*10 uday
      common/tear_out/psi_rsg(num_rs),psi_norm(num_rs), &
                      val_lambda(num_rs), &
                      delta_prime(num_rs),w_nc(num_rs), &
                      w_star1(num_rs),w_star2(num_rs), &
                      r_rsg(num_rs),w_sat(num_rs), &
                      m_dp(num_rs),n_dp(num_rs),num_rsd
      namelist/tear_dump/psi_rsg,psi_norm,val_lambda,delta_prime,w_nc, &
                      w_star1,w_star2,r_rsg,w_sat, &
                      m_dp,n_dp,num_rsd
!
!     Open teqdsk as output file and then write the tear_dump namelist
      if(ktear.eq.3)then
        idum=3
        ijtime=time(jtime)
        write (case(1),fmt="('  EFITD ')")
        write (case(2),fmt="('   ',a5)") mfvers(1)
        write (case(3),fmt="(a5,'   ')") mfvers(2)
        if (ishot.le.99999) then
          write (case(4),fmt="(' # ',i5)") ishot
        else
          write (case(4),fmt="(' #',i6)") ishot
        endif
        write (case(5),fmt="('  ',i4,'ms')") ijtime
        case(6)=' '
        let = 't'
        call getfnmu(itimeu,let,ishot,ijtime,eqdsk)
        open(unit=ntdsk,file=eqdsk,status='old',err=12932)
        close(unit=ntdsk,status='delete')
12932   continue
        open(unit=ntdsk,file=eqdsk,status='new')
        write (ntdsk,fmt='(6a8,3i4)') (case(i),i=1,6),idum,nw,nh
        write (ntdsk,tear_dump)
        close(unit=ntdsk)
      elseif((ktear.eq.2).or.(ktear.eq.1))then
!
!       Write data to 2 locations: a-file and pltout.out
!
!       Write to pltout.out
        ijtime=time(jtime)
         if (itek.ge.5) then
           m_write = 1
           iunit = 35
           if (kgraph.eq.0) then
             plotname='pltout.out'
           else
             let = 'p'
             call getfnmd(let,ishot,ijtime,plotname)
             if (istore .eq. 1)  &
                  plotname = store_dir(1:lstdir)//plotname
           endif
          if (m_write .eq. 1) then
            open(unit=iunit,file = plotname, status = 'old', &
                 access = 'append')
          elseif (m_write .eq. 0) then
            open (unit=iunit, file = plotname, form = 'unformatted', &
                  status = 'old',access = 'append')
          endif
         endif

        msg=1
        call init2d
        xphy=0.0
        yphy=0.0
        xabs=1.0
        yabs=8.0
        dyabs = 0.22
        write (text,8950) (mfvers(i),i=1,2)
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 33
        xpos(msg) = xabs
        ypos(msg) = yabs
        yabs = yabs - dyabs
        ht(msg) = 0.14
!
        call date(uday)
        write (text,fmt="(' date ran = ',a10)") uday
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 22
        xpos(msg) = xabs
        ypos(msg) = yabs
        yabs = yabs - dyabs
        ht(msg) = 0.14
!
        write (text,9000) ishot
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 25
        xpos(msg) = xabs
        ypos(msg) = yabs
        yabs = yabs - dyabs
        ht(msg) = 0.14
!
        write (text,9020) itime,itimeu
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 25
        xpos(msg) = xabs
        ypos(msg) = yabs
        yabs = yabs - dyabs
        ht(msg) = 0.14
!
        write (text,fmt="(2x,'Conventional Tearing')")
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 22
        xpos(msg) = xabs
        ypos(msg) = yabs-0.07
        yabs = yabs - dyabs-0.07
        ht(msg) = 0.14
!
        write(text,1001)
  1001 format(1x,'m',2x,'n', &
         1x,'    psi     ', &
         1x,'   r        ', &
         1x,'   lambda   ', &
         1x,'   deltaprime', &
         1x,'   wsat    ')
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 70
        xpos(msg) = xabs
        ypos(msg) = yabs
        yabs = yabs - dyabs
        ht(msg) = 0.14
!
        do i=1,num_rsd
          write(text,fmt='(i2,1x,i2,5(1x,e12.5))') &
            m_dp(i),n_dp(i),psi_rsg(i), &
            r_rsg(i), &
            val_lambda(i), &
            delta_prime(i), &
            w_sat(i)
          msg = msg + 1
          note(msg) = 1
          lmes(msg) = text
          imes(msg) = 70
          xpos(msg) = xabs
          ypos(msg) = yabs
          yabs = yabs - dyabs
          ht(msg) = 0.14
        enddo
!
        write (text,fmt="(2x,'Neoclassical Tearing')")
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 22
        xpos(msg) = xabs
        ypos(msg) = yabs-0.07
        yabs = yabs - dyabs-0.07
        ht(msg) = 0.14
!
        write(text,1002)
1002    format(1x,'m',2x,'n', &
          1x,'   wnc-press', &
          1x,'   wnc-n    ', &
          1x,'   wnc-Te   ', &
          1x,'   wnc-Ti   ')
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 57
        xpos(msg) = xabs
        ypos(msg) = yabs
        yabs = yabs - dyabs
        ht(msg) = 0.14
!
        do i=1,num_rsd
          write(text,fmt='(i2,1x,i2,4(1x,e12.5))') &
                              m_dp(i),n_dp(i),w_nc(i) 
          msg = msg + 1
          note(msg) = 1
          lmes(msg) = text
          imes(msg) = 57
          xpos(msg) = xabs
          ypos(msg) = yabs
          yabs = yabs - dyabs
          ht(msg) = 0.14
        enddo
!
        write (text,fmt='(2x,20hPolarization Tearing)')
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 22
        xpos(msg) = xabs
        ypos(msg) = yabs-0.07
        yabs = yabs - dyabs-0.07
        ht(msg) = 0.14
!
        write(text,1003)
 1003   format(1x,'m',1x,'n', &
                                      1x,'   w-star1  ', &
                                      1x,'   w-star2  ')
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 30
        xpos(msg) = xabs
        ypos(msg) = yabs
        yabs = yabs - dyabs
        ht(msg) = 0.14
!
        do i=1,num_rsd
          write(text,fmt='(i2,1x,i2,4(1x,e12.5))') &
                              m_dp(i),n_dp(i),w_star1(i), &
                              w_star2(i) 
          msg = msg + 1
          note(msg) = 1
          lmes(msg) = text
          imes(msg) = 62
          xpos(msg) = xabs
          ypos(msg) = yabs
          yabs = yabs - dyabs
          ht(msg) = 0.14
        enddo
!
        ncurve=0
        xlen=0.001
        ylen=0.001
        npltlen=100
        nxlen=100
        nylen=100
        ixnon=1
        iynon=1
        xorg=0.0
        xstp=0.00001
        xmax=0.00001
        yorg=0.0
        ystp=0.00001
        ymax=0.00001
        call curve2d (ncurve, ipag, ibrdr, grce, xphy, yphy, iorel, &
          xorl, yorl,hight, bngle, bshft, ptitle, npltlen,  xtitle, &
          nxlen, ytitle, nylen, xlen, ylen, xorg, xstp, xmax, yorg, &
          ystp, ymax, iaxis, xtck, ytck, ixnon, iynon, intax, intay, &
          isaxs, sorg, stp, smax, slen, sname, nslen,xpos,ypos, &
          igridx, igridy, idash, idot, ichdsh, ichdot, &
          thcrv, sclpc, dashme, dotme, chdhme, chdtme, markme, &
          clearx, mrc, tlen, nmrk, rat, x, y, nplt, ncnct, &
          icont, nword, zmat, ix, iy, zinc, line, mode, &
          lbflg, ithk, ipri, nline, draw, &
          nshd, sx, sy, nsxy, &
          sangle, sgap, ngaps, nvec, xfm, yfm, xto, yto, ivec, &
          msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, &
          iexit)
        close(unit=iunit)
!
!       Write to the a-file
        if (keqdsk.ge.1) then
          wform='formatted'
        else
          wform='unformatted'
        endif
!
        let = 'a'
        ijtime=time(jtime)
        call getfnmu(itimeu,let,ishot,itime,eqdsk)
!----------------------------------------------------------------------
!--        If (ISTORE = 1) Then                                      --
!--        Central directory to collect EFIT results is store_dir    --
!----------------------------------------------------------------------
          if (istore .eq. 1) then
             eqdsk = store_dir(1:lstdir)//eqdsk
          endif
      open(unit=neqdsk,file=eqdsk,status='old', &
                  form=wform,access='append')
!
        if (keqdsk.ge.1) then
          write (neqdsk,1040)num_rsd
          write (neqdsk,1040)(m_dp(i),i=1,num_rsd)
          write (neqdsk,1040)(n_dp(i),i=1,num_rsd)
          write (neqdsk,1050)(val_lambda(i),i=1,num_rsd)
          write (neqdsk,1050)(delta_prime(i),i=1,num_rsd)
          write (neqdsk,1050)(w_nc(i),i=1,num_rsd)
          write (neqdsk,1050)(psi_rsg(i),i=1,num_rsd)
          write (neqdsk,1050)(r_rsg(i),i=1,num_rsd)
          write (neqdsk,1050)(w_star1(i),i=1,num_rsd)
          write (neqdsk,1050)(w_star2(i),i=1,num_rsd)
          write (neqdsk,1050)(w_sat(i),i=1,num_rsd)
  
1040      format(1x,4i4)
1050      format (1x,4e16.9)
        else
          write (neqdsk)num_rsd
          write (neqdsk)(m_dp(i),i=1,num_rsd)
          write (neqdsk)(n_dp(i),i=1,num_rsd)
          write (neqdsk)(val_lambda(i),i=1,num_rsd)
          write (neqdsk)(delta_prime(i),i=1,num_rsd)
          write (neqdsk)(w_nc(i),i=1,num_rsd)
          write (neqdsk)(psi_rsg(i),i=1,num_rsd)
          write (neqdsk)(r_rsg(i),i=1,num_rsd)
          write (neqdsk)(w_star1(i),i=1,num_rsd)
          write (neqdsk)(w_star2(i),i=1,num_rsd)
          write (neqdsk)(w_sat(i),i=1,num_rsd)
        endif
        close(neqdsk)
 
      endif
 8950 format(1x,1('*'),' TEARING STABILITY  ',2a5,1('*'))
 9000 format (' shot #   = ',i10)
 9020 format (' t(ms,us) = ',i6,1x,i6)
      return
      end
