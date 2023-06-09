#include "config.f"
!**********************************************************************
!>
!!    writes output data to HDF5 file with OMAS implementation of the
!!      IMAS format
!!
!!    @param jtime : time index
!!    @param ktime : number of times (per cpu)
!!
!**********************************************************************
      subroutine write_omas_output(jtime,ktime)
      use set_kinds
      use commonblocks,only: c,wk,bkx,bky
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      real*8 seval,speval
      integer*4, intent(in) :: jtime,ktime
      integer*4 i,j,kk,ier
      real*8 ssibry,ssimag,xdiff,xdim,zdiff,zdim
      real*8 pds(6),gaps(4),rsps(4),zsps(4)
      real*8, dimension(nw) :: vprime,dvdrho,bwork,cwork,dwork
      real*8, dimension(nw,nh) :: psirz,pcurrz,brrz,bzrz,btrz
      logical group_exists
      character fname*15
      character*10 tindex,probeind,names(4)

      ! Setup flux scalars
      if (ipmeas(jtime).gt.0.0) then
        ssimag=-simag*twopi
        ssibry=-psibry*twopi
      else
        ssimag=simag*twopi
        ssibry=psibry*twopi
      endif

      ! Unroll 2D psi (needs transpose to match python) 
      psirz=0.0
      do i=1,nw
        do j=1,nh
          kk=(i-1)*nh+j
          if (ivacum.eq.0) then
            if (ipmeas(jtime).gt.0.0) then
              psirz(j,i)=-psi(kk)*twopi
            else
              psirz(j,i)=psi(kk)*twopi
            endif
          else
            psirz(j,i)=-psi(kk)*twopi
          endif
        enddo
      enddo

      ! Unroll and un-normalize 2D jtor (needs transpose)
      pcurrz=0.0
      if (ivacum.eq.0) then
        xdim=rgrid(nw)-rgrid(1)
        zdim=zgrid(nh)-zgrid(1)
        xdiff=xdim/(nw-1)
        zdiff=zdim/(nh-1)
        do i=1,nw
          do j=1,nh
            kk=(i-1)*nh+j
                pcurrz(j,i)=pcurrt(kk)/xdiff/zdiff
          enddo
        enddo
      endif

      ! Compute 2D Br, Bz, and Btor (probably not the most efficient way...)
      call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
      do i=1,nw
        do j=1,nh
          call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n333)
          brrz(j,i)=-pds(3)/rgrid(i)
          bzrz(j,i)=pds(2)/rgrid(i)
          btrz(j,i)=seval(nw,psirz(j,i),sigrid,fpol,bbfpol, &
                          ccfpol,ddfpol)/rgrid(i)
        enddo
      enddo

      ! Compute volume derivatives
      call zpline(nw,sigrid,volp,bwork,cwork,dwork)
      do i=1,nw
        vprime(i)=speval(nw,sigrid(i),sigrid,volp, &
                         bwork,cwork,dwork)/twopi
      enddo
      call zpline(nw,rhovn,volp,bwork,cwork,dwork)
      do i=1,nw
        dvdrho(i)=speval(nw,rhovn(i),rhovn,volp,bwork,cwork,dwork)
      enddo

      ! Set name and open file
      call setfnmomas(ishot,fname)
!----------------------------------------------------------------------
!--   If (ISTORE = 1) Then                                           --
!--   Central directory to collect EFIT results is store_dir         --
!----------------------------------------------------------------------
      if(istore .eq. 1) fname = store_dir(1:lstdir)//fname
      call fch5init
      call open_h5file(trim(fname),fileid, &
                  "EFIT output file in OMAS format",rootgid,h5in,h5err)

      ! setup file only with the first timeslice
      if (jtime.eq.1 .and. rank.eq.0) then
        call make_group(rootgid,"dataset_description",fid, &
                        group_exists,h5err)
        if (group_exists) then
          write(*,*) &
            "dataset_description already exists in file, not writing"
        else
          call make_group(fid,"data_entry",sid,group_exists,h5err)
          call dump_h5(sid,"machine",device,h5in,h5err)
          call dump_h5(sid,"pulse",ishot,h5in,h5err)
          call close_group("data_entry",sid,h5err)
          call close_group("dataset_description",fid,h5err)
        endif

        ! write limiter
        call make_group(rootgid,"wall",eqid,group_exists,h5err)
        if (group_exists) then
          write(*,*) "wall already exists in file, not writing"
        else
          call make_group(eqid,"description_2d",cid,group_exists,h5err)
          call make_group(cid,"0",pid,group_exists,h5err)
          call make_group(pid,"limiter",tid,group_exists,h5err)
          call make_group(tid,"type",sid,group_exists,h5err)
          call dump_h5(sid,"description","first wall",h5in,h5err)
          call dump_h5(sid,"name","first wall",h5in,h5err)
          call dump_h5(sid,"index",0,h5in,h5err)
          call close_group("type",sid,h5err)
          call make_group(tid,"unit",sid,group_exists,h5err)
          call make_group(sid,"0",fid,group_exists,h5err)
          call make_group(fid,"outline",nid,group_exists,h5err)
          call dump_h5(nid,"r",xlim(1:limitr),h5in,h5err)
          call dump_h5(nid,"z",ylim(1:limitr),h5in,h5err)
          call close_group("outline",nid,h5err)
          call close_group("0",fid,h5err)
          call close_group("unit",sid,h5err)
          call close_group("limiter",tid,h5err)
          call close_group("0",pid,h5err)
          call close_group("description_2d",cid,h5err)
          call close_group("wall",eqid,h5err)
        endif

        ! write extra coils used in probe corrections (DIII-D only)
        if (device == "DIII-D") then
          call make_group(rootgid,"coils_non_axisymmetric",fid, &
                          group_exists,h5err)
          if (group_exists) then
            write(*,*) &
              "coils_non_asym already exists in file, not writing"
          else
            call make_group(fid,"coil",sid,group_exists,h5err)
            call make_group(sid,"0",nid,group_exists,h5err)
            call dump_h5(nid,"identifier","C79",h5in,h5err)
            call make_group(nid,"current",cid,group_exists,h5err)
            call dump_h5(cid,"data",curc79,h5in,h5err)
            call close_group("current",cid,h5err)
            call close_group("0",nid,h5err)
            call make_group(sid,"1",nid,group_exists,h5err)
            call dump_h5(nid,"identifier","C139",h5in,h5err)
            call make_group(nid,"current",cid,group_exists,h5err)
            call dump_h5(cid,"data",curc139,h5in,h5err)
            call close_group("current",cid,h5err)
            call close_group("1",nid,h5err)
            call make_group(sid,"2",nid,group_exists,h5err)
            call dump_h5(nid,"identifier","C199",h5in,h5err)
            call make_group(nid,"current",cid,group_exists,h5err)
            call dump_h5(cid,"data",curc199,h5in,h5err)
            call close_group("current",cid,h5err)
            call close_group("1",nid,h5err)
            call close_group("coil",sid,h5err)
            call close_group("coils_non_axisymmetric",fid,h5err)
          endif
        endif
        ! Setup equilibrium properties
        call open_group(rootgid,"equilibrium",eqid,h5err)
        call make_group(eqid,"ids_properties",fid,group_exists,h5err)
        if (group_exists) then
          write(*,*) &
            "ids_properties already exists in file, not writing"
        else
          call dump_h5(fid,"comment","EFIT",h5in,h5err)
          call close_group("ids_properties",fid,h5err)
        endif
        call make_group(eqid,"vacuum_toroidal_field",fid,group_exists, &
                        h5err)
        if (group_exists) then
          write(*,*)"vacuum field already exists in file, expect issues"
        else
          call dump_h5(fid,"r0",rcentr,h5in,h5err)
          call close_group("vacuum_toroidal_field",fid,h5err)
        endif
        call close_group("equilibrium",eqid,h5err)
      endif

      ! start of equilibrium write
      call open_group(rootgid,"equilibrium",eqid,h5err)
      call open_group(eqid,"time_slice",tid,h5err)
      write(tindex,"(I0)") jtime-1+rank*ktime
      call make_group(tid,trim(tindex),sid,group_exists,h5err)
      if (group_exists) then
        write(*,*) "equilibrium index ",trim(tindex), &
                   " already exists in file, not writing"
        call close_group("equilibrium",eqid,h5err)
        call close_h5file(fileid,rootgid,h5err)
        return
      endif
      call dump_h5(sid,"time",time(jtime)/1000,h5in,h5err)

      ! write 2D variables
      call make_group(sid,"profiles_2d",cid,group_exists,h5err)
      call make_group(cid,"0",nid,group_exists,h5err)
      call dump_h5(nid,"psi",psirz,h5in,h5err)
      call dump_h5(nid,"j_tor",pcurrz,h5in,h5err)
      call dump_h5(nid,"b_field_r",brrz,h5in,h5err)
      call dump_h5(nid,"b_field_z",bzrz,h5in,h5err)
      call dump_h5(nid,"b_field_tor",btrz,h5in,h5err)
      call make_group(nid,"grid_type",fid,group_exists,h5err)
      call dump_h5(fid,"index",1,h5in,h5err)
      call close_group("grid_type",fid,h5err)
      call make_group(nid,"grid",fid,group_exists,h5err)
      call dump_h5(fid,"dim1",rgrid,h5in,h5err)
      call dump_h5(fid,"dim2",zgrid,h5in,h5err)
      call close_group("grid",fid,h5err)
      call close_group("0",nid,h5err)
      call close_group("profiles_2d",cid,h5err)

      ! write boundary
      call make_group(sid,"boundary_separatrix",cid,group_exists,h5err)
      call make_group(cid,"outline",nid,group_exists,h5err)
      call dump_h5(nid,"r",rbbbs(1:nbbbs),h5in,h5err)
      call dump_h5(nid,"z",zbbbs(1:nbbbs),h5in,h5err)
      call close_group("outline",nid,h5err)
      call make_group(cid,"x_point",nid,group_exists,h5err)
      call make_group(nid,"0",fid,group_exists,h5err)
      call dump_h5(fid,"r",rseps(1,jtime)/100.,h5in,h5err)
      call dump_h5(fid,"z",zseps(1,jtime)/100.,h5in,h5err)
      call close_group("0",fid,h5err)
      call make_group(nid,"1",fid,group_exists,h5err)
      call dump_h5(fid,"r",rseps(2,jtime)/100.,h5in,h5err)
      call dump_h5(fid,"z",zseps(2,jtime)/100.,h5in,h5err)
      call close_group("1",fid,h5err)
      call close_group("x_point",nid,h5err)
      call make_group(cid,"geometric_axis",nid,group_exists,h5err)
      call dump_h5(nid,"r",rout(jtime)/100.,h5in,h5err)
      call dump_h5(nid,"z",zout(jtime)/100.,h5in,h5err)
      call close_group("geometric_axis",nid,h5err)
      call make_group(cid,"closest_wall_point",nid,group_exists,h5err)
      call dump_h5(nid,"distance",dsep(jtime)/100.,h5in,h5err)
      !call dump_h5(nid,"r",,h5in,h5err)
      !call dump_h5(nid,"z",,h5in,h5err)
      call close_group("closest_wall_point",nid,h5err)
      call make_group(cid,"gap",nid,group_exists,h5err)
      gaps(1)=gapin(jtime); names(1)="inboard"
      gaps(2)=gapout(jtime); names(2)="outboard"
      gaps(3)=gaptop(jtime); names(3)="top"
      gaps(4)=gapbot(jtime); names(4)="bottom"
      j=0
      do i=1,4
        if(gaps(i).gt.1.0e+9_dp) cycle
        if(gaps(i).le.0.1_dp) cycle
        write(probeind,"(I0)") j
        call make_group(nid,trim(probeind),fid,group_exists,h5err)
        call dump_h5(fid,"name",names(i),h5in,h5err)
        call dump_h5(fid,"value",gaps(i)/100.,h5in,h5err)
        !call dump_h5(fid,"angle",,h5in,h5err)
        !call dump_h5(fid,"r",,h5in,h5err)
        !call dump_h5(fid,"z",,h5in,h5err)
        call close_group(trim(probeind),fid,h5err)
        j=j+1
      enddo
      call close_group("gap",nid,h5err)
      if (dsep(jtime).le.0.1_dp) then
        call dump_h5(cid,"type",0,h5in,h5err) ! limited
      else
        call dump_h5(cid,"type",1,h5in,h5err) ! diverted
      endif
      call make_group(cid,"strike_point",nid,group_exists,h5err)
      rsps(1)=rvsid(jtime); zsps(1)=zvsid(jtime) 
      rsps(2)=rvsod(jtime); zsps(2)=zvsod(jtime)
      rsps(3)=rvsiu(jtime); zsps(3)=zvsiu(jtime) 
      rsps(4)=rvsou(jtime); zsps(4)=zvsou(jtime) 
      j=0
      do i=1,4
        if(rsps(i).lt.0.1_dp) cycle
        if(zsps(i).lt.-80.) cycle
        if(abs(zsps(i)).le.0.1_dp) cycle
        write(probeind,"(I0)") j
        call make_group(nid,trim(probeind),fid,group_exists,h5err)
        call dump_h5(fid,"r",rsps(i)/100.,h5in,h5err)
        call dump_h5(fid,"z",zsps(i)/100.,h5in,h5err)
        call close_group(trim(probeind),fid,h5err)
        j=j+1
      enddo
      call close_group("strike_point",nid,h5err)
      call dump_h5(cid,"triangularity",atri(jtime),h5in,h5err)
      call dump_h5(cid,"triangularity_lower",ltri(jtime),h5in,h5err)
      call dump_h5(cid,"triangularity_upper",utri(jtime),h5in,h5err)
      call dump_h5(cid,"elongation",elong(jtime),h5in,h5err)
      call dump_h5(cid,"minor_radius",aminor(jtime)/100.,h5in,h5err)
      call dump_h5(cid,"psi",ssibry,h5in,h5err)
      call close_group("boundary_separatrix",cid,h5err)

      ! write 1D variables
      call make_group(sid,"profiles_1d",nid,group_exists,h5err)
      call dump_h5(nid,"psi",-(simag-sigrid*sidif)*twopi,h5in,h5err)
      if (ipmeas(jtime).gt.0.0) then
        call dump_h5(nid,"dpressure_dpsi",-pprime/twopi,h5in,h5err)
        call dump_h5(nid,"f_df_dpsi",-ffprim/twopi,h5in,h5err)
      else
        call dump_h5(nid,"dpressure_dpsi",pprime/twopi,h5in,h5err)
        call dump_h5(nid,"f_df_dpsi",ffprim/twopi,h5in,h5err)
      endif
      call dump_h5(nid,"q",qpsi,h5in,h5err)
      call dump_h5(nid,"f",fpol,h5in,h5err)
      call dump_h5(nid,"pressure",pres,h5in,h5err)
      call dump_h5(nid,"rho_tor_norm",rhovn,h5in,h5err)
      call dump_h5(nid,"volume",volp,h5in,h5err)
      call dump_h5(nid,"dvolume_dpsi",vprime,h5in,h5err)
      call dump_h5(nid,"dvolume_drho_tor",dvdrho,h5in,h5err)

      ! TODO need to compute and add these:
      ! j_tor [flux surface averaged]
      ! elongation
      ! triangularity_upper
      ! triangularity_lower
      ! area
      ! surface [r1surf, r2surf, or r2surg??]
      ! b_field_average
      ! b_field_max
      ! b_field_min
      ! centroid (r, r_max, r_min, z)
      ! darea_dpsi
      ! darea_drho_tor
      ! dpsi_drho_tor
      ! geometric_axis (r, z)
      ! gm1
      ! gm2
      ! gm5
      ! gm8
      ! gm9
      ! phi
      ! r_inboard
      ! r_outboard
      ! rho_tor
      ! squareness_lower_inner
      ! squareness_lower_outer
      ! squareness_upper_inner
      ! squareness_upper_outer
      ! trapped_fraction

      call close_group("profiles_1d",nid,h5err)

      ! write scalars
      call make_group(sid,"global_quantities",nid,group_exists,h5err)
      call dump_h5(nid,"psi_axis",ssimag,h5in,h5err)
      call dump_h5(nid,"psi_boundary",ssibry,h5in,h5err)
      call dump_h5(nid,"ip",ipmhd(jtime),h5in,h5err)
      call dump_h5(nid,"beta_pol",betap(jtime),h5in,h5err)
      call dump_h5(nid,"beta_tor",betat(jtime)/100,h5in,h5err)
      call dump_h5(nid,"beta_normal",abs(betan(jtime))/100,h5in,h5err)
      call dump_h5(nid,"li_3",li3(jtime),h5in,h5err)
      call dump_h5(nid,"volume",volume(jtime),h5in,h5err)
      call dump_h5(nid,"area",area(jtime)/1.0e+04_dp,h5in,h5err)
      call dump_h5(nid,"surface",psurfa(jtime),h5in,h5err)
      call make_group(nid,"magnetic_axis",cid,group_exists,h5err)
      call dump_h5(cid,"b_field_tor",btaxp(jtime),h5in,h5err)
      call dump_h5(cid,"r",rmaxis,h5in,h5err)
      call dump_h5(cid,"z",zmaxis,h5in,h5err)
      call close_group("magnetic_axis",cid,h5err)
      call dump_h5(nid,"q_axis",qm(jtime),h5in,h5err)
      call dump_h5(nid,"q_95",q95(jtime),h5in,h5err)
      call make_group(nid,"q_min",cid,group_exists,h5err)
      call dump_h5(cid,"rho_tor_norm",rhoqmin,h5in,h5err) ! TODO: need to compute tor instead of pol?
      call dump_h5(cid,"value",qmin,h5in,h5err)
      call close_group("q_min",cid,h5err)
      call close_group("global_quantities",nid,h5err)

      ! write convergence metrics
      call make_group(sid,"convergence",cid,group_exists,h5err)
      call make_group(cid,"grad_shafranov_deviation_expression",nid, &
                      group_exists,h5err)
      call dump_h5(nid,"description","Maximum absolute difference over &
                   &the plasma poloidal cross-section of the poloidal &
                   &flux between the current and preceding iteration, on &
                   &fixed grid points",h5in,h5err)
      call dump_h5(nid,"index",3,h5in,h5err)
      call dump_h5(nid,"name","max_absolute_psi_residual",h5in,h5err)
      call close_group("grad_shafranov_deviation_expression",nid,h5err)
      call dump_h5(cid,"grad_shafranov_deviation_value",terror(jtime), &
                   h5in,h5err)
      call dump_h5(cid,"iterations_n",nitera,h5in,h5err)
      call close_group("boundary",cid,h5err)

      ! write constraints
      call make_group(sid,"constraints",cid,group_exists,h5err)

      ! write total plasma current
      call make_group(cid,"ip",nid,group_exists,h5err)
      call dump_h5(nid,"exact",0,h5in,h5err)
      call dump_h5(nid,"measured",ipmeas(jtime),h5in,h5err)
      call dump_h5(nid,"measured_error_upper",sigcur,h5in,h5err)
      call dump_h5(nid,"weight",fwtcur,h5in,h5err)
      call dump_h5(nid,"reconstructed",ipmhd(jtime),h5in,h5err)
      call dump_h5(nid,"chi_squared",saiip,h5in,h5err)
      call close_group("ip",nid,h5err)

      ! write e, f and a coil currents
      call make_group(cid,"pf_current",nid,group_exists,h5err)
      do i=1,nesum
        write(probeind,"(I0)") i-1
        call make_group(nid,trim(probeind),fid,group_exists,h5err)
        call dump_h5(fid,"exact",0,h5in,h5err)
        call dump_h5(fid,"measured",eccurt(jtime,i),h5in,h5err)
        call dump_h5(fid,"measured_error_upper",sigec(i),h5in,h5err)
        call dump_h5(fid,"weight",fwtec(i),h5in,h5err)
        call dump_h5(fid,"reconstructed",cecurr(i),h5in,h5err)
        call dump_h5(fid,"chi_squared",saiec(i),h5in,h5err)
        call close_group(trim(probeind),fid,h5err)
      enddo
      do i=nesum+1,nesum+nfsum
        write(probeind,"(I0)") i-1
        call make_group(nid,trim(probeind),fid,group_exists,h5err)
        call dump_h5(fid,"exact",0,h5in,h5err)
        call dump_h5(fid,"measured",fccurt(jtime,i-nesum)/turnfc(i-nesum), &
                     h5in,h5err) ! convert A-t to A
        call dump_h5(fid,"measured_error_upper", &
                     sigfc(i-nesum)/turnfc(i-nesum),h5in,h5err) ! convert A-t to A
        call dump_h5(fid,"weight",fwtfc(i-nesum),h5in,h5err) ! non-dimensional
        call dump_h5(fid,"reconstructed",brsp(i-nesum)/turnfc(i-nesum), &
                     h5in,h5err) ! convert A-t to A
        call dump_h5(fid,"chi_squared",saifc(i-nesum),h5in,h5err) ! non-dimensional
        call close_group(trim(probeind),fid,h5err)
      enddo
      do i=nesum+nfsum+1,nesum+nfsum+nacoil
        write(probeind,"(I0)") i-1
        call make_group(nid,trim(probeind),fid,group_exists,h5err)
        call dump_h5(fid,"exact",0,h5in,h5err)
        call dump_h5(fid,"measured",accurt(jtime,i-nesum-nfsum),h5in,h5err)
        !call dump_h5(fid,"measured_error_upper",,h5in,h5err)
        !call dump_h5(fid,"weight",,h5in,h5err)
        call dump_h5(fid,"reconstructed",caccurt(jtime,1-nesum-nfsum),h5in,h5err)
        !call dump_h5(fid,"chi_squared",,h5in,h5err)
        call close_group(trim(probeind),fid,h5err)
      enddo
      call close_group("pf_current",nid,h5err)

      ! write magnetic probes
      call make_group(cid,"bpol_probe",nid,group_exists,h5err)
      do i=1,magpri
        write(probeind,"(I0)") i-1
        call make_group(nid,trim(probeind),fid,group_exists,h5err)
        call dump_h5(fid,"exact",0,h5in,h5err)
        call dump_h5(fid,"measured",expmpi(jtime,i),h5in,h5err)
        call dump_h5(fid,"measured_error_upper",sigmp2(i),h5in,h5err)
        call dump_h5(fid,"weight",fwtmp2(i),h5in,h5err)
        call dump_h5(fid,"reconstructed",cmpr2(i,jtime),h5in,h5err)
        call dump_h5(fid,"chi_squared",saimpi(i),h5in,h5err)
        call close_group(trim(probeind),fid,h5err)
      enddo
      call close_group("bpol_probe",nid,h5err)

      ! write flux loops
      call make_group(cid,"flux_loop",nid,group_exists,h5err)
      do i=1,nsilop
        write(probeind,"(I0)") i-1
        call make_group(nid,trim(probeind),fid,group_exists,h5err)
        call dump_h5(fid,"exact",0,h5in,h5err)
        call dump_h5(fid,"measured",silopt(jtime,i)*twopi,h5in,h5err)
        call dump_h5(fid,"measured_error_upper",sigsi(i)*twopi,h5in,h5err)
        call dump_h5(fid,"weight",fwtsi(i),h5in,h5err)
        call dump_h5(fid,"reconstructed",csilop(i,jtime)*twopi,h5in,h5err)
        call dump_h5(fid,"chi_squared",saisil(i),h5in,h5err)
        call close_group(trim(probeind),fid,h5err)
      enddo
      call close_group("flux_loop",nid,h5err)

      ! write MSE signals
      call make_group(cid,"mse_polarisation_angle",nid,group_exists, &
                      h5err)
      do i=1,nstark
        write(probeind,"(I0)") i-1
        call make_group(nid,trim(probeind),fid,group_exists,h5err)
        call dump_h5(fid,"exact",0,h5in,h5err)
        call dump_h5(fid,"measured",atan(tangam(jtime,i)),h5in,h5err)
        call dump_h5(fid,"measured_error_upper",atan(siggam(jtime,i)),h5in, &
                     h5err)
        call dump_h5(fid,"weight",fwtgam(i),h5in,h5err)
        call dump_h5(fid,"reconstructed",atan(cmgam(i,jtime)),h5in,h5err)
        call dump_h5(fid,"chi_squared",chigam(i),h5in,h5err)
        call close_group(trim(probeind),fid,h5err)
      enddo
      call close_group("mse_polarisation_angle",nid,h5err)
  
      ! write pressure constraints
      if (npress.gt.0) then
        call make_group(cid,"pressure",nid,group_exists,h5err)
        do i=1,npress
          write(probeind,"(I0)") i-1
          call make_group(nid,trim(probeind),fid,group_exists,h5err)
          call dump_h5(fid,"exact",0,h5in,h5err)
          call make_group(fid,"position",pid,group_exists,h5err)
          if (rpress(i).le.0) then
            call dump_h5(pid,"psi",abs(rpress(i)),h5in,h5err) ! dimensionless...
          else
            call dump_h5(pid,"r",rpress(i),h5in,h5err)
            call dump_h5(pid,"z",zpress(i),h5in,h5err)
          endif
          call close_group("position",pid,h5err)
          call dump_h5(fid,"measured",premea(i),h5in,h5err)
          call dump_h5(fid,"measured_error_upper",sigpre(i),h5in,h5err)
          call dump_h5(fid,"weight",fwtpre(i),h5in,h5err)
          call dump_h5(fid,"reconstructed",pressr(i),h5in,h5err)
          call dump_h5(fid,"chi_squared",saipre2(i),h5in,h5err)
          call close_group(trim(probeind),fid,h5err)
        enddo
        call close_group("pressure",nid,h5err)
      endif

      ! write rotation pressure constraints
      if (npress.gt.0) then
        call make_group(cid,"pressure_rotational",nid,group_exists,h5err)
        do i=1,npresw
          write(probeind,"(I0)") i-1
          call make_group(nid,trim(probeind),fid,group_exists,h5err)
          call dump_h5(fid,"exact",0,h5in,h5err)
          call make_group(fid,"position",pid,group_exists,h5err)
          if (rpresw(i).le.0) then
            call dump_h5(pid,"psi",abs(rpresw(i)),h5in,h5err) ! dimensionless...
          else
            call dump_h5(pid,"r",rpresw(i),h5in,h5err)
            call dump_h5(pid,"z",zpresw(i),h5in,h5err)
          endif
          call close_group("position",pid,h5err)
          call dump_h5(fid,"measured",premew(i),h5in,h5err)
          call dump_h5(fid,"measured_error_upper",sigprw(i),h5in,h5err)
          call dump_h5(fid,"weight",fwtprw(i),h5in,h5err)
          call dump_h5(fid,"reconstructed",presw(i),h5in,h5err)
          call dump_h5(fid,"chi_squared",saiprw(i),h5in,h5err)
          call close_group(trim(probeind),fid,h5err)
        enddo
        call close_group("pressure_rotational",nid,h5err)
      endif

      ! write current constraint
      if (npress.gt.0) then
        call make_group(cid,"j_tor",nid,group_exists,h5err)
        do i=1,kzeroj
          write(probeind,"(I0)") i-1
          call make_group(nid,trim(probeind),fid,group_exists,h5err)
          call dump_h5(fid,"exact",0,h5in,h5err)
          call make_group(fid,"position",pid,group_exists,h5err)
          if (rzeroj(i).eq.0.0) then
            if (i.eq.1 .and. sizeroj(i).lt.0.0) then
              call dump_h5(pid,"psi",psiwant,h5in,h5err) ! dimensionless...
            else
              call dump_h5(pid,"psi",sizeroj(i),h5in,h5err) ! dimensionless...
            endif
          elseif (rzeroj(i).gt.0.0) then
            call dump_h5(pid,"r",rzeroj(i),h5in,h5err)
            call dump_h5(pid,"z",sizeroj(i),h5in,h5err)
          endif
          call close_group("position",pid,h5err)
          call dump_h5(fid,"measured",vzeroj(i)*ipmeas(jtime)/carea, &
                       h5in,h5err) ! doesn't exactly match IMAS convention...
          !call dump_h5(fid,"measured_error_upper",,h5in,h5err)
          !call dump_h5(fid,"weight",,h5in,h5err)
          !call dump_h5(fid,"reconstructed",,h5in,h5err)
          !call dump_h5(fid,"chi_squared",,h5in,h5err)
          call close_group(trim(probeind),fid,h5err)
        enddo
        call close_group("j_tor",nid,h5err)
      endif

      ! write diamagnetic flux
      call make_group(cid,"diamagnetic_flux",nid,group_exists,h5err)
      call dump_h5(nid,"exact",0,h5in,h5err)
      call dump_h5(nid,"measured",diamag(jtime),h5in,h5err)
      call dump_h5(nid,"measured_error_upper",sigdia(jtime),h5in,h5err)
      call dump_h5(nid,"weight",fwtdlc,h5in,h5err)
      call dump_h5(nid,"reconstructed",cdflux(jtime),h5in,h5err)
      call dump_h5(nid,"chi_squared",chidlc,h5in,h5err)
      call close_group("diamagnetic_flux",nid,h5err)

      ! write total chi squared *These are not in the IMAS standard yet!
      call dump_h5(cid,"chi_squared_total",chifin,h5in,h5err)
      call dump_h5(cid,"chi_squared_extended",chitot,h5in,h5err)

      call close_group("constraints",cid,h5err)
      call close_group(trim(tindex),sid,h5err)
      call close_group("time_slice",tid,h5err)

      ! append time and vacuum field arrays (hopefully this order matches time-slices)
      call add_h5(eqid,"time",time(jtime)/1000_dp,h5in,h5err)
      call open_group(eqid,"vacuum_toroidal_field",fid,h5err)
      call add_h5(fid,"b0",bcentr(jtime),h5in,h5err)
      call close_group("vacuum_toroidal_field",fid,h5err)

      call close_group("equilibrium",eqid,h5err)
      call close_h5file(fileid,rootgid,h5err)
      return
      end subroutine write_omas_output

!**********************************************************************
!>
!!    writes input data to HDF5 file in OMAS dictionary format
!!      note: this does not match the IMAS XML schema, but can be
!!            converted to it
!! https://gafusion.github.io/omas/auto_examples/parse_codeparameters.html
!!
!!    @param jtime : time index
!!    @param ktime : number of times (per cpu)
!!
!**********************************************************************
      subroutine write_omas_input(jtime,ktime)
      use commonblocks,only: c,wk,bkx,bky,wgridpc,rfcpc
      use set_kinds
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none

      integer*4, intent(in) :: jtime,ktime
      integer*4 nbdryss
      real*8 xltype,xltype_180,timeus,timems
      real*8, dimension(mbdry) :: rbdryss,zbdryss
      logical group_exists
      character fname*15
      character*10 tindex

      itime=time(jtime)
      timems=itime
      timeus=(time(jtime)-timems)*1000.
      itimeu=timeus
      nbdryss=nbdry
      rbdryss=rbdry
      zbdryss=zbdry
      rbdry=0.0
      zbdry=0.0
      if (fitzts.eq.'te'.and.ztserr(jtime)) then
        nbdry=1
        rbdry(1)=1.94_dp
        zbdry(1)=ztssym(jtime)+0.5_dp*ztswid(jtime)
      endif
      n1coilt=n1coil
      n1coil=n1coils
      if(device=="DIII-D" .and. ishot.gt.108281) n1coil = 0
!----------------------------------------------------------------------
!--   get xltype, xltype_180 without reading limiters                --
!----------------------------------------------------------------------
      call getlim(0,xltype,xltype_180,.false.)
!-----------------------------------------------------------------------
!--   Write file                                                    --
!-----------------------------------------------------------------------
      call setfnmomas(ishot,fname)
!----------------------------------------------------------------------
!--   If (ISTORE = 1) Then                                           --
!--   Central directory to collect EFIT results is store_dir         --
!----------------------------------------------------------------------
      if(istore .eq. 1) fname = store_dir(1:lstdir)//fname
      call fch5init
      call open_h5file(trim(fname),fileid, &
                       "EFIT input file in OMAS format",rootgid,h5in,h5err)
      call open_group(rootgid,"equilibrium",eqid,h5err)
      call open_group(eqid,"code",cid,h5err)
      call open_group(cid,"parameters",pid,h5err)
      call open_group(pid,"time_slice",tid,h5err)
      write(tindex,"(I0)") jtime-1+rank*ktime
      call make_group(tid,trim(tindex),sid,group_exists,h5err)
      if (group_exists) then
        write(*,*) "code parameters index ",trim(tindex), &
                   " already exists in file, not writing"
        call close_group("time_slice",tid,h5err)
        call close_group("parameters",pid,h5err)
        call close_group("code",cid,h5err)
        call close_group("equilibrium",eqid,h5err)
        call close_h5file(fileid,rootgid,h5err)
        return
      endif

      call make_group(sid,"in1",nid,group_exists,h5err)
      call dump_h5(nid,"ishot",ishot,h5in,h5err)
      call dump_h5(nid,"itime",itime,h5in,h5err)
      call dump_h5(nid,"plasma",ipmeas(jtime),h5in,h5err)
      call dump_h5(nid,"itek",iteks,h5in,h5err)
      call dump_h5(nid,"itrace",itrace,h5in,h5err)
      call dump_h5(nid,"nxiter",nxiter,h5in,h5err)
      call dump_h5(nid,"fwtcur",swtcur,h5in,h5err)
      call dump_h5(nid,"kffcur",kffcurs,h5in,h5err)
      call dump_h5(nid,"coils",silopt(jtime,:),h5in,h5err)
      call dump_h5(nid,"fwtsi",swtsi,h5in,h5err)
      call dump_h5(nid,"expmp2",expmpi(jtime,:),h5in,h5err)
      call dump_h5(nid,"fwtmp2",swtmp2,h5in,h5err)
      call dump_h5(nid,"kppcur",kppcurs,h5in,h5err)
      call dump_h5(nid,"mxiter",mxiters,h5in,h5err)
      call dump_h5(nid,"ierchk",ierchk,h5in,h5err)
      call dump_h5(nid,"fwtqa",fwtqa,h5in,h5err)
      call dump_h5(nid,"qemp",qemp,h5in,h5err)
      call dump_h5(nid,"error",error,h5in,h5err)
      call dump_h5(nid,"limitr",-limid,h5in,h5err)
      call dump_h5(nid,"xlim",xlim,h5in,h5err)
      call dump_h5(nid,"ylim",ylim,h5in,h5err)
      call dump_h5(nid,"serror",serror,h5in,h5err)
      call dump_h5(nid,"nbdry",nbdry,h5in,h5err)
      call dump_h5(nid,"rbdry",rbdry,h5in,h5err)
      call dump_h5(nid,"zbdry",zbdry,h5in,h5err)
      call dump_h5(nid,"psibry",psibry,h5in,h5err)
      call dump_h5(nid,"nslref",nslref,h5in,h5err)
      call dump_h5(nid,"ibunmn",ibunmn,h5in,h5err)
      call dump_h5(nid,"btor",bcentr(jtime),h5in,h5err)
      call dump_h5(nid,"psibit",psibit,h5in,h5err)
      call dump_h5(nid,"bitmpi",bitmpi,h5in,h5err)
      call dump_h5(nid,"bitip",bitip,h5in,h5err)
      call dump_h5(nid,"icurrt",icurrt,h5in,h5err)
      call dump_h5(nid,"icinit",icinit,h5in,h5err)
      call dump_h5(nid,"brsp",fccurt(jtime,:),h5in,h5err)
      call dump_h5(nid,"iweigh",iweigh,h5in,h5err)
      call dump_h5(nid,"qenp",qenp,h5in,h5err)
      call dump_h5(nid,"fwtbp",fwtbp,h5in,h5err)
      call dump_h5(nid,"relip",relip,h5in,h5err)
      call dump_h5(nid,"zelip",zelipss,h5in,h5err)
      call dump_h5(nid,"aelip",aelip,h5in,h5err)
      call dump_h5(nid,"eelip",eelip,h5in,h5err)
      call dump_h5(nid,"qvfit",qvfit,h5in,h5err)
      call dump_h5(nid,"fwtdlc",fwtdlc,h5in,h5err)
      call dump_h5(nid,"betap0",betap0,h5in,h5err)
      call dump_h5(nid,"emp",emp,h5in,h5err)
      call dump_h5(nid,"enp",enp,h5in,h5err)
      call dump_h5(nid,"iconvr",iconvr,h5in,h5err)
      call dump_h5(nid,"icprof",icprof,h5in,h5err)
      call dump_h5(nid,"nextra",nextra,h5in,h5err)
      call dump_h5(nid,"ixstrt",ixstrt,h5in,h5err)
      call dump_h5(nid,"scrape",scrape,h5in,h5err)
      call dump_h5(nid,"errmin",errmin,h5in,h5err)
      call dump_h5(nid,"rbound",rbound,h5in,h5err)
      call dump_h5(nid,"npnef",npnef,h5in,h5err)
      call dump_h5(nid,"nptef",nptef,h5in,h5err)
      call dump_h5(nid,"fwacoil",fwacoil,h5in,h5err)
      call dump_h5(nid,"itimeu",itimeu,h5in,h5err)
      call dump_h5(nid,"rcentr",rcentr,h5in,h5err)
      call dump_h5(nid,"rzero",rzero,h5in,h5err)
      call dump_h5(nid,"gammap",gammap,h5in,h5err)
      call dump_h5(nid,"cfcoil",cfcoil,h5in,h5err)
      call dump_h5(nid,"fczero",fczero,h5in,h5err)
      call dump_h5(nid,"fcsum",fcsum,h5in,h5err)
      call dump_h5(nid,"islve",islve,h5in,h5err)
      call dump_h5(nid,"icntour",icntour,h5in,h5err)
      call dump_h5(nid,"iprobe",iprobe,h5in,h5err)
      call dump_h5(nid,"salpha",salpha,h5in,h5err)
      call dump_h5(nid,"srm",srm,h5in,h5err)
      call dump_h5(nid,"sbeta",sbeta,h5in,h5err)
      call dump_h5(nid,"ifref",ifref,h5in,h5err)
      call dump_h5(nid,"isumip",isumip,h5in,h5err)
      call dump_h5(nid,"n1coil",n1coil,h5in,h5err)
      call dump_h5(nid,"ifcurr",ifcurr,h5in,h5err)
      call dump_h5(nid,"iecurr",iecurr,h5in,h5err)
      call dump_h5(nid,"ecurrt",eccurt(jtime,:),h5in,h5err)
      call dump_h5(nid,"iecoil",iecoil,h5in,h5err)
      call dump_h5(nid,"vcurrt",vcurrt,h5in,h5err)
      call dump_h5(nid,"dflux",1.0e+03_dp*diamag(jtime),h5in,h5err)
      call dump_h5(nid,"sigdlc",1.0e+03_dp*sigdia(jtime),h5in,h5err)
      call dump_h5(nid,"iplim",iplim,h5in,h5err)
      call dump_h5(nid,"kinput",kinput,h5in,h5err)
      call dump_h5(nid,"limfag",limfag,h5in,h5err)
      call dump_h5(nid,"sigprebi",sigprebi,h5in,h5err)
      call dump_h5(nid,"fwtxx",fwtxx,h5in,h5err)
      call dump_h5(nid,"kprfit",kprfit,h5in,h5err)
      call dump_h5(nid,"pressr",pressr,h5in,h5err)
      call dump_h5(nid,"rpress",rpress,h5in,h5err)
      call dump_h5(nid,"zpress",zpress,h5in,h5err)
      call dump_h5(nid,"sigpre",sigpre,h5in,h5err)
      call dump_h5(nid,"npress",npress,h5in,h5err)
      call dump_h5(nid,"tethom",tethom,h5in,h5err)
      call dump_h5(nid,"rteth",rteth,h5in,h5err)
      call dump_h5(nid,"keqdsk",keqdsk,h5in,h5err)
      call dump_h5(nid,"zteth",zteth,h5in,h5err)
      call dump_h5(nid,"sgteth",sgteth,h5in,h5err)
      call dump_h5(nid,"npteth",npteth,h5in,h5err)
      call dump_h5(nid,"tionex",tionex,h5in,h5err)
      call dump_h5(nid,"rion",rion,h5in,h5err)
      call dump_h5(nid,"zion",zion,h5in,h5err)
      call dump_h5(nid,"sigti",sigti,h5in,h5err)
      call dump_h5(nid,"nption",nption,h5in,h5err)
      call dump_h5(nid,"dnethom",dnethom,h5in,h5err)
      call dump_h5(nid,"zeffvs",zeffvs,h5in,h5err)
      call dump_h5(nid,"rneth",rneth,h5in,h5err)
      call dump_h5(nid,"zneth",zneth,h5in,h5err)
      call dump_h5(nid,"sgneth",sgneth,h5in,h5err)
      call dump_h5(nid,"npneth",npneth,h5in,h5err)
      call dump_h5(nid,"pbeam",pbeam,h5in,h5err)
      call dump_h5(nid,"sibeam",sibeam,h5in,h5err)
      call dump_h5(nid,"nbeam",nbeam,h5in,h5err)
      call dump_h5(nid,"rzeroj",rzeroj,h5in,h5err)
      call dump_h5(nid,"xalpa",xalpa,h5in,h5err)
      call dump_h5(nid,"cgama",cgama,h5in,h5err)
      call dump_h5(nid,"ivesel",ivesel,h5in,h5err)
      call dump_h5(nid,"iexcal",iexcal,h5in,h5err)
      call dump_h5(nid,"iconsi",iconsi,h5in,h5err)
      call dump_h5(nid,"fwtfc",swtfc,h5in,h5err)
      call dump_h5(nid,"xltype",xltype,h5in,h5err)
      call dump_h5(nid,"kcalpa",kcalpa,h5in,h5err)
      call dump_h5(nid,"kcgama",kcgama,h5in,h5err)
      call dump_h5(nid,"calpa",calpa,h5in,h5err)
      call dump_h5(nid,"iacoil",iacoil,h5in,h5err)
      call dump_h5(nid,"limid",limid,h5in,h5err)
      call dump_h5(nid,"irfila",irfila,h5in,h5err)
      call dump_h5(nid,"jzfila",jzfila,h5in,h5err)
      call dump_h5(nid,"vloop",vloopt(jtime),h5in,h5err)
      call dump_h5(nid,"iqplot",iqplot,h5in,h5err)
      call dump_h5(nid,"siref",psiref(jtime),h5in,h5err)
      call dump_h5(nid,"denr",denrt(jtime,:),h5in,h5err)
      call dump_h5(nid,"denv",denvt(jtime,:),h5in,h5err)
      call dump_h5(nid,"xgama",xgama,h5in,h5err)
      call dump_h5(nid,"nptionf",nptionf,h5in,h5err)
      call dump_h5(nid,"currn1",curtn1(jtime),h5in,h5err)
      call dump_h5(nid,"ifitvs",ifitvs,h5in,h5err)
      call dump_h5(nid,"bitfc",bitfc,h5in,h5err)
      call dump_h5(nid,"relax",relax,h5in,h5err)
      call dump_h5(nid,"saimin",saimin,h5in,h5err)
      call dump_h5(nid,"icutfp",icutfp,h5in,h5err)
      call dump_h5(nid,"cutip",cutip,h5in,h5err)
      call dump_h5(nid,"iavem",iavem,h5in,h5err)
      call dump_h5(nid,"pnbeam",pbinj(jtime),h5in,h5err)
      call dump_h5(nid,"xltype_180",xltype_180,h5in,h5err)
      call dump_h5(nid,"sgprmin",sgprmin,h5in,h5err)
      call dump_h5(nid,"elomin",elomin,h5in,h5err)
      call dump_h5(nid,"fcurbd",fcurbd,h5in,h5err)
      call dump_h5(nid,"pcurbd",pcurbd,h5in,h5err)
      call dump_h5(nid,"prbdry",prbdry,h5in,h5err)
      call dump_h5(nid,"ndokin",ndokin,h5in,h5err)
      call dump_h5(nid,"zlowimp",zlowimp,h5in,h5err)
      call dump_h5(nid,"kskipvs",kskipvs,h5in,h5err)
      call dump_h5(nid,"limvs",limvs,h5in,h5err)
      call dump_h5(nid,"vcurfb",vcurfb,h5in,h5err)
      call dump_h5(nid,"kpressb",kpressb,h5in,h5err)
      call dump_h5(nid,"pressbi",pressbi,h5in,h5err)
      call dump_h5(nid,"prespb",prespb,h5in,h5err)
      call dump_h5(nid,"sigppb",sigppb,h5in,h5err)
      call dump_h5(nid,"kzeroj",kzeroj,h5in,h5err)
      call dump_h5(nid,"rminvs",rminvs,h5in,h5err)
      call dump_h5(nid,"rmaxvs",rmaxvs,h5in,h5err)
      call dump_h5(nid,"errbry",errbry,h5in,h5err)
      call dump_h5(nid,"fwtpre",fwtpre,h5in,h5err)
      call dump_h5(nid,"ibtcomp",ibtcomp,h5in,h5err)
      call dump_h5(nid,"klabel",klabel,h5in,h5err)
      call dump_h5(nid,"zmaxvs",zmaxvs,h5in,h5err)
      call dump_h5(nid,"dnbeam",dnbeam,h5in,h5err)
      call dump_h5(nid,"dmass",dmass,h5in,h5err)
      call dump_h5(nid,"nmass",nmass,h5in,h5err)
      call dump_h5(nid,"condin",condin,h5in,h5err)
      call dump_h5(nid,"iaveus",iaveus,h5in,h5err)
      call dump_h5(nid,"sgtimin",sgtimin,h5in,h5err)
      call dump_h5(nid,"kwripre",kwripre,h5in,h5err)
      call dump_h5(nid,"kbound",kbound,h5in,h5err)
      call dump_h5(nid,"alphafp",alphafp,h5in,h5err)
      call dump_h5(nid,"kframe",kframe,h5in,h5err)
      call dump_h5(nid,"zbound",zbound,h5in,h5err)
      call dump_h5(nid,"vsdamp",vsdamp,h5in,h5err)
      call dump_h5(nid,"zminvs",zminvs,h5in,h5err)
      call dump_h5(nid,"saicon",saicon,h5in,h5err)
      call dump_h5(nid,"kppfnc",kppfnc,h5in,h5err)
      call dump_h5(nid,"kppknt",kppknt,h5in,h5err)
      call dump_h5(nid,"ppknt",ppknt,h5in,h5err)
      call dump_h5(nid,"pptens",pptens,h5in,h5err)
      call dump_h5(nid,"kfffnc",kfffnc,h5in,h5err)
      call dump_h5(nid,"kffknt",kffknt,h5in,h5err)
      call dump_h5(nid,"ffknt",ffknt,h5in,h5err)
      call dump_h5(nid,"fftens",fftens,h5in,h5err)
      call dump_h5(nid,"fwtbdry",fwtbdry,h5in,h5err)
      call dump_h5(nid,"kwwfnc",kwwfnc,h5in,h5err)
      call dump_h5(nid,"kwwknt",kwwknt,h5in,h5err)
      call dump_h5(nid,"wwknt",wwknt,h5in,h5err)
      call dump_h5(nid,"wwtens",wwtens,h5in,h5err)
      call dump_h5(nid,"fwtec",swtec,h5in,h5err)
      call dump_h5(nid,"fitsiref",fitsiref,h5in,h5err)
      call dump_h5(nid,"bitec",bitec,h5in,h5err)
      call dump_h5(nid,"scalepr",scalepr,h5in,h5err)
      call dump_h5(nid,"scalesir",scalesir,h5in,h5err)
      call dump_h5(nid,"ppbdry",ppbdry,h5in,h5err)
      call dump_h5(nid,"kppbdry",kppbdry,h5in,h5err)
      call dump_h5(nid,"pp2bdry",pp2bdry,h5in,h5err)
      call dump_h5(nid,"kpp2bdry",kpp2bdry,h5in,h5err)
      call dump_h5(nid,"scalea",scalea,h5in,h5err)
      call dump_h5(nid,"sigrbd",sigrbd,h5in,h5err)
      call dump_h5(nid,"sigzbd",sigzbd,h5in,h5err)
      call dump_h5(nid,"nbskip",nbskip,h5in,h5err)
      call dump_h5(nid,"ffbdry",ffbdry,h5in,h5err)
      call dump_h5(nid,"kffbdry",kffbdry,h5in,h5err)
      call dump_h5(nid,"ff2bdry",ff2bdry,h5in,h5err)
      call dump_h5(nid,"kff2bdry",kff2bdry,h5in,h5err)
      call dump_h5(nid,"errsil",errsil,h5in,h5err)
      call dump_h5(nid,"vbit",vbit,h5in,h5err)
      call dump_h5(nid,"wwbdry",wwbdry,h5in,h5err)
      call dump_h5(nid,"kwwbdry",kwwbdry,h5in,h5err)
      call dump_h5(nid,"ww2bdry",ww2bdry,h5in,h5err)
      call dump_h5(nid,"kww2bdry",kww2bdry,h5in,h5err)
      call dump_h5(nid,"f2edge",f2edge,h5in,h5err)
      call dump_h5(nid,"fe_width",fe_width,h5in,h5err)
      call dump_h5(nid,"fe_psin",fe_psin,h5in,h5err)
      call dump_h5(nid,"kedgef",kedgef,h5in,h5err)
      call dump_h5(nid,"kersil",kersil,h5in,h5err)
      call dump_h5(nid,"iout",iout,h5in,h5err)
      call dump_h5(nid,"ixray",ixray,h5in,h5err)
      call dump_h5(nid,"pedge",pedge,h5in,h5err)
      call dump_h5(nid,"kedgep",kedgep,h5in,h5err)
      call dump_h5(nid,"pe_width",pe_width,h5in,h5err)
      call dump_h5(nid,"pe_psin",pe_psin,h5in,h5err)
      call dump_h5(nid,"table_dir",table_dir,h5in,h5err)
      call dump_h5(nid,"input_dir",input_dir,h5in,h5err)
      call dump_h5(nid,"store_dir",store_dir,h5in,h5err)
      call dump_h5(nid,"kautoknt",kautoknt,h5in,h5err)
      call dump_h5(nid,"akchiwt",akchiwt,h5in,h5err)
      call dump_h5(nid,"akerrwt",akerrwt,h5in,h5err)
      call dump_h5(nid,"kakloop",kakloop,h5in,h5err)
      call dump_h5(nid,"aktol",aktol,h5in,h5err)
      call dump_h5(nid,"kakiter",kakiter,h5in,h5err)
      call dump_h5(nid,"akgamwt",akgamwt,h5in,h5err)
      call dump_h5(nid,"akprewt",akprewt,h5in,h5err)
      call dump_h5(nid,"kpphord",kpphord,h5in,h5err)
      call dump_h5(nid,"kffhord",kffhord,h5in,h5err)
      call dump_h5(nid,"keehord",keehord,h5in,h5err)
      call dump_h5(nid,"psiecn",psiecn,h5in,h5err)
      call dump_h5(nid,"dpsiecn",dpsiecn,h5in,h5err)
      call dump_h5(nid,"fitzts",fitzts,h5in,h5err)
      call dump_h5(nid,"isolve",isolve,h5in,h5err)
      call dump_h5(nid,"iplcout",iplcout,h5in,h5err)
      call dump_h5(nid,"fitfcsum",fitfcsum,h5in,h5err)
      call dump_h5(nid,"fwtfcsum",fwtfcsum,h5in,h5err)
      call dump_h5(nid,"appendsnap",appendsnap,h5in,h5err)
      call dump_h5(nid,"idebug",idebug,h5in,h5err)
      call dump_h5(nid,"nbdrymx",nbdrymx,h5in,h5err)
      call dump_h5(nid,"nsol",nsol,h5in,h5err)
      call dump_h5(nid,"rsol",rsol,h5in,h5err)
      call dump_h5(nid,"zsol",zsol,h5in,h5err)
      call dump_h5(nid,"fwtsol",fwtsol,h5in,h5err)
      call dump_h5(nid,"efitversion",efitversion,h5in,h5err)
      call dump_h5(nid,"kbetapr",kbetapr,h5in,h5err)
      call dump_h5(nid,"nbdryp",nbdryp,h5in,h5err)
      call dump_h5(nid,"jdebug",jdebug,h5in,h5err)
      call dump_h5(nid,"ifindopt",ifindopt,h5in,h5err)
      call dump_h5(nid,"tolbndpsi",tolbndpsi,h5in,h5err)
      call dump_h5(nid,"siloplim",siloplim,h5in,h5err)
      call dump_h5(nid,"use_previous",use_previous,h5in,h5err)
      call close_group("in1",nid,h5err)
   
      call make_group(sid,"ink",nid,group_exists,h5err)
      call dump_h5(nid,"isetfb",isetfb,h5in,h5err)
      call dump_h5(nid,"ioffr",ioffr,h5in,h5err)
      call dump_h5(nid,"ioffz",ioffz,h5in,h5err)
      call dump_h5(nid,"ishiftz",ishiftz,h5in,h5err)
      call dump_h5(nid,"gain",gain,h5in,h5err)
      call dump_h5(nid,"gainp",gainp,h5in,h5err)
      call dump_h5(nid,"idplace",idplace,h5in,h5err)
      call dump_h5(nid,"symmetrize",symmetrize,h5in,h5err)
      call dump_h5(nid,"backaverage",backaverage,h5in,h5err)
      call dump_h5(nid,"lring",lring,h5in,h5err)
      call dump_h5(nid,"cupdown",cupdown,h5in,h5err)
      call close_group("ink",nid,h5err)
   
      call make_group(sid,"ins",nid,group_exists,h5err)
      call dump_h5(nid,"tgamma",tangam(jtime,1:nmselp),h5in,h5err)
      call dump_h5(nid,"sgamma",siggam(jtime,1:nmselp),h5in,h5err)
      call dump_h5(nid,"fwtgam",swtgam,h5in,h5err)
      call dump_h5(nid,"rrrgam",rrgam(jtime,1:nmselp),h5in,h5err)
      call dump_h5(nid,"zzzgam",zzgam(jtime,1:nmselp),h5in,h5err)
      call dump_h5(nid,"aa1gam",a1gam(jtime,1:nmselp),h5in,h5err)
      call dump_h5(nid,"aa2gam",a2gam(jtime,1:nmselp),h5in,h5err)
      call dump_h5(nid,"aa3gam",a3gam(jtime,1:nmselp),h5in,h5err)
      call dump_h5(nid,"aa4gam",a4gam(jtime,1:nmselp),h5in,h5err)
      call dump_h5(nid,"aa5gam",a5gam(jtime,1:nmselp),h5in,h5err)
      call dump_h5(nid,"aa6gam",a6gam(jtime,1:nmselp),h5in,h5err)
      call dump_h5(nid,"aa7gam",a7gam(jtime,1:nmselp),h5in,h5err)
      call dump_h5(nid,"iplots",iplots,h5in,h5err)
      call dump_h5(nid,"kdomse",kdomse,h5in,h5err)
      call dump_h5(nid,"msebkp",msebkp,h5in,h5err)
      call dump_h5(nid,"msefitfun",msefitfun,h5in,h5err)
      call dump_h5(nid,"mse_quiet",mse_quiet,h5in,h5err)
      call dump_h5(nid,"mse_spave_on",mse_spave_on,h5in,h5err)
      call dump_h5(nid,"kwaitmse",kwaitmse,h5in,h5err)
      call dump_h5(nid,"dtmsefull",dtmsefull,h5in,h5err)
      call dump_h5(nid,"mse_strict",mse_strict,h5in,h5err)
      call dump_h5(nid,"t_max_beam_off",t_max_beam_off,h5in,h5err)
      call dump_h5(nid,"ok_30rt",ok_30rt,h5in,h5err)
      call dump_h5(nid,"ok_210lt",ok_210lt,h5in,h5err)
      call dump_h5(nid,"mse_usecer",mse_usecer,h5in,h5err)
      call dump_h5(nid,"mse_certree",mse_certree,h5in,h5err)
      call dump_h5(nid,"mse_use_cer330",mse_use_cer330,h5in,h5err)
      call dump_h5(nid,"mse_use_cer210",mse_use_cer210,h5in,h5err)
      call dump_h5(nid,"tgammauncor",tangam_uncor(jtime,1:nmselp),h5in,h5err)
      call dump_h5(nid,"v30lt",v30lt,h5in,h5err)
      call dump_h5(nid,"v30rt",v30rt,h5in,h5err)
      call dump_h5(nid,"v210lt",v210lt,h5in,h5err)
      call dump_h5(nid,"v210rt",v210rt,h5in,h5err)
      call close_group("ins",nid,h5err)
   
      call make_group(sid,"in_msels",nid,group_exists,h5err)
      call dump_h5(nid,"kdomsels",kdomsels,h5in,h5err)
      call dump_h5(nid,"fmlscut",fmlscut,h5in,h5err)
      call dump_h5(nid,"synmsels",synmsels,h5in,h5err)
      call dump_h5(nid,"avemsels",avemsels,h5in,h5err)
      call close_group("in_msels",nid,h5err)

      call make_group(sid,"ina",nid,group_exists,h5err)
      call dump_h5(nid,"spatial_avg_gam",spatial_avg_gam,h5in,h5err)
      call close_group("ina",nid,h5err)

      call make_group(sid,"inece",nid,group_exists,h5err)
      call dump_h5(nid,"necein",necein,h5in,h5err)
      call dump_h5(nid,"teecein0",teecein0,h5in,h5err)
      call dump_h5(nid,"feece0",feece0,h5in,h5err)
      call dump_h5(nid,"errorece0",errorece0,h5in,h5err)
      call dump_h5(nid,"fwtece0",fwtece0,h5in,h5err)
      call dump_h5(nid,"fwtecebz0",fwtecebz0,h5in,h5err)
      call dump_h5(nid,"ecefit",ecefit,h5in,h5err)
      call dump_h5(nid,"ecebzfit",ecebzfit,h5in,h5err)
      call dump_h5(nid,"kfitece",kfitece,h5in,h5err)
      call dump_h5(nid,"kinputece",kinputece,h5in,h5err)
      call dump_h5(nid,"kcallece",kcallece,h5in,h5err)
      call dump_h5(nid,"nharm",nharm,h5in,h5err)
      call dump_h5(nid,"kfixro",kfixro,h5in,h5err)
      call dump_h5(nid,"rteo",rteo,h5in,h5err)
      call dump_h5(nid,"zteo",zteo,h5in,h5err)
      call dump_h5(nid,"kfixrece",kfixrece,h5in,h5err)
      call dump_h5(nid,"rtep",rtep,h5in,h5err)
      call dump_h5(nid,"rtem",rtem,h5in,h5err)
      call dump_h5(nid,"rpbit",rpbit,h5in,h5err)
      call dump_h5(nid,"rmbit",rmbit,h5in,h5err)
      call dump_h5(nid,"robit",robit,h5in,h5err)
      call dump_h5(nid,"nfit",nfit,h5in,h5err)
      call dump_h5(nid,"kcmin",kcmin,h5in,h5err)
      call dump_h5(nid,"mtxece",mtxece,h5in,h5err)
      call dump_h5(nid,"nconstr",nconstr,h5in,h5err)
      call dump_h5(nid,"eceiter",eceiter,h5in,h5err)
      call dump_h5(nid,"eceerror",eceerror,h5in,h5err)
      call close_group("inece",nid,h5err)

      call make_group(sid,"edgep",nid,group_exists,h5err)
      call dump_h5(nid,"symmetrize",symmetrize,h5in,h5err)
      call dump_h5(nid,"rpress",rpress,h5in,h5err)
      call dump_h5(nid,"pressr",pressr,h5in,h5err)
      call dump_h5(nid,"sigpre",sigpre,h5in,h5err)
      call dump_h5(nid,"npress",npress,h5in,h5err)
      call dump_h5(nid,"kprfit",kprfit,h5in,h5err)
      call dump_h5(nid,"kpressb",kpressb,h5in,h5err)
      call dump_h5(nid,"ndokin",ndokin,h5in,h5err)
      call dump_h5(nid,"kppfnc",kppfnc,h5in,h5err)
      call dump_h5(nid,"kfffnc",kfffnc,h5in,h5err)
      call dump_h5(nid,"kffcur",kffcurs,h5in,h5err)
      call dump_h5(nid,"kppcur",kppcurs,h5in,h5err)
      call dump_h5(nid,"mxiter",mxiters,h5in,h5err)
      call dump_h5(nid,"error",error,h5in,h5err)
      call dump_h5(nid,"errmin",errmin,h5in,h5err)
      call dump_h5(nid,"keecur",keecur,h5in,h5err)
      call close_group("edgep",nid,h5err)

      call make_group(sid,"iner",nid,group_exists,h5err)
      call dump_h5(nid,"keecur",keecur,h5in,h5err)
      call dump_h5(nid,"ecurbd",ecurbd,h5in,h5err)
      call dump_h5(nid,"keefnc",keefnc,h5in,h5err)
      call dump_h5(nid,"eetens",eetens,h5in,h5err)
      call dump_h5(nid,"keebdry",keebdry,h5in,h5err)
      call dump_h5(nid,"kee2bdry",kee2bdry,h5in,h5err)
      call dump_h5(nid,"eebdry",eebdry,h5in,h5err)
      call dump_h5(nid,"ee2bdry",ee2bdry,h5in,h5err)
      call dump_h5(nid,"eeknt",eeknt,h5in,h5err)
      call dump_h5(nid,"keeknt",keeknt,h5in,h5err)
      call dump_h5(nid,"keehord",keehord,h5in,h5err)
      call close_group("iner",nid,h5err)

      call make_group(sid,"insxr",nid,group_exists,h5err)
      call dump_h5(nid,"ksxr0",ksxr0,h5in,h5err)
      call dump_h5(nid,"ksxr2",ksxr2,h5in,h5err)
      call dump_h5(nid,"idosxr",idosxr,h5in,h5err)
      call close_group("insxr",nid,h5err)

      call make_group(sid,"inwant",nid,group_exists,h5err)
      call dump_h5(nid,"psiwant",psiwant,h5in,h5err)
      call dump_h5(nid,"vzeroj",vzeroj,h5in,h5err)
      call dump_h5(nid,"fwtxxj",fwtxxj,h5in,h5err)
      call dump_h5(nid,"fbetap",fbetap,h5in,h5err)
      call dump_h5(nid,"fbetan",fbetan,h5in,h5err)
      call dump_h5(nid,"fli",fli,h5in,h5err)
      call dump_h5(nid,"fqsiw",fqsiw,h5in,h5err)
      call dump_h5(nid,"jbeta",jbeta,h5in,h5err)
      call dump_h5(nid,"jli",jli,h5in,h5err)
      call dump_h5(nid,"alpax",alpax,h5in,h5err)
      call dump_h5(nid,"gamax",gamax,h5in,h5err)
      call dump_h5(nid,"jwantm",jwantm,h5in,h5err)
      call dump_h5(nid,"fwtxxq",fwtxxq,h5in,h5err)
      call dump_h5(nid,"fwtxxb",fwtxxb,h5in,h5err)
      call dump_h5(nid,"fwtxli",fwtxli,h5in,h5err)
      call dump_h5(nid,"znose",znose,h5in,h5err)
      call dump_h5(nid,"fwtbdry",fwtbdry,h5in,h5err)
      call dump_h5(nid,"nqwant",nqwant,h5in,h5err)
      call dump_h5(nid,"siwantq",siwantq,h5in,h5err)
      call dump_h5(nid,"kccoils",kccoils,h5in,h5err)
      call dump_h5(nid,"ccoils",ccoils,h5in,h5err)
      call dump_h5(nid,"rexpan",rexpan,h5in,h5err)
      call dump_h5(nid,"xcoils",xcoils,h5in,h5err)
      call dump_h5(nid,"kcloops",kcloops,h5in,h5err)
      call dump_h5(nid,"cloops",cloops,h5in,h5err)
      call dump_h5(nid,"xloops",xloops,h5in,h5err)
      call dump_h5(nid,"nccoil",nccoil,h5in,h5err)
      call dump_h5(nid,"sizeroj",sizeroj,h5in,h5err)
      call dump_h5(nid,"fitdelz",fitdelz,h5in,h5err)
      call dump_h5(nid,"relaxdz",relaxdz,h5in,h5err)
      call dump_h5(nid,"stabdz",stabdz,h5in,h5err)
      call dump_h5(nid,"table_dir",table_dir,h5in,h5err)
      call dump_h5(nid,"errdelz",errdelz,h5in,h5err)
      call dump_h5(nid,"oldccomp",oldccomp,h5in,h5err)
      call dump_h5(nid,"nicoil",nicoil,h5in,h5err)
      call dump_h5(nid,"oldcomp",oldcomp,h5in,h5err)
      call dump_h5(nid,"currc79",curc79(jtime),h5in,h5err)
      call dump_h5(nid,"currc139",curc139(jtime),h5in,h5err)
      call dump_h5(nid,"currc199",curc199(jtime),h5in,h5err)
      call dump_h5(nid,"curriu30",curiu30(jtime),h5in,h5err)
      call dump_h5(nid,"curriu90",curiu90(jtime),h5in,h5err)
      call dump_h5(nid,"curriu150",curiu150(jtime),h5in,h5err)
      call dump_h5(nid,"curril30",curil30(jtime),h5in,h5err)
      call dump_h5(nid,"curril90",curil90(jtime),h5in,h5err)
      call dump_h5(nid,"curril150",curil150(jtime),h5in,h5err)
      call dump_h5(nid,"ifitdelz",ifitdelz,h5in,h5err)
      call dump_h5(nid,"scaledz",scaledz,h5in,h5err)
      call close_group("inwant",nid,h5err)
   
      call make_group(sid,"invt",nid,group_exists,h5err)
      call dump_h5(nid,"omegat",omegat,h5in,h5err)
      call dump_h5(nid,"nomegat",nomegat,h5in,h5err)
      call dump_h5(nid,"enw",enw,h5in,h5err)
      call dump_h5(nid,"emw",emw,h5in,h5err)
      call dump_h5(nid,"betapw0",betapw0,h5in,h5err)
      call dump_h5(nid,"kdovt",kdovt,h5in,h5err)
      call dump_h5(nid,"kvtor",kvtor,h5in,h5err)
      call dump_h5(nid,"kwwcur",kwwcur,h5in,h5err)
      call dump_h5(nid,"rvtor",rvtor,h5in,h5err)
      call dump_h5(nid,"wcurbd",wcurbd,h5in,h5err)
      call dump_h5(nid,"preswb",preswb,h5in,h5err)
      call dump_h5(nid,"fwtprw",fwtprw,h5in,h5err)
      call dump_h5(nid,"npresw",npresw,h5in,h5err)
      call dump_h5(nid,"presw",presw,h5in,h5err)
      call dump_h5(nid,"sigprw",sigprw,h5in,h5err)
      call dump_h5(nid,"rpresw",rpresw,h5in,h5err)
      call dump_h5(nid,"zpresw",zpresw,h5in,h5err)
      call dump_h5(nid,"kplotp",kplotp,h5in,h5err)
      call dump_h5(nid,"nsplot",nsplot,h5in,h5err)
      call dump_h5(nid,"sbetaw",sbetaw,h5in,h5err)
      call dump_h5(nid,"comega",comega,h5in,h5err)
      call dump_h5(nid,"kcomega",kcomega,h5in,h5err)
      call dump_h5(nid,"xomega",xomega,h5in,h5err)
      call dump_h5(nid,"romegat",romegat,h5in,h5err)
      call dump_h5(nid,"zomegat",zomegat,h5in,h5err)
      call dump_h5(nid,"sigome",sigome,h5in,h5err)
      call dump_h5(nid,"scalepw",scalepw,h5in,h5err)
      call dump_h5(nid,"kwwfnc",kwwfnc,h5in,h5err)
      call dump_h5(nid,"kwwknt",kwwknt,h5in,h5err)
      call dump_h5(nid,"wwknt",wwknt,h5in,h5err)
      call dump_h5(nid,"wwtens",wwtens,h5in,h5err)
      call close_group("invt",nid,h5err)

      call close_group(trim(tindex),sid,h5err)
      call close_group("time_slice",tid,h5err)
      call close_group("parameters",pid,h5err)
      call close_group("code",cid,h5err)
      call close_group("equilibrium",eqid,h5err)
      call close_h5file(fileid,rootgid,h5err)
!----------------------------------------------------------------------
!--   Restore variables                                              --
!----------------------------------------------------------------------
      rbdry=rbdryss
      zbdry=zbdryss
      nbdry=nbdryss
      n1coil=n1coilt

      return
      end subroutine write_omas_input
