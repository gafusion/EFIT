#include "config.f"
! set_kinds
      module set_kinds
!**     set the type of variables like integer, real, etc...
        !HP integer, parameter :: rprec = selected_real_kind(20,100)
        !integer*4, parameter :: rprec = selected_real_kind(13,307)
        !integer*4, parameter :: iprec = selected_real_kind(4)
        !integer*4, parameter :: dp=rprec
        integer*4, parameter :: dp=selected_real_kind(15,307) ! REAL*8
        integer*4, parameter :: i4=selected_int_kind(9)
        integer*4, parameter :: i8=selected_int_kind(18)
        integer*4, parameter :: r4=selected_real_kind(6,37)
        !integer*4, parameter :: r8=selected_real_kind(13,307)
      end module set_kinds
!extvars
      module extvars
        public
!       table_dir    area where green tables are stored
!       input_dir    area where other input files are stored
!       store_dir    central directory to collect EFIT results
        character(256) table_dir,input_dir,store_dir,table_di2, &
                       link_efit,link_store
        integer*4 ltbdir,lindir,lstdir,ltbdi2
        character(7) :: efitvers ! git hash
        character(10) :: efitdate,efitdatealt ! commit dates (legacy formats)
        integer*4 :: efitversion ! commit date (legacy form)
        character(512) :: giturl,gitbranch,githash,gitdate,gitdatealt,lmods
        character(512) :: fc,fc_id,fcver
        character(5096) :: fcflags
        character(256) :: arch_type,hostname
      end module extvars
!eparm
      module eparm
        public
        ! Experiment dependant parameters
        character(10) :: device !< machine name
        integer*4 :: nsilop !< number of flux loops
        integer*4 :: nfsum  !< number of p.f. f-coil groups
        integer*4 :: nrogow !< number of rogowski coils
        integer*4 :: nacoil !< number of advance divertor coils
        integer*4 :: nfcoil !< number of p.f. f-coils
        integer*4 :: necoil !< number of p.f e-coils
        integer*4 :: mpress !< number of pressure data points
        integer*4 :: nvesel !< number of vessel segements
        integer*4 :: nesum  !< number of p.f. e-coil groups
        integer*4 :: magpri !< number of magnetic detectors
        integer*4 :: libim  !< number of lithium ion beam measurements
        integer*4 :: nmselp !< number of mse-lp channels
        integer*4 :: nstark !< number of li beam and mse-lp measurements
        integer*4 :: nmsels !< number of mse-ls channels
        integer*4 :: nnece  !< number of ece channels
        integer*4 :: nnecein,neceo,nnnte
        integer*4 :: ngam_vars,ngam_u,ngam_w !< dimensions of mse spatial averaging data
        integer*4 :: nlimit !< maximum number of limiter points
        integer*4 :: nlimbd !< maximum number of 'outer' limiter points
        integer*4 :: nangle !< nangle dimension of poloidal sxr, first part of xangle,zxray,rxray
        integer*4 :: ntangle !< dimension of toroidal xray, last part of xangle,zxray,rxray
        integer*4 :: nfbcoil !< nfbcoil (obsolete)
        integer*4 :: mccoil,micoil
        ! General parameters
        integer*4 :: ntime !< number of time slices
        integer*4 :: ndata
        integer*4 :: nwwcur
        integer*4 :: nffcur,nppcur,npcurn,necur, &
                     mfnpcr,nercur,nrsmat,nwcurn
        integer*4 :: msbdr2
        integer*4 :: ndim_crv
        integer*4 :: ndim,kxiter,kxiter_save,mqwant
        integer*4 :: nw,nh,nw_sub,nh_sub
        integer*4 :: nwnh,nh2,nwork,nwwf, &
                     nwf,lubicx,lubicy,kujunk,boundary_count, &
                     lr0,lz0,npoint,nxtrap
        integer*4 :: kubicx,kubicy
        integer*4 :: ncurrt
        integer*4 :: mbdry, mbdry1
        integer*4 :: msbdry
        integer*4 :: nxtram
        integer*4 :: nxtlim,nco2v,nco2r
        integer*4 :: modef,modep,modew
        integer*4 :: icycred_loopmax
        integer*4 :: nfourier !< nfourier number Fourier components of vessel current
      end module eparm
!errlims
      module errlims
        ! Experiment dependant checks on the solution error
        real*8 li_max,li_min,betap_max,plasma_diff, &
               aminor_max,aminor_min,elong_max,elong_min, &
               rout_max,rout_min,zout_max,zout_min, &
               rcurrt_max,rcurrt_min,zcurrt_max,zcurrt_min, &
               qstar_max,qstar_min,betat_max, &
               gapin_min,gapout_min,gaptop_min, &
               sepin_check,qout_max,qout_min, &
               dbpli_diff,delbp_diff
      end module errlims
! global_constants
      ! Calculate and store global constants like pi, e, gravity, etc.
      module global_constants
        use set_kinds, only: dp
        public
        real*8 :: pi,twopi,radeg,tmu0,tmu02
        real*8, parameter :: tmu = 2.0e-07_dp
      end module global_constants
!var_nio
      module var_nio
#if defined(USE_HDF5)
        use hdf5_api
#endif
        integer*4 nsnapf
        character*2 appendsnap
        character*86 snapextin
        character*100 snapfile,tmpdata
        integer*4, parameter :: nin=11,nout=10,ntty=5,nrsppc=25,neqdsk=38,nffile=40, &
                                nsave=49,nttyo=6
#if defined(USE_HDF5)
        type(hdf5ErrorType) :: h5err
        type(hdf5InOpts), save :: h5in
        integer(HID_T) :: fileid,rootgid,eqid,cid,pid,tid,sid,fid,nid
        contains
!       subprogram 1. fch5init.
!       initialize fcio h5 information for writes and reads.
        subroutine fch5init
          logical, save :: h5init=.false.
          if (.not.h5init) then
            call vshdf5_fcinit
            call vshdf5_inith5vars(h5in, h5err)
            h5in%verbose=.false.
            h5in%debug=.false.
            h5init=.true.
          endif
        end subroutine fch5init
!       subprogram 2. fch5resetvars.
!       reset h5 information for writes and reads.
        subroutine fch5resetvars
          h5in%mesh =  " "
          h5in%vsAxisLabels = " " 
          h5in%units =  " "
          h5in%vsCentering =  " "
          h5in%vsMD =  " "
          h5in%vsIndexOrder = " "
          h5in%vsLabels = " "
        end subroutine fch5resetvars
#endif
      end module var_nio
! error_control
      ! Error control and write out error messages consistently.
      module error_control
        use set_kinds, only: dp
        integer*4 :: currrank=-1 ! some invalid value
        real(dp) :: currtime=-1.0 ! some invalid value
        public :: errctrl_msg

      contains
        subroutine errctrl_setstate(currrank0,currtime0)
          use set_kinds, only: dp
          integer*4, intent(in) :: currrank0
          real(dp), intent(in) :: currtime0
          currrank = currrank0
          currtime = currtime0
        end subroutine errctrl_setstate

        subroutine errctrl_msg(subrstr,msgstr,mtype0)
          use var_nio, only: nttyo
          use set_kinds, only: dp
          character(len=*), intent(in) :: subrstr,msgstr
          integer*4, optional :: mtype0
          integer*4 :: mtype
          character(len=16) :: labelstr

          mtype = 1
          if (present(mtype0)) mtype = mtype0

          select case(mtype)
          case (2)
            labelstr = 'WARNING'
          case (3)
            labelstr = 'INFO'
          case default
            labelstr = 'ERROR'
          end select

          write(nttyo, '(a,a,a,a,i3,a,i6,a,a)') trim(labelstr),' in ', &
            subrstr,' at r=',currrank,', t=',int(currtime+1.e-8_dp),': ',msgstr

!         TODO: can cause race condition with MPI
!          open(unit=40,file='errfil.out',status='unknown', &
!               position='append')
!          write(40,'(a,a,a,a,i3,a,i6,a,a)') trim(labelstr),' in ', &
!            subrstr,' at r=',currrank,', t=',int(currtime),': ',msgstr
!          close(unit=40)

        end subroutine errctrl_msg
      end module error_control
!var_filech
      module var_filech
        character*4 ch1,ch2
      end module var_filech
!var_outp1
      module var_outp1
        real*8 sbpp,sbppa
      end module var_outp1
!var_iopen
      module var_iopen
        integer*4 ifcurr 
        real*8 errbry,relax
      end module var_iopen
!var_zcntrl
      module var_zcntrl
        integer*4 isetfb,ishiftz,ioffz,ioffr,idplace,lring 
        real*8 gain,gainp,delzmm
      end module var_zcntrl
!var_updown
      module var_updown
        real*8 dpsip,dpsip_last
        logical symmetrize,backaverage
      end module var_updown
!var_test
      module var_test
        real*8 zcontr,zcurnow
      end module var_test
!var_fitsiref
       module var_fitsiref
         real*8 csiref,saisref,fwtref,scalesir
         logical fitsiref
       end module var_fitsiref
!var_cnnn
       module var_cnnn
         integer*4, parameter :: n111=1,n222=2,n333=3,n444=4,n555=5,n666=6
         real*8, parameter :: x000=0.0,x111=1.0
       end module var_cnnn
!var_pcsys
       module var_pcsys
         integer*4 use_alternate_pointnames
         character*80 alternate_pointname_file
         logical*4 do_spline_fit
       end module var_pcsys
!var_pfedge
      module var_pfedge
        integer*4 kedgep,kedgef 
        real*8 pedge,pe_psin,pe_width, &
          f2edge,fe_psin,fe_width,constf2, &
          tpedge,tfedge,rdlcfe,rqape,rqafe,betped,betnped
      end module var_pfedge
!var_sxpoint
       module var_sxpoint
         real*8 sissep,rssep,zssep,sifsep,rfsep,zfsep,rminss, &
                zminss,rmaxss,zmaxss,rminfs,rmaxfs,zminfs,zmaxfs
       end module var_sxpoint
!var_consta
      module var_consta
        public
        real*8 :: tmu2,errcut
        integer*4 :: ibunmn,kinput,kcaldia,ibunmns
      end module var_consta
!var_rcfact
      module var_rcfact
        integer*4 ircfact
      end module var_rcfact
!var_curpo
      module var_curpro
        real*8 emf,emp,enf,enp,rbetap,rzero,pbetap,qenp,qemp,qenf
      end module var_curpro
!var_pfterm
      module var_pfterm
        integer*4 kffcur,kppcur,kpcurn,icalbet,kffcurs,kppcurs 
        real*8 chidflux,gammaf,gammap,cstab0, &
               vbtot2,vbtvac2,vbtor2,vbtvac,vbeta0, &
               vbtmag,btvvac2,btvtor2,btvtot2
      end module var_pfterm
!var_cfit
      module var_cfit
        integer*4 mxiter,nitera,nxiter,ixnn,isolve 
        real*8 error,errorm,errmin,delerr,delerb
      end module var_cfit 
!var_cgrid
      module var_cgrid
        real*8 darea,drgrid,dzgrid,qmaxis,cratio,dfsqe,cratiof
      end module var_cgrid
!var_extra
      module var_extra
        real*8 scrape,tolbndpsi
        integer*4 nextra,ixstrt,iprobe,iecoil,iexcal,iconsi,iqplot, &
                  klabel,ifindopt,iexcals
      end module var_extra
!var_conveg
      module var_conveg
        real*8 omega,relip,zelip,aelip,eelip,omecur 
      end module var_conveg
!var_limmm
      module var_limmm
        real*8 xlmin,xlmax,ylmin,ylmax
        integer*4 limid,limup,limbot
      end module var_limmm
!var_inaver
      module var_inaver
        integer*4 iavem,iaved,iavev,iaveus
      end module var_inaver
!var_vessel
      module var_vessel
        real*8,dimension(:),allocatable :: volecs,volecc,rsisec
        real*8,dimension(:),allocatable :: volfcs,volfcc 
        real*8,dimension(:),allocatable :: rvs,zvs,hvs,wvs,avs,avs2,rsisvs 
        real*8 powvs,pvscur,pscurn,ppscur,efreq
        integer*4 ivesel
      end module var_vessel
!var_cyclic_red
      module var_cyclic_red
        real*8,dimension(:,:),allocatable :: beti,abeti,wk1
        real*8,dimension(:),allocatable :: alphab,diag1 
        real*8,dimension(:),allocatable :: rhsdumy1
        real*8,dimension(:),allocatable :: phi,v,wk2,diagl,diagu
        real*8,dimension(:),allocatable :: tempgrid,tempgrid2
        real*8 diag,rhs_a_dumy,rhs_b_dumy
        integer*4 nhpwr
      end module var_cyclic_red
!var_scalem
      module var_scalem
        !use eparm,only:nrsmat,mfnpcr
        integer*4 infosc
        real*8,dimension(:), allocatable :: rowscale
        real*8,dimension(:), allocatable :: colscale
        real*8 rowcnd,colcnd,arspmax
        logical scalea
      end module var_scalem

      module var_solove
        integer*4 islve
        real*8 salpha,sbeta,srm,scc1,seee,saaa,srma, sbetaw
      end module var_solove 

      module var_buneman
        integer*4 :: mno,m,n
        real*8    :: drdz2,rgrid1,delrgrid,delz
        real*8    :: s,shift,dr,dz
      end module var_buneman

!var_contor
      module var_contor
        real*8,dimension(:),allocatable :: s1,s2,s3,bpolav
      end module var_contor
!var_mfield
      module var_mfield
      !use eparm,only:npoint
        real*8,dimension(:),allocatable :: bpol,plengt,bpolz
        real*8 siar,siaz
      end module var_mfield
!var_hist
      module var_hist
        integer*4,dimension(:), allocatable :: jerror
        real*8,dimension(:), allocatable :: elong,rout,zout,utri, &
          ltri,aminor,volume,betat,gaptop, &
          betap,li,gapin,gapout,qstar, &
          rcurrt,zcurrt,qout,sepin, &
          sepout,septop,sibdry,area, &
          wmhd,elongm,qm,terror, &
          rm,zm,gapbot,sepbot, &
          alpha,rttt,dbpli,delbp,oring, &
          sepexp,shearb, &
          xtch,ytch,q95,vertn,aaq1, &
          aaq2,aaq3,btaxp,btaxv, &
          psim,dsep,peak, &
          wbpol,taumhd,betapd,betatd, &
          wdia,taudia,wbpold, &
          qmerci,slantu,slantl,zeff, &
          zeffr,tave,rvsin,zvsin, &
          rvsout,zvsout,wpdot,wbdot, &
          vsurfa,cjor95,pp95,drsep, &
          yyy2,xnnc,wtherm,wfbeam,taujd3,tauthn, &
          li3,tflux,twagap
        real*8,dimension(:,:), allocatable :: rseps,zseps
        real*8 pasman,betatn,psiq1,betat2
        integer*4 jtwagap
      end module var_hist
!var_hist2
      module var_hist2
        real*8,dimension(:), allocatable :: qsiwant,cjorsw,cjor0, &
          ssiwant,ssi95,cjor99,cj1ave, &
          rmidin,rmidout,psurfa
        real*8 psiwant,rexpan,fexpan,qmin,fexpvs,shearc, &
          sepnose,ssi01,znose,rhoqmin
      end module var_hist2
!var_cshape
      module var_cshape
        real*8,dimension(:),allocatable :: xout,yout
        real*8 dpsi,rymin,rymax, &
          zxmin,zxmax,xmin,xmax,ymin,ymax,rmaxis,zmaxis,emaxis, &
          rminzm,rmaxzm,delrmax1,delrmax2
        integer*4 nfound
      end module var_cshape
!var_divdis
      module var_divdis
        real*8,dimension(:), allocatable :: dolubaf,dolubafm,diludom, &
                                            diludomm,dminux,dminlx, &
                                            ratsol,rvsiu,zvsiu,rvsou, &
                                            zvsou,rvsid,zvsid,rvsod, &
                                            zvsod
      end module var_divdis
!var_cpsi
      module var_cpsi
        real*8,dimension(:),allocatable :: psi,xpsi,vfbrrt,psipla
        real*8 vcurfb(3)
        real*8 psibry,simag,sidif,eouter,zplasm,zpwant,difpsi, &
               cupdown
      end module var_cpsi
!var_cvalue
      module var_cvalue
        real*8,dimension(:,:), allocatable :: csilop,csilopv 
        real*8,dimension(:,:), allocatable :: crogow
        real*8,dimension(:,:), allocatable :: cmpr2,cmpr2v 
        real*8,dimension(:), allocatable :: ipmhd,indent 
        real*8,dimension(:,:), allocatable :: ccbrsp
        real*8,dimension(:,:), allocatable :: caccurt
        real*8 cli,cqqxis,ci0,cipmp2 
      end module var_cvalue
!var_gtable
      module var_gtable
        real*8,dimension(:,:),allocatable,save :: gridfc
        real*8,dimension(:),allocatable,save :: rgrid
        real*8,dimension(:),allocatable,save :: zgrid
        real*8,dimension(:,:),allocatable,save :: rsilfc
        real*8,dimension(:,:),allocatable,save :: rmp2fc
        real*8,dimension(:,:),allocatable,save :: gridpc
        real*8,dimension(:,:),allocatable,save :: gsilpc
        real*8,dimension(:,:),allocatable,save :: gmp2pc
        real*8,dimension(:,:),allocatable,save :: rfcfc
        integer*4 iallocate_stat
      end module var_gtable
! jm.s
! NOTE : array sizes are grid size-dependent so they should be dynamically allocated, but
!        we cannot make them dynamic since they are included in a namelist
! NOTE : the largest possible grid size with Buneman's algorithm is 2049
!        (will not be the case with pefit)
! NOTE : npsi_ext (actual dimension of _ext arrays) used in code logic
      module profile_ext_mod
        integer*4 :: npsi_ext,nw_ext,nh_ext,nbdry_ext,limitr_ext
        real*8,dimension(2049) :: psin_ext
        real*8,dimension(2049) :: bpp_ext,cpp_ext,dpp_ext
        real*8,dimension(2049) :: bfp_ext,cfp_ext,dfp_ext
        real*8,dimension(:),allocatable :: rbdry_ext,zbdry_ext,xlim_ext, &
                                           ylim_ext,psirz_ext,pprime_ext,&
                                           ffprim_ext,qpsi_ext,fcoil_ext
        real*8 :: sign_ext,scalepp_ext,scaleffp_ext,cratio_ext, &
                  cratiop_ext,cratiof_ext,simag_ext,psibry_ext
        character*80 :: geqdsk_ext
        logical :: fixpp
      end module profile_ext_mod

! NOTE : keep track of times for which BCOIL and ECOIL data exist (see measurents.F90)
      module vtime_mod
        integer*4 :: nvtime
        real*8,dimension(:), allocatable :: vtime
        integer*4 :: ntims
        integer*4, parameter :: npmax=262144 ! sufficient for ms data from 262s shot... needs to be larger than ntims and match npmax in getdat.F90
      end module vtime_mod
