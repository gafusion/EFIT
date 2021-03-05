!exparm
      module exparm
      public

!     2000/03/11 00:01:29 lao
!     @(#)exparm.inc,v 1.6
!
!-------------------------------------------------------------------
!--  New E-coil connections                       LLao, 95/07/11  --
!--  Add 8 new probes for radiative divertor      LLao, 97/03/17  --
!--  Update MSE to 35 channels                    LLao, 97/04/21  --
!--  Separate machine dependent configuration                     --
!--   parameters from eparmdx.for    QPeng,97/09/24  --
!-- added ntangle for toroidal x-ray   QPeng,98/05/12  --
!--  Increase MSE channels from 35 to 36                98/12/01  --
!-------------------------------------------------------------------
!
!     magpri67 number of magnetic detectors at toroidal angle "1"
!     magpri322 number of magnetic detectors at toroidal angle "2"
!     magprirdp number of magnetic detectors for radiative divertor
!     magpri total number of magnetic detectors
!     mpress number of pressure data points
!     mse315
!     mse45
!     mse15
!     mse210
!     nstark total number of mse channels
!     ngam_vars, ngam_u, ngam_w dimensions of mse spatial averaging data
!heng necein    total number of ece channels
!     nacoil number of advance divertor coils
!     nangle dimension of poloidal sxr, first part of xangle,zxray,rxray
!     ntangle dimension of toroidal xray, last part of xangle,zxray,rxray
!     necoil number of ohmic heating coils
!     nesum number of p.f. coil groups
!     nfbcoil (obsolete)
!     nfcoil number of p.f. coils
!     nlimbd number of 'outer' limiter points
!     nlimit maximum number of limiter points
!     nsilop number of flux loops
!     nvesel number of vessel segements
!
!
      parameter (nsilds=3,nsilol=41)
      parameter (nsilop=nsilds+nsilol)
      parameter (nfcoil=18,nrogow=1,nacoil=1)
      parameter (mfcoil=18)
      parameter (necoil=122,nvesel=24,mpress=201)
      parameter (nesum=6)
      parameter (magpri67=29,magpri322=31,magprirdp=8,magudom=5)
      parameter (maglds=3)
      parameter (magpol=magpri67+magpri322+magprirdp+magudom)
      parameter (magpri=magpol+maglds)
      parameter (mse315=11,mse45=15,mse15=10,mse1h=4,mse315_2=5)
      parameter (mse210=24)
      parameter (libim=32)
      parameter (nmtark=mse315+mse45+mse15+mse1h+mse315_2+mse210)
      parameter (nstark=nmtark+libim)
      parameter (nmsels=16)
      parameter (nnece=40,nnecein=80,neceo=1,nnnte=801)
      parameter (ngam_vars=9,ngam_u=5,ngam_w=3)
      parameter (nlimit=160,nlimbd=6)
      parameter (nangle=64,ntangle=12)
      parameter (nfbcoil=12)
      parameter (mccoil=6,micoil=12)

      end module exparm

      module set_kinds
!**     set the type of variables like integer, real, etc...
        !HP integer, parameter :: rprec = selected_real_kind(20,100)
        !integer, parameter :: rprec = selected_real_kind(13,307)
        !integer, parameter :: iprec = selected_real_kind(4)
        !integer, parameter :: dp=rprec
        integer, parameter :: dp = selected_real_kind(15,307) ! REAL*8

       end module set_kinds

!expath
      module expath
      public

!     $Date: 2009/02/24 01:37:39 $ $Author: lao $
!     @(#)$RCSfile: modules-efitx.f90,v $ $Revision: 1.1.2.4 $
!
!     table_dir    area where green tables are stored, including
!   re*,rv*,fc*,ra*,ef*,fc*,fm*,ec*,rf*,ep*
!   rs*, dprobe.dat, pcsnames.dat,
!   (-fitweight.dat., -sxr94.dat, -efithelp.txt)
!     input_dir    area where other input files are stored, they include
!   n1coil.d3d
!   ccoil.d3d
!   btcomp.d3d
!   dcoef.dat
!   limiter files
!   snap files
!   + efithelp.txt
!   + fitweight.dat
!   + sxr94.dat
!     store_dir    central directory to collect EFIT results
!

!vas      common/datafile/table_dir,input_dir,store_dir, &
!vas              ltbdir,lindir,lstdir
      character*82 table_dir,input_dir,store_dir,table_di2,link_efitx,link_storex
      integer*4 ltbdir,lindir,lstdir,ltbdi2
!      data table_dir /'/link/efit/new_table/'/
!      data table_dir /'/task/efit/lao/efits/CVS/p/2006/'/
!sri-mpi
!      data table_dir /'/link/efit/2006/'/
      data table_dir /'link_efitx/2006/'/
!org      data table_dir /'/task/efit/lao/efits/CVS/p/2006/'/
!vas      data table_dir /'/u4/radhakri/Vas-green/2006/'/
!vas      data table_dir /'/u4/radhakri/Vas-green/2008/'/
!old      data table_dir /'/u4/radhakri/Renq-green/hydra/2008/'/
!unified table      data table_dir /'/task/imd/radhakri/Run/2008/'/
!lf95      data table_dir /'/u4/radhakri/Renq-green/lf95/2008/'/
!vas for linux      data table_dir /'/u4/radhakri/green/2008/'/
      data input_dir /'link_efitx/'/
      data store_dir /'link_storex/'/

      end module expath
!eparmdud129
      module eparmdud129
      public
!vas      implicit integer*4 (i-n), real*8 (a-h, o-z)
!
!   1997/10/09 00:01:35 peng
!
!  @(#)eparmdx.for,v 4.19
!  

!     2000/03/11 00:01:29 lao
!     @(#)exparm.inc,v 1.6
!
!-------------------------------------------------------------------
!--  New E-coil connections                       LLao, 95/07/11  --
!--  Add 8 new probes for radiative divertor      LLao, 97/03/17  --
!--  Update MSE to 35 channels                    LLao, 97/04/21  --
!--  Separate machine dependent configuration                     --
!--   parameters from eparmdx.for    QPeng,97/09/24  --
!-- added ntangle for toroidal x-ray   QPeng,98/05/12  --
!--  Increase MSE channels from 35 to 36                98/12/01  --
!-------------------------------------------------------------------
!
!     magpri67 number of magnetic detectors at toroidal angle "1"
!     magpri322 number of magnetic detectors at toroidal angle "2"
!     magprirdp number of magnetic detectors for radiative divertor
!     magpri total number of magnetic detectors
!     mpress number of pressure data points
!     mse315
!     mse45
!     mse15
!     mse210
!     nstark total number of mse channels
!     ngam_vars, ngam_u, ngam_w dimensions of mse spatial averaging data
!heng necein    total number of ece channels
!     nacoil number of advance divertor coils
!     nangle dimension of poloidal sxr, first part of xangle,zxray,rxray
!     ntangle dimension of toroidal xray, last part of xangle,zxray,rxray
!     necoil number of ohmic heating coils
!     nesum number of p.f. coil groups
!     nfbcoil (obsolete)
!     nfcoil number of p.f. coils
!     nlimbd number of 'outer' limiter points
!     nlimit maximum number of limiter points
!     nsilop number of flux loops
!     nvesel number of vessel segements
!
!
      parameter (nsilds=3,nsilol=41)
      parameter (nsilop=nsilds+nsilol)
      parameter (nfcoil=18,nrogow=1,nacoil=1)
      parameter (mfcoil=18)
      parameter (necoil=122,nvesel=24,mpress=201)
      parameter (nesum=6)
      parameter (magpri67=29,magpri322=31,magprirdp=8,magudom=5)
      parameter (maglds=3)
      parameter (magpol=magpri67+magpri322+magprirdp+magudom)
      parameter (magpri=magpol+maglds)
      parameter (mse315=11,mse45=15,mse15=10,mse1h=4,mse315_2=5)
      parameter (mse210=24)
      parameter (libim=32)
      parameter (nmtark=mse315+mse45+mse15+mse1h+mse315_2+mse210)
      parameter (nstark=nmtark+libim)
      parameter (nmsels=16)
      parameter (nnece=40,nnecein=80,neceo=1,nnnte=801)
      parameter (ngam_vars=9,ngam_u=5,ngam_w=3)
      parameter (nlimit=160,nlimbd=6)
      parameter (nangle=64,ntangle=12)
      parameter (nfbcoil=12)
      parameter (mccoil=6,micoil=12)

!
! --- experiment dependant parameters
!
!vas      include 'exparm.inc'
!
! --- general parameters
!
      parameter (ndata=61)
      parameter (nwwcur=32)
      parameter (nffcur=32,nppcur=32,npcurn=nffcur+nppcur &
           ,nercur=32, necur2=nercur*2 &
           ,mfnpcr=nfcoil+npcurn+nvesel+nwwcur+nesum+nfcoil &
                   + nercur &
           ,npcur2=npcurn*2 &
           ,nrsmat=nsilop+magpri+nrogow+nffcur+1+npcurn+nwwcur+ &
            mpress+nfcoil+nstark+nnece+neceo &
           ,nwcurn=nwwcur+npcurn,npcur3=npcurn*2 &
           ,nwcur2=nwcurn*2)
      parameter (ntime=1001)
!!257g      parameter (nw=257,nh=257,nwnh=nw*nh)
!!513g      parameter (nw=513,nh=513,nwnh=nw*nh)
!!129g      parameter (nw=129,nh=129,nwnh=nw*nh)
!!65g      parameter (nw=65,nh=65,nwnh=nw*nh)
!!lowg      parameter (npoint=800)
!!lowg      parameter(ndim = 700)
!!lowg      parameter (kxiter=250,mqwant=30)
!!higg      parameter (npoint=3200)
!!higg      parameter(ndim = 3200)
!!higg      parameter (kxiter=nw+2,mqwant=30)
      parameter (ndim=3200,kxiter=515,mqwant=30)
      integer :: nw,nh,nwnh,nh2,nwork,nwwf, &
           nwf,lubicx,lubicy,kujunk,boundary_count, &
           lr0,lz0,kubicx,kubicy,npoint,nxtrap
!      parameter (nh2=2*nh,nwrk=2*(nw+1)*nh)
      parameter (ncurrt=nvesel+nesum+nfcoil)
      parameter (mbdry=300, mbdry1=110)
      parameter (nbwork=nsilop)
      parameter (msbdry=mbdry+nsilop+nfcoil+1,msbdr2=2*msbdry)
      parameter (nrsma2=2*nrsmat)
!      parameter (nwwf=2*nw)
!      parameter (nwf=nwwf)
!      parameter (nxtram=10,nxtrap=npoint)
      parameter (nxtram=10)
      parameter (nxtlim=9,nco2v=3,nco2r=2)
!      parameter (kubicx = 4, kubicy = 4, lubicx = nw - kubicx + 1, &
!                 lubicy = nh - kubicy + 1, &
!                 kujunk = kubicx*kubicy*lubicx*lubicy)
      parameter (modef=4, modep=4, modew=4 , kubics=4 )
      parameter (icycred_loopmax=1290)
!      parameter (boundary_count=2*nh+2*(nw-2))
! added by Qian for Fourier components of vessel current
      parameter (nfourier=5)
!
!
! --- common block for areas that store green tables and general inputs
!
!vas      include 'expath.inc'

!     $Date: 2009/02/24 01:37:39 $ $Author: lao $
!     @(#)$RCSfile: modules-efitx.f90,v $ $Revision: 1.1.2.4 $
!
!     table_dir    area where green tables are stored, including
!   re*,rv*,fc*,ra*,ef*,fc*,fm*,ec*,rf*,ep*
!   rs*, dprobe.dat, pcsnames.dat,
!   (-fitweight.dat., -sxr94.dat, -efithelp.txt)
!     input_dir    area where other input files are stored, they include
!   n1coil.d3d
!   ccoil.d3d
!   btcomp.d3d
!   dcoef.dat
!   limiter files
!   snap files
!   + efithelp.txt
!   + fitweight.dat
!   + sxr94.dat
!     store_dir    central directory to collect EFIT results
!

!vas      common/datafile/table_dir,input_dir,store_dir, &
!vas              ltbdir,lindir,lstdir
!test      character*42 table_dir,input_dir,store_dir
!test      integer ltbdir,lindir,lstdir
!      data table_dir /'/link/efit/new_table/'/
!test      data table_dir /'/task/efit/lao/efits/CVS/p/2006/'/
!test      data input_dir /'/link/efit/'/
!test      data store_dir /'/link/store/'/

      end module eparmdud129

      ! Calculate and store global constants like pi, e, gravity, etc.
      module global_constants
        use set_kinds
        public
        real*8 :: pi=0,twopi=0,tmu=0,radeg=0
      contains
        subroutine set_constants()
          pi = 4.0_dp*atan(1.0_dp) ! calculate pi to machine precision
          twopi = 2.0_dp*pi
          radeg = pi/180.0_dp
          tmu = 2.0e-07_dp
        end subroutine
     end module global_constants

      ! Error control and write out error messages consistently.
      module error_control
        use set_kinds
        integer :: currrank=-1 ! some invalid value
        real(dp) :: currtime=-1.0 ! some invalid value
        public :: errctrl_msg

      contains
        subroutine errctrl_setstate(currrank0,currtime0)
          integer, intent(in) :: currrank0
          real(dp), intent(in) :: currtime0
          currrank = currrank0
          currtime = currtime0
        end subroutine

        subroutine errctrl_msg(subrstr,msgstr,mtype0)
          character(len=*), intent(in) :: subrstr,msgstr
          integer, optional :: mtype0
          integer :: mtype
          character(len=16) :: labelstr

          mtype = 1
          if (present(mtype0)) mtype = mtype0

          select case(mtype)
          case (1)
            labelstr = 'ERROR'
          case (2)
            labelstr = 'WARNING'
          case (3)
            labelstr = 'INFO'
          case default
            labelstr = 'ERROR'
          end select

          write(*, '(a,a,a,a,i3,a,i6,a,a)') trim(labelstr),' in ',subrstr,' at r=',currrank,', t=',int(currtime),': ',msgstr

          open(unit=40,file='errfil.out',status='unknown',access='append')
          write(40,'(a,a,a,a,i3,a,i6,a,a)') trim(labelstr),' in ',subrstr,' at r=',currrank,', t=',int(currtime),': ',msgstr
          close(unit=40)

        end subroutine
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
      use set_kinds
      integer*4 iopen,ifread,itcoil,ifcurr 
      real*8 errbry,xncoil,relax
      data relax/1.00/,errbry/1.0e-04_dp/
      data iopen/0/
      end module var_iopen
!var_zcntrl
       module var_zcntrl
       use set_kinds
       integer*4 isetfb,ishiftz,ioffz,ioffr,idplace,lring 
       real*8 gain,gainp,delzmm
       data gain/.22_dp/,gainp/0.75_dp/,ishiftz/0/,ioffz/7/, &
           ioffr/-7/,lring/0/
       end module var_zcntrl
!var_updown
       module var_updown
       real*8 dpsip,dpsip_last
       logical vfeed,symmetrize,backaverage
       data vfeed/.false./
       end module var_updown

!var_test
      module var_test
      real*8 zcontr,zcurnow
      end module var_test
!var_graphic
      module var_graphic
      integer*4 ivnit,n_write
      data ivnit/35/,n_write/1/
      end module var_graphic

!var_errslop
       module var_errslop
       use set_kinds
       real*8 aaslop,drslop
       data aaslop/0.6_dp/,drslop/0.003_dp/
       end module var_errslop

!var_fitsiref
       module var_fitsiref
       use set_kinds
       real*8 csiref,saisref,fwtref,scalesir
       logical fitsiref
       data scalesir/1.0e-3_dp/,fitsiref/.false./
       end module var_fitsiref

!var_cnnn
       module var_cnnn
       integer*4 n111,n222,n333,n444,n555,n666
       real*8 x000,x111
       data n111/1/,n222/2/,x111/1.0/,x000/0.0/,n333/3/,n444/4/,n555/5/ &
           ,n666/6/
       end module var_cnnn

!var_pcsys
       module var_pcsys
       integer*4 use_alternate_pointnames
       character*80 alternate_pointname_file
        logical*4          do_spline_fit
      data do_spline_fit/.true./
      data use_alternate_pointnames/0/
      data alternate_pointname_file/'/link/efit/pcsnames.dat'/

       end module var_pcsys


!var_pfedge
      module var_pfedge
      use set_kinds
      integer*4 kedgep,kedgef 
      real*8 pedge,pe_psin,pe_width, &
      f2edge,fe_psin,fe_width,constf2, &
      tpedge,tfedge,rdlcfe,rqape,rqafe,betped,betnped
      data kedgep/0/,pe_width/0.02_dp/,pe_psin/0.98_dp/,pedge/0.0/, &
           kedgef/0/,fe_width/0.02_dp/,fe_psin/0.98_dp/,f2edge/0.0/

      end module var_pfedge



!var_sxpoint
       module var_sxpoint
       real*8 sissep,rssep,zssep,sifsep,rfsep,zfsep,rminss, &
              zminss,rmaxss,zmaxss,rminfs,rmaxfs,zminfs,zmaxfs
       end module var_sxpoint

!var_consta
      module var_consta
        use set_kinds
        public
        real*8 :: tmu2,errcut,tmu0,tmu02
        integer*4 :: ibunmn,kinput,kcaldia=0
     end module var_consta

!var_rcfact
     module var_rcfact
     use set_kinds
     integer*4 ircfact
     data ircfact/0/
     end module var_rcfact
!var_curpo
     module var_curpro
     use set_kinds
     real*8 emf,emp,enf,enp,rbetap,rzero,pbetap,qenp,qemp,qenf

!vas      common/curpro/emf,emp,enf,enp,rbetap,rzero,pbetap,qenp,qemp,qenf
      data rzero/1.6955_dp/
     end module var_curpro

!var_pfterm
     module var_pfterm

     integer*4 kffcur,kppcur,kpcurn,icalbet 
     real*8    chidlc,gammaf,gammap,cstabz ,cstab0 &
            ,vbtot2,vbtvac2,vbtor2,vbtvac,vbeta0 &
            ,vbtmag,btvvac2,btvtor2,btvtot2
     data cstabz/0.0e-13/
     data icalbet/1/
     end module var_pfterm

!var_nio
     module var_nio
     integer*4 nin,nout,ntty,nrsppc,nrspfc,nttyo,neqdsk,nffile,nsave
     integer*4 nsnapf
     character*2 appendsnap
     character*100 snapfile, tmpdata, snapextin
     data nin/11/,nout/10/,ntty/5/nrsppc/25/,nrspfc/26/, &
     neqdsk/38/,nffile/40/,nsave/49/,nttyo/6/
     end module var_nio

!var_cfit
     module var_cfit
     use set_kinds
     integer*4 mxiter,idone,nitera,nxiter,ixnn,isolve  
     real*8   error,errorm,errmin,delerr,delerb
     data errmin/0.010_dp/,errorm/10./
     end module var_cfit 
!var_cgrid
     module var_cgrid
     real*8 darea,drgrid,dzgrid,qmaxis,cratio,dfsqe,cratiof
     data dfsqe/0.0/
     end module var_cgrid
!var_extra
     module var_extra
      real*8 scrape,tolbndpsi
      integer*4 nextra,ixstrt,iextra,iprobe,ifcoil,iecoil &
           ,iexcal,iconsi,iqplot,klabel,kthkcrv,ifindopt
      data ifcoil/1/,kthkcrv/0/,klabel/0/
     end module var_extra
!var_conveg
     module var_conveg
     use set_kinds
     real*8 omega,relip,zelip,aelip,eelip,errorq,omecur 
     integer*4 jjmax
     data relip/1.68_dp/,zelip/0.0/,aelip/0.60_dp/eelip/1.2_dp/,&
     omega/1.0/,errorq/1.0e-03_dp/
     data jjmax/1/
     end module var_conveg

!var_limmm
     module var_limmm
     real*8 xlmin,xlmax,ylmin,ylmax
     integer*4 limid,limup,limbot
     data limid/33/
     end module var_limmm

!var_inaver
     module var_inaver
     
     integer*4 iavem,iaved,iavev,iaveus
     data iaveus/0/

     end module var_inaver
!
!var_vessel
     module var_vessel
     use eparmdud129,only:nesum,nvesel,nfcoil
     real*8,dimension(nesum) :: volecs,volecc,rsisec
     real*8,dimension(nfcoil) :: volfcs,volfcc 
     real*8,dimension(nvesel) :: rvs,zvs,hvs,wvs,avs,avs2,rsisvs 
     real*8 powvs,pvscur,pscurn,ppscur,efreq,sumvs0
     integer*4 ivesel

!      common/vessel/ivesel,volecs(nesum),volecc(nesum),rsisvs(nvesel) & 
!        ,efreq,sumvs0,volfcs(nfcoil),volfcc(nfcoil) &  
!        ,rsisec(nesum),powvs,pvscur,pscurn,ppscur &  
!        ,rvs(nvesel),zvs(nvesel),hvs(nvesel),wvs(nvesel), & 
!        avs(nvesel),avs2(nvesel) 
     end module var_vessel
!var_cyclic_red
     module var_cyclic_red
     use eparmdud129,only:nw,nwnh,ncurrt,icycred_loopmax
     real*8,dimension(:,:),allocatable :: beti,abeti,wk1
     real*8,dimension(icycred_loopmax) :: alphab,diag1 
     real*8,dimension(:),allocatable :: rhsdumy1
     real*8,dimension(:),allocatable :: phi,v,wk2,diagl,diagu
     real*8,dimension(ncurrt) :: tempgrid,tempgrid2
     real*8 diag,rhs_a_dumy,rhs_b_dumy
     integer*4 nhpwr
       
!vas      common/cyclic_red/beti(icycred_loopmax,nw-2), &
!vas       abeti(icycred_loopmax,nw-2),wk1(icycred_loopmax,nw-2), &
!vas       alphab(icycred_loopmax),diag1(icycred_loopmax), &
!vas       rhsdumy1(nwnh),phi(nw),v(nw),wk2(nw),diagl(nw),diagu(nw), &
!vas       tempgrid(ncurrt),tempgrid2(ncurrt),diag,rhs_a_dumy,rhs_b_dumy, &
     end module var_cyclic_red

!var_scalem
     module var_scalem
     use eparmdud129,only:nrsmat,mfnpcr
     integer*4 infosc
!vas     real*8,dimension(:), allocatable :: rowscale
!vas     real*8,dimension(:), allocatable :: colscale
     real*8 rowcnd,colcnd,arspmax
     real*8,dimension(nrsmat) :: rowscale
     real*8,dimension(mfnpcr) :: colscale
     logical scalea
     data scalea/.false./
     end module var_scalem

     module var_solove

      integer*4 islve
      real*8 salpha,sbeta,srm,scc1,seee,saaa,srma, &
                    sbetaw

     end module var_solove 
     module var_bunemn

      integer*4 :: mno,m,n
      integer*4 :: nbmdim,nww,nhh
      real*8    :: drdz2,rgrid1,delrgrid,delz
      real*8    :: s,shift,dr,dz
!vas      real(kind=dp)    :: drdz2,rgrid1,delrgrid,delz
!vas      real(kind=dp)    :: s,shift,dr,dz

!      mno = nbmdim
!      m = nww
!      n = nhh
!      s = drdz2
!      shift = rgrid1
!      dr = delrgrid
!      dz = delz

      end module var_bunemn

!------ put all the remining common blocks into modules here
!var_contor
      module var_contor
      use eparmdud129,only:ntime
      real*8,dimension(ntime) :: s1,s2,s3,bpolav
!vas      common/contor/s1(ntime),s2(ntime),s3(ntime),bpolav(ntime)
      end module var_contor
!var_mfield
      module var_mfield
      use eparmdud129,only:npoint
      real*8,dimension(:),allocatable :: bpol,plengt,bpolz
      real*8 siar,siaz
!vas      common/mfield/bpol(npoint),plengt(npoint),bpolz(npoint),siar,siaz
      end module var_mfield
!var_hist
      module var_hist
      use eparmdud129,only:ntime
      integer, dimension(ntime) :: jerror
      real*8, dimension(ntime) :: eout,rout,zout,doutu &
        ,doutl,aout,vout,betat,otop &
        ,betap,ali,oleft,oright,qsta &
        ,rcurrt,zcurrt,qout,olefs &
        ,orighs,otops,sibdry,areao &
        ,wplasm,elongm,qqmagx,terror &
        ,rmagx,zmagx,obott,obots &
        ,alpha,rttt,dbpli,delbp,oring &
        , sepexp,shearb &
        ,xtch,ytch,qpsib,vertn,aaq1 &
        ,aaq2,aaq3,btaxp,btaxv &
        ,simagx,seplim,peak &
        ,wbpol,taumhd,betapd,betatd &
        ,alid,wplasmd,taudia,wbpold &
        ,qmerci,slantu,slantl,zeff, &
        zeffr,tave,rvsin,zvsin, &
        rvsout,zvsout,wpdot,wbdot, &
        vsurfa,cjor95,pp95,ssep, &
        yyy2,xnnc,wtherm,wfbeam,taujd3,tauthn, &
        ali3,tflux,twagap
        real*8,dimension(2,ntime) :: rseps,zseps
        real*8 pasman,betatn,psiq1,betat2
        integer jtwagap

!vas      common/hist/eout(ntime),rout(ntime),zout(ntime),doutu(ntime) &
!vas        ,doutl(ntime),aout(ntime),vout(ntime),betat(ntime),otop(ntime) &
!vas        ,betap(ntime),ali(ntime),oleft(ntime),oright(ntime),qsta(ntime) &
!vas        ,rcurrt(ntime),zcurrt(ntime),qout(ntime),olefs(ntime) &
!vas        ,orighs(ntime),otops(ntime),sibdry(ntime),areao(ntime) &
!vas        ,wplasm(ntime),elongm(ntime),qqmagx(ntime),terror(ntime) &
!vas        ,rmagx(ntime),zmagx(ntime),obott(ntime),obots(ntime) &
!vas        ,alpha(ntime),rttt(ntime),dbpli(ntime),delbp(ntime),oring(ntime) &
!vas        ,rseps(2,ntime),zseps(2,ntime),sepexp(ntime),shearb(ntime) &
!vas        ,xtch(ntime),ytch(ntime),qpsib(ntime),vertn(ntime),aaq1(ntime) &
!vas        ,aaq2(ntime),aaq3(ntime),btaxp(ntime),btaxv(ntime) &
!vas        ,simagx(ntime),jerror(ntime),seplim(ntime),peak(ntime) &
!vas        ,wbpol(ntime),taumhd(ntime),betapd(ntime),betatd(ntime) &
!vas        ,alid(ntime),wplasmd(ntime),taudia(ntime),wbpold(ntime) &
!vas        ,qmerci(ntime),slantu(ntime),slantl(ntime),zeff(ntime), &
!vas        zeffr(ntime),tave(ntime),rvsin(ntime),zvsin(ntime), &
!vas        rvsout(ntime),zvsout(ntime),wpdot(ntime),wbdot(ntime), &
!vas        vsurfa(ntime),cjor95(ntime),pp95(ntime),ssep(ntime), &
!vas        yyy2(ntime),xnnc(ntime),pasman,betatn,psiq1,betat2, &
!vas        wtherm(ntime),wfbeam(ntime),taujd3(ntime),tauthn(ntime), &
!vas        ali3(ntime)
      end module var_hist
!var_hist2
      module var_hist2
      use set_kinds
      use eparmdud129,only:ntime
      real*8,dimension(ntime) :: qsiwant,cjorsw,cjor0, &
        ssiwant,ssi95,cjor99,cj1ave &
        ,rmidin,rmidout,psurfa
      real*8 psiwant,rexpan,fexpan,qqmin,fexpvs,shearc &
        ,sepnose,ssi01,znose,rqqmin
      data psiwant/1.0/,rexpan/0.010_dp/,znose/-1.276_dp/
!vas      common/hist2/psiwant,qsiwant(ntime),cjorsw(ntime),cjor0(ntime), &
!vas        ssiwant(ntime),ssi95(ntime),rexpan,fexpan,qqmin,fexpvs,shearc &
!vas        ,sepnose,ssi01,znose,rqqmin,cjor99(ntime),cj1ave(ntime) &
!vas        ,rmidin(ntime),rmidout(ntime),psurfa(ntime)
      end module var_hist2
!var_cshape
      module var_cshape
      use set_kinds
      use eparmdud129,only:npoint
      real*8,dimension(:),allocatable :: xout,yout
      real*8 dpsi,rymin,rymax, &
        zxmin,zxmax,xmin,xmax,ymin,ymax,rmaxis,zmaxis, emaxis, &
        rminzm,rmaxzm,dismins,simins,delrmax1,delrmax2
      integer*4 nfound
      data emaxis/1.3_dp/
      end module var_cshape
!var_divdis
      module var_divdis
      use set_kinds
      use eparmdud129,only:ntime
      real*8,dimension(ntime) :: dolubaf,dolubafm,diludom, &
        diludomm,dminux,dminlx, &
        ratsol,rvsiu,zvsiu,rvsou, &
        zvsou,rvsid,zvsid,rvsod, &
        zvsod
        real*8 rubaf,zubaf,rlbaf,zlbaf,rudom,zudom
      data rubaf/1.372_dp/,rudom/1.0420_dp/,rlbaf/1.6810_dp/
      data zubaf/1.310_dp/,zudom/1.1624_dp/,zlbaf/-1.339_dp/

!vas      common/divdis/dolubaf(ntime),dolubafm(ntime),diludom(ntime), &
!vas        diludomm(ntime),dminux(ntime),dminlx(ntime), &
!vas        ratsol(ntime),rvsiu(ntime),zvsiu(ntime),rvsou(ntime), &
!vas        zvsou(ntime),rvsid(ntime),zvsid(ntime),rvsod(ntime), &
!vas        zvsod(ntime), &
!vas        rubaf,zubaf,rlbaf,zlbaf,rudom,zudom
      end module var_divdis
!var_cpsi
      module var_cpsi
      use eparmdud129,only:nwnh
      real*8,dimension(:),allocatable :: psi,xpsi,vfbrrt,psipla
      real*8 vcurfb(3)
      real*8 psibry,simag,sidif,eouter,zplasm,zpwant,vertfb,difpsi &
             ,cupdown
      data vertfb/0./,cupdown/-100000./ 
!vas      common/cpsi/psi(nwnh),psibry,simag,sidif,xpsi(nwnh),eouter, &
!vas                  zplasm,zpwant,vertfb,vcurfb(3),vfbrrt(nwnh),difpsi &
!vas                 ,cupdown,psipla(nwnh)
      end module var_cpsi
!var_cvalue
      module var_cvalue
      use eparmdud129,only:nsilop,ntime,nrogow,magpri,nfcoil,nacoil
      real*8,dimension(nsilop,ntime) :: csilop,csilopv 
      real*8,dimension(nrogow,ntime) :: crogow
      real*8,dimension(magpri,ntime) :: cmpr2,cmpr2v 
      real*8,dimension(ntime) :: cpasma,xndnt 
      real*8,dimension(nfcoil,ntime) :: ccbrsp
      real*8,dimension(ntime,nacoil) :: caccurt
      real*8 cbetap,cli,cqqxis,cbetat,ci0,cipmp2 
!vas      common/cvalue/csilop(nsilop,ntime),crogow(nrogow,ntime), &
!vas        cmpr2(magpri,ntime),cpasma(ntime),xndnt(ntime) &
!vas       ,cbetap,cli,cqqxis,cbetat,ci0,cipmp2 &
!vas       ,ccbrsp(nfcoil,ntime),caccurt(ntime,nacoil) &
!vas       ,csilopv(nsilop,ntime),cmpr2v(magpri,ntime)
      end module var_cvalue



!-- modules from ecomdu2
!var_gtable
      module var_gtable
      use eparmdud129,only:nsilop,magpri,nfcoil,nw,nh,nwnh,nbwork,mbdry
!vas      real*8,dimension(nwnh,nfcoil) :: gridfc
      real*8,dimension(:,:),allocatable,save :: gridfc
!vas      real*8,dimension(nw) :: rgrid
!vas      real*8,dimension(nh) :: zgrid
      real*8,dimension(:),allocatable,save :: rgrid
      real*8,dimension(:),allocatable,save :: zgrid
!vas      real*8,dimension(nsilop,nfcoil) :: rsilfc
!vas      real*8,dimension(magpri,nfcoil) :: rmp2fc
      real*8,dimension(:,:),allocatable,save :: rsilfc
      real*8,dimension(:,:),allocatable,save :: rmp2fc
!vas      real*8,dimension(nwnh,nw) :: gridpc
      real*8,dimension(:,:),allocatable,save :: gridpc
!vas      real*8,dimension(nsilop,nwnh) :: gsilpc
!vas      real*8,dimension(magpri,nwnh)  :: gmp2pc
      real*8,dimension(:,:),allocatable,save :: gsilpc
      real*8,dimension(:,:),allocatable,save :: gmp2pc
!rem by vas      real*8,dimension(nbwork,nwnh)  :: gwork
      real*8,dimension(nfcoil,nfcoil) :: rfcfc
!      real*8,dimension(mbdry*nwnh-magpri*nwnh-nbwork*nwnh-nfcoil*nfcoil) &
!      :: cjfgtable
!      real*8,dimension(:) :: cjfgtable
      integer*4 iallocate_stat
!vas      common/gtable/gridfc(nwnh,nfcoil),rgrid(nw),zgrid(nh) &
!vas           ,rsilfc(nsilop,nfcoil),rmp2fc(magpri,nfcoil) &
!vas           ,gridpc(nwnh,nw),gsilpc(nsilop,nwnh) &
!vas           ,gmp2pc(magpri,nwnh),gwork(nbwork,nwnh) &
!vas           ,rfcfc(nfcoil,nfcoil) &
!vas           ,cjfgtable(mbdry*nwnh-magpri*nwnh-nbwork*nwnh-nfcoil*nfcoil)
      end module var_gtable

! jm.s
! NOTE : array sizes are grid size-dependent so they should be dynamically allocated, but
!        we cannot make them dynamic since they are included in a namelist
! NOTE : we assume the largest usefule grid size is 1025
! NOTE : npsi_ext (actual dimension of _ext arrays) used in code logic and intentially set
!        to default value of -1
      module profile_ext_mod
        use eparmdud129,only:mbdry,nlimit
        integer :: npsi_ext=-1,nbdry_ext,limitr_ext
        real*8,dimension(1025) :: pprime_ext
        real*8,dimension(1025) :: ffprim_ext
        real*8,dimension(1025) :: psin_ext
        real*8,dimension(1025) :: qpsi_ext
        real*8,dimension(1025) :: bpp_ext, cpp_ext, dpp_ext
        real*8,dimension(1025) :: bfp_ext, cfp_ext, dfp_ext
        real*8 :: rbdry_ext(mbdry),zbdry_ext(mbdry),xlim_ext(nlimit),ylim_ext(nlimit)
        real*8 :: sign_ext, scalepp_ext, scaleffp_ext, cratio_ext, cratiop_ext, &
                  cratiof_ext,psirz_ext
        character*80 :: geqdsk_ext
        logical :: fixpp = .false.
      end module profile_ext_mod
! jm.e

! NOTE : keep track of times for which BCOIL and ECOIL data exist (see getecdud129.f90)
      module vtime_mod
         use eparmdud129,only:ntime
         integer :: nvtime = -1
         real*8,dimension(ntime) :: vtime
      end module vtime_mod
