      module var_acoilrz
      real*8,dimension(:),allocatable  :: racoil,zacoil,wacoil,hacoil
      end module var_acoilrz

      module var_ccurrn
      !use eparm,only:nwnh
      real*8 fconst,betap0
      real*8,dimension(:),allocatable :: pcurrt,pcurrw,pcurrtpp
      integer*4 icurrt,icinit,isicinit,icprof
      end module var_ccurrn

      module var_czero
      real*8,dimension(:),allocatable :: zero,www
      integer*4 iweigh
      end module var_czero

      module var_parame
      real*8,dimension(:),allocatable :: volp,pprime,pres,ffprim,fpol, &
        qpsi,r2surf,rpres,curmid,cjor,bpolss, &
        rmajz0,bprmaj,btrmaj,r1surf,r2surg, &
        bbfpol,ccfpol,ddfpol
      real*8 r1bdry,r2bdry,carea,jwantm,r2wdry,r4wdry,rpwdry,rp2wdry 
      real*8,dimension(:),allocatable :: r1sdry,r2sdry,vzeroj,sizeroj 
      end module var_parame

      module var_cqfit
      real*8 fwtqa,qvfit,fbrdy,fcentr,cjmaxi, &
             rqajtor,rqaftor,rqapetor,rqafetor
      integer*4 ivacum 
      real*8,dimension(:),allocatable :: rqajx,rjjjx
      real*8,dimension(:),allocatable :: rqafx,rjjfx
      real*8,dimension(:),allocatable :: rqawx,rjjwx
      end module var_cqfit

      module var_rmatri
      real*8,dimension(:),allocatable :: brsp
      real*8,dimension(:),allocatable :: chisq
      real*8 cond,tchifcc
      integer*4 nparam,nfnpcr,nfnwcr,nbase
      end module var_rmatri
  
      module var_chimag
      real*8,dimension(:),allocatable :: chi2rm
      end module var_chimag

      module var_exdata
      integer*4 ishot,itime,ifitvs,iacoil,itimeu,kersil
      real*8 serror,fwtcur,elomin,fwtbp,fwtdlc,errsil,rcentr
      real*8,dimension(:),allocatable :: fwtsi,rwtsi
      real*8,dimension(:),allocatable :: fwtmp2,rwtmp2
      real*8,dimension(:),allocatable :: fwacoil
      real*8,dimension(:),allocatable :: fwtfc
      real*8,dimension(:),allocatable :: cbrspv
      real*8,dimension(:),allocatable :: bcentr
      end module var_exdata

      module var_texdat
      integer*4 ierpla,ico2,ierrdi 
      real*8 bitip
      real*8,dimension(:,:),allocatable :: silopt
      real*8,dimension(:,:),allocatable :: expmpi
      real*8,dimension(:,:),allocatable :: accurt
      real*8,dimension(:,:),allocatable :: fccurt
      real*8,dimension(:,:),allocatable :: eccurt
      real*8,dimension(:,:),allocatable :: denvt
      real*8,dimension(:,:),allocatable :: denrt
      integer*4,dimension(:,:),allocatable :: iccurt
      real*8,dimension(:),allocatable :: ipmeas,time,pbinj 
      integer*4,dimension(:),allocatable :: ierpsi
      integer*4,dimension(:),allocatable :: ierpr
      integer*4,dimension(:),allocatable :: iermpi
      real*8,dimension(:),allocatable :: psibit
      real*8,dimension(:),allocatable :: prbit
      real*8,dimension(:),allocatable :: bitmpi
      real*8 :: vbit
           
      real*8,dimension(:),allocatable :: vloopt,psiref,diamag,sigdia
      real*8,dimension(:),allocatable :: bitfc
      integer*4,dimension(:),allocatable :: ierfc
      real*8,dimension(:),allocatable :: bitec
      integer*4,dimension(:),allocatable :: ierec
      real*8,dimension(:),allocatable ::bitic
      integer*4,dimension(:),allocatable ::ieric
      integer*4 ierdia(3)
      end module var_texdat

      module var_savfit
      real*8,dimension(:),allocatable ::  swtsi
      real*8,dimension(:),allocatable ::  swtmp2
      real*8,dimension(:),allocatable ::  swtfc
      real*8,dimension(:),allocatable ::  swtec
      real*8 swtcur,swtdlc
      end module var_savfit

      module var_cstark
      real*8,dimension(:,:),allocatable :: tangam,siggam,a1gam,a2gam, &
                                           a3gam,a4gam,tangam_uncor
      real*8,dimension(:),allocatable  :: fwtgam,chigam,swtgam
      real*8 v30lt,v30rt,v210rt,v210lt
      integer*4,dimension(:),allocatable :: iergam
      integer*4,dimension(:),allocatable :: mseport,mse_spave_on
      integer*4  kstark,iplots,kmtark,klibim,kdomse, &
                 msebkp,msefitfun,kwaitmse,mse_quiet, &
                 mse_strict,ok_30rt,ok_210lt,mse_usecer, &
                 mse_certree,mse_use_cer330,mse_use_cer210
      real*8 chigamt,chilibt,dtmsefull,t_max_beam_off 
      real*8,dimension(:,:),allocatable :: rbrpc,rbzpc,rgampc
      real*8,dimension(:,:),allocatable :: rbrfc,rbzfc
      real*8,dimension(:,:),allocatable :: rbrec,rbzec,rgamec
      real*8,dimension(:,:),allocatable :: rbrvs,rbzvs,rgamvs
      real*8,dimension(:,:),allocatable :: rgamfc
      real*8,dimension(:,:),allocatable :: rhsgam,rrgam,zzgam, &
        starkar,starkaz,a5gam,a6gam,a7gam,a8gam
      real*8,dimension(:,:),allocatable :: cmgam,spatial_fix
      real*8,dimension(:,:),allocatable :: rgamac,rbrac,rbzac
      real*8,dimension(:),allocatable :: btgam,sistark,qstark, &
        rmse_gain,rmse_slope,rmse_scale,rmse_offset,sigam, &
        bzmse,bzmsec,cjmse,cjmsec,rhogam
      real*8,dimension(:,:,:,:),allocatable :: spatial_avg_gam
      end module var_cstark

      module var_msels
      real*8,dimension(:,:),allocatable :: bmselt,sbmselt,fwtbmselt, &
        rrmselt,zzmselt,l1mselt,l2mselt,l4mselt,emselt,semselt, &
        fwtemselt,swtbmselt,swtemselt,cmmls,cmels,rhsmls,rhsels, &
        l3mselt,cmmls2,cmmlsv
      real*8,dimension(:),allocatable :: swtbmsels,swtemsels,chimls,chiels
      real*8 :: avemsels, fmlscut, tchimls, tchiels
      real*8,dimension(:,:),allocatable :: rmlspc
      real*8,dimension(:,:),allocatable :: rmlsec,relsec
      real*8,dimension(:,:),allocatable :: rmlsvs,relsvs
      real*8,dimension(:),allocatable ::  simls,sinmls,btmls,brmls,bzmls
      real*8,dimension(:),allocatable ::  chi2mls, chi2gamt

      integer*4,dimension(:,:),allocatable :: iermselt
      integer*4,dimension(:),allocatable :: kerrot
      integer*4 kdomsels, mmbmsels
      character*3 :: synmsels
      end module var_msels

      module var_rpedge
      real*8,dimension(:),allocatable :: rsilpe
      real*8,dimension(:),allocatable :: rmp2pe
      real*8 fgowpe
      real*8,dimension(:),allocatable :: rprepe
      real*8,dimension(:),allocatable ::  gbdrpe
      real*8,dimension(:),allocatable :: rgampe,rbrpe,rbzpe
      real*8,dimension(:),allocatable ::  rmlspe,relspe
      end module var_rpedge

      module var_rfedge
      real*8,dimension(:),allocatable :: rsilfe
      real*8,dimension(:),allocatable :: rmp2fe
      real*8 fgowfe
      real*8,dimension(:),allocatable :: gbdrfe
      real*8,dimension(:),allocatable :: rgamfe,rbrfe,rbzfe
      real*8,dimension(:),allocatable :: rmlsfe,relsfe
      end module var_rfedge

      module var_cece
      integer*4 ierecebz,necein,jo,kfixro,nece,kece,kecebz,mecein, &
                kfitece,kinputece,kcallece,nharm,kdoece, &
                mtxece,kfixrece,nfit,kcmin,nconstr,iwo
      character*4 eceiter
      real*8 receo,fwtecebz,fwtecebz0,rteo,zteo,robit,swtecebz,chiecebz, &
             ecebzfit,ecebzbit,fwtnow,zeceo,chisqfit,xfit(10),tchiece, &
             recebzdz,gecebzdz,eceerror
      real*8,dimension(:),allocatable :: teecein,feece,errorece, &
                                         teecein0,feece0,errorece0, &
                                         becein,recein,teeceinr
      real*8,dimension(:),allocatable :: recem,recep,fwtece,fwtece0, &
                                         swtece,chiece,ecefit,ecebit
      real*8,dimension(:),allocatable :: recebzfc
      real*8,dimension(:),allocatable :: gecebzpc
      real*8,dimension(:),allocatable :: recebzec
      real*8,dimension(:,:),allocatable  :: recefc
      real*8,dimension(:,:),allocatable :: gecepc,geceppc,gecempc
      real*8,dimension(:,:),allocatable  :: receec
      integer*4,dimension(:),allocatable :: ierece,iwp,iwm
      real*8,dimension(:,:),allocatable :: recepc
      real*8,dimension(:,:),allocatable :: brspece
      real*8,dimension(:),allocatable :: recebzpc
      real*8,dimension(:),allocatable :: brspecebz,cmecebz
      real*8,dimension(:,:),allocatable :: cmece
      real*8,dimension(:,:),allocatable :: recevs
      real*8,dimension(:),allocatable :: recebzvs
      real*8,dimension(:),allocatable :: recedz,gecedz,rtep,rtem,rpbit,rmbit
      real*8,dimension(:,:),allocatable :: receac
      real*8,dimension(:),allocatable  :: recebzac
      real*8,dimension(:),allocatable  :: teecer,rrr,bbf,teeceb
      real*8,dimension(:),allocatable :: receoi
      real*8,dimension(:,:),allocatable :: recemi,recepi
      end module var_cece

      module var_calchi
      real*8 chipre,chifin,chitot,chipasma
      real*8,dimension(:),allocatable :: saisil
      real*8,dimension(:),allocatable :: saimpi
      real*8,dimension(:),allocatable :: saipre
      real*8,dimension(:),allocatable :: chifcc
      end module var_calchi

      module var_dlc
      real*8 sigdlc
      real*8,dimension(:),allocatable :: dfluxc,cdflux,edflux
      real*8,dimension(:),allocatable :: rspdlc
      end module var_dlc
    
      module var_comco2
      real*8,dimension(:,:),allocatable :: rco2r
      real*8,dimension(:,:),allocatable :: rco2v
      real*8,dimension(:),allocatable :: chordv
      real*8,dimension(:),allocatable :: chordr
      real*8,dimension(:,:),allocatable :: dco2r
      real*8,dimension(:,:),allocatable :: dco2v
      end module var_comco2

      module var_check
      !use eparm,only:ntime
      integer*4,dimension(:,:),allocatable ::erflag
      integer*4 lflag
      integer*4, parameter :: nflag=21
      end module var_check

      module var_consum
      real*8 condno,condin 
      real*8,dimension(:),allocatable :: cerror,csibry,csimag, &
        cvolp,crmaxi,czmaxi,cemaxi,cqmaxi,cchisq,brfbc,tvfbrt,cdelz
      end module var_consum

      module var_cxray
      real*8,dimension(:),allocatable :: rxray,zxray,xangle
      integer*4 ksxr0(10),ksxr2(10),idosxr
      end module var_cxray

      module var_mercie
      integer*4 imerci
      real*8,dimension(:),allocatable :: rzzmax,zzmax
      end module var_mercie

      module opt_input
         integer*4 mode_in
         character cmdfile_in*15, shotfile_in*15
         integer*4 shot_in
         real*8 starttime_in
         real*8 deltatime_in
         integer*4 steps_in
         logical use_opt_input
         character snapext_in*86
         character(80),dimension(:),allocatable :: inpfile
      end module opt_input

      module mpi_info
         integer*4 :: rank, nproc, ierr
         integer*4,dimension(:),allocatable :: dist_data
         integer*4,dimension(:),allocatable :: dist_data_displs
      end module mpi_info

      module var_input
      logical write_kfile,fitfcsum,use_previous, &
              req_valid,req_valid_prior
      integer*4 icondn,itek,kdata,itrace,ierchk,ierchk_prior,iconvr, &
                ixray,itell,kprfit,ibound,ibatch,idite,ilaser,lookfw, &
                kdot,icutfp,keqdsk,kdofit,kbetapr,kplotpr,kpressb, &
                kfcurb,kpcurb,kzeroj,ncstne,ncstte, &
                kwripre,negcur,kframe,kskipvs,icntour,iavdpl, &
                limvs,kbound,kgraph,istore,iout,iout_prior,kdopre, &
                iishot,kktime,iplcout,iplcout_prior,ksigma,kwritime
      integer*4 iteks,mxiters,n1coils,ierchks
      integer*4 itekt,mxitert,n1coilt
      real*8 zelipss,zeliptt
      real*8 cutip,dtdot,xpsimin,fcurbd,pcurbd,prbdry,sgprmin, &
             prespb,tipbry,tepbry,dipbry,depbry,pbimpb,sigppb,sigpreb, &
             sigtipb,sigtepb,sigdipb,sigdepb,fwtpd,fwtfd,cstabte,cstabne, &
             dp1dxf,sgtimin,alphafp,xpsialp,vsdamp,siloplim, &
             rminvs,rmaxvs,zminvs,zmaxvs,zbound,yvs2,saimin, &
             fztor,fzpol,tcurrp,fpolvs,rbound,sigprebi,pressbi, &
             alphamu,saicon,rsepex,timeb,dtime
      real*8,dimension(:),allocatable :: rzeroj
      real*8,dimension(:),allocatable :: fwtpre
      real*8,dimension(:),allocatable :: vforcep,vforcet
      real*8,dimension(:),allocatable :: fwtfcsum
      character*82 snap_file
      character(80),dimension(:),allocatable :: ifname
      end module var_input

      module var_inputc
      character*12 mfitpop
      character(4),dimension(:),allocatable :: limloc
      character(10),dimension(:),allocatable :: vsname
      character(10),dimension(:),allocatable :: mpnam2
      character(10),dimension(:),allocatable :: lpname
      character filimt*100,cshot*6,jdebug*4
      integer*4 idebug
      end module var_inputc

      module var_switch
      integer*4 nqaxis,isumip,jbeta,jli,nqwant
      real*8 sumip,fbetap,fli,fbetat,fbetan 
      real*8,dimension(:),allocatable :: qsiw, pasmsw,fqsiw,siwantq
      real*8,dimension(:,:),allocatable :: fgowsw
      end module var_switch

      module var_siloop
      real*8,dimension(:),allocatable :: rsi,zsi,wsi,hsi
      end module var_siloop

      module var_slname
      real*8,dimension(:),allocatable :: as,as2
      integer*4 nslref
      end module var_slname

      module var_ecoil 
      !use eparm,only:necoil
      real*8,dimension(:),allocatable :: re,ze,we,he
      integer*4,dimension(:),allocatable :: ecid
      end module var_ecoil 

      module var_fcoil
      real*8,dimension(:),allocatable :: rf,zf,wf,hf,af,af2, &
                                         turnfc,fcturn
      integer*4,dimension(:),allocatable :: fcid
      end module var_fcoil 

      module var_mprobe
      real*8,dimension(:),allocatable :: xmp2,ymp2,amp2,smp2,patmp2
      end module var_mprobe

      module var_limite
      integer*4 limitr,iplim,limfag,limitr_180
      real*8,dimension(:),allocatable :: xlim,ylim,xlim_180,ylim_180
      real*8 yltype,yltype_180
      end module var_limite

      module var_mimite
      !use eparm,only:nlimit
      integer*4 mimitr,mimitr_180
      real*8,dimension(:),allocatable :: xmim,ymim,xmim_180,ymim_180
      end module var_mimite

      module var_udata
      integer*4,dimension(:),allocatable :: ipsi,irogw,imag2,iplasm, &
                                            idlopc,ifc,iec
      end module var_udata

      module var_morsum
      real*8,dimension(:),allocatable :: csumip,tratio,aveerr
      integer*4,dimension(:),allocatable :: iermax,jermax
      end module var_morsum

      module var_bdsend
      integer*4 nbbbs,nbskip,nbdrymx,nbdryp
      real*8,dimension(:),allocatable :: rbbbs,zbbbs
      end module var_bdsend

      module var_fxbry
      use set_kinds, only: dp
      integer*4 nbdry,nsol
      real*8 wsisol
      logical fitts 
      real*8,dimension(:),allocatable :: rbdry,zbdry,fwtbdry,fwtbry, &
                                         sigrbd,sigzbd,rbdry0,zbdry0
      real*8,dimension(:),allocatable :: rsol,zsol,fwtsol,fwtsolw
      real*8,dimension(:,:),allocatable :: rbdrfc, rsolfc
      real*8,dimension(:,:),allocatable :: rbdrac
      real*8,dimension(:,:),allocatable :: rbdrpc,rsolpc
      real*8,dimension(:,:),allocatable :: gbdrpc
      real*8, parameter :: dselsum=1.e-12_dp
      end module var_fxbry

      module var_fwtdz
      logical fitdelz
      integer*4 ndelzon,ifitdelz 
      real*8 errdelz,fgowdz,scaledz,stabdz,relaxdz,cdeljsum
      real*8,dimension(:),allocatable :: gmp2dz
      real*8,dimension(:),allocatable :: gsildz
      real*8,dimension(:),allocatable :: gbrdz,gbzdz,rgamdz
      real*8,dimension(:),allocatable :: rmlsdz,relsdz
      real*8,dimension(:),allocatable :: rpredz,rprwdz
      real*8,dimension(:),allocatable :: gbdrdz
      real*8,dimension(:),allocatable :: rdjdz
      end module var_fwtdz
      
      module var_combry
      real*8,dimension(:),allocatable ::erbloc, erbsloc
      real*8 erbmax,erbave,erbsmax
      end module var_combry

      module var_fbysta
      real*8 cfcoil, psibry0
      integer*4 ifref
      real*8,dimension(:),allocatable :: fcsum,fczero
      end module var_fbysta

      module var_prdata
      integer*4 npress,npteth,nption,npneth,nbeam,ndokin,nbrmcrd, &
                nptef,npnef,nptionf,nmass 
      real*8 pbeamb,pressb,zeffvs,zlowimp
      integer*4,dimension(:),allocatable :: ivbcuse
      real*8,dimension(:),allocatable :: pressr,rpress,zpress,tethom, &
        rteth,zteth,tionex,rion,zion,dnethom,rneth,zneth, &
        pbeam,sibeam,sgteth,sigti,sigpre,sgneth,precal, &
        dnbeam,dmass,scalepr,premea,saipre2
      real*8,dimension(:),allocatable :: pbimth,bremin,bremsig,brmrtan, &
                                         brmzelev,dnbthom
      end module var_prdata

      module var_cerfit
      integer*4 keecur,keefnc,keeknt,needer,keehord
      integer*4,dimension(:),allocatable :: keebdry,kee2bdry
      real*8 ecurbd,eetens
      real*8,dimension(:,:),allocatable :: rgamer
      real*8,dimension(:,:),allocatable :: rmlser,relser
      real*8,dimension(:),allocatable ::eeknt,eebdry,ee2bdry,cerer
      real*8,dimension(:),allocatable :: ermse
      real*8,dimension(:,:),allocatable :: e1rbz,e2rbz,e3rbr
      real*8,dimension(:),allocatable :: ermid,eshear,epoten,rhovn, &
                                         rpmid,xmid,sigrid,sipmid, &
                                         brhovn,crhovn,drhovn,rhopmid
      end module var_cerfit

      module var_ccgama
      integer*4 kcalpa,kcgama,kcomega
      real*8 fwtxx,fwtxxj,fwtxxq,fwtxxb,fwtxli 
      real*8,dimension(:,:),allocatable ::  calpa
      real*8,dimension(:),allocatable :: xalpa, alpax
      real*8,dimension(:,:),allocatable ::  cgama
      real*8,dimension(:),allocatable :: xgama, gamax
      real*8,dimension(:,:),allocatable ::  comega
      real*8,dimension(:),allocatable :: xomega
      end module var_ccgama

      module var_cccoils
      integer*4 kccoils,kcloops
      real*8,dimension(:,:),allocatable :: ccoils
      real*8,dimension(:),allocatable :: xcoils
      real*8,dimension(:,:),allocatable :: cloops
      real*8,dimension(:),allocatable :: xloops
      end module var_cccoils

      module var_tionfit
      real*8 chisqti,tibdry,stibdry
      real*8,dimension(:),allocatable :: tifit
      real*8,dimension(:),allocatable :: stitho,xsiion,tithom
      end module var_tionfit

      module var_telnfit
      real*8 chisqte,tebdry,stebdry
      real*8,dimension(:),allocatable :: tefit
      end module var_telnfit

      module var_dionfit
      real*8 dibdry,sdibdry
      real*8,dimension(:),allocatable :: snitho,dnitho
      end module var_dionfit

      module var_delnfit
      real*8,dimension(:),allocatable :: defit
      real*8 chisqne,debdry,sdebdry,fco2ne
      end module var_delnfit

      module var_tsrz
      use set_kinds, only: dp
      real*8,dimension(:),allocatable :: zuperts,rlibim
      real*8 zlowerts
      real*8, parameter :: rmajts=1.94_dp
      end module var_tsrz

      module var_climxx
      real*8,dimension(:),allocatable :: xlimbd,ylimbd
      end module var_climxx

      module var_fdbkgr
      real*8 gsum
      real*8,dimension(:,:),allocatable :: grdfdb
      end module var_fdbkgr

      module var_fdbkcl
      real*8 brfb(2)
      integer*4 kct1,kct2,kct3,kct4
      real*8,dimension(:),allocatable :: fb_plasma
      end module var_fdbkcl

      module var_coiln1
      !use eparm,only:ntime,magpri
      integer*4 n1coil,iern1
      real*8,dimension(:),allocatable :: signn1
      real*8,dimension(:),allocatable :: curtn1
      end module var_coiln1

      module var_coilcc
      !use eparm,only:ntime,mccoil,micoil
      integer*4 nccoil,iercc,nicoil
      logical oldccomp,oldcomp 
      real*8,dimension(:),allocatable :: curc139,curc79,curc199,curiu30, &
        curiu90,curiu150,curil30,curil90,curil150
      real*8,dimension(:,:),allocatable :: curccoi
      real*8,dimension(:,:),allocatable :: curicoi
      end module var_coilcc

      module var_btcomp
      !use eparm,only:ntime,magpri
      integer*4 ibtcomp,ierbtc
      real*8,dimension(:),allocatable :: signbt
      end module var_btcomp

      module var_subic
      !use eparm,only:modef,modep,modew
      integer*4 nodef,nodep,nodew
      real*8,dimension(:),allocatable :: xnodef
      real*8,dimension(:),allocatable :: xnodep
      real*8,dimension(:),allocatable :: xnodew
      integer*4,dimension(:),allocatable :: kbasef
      integer*4,dimension(:),allocatable :: kbasep
      integer*4,dimension(:),allocatable :: kbasew
      end module var_subic

      module var_vtor
      integer*4 nomegat,kvtor,kwwcur,kwcurn,kplotp,npresw, &
                nsplot,kdovt
      real*8 betapw0,enw,emw,rvtor,wcurbd,gammaw,rbetaw, & 
             preswb,chiprw
      real*8,dimension(:),allocatable :: omegat,rpresw,zpresw,presw, &
                                         sigprw,rpresws,fwtprw,scalepw, &
                                         sigome,rpreswv,romegat,zomegat
      real*8,dimension(:,:),allocatable :: rprwpc
      real*8,dimension(:),allocatable :: betapw,betatw,wplasw
      real*8,dimension(:),allocatable :: rgrvt,pwprim,pressw,prwcal, &
                                         saiprw,rgsvt,saiprw2,premew
      real*8,dimension(:),allocatable :: presst
      end module var_vtor

      module var_fitec
      real*8,dimension(:),allocatable :: fwtec,cecurr,chiecc
      end module var_fitec

      module var_pflocal
      real*8 psiecn,dpsiecn,rkec,cjeccd
      real*8,dimension(:),allocatable :: cjorec,ffprec
      end module var_pflocal

      module var_ctanhts
      character*2 fitzts
      real*8,dimension(:),allocatable :: ztssym,ztswid,ptssym
      logical,dimension(:),allocatable :: ztserr
      end module var_ctanhts

      module var_qsurfac
      real*8,dimension(:),allocatable :: psin32,psin21,rq32in, &
                                         rq21top
      end module var_qsurfac
!----------------------------------------------------------------------
!--   Magnetic measurement uncertainties
!----------------------------------------------------------------------
      module var_initerror
      real*8 sigpasma,sigref
      real*8,dimension(:),allocatable :: sigfcc
      real*8,dimension(:),allocatable :: sigecc
      real*8,dimension(:),allocatable :: sigsil
      real*8,dimension(:),allocatable :: sigmpi
      end module var_initerror

      module var_magerror
      integer*4 imagsigma
      real*8 errmag,errmagb
      end module var_magerror
!
      module commonblocks
      real*8,allocatable :: c(:,:,:,:),wk(:),copy(:,:),bkx(:),bky(:), &
        cw(:,:,:,:),wkw(:),copyw(:,:),bwx(:),bwy(:), &
        sifprw(:),bwprw(:),cwprw(:),dwprw(:), &
        sfprw(:),sprwp(:),wgridpc(:),rfcpc(:,:),ct(:,:,:,:), &
        wkt(:),bkrt(:),bkzt(:),psiold(:),psipold(:), &
        work(:),sifpre(:),bwpre(:),cwpre(:),dwpre(:),sfpre(:), &
        sprep(:),worka(:),zeros(:),byringr(:),byringz(:), &
        cjrf(:),cj(:,:,:,:),wkj(:),copyj(:,:),bjx(:),bjy(:), &
        cv(:,:,:,:),wkv(:),copyv(:,:),bvx(:),bvy(:),xsisii(:), &
        f_a(:),f_b(:),f_c(:),f_d(:),pp_a(:),pp_b(:),pp_c(:), &
        pp_d(:),chi_c(:,:,:,:),chi_wk(:),chi_copy(:,:),chi_bkx(:), &
        chi_bky(:),wxin(:),wyin(:),wxout(:),wyout(:),xouts(:), &
        youts(:),bpoo(:),bpooz(:),bpooc(:),bfpol(:),cfpol(:), &
        dfpol(:),rsplt(:),zsplt(:),csplt(:),xxtra(:,:),yxtra(:,:), &
        bpxtra(:,:),flxtra(:,:),fpxtra(:,:),worka2(:)
      end module commonblocks

      module efit_bdata
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      integer*4 iunit,m_write
      real*8 xlims(5),ylims(5),xlmins
      end module efit_bdata
