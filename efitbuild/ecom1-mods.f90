      module var_acoilrz
      use eparmdud129,only:nacoil
      real*8,dimension(nacoil) :: racoil,zacoil,wacoil, &
                     hacoil
      end module var_acoilrz

      module var_ccurrn
      use eparmdud129,only:nwnh
      real*8,dimension(:),allocatable :: pcurrt,pcurrw,pcurrtpp
      integer*4 icurrt,icinit,islpfc,icprof
      data islpfc/0/

      real*8 fconst,betap0
      end module var_ccurrn

      module var_czero
      use eparmdud129,only:nwnh
      real*8,dimension(:),allocatable :: zero,www
      integer*4 iweigh
      end module var_czero

      module var_parame
      use eparmdud129,only:nw,ndata
      real*8,dimension(:),allocatable :: volp,pprime,pres,ffprim,fpol &
         ,qpsi,r2surf,rpres,curmid,cjor,bpolss &
         ,rmajz0,bprmaj,btrmaj,r1surf,r2surg &
         ,bbfpol,ccfpol,ddfpol
      real*8 r1bdry,r2bdry,carea,jwantm,r2wdry,r4wdry,rpwdry,rp2wdry 
      real*8,dimension(ndata) :: r1sdry,r2sdry,vzeroj,sizeroj 
      data vzeroj(1)/0.0/,sizeroj(1)/-1.0/,jwantm/3/
      end module var_parame

      module var_cqfit
      use eparmdud129,only:nppcur,nffcur,nwwcur
      real*8 fwtqa,qvfit,fbrdy,fcentr,cjmaxi &
           ,goli,gocur,goqq,gobeta,cbrsp &
           ,rqajtor,rqaftor,rqapetor,rqafetor
      data cbrsp/3.0e-5/
      integer*4 nqiter,ivacum 
      data nqiter/10/
      real*8,dimension(nppcur) :: rqajx,rjjjx
      real*8,dimension(nffcur) ::rqafx,rjjfx
      real*8,dimension(nwwcur) :: rqawx,rjjwx
      end module var_cqfit

      module var_rmatri
      use eparmdud129,only:nrsmat,ntime
      real*8,dimension(nrsmat) :: brsp,brspss
      real*8,dimension(ntime) :: tsaisq
      real*8 cond,tsaifc
      integer*4 nparam,nfnpcr,nfnwcr,nbase
      end module var_rmatri
  
      module var_chimag
      use eparmdud129,only:ntime
      real*8,dimension(ntime) :: chi2rm
      end module var_chimag

      module var_exdat2
      use eparmdud129,only:ntime
      real*8,dimension(ntime) :: bcentr
      real*8 rcentr
      data rcentr/1.6955/
      end module var_exdat2

      module var_exdata
      use eparmdud129,only:nsilop,magpri,nacoil,nfcoil,npcurn
!      public
      integer*4 ishot,itime,ifitvs,iacoil,itimeu,kersil
      data errsil/0.03/
      data kersil/2/
!      common/exdata/ishot,itime,serror,fwtsi(nsilop),fwtmp2(magpri), &
!           fwacoil(nacoil),fwtcur,elomin,fwtbp,fwtdlc,fwtfc(nfcoil), &
!           ifitvs,rwtsi(nsilop),rwtmp2(magpri),iacoil, &
!           cbrspv(npcurn),itimeu,kersil,errsil

      real*8 serror,fwtcur,elomin,fwtbp,fwtdlc,errsil 
      real*8,dimension(nsilop) :: fwtsi,rwtsi
      real*8,dimension(magpri) :: fwtmp2,rwtmp2
      real*8,dimension(nacoil) :: fwacoil
      real*8,dimension(nfcoil) :: fwtfc
      real*8,dimension(npcurn) :: cbrspv
      data elomin/0.90/
      data iacoil/0/,fwacoil/1*0./
      end module var_exdata

      module var_texdat
      use eparmdud129,only:ntime,nsilop,magpri,nacoil & 
      ,nfcoil,nesum,nco2v,nco2r,micoil,nrogow
      integer*4 ierpla,ico2,ierrdi 
      real*8 bitip
      real*8,dimension(ntime,nsilop) :: silopt
      real*8,dimension(ntime,magpri) :: expmpi
      real*8,dimension(ntime,nacoil) :: accurt
      real*8,dimension(ntime,nfcoil) ::fccurt
      real*8,dimension(ntime,nesum) :: eccurt
      real*8,dimension(ntime,nco2v) :: denvt
      real*8,dimension(ntime,nco2r) :: denrt
      integer*4,dimension(ntime,micoil) ::iccurt
      real*8,dimension(ntime) :: pasmat,time,pbinj 
      integer*4,dimension(nsilop) :: ierpsi
      integer*4,dimension(nrogow) :: ierpr
      integer*4,dimension(magpri) :: iermpi
      real*8,dimension(nsilop) :: psibit
      real*8,dimension(nrogow) :: prbit
      real*8,dimension(magpri) :: bitmpi
      real*8 :: vbit
           
      real*8,dimension(ntime) :: vloopt,psiref,diamag,sigdia,psirefs 
      real*8,dimension(nfcoil) :: bitfc
      integer*4,dimension(nfcoil) :: ierfc
      real*8,dimension(nesum) :: bitec
      integer*4,dimension(nesum) :: ierec
      real*8,dimension(micoil) ::bitic
      integer*4,dimension(micoil) ::ieric
      integer*4 ierdia(3)
!sri-mpi
      integer*4 nopbinj
      data vbit/10./
      end module var_texdat

      module var_savfit
      use eparmdud129,only:nsilop,magpri,nfcoil,nesum
      real*8,dimension(nsilop) :: swtsi
      real*8,dimension(magpri) :: swtmp2
      real*8,dimension(nfcoil) :: swtfc
      real*8,dimension(nesum) :: swtec
      real*8 swtcur,swtdlc
      end module var_savfit

      module var_cstark
      use eparmdud129,only:ntime,nstark,npcurn,nfcoil,nesum &
                 ,nvesel,nacoil,ngam_vars,ngam_u,ngam_w
      real*8,dimension(ntime,nstark):: tangam, siggam,a1gam &
           ,a2gam,a3gam,a4gam,tangam_uncor
      real*8,dimension(nstark) :: fwtgam,chigam,swtgam
      real*8 v30lt,v30rt,v210rt,v210lt
      integer*4,dimension(nstark) :: iergam
      integer*4,dimension(nstark) :: mseport,mse_spave_on
      integer*4  kstark,iplots,kmtark,klibim,kdomse &
                 ,msebkp,msefitfun,kwaitmse,mse_quiet &
                 ,mse_strict,ok_30rt,ok_210lt,mse_usecer &
		 ,mse_certree,mse_use_cer330,mse_use_cer210
      data iplots/1/
      data  msebkp/0/,msefitfun/1/
      data kwaitmse/0/,dtmsefull/0.0/
      data mse_strict/0/,t_max_beam_off/0.0/,ok_30rt/0/,ok_210lt/0/
      data kdomse/0/
      data mse_usecer/0/,mse_certree/0/
      data mse_use_cer210/0/,mse_use_cer330/0/
      data v30lt/0.0/,v30rt/0.0/
      data v210lt/0.0/,v210rt/0.0/

      real*8 chigamt,chilibt,dtmsefull,t_max_beam_off 
      real*8,dimension(nstark,npcurn) :: rbrpc,rbzpc,rgampc
      real*8,dimension(nstark,nfcoil) :: rbrfc,rbzfc
      real*8,dimension(nstark,nesum) :: rbrec,rbzec,rgamec
      real*8,dimension(nstark,nvesel) :: rbrvs,rbzvs,rgamvs
      real*8,dimension(nstark,nfcoil) :: rgamfc
      real*8,dimension(ntime,nstark) :: rhsgam,rrgam,zzgam &
           ,starkar,starkaz,a5gam,a6gam,a7gam,a8gam
      real*8,dimension(nstark,ntime) :: cmgam,spatial_fix
      real*8,dimension(nstark,nacoil) :: rgamac,rbrac,rbzac
      real*8,dimension(nstark) :: btgam,sistark,qstark &
           ,rmse_gain,rmse_slope,rmse_scale,rmse_offset,sigam &
           ,bzmse,bzmsec,cjmse,cjmsec,rhogam
      real*8,dimension(nstark,ngam_vars,ngam_u,ngam_w) :: spatial_avg_gam
      end module var_cstark

      module var_msels
      use eparmdud129,only:ntime,nmsels,npcurn,nfcoil,nesum &
                 ,nvesel,nacoil
      real*8,dimension(ntime,nmsels):: bmselt,sbmselt,fwtbmselt, &
           rrmselt,zzmselt,l1mselt,l2mselt,l4mselt,emselt,semselt, &
           fwtemselt,swtbmselt,swtemselt,cmmls,cmels,rhsmls,rhsels, &
           l3mselt,cmmls2,cmmlsv
      real*8,dimension(nmsels)::swtbmsels,swtemsels,chimls,chiels
      real*8 :: avemsels, fmlscut, tchimls, tchiels
      real*8,dimension(nmsels,npcurn) :: rmlspc
      real*8,dimension(nmsels,nesum) :: rmlsec,relsec
      real*8,dimension(nmsels,nvesel) :: rmlsvs,relsvs
      real*8,dimension(nmsels) :: simls,sinmls,btmls,brmls,bzmls
      real*8,dimension(ntime) :: chi2mls, chi2gamt

      integer*4,dimension(ntime,nmsels) :: iermselt
      integer*4,dimension(ntime) :: kerrot
      integer*4 kdomsels, mmbmsels
      character*3 :: synmsels
      data avemsels/10./, synmsels/'SYN'/, kdomsels/0/, fmlscut/1.e-6/
      end module var_msels

      module var_rpedge
      use eparmdud129,only:nsilop,magpri,mpress,mbdry,nstark,nmsels
      real*8,dimension(nsilop) :: rsilpe
      real*8,dimension(magpri) :: rmp2pe
      real*8 fgowpe
      real*8,dimension(mpress) :: rprepe
      real*8,dimension(mbdry)   :: gbdrpe
      real*8,dimension(nstark) :: rgampe,rbrpe,rbzpe
      real*8,dimension(nmsels) :: rmlspe,relspe
      end module var_rpedge

      module var_rfedge
      use eparmdud129,only:nsilop,magpri,mbdry,nstark,nmsels
      real*8,dimension(nsilop) :: rsilfe
      real*8,dimension(magpri) :: rmp2fe
      real*8 fgowfe
      real*8,dimension(mbdry)   :: gbdrfe
      real*8,dimension(nstark) :: rgamfe,rbrfe,rbzfe
      real*8,dimension(nmsels) :: rmlsfe,relsfe
      end module var_rfedge

      module var_cece
      use eparmdud129,only:nnecein,nwnh,nesum,nnece,nfcoil,npcurn &
                            ,ntime,nvesel,nacoil,kxiter,nnnte
      integer*4 ierecebz,necein,jo,kfixro,nece,kece,kecebz,mecein &
                ,kfitece,kinputece,kcallece,nharm,kdoece,kgeteceb & 
                ,mtxece,ksetece,kfixrece,nfit,kcmin,nconstr

      data ksetece/0/,kgeteceb/0/,kcallece/2/
      data kdoece/0/,nconstr/1/,mtxece/0/

      character*4 eceiter
      data eceiter/'pair'/

      real*8 receo,fwtecebz,fwtecebz0,rteo,zteo,robit,swtecebz,chiecebz &
             ,ecebzfit,ecebzbit,fwtnow,zeceo,chisqfit,xfit(10),tchiece &
             ,recebzdz,gecebzdz,eceerror

      data eceerror/0.03/

      real*8,dimension(nnecein) :: teecein,feece,errorece &
           ,teecein0,feece0,errorece0 &
           ,becein,recein,teeceinr
      real*8,dimension(nnece) :: recem,recep,fwtece,fwtece0 &
                                 ,swtece,chiece,ecefit,ecebit
      real*8,dimension(nfcoil) :: recebzfc
      real*8,dimension(:),allocatable :: gecebzpc
      real*8,dimension(nesum) :: recebzec
      real*8,dimension(nnece,nfcoil) :: recefc
      real*8,dimension(:,:),allocatable :: gecepc,geceppc,gecempc
      real*8,dimension(nnece,nesum) :: receec
      integer*4,dimension(nnece) :: ierece
      real*8,dimension(nnece,npcurn) :: recepc
      real*8,dimension(ntime,nnece) :: brspece
      real*8,dimension(npcurn) :: recebzpc
      real*8,dimension(ntime) :: brspecebz,cmecebz
      real*8,dimension(nnece,ntime) :: cmece
      real*8,dimension(nnece,nvesel) :: recevs
      real*8,dimension(nvesel) :: recebzvs
      real*8,dimension(nnece) :: recedz,gecedz,rtep,rtem,rpbit,rmbit
      real*8,dimension(nnece,nacoil) :: receac
      real*8,dimension(nacoil) :: recebzac
      real*8,dimension(nnnte) :: teecer,rrr,bbf,teeceb
      real*8,dimension(kxiter) :: receoi
      real*8,dimension(kxiter,nnece) :: recemi,recepi
      end module var_cece

      module var_calchi
      use eparmdud129,only:nfcoil,nsilop,magpri,mpress,nrogow
      real*8 chipre,chitot,saiip
      real*8,dimension(nsilop) :: saisil
      real*8,dimension(magpri) :: saimpi
      real*8,dimension(nrogow) :: saipr
      real*8,dimension(mpress) :: saipre
      real*8,dimension(nfcoil) :: saifc
      end module var_calchi

      module var_dlc
      use eparmdud129,only:ntime,nffcur
      real*8 sigdlc
      real*8,dimension(ntime) :: dfluxc,cdflux,edflux
      real*8,dimension(nffcur) :: rspdlc
      end module var_dlc
    
      module var_comco2
      use eparmdud129,only:ntime,nco2r,nco2v
      real*8,dimension(nco2r,ntime) :: rco2r
      real*8,dimension(nco2v,ntime) :: rco2v
      real*8,dimension(nco2v) :: chordv
      real*8,dimension(nco2r) :: chordr
      real*8 zcentr
      real*8,dimension(ntime,nco2r) :: dco2r
      real*8,dimension(ntime,nco2v) :: dco2v
      data chordv/1.486,1.945,2.098/
      data chordr/.00, 0.1524/,zcentr/0./
      end module var_comco2

      module var_check
      use eparmdud129,only:ntime
      integer*4,dimension(ntime,30) :: erflag
      integer*4 kflag(30),lflag,ktimeo
      end module var_check

      module var_consum
      use eparmdud129,only:kxiter
      real*8 condno,condin 
      data condin/1.0e-06/
      real*8,dimension(kxiter) :: cerror,csibry,csimag, &
                  cvolp,crmaxi,czmaxi,cemaxi,cqmaxi,cchisq &
                  ,brfbc,tvfbrt,cdelz
      end module var_consum

      module var_cxray
      use eparmdud129,only:nangle,ntangle
      real*8,dimension(nangle+ntangle) :: rxray,zxray,xangle
      integer*4 ksxr0(10),ksxr2(10),idosxr
      data ksxr0/10*0/,ksxr2/10*0/,idosxr/1/
      data xangle/    120.,124.,128.,132.,136.,140.,144.,148., &
                     152.,156.,160.,164.,168.,172.,176.,180.,180., &
                      184.,188.,192.,196.,200.,204.,208.,212.,216., &
                      220.,224.,228.,232.,236.,240., &
                   294.5,290.5,286.5,282.5,278.5,274.5,270.5,266.5, &
                   262.5,258.5,254.5,250.5,246.5,242.5,238.5,234.5, &
                   234.5,230.5,226.5,222.5,218.5,214.5,210.5,206.5, &
                   202.5,198.5,194.5,190.5,186.5,182.5,178.5,174.5, &
                   0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./
      data zxray/-10.7,-10.7,-10.7,-10.7,-10.7,-10.7, &
                    -10.7,-10.7,-10.7,-10.7,-10.7,-10.7,-10.7,-10.7, &
                    -10.7,-10.7,-14.7,-14.7,-14.7,-14.7,-14.7,-14.7, &
                    -14.7,-14.7,-14.7,-14.7,-14.7,-14.7,-14.7,-14.7, &
                    -14.7,-14.7, &
                    130.1,130.1,130.1,130.1,130.1,130.1,130.1,130.1, &
                    130.1,130.1,130.1,130.1,130.1,130.1,130.1,130.1, &
                    132.6,132.6,132.6,132.6,132.6,132.6,132.6,132.6, &
                    132.6,132.6,132.6,132.6,132.6,132.6,132.6,132.6, &
                   0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./
      data rxray/248.9,248.9,248.9,248.9,248.9,248.9,248.9, &
                248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9, &
                248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9, &
                248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9, &
                197.6,197.6,197.6,197.6,197.6,197.6,197.6,197.6, &
                197.6,197.6,197.6,197.6,197.6,197.6,197.6,197.6, &
                194.5,194.5,194.5,194.5,194.5,194.5,194.5,194.5, &
                194.5,194.5,194.5,194.5,194.5,194.5,194.5,194.5, &
                   0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./
!---D3D-----------------------D3D----------------------------D3D-----
      end module var_cxray

      module var_mercie
      use eparmdud129,only:nw
      integer*4 imerci
      real*8,dimension(:),allocatable :: rzzmax,zzmax
      end module var_mercie

      ! OPT_INPUT >>>
      module opt_input
         integer*4 mode_in
         character cmdfile_in*15, shotfile_in*15
         integer*4 shot_in
         real*8 starttime_in
         real*8 deltatime_in
         integer*4 steps_in
         logical use_opt_input
         character snapext_in*82
         character(80),dimension(:),allocatable :: inpfile_in
      end module opt_input
      ! OPT_INPUT <<<

      ! MPI >>>
      module mpi_info
         integer :: rank, nproc, ierr
         integer*4,dimension(:),allocatable :: dist_data
         integer*4,dimension(:),allocatable :: dist_data_displs
         double precision,dimension(:,:),allocatable :: fwtgam_mpi
      end module mpi_info
      ! MPI <<<

!jal 2/23/04 add iplcout=1 print plasma and pf currents to gfile
      module var_input1
      use eparmdud129,only:ndata,mpress,nvesel,nfcoil
      logical write_Kfile ,fitfcsum

      integer*4 icondn,itek,kdata,itrace,ierchk,iconvr,ixray,itell &
           ,kprfit,licalc,ibound,ibatch,idite,ilaser, lookfw &
            ,kdot, icutfp, keqdsk,kdofit,kbetapr,kplotpr,kpressb &
            ,kdoqn,kfcurb,kpcurb, kzeroj, ncstne,ncstte,ncstfp,ncstpp &
            ,kgrid,kwripre,negcur, kframe,kskipvs,icntour, iavdpl &
            ,jwake,limvs, kbound,kgraph,istore,iout,kdopre &
            ,iishot,kktime,iplcout,ksigma, kwritime 
      integer*4 iteks, mxiters, zelipss, n1coils
      integer*4 itekt, mxitert, zeliptt, n1coilt
      integer*4 kerrorq

      data kgrid/1/,kwripre/0/,kwritime/0/
      data licalc/1/
      data ksigma/0/
      data kdoqn/0/,kerrorq/0/
      data icntour/0/,cstabne/0.0/,cstabte/0.0/ &
           ,limvs/1/,ncstpp/1/,ncstfp/1/,kzeroj/0/
      data kbound/0/
      data kdofit/0/
      data kdopre/0/

      real*8 cutip, dtdot, xpsimin,fcurbd,pcurbd,prbdry,sgprmin &
           ,prespb,tipbry,tepbry,dipbry,depbry,pbimpb,sigppb,sigpreb &
           ,sigtipb,sigtepb,sigdipb,sigdepb,fwtpd,fwtfd,cstabte,cstabne &
           ,dp1dxf, sgtimin, alphafp,xpsialp, vsdamp &
           ,rminvs,rmaxvs,zminvs,zmaxvs,relbps,zbound,yvs2,saimin &
           ,fztor,fzpol,tcurrp,fpolvs, rbound, dnmin, sigprebi,pressbi &
           ,alphamu,saicon,rsepex, ttimeb,ddtime
      data alphafp/0./,sigppb/1000./
      data kframe/0/,rminvs/0/,rmaxvs/100./,zminvs/-100./,zmaxvs/100./
      data kskipvs/0/,vsdamp/0/,relbps/0.004/,zbound/0.0/,rbound/0.0/
      data dnmin/1.0/
      data saimin/60.0000/,saicon/60.0/

      real*8,dimension(ndata) :: rzeroj
      real*8,dimension(mpress) :: fwtpre
      real*8,dimension(nvesel) :: vforcep,vforcet
      real*8,dimension(nfcoil) :: fwtfcsum
      character*82  snap_file 
      end module var_input1

      module var_inputc
      use eparmdud129,only:ntime,nvesel,magpri,nsilop

      character*12 mfitpop
      character*5 mfvers(2)
      data mfvers(1)/'10/13'/,mfvers(2)/'/2020'/
      character(4),dimension(ntime) :: limloc
      character(10),dimension(nvesel) :: vsname
      character(10),dimension(magpri) :: mpnam2
      character(10),dimension(nsilop) :: lpname
      character  filimt*100,cshot*6,jdebug*4
      integer idebug,efitversion
      data idebug/0/,efitversion/20201013/
      data jdebug/'NONE'/
      end module var_inputc

      module var_input4
      use eparmdud129,only:ntime
      character(80),dimension(ntime) :: ifname
      end module var_input4

      module var_switch
      use eparmdud129,only:npcurn,mqwant
      integer*4 nqaxis,isumip,jbeta,jli,nqwant
      data nqwant/0/,jbeta/1/,jli/2/
      real*8 sumip,fbetap,fli,fbetat,fbetan 
      real*8,dimension(mqwant) :: qsiw, pasmsw,fqsiw,siwantq
      real*8,dimension(npcurn,mqwant) :: fgowsw
      end module var_switch

      module var_siloop
      use eparmdud129,only:nsilop
      real*8,dimension(nsilop) :: rsi,zsi,wsi,hsi
      end module var_siloop

      module var_slname
      use eparmdud129,only:nsilop
      real*8,dimension(nsilop) :: as,as2
      integer*4 nslref
      data nslref/1/
      end module var_slname

      module var_ecoil 
      use eparmdud129,only:necoil
      real*8,dimension(necoil) :: re,ze,we,he,ecid
      end module var_ecoil 

      module var_fcoil 
      use eparmdud129,only:mfcoil
      real*8,dimension(mfcoil) :: rf,zf,wf,hf, &
           af,af2,rsisfc,turnfc,fcid,fcturn
      data fcid/1. ,2. ,3. ,4. ,5. ,6. ,7. ,8. ,9. , &
                10.,11.,12.,13.,14.,15.,16.,17.,18./
      end module var_fcoil 

      module var_mprobe
      use eparmdud129,only:magpri
      real*8,dimension(magpri) :: xmp2,ymp2,amp2,smp2,patmp2
      integer*4 nsmp2
      end module var_mprobe

      module var_limite
      use eparmdud129,only:nlimit
      integer*4 limitr,iplim,limfag,limitr_180
      data limfag/2/
      real*8,dimension(nlimit) :: xlim,ylim, xlim_180, ylim_180
      real*8,dimension(2) :: rwstrip1, zwstrip1,rwstrip2, zwstrip2 
      real*8 yltype, yltype_180
      logical dowstrip
      data rwstrip1(1)/1.33/,zwstrip1(1)/-1.363/,rwstrip1(2)/1.38/,zwstrip1(2)/-1.363/
      data rwstrip2(1)/1.4075/,zwstrip2(1)/-1.250/,rwstrip2(2)/1.4575/,zwstrip2(2)/-1.250/
      data dowstrip/.F./  
      end module var_limite

      module var_mimite
      use eparmdud129,only:nlimit
      integer*4 mimitr,mimitr_180
      real*8,dimension(nlimit) ::  xmim,ymim,xmim_180,ymim_180
      end module var_mimite

      module var_udata
      use eparmdud129,only:ntime
      integer*4,dimension(ntime) :: ipsi,irogw,imag2,iplasm &
           ,idlopc,ifc,iec
      end module var_udata

      module var_morsum
      use eparmdud129,only:kxiter
      real*8,dimension(kxiter) :: csumip,tratio,aveerr
      integer*4,dimension(kxiter) :: iermax,jermax
      end module var_morsum

      module var_bdsend
      use eparmdud129,only:npoint
      integer*4 nbbbs,nbskip,nbdrymx, nbdryp
      data nbskip/2/,nbdrymx/110/
      real*8,dimension(:),allocatable :: rbbbs,zbbbs
      end module var_bdsend

      module var_fxbry
      use eparmdud129,only:mbdry,nfcoil,nacoil,npcurn,nwnh
      integer*4 nbdry,nbdryss,nsol
      data nsol/0/
      logical fitts 
      real*8 wsisol,dselsum/1.e-12/
      real*8,dimension(mbdry) :: rbdry,zbdry,fwtbdry,fwtbry,sigrbd &
                                 ,sigzbd,rbdry0, zbdry0,rbdryss,zbdryss
      real*8,dimension(mbdry) :: rsol, zsol, fwtsol, fwtsolw
      real*8,dimension(mbdry,nfcoil) :: rbdrfc, rsolfc
      real*8,dimension(mbdry,nacoil) :: rbdrac
      real*8,dimension(:,:),allocatable :: rbdrpc,rsolpc
      real*8,dimension(mbdry,npcurn) :: gbdrpc
      end module var_fxbry

      module var_fwtdz
      use eparmdud129,only:magpri,nsilop,mbdry,nstark,mpress,nwnh,nmsels
      logical fitdelz
      integer*4 ndelzon,ifitdelz 
      data fitdelz/.false./,ndelzon/999/,relaxdz/1.000/, &
           stabdz/-1.e-4/,scaledz/1.e-03/,ifitdelz/1/
      data errdelz/0.06/
      real*8 errdelz,fgowdz,scaledz,stabdz,relaxdz,cdeljsum
      real*8,dimension(magpri) :: gmp2dz
      real*8,dimension(nsilop) :: gsildz
      real*8,dimension(nstark) :: gbrdz,gbzdz,rgamdz
      real*8,dimension(nmsels) :: rmlsdz,relsdz
      real*8,dimension(mpress) :: rpredz,rprwdz
      real*8,dimension(mbdry) :: gbdrdz
      real*8,dimension(:),allocatable :: rdjdz
      end module var_fwtdz
      
      module var_combry
      use eparmdud129,only:msbdry
      real*8,dimension(msbdry) :: erbloc, erbsloc
      real*8 erbmax,erbave,erbcom,erbsmax,erbsave
      data erbcom/1.0e-2/
      end module var_combry

      module var_fbysta
      use eparmdud129,only:nfcoil
      real*8 cfcoil, psibry0
      integer*4 ifref
      real*8,dimension(nfcoil) :: fcsum,fczero
      end module var_fbysta

      module var_prdata
      use eparmdud129,only:mpress,ndata
      integer*4 npress,npteth,nption,npneth,nbeam,ndokin,nbrmcrd &
               ,nptef,npnef,nptionf,nmass 
      data ndokin/1/
      real*8 pbeamb,pressb,zeffvs,zlowimp
      integer*4,dimension(ndata) :: ivbcuse

      real*8,dimension(mpress) :: pressr,rpress,zpress,tethom,rteth &
                   ,zteth,tionex,rion,zion,dnethom,rneth,zneth &
                   ,pbeam,sibeam,sgteth,sigti,sigpre,sgneth,precal &
                   ,dnbeam,dmass,scalepr,premea,saipre2
      real*8,dimension(ndata) :: pbimth,bremin,bremsig,brmrtan, &
                                 brmzelev,dnbthom
      end module var_prdata

      module var_cerfit
      use eparmdud129,only:nstark,nercur,nw,nmsels
      integer*4 keecur,keefnc,keeknt,needer,keehord
      integer*4,dimension(nercur) :: keebdry,kee2bdry
      real*8 ecurbd,eetens
      real*8,dimension(nstark,nercur) :: rgamer
      real*8,dimension(nmsels,nercur) :: rmlser,relser
      real*8,dimension(nercur) :: eeknt,eebdry,ee2bdry,cerer
      real*8,dimension(nstark) :: ermse
      real*8,dimension(nstark,nercur) :: e1rbz,e2rbz,e3rbr
      real*8,dimension(:),allocatable :: ermid, eshear,epoten,rhovn, &
                    rpmid,xmid,sigrid,sipmid, &
                    brhovn,crhovn,drhovn,rhopmid
      data keecur/0/,ecurbd/0.0/,keefnc/0/,eetens/5.0/
      end module var_cerfit

      module var_ccgama
      use eparmdud129,only:nppcur,nffcur,nwcurn
      integer*4 kcalpa,kcgama,kcomega
      real*8 fwtxx,fwtxxj,fwtxxq,fwtxxb,fwtxli 
      data fwtxx/0.2/,fwtxxj/1./
      data  fwtxxq/1./,fwtxxb/1./,fwtxli/1./
      data kcalpa/0/,kcgama/0/,kcomega/0/

      real*8,dimension(nppcur,nppcur) :: calpa
      real*8,dimension(nppcur) :: xalpa, alpax
      real*8,dimension(nffcur,nffcur) :: cgama
      real*8,dimension(nffcur) :: xgama, gamax
      real*8,dimension(nwcurn,nwcurn) :: comega
      real*8,dimension(nwcurn) :: xomega
      end module var_ccgama

      module var_cccoils
      use eparmdud129,only:nfcoil,nsilop
      integer*4 kccoils,kcloops
      data kccoils/0/
      real*8,dimension(nfcoil,nfcoil) :: ccoils
      real*8,dimension(nfcoil) :: xcoils
      real*8,dimension(nsilop,nsilop) :: cloops
      real*8,dimension(nsilop) :: xloops
      end module var_cccoils

      module var_tionfit
      use eparmdud129,only:nppcur,ndata
      real*8 chisqti, tibdry,stibdry
      real*8,dimension(nppcur) :: tifit
      real*8,dimension(ndata) :: stitho,xsiion,tithom
      end module var_tionfit

      module var_telnfit
      use eparmdud129,only:nppcur
      real*8 chisqte,tebdry,stebdry
      real*8,dimension(nppcur) :: tefit
      end module var_telnfit

      module var_dionfit
      use eparmdud129,only:ndata
      real*8 dibdry,sdibdry
      real*8,dimension(ndata) :: snitho,dnitho
      end module var_dionfit

      module var_delnfit
      use eparmdud129,only:nppcur
      real*8,dimension(nppcur) :: defit
      real*8 chisqne,debdry,sdebdry,fco2ne
      data fco2ne/1.0/
      end module var_delnfit

      module var_tsrz
      use eparmdud129,only:ntime
      real*8,dimension(ntime) :: zuperts
      real*8 zlowerts,rmajts,rlibim(ntime),zlibim
      data rmajts/1.94/,zlibim/-0.127/
      end module var_tsrz

      module var_climxx
      use eparmdud129,only:nlimbd
      real*8,dimension(nlimbd) :: xlimbd,ylimbd
      end module var_climxx

      module var_fdbkgr
      use eparmdud129,only:nwnh,nfbcoil
      real*8 gsum
      real*8,dimension(:,:),allocatable :: grdfdb
      end module var_fdbkgr

      module var_fdbkcl
      use eparmdud129,only:ntime
      real*8 brfb(2)
      integer*4 kct1,kct2,kct3,kct4
      real*8,dimension(ntime) :: fb_plasma
      end module var_fdbkcl

      module var_coiln1
      use eparmdud129,only:ntime,magpri
      integer*4 n1coil,iern1
      data n1coil/0/
      real*8,dimension(magpri) :: signn1
      real*8,dimension(ntime) :: curtn1
      end module var_coiln1

!*****							 	!EJS(2014)
!  Need to update values -- wherever they are specified:
!	mccoil = 6   (probably is 3 now)
!	micoil = 12  (probably is 6 now)
!*****
      module var_coilcc
      use eparmdud129,only:ntime,mccoil,micoil
      integer*4 nccoil,iercc,nicoil
      data nccoil/1/,nicoil/1/
      logical oldccomp,oldcomp 
      data oldccomp/.true./,oldcomp/.false./
!*****								!EJS(2014)
!	Not sure why these are separate arrays.
!	They seem to be the same as the columns of curccoi and curicoi
!	Are they only used in magsigma?  If so, the following changes are
!	probably not needed until some future time when magsigma is updated.
!*****
      real*8,dimension(ntime) :: curc139,curc79,curc199,curiu30 &
                    ,curiu90,curiu150,curil30,curil90,curil150 &
                    ,curc259,curc319,curc19        &		!EJS(2014)
                    ,curiu210,curiu270,curiu330    &		!EJS(2014)
                    ,curil210,curil270,curil330			!EJS(2014)
      real*8,dimension(ntime,mccoil) :: curccoi
      real*8,dimension(ntime,micoil) :: curicoi
      end module var_coilcc

      module var_btcomp
      use eparmdud129,only:ntime,magpri
      integer*4 ibtcomp,ierbtc
      data ibtcomp/1/
      real*8,dimension(magpri) :: signbt
      real*8,dimension(ntime) :: bti322
      end module var_btcomp

      module var_subic
      use eparmdud129,only:modef,modep,modew
      integer*4 nodef,nodep,nodew
      real*8,dimension(modef) :: xnodef
      real*8,dimension(modep) :: xnodep
      real*8,dimension(modew) :: xnodew
      integer*4,dimension(modef) :: kbasef
      integer*4,dimension(modep) :: kbasep
      integer*4,dimension(modew) :: kbasew
      end module var_subic

      module var_vtor
      use eparmdud129,only:mpress,nwwcur,ntime,nw,nwnh
      integer*4 nomegat,kvtor,kwwcur,kwcurn,kplotp,npresw &
                ,nsplot,kdovt
      real*8 betapw0,enw,emw,rvtor,wcurbd,gammaw,rbetaw, & 
             preswb,chiprw
      data kvtor/0/,rvtor/1.70/,preswb/0.0/,wcurbd/0.0/,betapw0/0.0/ &
           ,kwwcur/2/,kplotp/1/,kdovt/0/,nsplot/4/ 
      real*8,dimension(mpress) :: omegat,rpresw,zpresw,presw, &
                       sigprw,rpresws,fwtprw,romegat,zomegat &
                       ,sigome,rpreswv,scalepw
      real*8,dimension(mpress,nwwcur) :: rprwpc
      real*8,dimension(ntime) :: betapw,betatw,wplasw
      real*8,dimension(:),allocatable :: rgrvt,pwprim,pressw,prwcal, &
                              saiprw,rgsvt
      real*8,dimension(:),allocatable :: presst
      end module var_vtor

      module var_fitec
      use eparmdud129,only:nesum
      real*8,dimension(nesum) :: fwtec,cecurr,saiec
      logical writepc
      data writepc/.false./
      end module var_fitec

      module var_pflocal
      use eparmdud129,only:nw
      real*8 psiecn,dpsiecn,rkec,cjeccd
      data psiecn/0.0/,dpsiecn/0.0/
      real*8,dimension(:),allocatable :: cjorec,ffprec
      end module var_pflocal

      module var_ctanhts
      use eparmdud129,only:ntime
      character*2 fitzts
      data fitzts/'no'/
      real*8,dimension(ntime) :: ztssym,ztswid,ptssym
      logical,dimension(ntime) :: ztserr
      end module var_ctanhts

      module var_qsurfac
      use eparmdud129,only:ntime
      real*8,dimension(ntime) :: psin32,psin21,rq32in, &
                     rq21top
      end module var_qsurfac
!----------------------------------------------------------------------
!-- New magnetic uncertainties commons                               --
!----------------------------------------------------------------------
      module var_initerror
      use eparmdud129,only:nfcoil,nesum,nsilop,magpri
      real*8 sigmaip0
      real*8,dimension(nfcoil) :: sigmaf0
      real*8,dimension(nesum) :: sigmae0
      real*8,dimension(nsilop) :: sigmafl0
      real*8,dimension(magpri) :: sigmamp0
      end module var_initerror

      module var_magerror
      use eparmdud129,only:ntime,nfcoil,nesum,nsilop,magpri
      integer*4 imagsigma,icountmagsigma
      data imagsigma/0/, errmag/1.0e-3/ &
           , errmagb/1.e-2/
      real*8 errmag,errmagb
      real*8,dimension(ntime,nfcoil) :: sigmaf
      real*8,dimension(ntime) :: sigmab,sigmaip
      real*8,dimension(ntime,nesum) :: sigmae
      real*8,dimension(ntime,nsilop) :: sigmafl,gradsfl
      real*8,dimension(ntime,magpri) :: sigmamp,gradsmp,bpermp
      end module var_magerror

      module var_psilopdat
      use eparmdud129,only:ntime,nsilop
      real*8,dimension(nsilop) :: psircg,psi_k,psi_rc,vrespsi,t0psi
      real*8,dimension(ntime,nsilop) :: devpsi,rnavpsi
      integer*4,dimension(ntime,nsilop) :: navpsi
      end module var_psilopdat

      module var_plasmacurrdat
      use eparmdud129,only:ntime
      real*8 prcg,p_k,p_rc,vresp,t0p
      real*8,dimension(ntime) :: devp,navp,rnavp
      end module var_plasmacurrdat

      module var_ccoilsdat
      use eparmdud129,only:ntime,mccoil
      real*8,dimension(mccoil) :: ccrcg,cc_k,cc_rc,vrescc,t0cc
      real*8,dimension(ntime,mccoil) :: devcc
      integer*4,dimension(ntime,mccoil) :: navcc
      end module var_ccoilsdat

      module var_icoilsdat
      use eparmdud129,only:ntime,micoil
      real*8,dimension(micoil) :: xicrcg,xic_k,xic_rc,vresxic,t0xic
      real*8,dimension(ntime,micoil) :: devxic
      integer*4,dimension(ntime,micoil) :: navxic
      end module var_icoilsdat

      module var_n1coildat
      use eparmdud129,only:ntime
      real*8 xn1rcg,xn1_k,xn1_rc,vresxn1,t0xn1
      real*8,dimension(ntime) :: devxn1
      integer*4,dimension(ntime) :: navxn1
      end module var_n1coildat

      module var_vloopdat
      use eparmdud129,only:ntime
      real*8 vlrcg,vl_k,vl_rc, vresvl,t0vl
      real*8,dimension(ntime) :: devvl
      integer*4,dimension(ntime) :: navvl
      end module var_vloopdat

      module var_diamdat
      use eparmdud129,only:ntime
      real*8 diamrcg,diam_k,diam_rc,vresdiam,t0diam
      real*8,dimension(ntime) :: devdiam
      integer*4,dimension(ntime) :: navdiam
      end module var_diamdat

      module var_denvdat
      use eparmdud129,only:ntime,nco2v
      real*8,dimension(nco2v) :: denvrcg,denv_k,denv_rc,vresdenv,t0denv
      real*8,dimension(ntime,nco2v) :: devdenv
      integer*4,dimension(ntime,nco2v) ::    navdenv
      end module var_denvdat

      module var_denrdat
      use eparmdud129,only:ntime,nco2r
      real*8,dimension(nco2r) :: denrrcg,denr_k,denr_rc,vresdenr,t0denr
      real*8,dimension(ntime,nco2r) :: devdenr
      integer*4,dimension(ntime,nco2r) :: navdenr
      end module var_denrdat

      module var_magprobdat
      use eparmdud129,only:ntime,magpri
      real*8,dimension(magpri) :: xmprcg,xmp_k,xmp_rc,vresxmp,t0xmp
      real*8,dimension(ntime,magpri) :: devxmp,rnavxmp
      integer*4,dimension(ntime,magpri) :: navxmp
      end module var_magprobdat

      module var_btcompdat
      use eparmdud129,only:ntime
      real*8 btrcg,bt_k,bt_rc,vresbt,t0bt
      real*8,dimension(ntime) :: devbt
      integer*4,dimension(ntime) :: navbt
      end module var_btcompdat

      module var_btordat
      use eparmdud129,only:ntime
      real*8 bcrcg,bc_k,bc_rc,vresbc,t0bc
      real*8,dimension(ntime) :: devbc, rnavbc
      integer*4,dimension(ntime) :: navbc
      end module var_btordat

      module var_fcoildat
      use eparmdud129,only:ntime,nfcoil
      real*8,dimension(nfcoil) :: fcrcg,fc_k,fc_rc,vresfc,t0fc
      real*8,dimension(ntime,nfcoil) :: devfc,rnavfc
      integer*4,dimension(ntime,nfcoil) :: navfc
      end module var_fcoildat

      module var_ecoildat
      use eparmdud129,only:ntime,nesum
      real*8,dimension(nesum) :: ercg,e_k,e_rc,vrese,t0e
      real*8,dimension(ntime,nesum) :: deve,rnavec
      integer*4,dimension(ntime,nesum) :: navec
      end module var_ecoildat

      module var_beamdat
      use eparmdud129,only:ntime
      real*8 beamrcg,beam_k,beam_rc,vresbeam,t0beam
      real*8,dimension(ntime) :: devbeam
      integer*4,dimension(ntime) :: navbeam
      end module var_beamdat
!
      module commonblocks
      real*8,allocatable :: c(:,:,:,:),wk(:),copy(:,:),bkx(:),bky(:), &
           cw(:,:,:,:),wkw(:),copyw(:,:),bwx(:),bwy(:), &
           sifprw(:),bwprw(:),cwprw(:),dwprw(:),psirz(:,:), &
           sfprw(:),sprwp(:),wgridpc(:),rfcpc(:,:),ct(:,:,:,:), &
           wkt(:),bkrt(:),bkzt(:),psiold(:),psipold(:),psipp(:), &
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
!
!vasorg      character*4 limloc
!vasorg      character*5 mfvers
!vasorg      character*10 vsname,mpnam2,lpname
!vasorg      character*60 ifname
!vasorg      character filimt*50,cshot*6,mfitpop*12,fitzts*2
!vasorg      character*80 alternate_pointname_file
!vasorg      double precision gsum,tvfbrt
!vasorg      logical vfeed,symmetrize,backaverage,fitdelz,fitsiref,writepc
!this should be there in ecomdu1      logical vfeed,symmetrize,backaverage
!vasorg      logical scalea,fitts,ztserr,oldccomp,oldcomp
!vasorg      logical fitts,ztserr,oldccomp,oldcomp
!vasorg      logical*4 do_spline_fit
!vasorg      integer erflag
!vasorg      integer*4 use_alternate_pointnames
!vasorg      logical write_Kfile
!vasorg      logical fitfcsum
!vasorg      character*82 snap_file
