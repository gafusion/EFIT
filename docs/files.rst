Output files
============

g-files
-------
Output from file mode or snap(_ext) mode (1,2,3,7) runs when ``KEQDSK<2`` (from ``efitin`` or ``in1`` namelists).
Contains mostly arrays of calculated values such as flux on the grid, *ff'*, *p'*, and *q* profiles.
The detailed format for g-files can be found in the Fortran source code ``write_g.f90``. 
Briefly, a right-handed cylindrical coordinate system (*R*, :math:`{\phi}` , *Z*) is used. The g-file provides
profile information on on a uniform flux grid from the magnetic axis to the plasma boundary
and the poloidal flux function on the rectangular computation grid. Information on the plasma
boundary and the surrounding limiter contour are also provided.  The following Fortan code can
be used to read the files when they are in ASCII format (they can also be written in binary
with the same structure using single precision reals and integers):

.. code-block:: fortran

         character*10 case(6)
         integer*4 i,idum,nw,nh,nbbbs,limitr,kvtor,nmass,keecur, &
                ishot,itime, ! already defined: iplcout,nfcoil,nesum
         real*8 rdim,zdim,rcentr,rleft,zmid,rmaxis,zmaxis,simag, &
                sibry,bcentr,current,xdum,rvtor
         real*8 brsp(nfcoil),ecurrt(nesum)
         real*8,dimension(:),allocatable :: fpol,pres,ffprim,pprime, &
                qpsi,rbbbs,zbbbs,rlim,zlim,pressw,pwprim,dmion,rhovn, &
                epoten,rgrid,zgrid,brsp,ecurrt,pcurrt,pcurrz
         real*8,dimension(:,:),allocatable :: psirz

         read (neqdsk,2000) (case(i),i=1,6),idum,nw,nh
         allocate(psirz(nw,nh),fpol(nw),pres(nw),ffprim(nw), &
           pprime(nw),qpsi(nw),pressw(nw),pwprim(nw), &
                dmion(nw),rhovn(nw),epoten(nw),rgrid(nw),zgrid(nw), &
                pcurrt(nw*nh),pcurrz(nw,nh)
         read (neqdsk,2020) rdim,zdim,rcentr,rleft,zmid
         read (neqdsk,2020) rmaxis,zmaxis,simag,sibry,bcentr
         read (neqdsk,2020) current,simag,xdum,rmaxis,xdum
         read (neqdsk,2020) zmaxis,xdum,sibry,xdum,xdum
         read (neqdsk,2020) (fpol(i),i=1,nw)
         read (neqdsk,2020) (pres(i),i=1,nw)
         read (neqdsk,2020) (ffprim(i),i=1,nw)
         read (neqdsk,2020) (pprime(i),i=1,nw)
         read (neqdsk,2020) ((psirz(i,j),i=1,nw),j=1,nh)
         read (neqdsk,2020) (qpsi(i),i=1,nw)
         read (neqdsk,2022) nbbbs,limitr
         allocate(rbbbs(nbbbs),zbbbs(nbbbs),rlim(limitr),zlim(limitr))
         read (neqdsk,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)
         read (neqdsk,2020) (rlim(i),zlim(i),i=1,limitr)
         read (neqdsk,2024) kvtor,rvtor,nmass
         if (kvtor.gt.0) then
           read (neqdsk,2020) (pressw(i),i=1,nw)
           read (neqdsk,2020) (pwprim(i),i=1,nw)
         endif
         if (nmass.gt.0) then
           read (neqdsk,2020) (dmion(i),i=1,nw)
         endif
         read (neqdsk,2020) (rhovn(i),i=1,nw)
         read (neqdsk,2026) keecur
         if (keecur.gt.0) then
           read (neqdsk,2020) (epoten(i),i=1,nw)
         endif
         ! note: unlike the rest of the file, these optional extras have
         !       never been described completely with variables available here
         !       (since being added initially in 2004)
         if (iplcout.gt.0) then
           if (iplcout.eq.1) then
             if (ishot.le.99999) then
               read (neqdsk,3000) nw,nh,ishot,itime
             else
               read (neqdsk,3003) nw,nh,ishot,itime
             endif
             read (neqdsk,2020) rgrid(1),rgrid(nw),zgrid(1),zgrid(nh)
             read (neqdsk,2020) (brsp(i),i=1,nfcoil)  ! also in m and a-files
             read (neqdsk,2020) (ecurrt(i),i=1,nesum) ! also in m and a-files
             read (neqdsk,2020) (pcurrt(i),i=1,nw*nh)
           elseif (iplcout.eq.2) then
             read (neqdsk,2020) ((pcurrz(i),i=1,nw),j=1,nh)
           endif
         endif

    2000 format (6a8,3i4)
    2020 format (5e16.9)
    2022 format (2i5)
    2024 format (i5,e16.9,i5)
    2026 format (i5)
    3000 format (4i5)
    3003 format (2i5,i6,i5)

.. csv-table:: gEQDSK variables
   :file: tables/geqdsk.csv
   :widths: 20,80
   :header-rows: 1


a-files
-------

Output from file mode or snap(_ext) mode (1,2,3,7) when ``ICONVR>=0`` (from ``efitin``
or ``in1`` namelist).  Contains shape, convergence, and other global parameters .
Most values are scalar. The detailed format for a-files can be found in the
Fortran source code ``write_a.f90``.  The following Fortan code can be used to read
the files when they are in ASCII format, assuming it is called in a loop over
timeslices with variable sizes already defined (they can also be written in
binary with the same structure using single precision reals and integers)

.. code-block:: fortran

       character limloc*4,qmflag*3 
       character header*42,qmflag*3,fit_type*3
       integer*4 nlold,nlnew ! already defined: jj,magpri,magpri0,nsilop,nesum,ntime
       data nlold/40/,nlnew/41/
       integer*4 jflag(ntime),jerror(ntime),
       real*8 time(ntime),elong(ntime),rout(ntime),zout(ntime),utri(ntime), &
         ltri(ntime),aminor(ntime),volume(ntime),betat(ntime),gaptop(ntime), &
         betap(ntime),li(ntime),gapin(ntime),gapout(ntime),qsta(ntime), &
         rcurrt(ntime),zcurrt(ntime),qout(ntime),sepin(ntime), &
         sepout(ntime),septop(ntime),sibdry(ntime),area(ntime), &
         wmhd(ntime),elongm(ntime),qm(ntime),terror(ntime), &
         rm(ntime),zm(ntime),sepbot(ntime),sepbot(ntime), &
         alpha(ntime),rttt(ntime),dbpli(ntime),delbp(ntime),oring(ntime), &
         rseps(2,ntime),zseps(2,ntime),sepexp(ntime),shearb(ntime), &
         xtch(ntime),ytch(ntime),q95(ntime),vertn(ntime),aaq1(ntime), &
         aaq2(ntime),aaq3(ntime),btaxp(ntime),btaxv(ntime), &
         psim(ntime),dsep(ntime), &
         wbpol(ntime),taumhd(ntime),betapd(ntime),betatd(ntime), &
         li3(ntime),wdia(ntime),taudia(ntime),wbpold(ntime), &
         qmerci(ntime),slantu(ntime),slantl(ntime),zeff(ntime), &
         zeffr(ntime),tave(ntime),rvsin(ntime),zvsin(ntime), &
         rvsout(ntime),zvsout(ntime),wpdot(ntime),wbdot(ntime), &
         vsurfa(ntime),cjor95(ntime),pp95(ntime),drsep(ntime), &
         yyy2(ntime),xnnc(ntime),ipmeas,betatn,psiq1,betat2, &
         wtherm(ntime),wfbeam(ntime),taujd3(ntime),tauthn(ntime) &
         qsiwant(ntime),cjorsw(ntime),cjor0(ntime), &
         ssiwant(ntime),ssi95(ntime),rexpan,fexpan,qmin,fexpvs,shearc, &
         sepnose,ssi01,znose,rhoqmin,peak(ntime),dminux(ntime), &
         dminlx(ntime),dolubat(ntime),dolubafm(ntime),diludom(ntime), &
         diludomm(ntime),ratsol(ntime),rvsiu(ntime),zvsiu(ntime), &
         rvsid(ntime),zvsid(ntime),rvsou(ntime),zvsou(ntime), &
         rvsod(ntime),zvsod(ntime),condno(ntime),psin32(ntime), &
         psin21(ntime),rq32in(ntime),rq21top(ntime),chilibt(ntime), &
            xbetapr,tflux(ntime),tchimls,twagap(ntime)
       real*8 csilop(nsilop,ntime),cmpr2(magpri,ntime), &
         ccbrsp(nfcoil,ntime),eccurt(nesum,ntime)

       read (neqdsk,1060) time(jj),jflag(jj),lflag,limloc(jj), &
         mco2v,mco2r,qmflag,nlold,nlnew
       read (neqdsk,1040) chisq(jj),rcencm,bcentr(jj),ipmeas(jj)
       read (neqdsk,1040) ipmhd(jj),rout(jj),zout(jj),aminor(jj)
       read (neqdsk,1040) elong(jj),utri(jj),ltri(jj),volume(jj)
       read (neqdsk,1040) rcurrt(jj),zcurrt(jj),qsta(jj),betat(jj)
       read (neqdsk,1040) betap(jj),li(jj),gapin(jj),gapout(jj)
       read (neqdsk,1040) gaptop(jj),gapbot(jj),q95(jj),vertn(jj)
       read (neqdsk,1040) (rco2v(k,jj),k=1,mco2v)
       read (neqdsk,1040) (dco2v(jj,k),k=1,mco2v)
       read (neqdsk,1040) (rco2r(k,jj),k=1,mco2r)
       read (neqdsk,1040) (dco2r(jj,k),k=1,mco2r)
       read (neqdsk,1040) shearb(jj),bpolav(jj),s1(jj),s2(jj)
       read (neqdsk,1040) s3(jj),qout(jj),sepin(jj),sepout(jj)
       read (neqdsk,1040) septop(jj),sibdry(jj),area(jj),wmhd(jj)
       read (neqdsk,1040) terror(jj),elongm(jj),qm(jj),cdflux(jj)
       read (neqdsk,1040) alpha(jj),rttt(jj),psiref(jj),indent(jj)
       read (neqdsk,1040) rseps(1,jj),zseps(1,jj),rseps(2,jj), &
         zseps(2,jj)
       read (neqdsk,1040) sepexp(jj),sepbot(jj),btaxp(jj),btaxv(jj)
       read (neqdsk,1040) aaq1(jj),aaq2(jj),aaq3(jj),dsep(jj)
       read (neqdsk,1040) rm(jj),zm(jj),psim(jj),taumhd(jj)

       fluxx=diamag(jj)*1.0e-03
       read (neqdsk,1040) betapd(jj),betatd(jj),wdia(jj),fluxx
       read (neqdsk,1040) vloopt(jj),taudia(jj),qmerci(jj),tavem
       read (neqdsk, 1041) nsilop0,magpri0,nfcoil0,nesum0
       read (neqdsk,1040) (csilop(k,jj),k=1,nsilop0), &
         (cmpr2(k,jj),k=1,magpri0)
       read (neqdsk,1040) (ccbrsp(k,jj),k=1,nfcoil0)
       read (neqdsk,1040) (eccurt(jj,k),k=1,nesum0)
       read (neqdsk,1040) pbinj(jj),rvsin(jj),zvsin(jj),rvsout(jj)
       read (neqdsk,1040) zvsout(jj),vsurfa(jj),wpdot(jj),wbdot(jj)
       read (neqdsk,1040) slantu(jj),slantl(jj),zuperts(jj),chipre
       read (neqdsk,1040) cjor95(jj),pp95(jj),drsep(jj),yyy2(jj)
       read (neqdsk,1040) xnnc(jj),cprof,oring(jj),cjor0(jj)
       read (neqdsk,1040) fexpan,qmin,chigamt,ssi01
       read (neqdsk,1040) fexpvs,sepnose,ssi95(jj),rhoqmin
       read (neqdsk,1040) cjor99(jj),cj1ave(jj),rmidin(jj),rmidout(jj)
       read (neqdsk,1040) psurfa(jj), peak(jj),dminux(jj),dminlx(jj)
       read (neqdsk,1040) dolubaf(jj),dolubafm(jj),diludom(jj),diludomm(jj)
       read (neqdsk,1040) ratsol(jj),rvsiu(jj),zvsiu(jj),rvsid(jj)
       read (neqdsk,1040) zvsid(jj),rvsou(jj),zvsou(jj),rvsod(jj)
       read (neqdsk,1040) zvsod(jj),condno(jj),psin32(jj),psin21(jj)
       read (neqdsk,1040) rq32in(jj),rq21top(jj),chilibt(jj),xdum
       read (neqdsk,1040) xbetapr,tflux(jj),tchimls,twagap(jj)
       read (neqdsk,1042) header,fit_type

  1040 format (1x,4e16.9)
  1041 format (1x,4i5)
  1042 format (1x,a42,1x,a3)
  1060 format (1h*,f7.2,10x,i5,11x,i5,1x,a3,1x,i3,1x,i3,1x,a3,1x,2i5)

.. csv-table:: aEQDSK variables
   :file: tables/aeqdsk.csv
   :widths: 20,80
   :header-rows: 1


m-files 
------- 

Output from file mode or snap(_ext) mode (1,2,3,7) runs when specified by
``IOUT`` (see ``efitin`` or ``in1`` namelist).  Contains all diagnostic data,
uncertainties, and synthetic measurements which can be used as input to fitting
solutions and the quality of the fits (chi squared).  Also contains several
global quality of fit parameters, plasma coefficients and coil currents used for
correcting magnetic measurements.

.. csv-table:: mEQDSK variables
   :file: tables/meqdsk.csv
   :widths: 20,80
   :header-rows: 1
