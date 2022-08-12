.. highlight:: fortran

EFIT files
==========

gEQDSK
-----------------------------------------
Contains mainly arrays of calculated values such as flux on grid, ff', p', q profile.
Output from file mode or snap(_ext) mode (1,2,3,7) runs when KEQDSK<2 (from efitin or in1 namelists).  The detailed format for G EQDSK can be found in the Fortran source code
weq.f90. Briefly, a right-handed cylindrical coordinate system (R, f, Z) is used. The G EQDSK provides
information on the pressure, poloidal current function, q profile on a uniform flux
grid from the magnetic axis to the plasma boundary and the poloidal flux
function on the rectangular computation grid. Information on the plasma
boundary and the surrounding limiter contour in also provided :: 

	 character*10 case(6)
	 dimension psirz(nw,nh),fpol(1),pres(1),ffprim(1),
	 . pprime(1),qpsi(1),rbbbs(1),zbbbs(1),
	 . rlim(1),zlim(1)
	c
	 read (neqdsk,2000) (case(i),i=1,6),idum,nw,nh
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
	c
	 2000 format (6a8,3i4)
	 2020 format (5e16.9)
	 2022 format (2i5)
	 2024 format (i5,e16.9,i5)

.. csv-table:: gEQDSK variables
   :file: tables/geqdsk.csv
   :widths: 20,80
   :header-rows: 1


aeqdsk 
---------------------------------------

Contains shape, convergence, and other global parameters, mainly scalar values.  It is output from file mode or snap(_ext) mode (1,2,3,7) when ICONVR>=0 (from efitin or in1 namelist).  aeqdsk has the following fortran format :: 

	parameter (magpri67=29,magpri322=31,magprirdp=8)
	parameter (magpri=magpri67+magpri322+magprirdp)
	parameter (nsilop=41,nesum=6)
	parameter (ntime=400)
	character limloc*4,qmflag*3
	data nlold/35/,nlnew/39/
	dimension time(ntime),jflag(ntime),
	. eout(ntime),rout(ntime),zout(ntime),doutu(ntime),
	. doutl(ntime),aout(ntime),vout(ntime),betat(ntime),otop(ntime),
	. betap(ntime),ali(ntime),oleft(ntime),oright(ntime),qsta(ntime),
	. rcurrt(ntime),zcurrt(ntime),qout(ntime),olefs(ntime),
	. orighs(ntime),otops(ntime),sibdry(ntime),areao(ntime),
	. wplasm(ntime),elongm(ntime),qqmagx(ntime),terror(ntime),
	. rmagx(ntime),zmagx(ntime),obott(ntime),obots(ntime),
	. alpha(ntime),rttt(ntime),dbpli(ntime),delbp(ntime),oring(ntime),
	. rseps(2,ntime),zseps(2,ntime),sepexp(ntime),shearb(ntime),
	. xtch(ntime),ytch(ntime),qpsib(ntime),vertn(ntime),aaq1(ntime),
	. aaq2(ntime),aaq3(ntime),btaxp(ntime),btaxv(ntime),
	. simagx(ntime),jerror(ntime),seplim(ntime),
	. wbpol(ntime),taumhd(ntime),betapd(ntime),betatd(ntime),
	. alid(ntime),wplasmd(ntime),taudia(ntime),wbpold(ntime),
	. qmerci(ntime),slantu(ntime),slantl(ntime),zeff(ntime),
	. zeffr(ntime),tave(ntime),rvsin(ntime),zvsin(ntime),
	. rvsout(ntime),zvsout(ntime),wpdot(ntime),wbdot(ntime),
	. vsurfa(ntime),cjor95(ntime),pp95(ntime),ssep(ntime),
	. yyy2(ntime),xnnc(ntime),pasman,betatn,psiq1,betat2,
	. wtherm(ntime),wfbeam(ntime),taujd3(ntime),tauthn(ntime)
	. qsiwant(ntime),cjorsw(ntime),cjor0(ntime),
	. ssiwant(ntime),ssi95(ntime),rexpan,fexpan,qqmin,fexpvs,shearc,
	. sepnose,ssi01,znose,rqqmin,peak(ntime),dminux(ntime),
	. dminlx(ntime),dolubat(ntime),dolubafm(ntime),diludom(ntime),
	. diludomm(ntime),ratsol(ntime),rvsiu(ntime),zvsiu(ntime),
	. rvsid(ntime),zvsid(ntime),rvsou(ntime),zvsou(ntime),
	. rvsod(ntime),zvsod(ntime),condno(ntime),psin32(ntime),
	. psin21(ntime),rq32in(ntime),rq21top(ntime),chilibt(ntime)
	dimension csilop(nsilop,ntime),cmpr2(magpri,ntime)
	dimension ccbrsp(nfcoil,ntime),eccurt(nesum,ntime)
	c
	read (neqdsk,1060) time(jj),jflag(jj),lflag,limloc(jj),
	. mco2v,mco2r,qmflag,nlold,nlnew
	read (neqdsk,1040) tsaisq(jj),rcencm,bcentr(jj),pasmat(jj)
	read (neqdsk,1040) cpasma(jj),rout(jj),zout(jj),aout(jj)
	read (neqdsk,1040) eout(jj),doutu(jj),doutl(jj),vout(jj)
	read (neqdsk,1040) rcurrt(jj),zcurrt(jj),qsta(jj),betat(jj)
	read (neqdsk,1040) betap(jj),ali(jj),oleft(jj),oright(jj)
	read (neqdsk,1040) otop(jj),obott(jj),qpsib(jj),vertn(jj)
	read (neqdsk,1040) (rco2v(k,jj),k=1,mco2v)
	read (neqdsk,1040) (dco2v(jj,k),k=1,mco2v)
	read (neqdsk,1040) (rco2r(k,jj),k=1,mco2r)
	read (neqdsk,1040) (dco2r(jj,k),k=1,mco2r)
	read (neqdsk,1040) shearb(jj),bpolav(jj),s1(jj),s2(jj)
	read (neqdsk,1040) s3(jj),qout(jj),olefs(jj),orighs(jj)
	read (neqdsk,1040) otops(jj),sibdry(jj),areao(jj),wplasm(jj)
	read (neqdsk,1040) terror(jj),elongm(jj),qqmagx(jj),cdflux(jj)
	read (neqdsk,1040) alpha(jj),rttt(jj),psiref(jj),xndnt(jj)
	read (neqdsk,1040) rseps(1,jj),zseps(1,jj),rseps(2,jj)
	. ,zseps(2,jj)
	read (neqdsk,1040) sepexp(jj),obots(jj),btaxp(jj),btaxv(jj)
	read (neqdsk,1040) aaq1(jj),aaq2(jj),aaq3(jj),seplim(jj)
	read (neqdsk,1040) rmagx(jj),zmagx(jj),simagx(jj),taumhd(jj)
	fluxx=diamag(jj)*1.0e-03
	read (neqdsk,1040) betapd(jj),betatd(jj),wplasmd(jj),fluxx
	read (neqdsk,1040) vloopt(jj),taudia(jj),qmerci(jj),tavem
	read (neqdsk, 1041) nsilop0,magpri0,nfcoil0,nesum0
	read (neqdsk,1040) (csilop(k,jj),k=1,nsilop0),
	. (cmpr2(k,jj),k=1,magpri0)
	read (neqdsk,1040) (ccbrsp(k,jj),k=1,nfcoil0)
	read (neqdsk,1040) (eccurt(jj,k),k=1,nesum0)
	read (neqdsk,1040) pbinj(jj),rvsin(jj),zvsin(jj),rvsout(jj)
	read (neqdsk,1040) zvsout(jj),vsurfa(jj),wpdot(jj),wbdot(jj)
	read (neqdsk,1040) slantu(jj),slantl(jj),zuperts(jj),chipre
	read (neqdsk,1040) cjor95(jj),pp95(jj),ssep(jj),yyy2(jj)
	read (neqdsk,1040) xnnc(jj),cprof,oring(jj),cjor0(jj)
	read (neqdsk,1040) fexpan,qqmin,chigamt,ssi01
	read (neqdsk,1040) fexpvs,sepnose,ssi95(jj),rqqmin
	read (neqdsk,1040) cjor99(jj),cj1ave(jj),rmidin(jj),rmidout(jj)
	read (neqdsk,1040) psurfa(jj), peak(jj),dminux(jj),dminlx(jj)
	Lao 9/10/04
	read (neqdsk,1040) dolubaf(jj),dolubafm(jj),diludom(jj),diludomm(jj)
	read (neqdsk,1040) ratsol(jj),rvsiu(jj),zvsiu(jj),rvsid(jj)
	read (neqdsk,1040) zvsid(jj),rvsou(jj),zvsou(jj),rvsod(jj)
	read (neqdsk,1040) zvsod(jj),condno(jj),psin32(jj),psin21(jj)
	read (neqdsk,1040) rq32in(jj),rq21top(jj),chilibt(jj),xdum
	c
	1040 format (1x,4e16.9)
	1041 format (1x,4i5)
	1060 format (1h*,f7.2,10x,i5,11x,i5,1x,a3,1x,i3,1x,i3,1x,a3,1x,2i5)

.. csv-table:: aEQDSK variables
   :file: tables/aeqdsk.csv
   :widths: 20,80
   :header-rows: 1


mEQDSK
------
Output from file mode or snap(_ext) mode (1,2,3,7) runs when specified by IOUT (see efitin or in1 namelist).  Contains all diagnostic data, uncertainties, and synthetic measurements which can be used as input to fitting solutions and the quality of the fits (chi squared).  Also contains several global quality of fit parameters, plasma coefficients and coil currents used for correcting magnetic measurements.

.. csv-table:: mEQDSK variables
   :file: tables/meqdsk.csv
   :widths: 20,80
   :header-rows: 1



