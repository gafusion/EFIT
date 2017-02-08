!
! $Date: 2008/07/29 23:42:31 $ $Author: radhakri $
!
! @(#)$RCSfile: curve2d_var.inc,v $ $Revision: 1.1.2.2 $
!
	dimension x(ndim,ncrv), y(ndim,ncrv), bshft(2), thcrv(ncrv), &
      		sclpc(ncrv), sx(ncrv,ncrv), &
      		sy(ncrv,ncrv), sangle(ncrv), sgap(ncrv), &
      		xfm(ncrv), yfm(ncrv), xto(ncrv), yto(ncrv), rat(ndim), &
      		tlen(ncrv), zmat(ix,iy), &
                xmpos(mdim), ympos(mdim), anum(mdim), hgt(mdim)
	
	integer	nplt(ncrv), ncnct(ncrv), ncurve, &
         	pltlen, xnlen, ynlen, xtck, ytck, ixnon, iynon,  &
      		intax, intay, iaxis, igridx, igridy, ibrdr,  &
      		idash, idot, ichdot, ichdsh, &
      		markme(ncrv), dotme(ncrv), dashme(ncrv), chdtme(ncrv), &
      		chdhme(ncrv), nslen, isaxs, &
      		ngaps(ncrv), nsxy(ncrv), nshd, nvec, ivec(ncrv), iorel, &
      		icont(ncrv), nword, iy, ix, line, ithk, ipri, nline, &
      		mrc(ncrv), nmrk(ncrv), iexit, ipag

	character*20 ptitle, xtitle, ytitle, clearx(ncrv), sname
	character*20	mode, lbflg, draw

     	integer	imes(mdim), inum(mdim), note(mdim), iplce(mdim), m
	character*72	lmes(mdim)
	integer	iunit, m_write, m_read

