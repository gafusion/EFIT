!*************************************************************************
!**                                                                   **
!** INCLUDE FILE: curve2d.inc                                        **
!**                                                                   **
!**     FILE DESCRIPTION:                                        **
!**       Declares and defines dimensioned variables, and common  **
!**       block variables used by plot parameters   **
!**                                                                   **
!**     VARIABLE DESCRIPTION:      **
!**       ncurve  Number of curves    **
!**       ipag page dimension flag    **
!**       ibrdr   Page border flag                **
!**       grce >= 0 Enable grace margin            **
!**       xphy    X-coordinate of physical origin            **
!**       yphy    Y-coordinate of physical origin            **
!**       hight   Font size                    **
!**       bngle   Base rotation angle                **
!**       bshft   Base translation                **
!**       ptitle  Plot title                    **
!**       pltlen  Length of plot title                           **
!**       xtitle  X-axis name                                 **
!**       xnlen   Length of x-axis name                           **
!**       ytitle  Y-axis name                                   **
!**       ynlen   Length of y-axis name                **
!**       xlen    Length of x-axis legend                **
!**       ylen    Length of y-axis legend                         **
!**       iaxis   Axes flag                    **
!**       ixtck   X-axis tick marks                **
!**       iytck   Y-axis tick marks                **
!**       ixnon   X-axis tick marks or labels flag        **
!**       iynon   Y-axis tick marks or labels flag        **
!**       intax   Trailing zeroes on x-axis flag            **
!**       intay   Trailing zeroes on y-axis flag            **
!**       xomin   X-value at the physical origin            **
!**       xstp    X step size in units                **
!**       xomax   Value at x-axis limit                **
!**       yomin   X-value at the physical origin                    **
!**       ystp    X step size in units                **
!**       yomax   Value at y-axis limit                **
!**       iorel   New physical origin relative to current     **
!**       origin  flag                    **
!**       xorl    X-coordinate of relative origin            **
!**       yorl    Y-coordinate of relative origin            **
!**       igridx  Number of grid lines per step on x-axis             **
!**       igridy  Number of grid lines per step on y-axis             **
!**       idash   Dash grid lines flag                **
!**       idot    Dot grid lines flag                **
!**       ichdot  Chain dot grid lines flag                       **
!**       ichdsh  Chain dash grid lines flag                      **
!**       thcrv   Curve thickness        dim = ncrv             **
!**       sclpc   Scale curve marker     dim = ncrv            **
!**       ndshme  Dash curve flag        dim = ncrv        **
!**       ndotme  Dot curve flag        dim = ncrv        **
!**       nchdme  Chain dash curve flag    dim = ncrv        **
!**       nchtme  Chain ot curve flag    dim = ncrv        **
!**       markme  Curve marker             dim = ncrv             **
!**       clearx  Color            dim = ncrv        **
!**       xx      Array of x-coordinates                **
!**               dim = (nxy, ncrv)    **
!**       yy      Array of y-coordinates                **
!**               dim = (nxy, ncrv)    **
!**       nxy     Number of (x,y) points to be ploted         **
!**               dim = ncrv        **
!**       ncnct   Marker specification                **
!**       mrc     Custom interrupted line style            **
!**       tlen    overall length of the pattern            **
!**       nmrk    total number of marks and spaces        **
!**       rat     Array of ratios of marks and spaces to overall    **
!**               length                        **
!**       icont   Contour plotting flag                **
!**       nword   Number of words available in common block    **
!**       zmat    2-dimensional array containing Z data surface    **
!**               values                        **
!**       ix      X dimension of zmat                **
!**       iy      Y dimension of zmat                **
!**       zinc    Increment between z levels            **
!**       line    Index number in the group of curve         **
!**               characteristic                    **
!**       mode    'DOT' for dotted lines                **
!**               'DASH' for dashed lines              **
!**       lbflg   'NOLABELS' do not use labels             **
!**               'LABELS' use labels                **
!**       ithk    Line thickness                    **
!**       ipri    Line priority                    **
!**       draw    'DRAW' draw curves                **
!**               'NODRAW' do not draw curves            **
!**       nshd    Number of shades, shade area between 2 curves    **
!**       sxx     Array of x-coordinates                **
!**       syy     Array of y-coordinates                **
!**       nsxy    Number of (sx,sy) pairs                **
!**       sangle  angle of shading lines                **
!**       sgap    Array of shading gaps                **
!**       ngaps   Number of elements in sgaps            **
!**       nvec    Number of vectors                **
!**       xfm     X value of originating point    dim = nvec    **
!**       yfm     Y value of originating point    dim = nvec    **
!**       xto     X value of point pointed to    dim = nvec    **
!**       yto     Y value of point pointed to    dim = nvec    **
!**       ivec    Describes vector arrowhead    dim = nvec    **
!**       m       Number of text lines to be written        **
!**       note    Text string flag        dim = mdim    **
!**       lmes    Character string text        dim = mdim    **
!**       imes    Number of characters in lmes    dim = mdim    **
!**       anum    Real number text string        dim = mdim    **
!**       iplce   Number of decimal place        dim = mdim    **
!**       inum    Integer number text string    dim = mdim    **
!**       xmpos   X-position of text string    dim = mdim    **
!**       ympos   Y-position of text string    dim = mdim    **
!**       hgt     Font size of the text string    dim = mdim    **
!**       iexit   Subplot flag                    **
!**                                                                      **
!**     REFERENCES:                             **
!**          (1)  CA-DISSPLA manual                     **
!**          (2)                                                         **
!**                                                                      **
!**     RECORD OF MODIFICATION:                                          **
!**          03/26/93..........first created                             **
!**          05/10/93..........revised                      **
!**                                                                      **
!**
!**     $Date: 2009/02/12 22:43:48 $ $Author: radhakri $
!**
!**     @(#)$RCSfile: curve2d129.inc,v $ $Revision: 1.1.2.3 $
!**
!*************************************************************************

!sri-feb122009    parameter    (ndim = 700, ncrv=180, mdim = 200)
      parameter    (ncrv=180, mdim = 300)

      dimension xx(ndim,ncrv), yy(ndim,ncrv), thcrv(ncrv), sclpc(ncrv), bshft(2), &
        xpos(mdim), ypos(mdim), anum(mdim), ht(mdim), sangle(ncrv), sgap(ncrv), &
        sxx(ncrv,ncrv), syy(ncrv,ncrv), xfm(ncrv), yfm(ncrv), xto(ncrv), yto(ncrv), &
        tlen(ncrv), rat(ndim), ndotme(ncrv), ndshme(ncrv), ncdtme(ncrv), ncdhme(ncrv), &
        markme(ncrv), ncnct(ncrv), imes(mdim), note(mdim), inum(mdim), iplce(mdim), &
		ngaps(ncrv), nsxy(ncrv), ivec(ncrv), mrc(ncrv), nmrk(ncrv), nxy(ncrv)

      character*20 ptitle, xtitle, ytitle, sname, clearx(ncrv)
      character*72 lmes(mdim)
      character*20 lbflg, draw, mode, mxalf(mdim)

      common /var2d_page/ grce, xphy, yphy, xorl, yorl, hight, bngle, bshft, &
	    xlen, ylen, xomin, xstp, xomax, yomin, ystp, yomax, nn, ncurve, ipag, &
		ibrdr, nylen, iorel, nplen, nxlen, ptitle, xtitle, ytitle
      common /var2d_axis/ sorg, stp, smax, slen, xps, yps, iaxis, ixtck, iytck, &
	    ixnon, iynon, intax, intay, isaxs, nslen, igridx, igridy, idash, idot, ichdsh, ichdot, sname
      common /var2d_curve/ thcrv, sclpc, tlen, rat, ndshme, ndotme, ncdhme, ncdtme, &
        markme, mrc, nmrk, nxy, ncnct, clearx
      common /var2d_contur/ zinc, icont(ncrv), nword, ix, iy, line, &
        ithk, ipri, nline, lbflg, draw, mode
      common /var2d_shade/ sangle, sgap, xfm, yfm, xto, yto, nshd, &
	    nsxy, ngaps, nvec, ivec
      common /anote/ anum, xpos, ypos, ht, msg, note, imes, iplce, inum, iexit, lmes
