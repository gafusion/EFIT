!*************************************************************************
!**                                                                   **
!** MODULE FILE: curve2d_mod.f90                                        **
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
!**
!*************************************************************************
      module curve2d_mod
!      use eparm, only: ndim

      integer*4, parameter:: ndimc=3200 ! must be consistent with ndim in eparm
      integer*4, parameter:: ncrv=180
      integer*4, parameter:: mdim=300

      real*8 xx(ndimc,ncrv), yy(ndimc,ncrv), &
             sxx(ncrv,ncrv), syy(ncrv,ncrv)

      character*24 ptitle, xtitle, ytitle, sname, clearx(ncrv)
      character*72 lmes(mdim)
      character*24 lbflg, draw, mode, mxalf(mdim)

      integer*4 nn, ncurve, ipag, ibrdr, iorel, nplen, nxlen, nylen
      real*8 grce, xphy, yphy, xorl, yorl, hight, bngle, bshft(2), &
             xlen, ylen, xomin, xstp, xomax, yomin, ystp, yomax
      integer*4 iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
                ixaxs, nslen, igridx, igridy, idash, idot, ichdsh, ichdot
      real*8 sorg, stp, smax, slen, xps, yps
      integer*4 ndshme(ncrv), ndotme(ncrv), ncdhme(ncrv), ncdtme(ncrv), &
                markme(ncrv), mrc(ncrv), nmrk(ncrv), nxy(ncrv), ncnct(ncrv)
      real*8 thcrv(ncrv), sclpc(ncrv), tlen(ncrv)
      integer*4 icont(ncrv), nword, ix, iy, line, ithk, ipri, nline
      real*8 zinc
      integer*4 nshd, nsxy(ncrv), ngaps(ncrv), nvec, ivec(ncrv)
      real*8 sangle(ncrv), sgap(ncrv), xfm(ncrv), yfm(ncrv), xto(ncrv), &
             yto(ncrv)
      integer*4 msg, note(mdim), imes(mdim), iplce(mdim), inum(mdim), iexit
      real*8 anum(mdim), xpos(mdim), ypos(mdim), ht(mdim) 
      
      end module curve2d_mod
