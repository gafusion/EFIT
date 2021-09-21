!**********************************************************************
!! 
!>   block data routine to hold all data statements for
!!   variables that are in common blocks.
!!
!**********************************************************************
      block data efit_bdata
      ! not used here
!      use commonblocks,only: zeros,xouts,youts,bpoo,bpooz,bpooc
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*8 (i-n), real*8 (a-h,o-z)
      include 'env2d.inc'
      integer*8 limtrs
      common/wwork1/xlims(5),ylims(5),limtrs,xlmins
      data limtrs/5/
      data iunit/35/, m_write/1/, m_read/1/
      data out2d /'curve2d.dat'/, out2d_bin/'curve2d_bin.dat'/
      end

