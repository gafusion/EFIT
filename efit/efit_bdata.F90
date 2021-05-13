!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          block data routine to hold all data statements for      **
!**          variables that are in common blocks.                    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**                91..........first created, John Ferron            **
!**          11/11/93..........revised                               **
!**          10/01/97..........revised by Q.Peng,.ddd specific data  **
!**                            are moved to a different file         **
!**                                                                  **
!**********************************************************************
      block data efit_bdata
      use commonblocks,only: zeros,xouts,youts,bpoo,bpooz,bpooc
      include 'eparmdud129.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      include 'env2d.inc'
      common/wwork1/xlims(5),ylims(5),limtrs,xlmins
      data limtrs/5/
      data iunit/35/, m_write/1/, m_read/1/
      data out2d /'curve2d.dat'/, out2d_bin/'curve2d_bin.dat'/
      end

