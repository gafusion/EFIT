! iunit  logical unit number
! m_write  1 = Write plot parameters in ASCII format
!          0 = Write plot parameters in Binary format
! m_read   1 = Read plot parameters in ASCII format
!          0 = Read plot parameters in Binary format
!
! $Date: 2008/07/29 23:42:32 $ $Author: radhakri $
!
! @(#)$RCSfile: env2d.inc,v $ $Revision: 1.1.2.2 $ 
!

      common / flg2d / iunit, m_write, m_read, out2d, out2d_bin
      character*72 out2d, out2d_bin

