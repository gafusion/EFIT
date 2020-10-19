       
!---------------------------------------------------------------------
!--    Version 9-oct-91                                             --
!--    File subs.for                                                --
!--    Dummy subroutines for linking with efitd on the sun.         --
!---------------------------------------------------------------------
      SUBROUTINE DATE
      RETURN
      END
!     subroutine getpts
!     return
!     end
!     SUBROUTINE PTDATA
!     RETURN
!     END
!     subroutine getdia
!     return
!     end
!     subroutine stark_multi
!     return
!     end
      subroutine getzeff
      return
      end
!
!      SUBROUTINE PLTOUT
!      RETURN
!      END
!
!      SUBROUTINE ALTPLT
!      RETURN
!      END
!
      SUBROUTINE DONEPL
      RETURN
      END
      subroutine initial
      return
      end
! This routine is in libd3
      SUBROUTINE LIB$WAIT(TIME)
      RETURN
      END
!
! This routine is in libd3
       SUBROUTINE STR$UPCASE(OUT,IN)
       CHARACTER*(*) IN,OUT
       OUT=IN
       RETURN
       END
!
      SUBROUTINE LIB$DATE_TIME(STRING)
      CHARACTER*(*) STRING
      CHARACTER*24 DATETIME,FDATE
!C
      DATETIME=FDATE()
      STRING=DATETIME(1:23)
!C
      RETURN
      END
!      
!----------------------------------------------------------------------
! A couple of routines extracted from subsflow.for
     subroutine db_header(ishot,itime,header)

     character*(*) header
     header = ' '
     return
     end
!======================================================================
! This is actually a D3 routine.
      subroutine laserout(filelaser,num,qqqname)
      character*(*) filelaser,qqqname
      write(6,100)
  100 format(' Laserout is not available on this computer.')
      return
      end
!======================================================================
! DISSPLA dummy subroutines.
!
      SUBROUTINE AXSPLT
      RETURN
      END
!
      SUBROUTINE YGRAXS
      RETURN
      END
!
      SUBROUTINE XNONUM
      RETURN
      END
!
      SUBROUTINE YTICKS
      RETURN
      END
!
      SUBROUTINE GRAF
      RETURN
      END
!
      SUBROUTINE ENDGR
      RETURN
      END
!
      SUBROUTINE GRACE
      RETURN
      END
!
      SUBROUTINE YLOG
      RETURN
      END
!
      SUBROUTINE XINTAX
      RETURN
      END
!
      SUBROUTINE MESSAG
      RETURN
      END
!
      SUBROUTINE INTNO
      RETURN
      END
!
      SUBROUTINE CONTHN
      RETURN
      END
!
      SUBROUTINE BCOMON
      RETURN
      END
!
      SUBROUTINE TKF
      RETURN
      END
!
      SUBROUTINE QMS
      RETURN
      END
!
      SUBROUTINE ENDPL
      RETURN
      END
!
      SUBROUTINE RLMESS
      RETURN
      END
!
      SUBROUTINE SHADE
      RETURN
      END
!
      SUBROUTINE RLVEC
      RETURN
      END
!
      SUBROUTINE RLINT
      RETURN
      END
!
      SUBROUTINE CHNDOT
      RETURN
      END
!
      SUBROUTINE MIXALF
      RETURN
      END
!
      SUBROUTINE CONMAK
      RETURN
      END
!
      SUBROUTINE CONDIG
      RETURN
      END
!
      SUBROUTINE CONLIN
      RETURN
      END
!
      SUBROUTINE COMPRS
      RETURN
      END
!
      SUBROUTINE POPNAM
      RETURN
      END
!
      SUBROUTINE GAIMA300
      RETURN
      END
!
      SUBROUTINE TEKALL
      RETURN
      END
!
      SUBROUTINE BGNPL
      RETURN
      END
!
      SUBROUTINE SETDEV
      RETURN
      END
!
      SUBROUTINE TSEND
      RETURN
      END
!
      SUBROUTINE SCLPIC
      RETURN
      END
!
      SUBROUTINE MARKER
      RETURN
      END
!
      SUBROUTINE DOT
      RETURN
      END
!
      SUBROUTINE RLREAL
      RETURN
      END
!
      SUBROUTINE CONTUR
      RETURN
      END
!
      SUBROUTINE CONANG
      RETURN
      END
!
      SUBROUTINE NOBRDR
      RETURN
      END
!
      SUBROUTINE HEIGHT
      RETURN
      END
!
      SUBROUTINE PHYSOR
      RETURN
      END
!
      SUBROUTINE BANGLE
      RETURN
      END
!
      SUBROUTINE BSHIFT
      RETURN
      END
!
      SUBROUTINE TITLE
      RETURN
      END
!
      SUBROUTINE GRID
      RETURN
      END
!
      SUBROUTINE THKCRV
      RETURN
      END
!
      SUBROUTINE DASH
      RETURN
      END
!
      SUBROUTINE SETCLR
      RETURN
      END
!
      SUBROUTINE CURVE
      RETURN
      END
!
      SUBROUTINE RESET
      RETURN
      END
!
!
      SUBROUTINE INTAXS
      RETURN
      END
!
!
      SUBROUTINE DEVT
      RETURN
      END
!
!
      SUBROUTINE MRSCOD
      RETURN
      END
!
!
      SUBROUTINE YGRAF
      RETURN
      END
!
!
      SUBROUTINE PAGE
      RETURN
      END
!
!
      SUBROUTINE XTICKS
      RETURN
      END
!
      SUBROUTINE OREL
      RETURN
      END
!
      SUBROUTINE FRAME
      RETURN
      END
      SUBROUTINE LIMPOS(NSHOT,RLIM,RLIM180,IER)
      RETURN
      END
!
!   This routine is required if the CVS revision numbers are to
!   survive an optimization.
!
!
!   2008/01/10 23:35:07 osborne
!
      subroutine subs_d_pt_rev(i)
      CHARACTER*100 opt
      character*10 s
      if( i .eq. 0) s =  &
      '@(#)subs_d_pt.for,v 4.16\000'
      return
      end
