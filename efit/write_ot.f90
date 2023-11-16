!**********************************************************************
!>
!!    writes out time history of plasma parameters in o-file format
!!    
!!
!!    @param kwtime : number of time elements 
!!
!********************************************************************** 
      subroutine write_ot(kwtime)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      character*30 ofname, chname
      Character*28 xxtitle,yytitle,zztitle
      data czero/0.0/
!
      if(kwtime.le.1) return
      if(kwritime.ne.2) return
!
      xdum=0.0
      ydum=0.0
      itime00=time(1)
      call setfnmt('ot',ishot,itime00,ofname)
!
      ofname=ofname(1:14)//'_chi2mag'
      open(unit=74,status='old',file=ofname,iostat=ioerr)
      if(ioerr.eq.0) close(unit=74,status='delete')
      open(unit=74,status='new',file=ofname)
      xxtitle='Time(ms)'
      yytitle='Magnetic Chi Square'
      zztitle='EFIT'
      write (74,93024) xxtitle
      write (74,93024) yytitle
      write (74,93024) zztitle
      do i=1,kwtime
        if(kerrot(i).eq.0) &
          write (74,92924) time(i),chisq(i),xdum,xdum
      enddo
      close(unit=74)
!
      if (kstark.gt.0) then
        ofname=ofname(1:14)//'_chi2mse'
        open(unit=74,status='old',file=ofname,iostat=ioerr)
        if(ioerr.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=ofname)
        xxtitle='Time(ms)'
        yytitle='MSE Chi Square'
        zztitle='EFIT'
        write (74,93024) xxtitle
        write (74,93024) yytitle
        write (74,93024) zztitle
        do i=1,kwtime
          if(kerrot(i).eq.0) &
            write (74,92924) time(i),chi2gamt(i),xdum,xdum
        enddo
        close(unit=74)
      endif
!
      if (mmbmsels.gt.0) then
        ofname=ofname(1:14)//'_chi2mls'
        open(unit=74,status='old',file=ofname,iostat=ioerr)
        if(ioerr.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=ofname)
        xxtitle='Time(ms)'
        yytitle='MSE-LS Chi Square'
        zztitle='EFIT'
        write (74,93024) xxtitle
        write (74,93024) yytitle
        write (74,93024) zztitle
        do i=1,kwtime
          if(kerrot(i).eq.0) &
            write (74,92924) time(i),chi2mls(i),xdum,xdum
        enddo
        close(unit=74)
!
        do j=1,nmsels
          if (sinmls(j).gt.0.0) then
            ichan00=j
            call setfnmt('ot',ishot,ichan00,chname)
            ofname=ofname(1:14)//'_echan'//chname(13:14)
            open(unit=74,status='old',file=ofname,iostat=ioerr)
            if(ioerr.eq.0) close(unit=74,status='delete')
            open(unit=74,status='new',file=ofname)
            xxtitle='Time(ms)'
            yytitle='Measured B MSE-LS (T)'
            zztitle='MSE-LS Channel '//chname(13:14)
            write (74,93024) xxtitle
            write (74,93024) yytitle
            write (74,93024) zztitle
            do i=1,kwtime
              if (kerrot(i).eq.0) then
                iges=i
                ydum=sbmselt(iges,j)
                if(swtbmselt(iges,j).gt.1.e-06_dp) ydum=ydum/ &
                                                       swtbmselt(iges,j)
                write (74,92924) time(i),bmselt(iges,j),xdum,ydum
              endif
            enddo
            close(unit=74)
          endif
        enddo
!
        do j=1,nmsels
          if (sinmls(j).gt.0.0) then
            ichan00=j
            call setfnmt('ot',ishot,ichan00,chname)
            ofname=ofname(1:14)//'_cchan'//chname(13:14)
            open(unit=74,status='old',file=ofname,iostat=ioerr)
            if(ioerr.eq.0) close(unit=74,status='delete')
            open(unit=74,status='new',file=ofname)
            xxtitle='Time(ms)'
            yytitle='Computed B MSE-LS (T)'
            zztitle='MSE-LS Channel '//chname(13:14)
            write (74,93024) xxtitle
            write (74,93024) yytitle
            write (74,93024) zztitle
            do i=1,kwtime
              if (kerrot(i).eq.0) then
                iges=i
                write (74,92924) time(i),cmmls(iges,j),xdum,xdum
              endif
            enddo
            close(unit=74)
          endif
        enddo
!
      endif
      return
92924 format (4(1pe12.5,1x))
93024 format (a28)
      end subroutine write_ot
