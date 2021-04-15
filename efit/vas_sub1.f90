      subroutine inp_file_ch(nw,nh,ch1,ch2)
        
        integer,     intent(in)  :: nw, nh 
        character*4, intent(out) :: ch1, ch2
        
        write(ch1,'(i4)') nw
        ch1 = adjustl(ch1)
        write(ch2,'(i4)') nh
        ch2 = adjustl(ch2)
        
        return
        
      end subroutine inp_file_ch
