!+ mds_interface contains a set of f77 routines which makes the         
!  actual calls to the mds routines.  Mds routines must be called       
!  with case sensitivity preserved.  Thus, compile this file with +U    
!  under HP unix.                                                       
!-                                                                      
                                                                        
      logical function mds_start (blah) 
      ! Function to connect to Mds database                             
      ! Returns true if the library was successfully opened;            
      ! Otherwise returns false.                                        
      ! We need to have an argument to make the compiler happy!         
      ! Written: 1/28/99 by rjg                                         
      implicit none 
      integer blah 
                                                                        
      ! If started is true, we need not call MdsConnect                 
      logical started / .false. / 
      integer MdsConnect 
      integer*4 status 
                                                                        
                                                                        
                                                                        
      ! Connect to mdsip service if connection not already established. 
      ! The C function MdsConnect is a void so we should not be checking
      ! the status.  Also, we pay a penalty if we open everytime through
      ! the loop.  Only open once in a program.                         
      mds_start = started 
      if (started) then 
         continue 
      else 
         status = MdsConnect('atlas' //CHAR(0)) 
         ! Assume success                                               
         mds_start = .true. 
         started = mds_start 
      endif 
                                                                        
      if (.not. started) then 
          write (6,fmt='(A,I10)')                                       &
     &     ' mds_start: ERROR connecting to MDS; status = ', status     
      endif 
                                                                        
             ! routine mds_start                                        
      END                                           
                                                                        
                                                                        
      ! Function to open an Mds tree for a specific shot number         
      ! Returns true if the tree was successfully opened;               
      ! Otherwise returns false.                                        
      ! Assumes that connection to Mds database has been made           
      ! Written: 1/28/99 by rjg                                         
      logical function mds_open_tree (tree, ishot) 
      implicit none 
                               ! IN: Name of tree                       
      character *(*) tree 
                               ! IN: The shot number                    
      integer ishot 
      integer MdsOpen 
      integer*4 status 
                               ! Function to check status from Mds+     
      logical good_status 
                                ! Function to get length of a character 
!      integer lencstr 
                                                                        
      ! open the tree                                                   
      status = MdsOpen(tree(1:len_trim(tree))//CHAR(0), %val(ishot)) 
                                                                        
      mds_open_tree = good_status(status) 
                                                                        
      if (.not. mds_open_tree) write (6,fmt='(A,A,A,I6)')               &
     &      ' Mds_Open_Tree: ERROR OPENING TREE: ',                     &
     &        tree(1:len_trim(tree)),                                    &
     &      ' for shot ', ishot                                         
                                                                        
                ! function mds_open_tree                                
      END                                           
                                                                        
                                                                        
                                                                        
      ! Function to close the connection to Mds database                
      ! Does no work because it is not clear if we have a fortran       
      ! interface to do this.                                           
      ! Written: 1/28/99 by rjg                                         
!      function mds_close                                               
!      implicit none                                                    
!      logical mds_close                                                
!                                                                       
!      mds_close = .false.                                              
!      if (.not. mds_close) write (6,fmt='(A)')                         
!     +      ' Mds_Close: ERROR disconnecting from MDS '                
!                                                                       
!      end    ! function mds_close                                      
                                                                        
                                                                        
                                                                        
      ! Routine to return an array of reals from MDS                    
      ! This routine assumes that the proper tree has been successfully 
      ! opened.                                                         
      ! Written: 1/27/99 by rjg                                         
      subroutine mds_real_array(node, signal, maxarr, rarray, narr, ok) 
      implicit none 
                                  ! IN: the proper node name            
      character*(*) node 
                                  ! IN: the proper signal name          
      character*(*) signal 
                                  ! IN: size of rarray                  
      integer maxarr 
                                  ! OUT: array with data                
      real    rarray(maxarr) 
                                  ! OUT: number of valid points in rarra
      integer narr 
                                  ! OUT: true means success; false means
      logical ok 
                                  !      there was a problem            
                                                                        
      ! Local variables                                                 
      integer MdsValue 
      integer status, descr 
                                 ! Function to check if status return   
      logical good_status 
                                 ! from MDS function calls is okay      
      character string*200 
      character expression*200 
                                  ! Function to get length of character 
!      integer lencstr 
      integer len1,len2 
                                                                        
      ok = .false. 
                                                                        
      len1 = len_trim(node) 
      len2 = len_trim(signal) 
                                                                        
      write (string,fmt='(A,A)') node(1:len1), signal(1:len2) 
                                                                        
      len1 = len_trim(string) 
      write (expression,fmt='(A,A,A,A)')                                &
     &       'LONG(SIZE(', string(1:len1), '))', CHAR(0)                
                                                                        
      status = MdsValue(expression, descr(8,narr,0),0) 
                                                                        
      ! Make sure that we had a good status return                      
      if (.not. good_status(status)) then 
         ok = .false. 
!         print *,' Bad status from MdsValue for size = ', status,      
!     +           ' for signal ', signal                                
      ! Make sure that we have some data                                
      else if (narr .le. 0) then 
         ok = .false. 
!         print *,' Number of available data values = ', narr           
      ! Make sure that our array is big enough                          
      else if (narr .gt. maxarr) then 
         ok = .false. 
!         print *,' Your array is too small to accept all data'         
      ! Success - get the data                                          
      else 
         write (expression,fmt='(A,A)')  string(1:len1), CHAR(0) 
         status = MdsValue(expression, descr(10,rarray,maxarr,0),0) 
         if (good_status(status)) then 
            ok = .true. 
         else 
            ok = .false. 
!            print *,' Bad status from MdsValue for data = ', status,   
!     +              ' for signal ', signal                             
         endif 
      endif 
                                                                        
      END                                           
                                                                        
                                                                        
                                                                        
      ! Routine to return a single float                                
      ! This routine assumes that the proper tree has been successfully 
      ! opened.                                                         
      ! Written: 2/18/99 by rjg                                         
      subroutine mds_float (node, signal, rarg, ok) 
      implicit none 
                                  ! IN: the proper node name            
      character*(*) node 
                                  ! IN: the proper signal name          
      character*(*) signal 
                                  ! OUT: the data                       
      real    rarg 
                                  ! OUT: true means success; false means
      logical ok 
                                  !      there was a problem            
                                                                        
      ! Local variables                                                 
      integer MdsValue 
      integer status, descr 
                                 ! Function to check if status return   
      logical good_status 
                                 ! from MDS function calls is okay      
      character expression*200 
                                  ! Function to get length of character 
!      integer lencstr 
      integer len1,len2 
                                                                        
      ok = .false. 
                                                                        
      len1 = len_trim(node) 
      len2 = len_trim(signal) 
                                                                        
      write (expression,fmt='(A,A,A)') node(1:len1), signal(1:len2),    &
     &          CHAR(0)                                                 
                                                                        
         status = MdsValue(expression, descr(10,rarg,0),0) 
         if (good_status(status)) then 
            ok = .true. 
         else 
            ok = .false. 
!            print *,' Bad status from MdsValue for data = ', status,   
!     +              ' for signal ', signal                             
         endif 
                                                                        
            ! mds_float                                                 
      END                                           
                                                                        
                                                                        
                                                                        
      ! Routine to return a single integer                              
      ! This routine assumes that the proper tree has been successfully 
      ! opened.                                                         
      ! Written: 2/14/99 by rjg                                         
      subroutine mds_int (node, signal, iarg, ok) 
      implicit none 
                                  ! IN: the proper node name            
      character*(*) node 
                                  ! IN: the proper signal name          
      character*(*) signal 
                                  ! OUT: the data                       
      integer iarg 
                                  ! OUT: true means success; false means
      logical ok 
                                  !      there was a problem            
                                                                        
      ! Local variables                                                 
      integer MdsValue 
      integer status, descr 
                                 ! Function to check if status return   
      logical good_status 
                                 ! from MDS function calls is okay      
      character expression*200 
                                  ! Function to get length of character 
!      integer lencstr 
      integer len1,len2 
                                                                        
      ok = .false. 
                                                                        
      len1 = len_trim(node) 
      len2 = len_trim(signal) 
                                                                        
      write (expression,fmt='(A,A,A)') node(1:len1), signal(1:len2),    &
     &          CHAR(0)                                                 
                                                                        
         status = MdsValue(expression, descr(10,iarg,0),0) 
         if (good_status(status)) then 
            ok = .true. 
         else 
            ok = .false. 
!            print *,' Bad status from MdsValue for data = ', status,   
!     +              ' for signal ', signal                             
         endif 
                                                                        
            ! mds_int                                                   
      END                                           
                                                                        
                                                                        
                                                                        
      ! Routine to return a character string                            
      ! This routine assumes that the proper tree has been successfully 
      ! opened.                                                         
      ! Written: 2/24/99 by rjg                                         
      subroutine mds_char (node, signal, carg, ok) 
      implicit none 
                                  ! IN: the proper node name            
      character*(*) node 
                                  ! IN: the proper signal name          
      character*(*) signal 
                                  ! OUT: the data                       
      character*(*) carg 
                                  ! OUT: true means success; false means
      logical ok 
                                  !      there was a problem            
                                                                        
      ! Local variables                                                 
      integer MdsValue 
      integer status, descr 
                                 ! Function to check if status return   
      logical good_status 
                                 ! from MDS function calls is okay      
      character expression*200 
                                  ! Function to get length of character 
!      integer lencstr 
      integer len1,len2 
                                                                        
      ok = .false. 
                                                                        
      len1 = len_trim(node) 
      len2 = len_trim(signal) 
                                                                        
      write (expression,fmt='(A,A,A)') node(1:len1), signal(1:len2),    &
     &          CHAR(0)                                                 
                                                                        
         status = MdsValue(expression, descr(14,carg,0),0) 
         if (good_status(status)) then 
            ok = .true. 
         else 
            ok = .false. 
!            print *,' Bad status from MdsValue for data = ', status,   
!     +              ' for signal ', signal                             
         endif 
                                                                        
            ! mds_char                                                  
      END                                           
                                                                        
                                                                        
                                                                        
      ! Function to interpret status return from MDS function calls.    
      ! Return TRUE if status is odd; else return FALSE.                
      ! Written: 1/26/99 by rjg                                         
      logical function good_status (status) 
      implicit none 
      integer*4 status 
                                                                        
      integer num 
                                                                        
      num = status/2 
      num = status - 2*num 
                                                                        
      if (num .eq. 1) then 
         good_status = .true. 
      else 
         good_status = .false. 
      endif 
                                                                        
            ! function good_status                                      
      END                                           
                                                                        
                                                                        
!      ! Function to get length of a character string                    
!      ! Written: 1/28/99 by rjg                                         
!      integer function lencstr(str) 
!      implicit none 
!                           ! IN/OUT the string                          
!      character *(*) str 
!      integer len 
!                                                                        
!!      call str$trim(str,str,len) 
!      len = len_trim(str)
!      lencstr = len 
!                     ! Function lencstr                                 
!      END                                           
!
