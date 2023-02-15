!-----------------------------------------------------------------------
!     file hdf_api
!     hdf_api module
!      Very generic module meant for writing HDF5 files with particular
!      attributes. In case we want to convert over to C, it should make
!      this easier.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     code organization for hdf_api.
!-----------------------------------------------------------------------
!      0. check_dims 
!      1. h5accessMethod
!      2. vshdf5_fcinit
!      3. vshdf5_inith5vars
!      4. open_h5file
!      5. open_oldh5file
!      6. open_newh5file
!      7. close_h5file
!      8. open_group
!      9. make_group
!     10. close_group
!     11. test_group
!     12. get_nmember
!     13. obj_exists
!     14. attr_exists
!     15. make reference
!     16. make_mesh_group
!     17. make_time_group
!     18. make_vec_group
!     19. make_limits_group
!     20. write_attribute_ch_sc
!     21. write_attribute_ch_vec
!     22. write_attribute_int_sc
!     23. write_attribute_int_vec
!     24. write_attribute_intl_sc
!     25. write_attribute_int_vec
!     26. write_attribute_rl_sc
!     27. write_attribute_rl_vec
!     28. write_attribute_rls_sc
!     29. write_attribute_log_sc
!     30. write_attribute_log_vec
!     31. dump_h5in_attributes
!     32. dump_ch_sc
!     33. dump_ch_vec
!     34. dump_int
!     35. dump_int_1d
!     36. dump_int_2d
!     37. dump_h5_int_dum
!     38. dump_intl
!     39. dump_intl_1d
!     40. dump_h5_intl_dum
!     41. dump_rl
!     42. dump_rls
!     43. dump_rl_1d
!     44. dump_rl_2d
!     45. dump_rl_3d
!     46. dump_rl_4d
!     47. dump_rl_5d
!     48. dump_h5_rl_dum
!     49. dump_rls_1d
!     50. dump_rls_2d
!     51. dump_rls_3d
!     52. dump_rls_4d
!     53. dump_h5_rls_dum   
!     54. add_h5_int
!     55. add_h5_dbl
!     56. add_h5_int_1d
!     57. add_h5_1d
!     58. add_h5_2d
!     59. add_h5_3d
!     60. add_h5_4d
!     61. dset_dims
!     62. read_dims
!     63. read_ndims
!     64. read_ch_sc
!     65. read_ch_vec
!     66. read_int
!     67. read_int_1d
!     68. read_int_2d
!     69. read_intl
!     70. read_intl_1d
!     71. read_rl
!     72. read_rl_1d
!     73. read_rl_2d
!     74. read_rl_3d
!     75. read_rl_4d
!     76. read_5d
!     77. read_attribute_int_sc
!     78. read_attribute_int_vec
!     79. read_attribute_intl_sc
!     80. read_attribute_intl_vec
!     81. read_attribute_rl_sc
!     82. read_attribute_rl_vec
!     83. read_attribute_ch_sc
!     84. read_attribute_ch_vec
!     85. read_attribute_log_sc
!     86. read_attribute_log_sc
!     87. read_ch_sc_ex
!     88. read_ch_vec_ex
!     89. read_int_ex
!     90. read_int_1d_ex
!     91. read_int_2d_ex
!     92. read_intl_ex
!     93. read_intl_1d_ex
!     94. read_rl_ex
!     95. read_rl_1d_ex
!     96. read_rl_2d_ex
!     97. read_rl_3d_ex
!     98. read_rl_4d_ex
!     99. read_rl_5d_ex
!    100. read_log_ex
!    101. read_int_2d_sq
!    102. read_rl_2d_sq
!    102. read_rl_1d_sq
!------------------------------------------------------------------
!     module hdf_api
!-----------------------------------------------------------------------
      module hdf5_api
      use hdf5
      implicit none
      character(5), parameter, private :: h5fortranapiversion="1.0"
      integer, parameter, private :: i4=selected_int_kind(9)
      integer, parameter, private :: i8=selected_int_kind(18)
      integer, parameter, private :: r4=selected_real_kind(6,37)
      integer, parameter, private :: r8=selected_real_kind(13,307)
!-----------------------------------------------------------------------
!     Input parameters to control the attributes and how written out
!     The method used to determine whether to write the character
!     variables are written out is something like the following:
!        IF(TRIM(h5in%vsMD)/="") THEN
!           CALL write_attribute(dset_id,"vsMD",h5in%vsMD,errval)
!        ENDIF
!     This needs to be paid attention to if writing out several
!     variables in a row.
!
!     MESH DISCUSSION:
!     The key to making the vsSchema work is to associate fields with
!     their meshes.  The way to point the mesh it use the mesh member
!     the derived type below.  It is also used to define the mesh by
!     prepending the mesh type from visSchemawith "mesh-"; e.g.,
!      h5in%mesh="mesh-structured" defines a mesh and
!      h5in%mesh="/coreMesh"       points to where the mesh is defined.
!     For valid types of meshes, details on how the multi-domain
!     specification works, and centering issues, see the visSchema wiki:
!      https://ice.txcorp.com/trac/vizschema/wiki/
!
!-----------------------------------------------------------------------
!     IMPORTANT:::
!     It is very important to have codes use this API use the initvars
!     When adding variables to this derived type, make sure they are 
!      initializing correctly.
!-----------------------------------------------------------------------
      type hdf5inopts
        integer(hid_t) :: wrd_type       ! depricated, use typeconvert=T
        integer :: data_xfer_mode        ! always F (serial IO)
        logical :: noWrite               ! disable writes for parallel operation
        logical :: dotranspose           ! whether to tranpose 2D arrays
        logical :: verbose               ! whether to write verbose output
        logical :: debug                 ! write even more verbose output for debugging
        logical :: pio                   ! whether parallel i/o is used
        logical :: wrvstime              ! whether to write time to attribute
        logical :: typeconvert           ! whether to demote the type
        logical :: unitconvert           ! whether to write vsunitcnv
        integer(i4) :: comm     ! Communicator associated w/ open file
        integer(i4) :: info
        character(len=30) :: mesh             ! See above
        character(len=30) :: units            ! Units
        character(len=30) :: vsAxisLabels     ! Axis labels
        character(len=30) :: vsCentering      ! How to center variables on mesh
        character(len=30) :: vsMD             ! Multidomain variable
        character(len=30) :: vsTimeGroup      ! Time group label
        character(len=30) :: vsIndexOrder     ! Data ordering
        character(len=10000) :: vsLabels      ! Labels for (ic) (nqty)
        integer(i4) :: vsSpatialIndices(3)    ! Spatial indicies (<0 if 1 or 2D)
        real(r8) :: vsTime                    ! Time
        integer(i4) :: vsStep=0            ! Step # associated with time
        real(r8) :: vsUnitCnv=1.              ! conversion factor
        character(len=30), dimension(3) :: vsAxis   ! For rectilinear meshes
      end type
!-----------------------------------------------------------------------
!     Example of how to use h5err type after an fcapi call.
!     if(h5err%errBool) WRITE(*,fmt='(/,a,/)') h5err%errorMsg
!-----------------------------------------------------------------------
      type hdf5errortype
        logical :: errbool
        character(64) :: errormsg
      end type
!-----------------------------------------------------------------------
!     subprogram name interfaces, types defined as
!     ch     -> character
!     int    -> 4 byte int  h5t_std_i32le
!     intl   -> 8 byte int  h5t_std_i64le 
!     rls    -> 4 byte real h5t_ieee_f32le
!     rl     -> 8 byte real h5t_ieee_f64le
!-----------------------------------------------------------------------
      interface write_attribute
        module procedure write_attribute_ch_sc,write_attribute_ch_vec, &
                         write_attribute_int_sc,write_attribute_int_vec, &
                         write_attribute_intl_sc,write_attribute_intl_vec, &
                         write_attribute_rls_sc, &
                         write_attribute_rl_sc,write_attribute_rl_vec, &
                         write_attribute_log_sc,write_attribute_log_vec
      end interface
      interface dump_h5
        module procedure dump_ch_sc,dump_ch_vec, &
                         dump_int,dump_int_1d,dump_int_2d, &
                         dump_intl,dump_intl_1d, &
                         dump_rls,dump_rls_1d,dump_rls_2d,dump_rls_3d,dump_rls_4d, &
                         dump_rl,dump_rl_1d,dump_rl_2d,dump_rl_3d,dump_rl_4d,dump_rl_5d
      end interface
      ! This is like dump but does an append
      interface add_h5
        module procedure  add_h5_int,add_h5_dbl, &
                          add_h5_int_1d, add_h5_1d,add_h5_2d, &
                          add_h5_3d,add_h5_4d
      end interface
      interface read_h5
        module procedure read_ch_sc,read_ch_vec, &
                         read_int,read_int_1d,read_int_2d, &
                         read_intl,read_intl_1d, &
                         read_rl,read_rl_1d,read_rl_2d,read_rl_3d, &
                         read_rl_4d,read_rl_5d
      end interface
      interface read_attribute
        module procedure read_attribute_ch_sc,read_attribute_ch_vec, &
                         read_attribute_int_sc,read_attribute_int_vec, &
                         read_attribute_intl_sc,read_attribute_intl_vec, &
                         read_attribute_rl_sc,read_attribute_rl_vec, &
                         read_attribute_log_sc,read_attribute_log_vec
      end interface
      interface read_h5_ex
        module procedure read_ch_sc_ex,read_ch_vec_ex, &
                         read_int_ex,read_int_1d_ex,read_int_2d_ex, &
                         read_intl_ex,read_intl_1d_ex, &
                         read_rl_ex,read_rl_1d_ex,read_rl_2d_ex, &
                         read_rl_3d_ex,read_rl_4d_ex,read_rl_5d_ex, &
                         read_log_ex
      end interface
      interface read_h5_sq
        module procedure read_int_2d_sq,read_rl_2d_sq,read_rl_1d_sq
      end interface
    
      contains
!-----------------------------------------------------------------------
!     subprogram 0. check_dims 
!     Write mismatched dimension errors
!-----------------------------------------------------------------------
      subroutine check_dims(dims, fdims, errval)
      integer(hsize_t), dimension(:), intent(in) :: dims, fdims
      type(hdf5errortype), intent(inout) :: errval
        integer i
        do i = 1, size(dims)
          if (dims(i) /= fdims(i)) then
            write(*, *) "error: dims (", dims, ")"
            write(*, *) "  /=  fdims (", fdims, ")"
            write(*, *) "fdims = dims in the file (use h5ls)"
            write(*, *) "dims  = dims allocated for array to read"
            errval%errormsg = 'error: dims /= fdims'
            errval%errbool = .true.
            return
          endif
        enddo
        errval%errbool = .false.
      end subroutine check_dims
!-----------------------------------------------------------------------
!     subprogram 1. h5accessMethod
!     Return something the correct hdf5 access method given a more
!     memoral names
!-----------------------------------------------------------------------
      function h5accessmethod(access_method)
      integer(hid_t) :: h5accessmethod
      character(*), intent(in) :: access_method
      select case(access_method)
      case("overwr")
        h5accessmethod=h5f_acc_trunc_f                   ! overwrite file
      case("rdwr")
        h5accessmethod=h5f_acc_rdwr_f                    ! read-write
      case("rdonly")
        h5accessmethod=h5f_acc_rdonly_f                  ! read only
      end select
      return
      end function h5accessmethod
!-----------------------------------------------------------------------
!     subprogram 2. vshdf5_fcinit
!     Open fortran hdf5 and set open/close parameters.
!-----------------------------------------------------------------------
      subroutine vshdf5_fcinit()
      integer :: error
      ! integer :: err
      ! write(*, *) "vshdf5_fcinit: entered"
      call h5dont_atexit_f(error)
      ! write(*, *) "vshdf5_fcinit: h5dont_atexit_f returned."
      call h5open_f(error)
      ! write(*, *) "vshdf5_fcinit: h5open_f returned."
      ! write(*, *) "vshdf5_fcinit: leaving."
      return
      end subroutine vshdf5_fcinit
!-----------------------------------------------------------------------
!     subprogram 3. vshdf5_inith5vars
!     Initialize these variables to default values.
!-----------------------------------------------------------------------
      subroutine vshdf5_inith5vars(h5in, h5err)
      type(hdf5inopts), intent(inout) :: h5in
      type(hdf5errortype), intent(inout) :: h5err
!-----------------------------------------------------------------------
!     Defaults
!-----------------------------------------------------------------------
      ! According to xlf:
      ! (E) Null literal string is not permitted.  A single blank is assumed.
      ! So these should be single blanks to avoid warnings
      h5in%wrd_type=h5t_ieee_f64le
      h5in%data_xfer_mode=0
      h5in%noWrite=.false.
      h5in%vsCentering=" "
      h5in%doTranspose=.false.
      h5in%verbose=.false.
      h5in%debug=.false.
      h5in%pIO=.false.
      h5in%wrVsTime=.false.
      h5in%typeConvert=.false.
      h5in%unitConvert=.false.
      h5in%mesh =  " "
      h5in%vsAxisLabels = " " 
      h5in%units =  " "
      h5in%vsCentering =  " "
      h5in%vsMD =  " "
      ! Data-ordering options:
      ! [ix][iy][iz][ic] compMinorC [iz][iy][ix][ic] compMinorF
      ! [ic][ix][iy][iz] compMajorC [ic][iz][iy][ix] compMajorF
      h5in%vsIndexOrder = " "
      h5in%vsLabels = " "
      h5in%vsSpatialIndices=-1 ! Don't write
      h5err%errBool = .false.
      h5err%errorMsg =  " "
      return
      end subroutine vshdf5_inith5vars
!-----------------------------------------------------------------------
!     subprogram 4. open_h5file
!     Open file for writing and write file attributes
!     This is just a nice wrapper for open_newh5file and open_oldh5file where
!     we just specify the openmethod.  
!     Separate subroutines kept because for simplicity in debugging
!-----------------------------------------------------------------------
      subroutine open_h5file(openmethod,fname,fileid,fdesc,rootgid,h5in,h5err)
      character(*), intent(in) :: fname,fdesc,openmethod
      integer(hid_t), intent(out) :: fileid,rootgid
      type(hdf5errortype), intent(inout) :: h5err
      type(hdf5inopts), intent(inout) :: h5in
      integer,parameter :: fail=-1
    
      if (h5in%verbose) then
        write(*, *) " open_h5file: entered."
      endif
!-----------------------------------------------------------------------
!     Create and open the file
!     The if statements seem to work better to give the expected
!      behavior for read-write for already open file.
!-----------------------------------------------------------------------
      select case(openmethod)
      case('overwr')
        call open_newh5file(fname,fileid,fdesc,rootgid,h5in,h5err)
      case('append')
        call open_oldh5file(fname,fileid,rootgid,h5in,h5err)
      case default
        h5err%errorMsg = "open_h5file called with incorrect openmethod"
        h5err%errBool = .true.
        return
      end select
      return
      end subroutine open_h5file
!-----------------------------------------------------------------------
!     subprogram 5. open_oldh5file
!     Open file for writing and write file attributes
!     Create the group for the independent variables at this stage
!-----------------------------------------------------------------------
      subroutine open_oldh5file(fname,fileid,rootgid,h5in,h5err,rw)
      character(*), intent(in) :: fname
      integer(hid_t), intent(out) :: fileid,rootgid
      type(hdf5errortype), intent(inout) :: h5err
      type(hdf5inopts), intent(inout) :: h5in
      integer,parameter :: fail=-1
      integer :: error
      logical, intent(in), optional :: rw
    
      logical :: file_exists,read_write
      read_write=.false.
      if (present(rw)) read_write=rw
      if (h5in%verbose) then
        write(*, *) " open_oldh5file: entered."
      endif
!-----------------------------------------------------------------------
!     Create and open the file
!     The if statements seem to work better to give the expected
!      behavior for read-write for already open file.
!-----------------------------------------------------------------------
      inquire(file=trim(fname),exist=file_exists)
      if (.not. file_exists) then
        h5err%errormsg = 'error: file does not exist: '//fname
        if (h5in%verbose) &
          write(*, *) 'error: file does not exist: ', fname
        h5err%errbool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Open file
!-----------------------------------------------------------------------
      if (h5in%verbose) then
        write(*, *) "open_oldh5file: calling h5fopen_f."
      endif
        if (read_write) then
          call h5fopen_f(fname, H5F_ACC_RDWR_F, fileId, error)
        else
          call h5fopen_f(fname, H5F_ACC_RDONLY_F, fileId, error)
        endif
      if (h5in%verbose) then
        write(*, *) "open_oldh5file: h5fopen_f returned."
      endif
!-----------------------------------------------------------------------
!     Grab the root group id which is created by default
!-----------------------------------------------------------------------
      call h5gopen_f(fileId,"/",rootGid,error)
      if (error==FAIL) then
        h5err%errorMsg = 'ERROR: Error grabbing root ID: '//fname
        h5err%errBool = .true.
        return
      endif
      h5err%errBool = .false.
      return
      end subroutine open_oldh5file
!-----------------------------------------------------------------------
!     subprogram 6. open_newh5file
!     Open file for writing and write file attributes
!     Create the group for the independent variables at this stage
!-----------------------------------------------------------------------
      subroutine open_newh5file(fname,fileId,fdesc,rootGid,h5in,h5err)
      character(*), intent(in) :: fname,fdesc
      integer(HID_T), intent(out) :: fileId,rootGid
      type(hdf5InOpts), intent(inout) :: h5in
      type(hdf5ErrorType), intent(inout) :: h5err
      integer :: access_mode
      integer,parameter :: FAIL=-1
      integer :: error
      LOGICAL :: file_exists
!-----------------------------------------------------------------------
!     Open file
!     The if statements seem to work better to give the expected
!      behavior for read-write for already open file.
!-----------------------------------------------------------------------
      INQUIRE(FILE=TRIM(fname),EXIST=file_exists)
      if (.NOT. file_exists) then
        ! Always create with over-write to avoid errors
        access_mode=h5accessMethod("overwr")
        call h5fcreate_f(TRIM(fname),access_mode,fileId,error)
      else
        OPEN(UNIT=999,FILE=TRIM(fname),FORM='UNFORMATTED', &
             POSITION='REWIND',STATUS='REPLACE')
        CLOSE(UNIT=999)
        access_mode=h5accessMethod("overwr")
        call h5fcreate_f(TRIM(fname),access_mode,fileId,error)
      endif
      if (error==FAIL) then
        h5err%errorMsg = 'ERROR: Error opening file: '//fname
        h5err%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Grab the root group id which is created by default
!-----------------------------------------------------------------------
      call h5gopen_f(fileId,"/",rootGid,error)
      if (error==FAIL) then
        h5err%errorMsg = 'ERROR: Error grabbing root ID: '//fname
        h5err%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Put in things which describe the api
!-----------------------------------------------------------------------
!      call write_attribute(rootGid,"vsFcVERSION","1.0",h5in,h5err)
!-----------------------------------------------------------------------
!     Put in file description when creating file
!-----------------------------------------------------------------------
      call write_attribute(rootGid,"Description",fdesc,h5in,h5err)
!-----------------------------------------------------------------------
      h5err%errBool = .false.
      return
      end subroutine open_newh5file
!-----------------------------------------------------------------------
!     subprogram 7. close_h5file
!         Close the file associated with fileId.
!-----------------------------------------------------------------------
      subroutine close_h5file(fileId,root_id,h5err)
      integer(HID_T), intent(in)         :: fileId
      integer(HID_T), intent(in)         :: root_id
      type(hdf5ErrorType), intent(inout) :: h5err
    
      integer,parameter :: FAIL=-1
      integer :: error
      integer(i8) :: nobjs
!-----------------------------------------------------------------------
!     Close the root group and file
!-----------------------------------------------------------------------
      call h5gclose_f(root_id, error)
      call h5fclose_f(fileId, error)
      if (error==FAIL) then
        h5err%errorMsg = 'ERROR: Error in close_h5file'
        h5err%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      h5err%errBool = .false.
      return
      end subroutine close_h5file
!-----------------------------------------------------------------------
!     subprogram 8. open_group
!         Open a group in a safe way
!-----------------------------------------------------------------------
      subroutine open_group(inid,gname,gid,errval)
      character(*), intent(in) :: gname
      integer(HID_T), intent(in) :: inid
      integer(HID_T), intent(out) :: gid
      type(hdf5ErrorType), intent(inout) :: errval
      integer,parameter :: FAIL=-1
      integer :: error
!-------------------------------------------------------------------
!     Open group
!-----------------------------------------------------------------------
      call h5gopen_f(inid,gname,gid,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Error opening group: '//gname
        errval%errBool = .true.
        return
      endif
      errval%errBool = .false.
      return
      end subroutine open_group
!-----------------------------------------------------------------------
!     subprogram 9. make_group
!     Create a group in a safe way
!-----------------------------------------------------------------------
      subroutine make_group(inid,gname,gid,h5in,errval)
      character(*), intent(in) :: gname
      integer(HID_T), intent(in) :: inid
      integer(HID_T), intent(out) :: gid
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType) :: errval
      integer,parameter :: FAIL=-1
      integer :: error
      integer :: i
      character(1) :: alabel
      character(7) :: axisstr
!-----------------------------------------------------------------------
!     Create group
!-----------------------------------------------------------------------
      !call h5gopen_f(inid,gname,gid,error)
      ! If it failed, most likely it doesn't exist
      call h5gcreate_f(inid,gname,gid,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Error opening group: '//gname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
!      call write_attribute(gid,'CLASS','GROUP',h5in,errval)
!      call write_attribute(gid,'VERSION','1.0',h5in,errval)
      if(len_trim(h5in%mesh)>0) then
       if(h5in%mesh(1:5)=="mesh-".and..not.h5in%mesh(1:9)=="mesh-stru") then
         if(h5in%debug) WRITE(*,*) 'Writing vsType attributes'
         call write_attribute(gid,'vsType',"mesh",h5in,errval)
         if(h5in%debug) WRITE(*,*) 'Writing vsMesh attributes',h5in%mesh(6:)
         call write_attribute(gid,'vsKind',h5in%mesh(6:),h5in,errval)
         do i=1,3
           if (len_trim(h5in%vsAxis(i))>0) then
             write(alabel,fmt='(i1.1)') i-1
             axisstr="vsAxis"//alabel
             call write_attribute(gid,axisstr,h5in%vsAxis(i),h5in,errval)
           endif
         enddo
!SEK: Not sure 
!       else
!         if(h5in%debug) WRITE(*,*) 'Writing vsType attributes'
!         call write_attribute(dset_id,'vsType',"variable",h5in,errval)
!         if(h5in%debug) WRITE(*,*) 'Writing vsMesh attributes',h5in%mesh
!         call write_attribute(dset_id,'vsMesh',h5in%mesh,h5in,errval)
!         if(len_trim(h5in%vsCentering)>0) then
!           if(h5in%debug) WRITE(*,*) 'Writing vsCentering attributes'
!           call write_attribute(dset_id,'vsCentering',h5in%vsCentering, &
!                                  h5in,errval)
       endif
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine make_group
!-----------------------------------------------------------------------
!     subprogram 10. close_group
!     Close a group in a safe way
!-----------------------------------------------------------------------
      subroutine close_group(gname,inid,errval)
      character(*), intent(in) :: gname
      integer(HID_T), intent(in) :: inid
      type(hdf5ErrorType) :: errval
      integer,parameter :: FAIL=-1
      integer :: error
!-----------------------------------------------------------------------
!     Close group
!-----------------------------------------------------------------------
      call h5gclose_f(inid,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Error closing group: '//TRIM(gname)
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine close_group
!-----------------------------------------------------------------------
!     subprogram 11. test_group
!     See if group exists
!-----------------------------------------------------------------------
      subroutine test_group(inid,gname,group_exists,errval)
      integer(HID_T), intent(in) :: inid
      character(*), intent(in) :: gname
      LOGICAL, intent(out) :: group_exists
      type(hdf5ErrorType), intent(inout) :: errval
      integer,parameter :: FAIL=-1
      integer :: error
      integer(HID_T) :: gid
!-----------------------------------------------------------------------
!     First test that there is an object with this name
!-----------------------------------------------------------------------
      if (.not. obj_exists(inid,gname,errval)) then
        group_exists=.FALSE.
        errval%errBool = .false.
        return
      endif
!-----------------------------------------------------------------------
!     Determine whether group exists by trying to opening it and testing
!      error message
!-----------------------------------------------------------------------
      call h5gopen_f(inid,gname,gid,error)
      if (error==FAIL) then
        group_exists=.FALSE.
      else
        group_exists=.TRUE.
        call h5gclose_f(gid, error)
      endif
!-------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine test_group
!-----------------------------------------------------------------------
!     subprogram 12. get_nmember
!     Get the number of members of a group
!-----------------------------------------------------------------------
      subroutine get_nmembers(inid,gname,nmembers,errval)
      character(*), intent(in) :: gname
      integer(HID_T), intent(in) :: inid
      integer, intent(out) :: nmembers
      type(hdf5ErrorType) :: errval
      integer,parameter :: FAIL=-1
      integer :: error
!-----------------------------------------------------------------------
!     Determine whether group exists by trying to opening it and testing
!      error message
!-----------------------------------------------------------------------
      call h5gn_members_f(inid,gname,nmembers,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Error in get_nmembers for'//gname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine get_nmembers
!-----------------------------------------------------------------------
!     subprogram 13. obj_exists
!     Determine whether object oname at location inid exists
!-----------------------------------------------------------------------
      function obj_exists(inid,oname,errval)
      integer(HID_T), intent(in) :: inid
      character(*), intent(in) :: oname
      type(hdf5ErrorType) :: errval
      logical :: obj_exists
      integer,parameter :: FAIL=-1
      integer :: error
!-----------------------------------------------------------------------
!     Determine whether object exists
!-----------------------------------------------------------------------
      call h5eset_auto_f(0_i4,error) ! disable errors
      call h5oexists_by_name_f(inid,oname,obj_exists,error)
      call h5eset_auto_f(1_i4,error) ! enable errors
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Error in obj_exists for'//oname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      end function obj_exists
!-----------------------------------------------------------------------
!     subprogram 14. attr_exists
!     Determine whether attribute aname at location inid exists
!-----------------------------------------------------------------------
      function attr_exists(inid,oname,aname,errval)
      integer(HID_T), intent(in) :: inid
      character(*), intent(in) :: oname
      character(*), intent(in) :: aname
      type(hdf5ErrorType) :: errval
      logical :: attr_exists
      integer,parameter :: FAIL=-1
      integer :: error
!-----------------------------------------------------------------------
!     Determine whether object exists
!-----------------------------------------------------------------------
      call h5eset_auto_f(0_i4,error) ! disable errors
      call h5aexists_by_name_f(inid,oname,aname,attr_exists,error)
      call h5eset_auto_f(1_i4,error) ! enable errors
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Error in attr_exists for'//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      end function attr_exists
!-----------------------------------------------------------------------
!     subprogram 15. make reference
!     Simplify referencing of one object (source) to another (target)
!-----------------------------------------------------------------------
      subroutine make_reference(inid,tname,errval)
      integer(HID_T), intent(in) :: inid
      character*(*), intent(in) :: tname
      type(hdf5ErrorType) :: errval
      integer :: error
    
      integer(HSIZE_T), dimension(1) :: dimsr= (/4/)
      integer :: trank = 1
      integer(HID_T) :: type_id      ! Attribute Dataspace identifier
      integer(HID_T) :: tgt_id,tgt_sid       ! Attribute identifier
!-----------------------------------------------------------------------
!     See refobjexample.f90
!-----------------------------------------------------------------------
      !
      ! Create dataspace and dataset to store references to the objects
      !
      call h5screate_simple_f(trank, dimsr, tgt_sid, error)
      call h5dcreate_f(inid,tname,H5T_STD_REF_OBJ,tgt_sid,tgt_id,error)
      !
      ! Create a datatype and store in the file
      !
      call h5tcopy_f(h5t_ieee_f32le, type_id, error)
      call h5tcommit_f(inid, "MyType", type_id, error)
    
      errval%errBool = .false.
      return
      end subroutine make_reference
!-----------------------------------------------------------------------
!     subprogram 16. make_mesh_group
!     Defines a mesh group that points to other variables that define
!      the actual mesh
!-----------------------------------------------------------------------
      subroutine make_mesh_group(gInId,gridId,h5in,meshName,&
                  meshKind,axis0,axis1,axis2,transform,trName,errval)
      integer(HID_T), intent(in) :: gInId
      integer(HID_T), intent(inout) :: gridId
      type(hdf5InOpts), intent(inout) :: h5in
      character*(*), intent(in) :: meshname,axis0,axis1,axis2
      character*(*), intent(in) :: meshKind,transform,trName
      type(hdf5ErrorType), intent(inout) :: errval
      integer,parameter :: FAIL=-1
!-----------------------------------------------------------------------
!     Open the group
!-----------------------------------------------------------------------
      call make_group(gInId, meshName, gridId,h5in,errval)
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call write_attribute(gridId,'vsType',"mesh",h5in,errval)
      call write_attribute(gridId,'vsKind',meshKind,h5in,errval)
      call write_attribute(gridId,'vsAxis0',axis0,h5in,errval)
      call write_attribute(gridId,'vsAxis1',axis1,h5in,errval)
      call write_attribute(gridId,'vsAxis2',axis2,h5in,errval)
      h5in%vsAxisLabels=trim(axis0)//", "//trim(axis1)
      h5in%vsAxisLabels=h5in%vsAxisLabels//", "//trim(axis2)
      call write_attribute(gridId,'vsAxisLabels',h5in%vsAxisLabels,h5in,errval)
      if (len_trim(transform)>0) then
        call write_attribute(gridId,'vsTransform',transform,h5in,errval)
        call write_attribute(gridId,'vsTransformedMesh',trName,h5in,errval)
      endif
      if (len_trim(h5in%vsCentering)>0) then
        call write_attribute(gridId,'vsCentering',h5in%vsCentering,&
                             h5in,errval)
      endif
!-----------------------------------------------------------------------
!     vsMD: Multidomain cabilities
!-----------------------------------------------------------------------
      if (len_trim(h5in%vsMD)>0) then
        call write_attribute(gridId,"vsMD",h5in%vsMD,h5in,errval)
      endif
!-----------------------------------------------------------------------
      return
      end subroutine make_mesh_group
!-----------------------------------------------------------------------
!     subprogram 17. make_time_group
!     Make a group that contains the time data.  See:
!        https://ice.txcorp.com/trac/vizschema/wiki/OtherMetaData
!-----------------------------------------------------------------------
      subroutine make_time_group(gInId,h5in,h5err)
      integer(HID_T), intent(in) :: gInId
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: h5err
      integer,parameter :: FAIL=-1
      integer(HID_T) :: timeId
!-----------------------------------------------------------------------
!     Open the group
!-----------------------------------------------------------------------
      call make_group(gInId, h5in%vstimegroup, timeId,h5in,h5err)
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call write_attribute(timeId,'vsType',"time",h5in,h5err)
      write(*,*) h5in%vsStep
      call write_attribute(timeId,'vsStep',h5in%vsStep,h5in,h5err)
      call write_attribute(timeId,'vsTime',h5in%vsTime,h5in,h5err)
      if (len_trim(h5in%units)>0) then
        if(h5in%debug) WRITE(*,*) 'Writing time units',h5in%units
        call write_attribute(timeId,"units",h5in%units,h5in,h5err)
      endif
      call close_group(h5in%vstimegroup,timeId,h5err)
      return
      end subroutine make_time_group
!-----------------------------------------------------------------------
!     subprogram 18. make_vec_group
!     Make a group that defines a vector.  See:
!        https://ice.txcorp.com/trac/vizschema/wiki/OtherMetaData
!-----------------------------------------------------------------------
      subroutine make_vec_group(gInId,grName,veclabel,h5in,h5err)
      integer(HID_T), INTENT(IN) :: gInId
      character*(*), INTENT(IN) :: grName,veclabel
      type(hdf5InOpts), INTENT(IN) :: h5in
      type(hdf5ErrorType), INTENT(INOUT) :: h5err
      integer(HID_T) :: vecId
!-----------------------------------------------------------------------
!     Open the group
!-----------------------------------------------------------------------
      call make_group(gInId, grName, vecId,h5in,h5err)
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call write_attribute(vecId,'vsType',"vsVars",h5in,h5err)
      call write_attribute(vecId,grName,veclabel,h5in,h5err)
      call close_group(grName,vecId,h5err)
      return
      end subroutine make_vec_group
!-----------------------------------------------------------------------
!     subprogram 19. make_limits_group
!     Make a group that contains the visualization region data.  See:
!        https://ice.txcorp.com/trac/vizschema/wiki/OtherMetaData
!-----------------------------------------------------------------------
      SUBROUTINE make_limits_group(gInId,grName,vsKind,lowerBound,      &
        upperBound,h5in,h5err)
      integer(HID_T), intent(in) :: gInId
      character*(*), intent(in) :: grName,vsKind
      real(r8), dimension(:), intent(in) :: lowerBound, upperBound      
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: h5err
      integer(HID_T) :: limitId
!-----------------------------------------------------------------------
!     Open the group
!-----------------------------------------------------------------------
      call make_group(gInId, grName, limitId,h5in,h5err)
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call write_attribute(limitId,'vsType',"region",h5in,h5err)
      call write_attribute(limitId,'vsKind',vsKind,h5in,h5err)
      call write_attribute(limitId,'vsLowerBound',lowerBound,h5in,h5err)
      call write_attribute(limitId,'vsUpperBound',upperBound,h5in,h5err)
      call close_group(grName,limitId,h5err)
      return
      end subroutine make_limits_group
!-----------------------------------------------------------------------
!     subprogram 20. write_attribute_ch_sc
!     Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
      subroutine write_attribute_ch_sc(inid,aname,attribute,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character*(*), intent(in) :: aname,attribute
      type(hdf5InOpts), INTENT(IN) :: h5in
      type(hdf5ErrorType) :: errval
    
      integer(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      integer(HID_T) :: atype_id      ! Attribute Dataspace identifier
      integer(HID_T) :: attr_id       ! Attribute identifier
      integer(SIZE_T) :: attrlen    ! Length of the attribute string
      integer(HSIZE_T), dimension(1) :: data_dims
      integer,parameter :: FAIL=-1
      integer :: error
!-----------------------------------------------------------------------
!     If it is a null value then no need to write it out.
!-----------------------------------------------------------------------
      if (len_trim(attribute)==0) then
        errval%errBool = .false.
        return
      endif
!-----------------------------------------------------------------------
!     See attrexample.f90
!-----------------------------------------------------------------------
      data_dims(1) = 1
      attrlen=len_trim(attribute)
    
      ! Create the data space for the attribute.
      call h5screate_f(H5S_SCALAR_F, aspace_id, error)
    
      ! Create datatype for the attribute.
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attrlen, error)
    
      ! Create dataset attribute for the group
      call h5acreate_f(inid, aname, atype_id, aspace_id, attr_id, error)
    
      if (error==FAIL) then
        errval%errorMsg = 'Cannot create attribute '//aname//attribute
        errval%errBool = .true.
        return
      else
        ! Write the attribute data.
        if (.not.h5in%noWrite) &
          call h5awrite_f(attr_id, atype_id, trim(attribute), data_dims, error)
        ! Close the attribute.
        call h5aclose_f(attr_id, error)
      endif
    
      ! Close the dataspace.
      call h5sclose_f(aspace_id,error)
    
      errval%errBool = .false.
      return
      end subroutine write_attribute_ch_sc
!-----------------------------------------------------------------------
!     subprogram 21. write_attribute_ch_vec
!     Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
      subroutine write_attribute_ch_vec(inid,aname,attribute,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character*(*), intent(in) :: aname
      character*(*), dimension(:), intent(in) :: attribute
      type(hdf5InOpts), INTENT(IN) :: h5in
      type(hdf5ErrorType) :: errval
    
      integer :: arank=1
      integer(i4) :: i
      integer(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      integer(HID_T) :: atype_id      ! Attribute Dataspace identifier
      integer(HID_T) :: attr_id       ! Attribute identifier
      integer(HSIZE_T), dimension(1) :: adims  ! Attribute dimension
      integer(SIZE_T) :: attrlen    ! Length of the attribute string
      integer(HSIZE_T), dimension(1) :: data_dims
      integer,parameter :: FAIL=-1
      integer :: error
      integer(SIZE_T) :: ati
!-----------------------------------------------------------------------
!     If it is a null value then no need to write it out.
!-----------------------------------------------------------------------
      if (len_trim(attribute(1))==0) then
        errval%errBool = .false.
        return
      endif
!-----------------------------------------------------------------------
!     See attrexample.f90
!-----------------------------------------------------------------------
      data_dims(1) = size(attribute)
      adims(:) = (/data_dims(1)/)
      attrlen=0
      do i=1,data_dims(1)
        ati = len_trim(attribute(i))
        attrlen = max(attrlen,ati)
      enddo
    
      ! Create the data space for the attribute.
      call h5screate_simple_f(arank, adims, aspace_id, error)
    
      ! Create datatype for the attribute.
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attrlen, error)
    
      ! Create dataset attribute for the group
      call h5acreate_f(inid, aname, atype_id, aspace_id,attr_id, error)
    
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Can not create attribute '//aname
        errval%errBool = .true.
        return
      else
        ! Write the attribute data.
        if (.not.h5in%noWrite) &
          call h5awrite_f(attr_id, atype_id, attribute(1:attrlen), data_dims, error)
        ! Close the attribute.
        call h5aclose_f(attr_id, error)
      endif
    
      ! Close the dataspace.
      call h5sclose_f(aspace_id,error)
    
      errval%errBool = .false.
      return
      end subroutine write_attribute_ch_vec
!-----------------------------------------------------------------------
!     subprogram 22. write_attribute_int_sc
!     Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
      subroutine write_attribute_int_sc(inid,aname,attribute,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character*(*), intent(in) :: aname
      integer(i4), intent(in) :: attribute
      type(hdf5InOpts), INTENT(IN) :: h5in
      type(hdf5ErrorType) :: errval
    
      integer :: arank=1
      integer(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      integer(HID_T) :: atype_id      ! Attribute Dataspace identifier
      integer(HID_T) :: attr_id       ! Attribute identifier
      integer(HSIZE_T), dimension(1) :: adims = (/1/) ! Attribute dimension
      integer(HSIZE_T), dimension(1) :: data_dims
      integer,parameter :: FAIL=-1
      integer :: error
!-----------------------------------------------------------------------
!     See attrexample.f90
!-----------------------------------------------------------------------
      data_dims(1) = 1
    
      ! Create the data space for the attribute.
      call h5screate_simple_f(arank, adims, aspace_id, error)
    
      ! Create dataset attribute for the group
      call h5tcopy_f(h5t_std_i32le, atype_id, error)
      call h5acreate_f(inid,aname,atype_id,aspace_id,attr_id,error)
    
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Can not create attribute '//aname
        errval%errBool = .true.
        return
      else
        ! Write the attribute data.
        if (.not.h5in%noWrite) &
          call h5awrite_f(attr_id, atype_id, attribute, data_dims, error)
        ! Close the attribute.
        call h5aclose_f(attr_id, error)
      endif
    
      ! Close the dataspace.
      call h5sclose_f(aspace_id,error)
    
      errval%errBool = .false.
      return
      end subroutine write_attribute_int_sc
!-----------------------------------------------------------------------
!     subprogram 23. write_attribute_int_vec
!     Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
      subroutine write_attribute_int_vec(inid,aname,attribute,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character*(*), intent(in) :: aname
      integer(i4), dimension(:), intent(in) :: attribute
      type(hdf5InOpts), INTENT(IN) :: h5in
      type(hdf5ErrorType) :: errval
    
      integer :: arank=1
      integer(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      integer(HID_T) :: atype_id      ! Attribute Dataspace identifier
      integer(HID_T) :: attr_id       ! Attribute identifier
      integer(HSIZE_T), dimension(1) :: adims  ! Attribute dimension
      integer(HSIZE_T), dimension(1) :: data_dims
      integer,parameter :: FAIL=-1
      integer :: error
!-----------------------------------------------------------------------
!     See attrexample.f90
!-----------------------------------------------------------------------
      data_dims(1) = size(attribute)
      adims(:) = (/data_dims(1)/)
    
      ! Create the data space for the attribute.
      call h5screate_simple_f(arank, adims, aspace_id, error)
    
      ! Create dataset attribute for the group
      call h5tcopy_f(h5t_std_i32le, atype_id, error)
      call h5acreate_f(inid,aname,atype_id,aspace_id,attr_id,error)
    
    
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Can not create attribute '//aname
        errval%errBool = .true.
        return
      else
        ! Write the attribute data.
        if (.not.h5in%noWrite) &
          call h5awrite_f(attr_id, atype_id, attribute, data_dims, error)
        ! Close the attribute.
        call h5aclose_f(attr_id, error)
      endif
    
      ! Close the dataspace.
      call h5sclose_f(aspace_id,error)
    
      errval%errBool = .false.
      return
      end subroutine write_attribute_int_vec
!-----------------------------------------------------------------------
!     subprogram 24. write_attribute_intl_sc
!     Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
      subroutine write_attribute_intl_sc(inid,aname,attribute,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character*(*), intent(in) :: aname
      integer(i8), intent(in) :: attribute
      type(hdf5InOpts), INTENT(IN) :: h5in
      type(hdf5ErrorType) :: errval
    
      integer :: arank=1
      integer(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      integer(HID_T) :: atype_id      ! Attribute Dataspace identifier
      integer(HID_T) :: attr_id       ! Attribute identifier
      integer(HSIZE_T), dimension(1) :: adims = (/1/) ! Attribute dimension
      integer(HSIZE_T), dimension(1) :: data_dims
      integer,parameter :: FAIL=-1
      integer :: error
!-----------------------------------------------------------------------
!     See attrexample.f90
!-----------------------------------------------------------------------
      data_dims(1) = 1
    
      ! Create the data space for the attribute.
      call h5screate_simple_f(arank, adims, aspace_id, error)
    
      ! Create dataset attribute for the group
      call h5tcopy_f(h5t_std_i64le, atype_id, error)
      call h5acreate_f(inid,aname,atype_id,aspace_id,attr_id,error)
    
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Can not create attribute '//aname
        errval%errBool = .true.
        return
      else
        ! Write the attribute data.
        if (.not.h5in%noWrite) &
          call h5awrite_f(attr_id, atype_id, int(attribute,i4), data_dims, error)
        ! Close the attribute.
        call h5aclose_f(attr_id, error)
      endif
    
      ! Close the dataspace.
      call h5sclose_f(aspace_id,error)
    
      errval%errBool = .false.
      return
      end subroutine write_attribute_intl_sc
!-----------------------------------------------------------------------
!     subprogram 25. write_attribute_int_vec
!     Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
      subroutine write_attribute_intl_vec(inid,aname,attribute,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character*(*), intent(in) :: aname
      integer(i8), dimension(:), intent(in) :: attribute
      type(hdf5InOpts), INTENT(IN) :: h5in
      type(hdf5ErrorType) :: errval
    
      integer :: arank=1
      integer(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      integer(HID_T) :: atype_id      ! Attribute Dataspace identifier
      integer(HID_T) :: attr_id       ! Attribute identifier
      integer(HSIZE_T), dimension(1) :: adims  ! Attribute dimension
      integer(HSIZE_T), dimension(1) :: data_dims
      integer,parameter :: FAIL=-1
      integer :: error
!-----------------------------------------------------------------------
!     See attrexample.f90
!-----------------------------------------------------------------------
      data_dims(1) = size(attribute)
      adims(:) = (/data_dims(1)/)
    
      ! Create the data space for the attribute.
      call h5screate_simple_f(arank, adims, aspace_id, error)
    
      ! Create dataset attribute for the group
      call h5tcopy_f(h5t_std_i64le, atype_id, error)
      call h5acreate_f(inid,aname,atype_id,aspace_id,attr_id,error)
    
    
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Can not create attribute '//aname
        errval%errBool = .true.
        return
      else
        ! Write the attribute data.
        if (.not.h5in%noWrite) &
          call h5awrite_f(attr_id, atype_id, int(attribute,i4), data_dims, error)
        ! Close the attribute.
        call h5aclose_f(attr_id, error)
      endif
    
      ! Close the dataspace.
      call h5sclose_f(aspace_id,error)
    
      errval%errBool = .false.
      return
      end subroutine write_attribute_intl_vec
!-----------------------------------------------------------------------
!     subprogram 26. write_attribute_rl_sc
!-----------------------------------------------------------------------
      subroutine write_attribute_rl_sc(inid,aname,attribute,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character*(*), intent(in) :: aname
      real(r8), intent(in) :: attribute
      type(hdf5InOpts), INTENT(IN) :: h5in
      type(hdf5ErrorType), intent(out) :: errval
    
      integer(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      integer(HID_T) :: atype_id      ! Attribute Dataspace identifier
      integer(HID_T) :: attr_id       ! Attribute identifier
      integer(HSIZE_T), dimension(1) :: data_dims
      integer,parameter :: FAIL=-1
      integer :: error
!-----------------------------------------------------------------------
!     See attrexample.f90
!-----------------------------------------------------------------------
      data_dims(1) = 1
    
      ! Create the data space for the attribute.
      call h5screate_f(H5S_SCALAR_F, aspace_id, error)
    
      ! Create dataset attribute for the group
    
      call h5tcopy_f(h5t_ieee_f64le, atype_id, error)
      call h5acreate_f(inid,aname,atype_id,aspace_id,attr_id,error)
    
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Can not create attribute '//aname
        errval%errBool = .true.
        return
      else
        ! Write the attribute data.
!-PRE        if (h5in%typeConvert) then
!-PRE          call h5awrite_f(attr_id, atype_id, real(attribute,r4), data_dims, error)
!-PRE        else
        if (.not.h5in%noWrite) &
          call h5awrite_f(attr_id, atype_id, attribute, data_dims, error)
!-PRE        endif
        ! Close the attribute.
        call h5aclose_f(attr_id, error)
      endif
    
      ! Close the dataspace.
      call h5sclose_f(aspace_id,error)
    
      errval%errBool = .false.
      return
      end subroutine write_attribute_rl_sc
!-----------------------------------------------------------------------
!     subprogram 27. write_attribute_rl_vec
!     Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
      subroutine write_attribute_rl_vec(inid,aname,attribute,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character*(*), intent(in) :: aname
      real(r8), dimension(:), intent(in) :: attribute
      type(hdf5InOpts), INTENT(IN) :: h5in
      type(hdf5ErrorType) :: errval
    
      integer :: arank=1
      integer(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      integer(HID_T) :: atype_id      ! Attribute Dataspace identifier
      integer(HID_T) :: attr_id       ! Attribute identifier
      integer(HSIZE_T), dimension(1) :: adims  ! Attribute dimension
      integer(HSIZE_T), dimension(1) :: data_dims
      integer,parameter :: FAIL=-1
      integer :: error
!-------------------------------------------------------------------
!     See attrexample.f90
!-------------------------------------------------------------------
      data_dims(1) = size(attribute)
      adims(:) = (/data_dims(1)/)
    
      ! Create the data space for the attribute.
      call h5screate_simple_f(arank, adims, aspace_id, error)
    
      ! Create dataset attribute for the group
      call h5tcopy_f(h5t_ieee_f64le, atype_id, error)
      call h5acreate_f(inid,aname,atype_id,aspace_id,attr_id,error)
    
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Can not create attribute '//aname
        errval%errBool = .true.
        return
      else
        ! Write the attribute data.
!-PRE        if(h5in%typeConvert) then
!-PRE          call h5awrite_f(attr_id, atype_id, real(attribute,r4), data_dims, error)
!-PRE        else
        if (.not.h5in%noWrite) &
          call h5awrite_f(attr_id, atype_id, attribute, data_dims, error)
!-PRE        endif
        ! Close the attribute.
        call h5aclose_f(attr_id, error)
      endif
    
      ! Close the dataspace.
      CALL h5sclose_f(aspace_id,error)
    
      errval%errBool = .false.
      return
      end subroutine write_attribute_rl_vec
!-----------------------------------------------------------------------
!     subprogram 28. write_attribute_rls_sc
!-----------------------------------------------------------------------
      subroutine write_attribute_rls_sc(inid,aname,attribute,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character*(*), intent(in) :: aname
      real(r4), intent(in) :: attribute
      type(hdf5InOpts), INTENT(IN) :: h5in
      type(hdf5ErrorType), intent(out) :: errval
    
      integer(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      integer(HID_T) :: atype_id      ! Attribute Dataspace identifier
      integer(HID_T) :: attr_id       ! Attribute identifier
      integer(HSIZE_T), dimension(1) :: data_dims
      integer,parameter :: FAIL=-1
      integer :: error
!-----------------------------------------------------------------------
!     See attrexample.f90
!-----------------------------------------------------------------------
      data_dims(1) = 1
    
      ! Create the data space for the attribute.
      call h5screate_f(H5S_SCALAR_F, aspace_id, error)
    
      ! Create dataset attribute for the group
    
      call h5tcopy_f(h5t_ieee_f32le, atype_id, error)
      call h5acreate_f(inid,aname,atype_id,aspace_id,attr_id,error)
    
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Can not create attribute '//aname
        errval%errBool = .true.
        return
      else
        ! Write the attribute data.
        if (.not.h5in%noWrite) &
          call h5awrite_f(attr_id, atype_id, attribute, data_dims, error)
        ! Close the attribute.
        call h5aclose_f(attr_id, error)
      endif
    
      ! Close the dataspace.
      call h5sclose_f(aspace_id,error)
    
      errval%errBool = .false.
      return
      end subroutine write_attribute_rls_sc
!-----------------------------------------------------------------------
!     subprogram 29. write_attribute_log_sc
!-----------------------------------------------------------------------
      subroutine write_attribute_log_sc(inid,aname,attribute,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character*(*), intent(in) :: aname
      logical, intent(in) :: attribute
      type(hdf5InOpts), INTENT(IN) :: h5in
      type(hdf5ErrorType) :: errval
    
      integer :: arank=1
      integer(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      integer(HID_T) :: atype_id      ! Attribute Dataspace identifier
      integer(HID_T) :: attr_id       ! Attribute identifier
      integer(HSIZE_T), dimension(1) :: adims = (/1/) ! Attribute dimension
      integer(HSIZE_T), dimension(1) :: data_dims
      integer(i4) :: iattribute
      integer,parameter :: FAIL=-1
      integer :: error
!-----------------------------------------------------------------------
!     See attrexample.f90
!-----------------------------------------------------------------------
      data_dims(1) = 1
    
      if (attribute) then
        iattribute=1
      else
        iattribute=0
      endif
    
      ! Create the data space for the attribute.
      call h5screate_simple_f(arank, adims, aspace_id, error)
    
      ! Create dataset attribute for the group
      call h5tcopy_f(h5t_std_i32le, atype_id, error)
      call h5acreate_f(inid,aname,atype_id,aspace_id,attr_id,error)
    
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Can not create attribute '//aname
        errval%errBool = .true.
        return
      else
        ! Write the attribute data.
        if (.not.h5in%noWrite) &
          call h5awrite_f(attr_id, atype_id, iattribute, data_dims, error)
        ! Close the attribute.
        call h5aclose_f(attr_id, error)
      endif
    
      ! Close the dataspace.
      call h5sclose_f(aspace_id,error)
    
      errval%errBool = .false.
      return
      end subroutine write_attribute_log_sc
!-----------------------------------------------------------------------
!     subprogram 30. write_attribute_log_vec
!     Create a group for the independent vars (aka dimensions, scales)
!-----------------------------------------------------------------------
      subroutine write_attribute_log_vec(inid,aname,attribute,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character*(*), intent(in) :: aname
      logical, dimension(:), intent(in) :: attribute
      type(hdf5InOpts), INTENT(IN) :: h5in
      type(hdf5ErrorType) :: errval
    
      integer :: arank=1
      integer(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      integer(HID_T) :: atype_id      ! Attribute Dataspace identifier
      integer(HID_T) :: attr_id       ! Attribute identifier
      integer(HSIZE_T), dimension(1) :: adims  ! Attribute dimension
      integer(HSIZE_T), dimension(1) :: data_dims
      integer(i4), allocatable, dimension(:) :: iattribute
      integer,parameter :: FAIL=-1
      integer :: error,ii
!-----------------------------------------------------------------------
!     See attrexample.f90
!-----------------------------------------------------------------------
      data_dims(1) = size(attribute)
      allocate(iattribute(data_dims(1)))
      do ii=1,data_dims(1)
        if (attribute(ii)) then
          iattribute(ii)=1
        else
          iattribute(ii)=0
        endif
      enddo
    
      adims(:) = (/data_dims(1)/)
    
      ! Create the data space for the attribute.
      call h5screate_simple_f(arank, adims, aspace_id, error)
    
      ! Create dataset attribute for the group
      call h5tcopy_f(h5t_std_i32le, atype_id, error)
      call h5acreate_f(inid,aname,atype_id,aspace_id,attr_id,error)
    
    
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Can not create attribute '//aname
        errval%errBool = .true.
        return
      else
        ! Write the attribute data.
        if (.not.h5in%noWrite) &
          call h5awrite_f(attr_id, atype_id, iattribute, data_dims, error)
        ! Close the attribute.
        call h5aclose_f(attr_id, error)
      endif
      deallocate(iattribute)
    
      ! Close the dataspace.
      call h5sclose_f(aspace_id,error)
    
      errval%errBool = .false.
      return
      end subroutine write_attribute_log_vec
!-----------------------------------------------------------------------
!     subprogram 31. dump_h5in_attributes
!     Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
      subroutine dump_h5in_attributes(dset_id,h5in,h5err)
      integer(HID_T), intent(in) :: dset_id
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: h5err
      integer,parameter :: FAIL=-1
      integer(i4) :: nspat
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      if (len_trim(h5in%mesh)>0) then
       if (h5in%mesh(1:5)=="mesh-") then
         if (h5in%debug) WRITE(*,*) 'Writing vsType attribute'
         call write_attribute(dset_id,'vsType',"mesh",h5in,h5err)
         if (h5in%debug) WRITE(*,*) 'Writing vsMesh attribute ',h5in%mesh(6:)
         call write_attribute(dset_id,'vsKind',h5in%mesh(6:),h5in,h5err)
         if (len_trim(h5in%vsAxisLabels)>0) then
           if (h5in%debug) WRITE(*,*) 'Writing vsAxisLabels'
           call write_attribute(dset_id,'vsAxisLabels',&
                                  h5in%vsAxisLabels,h5in,h5err)
         endif
       elseif (h5in%mesh(1:16)=="variableWithMesh") then
         if (h5in%debug) WRITE(*,*) 'Writing vsType attribute'
         call write_attribute(dset_id,'vsType',"variableWithMesh",h5in,h5err)
         if (h5in%vsSpatialIndices(1)>=0) then
           nspat=1
           if (h5in%vsSpatialIndices(2)>=0) nspat=2
           if (h5in%vsSpatialIndices(3)>=0) nspat=3
           if (h5in%debug) &
             WRITE(*,*) 'Writing vsSpatialIndices attribute'
           call write_attribute(dset_id,'vsSpatialIndices', &
                                h5in%vsSpatialIndices(1:nspat),h5in,h5err)
         endif
         if (len_trim(h5in%vsLabels)>0) then
           if (h5in%debug) WRITE(*,*) 'Writing vsLabels attribute'
           call write_attribute(dset_id,'vsLabels',h5in%vsLabels,h5in,h5err)
         endif
         if (len_trim(h5in%vstimegroup)>0 .and. h5in%wrvstime) then
           if (h5in%debug) WRITE(*,*) 'Writing vsTimeGroup attribute'
           call write_attribute(dset_id,'vsTimeGroup',h5in%vstimegroup,h5in,h5err)
         endif
       else
         if (h5in%debug) WRITE(*,*) 'Writing vsType attribute'
         call write_attribute(dset_id,'vsType',"variable",h5in,h5err)
         if (h5in%debug) WRITE(*,*) 'Writing vsMesh attribute ',h5in%mesh
         call write_attribute(dset_id,'vsMesh',h5in%mesh,h5in,h5err)
         if (len_trim(h5in%vsCentering)>0) then
           if (h5in%debug) WRITE(*,*) 'Writing vsCentering attribute'
           call write_attribute(dset_id,'vsCentering',h5in%vsCentering, &
                                  h5in,h5err)
         endif
         if (len_trim(h5in%vsLabels)>0) then
           if (h5in%debug) WRITE(*,*) 'Writing vsLabels attribute'
           call write_attribute(dset_id,'vsLabels',h5in%vsLabels,h5in,h5err)
         endif
         if (len_trim(h5in%vstimegroup)>0 .and. h5in%wrvstime) then
           if (h5in%debug) WRITE(*,*) 'Writing vsTimeGroup attribute'
           call write_attribute(dset_id,'vsTimeGroup',h5in%vstimegroup,h5in,h5err)
         endif
       endif
      endif
!-----------------------------------------------------------------------
!     vsIndexOrder: set ordering (defined in initialization)
!-----------------------------------------------------------------------
      if (len_trim(h5in%vsIndexOrder)>0) then
        if (h5in%debug) WRITE(*,*) 'Writing vsIndexOrder attribute'
        call write_attribute(dset_id,'vsIndexOrder',h5in%vsIndexOrder,h5in,h5err)
      endif
!-----------------------------------------------------------------------
!     vsMD: Multidomain cabilities
!-----------------------------------------------------------------------
      if (len_trim(h5in%vsMD)>0) then
        if (h5in%debug) WRITE(*,*) 'Writing vsMD attribute ',h5in%vsMD
        call write_attribute(dset_id,"vsMD",h5in%vsMD,h5in,h5err)
      endif
!-----------------------------------------------------------------------
!     Label the units
!-----------------------------------------------------------------------
      if (len_trim(h5in%units)>0) then
        if (h5in%debug) WRITE(*,*) 'Writing units attribute ',h5in%units
        call write_attribute(dset_id,"units",h5in%units,h5in,h5err)
      endif
!-----------------------------------------------------------------------
!     If we have the ability to write the conversion factors than do so
!-----------------------------------------------------------------------
      if (h5in%unitConvert) then
        call write_attribute(dset_id,'vsUnitConvert',h5in%vsUnitCnv,h5in,h5err)
      endif
!-----------------------------------------------------------------------
      if (h5in%debug) WRITE(*,*) "Returning from h5in_attributes"
      return
      end subroutine dump_h5in_attributes
!-----------------------------------------------------------------------
!     subprogram 32. dump_ch_sc
!     Write an hdf5 string
      ! H5: needs to be tested!
!-----------------------------------------------------------------------
      subroutine dump_ch_sc(inid,aname,dataset,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character*(*), intent(in) :: aname,dataset
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: errval
      integer,parameter :: FAIL=-1
      integer :: error
    
      integer(HID_T) :: dspace_id     ! Dataspace identifier
      integer(HID_T) :: dtype_id      ! Datatype identifier
      integer(HID_T) :: dset_id       ! Dataset identifier
      integer(SIZE_T) :: dsetlen    ! Length of the dataset
      integer(HSIZE_T), dimension(1) :: dims
!-----------------------------------------------------------------------
!     If it is a null value then no need to write it out.
!-----------------------------------------------------------------------
      if (len_trim(dataset)==0) then
        errval%errBool = .false.
        return
      endif
!-------------------------------------------------------------------
      if(h5in%verbose) WRITE(*,*) 'Writing ', aname
      dims(1) = 1
      dsetlen = len_trim(dataset)
!-------------------------------------------------------------------
!     Create the data space.
!-------------------------------------------------------------------
      call h5screate_f(H5S_SCALAR_F,dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-------------------------------------------------------------------
!     Create the datatype.
!-------------------------------------------------------------------
      call h5tcopy_f(h5t_native_character,dtype_id,error)
      call h5tset_size_f(dtype_id,dsetlen,error)
!-------------------------------------------------------------------
!     Create the data set.
!-------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dcreate_f(inid,aname,dtype_id,dspace_id,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
      call h5dwrite_f(dset_id,dtype_id,trim(dataset),dims,error)
      ! H5: clean-up
!      call h5dwrite_f(dset_id,h5t_ieee_f64le,dataset,dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Data set write failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine dump_ch_sc
!-----------------------------------------------------------------------
!     subprogram 33. dump_ch_vec
!     Write an hdf5 character array
      ! H5: needs to be tested!
!-----------------------------------------------------------------------
      subroutine dump_ch_vec(inid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character*(*), intent(in) :: aname
      character(*), dimension(:), intent(in) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: errval
      integer,parameter :: FAIL=-1
      integer :: error
    
      integer :: arank=1
      integer(i4) :: i
      integer(HID_T) :: dspace_id     ! Dataspace identifier
      integer(HID_T) :: dtype_id      ! Datatype identifier
      integer(HID_T) :: dset_id       ! Dataset identifier
      integer(HSIZE_T), dimension(1) :: adims  ! Array dimension
      integer(SIZE_T) :: arrlen    ! Length of the character array
      integer(HSIZE_T), dimension(1) :: dims
      integer(SIZE_T) :: ari
!-----------------------------------------------------------------------
!     If it is a null value then no need to write it out.
!-----------------------------------------------------------------------
      if (len_trim(array(1))==0) then
        errval%errBool = .false.
        return
      endif
!-------------------------------------------------------------------
      if(h5in%verbose) WRITE(*,*) 'Writing ', aname
      dims(1) = size(array)
      adims(:) = (/dims(1)/)
      arrlen=0
      do i=1,adims(1)
        ari = len_trim(array(i))
        arrlen = max(arrlen,ari)
      enddo
!-------------------------------------------------------------------
!     Create the data space.
!-------------------------------------------------------------------
      call h5screate_simple_f(arank,adims,dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-------------------------------------------------------------------
!     Create the datatype.
!-------------------------------------------------------------------
      call h5tcopy_f(h5t_native_character,dtype_id,error)
      call h5tset_size_f(dtype_id,arrlen,error)
!-------------------------------------------------------------------
!     Create the data set.
!-------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dcreate_f(inid,aname,dtype_id,dspace_id,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
      call h5dwrite_f(dset_id,dtype_id,array(1:arrlen),dims,error)
      ! H5: clean-up
!      call h5dwrite_f(dset_id,h5t_ieee_f64le,array,dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Data set write failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine dump_ch_vec
!-----------------------------------------------------------------------
!     subprogram 34. dump_int
!     Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
      subroutine dump_int(inid,aname,val,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character*(*), intent(in) :: aname
      integer(i4), intent(in) :: val
      type(hdf5ErrorType) :: errval
      type(hdf5InOpts), intent(in) :: h5in
      integer,parameter :: FAIL=-1
      integer :: error
      integer(HID_T) :: dspace_id
      integer(HID_T) :: dset_id
      integer(HSIZE_T), dimension(1) :: dims=0
!-----------------------------------------------------------------------
      if(h5in%verbose) WRITE(*,*) 'Writing ', aname
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      call h5screate_f(H5S_SCALAR_F, dspace_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Create the data set.
!-----------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dcreate_f(inid,aname, &
                       h5t_std_i32le,dspace_id,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
      call h5dwrite_f(dset_id,h5t_std_i32le,val,dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Data set write failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine dump_int
!-----------------------------------------------------------------------
!     subprogram 35. dump_int_1d
!     Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
      subroutine dump_int_1d(inid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character*(*), intent(in) :: aname
      integer(i4), dimension(:), intent(in) :: array
      type(hdf5ErrorType) :: errval
      type(hdf5InOpts), intent(in) :: h5in
      integer,parameter :: FAIL=-1
      integer :: error
      integer :: rank
      integer(HID_T) :: dspace_id, dset_id
      integer(HSIZE_T), dimension(1) :: dims
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      if(h5in%verbose) WRITE(*,*) 'Writing ', aname
      rank = 1;            dims(:) = (/size(array,1)/)
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      call h5screate_simple_f(rank,dims,dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Create the data set.
!-----------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dcreate_f(inid,aname,h5t_std_i32le,dspace_id,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
      if (h5in%typeConvert) then
        call h5dwrite_f(dset_id,h5t_std_i32le,real(array,r4),dims, &
                        error)
      else
        call h5dwrite_f(dset_id,h5t_std_i32le,array,dims,error)
      endif
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Data set write failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine dump_int_1d
!-----------------------------------------------------------------------
!     subprogram 36. dump_int_2d
!     Create a "simple dataset" and write it out.
!-----------------------------------------------------------------------
      subroutine dump_int_2d(inid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character(*), intent(in) :: aname
      integer(i4), dimension(:,:), intent(in) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType) :: errval
      integer,parameter :: FAIL=-1
      integer :: error
    
      integer(HID_T) :: dspace_id, dset_id
      integer :: rank
      integer(HSIZE_T) :: dims(2)
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) 'Writing ', aname
      rank = 2
      if (h5in%doTranspose) then
        dims(2) = size(array,1);  dims(1) = size(array,2)
      else
        dims(1) = size(array,1);  dims(2) = size(array,2)
      endif
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      call h5screate_simple_f(rank,dims,dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Create the data set.
!-----------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dcreate_f(inid,aname,h5t_std_i32le,dspace_id,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
      if (h5in%doTranspose) then
        call h5dwrite_f(dset_id,h5t_std_i32le,transpose(array),dims, &
                        error)
      else
        call h5dwrite_f(dset_id,h5t_std_i32le,array,dims,error)
      endif
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Writing data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine dump_int_2d
!-----------------------------------------------------------------------
!     subprogram 37. dump_h5_int_dum
!-----------------------------------------------------------------------
      subroutine dump_h5_int_dum(inid,aname,dims,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character(*), intent(in) :: aname
      integer(i4), dimension(:), intent(in) :: dims
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType) :: errval
    
      integer,parameter :: FAIL=-1
      integer :: error
      integer(i4) :: ii
      integer(HID_T) :: dspace_id, dset_id
      integer :: rank
      integer(HSIZE_T), allocatable :: tdims(:)
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      rank=size(dims)
      allocate(tdims(rank))
      if (h5in%doTranspose) then
        do ii=1,size(dims)
          tdims(size(dims)-ii+1)=dims(ii)
        enddo
      else
        tdims=dims
      endif
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      call h5screate_simple_f(rank,tdims,dspace_id,error)
      deallocate(tdims)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Create the data set.
!-----------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dcreate_f(inid,aname,h5t_std_i32le,dspace_id, &
                       dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine dump_h5_int_dum
!-----------------------------------------------------------------------
!     subprogram 38. dump_intl
!     Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
      subroutine dump_intl(inid,aname,val,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character*(*), intent(in) :: aname
      integer(i8), intent(in) :: val
      type(hdf5ErrorType) :: errval
      type(hdf5InOpts), intent(in) :: h5in
      integer,parameter :: FAIL=-1
      integer :: error
      integer(HID_T) :: dspace_id
      integer(HID_T) :: dset_id
      integer(HSIZE_T), dimension(1) :: dims=0
!-----------------------------------------------------------------------
      if(h5in%verbose) WRITE(*,*) 'Writing ', aname
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      call h5screate_f(H5S_SCALAR_F, dspace_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Create the data set.
!-----------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dcreate_f(inid,aname, &
                       h5t_std_i64le,dspace_id,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
      call h5dwrite_f(dset_id,h5t_std_i32le,int(val,i4),dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Data set write failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine dump_intl
!-----------------------------------------------------------------------
!     subprogram 39. dump_intl_1d
!     Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
      subroutine dump_intl_1d(inid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character*(*), intent(in) :: aname
      integer(i8), dimension(:), intent(in) :: array
      integer(i4), dimension(:), allocatable :: intarray
      type(hdf5ErrorType) :: errval
      type(hdf5InOpts), intent(in) :: h5in
      integer,parameter :: FAIL=-1
      integer :: error
      integer :: rank
      integer(HID_T) :: dspace_id, dset_id
      integer(HSIZE_T), dimension(1) :: dims
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      if(h5in%verbose) WRITE(*,*) 'Writing ', aname
      rank = 1;            dims(:) = (/size(array,1)/)
      allocate(intarray(dims(1)))
      intarray=array
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      call h5screate_simple_f(rank,dims,dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Create the data set.
!-----------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dcreate_f(inid,aname,h5t_std_i64le,dspace_id,dset_id, &
                       error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
      call h5dwrite_f(dset_id,h5t_std_i64le,intarray,dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Data set write failed for '//aname
        errval%errBool = .true.
        return
      endif
      deallocate(intarray)
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine dump_intl_1d
!-----------------------------------------------------------------------
!     subprogram 40. dump_h5_intl_dum
!-----------------------------------------------------------------------
      subroutine dump_h5_intl_dum(inid,aname,dims,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character(*), intent(in) :: aname
      integer(i4), dimension(:), intent(in) :: dims
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType) :: errval
    
      integer,parameter :: FAIL=-1
      integer :: error
      integer(i4) :: ii
      integer(HID_T) :: dspace_id, dset_id
      integer :: rank
      integer(HSIZE_T), allocatable :: tdims(:)
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      rank=size(dims)
      allocate(tdims(rank))
      if (h5in%doTranspose) then
        do ii=1,size(dims)
          tdims(size(dims)-ii+1)=dims(ii)
        enddo
      else
        tdims=dims
      endif
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      call h5screate_simple_f(rank,tdims,dspace_id,error)
      deallocate(tdims)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Create the data set.
!-----------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dcreate_f(inid,aname,h5t_std_i64le,dspace_id, &
                       dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine dump_h5_intl_dum
!-----------------------------------------------------------------------
!     subprogram 41. dump_rl
!     Write an hdf5 array + references to independent vars
!-------------------------------------------------------------------
      subroutine dump_rl(inid,aname,val,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character*(*), intent(in) :: aname
      real(r8), intent(in) :: val
      type(hdf5ErrorType) :: errval
      type(hdf5InOpts), intent(in) :: h5in
      integer,parameter :: FAIL=-1
      integer :: error
      integer(HID_T) :: dspace_id, dset_id
      integer(HSIZE_T), dimension(1) :: dims=0
!-------------------------------------------------------------------
      if(h5in%verbose) WRITE(*,*) 'Writing ', aname
!-------------------------------------------------------------------
!     Create the data space.
!-------------------------------------------------------------------
      call h5screate_f(H5S_SCALAR_F, dspace_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-------------------------------------------------------------------
!     Create the data set.
!-------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dcreate_f(inid,aname,h5t_ieee_f64le,dspace_id,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
      if(h5in%typeConvert) then
        call h5dwrite_f(dset_id,h5t_ieee_f32le,real(val,r4),dims,error)
      else
        call h5dwrite_f(dset_id,h5t_ieee_f64le,val,dims,error)
      endif
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Data set write failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine dump_rl
!-------------------------------------------------------------------
!     subprogram 42. dump_rls
!     Write an hdf5 array + references to independent vars
!-------------------------------------------------------------------
      subroutine dump_rls(inid,aname,val,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character*(*), intent(in) :: aname
      real(r4), intent(in) :: val
      type(hdf5ErrorType) :: errval
      type(hdf5InOpts), intent(in) :: h5in
      integer,parameter :: FAIL=-1
      integer :: error
      integer(HID_T) :: dspace_id, dset_id
      integer(HSIZE_T), dimension(1) :: dims=0
!-------------------------------------------------------------------
      if(h5in%verbose) WRITE(*,*) 'Writing ', aname
!-------------------------------------------------------------------
!     Create the data space.
!-------------------------------------------------------------------
      call h5screate_f(H5S_SCALAR_F, dspace_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-------------------------------------------------------------------
!     Create the data set.
!-------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dcreate_f(inid,aname,h5t_ieee_f32le,dspace_id,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
      call h5dwrite_f(dset_id,h5t_ieee_f32le,val,dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Data set write failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine dump_rls
!-----------------------------------------------------------------------
!     subprogram 43. dump_rl_1d
!     Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
      subroutine dump_rl_1d(inid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character*(*), intent(in) :: aname
      real(r8), dimension(:), intent(in) :: array
      type(hdf5ErrorType) :: errval
      type(hdf5InOpts), intent(in) :: h5in
      integer,parameter :: FAIL=-1
      integer :: error
      integer :: rank
      integer(HID_T) :: dspace_id, dset_id
      integer(HSIZE_T), dimension(1) :: dims
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      if(h5in%verbose) WRITE(*,*) 'Writing ', aname
      rank = 1;            dims(:) = (/size(array,1)/)
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      call h5screate_simple_f(rank,dims,dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Create the data set.
!-----------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dcreate_f(inid,aname,h5t_ieee_f64le,dspace_id,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
      if(h5in%typeConvert) then
       call h5dwrite_f(dset_id,h5t_ieee_f32le,real(array,r4),dims,error)
      else
       call h5dwrite_f(dset_id,h5t_ieee_f64le,array,dims,error)
      endif
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Data set write failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine dump_rl_1d
!-----------------------------------------------------------------------
!     subprogram 44. dump_rl_2d
!     Create a "simple dataset" and write it out.
!-----------------------------------------------------------------------
      subroutine dump_rl_2d(inid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character(*), intent(in) :: aname
      real(r8), dimension(:,:), intent(in) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType) :: errval
      integer,parameter :: FAIL=-1
      integer :: error
    
      integer(HID_T) :: dspace_id, dset_id
      integer :: rank
      integer(HSIZE_T) :: dims(2)
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      if(h5in%verbose) WRITE(*,*) 'Writing ', aname
      rank = 2
      if (h5in%doTranspose) then
        dims(2) = size(array,1);  dims(1) = size(array,2)
      else
        dims(1) = size(array,1);  dims(2) = size(array,2)
      endif
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      call h5screate_simple_f(rank,dims,dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Create the data set.
!-----------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dcreate_f(inid,aname,h5t_ieee_f64le,dspace_id,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
      if (h5in%typeConvert) then
        if (h5in%doTranspose) then
          call h5dwrite_f(dset_id,h5t_ieee_f32le,transpose(real(array,r4)), &
                          dims,error)
        else
          call h5dwrite_f(dset_id,h5t_ieee_f32le,real(array,r4),dims,error)
        endif
      else
        if (h5in%doTranspose) then
          call h5dwrite_f(dset_id,h5t_ieee_f64le,transpose(array),dims, &
                          error)
        else
          call h5dwrite_f(dset_id,h5t_ieee_f64le,array,dims,error)
        endif
      endif
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Writing data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine dump_rl_2d
!-----------------------------------------------------------------------
!     subprogram 45. dump_rl_3d
!-----------------------------------------------------------------------
      subroutine dump_rl_3d(inid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character(*), intent(in) :: aname
      real(r8), dimension(:,:,:), intent(in) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType) :: errval
      integer,parameter :: FAIL=-1
      integer :: error
      integer :: i,j
      real(r8), dimension(:,:,:), allocatable :: tmparray
    
      integer(HID_T) :: dspace_id, dset_id
      integer :: rank
      integer(HSIZE_T) :: dims(3)
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) 'Writing ', aname
      rank = 3
      if (h5in%doTranspose) then
        dims(3)=size(array,1);dims(2)=size(array,2);dims(1)=size(array,3)
        allocate(tmparray(dims(1),dims(2),dims(3)))
        do i=1,dims(1); do j=1,dims(2)
          tmparray(i,j,:)=array(:,j,i)
        enddo; enddo
      else
        dims(1)=size(array,1);dims(2)=size(array,2);dims(3)=size(array,3)
      endif
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      call h5screate_simple_f(rank,dims,dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Create the data set.
!-----------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dcreate_f(inid,aname,h5t_ieee_f64le,dspace_id,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
      if (h5in%typeConvert) then
        if (h5in%doTranspose) then
          call h5dwrite_f(dset_id,h5t_ieee_f32le,real(tmparray,r4),dims,error)
          deallocate(tmparray)
        else
          call h5dwrite_f(dset_id,h5t_ieee_f32le,real(array,r4),dims,error)
        endif
      else
        if (h5in%doTranspose) then
          call h5dwrite_f(dset_id,h5t_ieee_f64le,tmparray,dims,error)
          deallocate(tmparray)
        else
          call h5dwrite_f(dset_id,h5t_ieee_f64le,array,dims,error)
        endif
      endif
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Writing data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      if(allocated(tmparray)) deallocate(tmparray)
      return
      end subroutine dump_rl_3d
!-----------------------------------------------------------------------
!     subprogram 46. dump_rl_4d
!-----------------------------------------------------------------------
      subroutine dump_rl_4d(inid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character(*), intent(in) :: aname
      real(r8), dimension(:,:,:,:), intent(in) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType) :: errval
      integer,parameter :: FAIL=-1
      integer :: error
      integer :: i,j,k
      real(r8), dimension(:,:,:,:), allocatable :: tmparray
    
      integer(HID_T) :: dspace_id, dset_id
      integer :: rank
      integer(HSIZE_T) :: dims(4)
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) 'Writing ', aname
      rank = 4
      if (h5in%doTranspose) then
        dims(4)=size(array,1);  dims(2)=size(array,3)
        dims(3)=size(array,2);  dims(1)=size(array,4)
        allocate(tmparray(dims(1),dims(2),dims(3),dims(4)))
        do i=1,dims(1); do j=1,dims(2); do k=1,dims(3)
          tmparray(i,j,k,:)=array(:,k,j,i)
        enddo; enddo; enddo
      else
        dims(1)=size(array,1);  dims(2)=size(array,2)
        dims(3)=size(array,3);  dims(4)=size(array,4)
      endif
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      call h5screate_simple_f(rank,dims,dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Create the data set.
!-----------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dcreate_f(inid,aname,h5t_ieee_f64le,dspace_id,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
      if (h5in%typeConvert) then
        if (h5in%doTranspose) then
          call h5dwrite_f(dset_id,h5t_ieee_f32le,real(tmparray,r4),dims, &
                          error)
          deallocate(tmparray)
        else
          call h5dwrite_f(dset_id,h5t_ieee_f32le,real(array,r4),dims, &
                          error)
        endif
      else
        if(h5in%doTranspose) then
          call h5dwrite_f(dset_id,h5t_ieee_f64le,tmparray,dims,error)
          deallocate(tmparray)
        else
          call h5dwrite_f(dset_id,h5t_ieee_f64le,array,dims,error)
        endif
      endif
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Writing data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      if (allocated(tmparray)) deallocate(tmparray)
      return
      end subroutine dump_rl_4d
!-----------------------------------------------------------------------
!     subprogram 47. dump_rl_5d
!-----------------------------------------------------------------------
      subroutine dump_rl_5d(inid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character(*), intent(in) :: aname
      real(r8), dimension(:,:,:,:,:), intent(in) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType) :: errval
      integer,parameter :: FAIL=-1
      integer :: error
      integer :: i,j,k,l
      real(r8), dimension(:,:,:,:,:), allocatable :: tmparray
    
      integer(HID_T) :: dspace_id, dset_id
      integer :: rank
      integer(HSIZE_T) :: dims(5)
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      if(h5in%verbose) WRITE(*,*) 'Writing ', aname
      rank = 5
      if (h5in%doTranspose) then
        dims(5)=size(array,1);
        dims(4)=size(array,2);  dims(2)=size(array,4)
        dims(3)=size(array,3);  dims(1)=size(array,5)
        allocate(tmparray(dims(1),dims(2),dims(3),dims(4),dims(5)))
        do i=1,dims(1); do j=1,dims(2); do k=1,dims(3); do l=1,dims(4)
          tmparray(i,j,k,l,:)=array(:,l,k,j,i)
        enddo; enddo; enddo; enddo
      else
        dims(1)=size(array,1);  dims(2)=size(array,2)
        dims(3)=size(array,3);  dims(4)=size(array,4)
        dims(5)=size(array,5);
      endif
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      call h5screate_simple_f(rank,dims,dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Create the data set.
!-----------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dcreate_f(inid,aname,h5t_ieee_f64le,dspace_id,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
      if (h5in%typeConvert) then
        if (h5in%doTranspose) then
          call h5dwrite_f(dset_id,h5t_ieee_f32le,real(tmparray,r4),dims, &
                          error)
          deallocate(tmparray)
        else
          call h5dwrite_f(dset_id,h5t_ieee_f32le,real(array,r4),dims, &
                          error)
        endif
      else
        if (h5in%doTranspose) then
          call h5dwrite_f(dset_id,h5t_ieee_f64le,tmparray,dims,error)
          deallocate(tmparray)
        else
          call h5dwrite_f(dset_id,h5t_ieee_f64le,array,dims,error)
        endif
      endif
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Writing data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      if (allocated(tmparray)) deallocate(tmparray)
      return
      end subroutine dump_rl_5d
!-----------------------------------------------------------------------
!     subprogram 48. dump_h5_rl_dum
!-----------------------------------------------------------------------
      subroutine dump_h5_rl_dum(inid,aname,dims,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character(*), intent(in) :: aname
      integer(i4), dimension(:), intent(in) :: dims
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType) :: errval
    
      integer,parameter :: FAIL=-1
      integer :: error
      integer(i4) :: ii
      integer(HID_T) :: dspace_id, dset_id
      integer :: rank
      integer(HSIZE_T), allocatable :: tdims(:)
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      rank=size(dims)
      allocate(tdims(rank))
      if (h5in%doTranspose) then
        do ii=1,size(dims)
          tdims(size(dims)-ii+1)=dims(ii)
        enddo
      else
        tdims=dims
      endif
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      call h5screate_simple_f(rank,tdims,dspace_id,error)
      deallocate(tdims)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Create the data set.
!-----------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dcreate_f(inid,aname,h5t_ieee_f64le,dspace_id,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine dump_h5_rl_dum
!-----------------------------------------------------------------------
!     subprogram 49. dump_rls_1d
!     Write an hdf5 array + references to independent vars
!-----------------------------------------------------------------------
      subroutine dump_rls_1d(inid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character*(*), intent(in) :: aname
      real(r4), dimension(:), intent(in) :: array
      type(hdf5ErrorType) :: errval
      type(hdf5InOpts), intent(in) :: h5in
      integer,parameter :: FAIL=-1
      integer :: error
      integer :: rank
      integer(HID_T) :: dspace_id, dset_id
      integer(HSIZE_T), dimension(1) :: dims
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      if(h5in%verbose) WRITE(*,*) 'Writing ', aname
      rank = 1;            dims(:) = (/size(array,1)/)
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      call h5screate_simple_f(rank,dims,dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Create the data set.
!-----------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dcreate_f(inid,aname,h5t_ieee_f32le,dspace_id,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
      call h5dwrite_f(dset_id,h5t_ieee_f32le,array,dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Data set write failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine dump_rls_1d
!-----------------------------------------------------------------------
!     subprogram 50. dump_rls_2d
!     Create a "simple dataset" and write it out.
!-----------------------------------------------------------------------
      subroutine dump_rls_2d(inid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character(*), intent(in) :: aname
      real(r4), dimension(:,:), intent(in) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType) :: errval
      integer,parameter :: FAIL=-1
      integer :: error
    
      integer(HID_T) :: dspace_id, dset_id
      integer :: rank
      integer(HSIZE_T) :: dims(2)
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) 'Writing ', aname
      rank = 2
      if (h5in%doTranspose) then
        dims(2) = size(array,1);  dims(1) = size(array,2)
      else
        dims(1) = size(array,1);  dims(2) = size(array,2)
      endif
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      call h5screate_simple_f(rank,dims,dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Create the data set.
!-----------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dcreate_f(inid,aname,h5t_ieee_f32le,dspace_id,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
      if (h5in%doTranspose) then
        call h5dwrite_f(dset_id,h5t_ieee_f32le,transpose(array),dims, &
                        error)
      else
        call h5dwrite_f(dset_id,h5t_ieee_f32le,array,dims,error)
      endif
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Writing data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine dump_rls_2d
!-----------------------------------------------------------------------
!     subprogram 51. dump_rls_3d
!-----------------------------------------------------------------------
      subroutine dump_rls_3d(inid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character(*), intent(in) :: aname
      real(r4), dimension(:,:,:), intent(in) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType) :: errval
      integer,parameter :: FAIL=-1
      integer :: error
      integer :: i,j
      real(r4), dimension(:,:,:), allocatable :: tmparray
    
      integer(HID_T) :: dspace_id, dset_id
      integer :: rank
      integer(HSIZE_T) :: dims(3)
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) 'Writing ', aname
      rank = 3
      if(h5in%doTranspose) then
        dims(3)=size(array,1);dims(2)=size(array,2);dims(1)=size(array,3)
        allocate(tmparray(dims(1),dims(2),dims(3)))
        do i=1,dims(1); do j=1,dims(2)
          tmparray(i,j,:)=array(:,j,i)
        enddo; enddo
      else
        dims(1)=size(array,1);dims(2)=size(array,2);dims(3)=size(array,3)
      endif
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      call h5screate_simple_f(rank,dims,dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Create the data set.
!-----------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dcreate_f(inid,aname,h5t_ieee_f32le,dspace_id,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
      if (h5in%doTranspose) then
        call h5dwrite_f(dset_id,h5t_ieee_f32le,tmparray,dims,error)
        deallocate(tmparray)
      else
        call h5dwrite_f(dset_id,h5t_ieee_f32le,array,dims,error)
      endif
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Writing data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine dump_rls_3d
!-----------------------------------------------------------------------
!     subprogram 52. dump_rls_4d
!-----------------------------------------------------------------------
      subroutine dump_rls_4d(inid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character(*), intent(in) :: aname
      real(r4), dimension(:,:,:,:), intent(in) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType) :: errval
      integer,parameter :: FAIL=-1
      integer :: error
      integer :: i,j,k
      real(r4), dimension(:,:,:,:), allocatable :: tmparray
    
      integer(HID_T) :: dspace_id, dset_id
      integer :: rank
      integer(HSIZE_T) :: dims(4)
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) 'Writing ', aname
      rank = 4
      if (h5in%doTranspose) then
        dims(4)=size(array,1);  dims(2)=size(array,3)
        dims(3)=size(array,2);  dims(1)=size(array,4)
        allocate(tmparray(dims(1),dims(2),dims(3),dims(4)))
        do i=1,dims(1); do j=1,dims(2); do k=1,dims(3)
          tmparray(i,j,k,:)=array(:,k,j,i)
        enddo; enddo; enddo
      else
        dims(1)=size(array,1);  dims(2)=size(array,2)
        dims(3)=size(array,3);  dims(4)=size(array,4)
      endif
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      call h5screate_simple_f(rank,dims,dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Create the data set.
!-----------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dcreate_f(inid,aname,h5t_ieee_f32le,dspace_id,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Write stored data to "name" data set.
!-----------------------------------------------------------------------
      if (h5in%doTranspose) then
        call h5dwrite_f(dset_id,h5t_ieee_f32le,tmparray,dims,error)
        deallocate(tmparray)
      else
        call h5dwrite_f(dset_id,h5t_ieee_f32le,array,dims,error)
      endif
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Writing data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine dump_rls_4d
!-----------------------------------------------------------------------
!     subprogram 53. dump_h5_rls_dum
!-----------------------------------------------------------------------
      subroutine dump_h5_rls_dum(inid,aname,dims,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character(*), intent(in) :: aname
      integer(i4), dimension(:), intent(in) :: dims
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType) :: errval
    
      integer,parameter :: FAIL=-1
      integer :: error
      integer(i4) :: ii
      integer(HID_T) :: dspace_id, dset_id  
      integer :: rank
      integer(HSIZE_T), allocatable :: tdims(:)
!-----------------------------------------------------------------------
!     Define the rank and dimensions of the data set to be created.
!-----------------------------------------------------------------------
      rank=size(dims)
      allocate(tdims(rank))
      if (h5in%doTranspose) then
        do ii=1,size(dims)
          tdims(size(dims)-ii+1)=dims(ii)
        enddo
      else
        tdims=dims
      endif
!-----------------------------------------------------------------------
!     Create the data space.
!-----------------------------------------------------------------------
      call h5screate_simple_f(rank,tdims,dspace_id,error)
      deallocate(tdims)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Create the data set.
!-----------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dcreate_f(inid,aname,h5t_ieee_f32le,dspace_id,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Create data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Add the VisSchema attributes
!-----------------------------------------------------------------------
      call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine dump_h5_rls_dum
!-----------------------------------------------------------------------
!     subprogram 54. add_h5_int
!     This adds data to an unlimited data space.  This is used for
!     writing things like time data
!     You pass in a scalar, and it adds it to a 1D array
!-----------------------------------------------------------------------
      subroutine add_h5_int(inid,aname,val,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character(*), intent(in) :: aname
      integer(i4), intent(in) :: val
      type(hdf5ErrorType) :: errval
      type(hdf5InOpts), intent(in) :: h5in
      integer,parameter :: FAIL=-1
      integer :: rank = 1
      integer :: error
      integer(HID_T) :: dspace_id, filespace
      integer(HID_T) :: dset_id
      integer(HSIZE_T), dimension(1) :: dims,maxdims,chunk_dims,extdims,offset
      integer(HSIZE_T), dimension(1) :: olddims,oldmaxdims
      integer(HID_T) :: cparms        !dataset creatation property identifier 
      LOGICAL(i4) :: dset_exists
!-----------------------------------------------------------------------
      if(h5in%verbose) WRITE(*,*) 'Writing ', aname
!-----------------------------------------------------------------------
!     First determine whether the dataset exists because different
!     whether data exists or not.
!-----------------------------------------------------------------------
      call h5lexists_f(inid, aname, dset_exists, error)
!-----------------------------------------------------------------------
!     Do the case of creating the data
!-----------------------------------------------------------------------
      if (.NOT. dset_exists) then
        !
        ! Create the data space with unlimited dimensions.
        maxdims = (/H5S_UNLIMITED_f/)
        dims = (/1/)
        call h5screate_simple_f(rank, dims, dspace_id, error, maxdims)
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Create data space failed for '//aname
          errval%errBool = .true.
          return
        endif
        ! Modify dataset creation properties, i.e. enable chunking
        call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, error)
        chunk_dims = (/1/)
        call h5pset_chunk_f(cparms, rank, chunk_dims, error)
    
        ! Create a new dataset within the file using cparms creation properties.
        call h5dcreate_f(inid,aname,h5t_std_i32le,dspace_id,dset_id,error,cparms)
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Create data set failed for '//aname
          errval%errBool = .true.
          return
        endif
        ! Write stored data to "name" data set.
        call h5dwrite_f(dset_id,h5t_std_i32le,val,dims,error)
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Data set write failed for '//aname
          errval%errBool = .true.
          return
        endif
        call h5pclose_f(cparms, error) !Close the property list.
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Close property list failed for '//aname
          errval%errBool = .true.
          return
        endif
        ! Add the VisSchema attributes on the first step
        call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Case for appending the dataset
!-----------------------------------------------------------------------
      else
        call h5dopen_f(inid,aname, dset_id, error)
        
        ! Open filespace in existing dataset and get existing shape and size of data
        call h5dget_space_f(dset_id, filespace, error)
    
        call H5Sget_simple_extent_dims_f(filespace, olddims, oldmaxdims,error)
        
        ! Extend the dataset. This call assures that dataset has the space
        dims = (/1/)
        extdims=dims
        extdims(1)=olddims(1)+dims(1)
        call h5dextend_f(dset_id, extdims, error)
    
        ! Define memory space
        call h5screate_simple_f(rank, dims, dspace_id, error)
    
        ! Open filespace in existing dataset and get existing shape and size of data
        call h5dget_space_f(dset_id, filespace, error)
    
        offset=0
        offset(1) = olddims(1)
        call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error) 
        
        ! Write the data to the hyperslab.
        call H5Dwrite_f(dset_id,h5t_std_i32le,val,dims,error,  &
                        file_space_id=filespace,mem_space_id=dspace_id)
        call h5sclose_f(filespace, error)
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine add_h5_int
!-----------------------------------------------------------------------
!     subprogram 55. add_h5_dbl
!     This adds data to an unlimited data space.  This is used for
!     writing things like time data
!     You pass in a scalar, and it adds it to a 1D array
!-----------------------------------------------------------------------
      subroutine add_h5_dbl(inid,aname,val,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character(*), intent(in) :: aname
      double precision, intent(in) :: val
      type(hdf5ErrorType) :: errval
      type(hdf5InOpts), intent(in) :: h5in
      integer,parameter :: FAIL=-1
      integer :: rank = 1
      integer :: error
      integer(HID_T) :: dspace_id=0, dset_id=0, filespace=0
      integer(HSIZE_T), dimension(1) :: dims=0,maxdims,extdims,offset
      integer(HSIZE_T), dimension(1) :: olddims,oldmaxdims
      integer(HID_T) :: cparms        !dataset creatation property identifier 
      LOGICAL(i4) :: dset_exists
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) 'Writing ', aname
!-----------------------------------------------------------------------
!     First determine whether the dataset exists because different
!     whether data exists or not.
!-----------------------------------------------------------------------
      call h5lexists_f(inid, aname, dset_exists, error)
!-----------------------------------------------------------------------
!     Do the case of creating the data
!-----------------------------------------------------------------------
      if (.NOT. dset_exists) then
        !
        ! Create the data space with unlimited dimensions.
        !
        maxdims = (/H5S_UNLIMITED_f/)
        dims = (/1/)
        call h5screate_simple_f(rank, dims, dspace_id, error, maxdims)
         if (error==FAIL) then
           errval%errorMsg = 'ERROR: Create data space failed for '//aname
           errval%errBool = .true.
           return
        endif
        ! Modify dataset creation properties, i.e. enable chunking
        call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, error)
        call h5pset_chunk_f(cparms, rank, dims, error)
    
        ! Create a new dataset within the file using cparms creation properties.
        call h5dcreate_f(inid,aname,h5t_ieee_f64le,dspace_id,dset_id,error,cparms)
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Create data set failed for '//aname
          errval%errBool = .true.
          return
        endif
        call h5pclose_f(cparms, error) !Close the property list.
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Close property list failed for '//aname
          errval%errBool = .true.
          return
        endif
    
        ! Write stored data to "name" data set.
        !call h5dwrite_f(dset_id,h5t_ieee_f64le,val,dims,error)
        call H5Dwrite_f(dset_id,h5t_ieee_f64le,val,dims,error)
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Data set write failed for '//aname
          errval%errBool = .true.
          return
        endif
        ! Add the VisSchema attributes on the first time step
        call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Case for appending the dataset
!-----------------------------------------------------------------------
      else
        call h5dopen_f(inid,aname, dset_id, error)
        
        ! Open filespace in existing dataset and get existing shape and size of data
        call h5dget_space_f(dset_id, filespace, error)
    
        call H5Sget_simple_extent_dims_f(filespace, olddims, oldmaxdims,error)
        
        ! Extend the dataset. This call assures that dataset has the space
        dims = (/1/)
        extdims=dims
        extdims(1)=olddims(1)+dims(1)
        call h5dextend_f(dset_id, extdims, error)
    
        ! Define memory space
        call h5screate_simple_f(rank, dims, dspace_id, error)
    
        ! Open filespace in existing dataset and get existing shape and size of data
        call h5dget_space_f(dset_id, filespace, error)
    
        offset=0
        offset(1) = olddims(1)
        call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error) 
        
        ! Write the data to the hyperslab.
        call H5Dwrite_f(dset_id,h5t_ieee_f64le,val,dims,error, &
                        file_space_id=filespace,mem_space_id=dspace_id)
        call h5sclose_f(filespace, error)
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine add_h5_dbl
!-----------------------------------------------------------------------
!     subprogram 56. add_h5_int_1d
!     This adds data to an unlimited data space.  This is used for
!     writing things like time data
!     You pass in a scalar, and it adds it to a array array
!-----------------------------------------------------------------------
      subroutine add_h5_int_1d(inid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character(*), intent(in) :: aname
      integer(i4), dimension(:), intent(in) :: array
      type(hdf5ErrorType) :: errval
      type(hdf5InOpts), intent(in) :: h5in
      integer,parameter :: FAIL=-1
      integer :: rank = 2
      integer :: asize
      integer :: error
      integer(HID_T) :: dspace_id, dset_id, filespace
      integer(HSIZE_T), dimension(2) :: dims,maxdims,chunk_dims,extdims,offset
      integer(HSIZE_T), dimension(2) :: olddims,oldmaxdims
      integer(HID_T) :: cparms        !dataset creatation property identifier 
      LOGICAL(i4) :: dset_exists
!-----------------------------------------------------------------------
      if(h5in%verbose) WRITE(*,*) 'Writing ', aname
!-----------------------------------------------------------------------
!     First determine whether the dataset exists because different
!     whether data exists or not.
!-----------------------------------------------------------------------
      call h5lexists_f(inid, aname, dset_exists, error)
      asize=size(array)
!-----------------------------------------------------------------------
!     Do the case of creating the data
!-----------------------------------------------------------------------
      if (.NOT. dset_exists) then
        !
        ! Create the data space with unlimited dimensions.
        maxdims = (/H5S_UNLIMITED_f, H5S_UNLIMITED_f/)
        ! For convenience, put the time step (extendible set, as the first index
        dims = (/1, asize/)
        call h5screate_simple_f(rank, dims, dspace_id, error, maxdims)
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Create data space failed for '//aname
          errval%errBool = .true.
          return
        endif
        ! Modify dataset creation properties, i.e. enable chunking
        call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, error)
        chunk_dims = dims
        call h5pset_chunk_f(cparms, rank, chunk_dims, error)
    
        ! Create a new dataset within the file using cparms creation properties.
        call h5dcreate_f(inid,aname,h5t_std_i32le,dspace_id,dset_id,error,cparms)
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Create data set failed for '//aname
          errval%errBool = .true.
          return
        endif
        ! Write stored data to "name" data set.
        call h5dwrite_f(dset_id,h5t_std_i32le,array,dims,error)
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Data set write failed for '//aname
          errval%errBool = .true.
          return
        endif
        call h5pclose_f(cparms, error) !Close the property list.
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Close property list failed for '//aname
          errval%errBool = .true.
          return
        endif
        ! Add the VisSchema attributes on the first step
        call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Case for appending the dataset
!-----------------------------------------------------------------------
      else
        call h5dopen_f(inid,aname, dset_id, error)
        
        ! Open filespace in existing dataset and get existing shape and size of data
        call h5dget_space_f(dset_id, filespace, error)
    
        call H5Sget_simple_extent_dims_f(filespace, olddims, oldmaxdims,error)
        
        ! Extend the dataset. This call assures that dataset has the space
        dims = (/1, asize/)
        extdims=dims
        extdims(1)=olddims(1)+dims(1)
        extdims(2)=asize
        call h5dextend_f(dset_id, extdims, error)
    
        ! Define memory space
        call h5screate_simple_f(rank, dims, dspace_id, error)
    
        ! Open filespace in existing dataset and get existing shape and size of data
        call h5dget_space_f(dset_id, filespace, error)
    
        offset=0
        offset(1) = olddims(1)
        call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error) 
        
        ! Write the data to the hyperslab.
        call H5Dwrite_f(dset_id,h5t_std_i32le,array,dims,error,  &
                        file_space_id=filespace,mem_space_id=dspace_id)
        call h5sclose_f(filespace, error)
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine add_h5_int_1d
!-----------------------------------------------------------------------
!     subprogram 57. add_h5_1d
!     This adds data to an unlimited data space.  This is used for
!     writing things like time data
!     You pass in a scalar, and it adds it to a array array
!-----------------------------------------------------------------------
      subroutine add_h5_1d(inid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character(*), intent(in) :: aname
      double precision, dimension(:), intent(in) :: array
      type(hdf5ErrorType) :: errval
      type(hdf5InOpts), intent(in) :: h5in
      integer,parameter :: FAIL=-1
      integer :: rank = 2
      integer :: asize
      integer :: error
      integer(HID_T) :: dspace_id, dset_id, filespace
      integer(HSIZE_T), dimension(2) :: dims,maxdims,chunk_dims,extdims,offset
      integer(HSIZE_T), dimension(2) :: olddims,oldmaxdims
      integer(HID_T) :: cparms        !dataset creatation property identifier 
      LOGICAL(i4) :: dset_exists
!-----------------------------------------------------------------------
      if(h5in%verbose) WRITE(*,*) 'Writing ', aname
!-----------------------------------------------------------------------
!     First determine whether the dataset exists because different
!     whether data exists or not.
!-----------------------------------------------------------------------
      call h5lexists_f(inid, aname, dset_exists, error)
      asize=size(array)

!-----------------------------------------------------------------------
!     Do the case of creating the data
!-----------------------------------------------------------------------
      if (.NOT. dset_exists) then
        !
        ! Create the data space with unlimited dimensions.
        maxdims = (/H5S_UNLIMITED_f, H5S_UNLIMITED_f/)
        ! For convenience, put the time step (extendible set, as the first index
        dims = (/1, asize/)
        call h5screate_simple_f(rank, dims, dspace_id, error, maxdims)
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Create data space failed for '//aname
          errval%errBool = .true.
          return
        endif
        ! Modify dataset creation properties, i.e. enable chunking
        call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, error)
        chunk_dims = dims
        call h5pset_chunk_f(cparms, rank, chunk_dims, error)
    
        ! Create a new dataset within the file using cparms creation properties.
        call h5dcreate_f(inid,aname,h5t_ieee_f64le,dspace_id,dset_id,error,cparms)
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Create data set failed for '//aname
          errval%errBool = .true.
          return
        endif
        ! Write stored data to "name" data set.
        call h5dwrite_f(dset_id,h5t_ieee_f64le,array,dims,error)
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Data set write failed for '//aname
          errval%errBool = .true.
          return
        endif
        call h5pclose_f(cparms, error) !Close the property list.
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Close property list failed for '//aname
          errval%errBool = .true.
          return
        endif
        ! Add the VisSchema attributes on the first step
        call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Case for appending the dataset
!-----------------------------------------------------------------------
      else
        call h5dopen_f(inid,aname, dset_id, error)
        
        ! Open filespace in existing dataset and get existing shape and size of data
        call h5dget_space_f(dset_id, filespace, error)
    
        call H5Sget_simple_extent_dims_f(filespace, olddims, oldmaxdims,error)
        
        ! Extend the dataset. This call assures that dataset has the space
        dims = (/1, asize/)
        extdims=dims
        extdims(1)=olddims(1)+dims(1)
        extdims(2)=asize
        call h5dextend_f(dset_id, extdims, error)
    
        ! Define memory space
        call h5screate_simple_f(rank, dims, dspace_id, error)
    
        ! Open filespace in existing dataset and get existing shape and size of data
        call h5dget_space_f(dset_id, filespace, error)
    
        offset=0
        offset(1) = olddims(1)
        call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error) 
        
        ! Write the data to the hyperslab.
        call H5Dwrite_f(dset_id,h5t_ieee_f64le,array,dims,error,  &
                        file_space_id=filespace,mem_space_id=dspace_id)
        call h5sclose_f(filespace, error)
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine add_h5_1d
!-----------------------------------------------------------------------
!     subprogram 58. add_h5_2d
!     This adds data to an unlimited data space.  This is used for
!     writing things like time data
!     You pass in a scalar, and it adds it to a array array
!-----------------------------------------------------------------------
      subroutine add_h5_2d(inid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character(*), intent(in) :: aname
      double precision, dimension(:,:), intent(in) :: array
      type(hdf5ErrorType) :: errval
      type(hdf5InOpts), intent(in) :: h5in
      integer,parameter :: FAIL=-1
      integer :: rank = 3
      integer :: error
      integer(HID_T) :: dspace_id, dset_id, filespace
      integer, dimension(2)          :: asize
      integer(HSIZE_T), dimension(3) :: dims,maxdims,chunk_dims,extdims,offset
      integer(HSIZE_T), dimension(3) :: olddims,oldmaxdims
      integer(HID_T) :: cparms        !dataset creatation property identifier 
      LOGICAL(i4) :: dset_exists
!-----------------------------------------------------------------------
      if(h5in%verbose) WRITE(*,*) 'Writing ', aname
!-----------------------------------------------------------------------
!     First determine whether the dataset exists because different
!     whether data exists or not.
!-----------------------------------------------------------------------
      call h5lexists_f(inid, aname, dset_exists, error)
      asize(1)=size(array,1)
      asize(2)=size(array,2)
      if (h5in%verbose) WRITE(*,*) 'Writing ', aname
    
      if (h5in%doTranspose) then
        dims = (/1, asize(2), asize(1)/)
      else
        dims = (/1, asize(1), asize(2)/)
      endif
!-----------------------------------------------------------------------
!     Do the case of creating the data
!-----------------------------------------------------------------------
      if (.NOT. dset_exists) then
        !
        ! Create the data space with unlimited dimensions.
        maxdims = (/H5S_UNLIMITED_f, H5S_UNLIMITED_f, H5S_UNLIMITED_f/)
        ! For convenience, put the time step (extendible set, as the first index
        call h5screate_simple_f(rank, dims, dspace_id, error, maxdims)
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Create data space failed for '//aname
          errval%errBool = .true.
          return
        endif
        ! Modify dataset creation properties, i.e. enable chunking
        call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, error)
        chunk_dims = dims
        call h5pset_chunk_f(cparms, rank, chunk_dims, error)
    
        ! Create a new dataset within the file using cparms creation properties.
        call h5dcreate_f(inid,aname,h5t_ieee_f64le,dspace_id,dset_id,error,cparms)
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Create data set failed for '//aname
          errval%errBool = .true.
          return
        endif
        ! Write stored data to "name" data set.
        if (h5in%doTranspose) then
          call h5dwrite_f(dset_id,h5t_ieee_f64le,transpose(array),dims,error)
        else
          call h5dwrite_f(dset_id,h5t_ieee_f64le,array,dims,error)
        endif
    
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Data set write failed for '//aname
          errval%errBool = .true.
          return
        endif
        call h5pclose_f(cparms, error) !Close the property list.
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Close property list failed for '//aname
          errval%errBool = .true.
          return
        endif
        ! Add the VisSchema attributes on the first step
        call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Case for appending the dataset
!-----------------------------------------------------------------------
      else
        call h5dopen_f(inid,aname, dset_id, error)
        
        ! Open filespace in existing dataset and get existing shape and size of data
        call h5dget_space_f(dset_id, filespace, error)
    
        call H5Sget_simple_extent_dims_f(filespace, olddims, oldmaxdims,error)
        
        ! Extend the dataset. This call assures that dataset has the space
        extdims=dims
        extdims(1)=olddims(1)+dims(1)
        !extdims(2)=asize(1)
        extdims(2)=dims(2)
        !extdims(3)=asize(2)
        extdims(3)=dims(3)
    
        call h5dextend_f(dset_id, extdims, error)
    
        ! Define memory space
        call h5screate_simple_f(rank, dims, dspace_id, error)
    
        ! Open filespace in existing dataset and get existing shape and size of data
        call h5dget_space_f(dset_id, filespace, error)
    
        offset=0
        offset(1) = olddims(1)
        call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error) 
        
        ! Write the data to the hyperslab.
        if (h5in%doTranspose) then
          call H5Dwrite_f(dset_id,h5t_ieee_f64le,transpose(array),dims,error,  &
                        file_space_id=filespace,mem_space_id=dspace_id)
        else
          call H5Dwrite_f(dset_id,h5t_ieee_f64le,array,dims,error,  &
                        file_space_id=filespace,mem_space_id=dspace_id)
        endif
    
        call h5sclose_f(filespace, error)
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine add_h5_2d
!-----------------------------------------------------------------------
!     subprogram 59. add_h5_3d
!     This adds data to an unlimited data space.  This is used for
!     writing things like time data
!     You pass in a scalar, and it adds it to a array array
!-----------------------------------------------------------------------
      subroutine add_h5_3d(inid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character(*), intent(in) :: aname
      double precision, dimension(:,:,:), intent(in) :: array
      type(hdf5ErrorType) :: errval
      type(hdf5InOpts), intent(in) :: h5in
      integer,parameter :: FAIL=-1
      integer :: rank = 4
      integer :: error
      integer(HID_T) :: dspace_id, dset_id, filespace
      integer, dimension(3)          :: asize
      integer(HSIZE_T), dimension(4) :: dims,maxdims,chunk_dims,extdims,offset
      integer(HSIZE_T), dimension(4) :: olddims,oldmaxdims
      integer(HID_T) :: cparms        !dataset creatation property identifier 
      LOGICAL(i4) :: dset_exists
      integer ::i,j
      real(r8), dimension(:,:,:), allocatable :: tmparray
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) 'Writing ', aname
!-----------------------------------------------------------------------
!     First determine whether the dataset exists because different
!     whether data exists or not.
!-----------------------------------------------------------------------
      call h5lexists_f(inid, aname, dset_exists, error)
      asize(1)=size(array,1)
      asize(2)=size(array,2)
      asize(3)=size(array,3)
      ! For convenience, put the time step (extendible set, as the first index
      if (h5in%verbose) WRITE(*,*) 'Writing ', aname
      if (h5in%doTranspose) then
        dims = (/1, asize(3), asize(2),asize(1)/)
        allocate(tmparray(dims(2),dims(3),dims(4)))
        do i=1,dims(2); do j=1,dims(3)
          tmparray(i,j,:)=array(:,j,i)
        enddo; enddo
      else
        dims = (/1, asize(1), asize(2),asize(3)/)
      endif
!-----------------------------------------------------------------------
!     Do the case of creating the data
!-----------------------------------------------------------------------
      if (.NOT. dset_exists) then
        !
        ! Create the data space with unlimited dimensions.
        maxdims = (/H5S_UNLIMITED_f, H5S_UNLIMITED_f, H5S_UNLIMITED_f, H5S_UNLIMITED_f/)
    
        call h5screate_simple_f(rank, dims, dspace_id, error, maxdims)
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Create data space failed for '//aname
          errval%errBool = .true.
          return
        endif
        ! Modify dataset creation properties, i.e. enable chunking
        call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, error)
        chunk_dims = dims
        call h5pset_chunk_f(cparms, rank, chunk_dims, error)
    
        ! Create a new dataset within the file using cparms creation properties.
        call h5dcreate_f(inid,aname,h5t_ieee_f64le,dspace_id,dset_id,error,cparms)
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Create data set failed for '//aname
          errval%errBool = .true.
          return
        endif
        ! Write stored data to "name" data set.
        if (h5in%doTranspose) then
          call h5dwrite_f(dset_id,h5t_ieee_f64le,tmparray,dims,error)
        else 
          call h5dwrite_f(dset_id,h5t_ieee_f64le,array,dims,error)
        endif
    
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Data set write failed for '//aname
          errval%errBool = .true.
          return
        endif
        call h5pclose_f(cparms, error) !Close the property list.
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Close property list failed for '//aname
          errval%errBool = .true.
          return
        endif
        ! Add the VisSchema attributes on the first step
        call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Case for appending the dataset
!-----------------------------------------------------------------------
      else
        call h5dopen_f(inid,aname, dset_id, error)
        
        ! Open filespace in existing dataset and get existing shape and size of data
        call h5dget_space_f(dset_id, filespace, error)
    
        call H5Sget_simple_extent_dims_f(filespace, olddims, oldmaxdims,error)
        
        ! Extend the dataset. This call assures that dataset has the space
        extdims=dims
        extdims(1)=olddims(1)+dims(1)
!        extdims(2)=asize(1)
!        extdims(3)=asize(2)
!        extdims(4)=asize(3)
        extdims(2)=dims(2)
        extdims(3)=dims(3)
        extdims(4)=dims(4)
    
        call h5dextend_f(dset_id, extdims, error)
    
        ! Define memory space
        call h5screate_simple_f(rank, dims, dspace_id, error)
    
        ! Open filespace in existing dataset and get existing shape and size of data
        call h5dget_space_f(dset_id, filespace, error)
    
        offset=0
        offset(1) = olddims(1)
        call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error) 
        
        ! Write the data to the hyperslab.
        if (h5in%doTranspose) then
          call H5Dwrite_f(dset_id,h5t_ieee_f64le,tmparray,dims,error,  &
                        file_space_id=filespace,mem_space_id=dspace_id)
        else 
          call H5Dwrite_f(dset_id,h5t_ieee_f64le,array,dims,error,  &
                        file_space_id=filespace,mem_space_id=dspace_id)
        endif
        call h5sclose_f(filespace, error)
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      if (allocated(tmparray)) deallocate(tmparray)
      return
      end subroutine add_h5_3d
!-----------------------------------------------------------------------
!     subprogram 60. add_h5_4d
!     This adds data to an unlimited data space.  This is used for
!     writing things like time data
!     You pass in a scalar, and it adds it to a array array
!-----------------------------------------------------------------------
      subroutine add_h5_4d(inid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: inid
      character(*), intent(in) :: aname
      double precision, dimension(:,:,:,:), intent(in) :: array
      type(hdf5ErrorType) :: errval
      type(hdf5InOpts), intent(in) :: h5in
      integer,parameter :: FAIL=-1
      integer :: rank = 5
      integer :: error
      integer(HID_T) :: dspace_id, dset_id, filespace
      integer, dimension(4)          :: asize
      integer(HSIZE_T), dimension(5) :: dims,maxdims,chunk_dims,extdims,offset
      integer(HSIZE_T), dimension(5) :: olddims,oldmaxdims
      integer(HID_T) :: cparms        !dataset creatation property identifier 
      LOGICAL(i4) :: dset_exists
      integer :: i,j,k
      real(r8), dimension(:,:,:,:), allocatable :: tmparray
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) 'Writing ', aname
!-----------------------------------------------------------------------
!     First determine whether the dataset exists because different
!     whether data exists or not.
!-----------------------------------------------------------------------
      call h5lexists_f(inid, aname, dset_exists, error)
      asize(1)=size(array,1)
      asize(2)=size(array,2)
      asize(3)=size(array,3)
      asize(4)=size(array,4)
    
      if (h5in%verbose) WRITE(*,*) 'Writing ', aname
      if (h5in%doTranspose) then
        dims = (/1, asize(4), asize(3),asize(2),asize(1)/)
        allocate(tmparray(dims(2),dims(3),dims(4),dims(5)))
        do i=1,dims(2); do j=1,dims(3); do k=1,dims(4)
          tmparray(i,j,k,:)=array(:,k,j,i)
        enddo; enddo; enddo
      else
        dims = (/1, asize(1), asize(2),asize(3),asize(4)/)
      endif
!-----------------------------------------------------------------------
!     Do the case of creating the data
!-----------------------------------------------------------------------
      if (.NOT. dset_exists) then
        !
        ! Create the data space with unlimited dimensions.
        maxdims = (/H5S_UNLIMITED_f, H5S_UNLIMITED_f, H5S_UNLIMITED_f, H5S_UNLIMITED_f, H5S_UNLIMITED_f/)
        ! For convenience, put the time step (extendible set, as the first index
        call h5screate_simple_f(rank, dims, dspace_id, error, maxdims)
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Create data space failed for '//aname
          errval%errBool = .true.
          return
        endif
        ! Modify dataset creation properties, i.e. enable chunking
        call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, error)
        chunk_dims = dims
        call h5pset_chunk_f(cparms, rank, chunk_dims, error)
    
        ! Create a new dataset within the file using cparms creation properties.
        call h5dcreate_f(inid,aname,h5t_ieee_f64le,dspace_id,dset_id,error,cparms)
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Create data set failed for '//aname
          errval%errBool = .true.
          return
        endif
        ! Write stored data to "name" data set.
        if (h5in%doTranspose) then
          call h5dwrite_f(dset_id,h5t_ieee_f64le,tmparray,dims,error)
        else
          call h5dwrite_f(dset_id,h5t_ieee_f64le,array,dims,error)
        endif
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Data set write failed for '//aname
          errval%errBool = .true.
          return
        endif
        call h5pclose_f(cparms, error) !Close the property list.
        if (error==FAIL) then
          errval%errorMsg = 'ERROR: Close property list failed for '//aname
          errval%errBool = .true.
          return
        endif
        ! Add the VisSchema attributes on the first step
        call dump_h5in_attributes(dset_id,h5in,errval)
!-----------------------------------------------------------------------
!     Case for appending the dataset
!-----------------------------------------------------------------------
      else
        call h5dopen_f(inid,aname, dset_id, error)
        
        ! Open filespace in existing dataset and get existing shape and size of data
        call h5dget_space_f(dset_id, filespace, error)
    
        call H5Sget_simple_extent_dims_f(filespace, olddims, oldmaxdims,error)
        
        ! Extend the dataset. This call assures that dataset has the space
        extdims=dims
        extdims(1)=olddims(1)+dims(1)
!        extdims(2)=asize(1)
!        extdims(3)=asize(2)
!        extdims(4)=asize(3)
!        extdims(5)=asize(4)
        extdims(2)=dims(2)
        extdims(3)=dims(3)
        extdims(4)=dims(4)
        extdims(5)=dims(5)
        call h5dextend_f(dset_id, extdims, error)
    
        ! Define memory space
        call h5screate_simple_f(rank, dims, dspace_id, error)
    
        ! Open filespace in existing dataset and get existing shape and size of data
        call h5dget_space_f(dset_id, filespace, error)
    
        offset=0
        offset(1) = olddims(1)
        call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error) 
        
        ! Write the data to the hyperslab.
        if (h5in%doTranspose) then
          call H5Dwrite_f(dset_id,h5t_ieee_f64le,tmparray,dims,error,  &
                        file_space_id=filespace,mem_space_id=dspace_id)
        else
          call H5Dwrite_f(dset_id,h5t_ieee_f64le,array,dims,error,  &
                        file_space_id=filespace,mem_space_id=dspace_id)
        endif
        call h5sclose_f(filespace, error)
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset and dataspace
!-----------------------------------------------------------------------
      call h5sclose_f(dspace_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data space failed for '//aname
        errval%errBool = .true.
        return
      endif
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      if (allocated(tmparray)) deallocate(tmparray)
      return
      end subroutine add_h5_4d
!-----------------------------------------------------------------------
!     subprogram 61. dset_dims
!     Read the dimensions of dataset associated (for internal use)
!-----------------------------------------------------------------------
      subroutine dset_dims(dset_id,dims,errval)
      integer(HID_T), intent(in) :: dset_id
      integer(HSIZE_T), dimension(:), intent(inout) :: dims
      type(hdf5ErrorType), intent(inout) :: errval
      integer,parameter :: FAIL=-1
      integer :: error
    
      integer(HSIZE_T), dimension(:), allocatable :: maxdims
      integer(HID_T) dspace_id
!-----------------------------------------------------------------------
!     Get dataset's dataspace handle.
!-----------------------------------------------------------------------
      allocate(maxdims(size(dims)))
      call h5dget_space_f(dset_id,dspace_id,error)
!-----------------------------------------------------------------------
!     Get dataspace's dimensinons.
!-----------------------------------------------------------------------
      call h5sget_simple_extent_dims_f(dspace_id,dims,maxdims,error)
!-----------------------------------------------------------------------
      errval%errBool = .false.
      deallocate(maxdims)
      return
      end subroutine dset_dims
!-----------------------------------------------------------------------
!     subprogram 62. read_dims
!     Read the dimensions of dataset associated with aname
!-----------------------------------------------------------------------
      subroutine read_dims(fid,aname,dims,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      integer(HSIZE_T), dimension(:), intent(inout) :: dims
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: errval
      integer,parameter :: FAIL=-1
      integer :: error
      integer(HID_T) :: dset_id

      integer(HSIZE_T), dimension(:), allocatable :: maxdims
      integer(HID_T) dspace_id
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading size of: ', aname
      call h5dopen_f(fid, aname,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Get dataset's dataspace handle.
!-----------------------------------------------------------------------
      allocate(maxdims(size(dims)))
      call h5dget_space_f(dset_id,dspace_id,error)
!-----------------------------------------------------------------------
!     Get dataspace's dimensinons.
!-----------------------------------------------------------------------
      call h5sget_simple_extent_dims_f(dspace_id,dims,maxdims,error)
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      deallocate(maxdims)
      return
      end subroutine read_dims
!-----------------------------------------------------------------------
!     subprogram 63. read_ndims
!     Read the number of dimensions of dataset associated with aname
!-----------------------------------------------------------------------
      subroutine read_ndims(fid,aname,rank,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      integer, intent(inout) :: rank
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: errval
      integer,parameter :: FAIL=-1
      integer :: error
      integer(HID_T) :: dset_id

      integer(HID_T) dspace_id
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading dimensions of: ', aname
      call h5dopen_f(fid, aname,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Get dataset's dataspace handle.
!-----------------------------------------------------------------------
      call h5dget_space_f(dset_id,dspace_id,error)
!-----------------------------------------------------------------------
!     Get dataspace's number of dimensinons.
!-----------------------------------------------------------------------
      call h5sget_simple_extent_ndims_f(dspace_id,rank,error)
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_ndims
!-----------------------------------------------------------------------
!     subprogram 64. read_ch_sc
!     Read simple data set
!-----------------------------------------------------------------------
      subroutine read_ch_sc(fid,aname,dataset,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      character*(*), intent(inout) :: dataset
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: errval
      integer,parameter :: FAIL=-1
      integer :: error

      integer(HID_T) :: dtype_id      ! Datatype identifier
      integer(HID_T) :: dset_id       ! Dataset identifier
      integer(SIZE_T) :: dsetlen    ! Length of the string
      integer(HSIZE_T), dimension(1) :: dims,fdims
!-----------------------------------------------------------------------
!     Open the dataset specified by aname.
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading character string: ', aname
      call h5dopen_f(fid,aname,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-------------------------------------------------------------------
!     Create the datatype.
!-------------------------------------------------------------------
      dims(1)=1
      dsetlen=len(dataset)
      call h5tcopy_f(h5t_native_character,dtype_id,error)
      call h5tset_size_f(dtype_id,dsetlen,error)
!-------------------------------------------------------------------
!     Read data set.
!-------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dread_f(dset_id,dtype_id,dataset,dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        call h5dclose_f(dset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset.
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_ch_sc
!-----------------------------------------------------------------------
!     subprogram 65. read_ch_vec
!     Read simple data set: 1d array
      ! H5: needs to be tested!
!-----------------------------------------------------------------------
      subroutine read_ch_vec(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      character*(*), dimension(:), intent(inout) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: errval
      integer,parameter :: FAIL=-1
      integer :: error
    
      integer(i4) :: i
      integer(HID_T) :: dtype_id      ! Datatype identifier
      integer(HID_T) :: dset_id       ! Dataset identifier
      integer(HSIZE_T), dimension(1) :: adims  ! Array dimension
      integer(SIZE_T) :: arrlen    ! Length of the array
      integer(HSIZE_T), dimension(1) :: dims,fdims
      integer(SIZE_T) :: ari
!-----------------------------------------------------------------------
!     Open the dataset specified by aname.
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading string array: ', aname
      call h5dopen_f(fid,aname,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-------------------------------------------------------------------
!     Check size and create the datatype.
!-------------------------------------------------------------------
      dims(1)=size(array)
      call dset_dims(dset_id,fdims,errval)
      call check_dims(dims,fdims,errval)
      if (errval%errBool) return
      adims(:) = (/dims(1)/)
      arrlen=0
      do i=1,adims(1)
        ari = len(array(i))
        arrlen = max(arrlen,ari)
      enddo
      call h5tcopy_f(h5t_native_character,dtype_id,error)
      call h5tset_size_f(dtype_id,arrlen,error)
!-------------------------------------------------------------------
!     Read data set.
!-------------------------------------------------------------------
      if(h5in%debug) write(*,*) 'creating dataspace ',aname,' with dims',dims
      call h5dread_f(dset_id,dtype_id,array(1:arrlen),dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        call h5dclose_f(dset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset.
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_ch_vec
!-----------------------------------------------------------------------
!     subprogram 66. read_int
!     Read simple data scalar
!-----------------------------------------------------------------------
      subroutine read_int(fid,aname,val,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      integer(i4), intent(inout) :: val
      type(hdf5ErrorType), intent(inout) :: errval
      integer(HSIZE_T), dimension(1) :: dims
      integer,parameter :: FAIL=-1
      integer :: error
      integer(HID_T) :: dset_id
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading integer value: ', aname
      call h5dopen_f(fid, aname, dset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Read data set
!-----------------------------------------------------------------------
      if (errval%errBool) return
      call h5dread_f(dset_id,h5t_std_i32le,val,dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        call h5dclose_f(dset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_int
!-----------------------------------------------------------------------
!     subprogram 67. read_int_1d
!     Read simple data set: 1d array
!-----------------------------------------------------------------------
      subroutine read_int_1d(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      integer(i4), dimension(:), intent(inout) :: array
      type(hdf5ErrorType), intent(inout) :: errval
      integer(HSIZE_T), dimension(1) :: dims, fdims
      integer,parameter :: FAIL=-1
      integer :: error
      integer(HID_T) :: dset_id
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading 1d i4 array: ', aname
      call h5dopen_f(fid, aname, dset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Read data set
!-----------------------------------------------------------------------
      dims(1)=size(array)
      call dset_dims(dset_id,fdims,errval)
      call check_dims(dims,fdims, errval)
      if (errval%errBool) return
      call h5dread_f(dset_id,h5t_std_i32le,array,dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        call h5dclose_f(dset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_int_1d
!-----------------------------------------------------------------------
!     subprogram 68. read_int_2d
!     Read simple data set: 2d array
!-----------------------------------------------------------------------
      subroutine read_int_2d(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      integer(i4), dimension(:,:), intent(inout) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: errval
      integer(HSIZE_T), dimension(2) :: dims, fdims
      integer,parameter :: FAIL=-1
      integer :: error 
      integer :: i
      real(i4), dimension(:,:), allocatable :: tmparray
    
      integer(HID_T) :: dset_id
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading 2d i4 array: ', aname
      call h5dopen_f(fid, aname, dset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Check dims
!-----------------------------------------------------------------------
      if (.NOT. h5in%doTranspose) then
        dims(1)=size(array,1); dims(2)=size(array,2)
      else
        dims(1)=size(array,2); dims(2)=size(array,1)
      endif
      call dset_dims(dset_id,fdims,errval)
      call check_dims(dims,fdims, errval)
      if (errval%errBool) then
        call h5dclose_f(dset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Read data set
!-----------------------------------------------------------------------
      if (.NOT. h5in%doTranspose) then
        call h5dread_f(dset_id,h5t_std_i32le,array,dims,error)
      else
        allocate(tmparray(dims(1),dims(2)))
        call h5dread_f(dset_id,h5t_std_i32le,tmparray,dims,error)
        do i=1,dims(1);      array(:,i)=tmparray(i,:);     enddo
        deallocate(tmparray)
      endif
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Read data set failed for '//aname
        errval%errBool = .true.
        call h5dclose_f(dset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_int_2d
!-----------------------------------------------------------------------
!     subprogram 69. read_intl
!     Read simple data set: scalar
!-----------------------------------------------------------------------
      subroutine read_intl(fid,aname,val,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      integer(i8), intent(inout) :: val
      type(hdf5ErrorType), intent(inout) :: errval
      integer(HSIZE_T), dimension(1) :: dims
      integer,parameter :: FAIL=-1
      integer :: error
      integer(HID_T) :: dset_id
    
      integer(i4) :: intval
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading integer value: ', aname
      call h5dopen_f(fid, aname, dset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Read data set
!-----------------------------------------------------------------------
      if (errval%errBool) return
      call h5dread_f(dset_id,h5t_std_i64le,intval,dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        call h5dclose_f(dset_id,error)
        return
      endif
      val = intval
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_intl
!-----------------------------------------------------------------------
!     subprogram 70. read_intl_1d
!     Read simple data set: 1d array
!-----------------------------------------------------------------------
      subroutine read_intl_1d(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      integer(i8), dimension(:), intent(inout) :: array
      type(hdf5ErrorType), intent(inout) :: errval
      integer(HSIZE_T), dimension(1) :: dims, fdims
      integer,parameter :: FAIL=-1
      integer :: error
      integer(HID_T) :: dset_id
    
      integer(i4), dimension(:), allocatable :: intarray
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading 1d i8 array: ', aname
      call h5dopen_f(fid, aname, dset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Read data set
!-----------------------------------------------------------------------
      dims(1)=size(array)
      allocate(intarray(dims(1)))
!      intarray=array  ! PRE-necessary?
      call dset_dims(dset_id,fdims,errval)
      call check_dims(dims,fdims, errval)
      if (errval%errBool) return
      call h5dread_f(dset_id,h5t_std_i64le,intarray,dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        call h5dclose_f(dset_id,error)
        return
      endif
      array = intarray
      deallocate(intarray)
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_intl_1d
!-----------------------------------------------------------------------
!     subprogram 71. read_rl
!     Read simple data set: scalar
!-----------------------------------------------------------------------
      subroutine read_rl(fid,aname,val,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      real(r8), intent(inout) :: val
      type(hdf5ErrorType), intent(inout) :: errval
      integer(HSIZE_T), dimension(1) :: dims
      integer,parameter :: FAIL=-1
      integer :: error
    
      integer(HID_T) :: dset_id
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading r8 val: ', aname
      call h5dopen_f(fid, aname, dset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Read data set
!-----------------------------------------------------------------------
      dims(1)=1
      if (errval%errBool) return
      call h5dread_f(dset_id,h5t_ieee_f64le,val,dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        call h5dclose_f(dset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_rl
!-----------------------------------------------------------------------
!     subprogram 72. read_rl_1d
!     Read simple data set: 1d array
!-----------------------------------------------------------------------
      subroutine read_rl_1d(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      real(r8), dimension(:), intent(inout) :: array
      type(hdf5ErrorType), intent(inout) :: errval
      integer(HSIZE_T), dimension(1) :: dims, fdims
      integer,parameter :: FAIL=-1
      integer :: error
    
      integer(HID_T) :: dset_id
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading 1d r8 array: ', aname
      call h5dopen_f(fid, aname, dset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Read data set
!-----------------------------------------------------------------------
      dims(1)=size(array)
      call dset_dims(dset_id,fdims,errval)
      call check_dims(dims,fdims, errval)
      if (errval%errBool) return
      call h5dread_f(dset_id,h5t_ieee_f64le,array,dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        call h5dclose_f(dset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_rl_1d
!-----------------------------------------------------------------------
!     subprogram 73. read_rl_2d
!     Read simple data set: 2d array
!-----------------------------------------------------------------------
      subroutine read_rl_2d(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      real(r8), dimension(:,:), intent(inout) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: errval
      integer(HSIZE_T), dimension(2) :: dims, fdims
      integer,parameter :: FAIL=-1
      integer :: error 
      integer :: i
      real(r8), dimension(:,:), allocatable :: tmparray
    
      integer(HID_T) :: dset_id
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading 2d r8 array: ', aname
      call h5dopen_f(fid, aname, dset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Check dims
!-----------------------------------------------------------------------
      if (.NOT. h5in%doTranspose) then
        dims(1)=size(array,1); dims(2)=size(array,2)
      else
        dims(1)=size(array,2); dims(2)=size(array,1)
      endif
      call dset_dims(dset_id,fdims,errval)
      call check_dims(dims,fdims, errval)
      if (errval%errBool) then
        call h5dclose_f(dset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Read data set
!-----------------------------------------------------------------------
      if (.NOT. h5in%doTranspose) then
        call h5dread_f(dset_id,h5t_ieee_f64le,array,dims,error)
      else
        allocate(tmparray(dims(1),dims(2)))
        call h5dread_f(dset_id,h5t_ieee_f64le,tmparray,dims,error)
        do i=1,dims(1);      array(:,i)=tmparray(i,:);     enddo
        deallocate(tmparray)
      endif
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Read data set failed for '//aname
        errval%errBool = .true.
        call h5dclose_f(dset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_rl_2d
!-----------------------------------------------------------------------
!     subprogram 74. read_rl_3d
!     Read simple data set: 3d array
!-----------------------------------------------------------------------
      subroutine read_rl_3d(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      real(r8), dimension(:,:,:), intent(inout) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: errval
      integer(HSIZE_T), dimension(3) :: dims, fdims
      integer,parameter :: FAIL=-1
      integer :: error
      integer :: i,j
      real(r8), dimension(:,:,:), allocatable :: tmparray
    
      integer(HID_T) :: dset_id
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if(h5in%verbose) WRITE(*,*) ' Reading 3d array: ', aname
      call h5dopen_f(fid, TRIM(aname), dset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Check dims
!-----------------------------------------------------------------------
      if (.NOT. h5in%doTranspose) then
        dims(1)=size(array,1)
        dims(2)=size(array,2)
        dims(3)=size(array,3)
      else
        dims(1)=size(array,3)
        dims(2)=size(array,2)
        dims(3)=size(array,1)
      endif
      call dset_dims(dset_id,fdims,errval)
      !call check_dims(dims,fdims, errval)
      if (errval%errBool) then
        call h5dclose_f(dset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Read data set
!-----------------------------------------------------------------------
      if (.NOT. h5in%doTranspose) then
        call h5dread_f(dset_id,h5t_ieee_f64le,array,dims,error)
      else
        allocate(tmparray(dims(1),dims(2),dims(3)))
        call h5dread_f(dset_id,h5t_ieee_f64le,tmparray,dims,error)
        do i=1,dims(1); do j=1,dims(2)
          array(:,j,i)=tmparray(i,j,:)
        enddo; enddo
        deallocate(tmparray)
      endif
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Read data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_rl_3d
!-----------------------------------------------------------------------
!     subprogram 75. read_rl_4d
!     Read simple data set: 4d array
!-----------------------------------------------------------------------
      subroutine read_rl_4d(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      real(r8), dimension(:,:,:,:), intent(inout) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: errval
      integer(HSIZE_T), dimension(4) :: dims, fdims
      integer,parameter :: FAIL=-1
      integer :: error
      integer :: i,j,k
      real(r8), dimension(:,:,:,:), allocatable :: tmparray
    
      integer(HID_T) :: dset_id
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading 4d array: ', aname
      call h5dopen_f(fid, TRIM(aname), dset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Check dims
!-----------------------------------------------------------------------
      if (.NOT. h5in%doTranspose) then
        dims(1)=size(array,1)
        dims(2)=size(array,2)
        dims(3)=size(array,3)
        dims(4)=size(array,4)
      else
        dims(1)=size(array,4)
        dims(2)=size(array,3)
        dims(3)=size(array,2)
        dims(4)=size(array,1)
      endif
      call dset_dims(dset_id,fdims,errval)
      !call check_dims(dims,fdims, errval)
      if (errval%errBool) then
        call h5dclose_f(dset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Read data set
!-----------------------------------------------------------------------
      if (.NOT. h5in%doTranspose) then
        call h5dread_f(dset_id,h5t_ieee_f64le,array,dims,error)
      else
        allocate(tmparray(dims(1),dims(2),dims(3),dims(4)))
        call h5dread_f(dset_id,h5t_ieee_f64le,tmparray,dims,error)
        do i=1,dims(1); do j=1,dims(2); do k=1,dims(3)
          array(:,k,j,i)=tmparray(i,j,k,:)
        enddo; enddo; enddo
        deallocate(tmparray)
      endif
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Read data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_rl_4d
!-----------------------------------------------------------------------
!     subprogram 76. read_5d
!     Read simple data set: 5d array
!-----------------------------------------------------------------------
      subroutine read_rl_5d(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      real(r8), dimension(:,:,:,:,:), intent(inout) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: errval
      integer(HSIZE_T), dimension(5) :: dims, fdims
      integer,parameter :: FAIL=-1
      integer :: error
      integer :: i,j,k,l
      real(r8), dimension(:,:,:,:,:), allocatable :: tmparray
    
      integer(HID_T) :: dset_id
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading 5d array: ', aname
      call h5dopen_f(fid, TRIM(aname), dset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Check dims
!-----------------------------------------------------------------------
      if (.NOT. h5in%doTranspose) then
        dims(1)=size(array,1)
        dims(2)=size(array,2)
        dims(3)=size(array,3)
        dims(4)=size(array,4)
        dims(5)=size(array,5)
      else
        dims(1)=size(array,5)
        dims(2)=size(array,4)
        dims(3)=size(array,3)
        dims(4)=size(array,2)
        dims(5)=size(array,1)
      endif
      call dset_dims(dset_id,fdims,errval)
      !call check_dims(dims,fdims, errval)
      if (errval%errBool) then
        call h5dclose_f(dset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Read data set
!-----------------------------------------------------------------------
      if (.NOT. h5in%doTranspose) then
        call h5dread_f(dset_id,h5t_ieee_f64le,array,dims,error)
      else
        allocate(tmparray(dims(1),dims(2),dims(3),dims(4),dims(5)))
        call h5dread_f(dset_id,h5t_ieee_f64le,tmparray,dims,error)
        do i=1,dims(1); do j=1,dims(2); do k=1,dims(3); do l=1,dims(4)
          array(:,l,k,j,i)=tmparray(i,j,k,l,:)
        enddo; enddo; enddo; enddo
        deallocate(tmparray)
      endif
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Read data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_rl_5d
!-----------------------------------------------------------------------
!     subprogram 77. read_attribute_int_sc
!     Read real scalar attribute
!-----------------------------------------------------------------------
      subroutine read_attribute_int_sc(fid,aname,val,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      integer(i4), intent(inout) :: val
      type(hdf5ErrorType), intent(inout) :: errval
      integer(HSIZE_T), dimension(1) :: dims
      integer,parameter :: FAIL=-1
      integer :: error
      integer(HID_T) aset_id
    
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading attribute: ', aname
      call h5aopen_name_f(fid, aname, aset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Read attribute
!-----------------------------------------------------------------------
      dims(1)=1
      call h5aread_f(aset_id,h5t_std_i32le,val,dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        call h5aclose_f(aset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5aclose_f(aset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_attribute_int_sc
!-----------------------------------------------------------------------
!     subprogram 78. read_attribute_int_vec
!     Read integer vector attribute
!-----------------------------------------------------------------------
      subroutine read_attribute_int_vec(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      integer(i4), dimension(:), intent(inout) :: array
      type(hdf5ErrorType), intent(inout) :: errval
      integer(HSIZE_T), dimension(1) :: dims
      integer,parameter :: FAIL=-1
      integer :: error
      integer(HID_T) aset_id
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading attribute: ', aname
      call h5aopen_name_f(fid, aname, aset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Read attribute
!-----------------------------------------------------------------------
      dims(1)=size(array)
      call h5aread_f(aset_id,h5t_std_i32le,array,dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        call h5aclose_f(aset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5aclose_f(aset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_attribute_int_vec
!-----------------------------------------------------------------------
!     subprogram 79. read_attribute_intl_sc
!     Read real scalar attribute
!-----------------------------------------------------------------------
      subroutine read_attribute_intl_sc(fid,aname,val,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      integer(i8), intent(inout) :: val
      type(hdf5ErrorType), intent(inout) :: errval
      integer(HSIZE_T), dimension(1) :: dims
      integer,parameter :: FAIL=-1
      integer :: error
      integer(HID_T) aset_id
    
      integer(i4) :: intval
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading attribute: ', aname
      call h5aopen_name_f(fid, aname, aset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Read attribute
!-----------------------------------------------------------------------
      dims(1)=1
      call h5aread_f(aset_id,h5t_std_i64le,intval,dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        call h5aclose_f(aset_id,error)
        return
      endif
      val = intval
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5aclose_f(aset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_attribute_intl_sc
!-----------------------------------------------------------------------
!     subprogram 80. read_attribute_intl_vec
!     Read integer vector attribute
!-----------------------------------------------------------------------
      subroutine read_attribute_intl_vec(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      integer(i8), dimension(:), intent(inout) :: array
      type(hdf5ErrorType), intent(inout) :: errval
      integer(HSIZE_T), dimension(1) :: dims
      integer,parameter :: FAIL=-1
      integer :: error
      integer(HID_T) aset_id
    
      integer(i4), dimension(:), allocatable :: intarray
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading attribute: ', aname
      call h5aopen_name_f(fid, aname, aset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Read attribute
!-----------------------------------------------------------------------
      dims(1)=size(array)
      allocate(intarray(dims(1)))
      intarray = array
      call h5aread_f(aset_id,h5t_std_i64le,intarray,dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        call h5aclose_f(aset_id,error)
        return
      endif
      array = intarray
      deallocate(intarray)
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5aclose_f(aset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_attribute_intl_vec
!-----------------------------------------------------------------------
!     subprogram 81. read_attribute_rl_sc
!     Read real scalar attribute
!-----------------------------------------------------------------------
      subroutine read_attribute_rl_sc(fid,aname,val,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      real(r8), intent(inout) :: val
      type(hdf5ErrorType), intent(inout) :: errval
      integer(HSIZE_T), dimension(1) :: dims
      integer,parameter :: FAIL=-1
      integer :: error
    
      integer(HID_T) aset_id
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading attribute: ', aname
      call h5aopen_name_f(fid, aname, aset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Read attribute
!-----------------------------------------------------------------------
      dims(1)=1
      call h5aread_f(aset_id,h5t_ieee_f64le,val,dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        call h5aclose_f(aset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5aclose_f(aset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_attribute_rl_sc
!-----------------------------------------------------------------------
!     subprogram 82. read_attribute_rl_vec
!     Read integer vector attribute
!-----------------------------------------------------------------------
      subroutine read_attribute_rl_vec(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      real(r8), dimension(:), intent(inout) :: array
      type(hdf5ErrorType), intent(inout) :: errval
      integer(HSIZE_T), dimension(1) :: dims
      integer,parameter :: FAIL=-1
      integer :: error
      integer(HID_T) aset_id
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading attribute: ', aname
      call h5aopen_name_f(fid, aname, aset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Read attribute
!-----------------------------------------------------------------------
      dims(1)=size(array)
      call h5aread_f(aset_id,h5t_ieee_f64le,array,dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        call h5aclose_f(aset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5aclose_f(aset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_attribute_rl_vec
!-----------------------------------------------------------------------
!     subprogram 83. read_attribute_ch_sc
!     Read character scalar attribute
!-----------------------------------------------------------------------
      subroutine read_attribute_ch_sc(fid,aname,val,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      character*(*), intent(inout) :: val
      type(hdf5ErrorType), intent(inout) :: errval
      integer(HSIZE_T), dimension(1) :: dims
      integer,parameter :: FAIL=-1
      integer :: error
    
      integer(SIZE_T) :: attrlen    ! Length of the attribute string
      integer(HID_T) :: atype_id      ! Attribute Dataspace identifier
      integer(HID_T) aset_id
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading attribute: ', aname
      call h5aopen_name_f(fid, aname, aset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Read attribute
!-----------------------------------------------------------------------
      dims(1)=1
      attrlen = len(val)
    
      ! Create datatype for the attribute.
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attrlen, error)
    
      call h5aread_f(aset_id, atype_id, val, dims, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        call h5aclose_f(aset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5aclose_f(aset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_attribute_ch_sc
!-----------------------------------------------------------------------
!     subprogram 84. read_attribute_ch_vec
!     Read integer vector attribute
!-----------------------------------------------------------------------
      subroutine read_attribute_ch_vec(fid,aname,attribute,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      character*(*), dimension(:), intent(inout) :: attribute
      type(hdf5ErrorType), intent(inout) :: errval
    
      integer(HSIZE_T), dimension(1) :: adims
      integer,parameter :: FAIL=-1
      integer :: error
      integer(HID_T) aset_id
      integer(HID_T) :: atype_id      ! Attribute Dataspace identifier
      integer(SIZE_T) :: attrlen    ! Length of the attribute string
      integer(HSIZE_T), dimension(1) :: data_dims
      integer(SIZE_T) :: ati,i
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading attribute: ', aname
      call h5aopen_name_f(fid, aname, aset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Read attribute
!-----------------------------------------------------------------------
      data_dims(1) = size(attribute)
      adims(:) = (/data_dims(1)/)
      attrlen=0
      do i=1,data_dims(1)
        ! attrlen=max(attrlen,len(attribute(i)))
        ati = len(attribute(i))
        attrlen=max(attrlen,ati)
      enddo
    
      ! Create datatype for the attribute.
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attrlen, error)
    
      call h5aread_f(aset_id, atype_id, attribute, adims, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        call h5aclose_f(aset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5aclose_f(aset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_attribute_ch_vec
!-----------------------------------------------------------------------
!     subprogram 85. read_attribute_log_sc
!     Read real scalar attribute
!-----------------------------------------------------------------------
      subroutine read_attribute_log_sc(fid,aname,val,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      logical, intent(inout) :: val
      type(hdf5ErrorType), intent(inout) :: errval
      integer(HSIZE_T), dimension(1) :: dims
      integer,parameter :: FAIL=-1
      integer :: error
      integer(HID_T) aset_id
      integer(i4) :: ival
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading attribute: ', aname
      call h5aopen_name_f(fid, aname, aset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Read attribute
!-----------------------------------------------------------------------
      dims(1)=1
      call h5aread_f(aset_id,h5t_std_i32le,ival,dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        call h5aclose_f(aset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5aclose_f(aset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      val = .true.
      if (ival==0) val = .false.
      return
      end subroutine read_attribute_log_sc
!-----------------------------------------------------------------------
!     subprogram 86. read_attribute_log_vec
!     Read integer vector attribute
!-----------------------------------------------------------------------
      subroutine read_attribute_log_vec(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      logical, dimension(:), intent(inout) :: array
      type(hdf5ErrorType), intent(inout) :: errval
      integer(HSIZE_T), dimension(1) :: dims
      integer,parameter :: FAIL=-1
      integer :: error,ii
      integer(HID_T) aset_id
    
      integer(i4), dimension(:), allocatable :: intarray
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading attribute: ', aname
      call h5aopen_name_f(fid, aname, aset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Read attribute
!-----------------------------------------------------------------------
      dims(1)=size(array)
      allocate(intarray(dims(1)))
      call h5aread_f(aset_id,h5t_std_i32le,intarray,dims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        call h5aclose_f(aset_id,error)
        return
      endif
      do ii=1,dims(1)
        if (intarray(ii)==0) then
          array(ii) = .false.
        else
          array(ii) = .true.
        endif
      enddo
      deallocate(intarray)
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5aclose_f(aset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_attribute_log_vec
!-----------------------------------------------------------------------
!     subprogram 87. read_ch_sc_ex
!     Read character scalar if the object exists, but leave the input
!       unchanged if it doesn't
!-----------------------------------------------------------------------
      subroutine read_ch_sc_ex(fid,aname,dataset,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      character*(*), intent(inout) :: dataset
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: errval

      if (obj_exists(fid,aname,errval)) &
        call read_ch_sc(fid,aname,dataset,h5in,errval)

      return
      end subroutine read_ch_sc_ex
!-----------------------------------------------------------------------
!     subprogram 88. read_ch_vec_ex
!     Read character vector if the object exists, but leave the input
!       unchanged if it doesn't
!     Only read the number of elements included in the file
      ! H5: needs to be tested!
!-----------------------------------------------------------------------
      subroutine read_ch_vec_ex(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      character*(*), dimension(:), intent(inout) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: errval
      integer,parameter :: FAIL=-1
      integer :: error
    
      integer(i4) :: i
      integer(HID_T) :: dtype_id      ! Datatype identifier
      integer(HID_T) :: dset_id       ! Dataset identifier
      integer(HSIZE_T), dimension(1) :: adims  ! Array dimension
      integer(SIZE_T) :: arrlen    ! Length of the array
      integer(i4) :: ndims     ! Dataset rank
      integer(SIZE_T) :: dims    ! Dataset length
      integer(HSIZE_T), dimension(1) :: fdims
      integer(SIZE_T) :: ari
!-----------------------------------------------------------------------
!     Check that the object exists.
!-----------------------------------------------------------------------
      if (.NOT.obj_exists(fid,aname,errval)) return
!-----------------------------------------------------------------------
!     Read dataset as a scalar if it is not an array.
!-----------------------------------------------------------------------
      call read_ndims(fid,aname,ndims,h5in,errval)
      if (ndims < 1) then
        call read_ch_sc(fid,aname,array(1),h5in,errval)
        return
      elseif (ndims > 1) then
        write(*, *) "error: dataset is not 1d"
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Open the dataset specified by aname.
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading string array: ', aname
      call h5dopen_f(fid,aname,dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-------------------------------------------------------------------
!     Check size and create the datatype.
!-------------------------------------------------------------------
      dims=size(array,1)
      call dset_dims(dset_id,fdims,errval)
      if (errval%errBool) return
      if (fdims(1) > dims) then
        write(*, *) "error: reading ",aname," from hdf5"
        write(*, *) "fdims (", fdims(1), ")  >  dims (", dims, ")"
        write(*, *) "fdims = dims in the file (use h5ls)"
        write(*, *) "dims  = dims allocated for array to read"
        errval%errormsg = 'error: fdims > dims'
        errval%errbool = .true.
        return
      endif
      adims(:) = (/fdims(1)/)
      arrlen=0
      do i=1,adims(1)
        ari = len(array(i))
        arrlen = max(arrlen,ari)
      enddo
      call h5tcopy_f(h5t_native_character,dtype_id,error)
      call h5tset_size_f(dtype_id,arrlen,error)
!-------------------------------------------------------------------
!     Read data set.
!-------------------------------------------------------------------
      if (h5in%debug) write(*,*) 'creating dataspace ',aname, &
                                 ' with dims ',fdims(1)
      call h5dread_f(dset_id,dtype_id,array(1:arrlen),fdims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        call h5dclose_f(dset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset.
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_ch_vec_ex
!-----------------------------------------------------------------------
!     subprogram 89. read_int_ex
!     Read integer scalar if the object exists, but leave the input
!       unchanged if it doesn't
!-----------------------------------------------------------------------
      subroutine read_int_ex(fid,aname,val,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      integer(i4), intent(inout) :: val
      type(hdf5ErrorType), intent(inout) :: errval
    
      if (obj_exists(fid,aname,errval)) &
        call read_int(fid,aname,val,h5in,errval)

      return
      end subroutine read_int_ex
!-----------------------------------------------------------------------
!     subprogram 90. read_int_1d_ex
!     Read integer vector if the object exists, but leave the input
!       unchanged if it doesn't
!     Only read the number of elements included in the file
!-----------------------------------------------------------------------
      subroutine read_int_1d_ex(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      integer(i4), dimension(:), intent(inout) :: array
      type(hdf5ErrorType), intent(inout) :: errval
      integer(i4) :: ndims
      integer(HSIZE_T) :: dims
      integer(HSIZE_T), dimension(1) :: fdims
      integer,parameter :: FAIL=-1
      integer :: error
      integer(HID_T) :: dset_id
!-----------------------------------------------------------------------
!     Check that the object exists.
!-----------------------------------------------------------------------
      if (.NOT.obj_exists(fid,aname,errval)) return
!-----------------------------------------------------------------------
!     Read dataset as a scalar if it is not an array.
!-----------------------------------------------------------------------
      call read_ndims(fid,aname,ndims,h5in,errval)
      if (ndims < 1) then
        call read_int(fid,aname,array(1),h5in,errval)
        return
      elseif (ndims > 1) then
        write(*, *) "error: dataset is not 1d"
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading 1d i4 array: ', aname
      call h5dopen_f(fid, aname, dset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-------------------------------------------------------------------
!     Check size.
!-------------------------------------------------------------------
      dims=size(array)
      call dset_dims(dset_id,fdims,errval)
      if (errval%errBool) return
      if (fdims(1) > dims) then
        write(*, *) "error: reading ",aname," from hdf5"
        write(*, *) "fdims (", fdims(1), ")  >  dims (", dims, ")"
        write(*, *) "dims  = dims allocated for array to read"
        errval%errormsg = 'error: fdims > dims'
        errval%errbool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Read data set.
!-----------------------------------------------------------------------
      call h5dread_f(dset_id,h5t_std_i32le,array(1:fdims(1)),fdims, &
                     error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        call h5dclose_f(dset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_int_1d_ex
!-----------------------------------------------------------------------
!     subprogram 91. read_int_2d_ex
!     Read integer 2d array if the object exists, but leave the input
!       unchanged if it doesn't
!     Also read the array if it has been flattened into 1d with
!       consistent length
!-----------------------------------------------------------------------
      subroutine read_int_2d_ex(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      integer(i4), dimension(:,:), intent(inout) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: errval
      integer(i4) :: ndims
      integer(HSIZE_T), dimension(1) :: fdims1
      integer(HSIZE_T), dimension(2) :: dims,fdims
      integer(i4), allocatable, dimension(:,:) :: flat
      integer(i4), allocatable, dimension(:)  :: vec
!-----------------------------------------------------------------------
!     Check that the object exists.
!-----------------------------------------------------------------------
      if (.NOT.obj_exists(fid,aname,errval)) return
!-----------------------------------------------------------------------
!     Fit the dataset to the specified dimensions if possible.
!-----------------------------------------------------------------------
      dims(1)=size(array,1)
      dims(2)=size(array,2)
      call read_ndims(fid,aname,ndims,h5in,errval)
      if (ndims.eq.2) then
        call read_dims(fid,aname,fdims,h5in,errval)
        if ((fdims(1).le.dims(1)).and.(fdims(2).le.dims(2))) then
          call read_int_2d(fid,aname,array(1:fdims(1),1:fdims(2)),h5in,errval)
        elseif ((fdims(1).eq.dims(2)).and.(fdims(2).eq.dims(1))) then
          allocate(flat(dims(2),dims(1)))
          call read_int_2d(fid,aname,flat,h5in,errval)
          array=transpose(flat)
          deallocate(flat)
        elseif ((fdims(1).eq.1).and.(fdims(2).eq.dims(1)*dims(2))) then
          allocate(flat(1,dims(1)*dims(2)))
          call read_int_2d(fid,aname,flat,h5in,errval)
          array=reshape(flat,(/dims(1),dims(2)/))
          deallocate(flat)
        elseif ((fdims(1).eq.dims(1)*dims(2)).and.(fdims(2).eq.1)) then
          allocate(flat(dims(1)*dims(2),1))
          call read_int_2d(fid,aname,flat,h5in,errval)
          array=reshape(flat,(/dims(1),dims(2)/))
          deallocate(flat)
        else
          write(*, *) "error: reading ",aname," from hdf5"
          write(*, *) "fdims (", fdims(1), ",", fdims(2), ")"
          write(*, *) "inconsistent with  dims (", dims(1), &
                      ",", dims(2), ")"
          write(*, *) "fdims = dims in the file (use h5ls)"
          write(*, *) "dims  = dims specified for array to read"
          errval%errormsg = 'error: fdims inconsistent with dims'
          errval%errbool = .true.
        endif
      elseif (ndims.eq.1) then
        call read_dims(fid,aname,fdims1,h5in,errval)
        if (fdims1(1).eq.dims(1)*dims(2)) then
          allocate(vec(fdims1(1)))
          call read_int_1d(fid,aname,vec,h5in,errval)
          array=reshape(vec,(/dims(1),dims(2)/))
          deallocate(vec)
        else
          write(*, *) "error: reading ",aname," from hdf5"
          write(*, *) "fdims (", fdims1(1), ")"
          write(*, *) "inconsistent with  dims (", dims(1), &
                      ",", dims(2), ")"
          write(*, *) "fdims = dims in the file (use h5ls)"
          write(*, *) "dims  = dims specified for array to read"
          errval%errormsg = 'error: fdims inconsistent with dims'
          errval%errbool = .true.
        endif
      else
        write(*, *) "error: reading ",aname," from hdf5"
        write(*, *) "ndims (", ndims, ") inconsistent"
        write(*, *) "with 2d array"
        write(*, *) "ndims = #dims in the file (use h5ls)"
        errval%errormsg = 'error: ndims /= 2'
        errval%errbool = .true.
      endif
      return
      end subroutine read_int_2d_ex
!-----------------------------------------------------------------------
!     subprogram 92. read_intl_ex
!     Read integer scalar if the object exists, but leave the input
!       unchanged if it doesn't
!-----------------------------------------------------------------------
      subroutine read_intl_ex(fid,aname,val,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      integer(i8), intent(inout) :: val
      type(hdf5ErrorType), intent(inout) :: errval
    
      if (obj_exists(fid,aname,errval)) &
        call read_intl(fid,aname,val,h5in,errval)

      return
      end subroutine read_intl_ex
!-----------------------------------------------------------------------
!     subprogram 93. read_intl_1d_ex
!     Read integer vector if the object exists, but leave the input
!       unchanged if it doesn't
!     Only read the number of elements included in the file
!-----------------------------------------------------------------------
      subroutine read_intl_1d_ex(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      integer(i8), dimension(:), intent(inout) :: array
      type(hdf5ErrorType), intent(inout) :: errval
      integer(i4) :: ndims
      integer(HSIZE_T) :: dims
      integer(HSIZE_T), dimension(1) :: fdims
      integer,parameter :: FAIL=-1
      integer :: error
      integer(HID_T) :: dset_id
    
      integer(i4), dimension(:), allocatable :: intarray
!-----------------------------------------------------------------------
!     Check that the object exists.
!-----------------------------------------------------------------------
      if (.NOT.obj_exists(fid,aname,errval)) return
!-----------------------------------------------------------------------
!     Read dataset as a scalar if it is not an array.
!-----------------------------------------------------------------------
      call read_ndims(fid,aname,ndims,h5in,errval)
      if (ndims < 1) then
        call read_intl(fid,aname,array(1),h5in,errval)
        return
      elseif (ndims > 1) then
        write(*, *) "error: dataset is not 1d"
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if (h5in%verbose) WRITE(*,*) ' Reading 1d i8 array: ', aname
      call h5dopen_f(fid, aname, dset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-------------------------------------------------------------------
!     Check size.
!-------------------------------------------------------------------
      dims=size(array)
      call dset_dims(dset_id,fdims,errval)
      allocate(intarray(fdims(1)))
!      intarray=array(1:fdims(1))
      if (errval%errBool) return
      if (fdims(1) > dims) then
        write(*, *) "error: reading ",aname," from hdf5"
        write(*, *) "fdims (", fdims(1), ")  >  dims (", dims, ")"
        write(*, *) "fdims = dims in the file (use h5ls)"
        write(*, *) "dims  = dims allocated for array to read"
        errval%errormsg = 'error: fdims > dims'
        errval%errbool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Read data set.
!-----------------------------------------------------------------------
      call h5dread_f(dset_id,h5t_std_i64le,intarray,fdims,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        call h5dclose_f(dset_id,error)
        return
      endif
      array(1:fdims(1)) = intarray
      deallocate(intarray)
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_intl_1d_ex
!-----------------------------------------------------------------------
!     subprogram 94. read_rl_ex
!     Read floating point scalar if the object exists, but leave the
!       input unchanged if it doesn't
!-----------------------------------------------------------------------
      subroutine read_rl_ex(fid,aname,val,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      real(r8), intent(inout) :: val
      type(hdf5ErrorType), intent(inout) :: errval
    
      if (obj_exists(fid,aname,errval)) &
        call read_rl(fid,aname,val,h5in,errval)

      return
      end subroutine read_rl_ex
!-----------------------------------------------------------------------
!     subprogram 95. read_rl_1d_ex
!     Read floating point vector if the object exists, but leave the
!       input unchanged if it doesn't
!     Only read the number of elements included in the file
!-----------------------------------------------------------------------
      subroutine read_rl_1d_ex(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      real(r8), dimension(:), intent(inout) :: array
      type(hdf5ErrorType), intent(inout) :: errval
      integer(i4) :: ndims
      integer(HSIZE_T) :: dims
      integer(HSIZE_T), dimension(1) :: fdims
      integer,parameter :: FAIL=-1
      integer :: error
    
      integer(HID_T) :: dset_id
!-----------------------------------------------------------------------
!     Check that the object exists.
!-----------------------------------------------------------------------
      if (.NOT.obj_exists(fid,aname,errval)) return
!-----------------------------------------------------------------------
!     Read dataset as a scalar if it is not an array.
!-----------------------------------------------------------------------
      call read_ndims(fid,aname,ndims,h5in,errval)
      if (ndims < 1) then
        call read_rl(fid,aname,array(1),h5in,errval)
        return
      elseif (ndims > 1) then
        write(*, *) "error: dataset is not 1d"
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Open the dataset specified by aname
!-----------------------------------------------------------------------
      if(h5in%verbose) WRITE(*,*) ' Reading 1d r8 array: ', aname
      call h5dopen_f(fid, aname, dset_id, error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Find data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-------------------------------------------------------------------
!     Check size.
!-------------------------------------------------------------------
      dims=size(array)
      call dset_dims(dset_id,fdims,errval)
      if (errval%errBool) return
      if (fdims(1) > dims) then
        write(*, *) "error: reading ",aname," from hdf5"
        write(*, *) "fdims (", fdims(1), ")  >  dims (", dims, ")"
        write(*, *) "fdims = dims in the file (use h5ls)"
        write(*, *) "dims  = dims allocated for array to read"
        errval%errormsg = 'error: fdims > dims'
        errval%errbool = .true.
        return
      endif
!-----------------------------------------------------------------------
!     Read data set.
!-----------------------------------------------------------------------
      call h5dread_f(dset_id,h5t_ieee_f64le,array(1:fdims(1)),fdims, &
                     error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Reading data set failed for '//aname
        errval%errBool = .true.
        call h5dclose_f(dset_id,error)
        return
      endif
!-----------------------------------------------------------------------
!     Terminate access to the dataset
!-----------------------------------------------------------------------
      call h5dclose_f(dset_id,error)
      if (error==FAIL) then
        errval%errorMsg = 'ERROR: Close data set failed for '//aname
        errval%errBool = .true.
        return
      endif
!-----------------------------------------------------------------------
      errval%errBool = .false.
      return
      end subroutine read_rl_1d_ex
!-----------------------------------------------------------------------
!     subprogram 96. read_rl_2d_ex
!     Read floating point 2d array if the object exists, but leave the
!       input unchanged if it doesn't
!     Also read the array if it has been flattened into 1d with
!       consistent length
!-----------------------------------------------------------------------
      subroutine read_rl_2d_ex(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      real(r8), dimension(:,:), intent(inout) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: errval
      integer(i4) :: ndims
      integer(HSIZE_T), dimension(1) :: fdims1
      integer(HSIZE_T), dimension(2) :: dims,fdims
      real(r8), allocatable, dimension(:,:) :: flat
      real(r8), allocatable, dimension(:)  :: vec
!-----------------------------------------------------------------------
!     Check that the object exists.
!-----------------------------------------------------------------------
      if (.NOT.obj_exists(fid,aname,errval)) return
!-----------------------------------------------------------------------
!     Fit the dataset to the specified dimensions if possible.
!-----------------------------------------------------------------------
      dims(1)=size(array,1)
      dims(2)=size(array,2)
      call read_ndims(fid,aname,ndims,h5in,errval)
      if (ndims.eq.2) then
        call read_dims(fid,aname,fdims,h5in,errval)
        if ((fdims(1).le.dims(1)).and.(fdims(2).le.dims(2))) then
          call read_rl_2d(fid,aname,array(1:fdims(1),1:fdims(2)),h5in,errval)
        elseif ((fdims(1).eq.dims(2)).and.(fdims(2).eq.dims(1))) then
          allocate(flat(dims(2),dims(1)))
          call read_rl_2d(fid,aname,flat,h5in,errval)
          array=transpose(flat)
          deallocate(flat)
        elseif ((fdims(1).eq.1).and.(fdims(2).eq.dims(1)*dims(2))) then
          allocate(flat(1,dims(1)*dims(2)))
          call read_rl_2d(fid,aname,flat,h5in,errval)
          array=reshape(flat,(/dims(1),dims(2)/))
          deallocate(flat)
        elseif ((fdims(1).eq.dims(1)*dims(2)).and.(fdims(2).eq.1)) then
          allocate(flat(dims(1)*dims(2),1))
          call read_rl_2d(fid,aname,flat,h5in,errval)
          array=reshape(flat,(/dims(1),dims(2)/))
          deallocate(flat)
        else
          write(*, *) "error: reading ",aname," from hdf5"
          write(*, *) "fdims (", fdims(1), ",", fdims(2), ")"
          write(*, *) "inconsistent with  dims (", dims(1), &
                      ",", dims(2), ")"
          write(*, *) "fdims = dims in the file (use h5ls)"
          write(*, *) "dims  = dims specified for array to read"
          errval%errormsg = 'error: fdims inconsistent with dims'
          errval%errbool = .true.
        endif
      elseif (ndims.eq.1) then
        call read_dims(fid,aname,fdims1,h5in,errval)
        if (fdims1(1).eq.dims(1)*dims(2)) then
          allocate(vec(fdims1(1)))
          call read_rl_1d(fid,aname,vec,h5in,errval)
          array=reshape(vec,(/dims(1),dims(2)/))
          deallocate(vec)
        else
          write(*, *) "error: reading ",aname," from hdf5"
          write(*, *) "fdims (", fdims1(1), ")"
          write(*, *) "inconsistent with  dims (", dims(1), &
                      ",", dims(2), ")"
          write(*, *) "fdims = dims in the file (use h5ls)"
          write(*, *) "dims  = dims specified for array to read"
          errval%errormsg = 'error: fdims inconsistent with dims'
          errval%errbool = .true.
        endif
      else
        write(*, *) "error: reading ",aname," from hdf5"
        write(*, *) "ndims (", ndims, ") inconsistent"
        write(*, *) "with 2d array"
        write(*, *) "ndims = #dims in the file (use h5ls)"
        errval%errormsg = 'error: ndims /= 2'
        errval%errbool = .true.
      endif
      return
      end subroutine read_rl_2d_ex
!-----------------------------------------------------------------------
!     subprogram 97. read_rl_3d_ex
!     Read floating point 3d array if the object exists, but leave the
!       input unchanged if it doesn't
!-----------------------------------------------------------------------
      subroutine read_rl_3d_ex(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      real(r8), dimension(:,:,:), intent(inout) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: errval
    
      if (obj_exists(fid,aname,errval)) &
        call read_rl_3d(fid,aname,array,h5in,errval)

      return
      end subroutine read_rl_3d_ex
!-----------------------------------------------------------------------
!     subprogram 98. read_rl_4d_ex
!     Read floating point 4d array if the object exists, but leave the
!       input unchanged if it doesn't
!-----------------------------------------------------------------------
      subroutine read_rl_4d_ex(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      real(r8), dimension(:,:,:,:), intent(inout) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: errval
    
      if (obj_exists(fid,aname,errval)) &
        call read_rl_4d(fid,aname,array,h5in,errval)

      end subroutine read_rl_4d_ex
!-----------------------------------------------------------------------
!     subprogram 99. read_5d_ex
!     Read floating point 5d array if the object exists, but leave the
!       input unchanged if it doesn't
!-----------------------------------------------------------------------
      subroutine read_rl_5d_ex(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      real(r8), dimension(:,:,:,:,:), intent(inout) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: errval
    
      if (obj_exists(fid,aname,errval)) &
        call read_rl_5d(fid,aname,array,h5in,errval)

      return
      end subroutine read_rl_5d_ex
!-----------------------------------------------------------------------
!     subprogram 100. read_log_ex
!     Read logical as an integer if the object exists, but leave the
!       input unchanged if it doesn't
!-----------------------------------------------------------------------
      subroutine read_log_ex(fid,aname,bool,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      type(hdf5InOpts), intent(in) :: h5in
      logical, intent(inout) :: bool
      type(hdf5ErrorType), intent(inout) :: errval
      integer(i4) :: val
!-----------------------------------------------------------------------
!     Check that the object exists.
!-----------------------------------------------------------------------
      if (.NOT.obj_exists(fid,aname,errval)) return
!-----------------------------------------------------------------------
!     Read data as an integer and convert to logical
!-----------------------------------------------------------------------
      call read_int(fid,aname,val,h5in,errval)
      if (errval%errBool) return
      bool=.false.
      if (val.ne.0) bool=.true.
      return
      end subroutine read_log_ex
!-----------------------------------------------------------------------
!     subprogram 101. read_int_2d_sq
!     Read integer 2d array if the object exists, but leave the input
!       unchanged if it doesn't
!     If the dataset is smaller than available space, ensure that it is
!       read into the subdominant square
!     If it is larger, only read in the subdominant square
!-----------------------------------------------------------------------
      subroutine read_int_2d_sq(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      integer(i4), dimension(:,:), intent(inout) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: errval
      integer(i4) :: ndims
      integer(HSIZE_T) :: dims1
      integer(HSIZE_T), dimension(1) :: fdims1,sqdims1
      integer(HSIZE_T), dimension(2) :: dims,fdims,sqdims
      integer(i4), allocatable, dimension(:,:) :: mat
      integer(i4), allocatable, dimension(:)  :: vec
!-----------------------------------------------------------------------
!     Check that the object exists.
!-----------------------------------------------------------------------
      if (.NOT.obj_exists(fid,aname,errval)) return
!-----------------------------------------------------------------------
!     Fit the dataset into a square if possible.
!-----------------------------------------------------------------------
      dims(1)=size(array,1)
      dims(2)=size(array,2)
      dims1=dims(1)*dims(2)
      call read_ndims(fid,aname,ndims,h5in,errval)
      if (ndims.eq.2) then
        call read_dims(fid,aname,fdims,h5in,errval)
        sqdims=int(sqrt(real(fdims)))
        if (fdims(1).eq.fdims(2)) then
          if ((fdims(1).le.dims(1)) .and.(fdims(2).le.dims(2))) then
            call read_int_2d(fid,aname,array(1:fdims(1),1:fdims(2)),h5in,errval)
          else
            allocate(mat(fdims(1),fdims(2)))
            call read_int_2d(fid,aname,mat,h5in,errval)
            array=mat(1:dims(1),1:dims(2))
            deallocate(mat)
          endif
        elseif ((fdims(1).eq.1).and.(fdims(2).eq.sqdims(2)**2)) then
          allocate(mat(1,fdims(2)))
          call read_int_2d(fid,aname,mat,h5in,errval)
          if ((sqdims(2).le.dims(1)).and.(sqdims(2).le.dims(2))) then
            array(1:sqdims(2),1:sqdims(2))=reshape(mat,(/sqdims(2),sqdims(2)/))
          else
            array=reshape(mat(1,1:dims1),(/dims(1),dims(2)/))
          endif
          deallocate(mat)
        elseif ((fdims(2).eq.1).and.(fdims(1).eq.sqdims(1)**2)) then
          allocate(mat(fdims(1),1))
          call read_int_2d(fid,aname,mat,h5in,errval)
          if ((sqdims(1).le.dims(1)).and.(sqdims(1).le.dims(2))) then
            array(1:sqdims(1),1:sqdims(1))=reshape(mat,(/sqdims(1),sqdims(1)/))
          else
            array=reshape(mat(1:dims1,1),(/dims(1),dims(2)/))
          endif
          deallocate(mat)
        else
          write(*, *) "error: reading ",aname," from hdf5"
          write(*, *) "fdims (", fdims(1), ",", fdims(2), ")"
          write(*, *) "cannot be made congruent to square"
          write(*, *) "dims (", dims(1), ",", dims(2), ")"
          write(*, *) "fdims = dims in the file (use h5ls)"
          write(*, *) "dims  = dims specified for array to read"
          errval%errormsg = 'error: fdims inconsistent with dims'
          errval%errbool = .true.
        endif
      elseif (ndims.eq.1) then
        call read_dims(fid,aname,fdims1,h5in,errval)
        sqdims1=int(sqrt(real(fdims1)))
        if ((fdims1(1).eq.sqdims1(1)**2)) then
          allocate(vec(fdims1(1)))
          call read_int_1d(fid,aname,vec,h5in,errval)
          if ((sqdims1(1).le.dims(1)).and.(sqdims1(1).le.dims(2))) then
            array(1:sqdims1(1),1:sqdims1(1))=reshape(vec,(/sqdims1(1),sqdims1(1)/))
          else
            allocate(mat(sqdims1(1),sqdims1(1)))
            mat=reshape(vec,(/sqdims1(1),sqdims1(1)/))
            array=mat(1:dims(1),1:dims(2))
            deallocate(mat)
          endif
          deallocate(vec)
        else
          write(*, *) "error: reading ",aname," from hdf5"
          write(*, *) "fdims (", fdims1(1), ")"
          write(*, *) "cannot be made congruent to square"
          write(*, *) "dims (", dims(1), ",", dims(2), ")"
          write(*, *) "fdims = dims in the file (use h5ls)"
          write(*, *) "dims  = dims specified for array to read"
          errval%errormsg = 'error: fdims inconsistent with dims'
          errval%errbool = .true.
        endif
      else
        write(*, *) "error: reading ",aname," from hdf5"
        write(*, *) "ndims (", ndims, ") inconsistent"
        write(*, *) "with 2d array"
        write(*, *) "ndims = #dims in the file (use h5ls)"
        errval%errormsg = 'error: ndims /= 2'
        errval%errbool = .true.
      endif
      return
      end subroutine read_int_2d_sq
!-----------------------------------------------------------------------
!     subprogram 102. read_rl_2d_sq
!     Read floating point 2d array if the object exists, but leave the
!       input unchanged if it doesn't
!     If the dataset is smaller than available space, ensure that it is
!       read into the subdominant square
!     If it is larger, only read in the subdominant square
!-----------------------------------------------------------------------
      subroutine read_rl_2d_sq(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      real(r8), dimension(:,:), intent(inout) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: errval
      integer(i4) :: ndims
      integer(HSIZE_T) :: dims1
      integer(HSIZE_T), dimension(1) :: fdims1,sqdims1
      integer(HSIZE_T), dimension(2) :: dims,fdims,sqdims
      real(r8), allocatable, dimension(:,:) :: mat
      real(r8), allocatable, dimension(:)  :: vec
!-----------------------------------------------------------------------
!     Check that the object exists.
!-----------------------------------------------------------------------
      if (.NOT.obj_exists(fid,aname,errval)) return
!-----------------------------------------------------------------------
!     Fit the dataset to the specified dimensions if possible.
!-----------------------------------------------------------------------
      dims(1)=size(array,1)
      dims(2)=size(array,2)
      dims1=dims(1)*dims(2)
      call read_ndims(fid,aname,ndims,h5in,errval)
      if (ndims.eq.2) then
        call read_dims(fid,aname,fdims,h5in,errval)
        sqdims=int(sqrt(real(fdims)))
        if (fdims(1).eq.fdims(2)) then
          if ((fdims(1).le.dims(1)) .and.(fdims(2).le.dims(2))) then
            call read_rl_2d(fid,aname,array(1:fdims(1),1:fdims(2)),h5in,errval)
          else
            allocate(mat(fdims(1),fdims(2)))
            call read_rl_2d(fid,aname,mat,h5in,errval)
            array=mat(1:dims(1),1:dims(2))
            deallocate(mat)
          endif
        elseif ((fdims(1).eq.1).and.(fdims(2).eq.sqdims(2)**2)) then
          allocate(mat(1,fdims(2)))
          call read_rl_2d(fid,aname,mat,h5in,errval)
          if ((sqdims(2).le.dims(1)).and.(sqdims(2).le.dims(2))) then
            array(1:sqdims(2),1:sqdims(2))=reshape(mat,(/sqdims(2),sqdims(2)/))
          else
            array=reshape(mat(1,1:dims1),(/dims(1),dims(2)/))
          endif
          deallocate(mat)
        elseif ((fdims(2).eq.1).and.(fdims(1).eq.sqdims(1)**2)) then
          allocate(mat(fdims(1),1))
          call read_rl_2d(fid,aname,mat,h5in,errval)
          if ((sqdims(1).le.dims(1)).and.(sqdims(1).le.dims(2))) then
            array(1:sqdims(1),1:sqdims(1))=reshape(mat,(/sqdims(1),sqdims(1)/))
          else
            array=reshape(mat(1:dims1,1),(/dims(1),dims(2)/))
          endif
          deallocate(mat)
        else
          write(*, *) "error: reading ",aname," from hdf5"
          write(*, *) "fdims (", fdims(1), ",", fdims(2), ")"
          write(*, *) "cannot be made congruent to square"
          write(*, *) "dims (", dims(1), ",", dims(2), ")"
          write(*, *) "fdims = dims in the file (use h5ls)"
          write(*, *) "dims  = dims specified for array to read"
          errval%errormsg = 'error: fdims inconsistent with dims'
          errval%errbool = .true.
        endif
      elseif (ndims.eq.1) then
        call read_dims(fid,aname,fdims1,h5in,errval)
        sqdims1=int(sqrt(real(fdims1)))
        if ((fdims1(1).eq.sqdims1(1)**2)) then
          allocate(vec(fdims1(1)))
          call read_rl_1d(fid,aname,vec,h5in,errval)
          if ((sqdims1(1).le.dims(1)).and.(sqdims1(1).le.dims(2))) then
            array(1:sqdims1(1),1:sqdims1(1))=reshape(vec,(/sqdims1(1),sqdims1(1)/))
          else
            allocate(mat(sqdims1(1),sqdims1(1)))
            mat=reshape(vec,(/sqdims1(1),sqdims1(1)/))
            array=mat(1:dims(1),1:dims(2))
            deallocate(mat)
          endif
          deallocate(vec)
        else
          write(*, *) "error: reading ",aname," from hdf5"
          write(*, *) "fdims (", fdims1(1), ")"
          write(*, *) "cannot be fit into square subset of"
          write(*, *) "dims (", dims(1), ",", dims(2), ")"
          write(*, *) "fdims = dims in the file (use h5ls)"
          write(*, *) "dims  = dims specified for array to read"
          errval%errormsg = 'error: fdims inconsistent with dims'
          errval%errbool = .true.
        endif
      else
        write(*, *) "error: reading ",aname," from hdf5"
        write(*, *) "ndims (", ndims, ") inconsistent"
        write(*, *) "with 2d array"
        write(*, *) "ndims = #dims in the file (use h5ls)"
        errval%errormsg = 'error: ndims /= 2'
        errval%errbool = .true.
      endif
      return
      end subroutine read_rl_2d_sq
!-----------------------------------------------------------------------
!     subprogram 103. read_rl_1d_sq
!     Read floating point 2d object into a 1d array if it exists, but
!       leave the input unchanged if it doesn't
!     If the dataset is smaller than available space, ensure that it is
!       read into the subdominant square
!     If it is larger, only read in the subdominant square
!-----------------------------------------------------------------------
      subroutine read_rl_1d_sq(fid,aname,array,h5in,errval)
      integer(HID_T), intent(in) :: fid
      character*(*), intent(in) :: aname
      real(r8), dimension(:), intent(inout) :: array
      type(hdf5InOpts), intent(in) :: h5in
      type(hdf5ErrorType), intent(inout) :: errval
      integer(i4) :: ndims
      integer(HSIZE_T) :: dims1,dims
      integer(HSIZE_T), dimension(1) :: fdims1,sqdims1
      integer(HSIZE_T), dimension(2) :: fdims,sqdims
      real(r8), allocatable, dimension(:,:) :: mat,arrsq
      real(r8), allocatable, dimension(:)  :: vec
!-----------------------------------------------------------------------
!     Check that the object exists.
!-----------------------------------------------------------------------
      if (.NOT.obj_exists(fid,aname,errval)) return
!-----------------------------------------------------------------------
!     Fit the dataset to the specified dimensions if possible.
!-----------------------------------------------------------------------
      dims1=size(array,1)
      dims=int(sqrt(real(dims1)))
      allocate(arrsq(dims,dims))
      call read_ndims(fid,aname,ndims,h5in,errval)
      if (ndims.eq.2) then
        call read_dims(fid,aname,fdims,h5in,errval)
        sqdims=int(sqrt(real(fdims)))
        if (fdims(1).eq.fdims(2)) then
          if (fdims(1).le.dims) then
            call read_rl_2d(fid,aname,arrsq(1:fdims(1),1:fdims(2)),h5in,errval)
            array=reshape(arrsq,(/dims1/))
          else
            allocate(mat(fdims(1),fdims(2)))
            call read_rl_2d(fid,aname,mat,h5in,errval)
            array=reshape(mat(1:dims,1:dims),(/dims1/))
            deallocate(mat)
          endif
        elseif ((fdims(1).eq.1).and.(fdims(2).eq.sqdims(2)**2)) then
          allocate(mat(1,fdims(2)))
          call read_rl_2d(fid,aname,mat,h5in,errval)
          if (sqdims(2).le.dims) then
            arrsq(1:sqdims(2),1:sqdims(2))=reshape(mat,(/sqdims(2),sqdims(2)/))
          else
            arrsq=reshape(mat(1,1:dims1),(/dims,dims/))
          endif
          array=reshape(arrsq,(/dims1/))
          deallocate(mat)
        elseif ((fdims(2).eq.1).and.(fdims(1).eq.sqdims(1)**2)) then
          allocate(mat(fdims(1),1))
          call read_rl_2d(fid,aname,mat,h5in,errval)
          if (sqdims(1).le.dims) then
            arrsq(1:sqdims(1),1:sqdims(1))=reshape(mat,(/sqdims(1),sqdims(1)/))
          else
            arrsq=reshape(mat(1:dims1,1),(/dims,dims/))
          endif
          array=reshape(arrsq,(/dims1/))
          deallocate(mat)
        else
          write(*, *) "error: reading ",aname," from hdf5"
          write(*, *) "fdims (", fdims(1), ",", fdims(2), ")"
          write(*, *) "cannot be made congruent to square"
          write(*, *) "dims (", dims, ",", dims, ")"
          write(*, *) "fdims = dims in the file (use h5ls)"
          write(*, *) "dims  = dims specified for array to read"
          errval%errormsg = 'error: fdims inconsistent with dims'
          errval%errbool = .true.
        endif
      elseif (ndims.eq.1) then
        call read_dims(fid,aname,fdims1,h5in,errval)
        sqdims1=int(sqrt(real(fdims1)))
        if ((fdims1(1).eq.sqdims1(1)**2)) then
          allocate(vec(fdims1(1)))
          call read_rl_1d(fid,aname,vec,h5in,errval)
          if (sqdims1(1).le.dims) then
            arrsq(1:sqdims1(1),1:sqdims1(1))=reshape(vec,(/sqdims1(1),sqdims1(1)/))
          else
            allocate(mat(sqdims1(1),sqdims1(1)))
            mat=reshape(vec,(/sqdims1(1),sqdims1(1)/))
            arrsq=mat(1:dims,1:dims)
            arrsq=reshape(vec(1:dims1),(/dims,dims/))
            deallocate(mat)
          endif
          array=reshape(arrsq,(/dims1/))
          deallocate(vec)
        else
          write(*, *) "error: reading ",aname," from hdf5"
          write(*, *) "fdims (", fdims1(1), ")"
          write(*, *) "cannot be fit into square subset of"
          write(*, *) "dims (", dims, ",", dims, ")"
          write(*, *) "fdims = dims in the file (use h5ls)"
          write(*, *) "dims  = dims specified for array to read"
          errval%errormsg = 'error: fdims inconsistent with dims'
          errval%errbool = .true.
        endif
      else
        write(*, *) "error: reading ",aname," from hdf5"
        write(*, *) "ndims (", ndims, ") inconsistent"
        write(*, *) "with 2d array"
        write(*, *) "ndims = #dims in the file (use h5ls)"
        errval%errormsg = 'error: ndims /= 2'
        errval%errbool = .true.
      endif
      deallocate(arrsq)
      return
      end subroutine read_rl_1d_sq
!-----------------------------------------------------------------------
!     close module.
!-----------------------------------------------------------------------
      end module hdf5_api
