!-------------------------------------------------------------------------------
!< this file is to be used on serial machines.  it contains two parts.
!  the first is a module of parameters and the second is a set of
!  dummy routines to replace mpi routines.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!> module for interface to F77 MPI include file
!  defines all MPI variables via parameter statements
!  use this module to define machine-specific MPI datatypes
!-------------------------------------------------------------------------------
MODULE mpi_efit
  IMPLICIT NONE
END MODULE mpi_efit

!-------------------------------------------------------------------------------
! stubs for MPI calls - use on serial machine to replace real MPI
! library so NIMROD will still compile
! except for mpi_comm_rank, mpi_comm_size and mpi_allreduce
! these are all no-operation rouines
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_testany(nrecv,recv_request,irecv,flag,status,ierror)
  IMPLICIT NONE

  REAL :: recv_request(nrecv)
  REAL :: status(1)
  INTEGER :: nrecv,irecv,ierror
  LOGICAL :: flag

  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE comm_create_(nsend,procsend,comm,nrecv,plan)
  IMPLICIT NONE

  INTEGER :: nsend,nrecv,procsend,comm
  DOUBLE PRECISION :: plan

  nrecv=nsend
  plan=REAL(nsend)
  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE comm_destroy_(plan)
  IMPLICIT NONE

  DOUBLE PRECISION :: plan

  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE comm_do_(plan,sendbuf,n,recvbuf)
  IMPLICIT NONE

  INTEGER :: n,sendbuf(*),recvbuf(*)
  DOUBLE PRECISION :: plan
  INTEGER :: i

  DO i=1,nint(plan)
    recvbuf(i)=sendbuf(i)
  ENDDO
  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_init(ierror)
  IMPLICIT NONE

  INTEGER :: ierror,nprocs

  nprocs=1
  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_finalize(ierror)
  IMPLICIT NONE

  INTEGER :: ierror

  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_init_thread(tlvlin,tlvlout,ierror)
  IMPLICIT NONE

  INTEGER :: tlvlin,tlvlout,ierror

  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_comm_rank(comm,rank,ierror)
  IMPLICIT NONE

  INTEGER :: comm,rank,ierror

  rank = 0
  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_comm_size(comm,size,ierror)
  IMPLICIT NONE

  INTEGER :: comm,size,ierror

  size = 1
  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_bcast(buffer,count,datatype,root,comm,ierror)
  IMPLICIT NONE

  INTEGER :: buffer,count,datatype,root,comm,ierror

  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_barrier(comm,ierror)
  IMPLICIT NONE

  INTEGER :: comm,ierror

  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_allreduce(sendbuf,recvbuf,count,datatype,op,comm,ierror)
  IMPLICIT NONE

  INTEGER :: count,sendbuf(count),recvbuf(count),datatype,op,comm,ierror
  INTEGER :: i

  DO i=1,count
    recvbuf(i)=sendbuf(i)
  ENDDO
  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_reduce(sendbuf,recvbuf,count,datatype,op,root,comm,ierror)
  IMPLICIT NONE

  INTEGER :: count,sendbuf(count),recvbuf(count),datatype,op,root,comm,ierror
  INTEGER :: i

  DO i=1,count
    recvbuf(i)=sendbuf(i)
  ENDDO
  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_reduce_scatter(sendbuf,recvbuf,recvcounts,datatype,op,comm,      &
                              ierror)
  IMPLICIT NONE

  INTEGER :: sendbuf,recvbuf,recvcounts,datatype,op,comm,ierror

  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_send(buf,count,datatype,dest,tag,comm,ierror)
  IMPLICIT NONE

  INTEGER :: buf,count,datatype,dest,tag,comm,ierror

  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_recv(buf,count,datatype,source,tag,comm,status,ierror)
  IMPLICIT NONE

  INTEGER :: buf,count,datatype,source,tag,comm,status,ierror

  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_irecv(buf,count,datatype,source,tag,comm,request,ierror)
  IMPLICIT NONE

  INTEGER :: buf,count,datatype,source,tag,comm,request,ierror

  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_isend(buf,count,datatype,source,tag,comm,request,ierror)
  IMPLICIT NONE

  INTEGER :: buf,count,datatype,source,tag,comm,request,ierror

  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_sendrecv(sendbuf,scounts,stypes,dest,sendtag,recvbuf,rcounts,    &
                        rtypes,source,recvtag,comm,status,ierror)
  IMPLICIT NONE

  INTEGER :: scounts,sendbuf(scounts),stypes,dest,sendtag,rcounts,     &
             recvbuf(rcounts),rtypes,source,recvtag,comm,status,ierror
  INTEGER :: i

  DO i=1,MIN(rcounts,scounts)
    recvbuf(i)=sendbuf(i)
  ENDDO
  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_test(request,flag,status,ierror)
  IMPLICIT NONE

  LOGICAL :: flag
  INTEGER :: request,status,ierror

  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_comm_split(datatype1,ilayer,ii,datatype2,ierror)
  IMPLICIT NONE

  INTEGER :: datatype1,datatype2,ilayer,ierror,ii

  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_wait(recv_request,status,ierror)
  IMPLICIT NONE

  REAL :: recv_reqest
  REAL :: recv_request(1)
  REAL :: status
  REAL :: statuses(1)
  INTEGER :: nrecv,ierror

  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_waitall(nrecv,recv_request,statuses,ierror)
  IMPLICIT NONE

  INTEGER :: nrecv,ierror,recv_request(nrecv),statuses(1,1)

  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_waitany(nrecv,recv_request,irecv,status,ierror)
  IMPLICIT NONE

  INTEGER :: nrecv,irecv,ierror,recv_request(nrecv),status(1)

  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_allgather(sendbuf,counts,datatypes,recvbuf,countr,displs,        &
                         datatyper,comm,ierror)
  IMPLICIT NONE

  INTEGER :: counts,sendbuf(counts),recvbuf(counts),countr,datatypes,datatyper, &
             displs,comm,ierror
  INTEGER :: i

  DO i=1,counts
    recvbuf(i)=sendbuf(i)
  ENDDO
  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_gather(sendbuf,counts,datatypes,recvbuf,countr,displs,datatyper, &
                      comm,ierror)
  IMPLICIT NONE

  INTEGER :: counts,sendbuf(counts),recvbuf(counts),countr,datatypes,datatyper, &
             displs,comm,ierror
  INTEGER :: i

  DO i=1,counts
    recvbuf(i)=sendbuf(i)
  ENDDO
  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_allgatherv(sendbuf,counts,datatypes,recvbuf,countr,displs,       &
                          datatyper,comm,ierror)
  IMPLICIT NONE

  INTEGER :: counts,sendbuf(counts),recvbuf(counts),countr,datatypes,datatyper, &
             displs,comm,ierror
  INTEGER :: i

  DO i=1,counts
    recvbuf(i)=sendbuf(i)
  ENDDO
  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_alltoallv(sendbuf,scounts,sdispls,datatypes,recvbuf,rcounts,     &
                         rdispls,datatyper,comm,ierror)
  IMPLICIT NONE

  INTEGER :: scounts,sendbuf(scounts),sdispls,rcounts,recvbuf(rcounts),rdispls, &
             datatypes,datatyper,comm,ierror
  INTEGER :: i

  DO i=1,MIN(rcounts,scounts)
    recvbuf(i)=sendbuf(i)
  ENDDO
  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_abort(comm,code,ierror)
  IMPLICIT NONE

  INTEGER :: comm,code,ierror

  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_get_address(datum,addrs,ierror)
  IMPLICIT NONE

  INTEGER :: datum,addrs,ierror

  addrs=0
  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_address(datum,addrs,ierror)
  IMPLICIT NONE

  INTEGER :: datum,addrs,ierror

  addrs=0
  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_pack(inbuf,incount,dtype,outbuf,outcount,pos,comm,ierror)
  IMPLICIT NONE

  INTEGER :: inbuf,incount,dtype,outbuf,outcount,pos,comm,ierror

  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_unpack(inbuf,insize,pos,outbuf,outcount,dtype,comm,ierror)
  IMPLICIT NONE

  INTEGER :: inbuf,insize,dtype,outbuf,outcount,pos,comm,ierror

  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_pack_size(incount,dtype,comm,insize,ierror)
  IMPLICIT NONE

  INTEGER :: incount,dtype,comm,insize,ierror

  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_type_struct(incount,blength,displ,atype,newtype,ierror)
  IMPLICIT NONE

  INTEGER :: incount,blength,displ,atype,newtype,ierror

  newtype=incount*blength
  RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!* serial mpi stub routine
!-------------------------------------------------------------------------------
SUBROUTINE mpi_type_commit(dtype,ierror)
  IMPLICIT NONE

  INTEGER :: dtype,ierror

  RETURN
END SUBROUTINE

