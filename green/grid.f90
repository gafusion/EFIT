      module grid

      use machine, only: nw,nh,nwnh
      implicit none

      integer*4 :: igrid
      real*8 :: rleft,rright,zbotto,ztop
      real*8 :: dr,dz
      real*8, dimension(:), allocatable :: rgrid
      real*8, dimension(:), allocatable :: zgrid
      real*8, dimension(:,:), allocatable :: brgridfc,bzgridfc

      contains
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          grid computes the green's functions at (r,z)            **
!**          due to plasma currents and f coils.                     **
!**                                                                  **
!**********************************************************************
      subroutine efund_grid
      use machine, only: nfcoil,nfsum,nw,nh,nwnh
      use utils, only: pi,tmu,psical,flux
      use fcoil
      use nio
      implicit none
      integer*4 i,ih,ii,iw,j,k,kk,ni,nj
      real*8 ab,cmut,dhc,drc,dwc,dzc,fridpc,rtmp,rsum,z,z1,z2,zdif
      real*8, dimension(:,:), allocatable :: gridfc,gridpc,ggridfc
      real*8 taf(nfcoil),taf2(nfcoil)
      integer*4, parameter :: isplit=17
      real*8, parameter :: aaa=0.0,tole=1.0e-10
!
      if (.not. allocated(gridfc)) then
        allocate(gridfc(nwnh,nfcoil))
        gridfc = 0.0
      endif
      if (.not. allocated(gridpc)) then
        allocate(gridpc(nwnh,nw))
        gridpc = 0.0
      endif
      if (.not. allocated(ggridfc)) then
        allocate(ggridfc(nwnh,nfsum))
        ggridfc = 0.0
      endif
!
      if(igrid.eq.0) return
!----------------------------------------------------------------------
!--  compute the green's functions at (r,z) due to f coils           --
!----------------------------------------------------------------------
      do ii=1,nfcoil
         taf(ii) = tan(af(ii)*pi/180.0)
         taf2(ii) = tan(af2(ii)*pi/180.0)
         do ni=1,nw
            do nj=1,nh
               kk = (ni-1)*nh + nj
               rsum = 0.
               dwc = wf(ii)/isplit
               dhc = hf(ii)/isplit
               if (af(ii) .eq. 0. .and. af2(ii) .eq. 0.) then
                  z = zf(ii)-.5*hf(ii)+.5*dhc
                  ab = rf(ii)-.5*wf(ii)+.5*dwc
                  do iw = 1,isplit
                     drc = ab+(iw-1)*dwc+iw*tole
                     do ih = 1,isplit
                        dzc = z+(ih-1)*dhc
                        rtmp = psical(drc,rgrid(ni),zgrid(nj)-dzc)
                        rsum = rsum+rtmp
                     enddo
                  enddo
               elseif (af(ii) .ne. 0.) then
                  z1 = zf(ii)-taf(ii)*(wf(ii)-dwc)/2.-.5*hf(ii)+.5*dhc
                  ab = rf(ii)-.5*wf(ii)+.5*dwc
                  do iw = 1,isplit
                     drc = ab+(iw-1)*dwc+iw*tole
                     z2 = z1+(iw-1)*taf(ii)*dwc
                     do ih = 1,isplit
                        dzc = z2+(ih-1)*dhc
                        rtmp = psical(drc,rgrid(ni),zgrid(nj)-dzc)
                        rsum = rsum+rtmp
                     enddo
                  enddo
               else
                  do ih = 1,isplit
                     dzc = zf(ii)-.5*hf(ii)+.5*dhc+dhc*(ih-1)
                     do iw = 1,isplit
                        drc = rf(ii)-.5*wf(ii)-.5*hf(ii)/taf2(ii) &
                              +.5*dwc+.5*dhc/taf2(ii) &
                              +dhc/taf2(ii)*(ih-1)+dwc*(iw-1)
                        rtmp = psical(drc,rgrid(ni),zgrid(nj)-dzc)
                        rsum = rsum+rtmp
                     enddo
                  enddo
               endif
               cmut = rsum*2.e-07/(isplit*isplit)
               gridfc(kk,ii) = cmut
            enddo
         enddo
      enddo
!----------------------------------------------------------------------
!--  compute the green's functions at (r,z) due to itself            --
!----------------------------------------------------------------------
      do i=1,nw
         do j=1,nh
            kk=(i-1)*nh+j
            do ni=1,nw
               if ((j.gt.1).or.(i.ne.ni)) then
                  zdif=(j-1)*dz
                  gridpc(kk,ni)=psical(rgrid(i),rgrid(ni),zdif)*tmu
               else
                  call flux(rgrid(ni),aaa,dr,dz,aaa,aaa,rgrid(ni),aaa, &
                            dr,dz,aaa,aaa,fridpc)
                  gridpc(kk,ni)=fridpc*0.5/pi
               endif
            enddo 
         enddo
      enddo
!----------------------------------------------------------------------
!--  store green's function table                                    --
!----------------------------------------------------------------------
      ggridfc=0.0
      do i=1,nfcoil
         ggridfc(:,fcid(i))=ggridfc(:,fcid(i))+fcturn(i)*gridfc(:,i)
      enddo
!
      print*,'file name : ','ec'//trim(ch1)//trim(ch2)//'.ddd' 
!
!vasorg      open(unit=ncontr,status='unknown',file='econto.dat', &
      open(unit=ncontr,status='unknown',file='ec'//trim(ch1)// &
           trim(ch2)//'.ddd',form='unformatted')
      write (ncontr) nw,nh
      write (ncontr) rgrid,zgrid
      write (ncontr) ggridfc
      write (ncontr) gridpc
!vas just for testing
!      open(35,file='test-ec1.dat',status='new')
!      write (35,*) nw,nh
!      write (35,1009) rgrid,zgrid
!      close(35)
!      open(35,file='test-ec2.dat',status='new')
!      write (35,1009) ggridfc
!      close(35)
!      open(35,file='test-ec3.dat',status='new')
!      write (35,1009) gridpc
!      close(35)
!1009  format(3(1x,e14.8))
      close(unit=ncontr)
!
      if (allocated(gridfc)) deallocate(gridfc)
      if (allocated(ggridfc)) deallocate(ggridfc)
      if (allocated(gridpc)) deallocate(gridpc)
!
      return
      end subroutine efund_grid

      end module grid
