      module machine

      implicit none
      public
      
      integer*4 nfcoil,nsilop,magpri,nrogow,necoil,nesum, &
                nfsum,nvsum,nvesel,nacoil

      integer*4 nw,nh,nwnh
      character*10 device


      integer*4 nnece,nnecein,neceo,mpress,libim,nmselp,nmsels,&
                nnnte,ngam_vars,ngam_u,ngam_w,nlimit,nlimbd,nangle,&
                ntangle,nfbcoil,mccoil,micoil,ndata,nwwcur,&
                nffcur,nppcur,nercur,ntime,ndim,kxiter,mqwant,& 
                mbdry,mbdry1,nxtram,nxtlim,nco2v,nco2r,modef,&
                modep,modew,nfourier

      end module machine
