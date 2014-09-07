c -*- Mode: Fortran; -*-
c----------------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2012
c     King Abdullah University of Science and Technology
c     All Rights Reserved
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c     recycling inflow velocity, from iprocxre to iprocx=1.
c     mirror method
c----------------------------------------------------------------------
      subroutine exchangeinflow_mirror(iprocxre,uinlett)
      use mesh_parms
      use mesh_common
      use InflowVelBC  
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      integer,parameter:: nfields=3
      integer:: iprocxre
      double precision:: uinlett(iylo:iyhi,izlo:izhi,3)
c      
#ifdef PARALLEL
      double precision, dimension
     &  (IYLO:IYHI,IZLO:IZHI,nfields) 
     &     :: zbuffer_send, zbuffer_recv
      integer:: XBUFFSIZE, YBUFFSIZE, ZBUFFSIZE
      integer:: ii,jj,mm,kk,l,idest
      integer:: iprocnum, iproc_re
      integer:: iipy, iipz, ipc_to, ipc_from
      integer:: msg_id_send_re
      integer:: msg_id_recv_re
c
c     Z DIRECTION COMMUNICATION, done for 3D simulations only
      ZBUFFSIZE = nfields*(IYHI-IYLO+1)*(IZHI-IZLO+1)    
      if(iprocx.eq. iprocxre) then
        iproc_re = iprocnum(1, yprocs+1-iprocy, iprocz)
        ipc_to   =  my_iproc(iproc_re)
        do l = 1,nfields
          do mm = 1-NGHOST,NZlocal+NGHOST
            do jj = 1-NGHOST,NYlocal+NGHOST
              zbuffer_send(jj,mm,l) = uinlett(jj,mm,l)
            enddo
          enddo
        enddo
        call MPI_Isend(zbuffer_send, ZBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        ipc_to, MSG_XCH_ZHI_TAG, comm3D, msg_id_send_re, ierr)
        call ErrorHandler(ierr,ERROR_SEND)             
      endif
c
      if(iprocx.eq. 1) then
        iproc_re = iprocnum(iprocxre, yprocs+1-iprocy, iprocz) 
        ipc_from = my_iproc(iproc_re)
        call MPI_Irecv(zbuffer_recv, ZBUFFSIZE, MPI_DOUBLE_PRECISION,
     &      ipc_from, MSG_XCH_ZHI_TAG, comm3D, msg_id_recv_re, ierr) 
        call ErrorHandler(ierr,ERROR_RECV)
      endif
c
      if(iprocx.eq. 1) then
        call MPI_Wait(msg_id_recv_re, status, ierr)
        call ErrorHandler(ierr,ERROR_WAIT)
        do jj=iylo, iyhi, 1
          do kk=izlo, izhi, 1
            Uin(1,jj,kk,1) =  zbuffer_recv(nylocal+1-jj,kk,1)
            Uin(1,jj,kk,2) = -zbuffer_recv(nylocal+1-jj,kk,2)
            Uin(1,jj,kk,3) =  zbuffer_recv(nylocal+1-jj,kk,3)
          enddo 
        enddo
      endif
c
      if(iprocx.eq. iprocxre) then
        call MPI_Wait(msg_id_send_re, status, ierr)
        call ErrorHandler(ierr,ERROR_WAIT)
      endif
c
#endif
c end PARALLEL
      return
      end subroutine exchangeinflow_mirror
c 
c----------------------------------------------------------------------
c     recycling inflow velocity, from iprocxre to iprocx=1.
c     mirror method
c----------------------------------------------------------------------
      subroutine exchangeinflow5_mirror(iprocxre,uinletn)
      use mesh_parms
      use mesh_common
      use InflowVelBC  
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      integer,parameter:: nfields=3
      integer:: iprocxre
      double precision:: uinletn(nghost+1,iylo:iyhi,izlo:izhi,3)
c      
#ifdef PARALLEL
      double precision, dimension
     &  (nghost+1,IYLO:IYHI,IZLO:IZHI,nfields) 
     &     :: zbuffer_send, zbuffer_recv
      integer:: XBUFFSIZE, YBUFFSIZE, ZBUFFSIZE
      integer:: ii,jj,mm,kk,l,idest, ng
      integer:: iprocnum, iproc_re
      integer:: iipy, iipz, ipc_to, ipc_from
      integer:: msg_id_send_re
      integer:: msg_id_recv_re
c
c     Z DIRECTION COMMUNICATION, done for 3D simulations only
      ZBUFFSIZE = (nghost+1)*nfields*(IYHI-IYLO+1)*(IZHI-IZLO+1)    
      if(iprocx.eq. iprocxre) then
        iproc_re = iprocnum(1, yprocs+1-iprocy, iprocz)
        ipc_to   =  my_iproc(iproc_re)
        do ng=1, nghost+1, 1
        do l = 1,nfields
          do mm = 1-NGHOST,NZlocal+NGHOST
            do jj = 1-NGHOST,NYlocal+NGHOST
              zbuffer_send(ng,jj,mm,l) = uinletn(ng,jj,mm,l)
            enddo
          enddo
        enddo
        enddo
        call MPI_Isend(zbuffer_send, ZBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        ipc_to, MSG_XCH_ZHI_TAG, comm3D, msg_id_send_re, ierr)
        call ErrorHandler(ierr,ERROR_SEND)             
      endif
c
      if(iprocx.eq. 1) then
        iproc_re = iprocnum(iprocxre, yprocs+1-iprocy, iprocz) 
        ipc_from = my_iproc(iproc_re)
        call MPI_Irecv(zbuffer_recv, ZBUFFSIZE, MPI_DOUBLE_PRECISION,
     &      ipc_from, MSG_XCH_ZHI_TAG, comm3D, msg_id_recv_re, ierr) 
        call ErrorHandler(ierr,ERROR_RECV)
      endif
c
      if(iprocx.eq. 1) then
        call MPI_Wait(msg_id_recv_re, status, ierr)
        call ErrorHandler(ierr,ERROR_WAIT)
        do ng=1, nghost+1, 1
        do jj=iylo, iyhi, 1
          do kk=izlo, izhi, 1
            Uin(ng,jj,kk,1) =  zbuffer_recv(ng,nylocal+1-jj,kk,1)
            Uin(ng,jj,kk,2) = -zbuffer_recv(ng,nylocal+1-jj,kk,2)
            Uin(ng,jj,kk,3) =  zbuffer_recv(ng,nylocal+1-jj,kk,3)
          enddo 
        enddo
        enddo
      endif
c
      if(iprocx.eq. iprocxre) then
        call MPI_Wait(msg_id_send_re, status, ierr)
        call ErrorHandler(ierr,ERROR_WAIT)
      endif
c
#endif
c end PARALLEL
      return
      end subroutine exchangeinflow5_mirror
c 
c----------------------------------------------------------------------
c
c     recycling inflow velocity, from iprocxre to iprocx=1.
c----------------------------------------------------------------------
      subroutine exchangeinflow(iprocxre,uinlett)
      use mesh_parms
      use mesh_common
      use InflowVelBC  
#ifdef PARALLEL
      use mpistuff
#endif 
      implicit none
      integer,parameter:: nfields=3
      integer:: iprocxre
      double precision:: uinlett(iylo:iyhi,izlo:izhi,3)
c      
#ifdef PARALLEL
      double precision, dimension
     &  (IYLO:IYHI,IZLO:IZHI,nfields) 
     &     :: zbuffer_send, zbuffer_recv
      integer:: XBUFFSIZE, YBUFFSIZE, ZBUFFSIZE
      integer:: ii,jj,mm,kk,l,idest
      integer:: iprocnum, iproc_re
      integer:: iipy, iipz, ipc_to, ipc_from
      integer:: msg_id_send_re
      integer:: msg_id_recv_re
c     
c     Z DIRECTION COMMUNICATION, done for 3D simulations only
      ZBUFFSIZE = nfields*(IYHI-IYLO+1)*(IZHI-IZLO+1)
      if(iprocx.eq. iprocxre) then
        iproc_re = iprocnum(1, iprocy, iprocz)
        ipc_to   =  my_iproc(iproc_re)
        do l = 1,nfields
          do mm = 1-NGHOST,NZlocal+NGHOST
            do jj = 1-NGHOST,NYlocal+NGHOST
              zbuffer_send(jj,mm,l) = uinlett(jj,mm,l)
            enddo
          enddo
        enddo
        call MPI_Isend(zbuffer_send, ZBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        ipc_to, MSG_XCH_ZHI_TAG, comm3D, msg_id_send_re, ierr)
        call ErrorHandler(ierr,ERROR_SEND)             
      endif
c
      if(iprocx.eq. 1) then
        iproc_re = iprocnum(iprocxre, iprocy, iprocz) 
        ipc_from = my_iproc(iproc_re)
        call MPI_Irecv(zbuffer_recv, ZBUFFSIZE, MPI_DOUBLE_PRECISION,
     &      ipc_from, MSG_XCH_ZHI_TAG, comm3D, msg_id_recv_re, ierr) 
        call ErrorHandler(ierr,ERROR_RECV)
      endif
c
      if(iprocx.eq. 1) then
        call MPI_Wait(msg_id_recv_re, status, ierr)
        call ErrorHandler(ierr,ERROR_WAIT)
        do jj=iylo, iyhi, 1
          do kk=izlo, izhi, 1
            Uin(1,jj,kk,1) =  zbuffer_recv(jj,kk,1)
            Uin(1,jj,kk,2) =  zbuffer_recv(jj,kk,2)
            Uin(1,jj,kk,3) =  zbuffer_recv(jj,kk,3)
          enddo 
        enddo
      endif
c
      if(iprocx.eq. iprocxre) then
        call MPI_Wait(msg_id_send_re, status, ierr)
        call ErrorHandler(ierr,ERROR_WAIT)
      endif
c
#endif
c end PARALLEL
      return
      end subroutine exchangeinflow
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c     recycling inflow velocity, from iprocxre to iprocx=1.
c----------------------------------------------------------------------
      subroutine exchangeinflow5(iprocxre,uinletn)
      use mesh_parms
      use mesh_common
      use InflowVelBC  
#ifdef PARALLEL
      use mpistuff
#endif 
      implicit none
      integer,parameter:: nfields=3
      integer:: iprocxre
      double precision:: uinletn(nghost+1,iylo:iyhi,izlo:izhi,3)
c      
#ifdef PARALLEL
      double precision, dimension
     &  (nghost+1,IYLO:IYHI,IZLO:IZHI,nfields) 
     &     :: zbuffer_send, zbuffer_recv
      integer:: XBUFFSIZE, YBUFFSIZE, ZBUFFSIZE
      integer:: ii,jj,mm,kk,l,idest, ng
      integer:: iprocnum, iproc_re
      integer:: iipy, iipz, ipc_to, ipc_from
      integer:: msg_id_send_re
      integer:: msg_id_recv_re
c     
c     Z DIRECTION COMMUNICATION, done for 3D simulations only
      ZBUFFSIZE = (nghost+1)*nfields*(IYHI-IYLO+1)*(IZHI-IZLO+1)
      if(iprocx.eq. iprocxre) then
        iproc_re = iprocnum(1, iprocy, iprocz)
        ipc_to   =  my_iproc(iproc_re)
        do ng= 1, nghost+1,1
        do l = 1,nfields
          do mm = 1-NGHOST,NZlocal+NGHOST
            do jj = 1-NGHOST,NYlocal+NGHOST
              zbuffer_send(ng,jj,mm,l) = uinletn(ng,jj,mm,l)
            enddo
          enddo
        enddo
        enddo
        call MPI_Isend(zbuffer_send, ZBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        ipc_to, MSG_XCH_ZHI_TAG, comm3D, msg_id_send_re, ierr)
        call ErrorHandler(ierr,ERROR_SEND)             
      endif
c
      if(iprocx.eq. 1) then
        iproc_re = iprocnum(iprocxre, iprocy, iprocz) 
        ipc_from = my_iproc(iproc_re)
        call MPI_Irecv(zbuffer_recv, ZBUFFSIZE, MPI_DOUBLE_PRECISION,
     &      ipc_from, MSG_XCH_ZHI_TAG, comm3D, msg_id_recv_re, ierr) 
        call ErrorHandler(ierr,ERROR_RECV)
      endif
c
      if(iprocx.eq. 1) then
        call MPI_Wait(msg_id_recv_re, status, ierr)
        call ErrorHandler(ierr,ERROR_WAIT)
        do ng=1, nghost+1, 1
        do jj=iylo, iyhi, 1
          do kk=izlo, izhi, 1
            Uin(ng,jj,kk,1) =  zbuffer_recv(ng,jj,kk,1)
            Uin(ng,jj,kk,2) =  zbuffer_recv(ng,jj,kk,2)
            Uin(ng,jj,kk,3) =  zbuffer_recv(ng,jj,kk,3)
          enddo 
        enddo
        enddo
      endif
c
      if(iprocx.eq. iprocxre) then
        call MPI_Wait(msg_id_send_re, status, ierr)
        call ErrorHandler(ierr,ERROR_WAIT)
      endif
c
#endif
c end PARALLEL
      return
      end subroutine exchangeinflow5
c
c----------------------------------------------------------------------

