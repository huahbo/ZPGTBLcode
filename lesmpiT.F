c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2010
c     King Abdullah University of Science and Technology
c     All Rights Reserved
c-----------------------------------------------------------------
c     $Log: bdryexchange.F,v $
c     Revision 1.1.1.1  2010/12/31 12:11:16  samtanr
c     Plasmoid reconnection MHD code. 
c
c-----------------------------------------------------------------
c
c     ================================================================
      
      subroutine exchangeTct(grid_array,nfields,
     & ITXLO,ITXHI,ITYLO,ITYHI,ITZLO,ITZHI )
c     
c     precondition:
c               grid_array holds the local data to be communicated
c               nfields holds the number of distinct fluid variables
c     postcondition:
c               grid_array holds the local data, post communication
c
c     Note: there are a number of #ifdef statements, meant to test 
c     for directional periodicity, depending on whether processes 
c     are on the domain boundary or not.
c     
      use mesh_parms
      use mesh_common
#ifdef PARALLEL
      use mpistuff
#endif
      integer,parameter:: nghostTensor=2
      integer:: ITXLO  != -1
      integer:: ITXHI  != nxlocal+2
      integer:: ITYLO  != -1
      integer:: ITYHI  != nylocal+2
      integer:: ITZLO  != -1
      integer:: ITZHI  != nzlocal+2     
c
      integer :: nfields
      double precision, dimension
     &   (ITXLO:ITXHI,ITYLO:ITYHI,ITZLO:ITZHI,nfields) 
     &     :: grid_array
      
#ifdef PARALLEL
      double precision, dimension
     &  (NGHOSTTENSOR,ITYLO:ITYHI,ITZLO:ITZHI,nfields) 
     &     :: xbuffer_send, xbuffer_recv
#ifndef ONE_D
      double precision, dimension
     &  (ITXLO:ITXHI,NGHOSTTENSOR,ITZLO:ITZHI,nfields) 
     &     :: ybuffer_send, ybuffer_recv
#ifndef TWO_D
      double precision, dimension
     &  (ITXLO:ITXHI,ITYLO:ITYHI,NGHOSTTENSOR,nfields) 
     &     :: zbuffer_send, zbuffer_recv
#endif
#endif

      integer:: XBUFFSIZE, YBUFFSIZE, ZBUFFSIZE
      integer:: ii,jj,mm,kk,l,idest
      integer:: iprocnum
      integer:: msg_id_send_x_low
      integer:: msg_id_send_x_hi
      integer:: msg_id_recv_x_low
      integer:: msg_id_recv_x_hi
      
      integer:: msg_id_send_y_low
      integer:: msg_id_send_y_hi
      integer:: msg_id_recv_y_low
      integer:: msg_id_recv_y_hi
      
      integer:: msg_id_send_z_low
      integer:: msg_id_send_z_hi
      integer:: msg_id_recv_z_low
      integer:: msg_id_recv_z_hi
c     
c
c     set buffer sizes
      XBUFFSIZE =nfields*(NGHOSTTENSOR)*(ITYHI-ITYLO+1)*(ITZHI-ITZLO+1)
      YBUFFSIZE =nfields*(ITXHI-ITXLO+1)*(NGHOSTTENSOR)*(ITZHI-ITZLO+1)
      ZBUFFSIZE =nfields*(ITXHI-ITXLO+1)*(ITYHI-ITYLO+1)*(NGHOSTTENSOR) 

      
c     -------X DIRECTION COMMUNICATION, done for all simulations
c
c     
c     Update x-low boundaries
c
#ifndef XPERIODIC
      if ( iprocx > 1 ) then
#endif
         call MPI_Irecv(xbuffer_recv, XBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        left, MSG_XCH_XHI_TAG, comm3D, msg_id_recv_x_hi, ierr)
         call ErrorHandler(ierr,ERROR_RECV)
#ifndef XPERIODIC
      endif
#endif


#ifndef XPERIODIC
      if ( iprocx < XPROCS ) then
#endif
         do l = 1,nfields
            do kk = ITZLO,ITZHI
               do jj = 1-NGHOSTTENSOR,NYlocal+NGHOSTTENSOR
                  do mm = 1,NGHOSTTENSOR
                     xbuffer_send(mm,jj,kk,l) = 
     &                    grid_array(NXlsize+1-mm,jj,kk,l)
                  enddo
               enddo
            enddo
         enddo
         
         call MPI_Isend(xbuffer_send, XBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        right, MSG_XCH_XHI_TAG, comm3D, msg_id_send_x_hi, ierr)
         call ErrorHandler(ierr,ERROR_SEND)
#ifndef XPERIODIC
      endif
#endif


#ifndef XPERIODIC
      if ( iprocx > 1 ) then
#endif
         call MPI_Wait(msg_id_recv_x_hi, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
         
         do l = 1,nfields
            do kk = ITZLO,ITZHI
               do jj = 1-NGHOSTTENSOR,NYlocal+NGHOSTTENSOR
                  do mm = 1,NGHOSTTENSOR
                     grid_array(1-mm,jj,kk,l) = 
     &                    xbuffer_recv(mm,jj,kk,l)
                  enddo
               enddo
            enddo
         enddo
#ifndef XPERIODIC
      endif
#endif
      

#ifndef XPERIODIC
      if ( iprocx < XPROCS ) then
#endif
         call MPI_Wait(msg_id_send_x_hi, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
#ifndef XPERIODIC
      endif
#endif
c
c      
c     update x-high boundaries
c
#ifndef XPERIODIC
      if ( iprocx < XPROCS ) then
#endif
         call MPI_Irecv(xbuffer_recv, XBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        right, MSG_XCH_XLOW_TAG, comm3D, msg_id_recv_x_low, ierr)
         call ErrorHandler(ierr,ERROR_RECV)
#ifndef XPERIODIC
      endif
#endif
      

#ifndef XPERIODIC
      if ( iprocx > 1 ) then
#endif
         do l = 1,nfields
            do kk = ITZLO,ITZHI
               do jj = 1-NGHOSTTENSOR,NYlocal+NGHOSTTENSOR
                  do mm = 1,NGHOSTTENSOR
                     xbuffer_send(mm,jj,kk,l) = 
     &                    grid_array(mm,jj,kk,l)
                  enddo
               enddo
            enddo
         enddo
         
         call MPI_Isend(xbuffer_send, XBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        left, MSG_XCH_XLOW_TAG, comm3D, msg_id_send_x_low, ierr)
         call ErrorHandler(ierr,ERROR_SEND)
#ifndef XPERIODIC
      endif
#endif

      
#ifndef XPERIODIC
      if ( iprocx < XPROCS ) then
#endif
         call MPI_Wait(msg_id_recv_x_low, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
         
         do l = 1,nfields
            do kk = ITZLO,ITZHI
               do jj = 1-NGHOSTTENSOR,NYlocal+NGHOSTTENSOR
                  do mm = 1,NGHOSTTENSOR
                     grid_array(NXlsize+mm,jj,kk,l) = 
     &                    xbuffer_recv(mm,jj,kk,l)
                  enddo
               enddo
            enddo
         enddo
#ifndef XPERIODIC
      endif
#endif

      
#ifndef XPERIODIC
      if ( iprocx > 1 ) then
#endif
         call MPI_Wait(msg_id_send_x_low, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
#ifndef XPERIODIC
      endif
#endif
c
c     
c     
c     -------Y DIRECTION COMMUNICATION, done for 2D, 2.5D, 3D runs
#ifndef ONE_D
c     
c     
c     update y-low boundaries
#ifndef YPERIODIC
      if ( iprocy > 1 ) then
#endif
         call  MPI_Irecv(ybuffer_recv, YBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        bottom, MSG_XCH_YHI_TAG, comm3D, msg_id_recv_y_hi, ierr)
         call ErrorHandler(ierr,ERROR_RECV)
#ifndef YPERIODIC
      endif
#endif
      

#ifndef YPERIODIC
      if ( iprocy < YPROCS ) then
#endif
         do l = 1,nfields
            do kk = ITZLO,ITZHI
               do mm = 1,NGHOSTTENSOR
                  do ii = 1-NGHOSTTENSOR,NXlocal+NGHOSTTENSOR
                     ybuffer_send(ii,mm,kk,l) = 
     &                    grid_array(ii,NYlsize+1-mm,kk,l)
                  enddo
               enddo
            enddo
         enddo
         
         call MPI_Isend(ybuffer_send, YBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        top, MSG_XCH_YHI_TAG, comm3D, msg_id_send_y_hi, ierr)
         call ErrorHandler(ierr,ERROR_SEND)
#ifndef YPERIODIC
      endif
#endif

      
#ifndef YPERIODIC
      if ( iprocy > 1 ) then
#endif
         call MPI_Wait(msg_id_recv_y_hi, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
         
         do l = 1,nfields
            do kk = ITZLO,ITZHI
               do mm = 1,NGHOSTTENSOR
                  do ii = 1-NGHOSTTENSOR,NXlocal+NGHOSTTENSOR
                     grid_array(ii,1-mm,kk,l) = 
     &                    ybuffer_recv(ii,mm,kk,l)
                  enddo
               enddo
            enddo
         enddo
#ifndef YPERIODIC
      endif
#endif
      

#ifndef YPERIODIC
      if ( iprocy < YPROCS ) then
#endif
         call MPI_Wait(msg_id_send_y_hi, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
#ifndef YPERIODIC
      endif
#endif
c
c      
c     update y-high boundaries
c
#ifndef YPERIODIC
      if ( iprocy < YPROCS ) then
#endif
         call MPI_Irecv(ybuffer_recv, YBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        top, MSG_XCH_YLOW_TAG, comm3D, msg_id_recv_y_low, ierr)
         call ErrorHandler(ierr,ERROR_RECV)
#ifndef YPERIODIC
      endif
#endif

      
#ifndef YPERIODIC
      if ( iprocy > 1 ) then
#endif
         do l = 1,nfields
            do kk = ITZLO,ITZHI
               do mm = 1,NGHOSTTENSOR
                  do ii = 1-NGHOSTTENSOR,NXlocal+NGHOSTTENSOR
                     ybuffer_send(ii,mm,kk,l) = 
     &                    grid_array(ii,mm,kk,l)
                  enddo
               enddo
            enddo
         enddo
         
         call MPI_Isend(ybuffer_send, YBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        bottom, MSG_XCH_YLOW_TAG, comm3D, msg_id_send_y_low, ierr)
         call ErrorHandler(ierr,ERROR_SEND)
#ifndef YPERIODIC
      endif
#endif

      
#ifndef YPERIODIC
      if ( iprocy < YPROCS ) then
#endif
         call MPI_Wait(msg_id_recv_y_low, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
         
         do l = 1,nfields
            do kk = ITZLO,ITZHI
               do mm = 1,NGHOSTTENSOR
                  do ii = 1-NGHOSTTENSOR,NXlocal+NGHOSTTENSOR
                     grid_array(ii,NYlsize+mm,kk,l) = 
     &                    ybuffer_recv(ii,mm,kk,l)
                  enddo
               enddo
            enddo
         enddo
#ifndef YPERIODIC
      endif
#endif

      
#ifndef YPERIODIC
      if ( iprocy > 1 ) then
#endif
         call MPI_Wait(msg_id_send_y_low, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
#ifndef YPERIODIC
      endif
#endif
c
c     
c     
c     -------Z DIRECTION COMMUNICATION, done for 3D simulations only
#ifndef TWO_D
c
c     
c     update z-low boundaries
#ifndef ZPERIODIC
      if ( iprocz > 1 ) then
#endif
         call MPI_Irecv(zbuffer_recv, ZBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        behind, MSG_XCH_ZHI_TAG, comm3D, msg_id_recv_z_hi, ierr) 
         call ErrorHandler(ierr,ERROR_RECV)
#ifndef ZPERIODIC
      endif
#endif

      
#ifndef ZPERIODIC
      if ( iprocz < ZPROCS ) then
#endif
         do l = 1,nfields
            do mm = 1,NGHOSTTENSOR
               do jj = 1-NGHOSTTENSOR,NYlocal+NGHOSTTENSOR
                  do ii = 1-NGHOSTTENSOR,NXlocal+NGHOSTTENSOR
                     zbuffer_send(ii,jj,mm,l) = 
     &                    grid_array(ii,jj,NZlsize+1-mm,l)
                  enddo
               enddo
            enddo
         enddo
         
         call MPI_Isend(zbuffer_send, ZBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        forward, MSG_XCH_ZHI_TAG, comm3D, msg_id_send_z_hi, ierr)
         call ErrorHandler(ierr,ERROR_SEND)
#ifndef ZPERIODIC
      endif
#endif
      

#ifndef ZPERIODIC
      if ( iprocz > 1 ) then
#endif
         call MPI_Wait(msg_id_recv_z_hi, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
         
         do l = 1,nfields
            do mm = 1,NGHOSTTENSOR
               do jj = 1-NGHOSTTENSOR,NYlocal+NGHOSTTENSOR
                  do ii = 1-NGHOSTTENSOR,NXlocal+NGHOSTTENSOR
                     grid_array(ii,jj,1-mm,l) = 
     &                    zbuffer_recv(ii,jj,mm,l)
                  enddo
               enddo
            enddo
         enddo
#ifndef ZPERIODIC
      endif
#endif
      

#ifndef ZPERIODIC
      if ( iprocz < ZPROCS ) then
#endif
         call MPI_Wait(msg_id_send_z_hi, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
#ifndef ZPERIODIC
      endif
#endif
c      
c     update z-high boundaries
c
#ifndef ZPERIODIC
      if ( iprocz < ZPROCS ) then
#endif
         call MPI_Irecv(zbuffer_recv, ZBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        forward, MSG_XCH_ZLOW_TAG, comm3D,msg_id_recv_z_low, ierr)
         call ErrorHandler(ierr,ERROR_RECV)
#ifndef ZPERIODIC
      endif
#endif

      
#ifndef ZPERIODIC
      if ( iprocz > 1 ) then
#endif
         do l = 1,nfields
            do mm = 1,NGHOSTTENSOR
               do jj = 1-NGHOSTTENSOR,NYlocal+NGHOSTTENSOR
                  do ii = 1-NGHOSTTENSOR,NXlocal+NGHOSTTENSOR
                     zbuffer_send(ii,jj,mm,l) = 
     &                    grid_array(ii,jj,mm,l)
                  enddo
               enddo
            enddo
         enddo
         
         call MPI_Isend(zbuffer_send, ZBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        behind, MSG_XCH_ZLOW_TAG, comm3D, msg_id_send_z_low, ierr)
         call ErrorHandler(ierr,ERROR_SEND)
#ifndef ZPERIODIC
      endif
#endif

      
#ifndef ZPERIODIC
      if ( iprocz < ZPROCS ) then
#endif
         call MPI_Wait(msg_id_recv_z_low, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
         
         do l = 1,nfields
            do mm = 1,NGHOSTTENSOR
               do jj = 1-NGHOSTTENSOR,NYlocal+NGHOSTTENSOR
                  do ii = 1-NGHOSTTENSOR,NXlocal+NGHOSTTENSOR
                     grid_array(ii,jj,NZlsize+mm,l) = 
     &                    zbuffer_recv(ii,jj,mm,l)
                  enddo
               enddo
            enddo
         enddo
#ifndef ZPERIODIC
      endif
#endif

      
#ifndef ZPERIODIC
      if ( iprocz > 1 ) then
#endif
         call MPI_Wait(msg_id_send_z_low, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
#ifndef ZPERIODIC
      endif
#endif

#endif
c end TWO_D
#endif
c end ONE_D
#endif
c end PARALLEL

      return
      end subroutine exchangeTct
c-----------------------------------------------------------------------
c
      subroutine ErrorHandlerT(mpierr,errortype)
      use mesh_common
#ifdef PARALLEL
      use mpistuff
#endif
      integer::mpierr,errortype
#ifdef PARALLEL
      if(mpierr.ne.MPI_SUCCESS) then
         write(0,*) 'FLUID: MPI RETURN VALUE',iproc_idx,mpierr,errortype
      endif
#endif
      return
      end subroutine ErrorHandlerT
c-----------------------------------------------------------------------

c----------------------------------------------------------------------
c exchange buffer region in Uwall
      subroutine exchangebcuwall
c     
      use mesh_parms
      use mesh_common
      use virtualwallBC
#ifdef PARALLEL
      use mpistuff
#endif
      integer,parameter:: nfields=3
c      
#ifdef PARALLEL
      double precision, dimension
     &  (ixlo:ixhi,NGHOST,nfields) 
     &     :: ybuffer_send, ybuffer_recv
      double precision, dimension
     &  (NGHOST,iylo:iyhi,nfields)
     &     :: xbuffer_send, xbuffer_recv

      integer:: YBUFFSIZE, XBUFFSIZE
      integer:: ii,jj,mm,kk,l,idest
      integer:: iprocnum
      
      integer:: msg_id_send_y_low
      integer:: msg_id_send_y_hi
      integer:: msg_id_recv_y_low
      integer:: msg_id_recv_y_hi
c
      integer:: msg_id_send_x_low
      integer:: msg_id_send_x_hi
      integer:: msg_id_recv_x_low
      integer:: msg_id_recv_x_hi     

c
c     set buffer sizes
      YBUFFSIZE = nfields*(ixhi-ixlo+1)*(NGHOST) 
c 
c     -------Y DIRECTION COMMUNICATION, done for 2D
c     
c     
             
c     update y-low boundaries
#ifndef YPERIODIC
      if ( iprocy > 1 ) then
#endif
         call  MPI_Irecv(ybuffer_recv, YBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        bottom, MSG_XCH_YHI_TAG, comm3D, msg_id_recv_y_hi, ierr)
         call ErrorHandler(ierr,ERROR_RECV)
#ifndef YPERIODIC
      endif
#endif
      

#ifndef YPERIODIC
      if ( iprocy < YPROCS ) then
#endif
         do l = 1,nfields
               do mm = 1,NGHOST
                  do ii = ixlo, ixhi
                     ybuffer_send(ii,mm,l) = 
     &                   bc_uwall(ii,nylocal+1-mm,l)
                  enddo
               enddo
          enddo
         
         call MPI_Isend(ybuffer_send, YBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        top, MSG_XCH_YHI_TAG, comm3D, msg_id_send_y_hi, ierr)
         call ErrorHandler(ierr,ERROR_SEND)
#ifndef YPERIODIC
      endif
#endif

      
#ifndef YPERIODIC
      if ( iprocy > 1 ) then
#endif
         call MPI_Wait(msg_id_recv_y_hi, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
         
         do l = 1,nfields
               do mm = 1,NGHOST
                  do ii = ixlo, ixhi
                     bc_uwall(ii,1-mm,l) = 
     &                    ybuffer_recv(ii,mm,l)
                  enddo
               enddo
         enddo
#ifndef YPERIODIC
      endif
#endif
      

#ifndef YPERIODIC
      if ( iprocy < YPROCS ) then
#endif
         call MPI_Wait(msg_id_send_y_hi, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
#ifndef YPERIODIC
      endif
#endif
c
c      
c     update y-high boundaries
c
#ifndef YPERIODIC
      if ( iprocy < YPROCS ) then
#endif
         call MPI_Irecv(ybuffer_recv, YBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        top, MSG_XCH_YLOW_TAG, comm3D, msg_id_recv_y_low, ierr)
         call ErrorHandler(ierr,ERROR_RECV)
#ifndef YPERIODIC
      endif
#endif

      
#ifndef YPERIODIC
      if ( iprocy > 1 ) then
#endif
         do l = 1,nfields
               do mm = 1,NGHOST
                  do ii = ixlo, ixhi
                     ybuffer_send(ii,mm,l) = 
     &                    bc_uwall(ii,mm,l)
                  enddo
               enddo
         enddo
         
         call MPI_Isend(ybuffer_send, YBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        bottom, MSG_XCH_YLOW_TAG, comm3D, msg_id_send_y_low, ierr)
         call ErrorHandler(ierr,ERROR_SEND)
#ifndef YPERIODIC
      endif
#endif

      
#ifndef YPERIODIC
      if ( iprocy < YPROCS ) then
#endif
         call MPI_Wait(msg_id_recv_y_low, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
         
         do l = 1,nfields
               do mm = 1,NGHOST
                  do ii = ixlo, ixhi
                     bc_uwall(ii,NYlsize+mm,l) = 
     &                    ybuffer_recv(ii,mm,l)
                  enddo
               enddo
         enddo
#ifndef YPERIODIC
      endif
#endif

      
#ifndef YPERIODIC
      if ( iprocy > 1 ) then
#endif
         call MPI_Wait(msg_id_send_y_low, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
#ifndef YPERIODIC
      endif
#endif
c
c
c
c     
      XBUFFSIZE = nfields*(iyhi-iylo+1)*(NGHOST) 
c     
c     -------X DIRECTION COMMUNICATION, done for 2D simulations only
c
c     
c     update x-low boundaries
#ifndef XPERIODIC
      if ( iprocx > 1 ) then
#endif
         call MPI_Irecv(xbuffer_recv, XBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        left, MSG_XCH_XHI_TAG, comm3D, msg_id_recv_x_hi, ierr) 
         call ErrorHandler(ierr,ERROR_RECV)
#ifndef XPERIODIC
      endif
#endif

      
#ifndef XPERIODIC
      if ( iprocx < XPROCS ) then
#endif
         do l = 1,nfields
            do mm = 1,NGHOST
               do jj = 1-NGHOST,NYlocal+NGHOST
                     xbuffer_send(mm,jj,l) = 
     &                    bc_uwall(NXlsize+1-mm,jj,l)
               enddo
            enddo
         enddo
         
         call MPI_Isend(xbuffer_send, XBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        right, MSG_XCH_XHI_TAG, comm3D, msg_id_send_x_hi, ierr)
         call ErrorHandler(ierr,ERROR_SEND)
#ifndef XPERIODIC
      endif
#endif
      

#ifndef XPERIODIC
      if ( iprocx > 1 ) then
#endif
         call MPI_Wait(msg_id_recv_x_hi, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
         
         do l = 1,nfields
            do mm = 1,NGHOST
               do jj = 1-NGHOST,NYlocal+NGHOST
                     bc_uwall(1-mm,jj,l) = 
     &                    xbuffer_recv(mm,jj,l)
               enddo
            enddo
         enddo
#ifndef XPERIODIC
      endif
#endif
      

#ifndef XPERIODIC
      if ( iprocx < XPROCS ) then
#endif
         call MPI_Wait(msg_id_send_x_hi, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
#ifndef XPERIODIC
      endif
#endif
c      
c     update x-high boundaries
c
#ifndef XPERIODIC
      if ( iprocx < XPROCS ) then
#endif
         call MPI_Irecv(xbuffer_recv, XBUFFSIZE, MPI_DOUBLE_PRECISION,
     &        right, MSG_XCH_XLOW_TAG, comm3D, msg_id_recv_x_low, ierr)
         call ErrorHandler(ierr,ERROR_RECV)
#ifndef XPERIODIC
      endif
#endif

      
#ifndef XPERIODIC
      if ( iprocx > 1 ) then
#endif
         do l = 1,nfields
            do mm = 1,NGHOST
               do jj = 1-NGHOST,NYlocal+NGHOST
                     xbuffer_send(mm,jj,l) = 
     &                    bc_uwall(mm,jj,l)
               enddo
            enddo
         enddo
         
         call MPI_Isend(xbuffer_send, XBUFFSIZE, MPI_DOUBLE_PRECISION,
     &         left, MSG_XCH_XLOW_TAG, comm3D, msg_id_send_x_low, ierr)
         call ErrorHandler(ierr,ERROR_SEND)
#ifndef XPERIODIC
      endif
#endif

     
#ifndef XPERIODIC
      if ( iprocx < XPROCS ) then
#endif
         call MPI_Wait(msg_id_recv_x_low, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
         
         do l = 1,nfields
            do mm = 1,NGHOST
               do jj = 1-NGHOST,NYlocal+NGHOST
                     bc_uwall(NXlsize+mm,jj,l) = 
     &                    xbuffer_recv(mm,jj,l)
               enddo
            enddo
         enddo
#ifndef XPERIODIC
      endif
#endif

      
#ifndef XPERIODIC
      if ( iprocx > 1 ) then
#endif
         call MPI_Wait(msg_id_send_x_low, status, ierr)
         call ErrorHandler(ierr,ERROR_WAIT)
#ifndef XPERIODIC
      endif
#endif
c
c     
#endif 
c
      return
      end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
