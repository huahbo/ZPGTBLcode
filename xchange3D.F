c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2010
c     King Abdullah University of Science and Technology
c     All Rights Reserved
c-----------------------------------------------------------------
c     $Log: xchange3D.F,v $
c     Revision 1.4  2011/10/19 15:14:49  samtanr
c     Added 2D option via TWO_D defines.
c
c     Revision 1.3  2011/10/08 13:39:00  samtanr
c     Fixed exchanges.
c
c     Revision 1.2  2011/10/07 13:21:09  samtanr
c     Fixed bug on size of array.
c
c     Revision 1.1  2011/07/25 11:32:01  samtanr
c     Original source.
c
c     Revision 1.2  2011/01/09 08:14:49  samtanr
c     Removed some tabs; and reformatted some code.
c
c     Revision 1.1.1.1  2010/12/31 12:11:16  samtanr
c     Plasmoid reconnection MHD code. 
c
c-----------------------------------------------------------------
        subroutine ExchangePsi3D(psi,imax,jmax,kmax,neqn,xChangeFlag)
c-----------------------------------------------------------------
c       Description: Exchange values at overlapping processor
c                    boundaries for psi
c       Input: psi - field with ghost which must be exchanged.
c              imax,jmax,kmax: mesh extent
c              xChangeFlag: 0 (no exchange in which case this is 
c                              an empty routine)
c                         : 1 (do the exchange)
c-----------------------------------------------------------------
        use mesh_uparms
        use mesh_common
#ifdef PARALLEL
        use mpistuff
#endif

        integer,parameter:: MSG_XCH_XLOW_TAG=1,MSG_XCH_XHI_TAG=2
        integer,parameter:: MSG_XCH_YLOW_TAG=3,MSG_XCH_YHI_TAG=4
        integer,parameter:: MSG_XCH_ZLOW_TAG=5,MSG_XCH_ZHI_TAG=6
        integer,parameter:: MSG_MAX_TAG=10
c
        integer:: imax,jmax,kmax,neqn,xChangeFlag
#ifdef TWO_D
        double precision:: psi(-2:imax+3,jmax,-2:kmax+3,neqn)
#else
        double precision:: psi(-2:imax+3,-2:jmax+3,-2:kmax+3,neqn)
#endif
        integer, parameter:: ng=3

#ifdef PARALLEL
#ifndef TWO_D
        double precision:: xbuffer_send(ng,jmax+1,kmax+1,neqn)
        double precision:: xbuffer_recv(ng,jmax+1,kmax+1,neqn)
        double precision:: ybuffer_send(imax+1,ng,kmax+1,neqn)
        double precision:: ybuffer_recv(imax+1,ng,kmax+1,neqn)
        double precision:: zbuffer_send(imax+1,jmax+1,ng,neqn)
        double precision:: zbuffer_recv(imax+1,jmax+1,ng,neqn)
#else
        double precision:: xbuffer_send(ng,jmax,kmax+1,neqn)
        double precision:: xbuffer_recv(ng,jmax,kmax+1,neqn)
        double precision:: zbuffer_send(imax+1,jmax,ng,neqn)
        double precision:: zbuffer_recv(imax+1,jmax,ng,neqn)
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
        integer:: ybufsize, zbufsize
        integer:: jhi
c       
#ifndef TWO_D
        jhi=jmax+1
#else
        jhi=jmax
#endif

        if(xChangeFlag.eq.0) return

C       -------X DIRECTION COMMUNICATION
c       Update x-low boundaries

         XBUFFSIZE=(ng) * ((jhi)*(kmax+1))*neqn
         YBUFFSIZE=(ng) * ((imax+1)*(kmax+1))*neqn
         ZBUFFSIZE=(ng) * ((imax+1)*(jhi))*neqn
#ifndef XPERIODIC
      if ( iprocx > 1 ) then
#endif
c           isource = iprocnum (iprocx-1, iprocy, iprocz)-1           
           call MPI_Irecv(xbuffer_recv, XBUFFSIZE, MPI_DOUBLE_PRECISION,
     1           left, MSG_XCH_XHI_TAG, comm3D, msg_id_recv_x_hi, ierr
     2           )        
           Call ErrorHandler(ierr,ERROR_RECV)
#ifndef XPERIODIC
        endif
#endif
       
#ifndef XPERIODIC
      if ( iprocx < XPROCS ) then
#endif
           xbuffer_send(1:ng,:,:,:)=psi(imax-ng+1:imax,
     &        1:jhi,1:kmax+1,:)

           call MPI_Isend(xbuffer_send, XBUFFSIZE, MPI_DOUBLE_PRECISION,
     1                right, MSG_XCH_XHI_TAG, comm3D, msg_id_send_x_hi
     &                , ierr)
           Call ErrorHandler(ierr,ERROR_SEND)
#ifndef XPERIODIC
        endif
#endif
       
#ifndef XPERIODIC
      if ( iprocx > 1 ) then
#endif

           call MPI_Wait(msg_id_recv_x_hi, status, ierr)
           Call ErrorHandler(ierr,ERROR_WAIT)

           psi(0-ng+1:0,1:jhi,1:kmax+1,:)=xbuffer_recv(:,:,:,:)
#ifndef XPERIODIC
        endif
#endif

#ifndef XPERIODIC
      if ( iprocx < XPROCS ) then
#endif
           call MPI_Wait(msg_id_send_x_hi, status, ierr)
           Call ErrorHandler(ierr,ERROR_WAIT)
#ifndef XPERIODIC
        endif
#endif

C        update x-high boundaries

#ifndef XPERIODIC
      if ( iprocx < XPROCS ) then
#endif
c           isource = iprocnum (iprocx+1, iprocy, iprocz)-1           
           call MPI_Irecv(xbuffer_recv, XBUFFSIZE, MPI_DOUBLE_PRECISION,
     1           right, MSG_XCH_XLOW_TAG, comm3D, msg_id_recv_x_low,
     2           ierr) 
           Call ErrorHandler(ierr,ERROR_RECV)
#ifndef XPERIODIC
        endif
#endif

#ifndef XPERIODIC
      if ( iprocx > 1 ) then
#endif
           xbuffer_send(:,:,:,:)=psi(1:ng,1:jhi,1:kmax+1,:)

           call MPI_Isend(xbuffer_send, XBUFFSIZE, MPI_DOUBLE_PRECISION,
     1                left, MSG_XCH_XLOW_TAG, comm3D, msg_id_send_x_low
     &               , ierr)
           Call ErrorHandler(ierr,ERROR_SEND)
#ifndef XPERIODIC
        endif
#endif

#ifndef XPERIODIC
      if ( iprocx < XPROCS ) then
#endif

           call MPI_Wait(msg_id_recv_x_low, status, ierr)
           Call ErrorHandler(ierr,ERROR_WAIT)

           psi(imax+1:imax+ng,1:jhi,1:kmax+1,:)=xbuffer_recv(:,:,:,:)
#ifndef XPERIODIC
        endif
#endif

#ifndef XPERIODIC
      if ( iprocx > 1 ) then
#endif
           call MPI_Wait(msg_id_send_x_low, status, ierr)
           Call ErrorHandler(ierr,ERROR_WAIT)
#ifndef XPERIODIC
        endif
#endif

C       -------Y DIRECTION COMMUNICATION

C        update y-low boundaries
c        ybufsize=(ihi-ilo+1)*(IZHI-IZLO+1)*REAL_SIZE*nfields*NGHOST

#ifndef TWO_D
#ifndef YPERIODIC
      if ( iprocy > 1 ) then
#endif
           call  MPI_Irecv(ybuffer_recv, YBUFFSIZE, MPI_DOUBLE_PRECISION
     1                 , bottom, MSG_XCH_YHI_TAG, comm3D, 
     &                 msg_id_recv_y_hi,ierr)
           Call ErrorHandler(ierr,ERROR_RECV)
#ifndef YPERIODIC
        endif
#endif
       
#ifndef YPERIODIC
      if ( iprocy < YPROCS ) then
#endif
        
           ybuffer_send(1:imax+1,1:ng,1:kmax+1,:)=
     &        psi(1:imax+1,jmax-ng+1:jmax,1:kmax+1,:)

           call MPI_Isend(ybuffer_send, YBUFFSIZE, MPI_DOUBLE_PRECISION
     1                , top, MSG_XCH_YHI_TAG, comm3D, msg_id_send_y_hi
     &                , ierr)

           Call ErrorHandler(ierr,ERROR_SEND)
#ifndef YPERIODIC
        endif
#endif
       
#ifndef YPERIODIC
      if ( iprocy > 1 ) then
#endif


           call MPI_Wait(msg_id_recv_y_hi, status, ierr)
           Call ErrorHandler(ierr,ERROR_WAIT)
           psi(1:imax+1,0-ng+1:0,1:kmax+1,:)=
     &          ybuffer_recv(1:imax+1,1:ng,1:kmax+1,:)
#ifndef YPERIODIC
        endif
#endif
        
#ifndef YPERIODIC
      if ( iprocy < YPROCS ) then
#endif

           call MPI_Wait(msg_id_send_y_hi, status, ierr)
           Call ErrorHandler(ierr,ERROR_WAIT)
#ifndef YPERIODIC
        endif
#endif

C        update y-high boundaries

#ifndef YPERIODIC
      if ( iprocy < YPROCS ) then
#endif
           call MPI_Irecv(ybuffer_recv, YBUFFSIZE, MPI_DOUBLE_PRECISION,
     1                top, MSG_XCH_YLOW_TAG, comm3D, msg_id_recv_y_low,
     2                ierr)
           Call ErrorHandler(ierr,ERROR_RECV)
#ifndef YPERIODIC
        endif
#endif
       
#ifndef YPERIODIC
      if ( iprocy > 1 ) then
#endif
           ybuffer_send(1:imax+1,1:ng,1:kmax+1,:)=
     &        psi(1:imax+1,1:ng,1:kmax+1,:)

           call MPI_Isend(ybuffer_send, YBUFFSIZE, MPI_DOUBLE_PRECISION,
     1                bottom, MSG_XCH_YLOW_TAG, comm3D,
     &                msg_id_send_y_low, ierr)  
           Call ErrorHandler(ierr,ERROR_SEND)
#ifndef YPERIODIC
        endif
#endif

#ifndef YPERIODIC
      if ( iprocy < YPROCS ) then
#endif
           call MPI_Wait(msg_id_recv_y_low, status, ierr)
           Call ErrorHandler(ierr,ERROR_WAIT)
           psi(1:imax+1,jmax+1:jmax+ng,1:kmax+1,:)=
     &          ybuffer_recv(1:imax+1,1:ng,1:kmax+1,:)
#ifndef YPERIODIC
        endif
#endif

#ifndef YPERIODIC
      if ( iprocy > 1 ) then
#endif
           call MPI_Wait(msg_id_send_y_low, status, ierr)
           Call ErrorHandler(ierr,ERROR_WAIT)
#ifndef YPERIODIC
        endif
#endif
c
c
#endif
C       -------Z DIRECTION COMMUNICATION
c       Z direction is periodic
C        update z-low boundaries
#ifndef ZPERIODIC
      if ( iprocz > 1 ) then
#endif
           call  MPI_Irecv(zbuffer_recv, ZBUFFSIZE, MPI_DOUBLE_PRECISION
     1                 , behind, MSG_XCH_ZHI_TAG, comm3D, 
     &                 msg_id_recv_z_hi, ierr)
           Call ErrorHandler(ierr,ERROR_RECV)
#ifndef ZPERIODIC
      endif
#endif
       
#ifndef ZPERIODIC
      if ( iprocz < ZPROCS ) then
#endif
        
           zbuffer_send(1:imax+1,1:jhi,1:ng,:)=
     &        psi(1:imax+1,1:jhi,kmax-ng+1:kmax,:)

           call MPI_Isend(zbuffer_send, ZBUFFSIZE, MPI_DOUBLE_PRECISION
     1           , forward, MSG_XCH_ZHI_TAG, comm3D, msg_id_send_z_hi
     &           , ierr)

           Call ErrorHandler(ierr,ERROR_SEND)
#ifndef ZPERIODIC
      endif
#endif
       
#ifndef ZPERIODIC
      if ( iprocz > 1 ) then
#endif
           call MPI_Wait(msg_id_recv_z_hi, status, ierr)
           Call ErrorHandler(ierr,ERROR_WAIT)
           psi(1:imax+1,1:jhi,0-ng+1:0,:)=
     &          zbuffer_recv(1:imax+1,1:jhi,1:ng,:)
#ifndef ZPERIODIC
      endif
#endif
        
#ifndef ZPERIODIC
      if ( iprocz < ZPROCS ) then
#endif
           call MPI_Wait(msg_id_send_z_hi, status, ierr)
           Call ErrorHandler(ierr,ERROR_WAIT)
#ifndef ZPERIODIC
      endif
#endif

C        update z-high boundaries

#ifndef ZPERIODIC
      if ( iprocz < ZPROCS ) then
#endif
           call MPI_Irecv(zbuffer_recv, ZBUFFSIZE, MPI_DOUBLE_PRECISION,
     1                forward, MSG_XCH_ZLOW_TAG, comm3D, 
     &                msg_id_recv_z_low, ierr)
           Call ErrorHandler(ierr,ERROR_RECV)
#ifndef ZPERIODIC
      endif
#endif
       
#ifndef ZPERIODIC
      if ( iprocz > 1 ) then
#endif
           zbuffer_send(1:imax+1,1:jhi,1:ng,:)=
     &        psi(1:imax+1,1:jhi,1:ng,:)

           call MPI_Isend(zbuffer_send, ZBUFFSIZE, MPI_DOUBLE_PRECISION,
     1                behind, MSG_XCH_ZLOW_TAG, comm3D, 
     &                msg_id_send_z_low,ierr)  
           Call ErrorHandler(ierr,ERROR_SEND)
#ifndef ZPERIODIC
      endif
#endif

#ifndef ZPERIODIC
      if ( iprocz < ZPROCS ) then
#endif
           call MPI_Wait(msg_id_recv_z_low, status, ierr)
           Call ErrorHandler(ierr,ERROR_WAIT)
           psi(1:imax+1,1:jhi,kmax+1:kmax+ng,:)=
     &          zbuffer_recv(1:imax+1,1:jhi,1:ng,:)
#ifndef ZPERIODIC
      endif
#endif


#ifndef ZPERIODIC
      if ( iprocz > 1 ) then
#endif
           call MPI_Wait(msg_id_send_z_low, status, ierr)
           Call ErrorHandler(ierr,ERROR_WAIT)
#ifndef ZPERIODIC
      endif
#endif
#endif

        return
        end
c
c-----------------------------------------------------------------------
        subroutine GatherPsi3D(psi,imax,jmax,kmax,
     &     psig,imaxg,jmaxg,kmaxg)
c----------------------------------------------------------------------
c       Description: All gather of psi - collected when 
c                    lMaxLocal is reached in multigrid
c                    boundaries for psi
c       Input: psi - field with ghost which must be gathered
c              imax,jmax,kmax: local mesh extent 
c              imaxg,jmaxg,kmaxg: global mesh extent 
c       Output: psig - gathered field with global mesh extent
c-----------------------------------------------------------------
        use mesh_uparms
        use mesh_common
#ifdef PARALLEL
        use mpistuff
#endif
        integer:: imax, jmax, kmax
        integer:: imaxg, jmaxg, kmaxg

        double precision:: psi(0:imax+1,0:jmax+1,0:kmax+1)
        double precision:: psig(0:imaxg+1,0:jmaxg+1,0:kmaxg+1)

#ifdef PARALLEL
        double precision:: buffer_send(imax,jmax,kmax)
        double precision:: buffer_recv(imax,jmax,kmax,
     &                                 XPROCS*YPROCS*ZPROCS)

        integer:: bufferSize
        integer:: recvBufferSize
        integer:: ii,jj,mm,kk,l,idest
        integer:: iprocnum
        integer:: msg_id_send
        integer:: msg_id_recv
        integer:: ipx,ipy,ipz,coords(3)


        bufferSize=imax*jmax*kmax
        buffer_send=psi(1:imax,1:jmax,1:kmax)
        
        call MPI_AllGather(buffer_send,bufferSize,MPI_DOUBLE_PRECISION,
     &   buffer_recv,bufferSize,MPI_DOUBLE_PRECISION,
     &               comm3D,
     &                ierr)
c
        do l=1,XPROCS*YPROCS*ZPROCS,1
         call MPI_Cart_Coords(comm3D, l-1, 3, coords, ierr)
        Call ErrorHandler(ierr,ERROR_CARTCOORDS)
         ipx = coords(1) + 1
         ipy = coords(2) + 1
         ipz = coords(3) + 1           
        psig((ipx-1)*imax+1:ipx*imax,(ipy-1)*jmax+1:ipy*jmax,
     &       (ipz-1)*kmax+1:ipz*kmax)
     &         =buffer_recv(:,:,:,l)
        enddo
#else
        psig=psi
#endif
        return
        end
c-----------------------------------------------------------------------
        subroutine GatherRes3D(res,imax,jmax,kmax,
     &                         resg,imaxg,jmaxg,kmaxg)
c-----------------------------------------------------------------
c       Description: All gather of psi - collected when 
c                    lMaxLocal is reached in multigrid
c                    boundaries for psi
c       Input: res - field with no ghost which must be gathered
c              imax,jmax,kmax: local mesh extent 
c              imaxg,jmaxg,kmaxg: global mesh extent 
c       Output: resg - gathered field with global mesh extent
c-----------------------------------------------------------------
        use mesh_uparms
        use mesh_common
#ifdef PARALLEL
        use mpistuff
#endif

        integer:: imax, jmax, kmax
        integer:: imaxg, jmaxg, kmaxg
        double precision:: res(imax,jmax,kmax)
        double precision:: resg(imaxg,jmaxg,kmaxg)

#ifdef PARALLEL
        double precision:: buffer_send(imax,jmax,kmax)
        double precision:: buffer_recv(imax,jmax,kmax,
     &                              XPROCS*YPROCS*ZPROCS)

        integer:: bufferSize
        integer:: ii,jj,mm,kk,l,idest, i,j
        integer:: iprocnum
        integer:: msg_id_send
        integer:: msg_id_recv
        integer:: ipx,ipy,ipz,coords(3)


        bufferSize=imax*jmax
        buffer_send=res(1:imax,1:jmax,1:kmax)
c
        call MPI_AllGather(buffer_send,bufferSize,MPI_DOUBLE_PRECISION,
     &   buffer_recv,bufferSize,MPI_DOUBLE_PRECISION,
     &               comm3D,
     &                ierr)
c
        do l=1,XPROCS*YPROCS*ZPROCS,1
         call MPI_Cart_Coords(comm3D, l-1, 3, coords, ierr)
        Call ErrorHandler(ierr,ERROR_CARTCOORDS)
         ipx = coords(1) + 1
         ipy = coords(2) + 1
         ipz = coords(3) + 1           
c
        resg((ipx-1)*imax+1:ipx*imax,(ipy-1)*jmax+1:ipy*jmax,
     &       (ipz-1)*kmax+1:ipz*kmax)
     &         =buffer_recv(:,:,:,l)
c
        enddo
#else
        resg=res
#endif
        return
        end
c-----------------------------------------------------------------------
        subroutine ScatterPsi3D(psi,imax,jmax,kmax,
     &                          psig,imaxg,jmaxg,kmaxg)
c-----------------------------------------------------------------
c       Description: All scatter of psi when going back to
c                    a level whence distributed multigrid is continued.
c       Input: psig - field with ghost which must be scattered
c              imax,jmax,kmax: local mesh extent 
c              imaxg,jmaxg,kmaxg: global mesh extent 
c       Output: psi - local field
c-----------------------------------------------------------------
        use mesh_common
        integer:: imax, jmax, kmax
        integer:: imaxg, jmaxg, kmaxg
        double precision:: psi(0:imax+1,0:jmax+1,0:kmax+1)
        double precision:: psig(0:imaxg+1,0:jmaxg+1,0:kmaxg+1)

        psi(1:imax,1:jmax,1:kmax)=
     &   psig((iprocx-1)*imax+1:iprocx*imax,
     &        (iprocy-1)*jmax+1:iprocy*jmax,
     &        (iprocz-1)*kmax+1:iprocz*kmax)
        return
        end
c-----------------------------------------------------------------------
        subroutine ScatterRes3D(res,imax,jmax,kmax,
     &                          resg,imaxg,jmaxg,kmaxg)
c-----------------------------------------------------------------
c       Description: All scatter of res when going back to
c                    a level whence distributed multigrid is continued.
c       Input: resg - field with no ghost which must be scattered
c              imax,jmax,kmax: local mesh extent 
c              imaxg,jmaxg,kmaxg: global mesh extent 
c       Output: res - local field
c-----------------------------------------------------------------
        use mesh_common
        integer:: imax, jmax, kmax
        integer:: imaxg, jmaxg, kmaxg
        double precision:: res(imax,jmax,kmax)
        double precision:: resg(imaxg,jmaxg,kmaxg)

        res(1:imax,1:jmax,1:kmax)=
     &   resg((iprocx-1)*imax+1:iprocx*imax,
     &        (iprocy-1)*jmax+1:iprocy*jmax,
     &        (iprocz-1)*kmax+1:iprocz*kmax)
        return
        end
c--------------------------------------------------------------
      subroutine MaxNorm3D(res,imax,jmax,kmax,neqn,resMax)
c-----------------------------------------------------------------
c       Description: Compute max norm of res
c       Input: res - field with no ghost 
c              imax,jmax,kmax: local mesh extent 
c       Output: resMax - max norm
c-----------------------------------------------------------------
        use mesh_uparms
        use mesh_common
#ifdef PARALLEL
        use mpistuff
#endif
        integer:: imax, jmax, kmax,neqn
        double precision:: res(imax,jmax,kmax,neqn)
        double precision:: resMax
        double precision:: resMaxLocal
#ifdef PARALLEL
        resMaxLocal=maxval(dabs(res))
        call MPI_Allreduce(resMaxLocal, resMax, 1, MPI_DOUBLE_PRECISION,
     &       MPI_MAX, comm3D, ierr)
        Call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
        resMax=maxval(dabs(res))
#endif
        return
        end

c--------------------------------------------------------------
      subroutine L2Norm3D(res,imax,jmax,kmax,neqn,resMax)
c-----------------------------------------------------------------
c       Description: Compute max norm of res
c       Input: res - field with no ghost 
c              imax,jmax,kmax: local mesh extent 
c       Output: resMax - max norm
c-----------------------------------------------------------------
        use mesh_uparms
        use mesh_common
#ifdef PARALLEL
        use mpistuff
#endif
        integer:: imax, jmax, kmax,neqn
#ifdef TWO_D
        double precision:: res(0:imax+1,0:jmax+1,kmax,neqn)
#else
        double precision:: res(0:imax+1,0:jmax+1,0:kmax+1,neqn)
#endif
        double precision:: resMax
        double precision:: resMaxLocal
#ifdef PARALLEL
        resMaxLocal=sum(res(1:imax,1:jmax,1:kmax,:)**2)
        call MPI_Allreduce(resMaxLocal, resMax, 1, MPI_DOUBLE_PRECISION,
     &       MPI_SUM, comm3D, ierr)
        Call ErrorHandler(ierr,ERROR_ALLREDUCE)
        resMax=dsqrt(resmax/(imax*jmax*kmax*xprocs*yprocs*zprocs*neqn))
#else
        resMaxLocal=sum(res(1:imax,1:jmax,1:kmax,:)**2)
        resMax=dsqrt(resmaxLocal/(imax*jmax*kmax*neqn))
#endif
        return
        end
c-----------------------------------------------------------------------
      double precision function InnerProd(psiT,psi,imax,jmax,kmax,neqn)
#ifdef PARALLEL
        use mpistuff
#endif
      implicit none
      integer:: imax,jmax,kmax,neqn
#ifdef TWO_D
      double precision:: psi(0:imax+1,0:jmax+1,1,neqn)
      double precision:: psiT(0:imax+1,0:jmax+1,1,neqn)
#else
      double precision:: psi(0:imax+1,0:jmax+1,0:kmax+1,neqn)
      double precision:: psiT(0:imax+1,0:jmax+1,0:kmax+1,neqn)
#endif
        double precision:: innerProdLocal, tmpInnerProd
#ifdef PARALLEL
        innerProdLocal=sum(psiT(1:imax,1:jmax,1:kmax,:)*
     &   psi(1:imax,1:jmax,1:kmax,:))
        call MPI_Allreduce(innerProdLocal, tmpInnerProd, 
     &       1, MPI_DOUBLE_PRECISION,
     &       MPI_SUM, comm3D, ierr)
        Call ErrorHandler(ierr,ERROR_ALLREDUCE)
        innerProd=tmpInnerProd
#else
      InnerProd=sum(psiT(1:imax,1:jmax,1:kmax,:)*
     &   psi(1:imax,1:jmax,1:kmax,:))   
#endif
      return
      end
