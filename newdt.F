c       Ravi Samtaney
c       Copyright 1997

       subroutine newdt (dt,ux)
c
       use mesh
       use mesh_common
       use properties
#ifdef PARALLEL
       use mpistuff
#endif
       implicit none
        double precision:: ux(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI, nvar)
       double precision:: xspmax, yspmax, zspmax
c       double precision:: u,v,w
       double precision:: dt, dt_local
       double precision:: cfl
c       integer:: i,j,k

c       compute maximum velocity and adjust time step to satisfy
c       cfl criterion

c       cfl=0.45
       cfl=0.5D0
c       cfl=0.225
c       cfl=0.1
c       dt_local = -1.0e0
c       xspmax=1.D16
       yspmax=1.D99
       zspmax=1.D99
       dt=1.D-04
c       return

c       do k = 1,nzlocal
c           do j = 1,nylocal
c              do i = 1,nxlocal+1
c                 u = ux(i,j,k,1)
c                 xspmax = dmin1(xspmax,dx/(max1(dabs(u),1.D-08)))
c              enddo
c           enddo
c        enddo
       xspmax=dx/
     &      max(maxval(ux(1:nxlocal+1,1:nylocal,1:nzlocal,1)),1.D-08)
#ifndef TWO_D
       yspmax=dz/
     &      max(maxval(ux(1:nxlocal,1:nylocal+1,1:nzlocal,2)),1.D-08)
#endif
       zspmax=dz/
     &      max(maxval(ux(1:nxlocal,1:nylocal,1:nzlocal+1,3)),1.D-08)

#ifndef TWO_D
       dt_local=min(xspmax,yspmax,zspmax)
#else
       dt_local=min(xspmax,zspmax)
#endif
c        write(6,*) 'dtx=',dt_local,xspmax,iprocx,iprocz
c     &      max(maxval(ux(1:nxlocal,1:nylocal,1:nzlocal+1,3)),1.D-08)
c        dt=dt_local
#ifdef PARALLEL
        call MPI_Allreduce(dt_local, dt, 1, MPI_DOUBLE_PRECISION,
     &       MPI_MIN, comm3D, ierr)
        dt = dt* cfl
        dt = 0.08d0 
#else   
       dt=dt_local*cfl
c        dt = 0.2d0
#endif  
       return
       end
