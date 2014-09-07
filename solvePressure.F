c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2011
c     King Abdullah University of Science and Technology
c     All Rights Reserved
c-----------------------------------------------------------------
c     $Log: solvePressure.F,v $
c     Revision 1.7  2011/10/31 12:08:58  samtanr
c     Extrpolated u,v,w at x=xouter.
c     For now set u,p direchlet bc on upper boundary.
c
c     Revision 1.6  2011/10/19 15:14:49  samtanr
c     Added 2D option via TWO_D defines.
c
c     Revision 1.5  2011/10/18 13:03:39  samtanr
c     Commented writing sum of rhsp.
c
c     Revision 1.4  2011/10/08 13:38:22  samtanr
c     Fixed bug in remove zero mode.
c
c     Revision 1.3  2011/10/07 13:20:43  samtanr
c     Cleaned up BC. Removezeromode subroutine added. Changed post/pre
c     smooth iterations to 10 from 4.
c
c     Revision 1.2  2011/10/03 08:18:29  samtanr
c     Intermediate checkin.
c
c     Revision 1.1  2011/07/25 11:32:01  samtanr
c     Original source.
c
c     Wind code: Incompressible Navier-Stokes LES
c
c-----------------------------------------------------------------
c
      module MGParamsp
      save
      integer, parameter:: lmax=3
      integer:: bchFlag
      integer:: xchangeFlag
      end module MGParamsp

      module PressParameters
      save
      type LevelData
        double precision:: delx, dely, delz
      end type levelData
      type(LevelData), allocatable:: ldata(:)
      double precision:: aa,bb,cc,dd
      end module PressParameters
c
c
      module PressureFieldSolver
      contains
       subroutine SolvePressure(press,rhsp)
       use mesh
       use mesh_common
       use properties
       use PressParameters
c       use PressBC
       use mgparamsp
       use iounits
#ifdef PARALLEL
        use mpistuff
#endif
c
       implicit none
c
c       integer:: neqn 
#ifndef TWO_D
       double precision:: press(-2:nxlocal+3,-2:nylocal+3,-2:nzlocal+3)
#else
       double precision:: press(-2:nxlocal+3,nylocal,-2:nzlocal+3)
#endif
       double precision:: rhsp(1:nxlocal,nylocal,nzlocal)
c
c       integer:: i,j,k,l
c
c
       call SetPressParameters(lmax)
c
c       call RemoveZeroMode(rhsp)
c         write(*,*) sum(rhsp),rhsp(3,1,7)
c       rhsp=rhsp-sum(rhsp)/(nxlocal*nylocal*nzlocal)
c          write(*,*)'RHS of p solver :', sum(rhsp)
c          call WriteHDF5Filedivu(rhsp,1)
c          stop

c       write(6,*) 'SUM RHSP=',sum(rhsp)
C     do j=1,nylocal,1
c          write(6,*) 'RHSP j=',j,maxval(rhsp(:,j,:)),
c     &         minval(rhsp(:,j,:))
c       enddo
          
c       rhsp=rhsp-sum(pBCx)*dy*dz-
c     &      sum(pBCy)*dx*dz-
c     &      sum(pBCz)*dx*dy

       call SolveEllipticEquation(press,rhsp,1)
c
       call XBC(press,nxlocal,nylocal,nzlocal,1,dx)
#ifndef TWO_D
       call YBCPeriodic(press,nxlocal,nylocal,nzlocal,1)
#endif
       call ZBC(press,nxlocal,nylocal,nzlocal,1,dz)
       deallocate(ldata)
       return
       end subroutine
c
c-----------------------------------------------------------------------
       subroutine RemoveZeroMode(rhsp)
       use mesh
       use mesh_common
#ifdef PARALLEL
        use mpistuff
#endif
        implicit none
       double precision:: rhsp(1:nxlocal,nylocal,nzlocal)
c
       double precision:: rhsAvg
       double precision:: rhsSumLocal
c       integer:: i,j,k,l
c
       rhsSumLocal=0.D0
#ifdef PARALLEL
        rhsSumLocal=sum(rhsp(1:nxlocal,1:nylocal,1:nzlocal))
        call MPI_Allreduce(rhsSumLocal, rhsAvg, 1, MPI_DOUBLE_PRECISION,
     &       MPI_SUM, comm3D, ierr)
        Call ErrorHandler(ierr,ERROR_ALLREDUCE)
        rhsAvg=rhsAvg/(nx*ny*nz)
#else
        rhsAvg=sum(rhsp(1:nxlocal,1:nylocal,1:nzlocal))
        rhsAvg=rhsAvg/(nx*ny*nz)
#endif        
c        write(6,*) 'RHSP',sum(rhsp)
        rhsp=rhsp-rhsAvg
c        write(6,*) 'RHSP',sum(rhsp)
       return
       end subroutine RemoveZeroMode
c-----------------------------------------------------------------------
      subroutine SetPressParameters(lmax)
c
      use PressParameters
      use mesh
      use properties
      implicit none
      integer:: lmax
c      integer:: i,j,k,imax,jmax,kmax
      integer:: l
c
c     Allocate lmax levels of data
      allocate(ldata(0:lmax))
c
      ldata(0)%delx=dxi
      ldata(0)%dely=dyi
      ldata(0)%delz=dzi
c
c      ldata(0)%delt=dt
      do l=1,lmax,1
         ldata(l)%delx=0.5D0*ldata(l-1)%delx
         ldata(l)%dely=0.5D0*ldata(l-1)%dely
         ldata(l)%delz=0.5D0*ldata(l-1)%delz
c   
      enddo
c
c      aa=16.D0/12.D0
c      bb=-30.D0/12.D0
c      cc=-1.D0/12.D0
c      aa=1.D0
c      bb=-2.D0
c      cc=0.D0
c      aa=1.D0/576.D0
c      bb=-54.D0/576.D0
c      cc=783.D0/576.D0
c      dd=-1460.D0/576.D0
c
      aa=1.D0/576.D0
      bb=-3.D0/32.D0
      cc=87.D0/64.D0
      dd=-365.D0/144.D0
      return
      end subroutine
c-----------------------------------------------------------------------
      subroutine LaplacePress(vx,rhs,imax,jmax,kmax,
     &     neqn,level,xChangeFlag)
c
      use properties
      use PressParameters
c
      implicit none
c
      integer:: imax,jmax,kmax,level,xChangeFlag 
      integer:: neqn
      double precision:: rhs(imax,jmax,kmax)
#ifndef TWO_D
      double precision:: vx(-2:imax+3,-2:jmax+3,-2:kmax+3)
#else
      double precision:: vx(-2:imax+3,jmax,-2:kmax+3)
#endif
c
      double precision:: dxi,dyi,dzi
c     
      integer:: i,j,k !,l
c
      dxi=ldata(level)%delx**2
      dyi=ldata(level)%dely**2
      dzi=ldata(level)%delz**2
c
      do k=1,kmax,1
         do j=1,jmax,1
            do i=1,imax,1
c$$$               rhs(i,j,k)=(cc*(vx(i+2,j,k)+vx(i-2,j,k))
c$$$     &                   +aa*(vx(i+1,j,k)+vx(i-1,j,k))
c$$$     &                   +bb*vx(i,j,k))*dxi+
c$$$     &              (cc*(vx(i,j+2,k)+vx(i,j-2,k))
c$$$     &                   +aa*(vx(i,j+1,k)+vx(i,j-1,k))
c$$$     &                   +bb*vx(i,j,k))*dyi+
c$$$     &              (cc*(vx(i,j,k+2)+vx(i,j,k-2))
c$$$     &                   +aa*(vx(i,j,k+1)+vx(i,j,k-1))
c$$$     &                   +bb*vx(i,j,k))*dzi
               rhs(i,j,k)=(aa*(vx(i+3,j,k)+vx(i-3,j,k))
     &              +bb*(vx(i+2,j,k)+vx(i-2,j,k))
     &              +cc*(vx(i+1,j,k)+vx(i-1,j,k))
     &              +dd*vx(i,j,k))*dxi+
#ifndef TWO_D
     &              (aa*(vx(i,j+3,k)+vx(i,j-3,k))
     &              +bb*(vx(i,j+2,k)+vx(i,j-2,k))
     &              +cc*(vx(i,j+1,k)+vx(i,j-1,k))
     &              +dd*vx(i,j,k))*dyi+
#endif
     &              (aa*(vx(i,j,k+3)+vx(i,j,k-3))
     &              +bb*(vx(i,j,k+2)+vx(i,j,k-2))
     &              +cc*(vx(i,j,k+1)+vx(i,j,k-1))
     &              +dd*vx(i,j,k))*dzi
            enddo
         enddo
      enddo
c
      return
      end subroutine
c     
c-----------------------------------------------------------------------
      subroutine ApplyOperator(press,res,imax,jmax,kmax,neqn,level,dx)
      use PressParameters
      use mgparamsp
      implicit none
      integer:: imax,jmax,kmax,level,neqn
      double precision:: dx(3)
#ifndef TWO_D
      double precision:: press(-2:imax+3,-2:jmax+3,-2:kmax+3)
      double precision:: res(-2:imax+3,-2:jmax+3,-2:kmax+3)
#else
      double precision:: press(-2:imax+3,jmax,-2:kmax+3)
      double precision:: res(-2:imax+3,jmax,-2:kmax+3)
#endif
c
c     Local vars
      double precision:: lapp(imax,jmax,kmax)
c     
      integer:: i,j,k
c
c      write(6,*) 'APPLY OP BEGIN',maxval(press),minval(press)
      if (bcHFlag.eq.0) then
         call BC3D(press,dx,imax,jmax,kmax,xChangeFlag)
      else
         call BCHomogeneous3D(press,imax,jmax,kmax,xChangeFlag)
      endif
      Call ExchangePsi3D(press,imax,jmax,kmax,neqn,xChangeFlag)
c
c      write(6,*) 'APPLY OP BC',maxval(press),minval(press)
      call LaplacePress(press,lapp,imax,jmax,kmax,neqn,level,
     &     xChangeFlag)
c
c      write(6,*) 'APPLY OP PRESS',maxval(press),minval(press)
         res=0.D0
c
         do k=1,kmax,1
            do j=1,jmax,1
               do i=1,imax,1
                  res(i,j,k)=lapp(i,j,k)
               enddo
            enddo
         enddo
c
         return
         end subroutine
c-----------------------------------------------------------------------
      subroutine BCHomogeneous3D(press,imax,jmax,kmax,xChangeFlag)
c-----------------------------------------------------------------------
c     Description: 
c     Homogeneous Dirichlet BCs for psi in x
c     Homogeneous Neumann BCs for psi in y
c     Homogeneous Neumann BCs for psi in z
c     Input: psi - solution
c            imax,jmax,kmax: mesh size at this level
c            xChangeFlag - whether we are doing this on 1 or many procs.
c     Output: psi: Modified at the ghost cells at physical boundaries
c     
c-----------------------------------------------------------------------
      use mesh_uparms
      use mesh_common
       implicit none
      integer:: imax,jmax,kmax
      integer:: xChangeFlag
#ifndef TWO_D
      double precision:: press(-2:imax+3,-2:jmax+3,-2:kmax+3,1)
#else
      double precision:: press(-2:imax+3,jmax,-2:kmax+3,1)
#endif

      call XBCH(press,imax,jmax,kmax,1)
#ifndef TWO_D
      call YBCPeriodic(press,imax,jmax,kmax,1)
#endif
      call ZBCH(press,imax,jmax,kmax,1)
      return
      end subroutine

c-----------------------------------------------------------------------
      subroutine BC3D(press,dx,imax,jmax,kmax,xChangeFlag)
c-----------------------------------------------------------------------
c     Homogeneous Dirichlet BCs for psi in x
c     Homogeneous Neumann BCs for psi in y
c     Important to make the distinction with the 
c     BCHomogeneous routine in case we have different BCs.
c     Input: psi - solution
c     imax,jmax: mesh size at this level
c     xChangeFlag - whether we are doing this on 1 or many procs.
c     Output: psi: Modified at the ghost cells at physical boundaries
c     
c-----------------------------------------------------------------------
      use mesh_uparms
      use mesh_common
       implicit none
      integer:: imax,jmax,kmax
      integer:: xChangeFlag
#ifndef TWO_D
      double precision:: press(-2:imax+3,-2:jmax+3,-2:kmax+3,1)
#else
      double precision:: press(-2:imax+3,jmax,-2:kmax+3,1)
#endif
      double precision:: dx(3)
c      integer:: i,j,k,l
c     
c     
      call XBC(press,imax,jmax,kmax,1,dx(1))
#ifndef TWO_D
      call YBCPeriodic(press,imax,jmax,kmax,1)
#endif
      call ZBC(press,imax,jmax,kmax,1,dx(3))
      return
      end subroutine

c-----------------------------------------------------------------------
       subroutine SolveEllipticEquation(psi,res,neqn)
c-----------------------------------------------------------------
c     Input: res - RHS of the equations L(phi)=res
c     Implicit Input: lMax and lMaxLocal which are determined 
c                     outside of this subroutine, are the
c                     maximum no of global and local levels, resp.
c     Output: phi
c     Multigrid implementation
c     Author: Ravi Samtaney
c     Assumed parameters: No of presmooth and postsmooth iterations 
c                         set to 10.
c                         Convergence tolerance set to 10^-12
c                         Max no. of V-cycles =10. 
c    
c     
c     
       use mesh
       use properties
       use mesh_common
       use iounits
#ifdef PARALLEL
      use mpistuff
#endif
      use MGParamsp
c     
       implicit none
c=======Declarations=========
c        double precision:: dt
c
        integer:: neqn
#ifndef TWO_D
        double precision:: psi(-2:nxlocal+3,-2:nylocal+3,-2:nzlocal+3)
#else
        double precision:: psi(-2:nxlocal+3,nylocal,-2:nzlocal+3)
#endif
        double precision:: res(1:nxlocal,nylocal,nzlocal)
c
c        double precision,allocatable:: psi(:,:,:,:)
c        double precision,allocatable:: psi0(:,:,:,:)
#ifndef TWO_D
        double precision:: psi0(-2:nxlocal+3,-2:nylocal+3,-2:nzlocal+3)
#else
        double precision:: psi0(-2:nxlocal+3,nylocal,-2:nzlocal+3)
#endif
        double precision:: resx(nxlocal,nylocal,nzlocal,neqn)
c
c        integer:: indx,indy
c        integer:: i,j,k,l
        integer:: maxOuterIter
        integer:: iter, doMG
c
        double precision:: resMax
        double precision:: tolerance
        double precision:: dxl(3)
c
        doMg=1
        maxOuterIter=10
cc
        bcHFlag=0
        xChangeFlag=1
        resMax=1.D16
        tolerance=1.D-7
        iter=1
c     Set initial mesh spacing here
        dxl(1)=dx
        dxl(2)=dy
        dxl(3)=dz
c
        resx=0.D0
c
c        write(6,*) 'Res',minval(res),maxval(res),neqn
        call Residual3D(psi,res,resx,nxlocal,nylocal,nzlocal,neqn,0,dxl)
        call MaxNorm3D(resx,nxlocal,nylocal,nzlocal,neqn,resMax)
c
c        if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
c           write(6,*) 'Initial Pressure Max Norm 3D Residual=',resMax
c        endif
c
        call BC3D(psi,dxl,nxlocal,nylocal,nzlocal,0)
        psi0=psi
c        write(6,*) 'Sum psi', sum(psi), sum(psi0)
c
        if(doMG.eq.1) then
        do while(iter.le.maxOuterIter.and.resMax.gt.tolerance)
           call MGSolve3D(psi,psi0,res,nxlocal,nylocal,nzlocal,neqn,
     &          0,dxl)
c           call BC3D(psi,dxl,nxlocal,nylocal,nzlocal,0)
c           write(6,*) 'Sum psi', sum(psi), sum(psi0)
           psi0=psi
           resx=0.D0
           call Residual3D(psi,res,resx,nxlocal,nylocal,nzlocal,neqn,
     &          0,dxl)
           call MaxNorm3D(resx,nxlocal,nylocal,nzlocal,neqn,resMax)
c           if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
c              write(6,*) 'Max Pressure Norm 3D Residual=',iter,resMax
c           endif
           iter=iter+1
        enddo
c     else use G-S
        else
           call GaussSeidel(psi0,psi,res,50000,nxlocal,nylocal,nzlocal,
     &          neqn, 0,dxl)
        endif
c
        call Residual3D(psi,res,resx,nxlocal,nylocal,nzlocal,neqn,0,dxl)
        call MaxNorm3D(resx,nxlocal,nylocal,nzlocal,neqn,resMax)
        if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
           write(ihis,*) 'Final Pressure Max Norm 3D Residual=',
     &          iter-1,resMax
c           write(6,*) 'Final Pressure Max Norm 3D Residual=',
c     &          iter-1,resMax
        endif
        return
        end subroutine
c     

c     
c-----------------------------------------------------------------------
      recursive subroutine MGSolve3D(psi,psiT,res,
     &     imax,jmax,kmax,neqn,
     &     level,dx)
      use MGParamsp
      use mesh_uparms
      use mesh_common
      implicit none
      integer:: imax,jmax,kmax,neqn,level
      integer:: imaxg,jmaxg,kmaxg
#ifndef TWO_D
      double precision:: psi(-2:imax+3,-2:jmax+3,-2:kmax+3,neqn)
      double precision:: psiT(-2:imax+3,-2:jmax+3,-2:kmax+3,neqn)
#else
      double precision:: psi(-2:imax+3,jmax,-2:kmax+3,neqn)
      double precision:: psiT(-2:imax+3,jmax,-2:kmax+3,neqn)
#endif
      double precision:: res(imax,jmax,kmax,neqn)
      double precision:: resCorr(imax,jmax,kmax,neqn)
      double precision, allocatable:: psiCoarse(:,:,:,:)
      double precision, allocatable:: psiTCoarse(:,:,:,:)
      double precision, allocatable:: resCorrCoarse(:,:,:,:)
      double precision, allocatable:: resg(:,:,:,:)
      double precision, allocatable:: psig(:,:,:,:), psiTg(:,:,:,:)
      double precision:: dx(3)
      double precision:: dxc(3)
      integer:: maxPreSmoothIterations
      integer:: maxPostSmoothIterations
      integer:: gammaK
      integer:: iterGamma
      integer:: imaxc,jmaxc,kmaxc
      maxPreSmoothIterations=10
      maxPostSmoothIterations=10
      gammaK=1
      
      
      if(level.eq.lMax) then
         if(level.eq.0) bcHFlag=0
         call BottomSolve3D(psiT,psi,res,imax,jmax,kmax,neqn,level,dx)
c
      else
         call Smooth3D(psiT,psi,res,maxPreSmoothIterations,
     &        imax,jmax,kmax,neqn,level,dx)
         call Residual3D(psi,res,resCorr,imax,jmax,kmax,neqn,level,dx)
         imaxg=(imax)
#ifndef TWO_D
         jmaxg=(jmax)
#else
         jmaxg=1
#endif
         kmaxg=(kmax)
         imaxc=imaxg/2
#ifndef TWO_D
         jmaxc=jmaxg/2
#else
         jmaxc=1
#endif
         kmaxc=kmaxg/2
c
#ifndef TWO_D
         allocate(psig(-2:imaxg+3,-2:jmaxg+3,-2:kmaxg+3,neqn))
         allocate(psiTg(-2:imaxg+3,-2:jmaxg+3,-2:kmaxg+3,neqn))
         allocate(resg(imaxg,jmaxg,kmaxg,neqn))
         allocate(resCorrCoarse(imaxc,jmaxc,kmaxc,neqn))
         allocate(psiCoarse(-2:imaxc+3,-2:jmaxc+3,-2:kmaxc+3,neqn))
         allocate(psiTCoarse(-2:imaxc+3,-2:jmaxc+3,-2:kmaxc+3,neqn))
#else
         allocate(psig(-2:imaxg+3,jmaxg,-2:kmaxg+3,neqn))
         allocate(psiTg(-2:imaxg+3,jmaxg,-2:kmaxg+3,neqn))
         allocate(resg(imaxg,jmaxg,kmaxg,neqn))
         allocate(resCorrCoarse(imaxc,jmaxc,kmaxc,neqn))
         allocate(psiCoarse(-2:imaxc+3,jmaxc,-2:kmaxc+3,neqn))
         allocate(psiTCoarse(-2:imaxc+3,jmaxc,-2:kmaxc+3,neqn))
#endif
c
         resg=resCorr
         psig=psi
         psiTg=psiT
         call Restrict3D(resg,resCorrCoarse,
     &           imaxg,jmaxg,kmaxg,neqn)
         psiTCoarse=0.D0
         psiCoarse=0.D0
         bcHFlag=1
         dxc=2.D0*dx
c
         do iterGamma=1,gammaK,1
               call MGSolve3D(psiCoarse,psiTCoarse,resCorrCoarse,
     &              imaxc,jmaxc,kmaxc,neqn,
     &              level+1,dxc)
               psiTCoarse=psiCoarse
         enddo
         resCorr=0.D0
         call Prolong3D(resg,imaxg,jmaxg,kmaxg,neqn, psiCoarse)
c
            resCorr=resg
            psi=psig
            psiT=psiTg
c
         psi(1:imax,1:jmax,1:kmax,:)=psi(1:imax,1:jmax,1:kmax,:)+
     &        resCorr(1:imax,1:jmax,1:kmax,:)
         bcHFlag=1
         psiT=psi
c
         if(level.eq.0) bcHFlag=0
         call Smooth3D(psi,psi,res,maxPostSmoothIterations,
     &        imax,jmax,kmax,neqn,level,dx)
c
         deallocate(psig)
         deallocate(psiTg)
         deallocate(resg)
         deallocate(resCorrCoarse)
         deallocate(psiCoarse)
         deallocate(psiTCoarse)
      endif
c
      return
      end subroutine
c-----------------------------------------------------------------------
      subroutine Smooth3D(psiT,psi,res,maxIterations,
     &     imax,jmax,kmax,neqn,
     &     level,dx)
c-----------------------------------------------------------------------
c     Input: psiT - psi_tilde
c     res  - residual
c     maxIterations - max no. of smoothing iterations
c     imax,jmax: mesh size
c     level - multi-grid level (not used here)
c     dx - mesh spacing
c     Implicit input: xChangeFlag  =1,0. Do (not) Exchange ghost 
c     information at this level if this is 1 (0).
c     bchflag = 0,1. Use (do not use)  homogeneous 
c     boundary conditions if this is  1 (0)
c     Output: psi - after smoothing of psiT
c     Description: Uses Gauss-Siedel red-black ordering
c-----------------------------------------------------------------------
      use mesh_uparms
      use mesh_common
      use MGParamsp
       implicit none
       integer:: imax, jmax,kmax,neqn,level,maxIterations
#ifndef TWO_D
      double precision:: psi(-2:imax+3,-2:jmax+3,-2:kmax+3,neqn)
      double precision:: psiT(-2:imax+3,-2:jmax+3,-2:kmax+3,neqn)
#else
      double precision:: psi(-2:imax+3,jmax,-2:kmax+3,neqn)
      double precision:: psiT(-2:imax+3,jmax,-2:kmax+3,neqn)
#endif
      double precision:: res(imax,jmax,kmax,neqn)
c     
c      integer:: i,j,k,iter
      double precision:: dx(3)
c     
      call GaussSeidel(psiT,psi,res,maxIterations,
     &        imax,jmax,kmax,neqn,level,dx)
      return
c
      end subroutine
c
c-----------------------------------------------------------------------
      subroutine BottomSolve3D(psiT,psi,res,
     &        imax,jmax,kmax,neqn,level,dx)
c-----------------------------------------------------------------------
c     Description: Used when mesh cannot be coarsened anymore. 
c                  Still uses Gauss-Siedel Red-Black for
c                  100 iterations. 
c     Issues: 100 is hard -coded but seems to work for the reconnection
c             problem. A better bottom solver (e.g. biconjugate gradient)
c             would be desirable.
c     Input: psiT - psi_tilde
c            res - residual
c            imax,jmax: local mesh size
c            level - multi-grid level
c            dx - local mesh spacing
c     Output: psi
c-----------------------------------------------------------------------
       implicit none
       integer:: imax,jmax,kmax,neqn,level
#ifndef TWO_D
      double precision:: psi(-2:imax+3,-2:jmax+3,-2:kmax+3,neqn)
      double precision:: psiT(-2:imax+3,-2:jmax+3,-2:kmax+3,neqn)
#else
      double precision:: psi(-2:imax+3,jmax,-2:kmax+3,neqn)
      double precision:: psiT(-2:imax+3,jmax,-2:kmax+3,neqn)
#endif
      double precision:: res(imax,jmax,kmax,neqn)
      double precision:: dx(3)
c      integer:: i,j,k
c
      call GaussSeidel(psiT,psi,res,1000,
     &        imax,jmax,kmax,neqn,level,dx)
      return
      end subroutine

c-----------------------------------------------------------------------
      Subroutine Residual3D(psi,res,resCorr,imax,jmax,kmax,neqn,
     &     level,dx)
c-----------------------------------------------------------------------
c     Description: Calculate the residual Residual=L(phi)- RHS
c     Input: psi 
c            res - rhs 
c            imax,jmax: local mesh size
c            dx - local mesh spacing
c     Output: resCorr  evualated as resCorr=L(psi)-res
c-----------------------------------------------------------------------
      implicit none
      integer:: imax,jmax,kmax,neqn,level
#ifndef TWO_D
      double precision:: psi(-2:imax+3,-2:jmax+3,-2:kmax+3,neqn)
#else
      double precision:: psi(-2:imax+3,jmax,-2:kmax+3,neqn)
#endif
      double precision:: resCorr(imax,jmax,kmax,neqn)
      double precision:: res(imax,jmax,kmax,neqn)
      double precision:: dx(3)
c
      integer:: i,j,k
#ifndef TWO_D
      double precision:: resTmp(-2:imax+3,-2:jmax+3,-2:kmax+3,neqn)
#else
      double precision:: resTmp(-2:imax+3,jmax,-2:kmax+3,neqn)
#endif
c     
c
c      write(6,*) 'RES',imax,jmax,kmax,neqn,dx
      resCorr=0.D0
      call ApplyOperator(psi,restmp,imax,jmax,kmax,neqn,level,dx)
      do k=1,kmax,1
      do j=1,jmax,1
         do i=1,imax,1
            resCorr(i,j,k,:)=res(i,j,k,:)
     &           -resTmp(i,j,k,:)
         enddo
      enddo
      enddo
c
      return
      end subroutine
c-----------------------------------------------------------------------
      Subroutine Restrict3D(psi,psiCoarse,imax,jmax,kmax,neqn)
c-----------------------------------------------------------------------
c     Description: Restrict the solution onto the next coarser level
c     Input: psi - solution at current level
c            imax,jmax: mesh size at the current level
c     Output: phiCoarse: solution on the next coarser level
c-----------------------------------------------------------------------
       implicit none
      integer:: imax,jmax,kmax,neqn
      double precision:: psi(imax,jmax,kmax,neqn)
#ifndef TWO_D
      double precision:: psiCoarse(imax/2,jmax/2,kmax/2,neqn)
#else
      double precision:: psiCoarse(imax/2,jmax,kmax/2,neqn)
#endif
c
      integer:: i,j,k
      integer:: ii,jj,kk
c
#ifndef TWO_D
      do k=1,kmax/2,1
         kk=(k-1)*2+1
         do j=1,jmax/2,1
            jj=(j-1)*2+1
            do i=1,imax/2,1
               ii=(i-1)*2+1
               psiCoarse(i,j,k,:)=1.D0/8.D0*(
     &              psi(ii,jj,kk,:)+psi(ii+1,jj,kk,:)
     &              +psi(ii+1,jj+1,kk,:)+psi(ii,jj+1,kk,:)+
     &              psi(ii,jj,kk+1,:)+psi(ii+1,jj,kk+1,:)
     &              +psi(ii+1,jj+1,kk+1,:)+psi(ii,jj+1,kk+1,:))
            enddo
         enddo
      enddo
#else
      jj=1
      j=1
      do k=1,kmax/2,1
         kk=(k-1)*2+1
            do i=1,imax/2,1
               ii=(i-1)*2+1
               psiCoarse(i,j,k,:)=1.D0/4.D0*(
     &              psi(ii,jj,kk,:)+psi(ii+1,jj,kk,:)
     &              +psi(ii,jj,kk+1,:)+psi(ii+1,jj,kk+1,:))
            enddo
      enddo
#endif

      return
      end subroutine
c-----------------------------------------------------------------------
      subroutine  Prolong3D(psi,imax,jmax,kmax, neqn,psiCoarse)
c-----------------------------------------------------------------------
c     Description: Prolong the solution onto the next finer level
c     Input: psiCoarse - solution at current level
c            imax,jmax: mesh size at the next finer level
c     Output: phi: solution prolonged on the next finer level
c     Issues:  Is this second order accurate? Do we care, because the 
c              entire MG solver is converging well?
c-----------------------------------------------------------------------
       implicit none
      integer:: imax,jmax,kmax,neqn
      double precision:: psi(imax,jmax,kmax,neqn)
#ifndef TWO_D
      double precision::
     &     psiCoarse(-2:imax/2+3,-2:jmax/2+3,-2:kmax/2+3,neqn)
#else
      double precision::
     &     psiCoarse(-2:imax/2+3,jmax,-2:kmax/2+3,neqn)
#endif
      integer:: i,j,k
      integer:: ii,jj,kk
#ifndef TWO_D
      do k=1,kmax/2,1
         kk=(k-1)*2+1
      do j=1,jmax/2,1
         jj=(j-1)*2+1
         do i=1,imax/2,1
            ii=(i-1)*2+1
            psi(ii,jj,kk,:)=psiCoarse(i,j,k,:)
            psi(ii+1,jj,kk,:)=psiCoarse(i,j,k,:)
            psi(ii+1,jj+1,kk,:)=psiCoarse(i,j,k,:)
            psi(ii,jj+1,kk,:)=psiCoarse(i,j,k,:)
            psi(ii,jj,kk+1,:)=psiCoarse(i,j,k,:)
            psi(ii+1,jj,kk+1,:)=psiCoarse(i,j,k,:)
            psi(ii+1,jj+1,kk+1,:)=psiCoarse(i,j,k,:)
            psi(ii,jj+1,kk+1,:)=psiCoarse(i,j,k,:)
         enddo
      enddo
      enddo
#else
      j=1
      jj=1
      do k=1,kmax/2,1
         kk=(k-1)*2+1
         do i=1,imax/2,1
            ii=(i-1)*2+1
            psi(ii,jj,kk,:)=psiCoarse(i,j,k,:)
            psi(ii+1,jj,kk,:)=psiCoarse(i,j,k,:)
            psi(ii,jj,kk+1,:)=psiCoarse(i,j,k,:)
            psi(ii+1,jj,kk+1,:)=psiCoarse(i,j,k,:)
         enddo
      enddo
#endif
      return
      end subroutine

c-----------------------------------------------------------------------
      Subroutine GaussSeidel(psiT,psi,res,maxIter,imax,jmax,kmax,neqn,
     &     level,dx)
c-----------------------------------------------------------------------
c     Description: Calculate the residual Residual=L(phi)- RHS
c     Input: psiT
c     Output: psi
c     res - rhs 
c     imax,jmax,kmax: local mesh size
c     dx - local mesh spacing
c     Output: resCorr  evaluated as resCorr=L(psi)-res
c-----------------------------------------------------------------------
      use mgparamsp
      use PressParameters
      implicit none
      integer:: imax,jmax,kmax, neqn,level
#ifndef TWO_D
      double precision:: psi(-2:imax+3,-2:jmax+3,-2:kmax+3,neqn)
      double precision:: psiT(-2:imax+3,-2:jmax+3,-2:kmax+3,neqn)
#else
      double precision:: psi(-2:imax+3,jmax,-2:kmax+3,neqn)
      double precision:: psiT(-2:imax+3,jmax,-2:kmax+3,neqn)
#endif
      double precision:: res(imax,jmax,kmax,neqn)
      double precision:: resx(imax,jmax,kmax,neqn)
      double precision:: resMax
c     
      double precision:: dx(3)
      double precision:: dxi,dyi,dzi
c     
      integer:: i,j,k
c     
c      double precision:: fac1,fac2
      double precision:: tolerance
      integer:: iter, MaxIter
c     
c     
      tolerance=1.D-9
      dxi=ldata(level)%delx**2
      dyi=ldata(level)%dely**2
      dzi=ldata(level)%delz**2
c      fac2=1.D0/(1.D0+(2.D0*dxi+2.D0*dyi+2.D0*dzi)*fac1)
      Call ExchangePsi3D(psiT,imax,jmax,kmax,neqn,xChangeFlag)
      psi=psiT
      Call ExchangePsi3D(psi,imax,jmax,kmax,neqn,xChangeFlag)
c     
      iter=1
      resMax=1.D16
      do while(iter.le.MaxIter.and.resMax.gt.tolerance)
c         call BC2D(psi,dx,imax,jmax,xChangeFlag)
         if (bcHFlag.eq.0) then
            call BC3D(psi,dx,imax,jmax,kmax,xChangeFlag)
         else
            call BCHomogeneous3D(psi,imax,jmax,kmax,xChangeFlag)
         endif
c
         do k=1,kmax,1
         do j=1,jmax,1
            do i=1,imax,1
#ifndef TWO_D
               psi(i,j,k,:)=(-res(i,j,k,:)+
     &              (aa*(psi(i+3,j,k,:)+psi(i-3,j,k,:))
     &              +bb*(psi(i+2,j,k,:)+psi(i-2,j,k,:))
     &              +cc*(psi(i+1,j,k,:)+psi(i-1,j,k,:))
     &              )*dxi+
     &              (aa*(psi(i,j+3,k,:)+psi(i,j-3,k,:))
     &              +bb*(psi(i,j+2,k,:)+psi(i,j-2,k,:))
     &              +cc*(psi(i,j+1,k,:)+psi(i,j-1,k,:))
     &              )*dyi+
     &              (aa*(psi(i,j,k+3,:)+psi(i,j,k-3,:))
     &              +bb*(psi(i,j,k+2,:)+psi(i,j,k-2,:))
     &              +cc*(psi(i,j,k+1,:)+psi(i,j,k-1,:))
     &              )*dzi)/
     &            (-dd*(dxi+dyi+dzi))
#else
               psi(i,j,k,:)=(-res(i,j,k,:)+
     &              (aa*(psi(i+3,j,k,:)+psi(i-3,j,k,:))
     &              +bb*(psi(i+2,j,k,:)+psi(i-2,j,k,:))
     &              +cc*(psi(i+1,j,k,:)+psi(i-1,j,k,:))
     &              )*dxi+
     &              (aa*(psi(i,j,k+3,:)+psi(i,j,k-3,:))
     &              +bb*(psi(i,j,k+2,:)+psi(i,j,k-2,:))
     &              +cc*(psi(i,j,k+1,:)+psi(i,j,k-1,:))
     &              )*dzi)/
     &            (-dd*(dxi+dzi))
#endif

            enddo
         enddo
         enddo
         Call ExchangePsi3D(psi,imax,jmax,kmax,neqn,xChangeFlag)
         call Residual3D(psi,res,resx,imax,jmax,kmax,neqn,level,dx)
         call MaxNorm3D(resx,imax,jmax,kmax,neqn,resMax)
c         write(6,*) 'Level',level,'Max Norm Residual=',iter,resMax
c
c     
         iter=iter+1
      enddo
c       do j=1,jmax,1
c          write(6,*) 'IN GS j=',j,minval(psi(:,j,:,1)),
c     &      maxval(psi(:,j,:,1))
c       enddo
      call Residual3D(psi,res,resx,imax,jmax,kmax,neqn,level,dx)
      call MaxNorm3D(resx,imax,jmax,kmax,neqn,resMax)
c      write(6,*) 'Level',level,'Max Norm Residual=',iter,resMax
c     
      return
      end subroutine
c-----------------------------------------------------------------------
        subroutine XBC(psi,imax,jmax,kmax,ncomp,dx)
        use mesh_parms
        use mesh_common
c        use PressBC
        implicit none
        integer:: imax,jmax,kmax,ncomp
#ifndef TWO_D
        double precision:: psi(-2:imax+3,-2:jmax+3,-2:kmax+3,ncomp)
#else
        double precision:: psi(-2:imax+3,jmax,-2:kmax+3,ncomp)
#endif
        double precision:: dx
        integer i,j,k,l
c       
c     When level=0 we call this
        if (iprocx .eq. 1) then
           i=0
           do k=1,kmax,1
              do j=1,jmax,1
                 do l=1,ncomp,1
                    psi(i,j,k,l)=psi(i+1,j,k,l)
                    psi(i-1,j,k,l)=psi(i+2,j,k,l)
                    psi(i-2,j,k,l)=psi(i+3,j,k,l)
c                    psi(i,j,k,l)=psi(i+1,j,k,l)-pBCx(j,k,1)*dx
c                    psi(i-1,j,k,l)=psi(i+2,j,k,l)-2.D0*pBCx(j,k,1)*dx
c                    psi(i-2,j,k,l)=psi(i+3,j,k,l)-3.D0*pBCx(j,k,1)*dx
                 enddo
              enddo
           enddo
        endif
c       
        if (iprocx .eq. xprocs) then
           i=imax+1
           do k=1,kmax,1
                 do j=1,jmax,1
                    do l=1,ncomp,1
c                       psi(i,j,k,l)=psi(i-1,j,k,l)+pBCx(j,k,2)*dx
c                       psi(i+1,j,k,l)=psi(i-2,j,k,l)+pBCx(j,k,2)*dx*2.D0
c                       psi(i+2,j,k,l)=psi(i-3,j,k,l)+pBCx(j,k,2)*dx*3.D0
                       psi(i,j,k,l)=psi(i-1,j,k,l)
                       psi(i+1,j,k,l)=psi(i-2,j,k,l)
                       psi(i+2,j,k,l)=psi(i-3,j,k,l)
                    enddo
                 enddo
              enddo
           endif
c       
        return
        end subroutine
c
c-----------------------------------------------------------------------
        subroutine XBCH(psi,imax,jmax,kmax,ncomp)
        use mesh_parms
        use mesh_common
        implicit none
        integer:: imax,jmax,kmax,ncomp
#ifndef TWO_D
        double precision:: psi(-2:imax+3,-2:jmax+3,-2:kmax+3,ncomp)
#else
        double precision:: psi(-2:imax+3,jmax,-2:kmax+3,ncomp)
#endif
        integer i,j,k,l
c       
c     When level=0 we call this
        if (iprocx .eq. 1) then
           i=0
           do k=1,kmax,1
              do j=1,jmax,1
                 do l=1,ncomp,1
                    psi(i,j,k,l)=psi(i+1,j,k,l)
                    psi(i-1,j,k,l)=psi(i+2,j,k,l)
                    psi(i-2,j,k,l)=psi(i+3,j,k,l)
                 enddo
              enddo
           enddo
        endif
c       
        if (iprocx .eq. xprocs) then
           i=imax+1
           do k=1,kmax,1
                 do j=1,jmax,1
                    do l=1,ncomp,1
                       psi(i,j,k,l)=psi(i-1,j,k,l)
                       psi(i+1,j,k,l)=psi(i-2,j,k,l)
                       psi(i+2,j,k,l)=psi(i-3,j,k,l)
                    enddo
                 enddo
              enddo
           endif
c       
        return
        end subroutine
c
c-----------------------------------------------------------------------
c       YBC periodic
        subroutine YBCPeriodic(psi,imax,jmax,kmax,ncomp)
        use mesh_parms
        use mesh_common
        implicit none
        integer:: imax,jmax,kmax,ncomp
      double precision:: psi(-2:imax+3,-2:jmax+3,-2:kmax+3,ncomp)
        integer i,j,k,l
c       
c     we rely on this for periodic bc if yprocs = 1 otherwise
c     the exchange routine will take care of it.
        if (yprocs .eq. 1) then
           j=0
           do k=1,kmax,1
              do i=1,imax,1
                 do l=1,ncomp,1
                    psi(i,j,k,l)=psi(i,jmax,k,l)
                    psi(i,j-1,k,l)=psi(i,jmax-1,k,l)
                    psi(i,j-2,k,l)=psi(i,jmax-2,k,l)
                 enddo
              enddo
           enddo
           j=jmax+1
           do k=1,kmax,1
              do i=1,imax,1
                 do l=1,ncomp,1
                    psi(i,j,k,l)=psi(i,1,k,l)
                    psi(i,j+1,k,l)=psi(i,2,k,l)
                    psi(i,j+2,k,l)=psi(i,3,k,l)
                 enddo
              enddo
           enddo
        endif
c       
        return
        end subroutine
c-----------------------------------------------------------------------
c       ZBC 
        subroutine ZBC(psi,imax,jmax,kmax,ncomp,dz)
        use mesh_parms
        use mesh_common
        use PressParameters
c        use PressBC
        implicit none
        integer:: imax,jmax,kmax,ncomp
#ifndef TWO_D
        double precision:: psi(-2:imax+3,-2:jmax+3,-2:kmax+3,ncomp)
#else
        double precision:: psi(-2:imax+3,jmax,-2:kmax+3,ncomp)
#endif
        double precision:: dz
        integer i,j,k,l
c       
        if (iprocz .eq. 1) then
           k=0
           do j=1,jmax,1
              do i=1,imax,1
                 do l=1,ncomp,1
c                    psi(i,j,k,l)=psi(i,j,k+1,l)-pBCz(i,j,1)*dz
c                    psi(i,j,k-1,l)=psi(i,j,k+2,l)-pBCz(i,j,1)*dz*2.D0
c                    psi(i,j,k-2,l)=psi(i,j,k+3,l)-pBCz(i,j,1)*dz*3.D0
                    psi(i,j,k,l)=psi(i,j,k+1,l)
                    psi(i,j,k-1,l)=psi(i,j,k+2,l)
                    psi(i,j,k-2,l)=psi(i,j,k+3,l)
                 enddo
              enddo
           enddo
        endif
c       
        if (iprocz .eq. zprocs) then
           k=kmax+1
           do j=1,jmax,1
              do i=1,imax,1
                 do l=1,ncomp,1
c                    psi(i,j,k,l)=psi(i,j,k-1,l)+pBCz(i,j,2)*dz
c                    psi(i,j,k+1,l)=psi(i,j,k-2,l)+pBCz(i,j,2)*dz*2.D0
c                    psi(i,j,k+2,l)=psi(i,j,k-3,l)+pBCz(i,j,2)*dz*3.D0
                    psi(i,j,k,l)=-psi(i,j,k-1,l)
                    psi(i,j,k+1,l)=-psi(i,j,k-2,l)
                    psi(i,j,k+2,l)=-psi(i,j,k-3,l)
                 enddo
              enddo
           enddo
        endif
c       
        return
        end subroutine
c
c
c-----------------------------------------------------------------------
c       ZBC 
        subroutine ZBCH(psi,imax,jmax,kmax,ncomp)
        use mesh_parms
        use mesh_common
        implicit none
        integer:: imax,jmax,kmax,ncomp
#ifndef TWO_D
        double precision:: psi(-2:imax+3,-2:jmax+3,-2:kmax+3,ncomp)
#else
        double precision:: psi(-2:imax+3,jmax,-2:kmax+3,ncomp)
#endif
        integer i,j,k,l
c       
        if (iprocz .eq. 1) then
           k=0
           do j=1,jmax,1
              do i=1,imax,1
                 do l=1,ncomp,1
                    psi(i,j,k,l)=psi(i,j,k+1,l)
                    psi(i,j,k-1,l)=psi(i,j,k+2,l)
                    psi(i,j,k-2,l)=psi(i,j,k+3,l)
                 enddo
              enddo
           enddo
        endif
c       
        if (iprocz .eq. zprocs) then
           k=kmax+1
           do j=1,jmax,1
              do i=1,imax,1
                 do l=1,ncomp,1
                    psi(i,j,k,l)=-psi(i,j,k-1,l)
                    psi(i,j,k+1,l)=-psi(i,j,k-2,l)
                    psi(i,j,k+2,l)=-psi(i,j,k-3,l)
                 enddo
              enddo
           enddo
        endif
c       
        return
        end subroutine
c
c
c-----------------------------------------------------------------------
       end module pressurefieldsolver
