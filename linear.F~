c-----------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2011
c     King Abdullah University of Science and Technology
c     All Rights Reserved
c-----------------------------------------------------------------
c     $Log: linear.F,v $
c     Revision 1.5  2011/11/05 13:55:06  chengw
c     *** empty log message ***
c
c     Revision 1.4  2011/10/19 15:14:49  samtanr
c     Added 2D option via TWO_D defines.
c
c     Revision 1.3  2011/10/07 13:19:45  samtanr
c     Argument uxp added to subroutine call.
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
       subroutine  ComputeLinearTerms(Lu,gradp,ux,uxp,Adu,px,
     &                                alpha,beta,dt,timestep,substep )
c-----------------------------------------------------------------
       use mesh
       use mesh_common
       use properties
#ifdef LES
       use virtualwallBC
#endif
#ifdef PARALLEL
       use mpistuff
#endif
       use VelocityFieldSolver
       use PressureFieldSolver
c
       implicit none
c=======Declarations=========
       double precision:: ux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
       double precision:: px(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI)
       double precision:: dt
c
#ifndef TWO_D
       double precision:: Adu(nxlocal+1,nylocal+1,nzlocal+1,nvar)
       double precision:: Lu(nxlocal+1,nylocal+1,nzlocal+1,nvar)
       double precision:: gradp(nxlocal+1,nylocal+1,nzlocal+1,nvar)
#else
       double precision:: Adu(nxlocal+1,nylocal,nzlocal+1,nvar)
       double precision:: Lu(nxlocal+1,nylocal,nzlocal+1,nvar)
       double precision:: gradp(nxlocal+1,nylocal,nzlocal+1,nvar)
#endif
c
       double precision:: uxp(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
       double precision:: utmp(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c
c     
       double precision:: alpha,beta
c
       double precision:: divu(nxlocal,nylocal,nzlocal)
c
#ifndef TWO_D
       double precision:: Rhsu(1:nxlocal+1,nylocal+1,nzlocal+1,3)
       double precision:: uu(-2:nxlocal+3,-2:nylocal+3,-2:nzlocal+3,3)
       double precision:: press(-2:nxlocal+3,-2:nylocal+3,-2:nzlocal+3)
#else
       double precision:: Rhsu(1:nxlocal+1,nylocal,nzlocal+1,3)
       double precision:: uu(-2:nxlocal+3,nylocal,-2:nzlocal+3,3)
       double precision:: press(-2:nxlocal+3,nylocal,-2:nzlocal+3)
#endif

c
c       double precision:: pbdry(nxlocal,nylocal)
c     
       integer:: i,j,k, timestep, substep
       double precision:: aa,bb,cc
       double precision:: lap
c
       aa=16.D0/12.D0
       bb=-30.D0/12.D0
       cc=-1.D0/12.D0
       uu=0.D0
       Rhsu = 0.d0
c
       call SetboundaryValuesST(ux)
       call exchange(ux,3)
c     Solve x-momentum equation
c     u-velocity rhs
       do k=1,nzlocal,1
          do j=1,nylocal,1
             do i=1,nxlocal+1,1
                lap=(cc*(ux(i-2,j,k,1)+ux(i+2,j,k,1))+
     &               aa*(ux(i-1,j,k,1)+ux(i+1,j,k,1))
     &               +bb*ux(i,j,k,1))*dxi*dxi+
#ifndef TWO_D
     &               (cc*(ux(i,j-2,k,1)+ux(i,j+2,k,1))+
     &               aa*(ux(i,j-1,k,1)+ux(i,j+1,k,1))
     &               +bb*ux(i,j,k,1))*dyi*dyi+
#endif
     &               (cc*(ux(i,j,k-2,1)+ux(i,j,k+2,1))+
     &               aa*(ux(i,j,k-1,1)+ux(i,j,k+1,1))
     &               +bb*ux(i,j,k,1))*dzi*dzi
                Rhsu(i,j,k,1)=ux(i,j,k,1)
     &               +dt*alpha/Reynolds*lap-dt*Adu(i,j,k,1)
c     Initial guess for uu
                uu(i,j,k,1)=ux(i,j,k,1)
             enddo
          enddo
       enddo
#ifndef TWO_D
c     v-velocity rhs
       do k=1,nzlocal,1
          do j=1,nylocal+1,1
             do i=1,nxlocal,1
                lap=(cc*(ux(i-2,j,k,2)+ux(i+2,j,k,2))+
     &               aa*(ux(i-1,j,k,2)+ux(i+1,j,k,2))
     &               +bb*ux(i,j,k,2))*dxi*dxi+
     &               (cc*(ux(i,j-2,k,2)+ux(i,j+2,k,2))+
     &               aa*(ux(i,j-1,k,2)+ux(i,j+1,k,2))
     &               +bb*ux(i,j,k,2))*dyi*dyi+
     &               (cc*(ux(i,j,k-2,2)+ux(i,j,k+2,2))+
     &               aa*(ux(i,j,k-1,2)+ux(i,j,k+1,2))
     &               +bb*ux(i,j,k,2))*dzi*dzi
                Rhsu(i,j,k,2)=ux(i,j,k,2)
     &               +dt*alpha/Reynolds*lap-dt*Adu(i,j,k,2)
c     Initial guess for uu
                uu(i,j,k,2)=ux(i,j,k,2)
             enddo
          enddo
       enddo
#endif
c     w-velocity rhs
       do k=1,nzlocal+1,1
          do j=1,nylocal,1
             do i=1,nxlocal,1
                lap=(cc*(ux(i-2,j,k,3)+ux(i+2,j,k,3))+
     &               aa*(ux(i-1,j,k,3)+ux(i+1,j,k,3))
     &               +bb*ux(i,j,k,3))*dxi*dxi+
#ifndef TWO_D
     &               (cc*(ux(i,j-2,k,3)+ux(i,j+2,k,3))+
     &               aa*(ux(i,j-1,k,3)+ux(i,j+1,k,3))
     &               +bb*ux(i,j,k,3))*dyi*dyi+
#endif
     &               (cc*(ux(i,j,k-2,3)+ux(i,j,k+2,3))+
     &               aa*(ux(i,j,k-1,3)+ux(i,j,k+1,3))
     &               +bb*ux(i,j,k,3))*dzi*dzi
                Rhsu(i,j,k,3)=ux(i,j,k,3)
     &               +dt*alpha/Reynolds*lap-dt*Adu(i,j,k,3)
c     Initial guess for uu
                uu(i,j,k,3)=ux(i,j,k,3)
             enddo
          enddo
       enddo
c
c       write(6,*) 'RHSU 2',maxval(Rhsu(:,:,:,2)),minval(Rhsu(:,:,:,2))
c       write(6,*) 'UU',maxval(uu(:,:,:,2))
c       write(6,*) 'ADU',maxval(adu(:,:,:,2))
c-----------------------------------les part
c setting boundary 
c#ifdef LES
#ifdef LES
       utmp = ux
       call les_virtualwall(utmp,timestep,substep)
#endif
c-------------------------------------------------
c
       call SolveHelmholtz(uu,Rhsu,beta*dt/Reynolds)
c        write(*,*) uu(1,1,1,1), uu(2,1,1,1), uu(3,1,1,1)
c
       uxp(1:nxlocal+1,1:nylocal,1:nzlocal,1)=
     &      uu(1:nxlocal+1,1:nylocal,1:nzlocal,1)
#ifndef TWO_D
       uxp(1:nxlocal,1:nylocal+1,1:nzlocal,2)=
     &      uu(1:nxlocal,1:nylocal+1,1:nzlocal,2)
#else
       uxp(1:nxlocal,1:nylocal,1:nzlocal,2)=0.D0
#endif
       uxp(1:nxlocal,1:nylocal,1:nzlocal+1,3)=
     &      uu(1:nxlocal,1:nylocal,1:nzlocal+1,3)

       call SetboundaryValuesST(uxp)
       call exchange(uxp,3)
c
       call computeLu(Lu,ux,uxp,alpha,beta)
c
c       call SetPressureBC(ux,uxp,Lu,Adu,dt)
c     
c     Now compute the divergence of uu and use that as rhs for
c     the pressure poisson solver
       call DivergenceVelocity(divu,uxp)
c       call WriteHDF5Filedivu(divu,1)
c
       press=0.D0
#ifndef TWO_D
c       press(-2:nxlocal+3,-2:nylocal+3,-2:nzlocal+3) =
c     &    px(-2:nxlocal+3,-2:nylocal+3,-2:nzlocal+3)
#else
c       press(-2:nxlocal+3,nylocal,-2:nzlocal+3) =
c     &    px(-2:nxlocal+3,nylocal,-2:nzlocal+3)
#endif

c
       call SolvePressure(press,divu)
c
#ifndef TWO_D
c          px(-2:nxlocal+3,-2:nylocal+3,-2:nzlocal+3) =
c     & press(-2:nxlocal+3,-2:nylocal+3,-2:nzlocal+3)
#else
c          px(-2:nxlocal+3,nylocal,-2:nzlocal+3) =
c     & press(-2:nxlocal+3,nylocal,-2:nzlocal+3)
#endif



c
#ifdef LES
       call les_dpdx(press)
#endif
       call CorrectVelocity(uxp,gradp,press)
c
       return
       end subroutine ComputeLinearTerms
c
c-----------------------------------------------------------------------
       subroutine CorrectVelocity(uu,gradp,press)
c-----------------------------------------------------------------
       use mesh
       use properties
#ifdef PARALLEL
       use mpistuff
#endif

c
       implicit none
c=======Declarations=========
c       double precision:: gradp(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c
c       double precision:: uu(-2:nxlocal+3,-2:nylocal+3,-2:nzlocal+3,3)
       double precision:: uu(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c
#ifndef TWO_D
       double precision:: gradp(nxlocal+1,nylocal+1,nzlocal+1,nvar)
       double precision:: press(-2:nxlocal+3,-2:nylocal+3,-2:nzlocal+3)
#else
       double precision:: gradp(nxlocal+1,nylocal,nzlocal+1,nvar)
       double precision:: press(-2:nxlocal+3,nylocal,-2:nzlocal+3)
#endif
c     
       double precision:: aa,bb
       integer:: i,j,k
c     
c     u is defined on i+1/2,j,k
c     v is defined on i, j+1/2,k
c     w is defined on i, j, k+1/2
c
       aa=9.D0/8.D0
       bb=-1.D0/8.D0
c
       do i=1,nxlocal+1,1
          do j=1,nylocal,1
             do k=1,nzlocal,1
                gradp(i,j,k,1)=aa*(press(i,j,k)-press(i-1,j,k))*dxi+
     &               bb*(press(i+1,j,k)-press(i-2,j,k))*dxi3
                uu(i,j,k,1)=uu(i,j,k,1)-gradp(i,j,k,1)
             enddo
          enddo
       enddo
#ifndef TWO_D
       do i=1,nxlocal,1
          do j=1,nylocal+1,1
             do k=1,nzlocal,1
                gradp(i,j,k,2)=aa*(press(i,j,k)-press(i,j-1,k))*dyi+
     &               bb*(press(i,j+1,k)-press(i,j-2,k))*dyi3
                uu(i,j,k,2)=uu(i,j,k,2)-gradp(i,j,k,2)
             enddo
          enddo
       enddo
#else
       gradp(:,:,:,2)=0.D0
#endif
c       return
       do i=1,nxlocal,1
          do j=1,nylocal,1
             do k=1,nzlocal+1,1
                gradp(i,j,k,3)=aa*(press(i,j,k)-press(i,j,k-1))*dzi+
     &               bb*(press(i,j,k+1)-press(i,j,k-2))*dzi3
                uu(i,j,k,3)=uu(i,j,k,3)-gradp(i,j,k,3)
             enddo
          enddo
       enddo
c
c       write(*,*) maxval(dabs(gradp(:,:,:,1))),
c     &   maxval(dabs(gradp(:,:,:,3)))
c
c           call WriteHDF5Filep(gradp,1)
c
       return
       end subroutine CorrectVelocity
c-----------------------------------------------------------------------
       subroutine DivergenceVelocity(divu,uu)
c-----------------------------------------------------------------
       use mesh
       use properties
#ifdef PARALLEL
       use mpistuff
#endif

c
       implicit none
c=======Declarations=========
       double precision:: divu(nxlocal,nylocal,nzlocal)
c
c       double precision:: uu(-2:nxlocal+3,-2:nylocal+3,-2:nzlocal+3,3)
       double precision:: uu(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c     
       double precision:: aa,bb
       integer:: i,j,k
c     
c     u is defined on i+1/2,j,k
c     v is defined on i, j+1/2,k
c     w is defined on i, j, k+1/2
c
       aa=9.D0/8.D0
       bb=-1.D0/8.D0
c
c     cell center i,j is bounded by (i-1/2 j-1/2) = (i,j), 
c                                   (i+1/2,j-1/2) = (i+1,j)
c                                   (i+1/2,j+1/2) = (i+1,j+1),
c                                   (i-1/2,j+1/2) = (i,j+1) 
c     consider d(u)/dx term - these are first evaluated at
c     cell centers 
c       write(6,*) 'IN DIVERG u2=',minval(uu(:,:,:,2)),
c     &      maxval(uu(:,:,:,2))
       do i=1,nxlocal,1
          do j=1,nylocal,1
             do k=1,nzlocal,1
                divu(i,j,k)=aa*(uu(i+1,j,k,1)-uu(i,j,k,1))*dxi+
     &               bb*(uu(i+2,j,k,1)-uu(i-1,j,k,1))*dxi3+
#ifndef TWO_D
     &               aa*(uu(i,j+1,k,2)-uu(i,j,k,2))*dyi+
     &               bb*(uu(i,j+2,k,2)-uu(i,j-1,k,2))*dyi3+
#endif
     &               aa*(uu(i,j,k+1,3)-uu(i,j,k,3))*dzi+
     &               bb*(uu(i,j,k+2,3)-uu(i,j,k-1,3))*dzi3
             enddo
          enddo
       enddo
c             stop
       return
c
      end subroutine DivergenceVelocity
c-----------------------------------------------------------------------
       subroutine  ComputeLu(Lu,ux,uu,alpha,beta)
c-----------------------------------------------------------------
       use mesh
       use mesh_common
       use properties
#ifdef PARALLEL
       use mpistuff
#endif
c
       implicit none
c=======Declarations=========
#ifndef TWO_D
       double precision:: Lu(nxlocal+1,nylocal+1,nzlocal+1,nvar)
#else
       double precision:: Lu(nxlocal+1,nylocal,nzlocal+1,nvar)
#endif
       double precision:: ux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
       double precision:: uu(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c       double precision:: uu(-2:nxlocal+3,-2:nylocal+3,-2:nzlocal+3,3)
c     
       double precision:: aa,bb,cc
       integer:: i,j,k
c
       double precision:: alpha,beta
c
       double precision:: lap1,lap2
c
       aa=16.D0/12.D0
       bb=-30.D0/12.D0
       cc=-1.D0/12.D0
       do k=1,nzlocal,1
          do j=1,nylocal,1
             do i=1,nxlocal+1,1
                lap1=(cc*(ux(i-2,j,k,1)+ux(i+2,j,k,1))+
     &               aa*(ux(i-1,j,k,1)+ux(i+1,j,k,1))
     &               +bb*ux(i,j,k,1))*dxi*dxi+
#ifndef TWO_D
     &               (cc*(ux(i,j-2,k,1)+ux(i,j+2,k,1))+
     &               aa*(ux(i,j-1,k,1)+ux(i,j+1,k,1))
     &               +bb*ux(i,j,k,1))*dyi*dyi+
#endif
     &               (cc*(ux(i,j,k-2,1)+ux(i,j,k+2,1))+
     &               aa*(ux(i,j,k-1,1)+ux(i,j,k+1,1))
     &               +bb*ux(i,j,k,1))*dzi*dzi
                lap2=(cc*(uu(i-2,j,k,1)+uu(i+2,j,k,1))+
     &               aa*(uu(i-1,j,k,1)+uu(i+1,j,k,1))
     &               +bb*uu(i,j,k,1))*dxi*dxi+
#ifndef TWO_D
     &               (cc*(uu(i,j-2,k,1)+uu(i,j+2,k,1))+
     &               aa*(uu(i,j-1,k,1)+uu(i,j+1,k,1))
     &               +bb*uu(i,j,k,1))*dyi*dyi+
#endif
     &               (cc*(uu(i,j,k-2,1)+uu(i,j,k+2,1))+
     &               aa*(uu(i,j,k-1,1)+uu(i,j,k+1,1))
     &               +bb*uu(i,j,k,1))*dzi*dzi
                Lu(i,j,k,1)=alpha*lap1+beta*lap2
             enddo
          enddo
       enddo
c     v-velocity 
#ifndef TWO_D
       do k=1,nzlocal,1
          do j=1,nylocal+1,1
             do i=1,nxlocal,1
                lap1=(cc*(ux(i-2,j,k,2)+ux(i+2,j,k,2))+
     &               aa*(ux(i-1,j,k,2)+ux(i+1,j,k,2))
     &               +bb*ux(i,j,k,2))*dxi*dxi+
     &               (cc*(ux(i,j-2,k,2)+ux(i,j+2,k,2))+
     &               aa*(ux(i,j-1,k,2)+ux(i,j+1,k,2))
     &               +bb*ux(i,j,k,2))*dyi*dyi+
     &               (cc*(ux(i,j,k-2,2)+ux(i,j,k+2,2))+
     &               aa*(ux(i,j,k-1,2)+ux(i,j,k+1,2))
     &               +bb*ux(i,j,k,2))*dzi*dzi
                lap2=(cc*(uu(i-2,j,k,2)+uu(i+2,j,k,2))+
     &               aa*(uu(i-1,j,k,2)+uu(i+1,j,k,2))
     &               +bb*uu(i,j,k,2))*dxi*dxi+
     &               (cc*(uu(i,j-2,k,2)+uu(i,j+2,k,2))+
     &               aa*(uu(i,j-1,k,2)+uu(i,j+1,k,2))
     &               +bb*uu(i,j,k,2))*dyi*dyi+
     &               (cc*(uu(i,j,k-2,2)+uu(i,j,k+2,2))+
     &               aa*(uu(i,j,k-1,2)+uu(i,j,k+1,2))
     &               +bb*uu(i,j,k,2))*dzi*dzi
                Lu(i,j,k,2)=alpha*lap1+beta*lap2
             enddo
          enddo
       enddo
#else
       Lu(:,:,:,2)=0.D0
#endif       
c     w-velocity rhs
       do k=1,nzlocal+1,1
          do j=1,nylocal,1
             do i=1,nxlocal,1
                lap1=(cc*(ux(i-2,j,k,3)+ux(i+2,j,k,3))+
     &               aa*(ux(i-1,j,k,3)+ux(i+1,j,k,3))
     &               +bb*ux(i,j,k,3))*dxi*dxi+
#ifndef TWO_D
     &               (cc*(ux(i,j-2,k,3)+ux(i,j+2,k,3))+
     &               aa*(ux(i,j-1,k,3)+ux(i,j+1,k,3))
     &               +bb*ux(i,j,k,3))*dyi*dyi+
#endif
     &               (cc*(ux(i,j,k-2,3)+ux(i,j,k+2,3))+
     &               aa*(ux(i,j,k-1,3)+ux(i,j,k+1,3))
     &               +bb*ux(i,j,k,3))*dzi*dzi  
              lap2=(cc*(uu(i-2,j,k,3)+uu(i+2,j,k,3))+
     &               aa*(uu(i-1,j,k,3)+uu(i+1,j,k,3))
     &               +bb*uu(i,j,k,3))*dxi*dxi+
#ifndef TWO_D
     &               (cc*(uu(i,j-2,k,3)+uu(i,j+2,k,3))+
     &               aa*(uu(i,j-1,k,3)+uu(i,j+1,k,3))
     &               +bb*uu(i,j,k,3))*dyi*dyi+
#endif
     &               (cc*(uu(i,j,k-2,3)+uu(i,j,k+2,3))+
     &               aa*(uu(i,j,k-1,3)+uu(i,j,k+1,3))
     &               +bb*uu(i,j,k,3))*dzi*dzi
                Lu(i,j,k,3)=alpha*lap1+beta*lap2
             enddo
          enddo
       enddo
c
       end subroutine ComputeLu
c-----------------------------------------------------------------------
c       subroutine  SetPressureBC(ux,uu,Lu,Adu,dt)
c-----------------------------------------------------------------
c       use mesh
c       use mesh_common
c       use properties
c#ifdef PARALLEL
c       use mpistuff
c#endif
c       use PressBC
c
c       implicit none
cc=======Declarations=========
c       double precision:: ux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c       double precision:: uu(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c       double precision:: Adu(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c       double precision:: Lu(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c#ifndef TWO_D
c       double precision:: Lu(nxlocal+1,nylocal+1,nzlocal+1,nvar)
c       double precision:: Adu(nxlocal+1,nylocal+1,nzlocal+1,nvar)
c#else
c       double precision:: Lu(nxlocal+1,nylocal,nzlocal+1,nvar)
c       double precision:: Adu(nxlocal+1,nylocal,nzlocal+1,nvar)
c#endif
c
c       double precision:: dt
c     Local vars
c       double precision:: dudt,dvdt,dwdt
c     
c       integer:: i,j,k
c
c
c       pBCx=0.D0
c       pBCy=0.D0
c       pBCz=0.D0
c       return
cc
c       
c     X-Boundary
c       do k=1,nzlocal,1
c          do j=1,nylocal,1
c             if(iprocx.eq.1) then
c     Left Boundary
c                i=1
c     Not sure whether dudt as follows is really the correct way
c             dudt=(uu(i,j,k,1)-ux(i,j,k,1))/dt
c     Temporarily set dudt=0
c                dudt=0.D0
c                pBCx(j,k,1)=-Adu(i,j,k,1)+Lu(i,j,k,1)+dudt
c             endif
c             if(iprocx.eq.XPROCS) then
c     Right Boundary
c                i=nxlocal+1
c             dudt=(uu(i,j,k,1)-ux(i,j,k,1))/dt
c                dudt=0.D0
c                pBCx(j,k,2)=-Adu(i,j,k,1)+Lu(i,j,k,1)+dudt
c             endif
c          enddo
c       enddo
c     Y-Boundary
c       do k=1,nzlocal,1
c          do i=1,nxlocal,1
c             if(iprocy.eq.1) then
cc     Front Boundary
c                j=1
c     Temporarily set dvdt=0
c                dvdt=0.D0
c                pBCy(i,k,1)=-Adu(i,j,k,2)+Lu(i,j,k,2)+dvdt
c             endif
c     Right Boundary
c             if(iprocy.eq.YPROCS) then
c                j=nylocal+1
c                dvdt=0.D0
c                pBCy(i,k,2)=-Adu(i,j,k,2)+Lu(i,j,k,2)+dvdt
c             endif
c          enddo
c       enddo
cc     Z-Boundary
c       do j=1,nylocal,1
c          do i=1,nxlocal,1
c             if(iprocz.eq.1) then
c     Front Boundary
c                k=1
c     Temporarily set dwdt=0
c                dwdt=0.D0
c                pBCz(i,j,1)=-Adu(i,j,k,3)+Lu(i,j,k,3)+dwdt
c             endif
c     Right Boundary
c             if(iprocz.eq.ZPROCS) then
c                k=nzlocal+1
c                dwdt=0.D0
c                pBCz(i,j,2)=-Adu(i,j,k,3)+Lu(i,j,k,3)+dwdt
c             endif
c          enddo
c       enddo
c
c
c       end subroutine SetPressureBC

