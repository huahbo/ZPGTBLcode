c-----------------------------------------------------------------
c     Ravi Samtaney
c     Mustafa M. Rahman
c     Copyright 2014
c     King Abdullah University of Science and Technology
c     All Rights Reserved
c-----------------------------------------------------------------
c
c     Revision 1.8  2014/05/22 15:30:00  rahmanmm
c     Original source.
c
c     Wind code: Incompressible Navier-Stokes LES
c
c-----------------------------------------------------------------
       subroutine SolidSolve(ux,px,timestep)  !!! Remove the px, ! Replace ux -> wux !!SolidSolve(wux,timestep)
c-----------------------------------------------------------------
       use mesh
       use mesh_common
       use time
       use properties
       use iounits
       use stats
#ifdef PARALLEL
       use mpistuff
#endif
       use RKCoefficients
c
       implicit none
c=======Declarations=========
       double precision:: ux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)	!ux -> !! wux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
       double precision:: px(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI)		!! Remove px
       integer:: timestep						! Can be done
c
c     RK intermediate stages solution
       double precision:: uxp(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)	!uxp -> wuxp(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
       double precision:: uxpp(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)	!uxpp -> wuxpp
       double precision:: utmp(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)	!utmp -> wutmp
c
       double precision:: wforce(nxlocal+1,nylocal+1,nzlocal+1,nvar) !Remove wind turbine force at various locations
c
#ifndef TWO_D
       double precision:: Adu(nxlocal+1,nylocal+1,nzlocal+1,nvar)	!Adu -> Adwu
       double precision:: Adup(nxlocal+1,nylocal+1,nzlocal+1,nvar)	!Adup -> Adwup
       double precision:: Ads(nxlocal+1,nylocal+1,nzlocal+1,nvar)	!Ads -> Ads
c
       double precision:: Lu(nxlocal+1,nylocal+1,nzlocal+1,nvar)	!! Remove 
       double precision:: Lup(nxlocal+1,nylocal+1,nzlocal+1,nvar)	!! Remove 
       double precision:: gradp(nxlocal+1,nylocal+1,nzlocal+1,nvar)	!! Remove 
#else
       double precision:: Adu(nxlocal+1,nylocal,nzlocal+1,nvar)		!Adu -> Adwu
       double precision:: Adup(nxlocal+1,nylocal,nzlocal+1,nvar)	!Adup -> Adwup	
       double precision:: Ads(nxlocal+1,nylocal,nzlocal+1,nvar)		!Ads -> Ads
c
       double precision:: Lu(nxlocal+1,nylocal,nzlocal+1,nvar)
       double precision:: Lup(nxlocal+1,nylocal,nzlocal+1,nvar)
       double precision:: gradp(nxlocal+1,nylocal,nzlocal+1,nvar)
#endif
c

c----------------- Stage before RK : Rescaling inflow-------------------
#ifdef CODEA
       if(ttot.lt. 250.d0) then
           Tave = 10.d0
       endif
       if(ttot.gt. 250.d0 .and. ttot.lt.350.d0)then
          Tave = 100.d0
       endif 
       if(ttot.gt.350.d0) then
          Tave = 100.d0 + ttot - 350.d0
       endif
       call setinflow(ux,timestep)  
#endif 
#ifdef CODEB
       if(iprocx.eq.1) call readinflowB(timestep)
#endif

c----------------- Stage  RK -------------------------------------------
       call SetBoundaryValuesST(ux)
       call exchange(ux,nvar)
       gradp=0.D0
       Lu=0.D0
       Lup=0.D0
       Adu=0.D0
       Adup=0.D0
       Ads=0.D0
       uxp=ux
       uxpp=ux
c
c       write(6,*) 'Fluid solve dt=',dt, beta(1), Reynolds
c
c----------------- Stage 1 of RK-----------------------------------------
       call ComputeAdvectionTerms(Adu,ux)            ! Adu is for intermediate v* i.e. ux     
#ifdef LES
c       call fftnonlinear(Adu)
       utmp = ux
       call sgsnonlinear(Adu,utmp,timestep)	!Remove: No subgridscale model
#endif
       Ads=Adu*gamma(1)
#ifdef TURBINE
       utmp = ux
       call windforce(utmp,wforce,timestep,1)	!Remove wind turbine at various locations
       Ads = Ads + wforce			! Remove this Wind force term by wind turbine.
#endif
       call outflow(ux,1)			! Will consider the outflow BCs for solid later
       call ComputeLinearTerms(Lu,gradp,ux,uxp,Ads,px,	!Remove
     &                         alpha(1),beta(1),dt,timestep,1) !As computes the Linear part of Eq(3)
       call SetBoundaryValuesST(uxp) 		!Later Bounadary Conditions: ground wall + other BCs
       call exchange(uxp,nvar)				!Later for MPI: Parallel computation
c
c
c      check
c       write(6,*)'Min/max Adu',maxval(Adu(:,:,:,1)),minval(Adu(:,:,:,1))
c       write(6,*) 'Min/max Lu',maxval(Lu(:,:,:,1)),minval(Lu(:,:,:,1))
c       write(6,*) 'Min/max gradp',maxval(gradp(:,:,:,1)),
c     &      minval(gradp(:,:,:,1)), maxval(dabs(gradp(:,:,:,3)))
c
c------------------Stage 2 of RK----------------------------------------
       gradp=0.D0
       call ComputeAdvectionTerms(Adup,uxp)    ! Adup is for intermediate v* i.e. uxp
#ifdef LES
c       call fftnonlinear(Adup)
       utmp = uxp
       call sgsnonlinear(Adup,utmp,timestep)	!Remove: No subgridscale model
#endif
       Ads=Adu*zeta(1)+Adup*gamma(2)
#ifdef TURBINE
       utmp = uxp
       call windforce(utmp,wforce,timestep,2)		!Remove wind turbine at various locations
       Ads = Ads + wforce			! Remove this Wind force term by wind turbine.
#endif
       call outflow(uxp,2)			! Will consider the outflow BCs for solid later
       
       call ComputeLinearTerms(Lup,gradp,uxp,uxpp,Ads,px, !Remove
     &                         alpha(2),beta(2),dt,timestep,2)  !As computes the Linear part of Eq(3)
       call SetBoundaryValuesST(uxpp)		!Later Bounadary Conditions: ground wall + other BCs
       call exchange(uxpp,nvar)			!Later for MPI: Parallel computation
c       write(6,*)'Min/max Adu',
c     &      maxval(Adup(:,:,:,1)),minval(Adup(:,:,:,1))
c       write(6,*) 'Min/max Lu',maxval(Lup(:,:,:,1)),minval(Lup(:,:,:,1))
c       write(6,*) 'Min/max gradp',maxval(gradp(:,:,:,1)),
c     &      minval(gradp(:,:,:,1)), maxval(dabs(gradp(:,:,:,3)))
c
c-------------------Stage 3 of RK---------------------------------------
       gradp = 0.d0
       call ComputeAdvectionTerms(Adu,uxpp)    ! Adu is for intermediate v* i.e. uxpp
#ifdef LES
c       call fftnonlinear(Adu)
       utmp = uxpp
       call sgsnonlinear(Adu,utmp,timestep) 	!Remove: No subgridscale model
#endif
       Ads=Adup*zeta(2)+Adu*gamma(3)  
#ifdef TURBINE
       utmp = uxpp
       call windforce(utmp,wforce,timestep,3) 	!Remove wind turbine at various locations
       Ads = Ads + wforce		! Will consider the outflow BCs for solid later
#endif
       call outflow(uxpp,3)
       call ComputeLinearTerms(Lu,gradp,uxpp,ux,Ads,px,  !Remove
     &                         alpha(3),beta(3),dt,timestep,3)  !As computes the Linear part of Eq(3)
       call SetBoundaryValuesST(ux) 		!Later Bounadary Conditions: ground wall + other BCs
       call exchange(ux,nvar)		!Later for MPI: Parallel computation

c       write(6,*)'Min/max Adu',maxval(Adu(:,:,:,1)),minval(Adu(:,:,:,1))
c       write(6,*) 'Min/max Lu',maxval(Lu(:,:,:,1)),minval(Lu(:,:,:,1))
c       write(6,*) 'Min/max gradp',maxval(gradp(:,:,:,1)),
c     &      minval(gradp(:,:,:,1)), maxval(dabs(gradp(:,:,:,3)))
c
c
       call checkdivergence(ux)
c
       return
       end
c     
c-----------------------------------------------------------------------
c      check divergence
c-----------------------------------------------------------------------
       subroutine checkdivergence(uu)
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
       aa=9.D0/8.D0
       bb=-1.D0/8.D0
c
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
c       write(*,*) 'Check divergence :', sum(divu)
       return
      end 
