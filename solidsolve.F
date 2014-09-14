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
       subroutine SolidSolve(sux,timestep)  !!! Remove the px, ! Replace ux -> wux !!SolidSolve(wux,timestep)
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
       double precision:: sux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)	!ux -> !! wux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
       !double precision:: px(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI)		!! Remove px
       integer:: timestep						! Can be done
c
c     RK intermediate stages solution
       double precision:: suxp(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)	!uxp -> wuxp(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
       double precision:: suxpp(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)	!uxpp -> wuxpp
       double precision:: sutmp(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)	!utmp -> wutmp
c
       !double precision:: wforce(nxlocal+1,nylocal+1,nzlocal+1,nvar) !Remove wind turbine force at various locations
c
#ifndef TWO_D
       double precision:: Adsu(nxlocal+1,nylocal+1,nzlocal+1,nvar)	!Adu -> Adwu
       double precision:: Adsup(nxlocal+1,nylocal+1,nzlocal+1,nvar)	!Adup -> Adwup
       double precision:: Adss(nxlocal+1,nylocal+1,nzlocal+1,nvar)	!Ads -> Ads
c
       !double precision:: Lu(nxlocal+1,nylocal+1,nzlocal+1,nvar)	!! Remove 
       !double precision:: Lup(nxlocal+1,nylocal+1,nzlocal+1,nvar)	!! Remove 
       !double precision:: gradp(nxlocal+1,nylocal+1,nzlocal+1,nvar)	!! Remove 
#else
       double precision:: Adsu(nxlocal+1,nylocal,nzlocal+1,nvar)		!Adu -> Adwu
       double precision:: Adsup(nxlocal+1,nylocal,nzlocal+1,nvar)	!Adup -> Adwup	
       double precision:: Adss(nxlocal+1,nylocal,nzlocal+1,nvar)		!Ads -> Ads
c
       !double precision:: Lu(nxlocal+1,nylocal,nzlocal+1,nvar)
       !double precision:: Lup(nxlocal+1,nylocal,nzlocal+1,nvar)
       !double precision:: gradp(nxlocal+1,nylocal,nzlocal+1,nvar)
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
       call SetBoundaryValuesST(sux)
       call exchange(sux,nvar)
       !gradp=0.D0
       !Lu=0.D0
       !Lup=0.D0
       Adsu=0.D0
       Adsup=0.D0
       Adss=0.D0
       suxp=sux
       suxpp=sux
c
c       write(6,*) 'Fluid solve dt=',dt, beta(1), Reynolds
c
c----------------- Stage 1 of RK-----------------------------------------
       call ComputeAdvectionTerms(Adu_s,ux_s)            ! Adu is for intermediate v* i.e. ux !!! call ComputeDragTerms(Adup_s,uxp_s)  inside ComputeAdvectionTerms(Adu_s,ux_s) but we need fluid velocity too...
           
#ifdef LES
c       call fftnonlinear(Adu)
       !sutmp = sux
       !call sgsnonlinear(Adsu,sutmp,timestep)	!Remove: No subgridscale model
#endif
       Ads_s=Adu_s*gamma(1)
#ifdef TURBINE
       !sutmp = sux
       !call windforce(sutmp,swforce,timestep,1)	!Remove wind turbine at various locations
       !Adss = Adss + swforce			! Remove this Wind force term by wind turbine.
#endif
       call outflow(sux,1)			! Will consider the outflow BCs for solid later
       call ComputeLinearTerms(Lu,gradp,ux,uxp,Ads,px,	! Edit it...
     &                         alpha(1),beta(1),dt,timestep,1) !Remove the compute of Linear part of Eq(3)
       call SetBoundaryValuesST(uxp_s) 		!Later Bounadary Conditions: ground wall + other BCs
       call exchange(uxp_s,nvar)				!Later for MPI: Parallel computation
c
c
c      check
c       write(6,*)'Min/max Adu',maxval(Adu(:,:,:,1)),minval(Adu(:,:,:,1))
c       write(6,*) 'Min/max Lu',maxval(Lu(:,:,:,1)),minval(Lu(:,:,:,1))
c       write(6,*) 'Min/max gradp',maxval(gradp(:,:,:,1)),
c     &      minval(gradp(:,:,:,1)), maxval(dabs(gradp(:,:,:,3)))
c
c------------------Stage 2 of RK----------------------------------------
       !gradp=0.D0
       call ComputeAdvectionTerms(Adup_s,uxp_s)    ! Adup is for intermediate v* i.e. uxp  !!! call ComputeDragTerms(Adup_s,uxp_s)  inside ComputeAdvectionTerms(Adu_s,ux_s) but we need fluid velocity too...
#ifdef LES
c       call fftnonlinear(Adup)
       !sutmp = suxp
       !call sgsnonlinear(Adup,utmp,timestep)	!Remove: No subgridscale model
#endif
       Ads_s=Adu_s*zeta(1)+Adup_s*gamma(2)
#ifdef TURBINE
       !utmp = uxp
       !call windforce(utmp,wforce,timestep,2)		!Remove wind turbine at various locations
       !Ads = Ads + wforce			! Remove this Wind force term by wind turbine.
#endif
       call outflow(uxp_s,2)			! Will consider the outflow BCs for solid later
       
       !call ComputeLinearTerms(Lup,gradp,uxp,uxpp,Ads,px, !Remove
     !&                         alpha(2),beta(2),dt,timestep,2)  !As computes the Linear part of Eq(3)
       call SetBoundaryValuesST(uxpp_s)		!Later Bounadary Conditions: ground wall + other BCs
       call exchange(uxpp_s,nvar)			!Later for MPI: Parallel computation
c       write(6,*)'Min/max Adu',
c     &      maxval(Adup(:,:,:,1)),minval(Adup(:,:,:,1))
c       write(6,*) 'Min/max Lu',maxval(Lup(:,:,:,1)),minval(Lup(:,:,:,1))
c       write(6,*) 'Min/max gradp',maxval(gradp(:,:,:,1)),
c     &      minval(gradp(:,:,:,1)), maxval(dabs(gradp(:,:,:,3)))
c
c-------------------Stage 3 of RK---------------------------------------
       !gradp = 0.d0
       call ComputeAdvectionTerms(Adu_s,uxpp_s)    ! Adu is for intermediate v* i.e. uxpp
#ifdef LES
c       call fftnonlinear(Adu_s)
       !utmp = uxpp
       !call sgsnonlinear(Adu,utmp,timestep) 	!Remove: No subgridscale model
#endif
       Ads_s=Adup_s*zeta(2)+Adu_s*gamma(3)  
#ifdef TURBINE
       !utmp = uxpp
       !call windforce(utmp,wforce,timestep,3) 	!Remove wind turbine at various locations
       !Ads = Ads + wforce		! Will consider the outflow BCs for solid later
#endif
       call outflow(uxpp_s,3)
       !call ComputeLinearTerms(Lu,gradp,uxpp,ux,Ads,px,  !Remove
     !&                         alpha(3),beta(3),dt,timestep,3)  !As computes the Linear part of Eq(3)
       call SetBoundaryValuesST(ux_s) 		!Later Bounadary Conditions: ground wall + other BCs
       call exchange(ux_s,nvar)		!Later for MPI: Parallel computation

c       write(6,*)'Min/max Adu',maxval(Adu(:,:,:,1)),minval(Adu(:,:,:,1))
c       write(6,*) 'Min/max Lu',maxval(Lu(:,:,:,1)),minval(Lu(:,:,:,1))
c       write(6,*) 'Min/max gradp',maxval(gradp(:,:,:,1)),
c     &      minval(gradp(:,:,:,1)), maxval(dabs(gradp(:,:,:,3)))
c
c
       call checkdivergence(ux_s)
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
