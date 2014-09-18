c     
c
c-----------------------------------------------------------------------------
c
       module outflowBC
       save
        double precision,allocatable:: bc_outflowu(:,:,:,:)
        double precision,allocatable:: bc_outflow_w_s(:,:,:) 
       end module outflowBC
c------------------------------------------------
       subroutine allocate_outflow
       use mesh
       use mesh_uparms
       use outflowBC  
c
        allocate(bc_outflowu(nghost,iylo:iyhi,izlo:izhi,3))    !!??? For bc_outflow_w_s required?
        !
        bc_outflowu = 1.d3
        
        allocate(bc_outflow_w_s(nghost,iylo:iyhi,izlo:izhi))    !!??? For bc_outflow_w_s required?
        !
        bc_outflow_w_s = 1.d3
        
       return
       end 
ccccccc---------------------------------
       subroutine init_outflow(un)
       use mesh
       use mesh_common
       use properties
       use outflowBC
c
       implicit none
       double precision:: un(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
       integer:: i       
c
        do i=1, nghost, 1
          bc_outflowu(i,iylo:iyhi,izlo:izhi,1)
     &    =un(nxlocal+1,iylo:iyhi,izlo:izhi,1)
          bc_outflowu(i,iylo:iyhi,izlo:izhi,2)
     &      =un(nxlocal,iylo:iyhi,izlo:izhi,2)
          bc_outflowu(i,iylo:iyhi,izlo:izhi,3)
     &      =un(nxlocal,iylo:iyhi,izlo:izhi,3)
        enddo
c
       return
       end
ccccccc---------------------------------
       subroutine init_outflow_w_s(wn_s)
       use mesh
       use mesh_common
       use properties
       use outflowBC
c
       implicit none
       double precision:: wn_s(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
       integer:: i       
c
        do i=1, nghost, 1
          bc_outflow_w_s(i,iylo:iyhi,izlo:izhi)
     &    =wn_s(nxlocal+1,iylo:iyhi,izlo:izhi)
        enddo
c
       return
       end       
c-----------------------------------------------------------------------------
c      Time marching of the outflow boundary condition for fluid
c-----------------------------------------------------------------------------
       subroutine outflow(un,substep)
       use mesh
       use mesh_common
       use properties
       use time
       use RKCoefficients
       use outflowBC
#ifdef PARALLEL
       use mpistuff
#endif
c
       double precision:: un(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
       double precision:: bc_outflow(IYLO:IYHI,IZLO:IZHI,nvar)
       double precision:: alpha_t, beta_t
       double precision:: convec_v(IZLO:IZHI),convec_w(IZLO:IZHI)
       integer:: i,j,k, substep
     
c%%%%% 
      bc_outflow(:,:,:)=bc_outflowu(1,:,:,:)       
c
c     first step: get the outflow velocity at the outlet
       alpha_t = alpha(substep)
       beta_t = beta(substep)
       convec_v = 0.d0
       convec_w = 0.d0
       i = nxlocal+1
       do k= izlo, izhi, 1
         do j= iylo, iyhi, 1  
           convec_v(k) = convec_v(k) + un(i,j,k,1)           
         enddo
         convec_v(k) = convec_v(k)/(iyhi-iylo+1) *(alpha_t+beta_t)
       enddo
       do k=izlo+1, izhi, 1
         convec_w(k)= (convec_v(k)+convec_v(k-1))/2.d0
       enddo
        convec_w(izlo)= convec_v(izlo)
c
       i = nxlocal
c      u component
       do j= iylo, iyhi, 1
         do k= izlo, izhi, 1
          bc_outflow(j,k,1) = bc_outflow(j,k,1) -
     &     convec_v(k)*dt*(bc_outflow(j,k,1)-un(i+1,j,k,1))/dx
         enddo
       enddo
c      v component
       do j =iylo, iyhi, 1
         do k= izlo, izhi, 1
          bc_outflow(j,k,2) = bc_outflow(j,k,2) -
     &     convec_v(k)*dt*(bc_outflow(j,k,2)-un(i,j,k,2))/dx
         enddo
       enddo
c       w component
       do j= iylo, iyhi, 1
         do k= izlo, izhi, 1
          bc_outflow(j,k,3) = bc_outflow(j,k,3) - 
     &     convec_w(k)*dt*(bc_outflow(j,k,3)-un(i,j,k,3))/dx
         enddo
       enddo
c      check output
c       write(*,*)'MAX/MIN outflow :',
c     &            maxval(bc_outflow(:,:,1)),minval(bc_outflow(:,:,1))

c%%%
c      get the velocity for ghost cells
c
      if (iprocx .eq. xprocs) then
         i=nxlocal+1
         do k=IZLO,IZHI,1
            do j=IYLO,IYHI,1
             bc_outflowu(1,j,k,1)=bc_outflow(j,k,1)
             bc_outflowu(2,j,k,1)=2.D0*bc_outflowu(1,j,k,1)-un(i,j,k,1)
             bc_outflowu(3,j,k,1)=
     &         4.d0*(un(i,j,k,1)+bc_outflowu(2,j,k,1)) 
     &              -un(i-1,j,k,1)-6.d0*bc_outflowu(1,j,k,1)
            enddo
         enddo     
         do j=IYLO,IYHI,1
           do k=IZLO,IZHI,1
             bc_outflowu(1,j,k,2:3)= bc_outflow(j,k,2:3)
             bc_outflowu(2,j,k,2:3)= 3.d0*bc_outflow(j,k,2:3)
     &        -3.d0*un(i-1,j,k,2:3) + 1.d0*un(i-2,j,k,2:3)
             bc_outflowu(3,j,k,2:3)= 
     &          2.d0*bc_outflowu(2,j,k,2:3)-bc_outflowu(1,j,k,2:3)
             bc_outflowu(4,j,k,2:3) = 
     &          2.d0*bc_outflowu(3,j,k,2:3)-bc_outflowu(2,j,k,2:3)
            enddo
         enddo
      endif

        j=1
       do k=izlo, izhi,1
c         k=5
c         write(*,*)k, un(nxlocal,j,k,3), un(nxlocal-1,j,k,3),
c     &                un(nxlocal-3,j,k,3)
c         write(*,*)k,bc_outflow(j,k,3), bc_outflowu(1,j,k,3),
c     &   bc_outflowu(2,j,k,3),bc_outflowu(3,j,k,3),bc_outflowu(4,j,k,3)
      enddo
c
       return
       end
c
c-----------------------------------------------------------------------------
c      Time marching of the outflow boundary condition for solid
c-----------------------------------------------------------------------------
       subroutine outflowSolid(w_s, ux_s,substep)
       use mesh
       use mesh_common
       use properties
       use time
       use RKCoefficients
       use outflowBC
#ifdef PARALLEL
       use mpistuff
#endif
c
       double precision:: ux_s(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar
     &                           ,nspecis_s)
       double precision:: w_s(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI
     &                           ,nspecis_s)
       integer :: snumb_s, substep
     
       do snumb_s =1,nspecis_s,1
          !call outflow(w_s(:,:,:,snumb),substep)  ?????
          call outflow_w_s(w_s(:,:,:,snumb_s),substep)  
          call outflow(ux_s(:,:,:,:,snumb_s),substep)
       enddo
c
       return
       end subroutine outflowSolid
c-----------------------------------------------------------------------------

c-----------------------------------------------------------------------------
       subroutine outflow_w_s(wn_s,substep)
       use mesh
       use mesh_common
       use properties
       use time
       use RKCoefficients
       use outflowBC
#ifdef PARALLEL
       use mpistuff
#endif
c
       double precision:: wn_s(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
       double precision:: bc_outflow_s(IYLO:IYHI,IZLO:IZHI)
       double precision:: alpha_t, beta_t
       double precision:: convec_v(IZLO:IZHI),convec_w(IZLO:IZHI)
       integer:: i,j,k, substep
     
c%%%%% 
      bc_outflow_s(:,:)=bc_outflow_w_s(1,:,:)     !!!! ??? define bc_outflow_w_s(1,:,:) in module outflowBC, UP(TOP)
c
c     first step: get the outflow velocity at the outlet
       alpha_t = alpha(substep)
       beta_t = beta(substep)
       convec_v = 0.d0
       convec_w = 0.d0
       i = nxlocal+1
       do k= izlo, izhi, 1
         do j= iylo, iyhi, 1  
           convec_v(k) = convec_v(k) + wn_s(i,j,k)       !???  define subroutine init_outflow(wn_s) above similar to un   
         enddo
         convec_v(k) = convec_v(k)/(iyhi-iylo+1) *(alpha_t+beta_t)
       enddo
       do k=izlo+1, izhi, 1
         convec_w(k)= (convec_v(k)+convec_v(k-1))/2.d0
       enddo
        convec_w(izlo)= convec_v(izlo)
c
       i = nxlocal
c      u component
       do j= iylo, iyhi, 1
         do k= izlo, izhi, 1
!          bc_outflow(j,k) = bc_outflow_s(j,k) !-
!     &     convec_v(k)*dt*( bc_outflow_s(j,k) )/dx            !! un->wn_s    -wn_s(i+1,j,k)
         enddo
       enddo
!! Not required       
!c      v component
!       do j =iylo, iyhi, 1
!         do k= izlo, izhi, 1
!          bc_outflow(j,k,2) = bc_outflow(j,k,2) -
!     &     convec_v(k)*dt*(bc_outflow(j,k,2)-un(i,j,k,2))/dx
!         enddo
!       enddo
!c       w component
!       do j= iylo, iyhi, 1
!         do k= izlo, izhi, 1
!          bc_outflow(j,k,3) = bc_outflow(j,k,3) - 
!     &     convec_w(k)*dt*(bc_outflow(j,k,3)-un(i,j,k,3))/dx
!         enddo
!       enddo
c      check output
c       write(*,*)'MAX/MIN outflow :',
c     &            maxval(bc_outflow(:,:,1)),minval(bc_outflow(:,:,1))

c%%%
c      get the velocity for ghost cells
c
      if (iprocx .eq. xprocs) then
         i=nxlocal+1
         do k=IZLO,IZHI,1
            do j=IYLO,IYHI,1
             bc_outflow_w_s(1,j,k)=bc_outflow_s(j,k)
             bc_outflow_w_s(2,j,k)=2.D0*bc_outflow_w_s(1,j,k)
     &                              -wn_s(i,j,k)  !un->wn_s
             bc_outflow_w_s(3,j,k)=
     &         4.d0*(wn_s(i,j,k)+bc_outflow_w_s(2,j,k)) 
     &              -wn_s(i-1,j,k)-6.d0*bc_outflow_w_s(1,j,k)
            enddo
         enddo 
!!! Consider above interpolation if not the above interpolation then this ???            
!         do j=IYLO,IYHI,1
!           do k=IZLO,IZHI,1
!             bc_outflow_w_s(1,j,k)= bc_outflow_s(j,k)
!             bc_outflow_w_s(2,j,k)= 3.d0*bc_outflow_s(j,k)   !un->wn_s
!     &        -3.d0*wn_s(i-1,j,k) + 1.d0*wn_s(i-2,j,k)
!             bc_outflow_w_s(3,j,k)= 
!     &          2.d0*bc_outflow_w_s(2,j,k)-bc_outflow_w_s(1,j,k)
!             bc_outflow_w_s(4,j,k) = 
!     &          2.d0*bc_outflow_w_s(3,j,k)-bc_outflow_w_s(2,j,k)
!            enddo
!         enddo
      endif

        j=1
       do k=izlo, izhi,1
c         k=5
c         write(*,*)k, un(nxlocal,j,k,3), un(nxlocal-1,j,k,3),
c     &                un(nxlocal-3,j,k,3)
c         write(*,*)k,bc_outflow(j,k,3), bc_outflowu(1,j,k,3),
c     &   bc_outflowu(2,j,k,3),bc_outflowu(3,j,k,3),bc_outflowu(4,j,k,3)
      enddo
c
       return
       end
