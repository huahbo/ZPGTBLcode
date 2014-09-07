c----------------------------------------------------------------------
c     Wan Cheng 
c     August 2012
c----------------------------------------------------------------------
c
c
c
c----------------------------------------------------------------------
c     init diagnostics 
c----------------------------------------------------------------------
      subroutine init_mean(ux)
      use mesh
      use mesh_common
#ifdef PARALLEL
      use mpistuff
#endif
c      use time
      use stats
      use properties
      use virtualwallBC
c
      implicit none
      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI, nvar)
c
      double precision::   utmp(IXLO:IXHI,IZLO:IZHI,nvar) 
      double precision::   uwall_t(1:nxlocal), guwall_t(1:nx)
      double precision:: g_del99t(1:nx),g_deltat(1:nx),g_thetat(1:nx)
      double precision:: g_utaut(1:nx)
      double precision:: g_u2(nz),g_u3(nz),g_u4(nz),g_Te(nz,6)
      double precision:: g_uv(nz),g_uw(nz),g_vv(nz),g_vw(nz),g_ww(nz)
      integer:: k, ilocal
c
      uwall_t  = 0.d0
      guwall_t = 0.d0

c
c     local mean velocity 
      call get_utmp(ux,utmp,uwall_t)
      st_umean = 0.d0
      st_umean(1:nxlocal+1,izlo:izhi,1:3)=
     &    utmp(1:nxlocal+1,izlo:izhi,1:3)
      call get_guwall(uwall_t,guwall_t)
      g_vwall = guwall_t
c      
c     utau
      call get_utau_mean(g_utaut)
      g_utau(1:nx) = g_utaut(1:nx)
      if(iprocx.lt.xprocs) then
      st_utau(1:nxlocal+1)=g_utau((iprocx-1)*nxlocal+1:iprocx*nxlocal+1)
      else
      st_utau(1:nxlocal)=g_utau((iprocx-1)*nxlocal+1:iprocx*nxlocal) 
      st_utau(nxlocal+1)=st_utau(nxlocal) ! the last point, give this value is ok   
      endif
c
c     fluctuation
      do ilocal=1, nxlocal, 1
        call get_ufluc(ilocal,ux,g_u2,g_u3,g_u4,
     &                 g_uv,g_uw,g_vv,g_vw,g_ww,g_Te)
        do k=1, nz
          g_fu2(ilocal,k)= g_u2(k) 
          g_fu3(ilocal,k)= g_u3(k) 
          g_fu4(ilocal,k)= g_u4(k) 
          g_fuv(ilocal,k)= g_uv(k) 
          g_fuw(ilocal,k)= g_uw(k) 
          g_fvv(ilocal,k)= g_vv(k) 
          g_fvw(ilocal,k)= g_vw(k) 
          g_fww(ilocal,k)= g_ww(k) 
          g_Ten(ilocal,k,1:6)=g_Te(k,1:6)
        enddo
      enddo
c
c     boundary layer thickness
      call get_blthickness(utmp,uwall_t,g_del99t,g_deltat,g_thetat)
      g_del99(1:nx)= g_del99t(1:nx) 
      g_delta(1:nx)= g_deltat(1:nx) 
      g_theta(1:nx)= g_thetat(1:nx) 
      if(iprocx.lt.xprocs) then
      st_del99(1:nxlocal+1)=
     &        g_del99((iprocx-1)*nxlocal+1:iprocx*nxlocal+1)
      st_delta(1:nxlocal+1)=
     &        g_delta((iprocx-1)*nxlocal+1:iprocx*nxlocal+1)
      st_theta(1:nxlocal+1)=
     &        g_theta((iprocx-1)*nxlocal+1:iprocx*nxlocal+1)
      else
      st_del99(1:nxlocal)=g_del99((iprocx-1)*nxlocal+1:iprocx*nxlocal)
      st_del99(nxlocal+1)=g_del99(nx)
      st_delta(1:nxlocal)=g_delta((iprocx-1)*nxlocal+1:iprocx*nxlocal)
      st_delta(nxlocal+1)=g_delta(nx)
      st_theta(1:nxlocal)=g_theta((iprocx-1)*nxlocal+1:iprocx*nxlocal)
      st_theta(nxlocal+1)=g_theta(nx)
      endif
c
c     
c  check output screen
c       if(iprocx.eq.1 .and. iprocy.eq.1 .and. iprocz.eq.1) then
c         do i=1, nx
c          write(*,*) i, g_del99(i), g_delta(i), g_theta(i)
c         enddo 
c       endif
c         if( iprocy.eq.1 .and. iprocz.eq.1) then
c          do i=1, nxlocal, 1
c            i= 96
c         write(*,*)'blthickness :',i,g_del99(i),g_delta(i),g_theta(i),
c     &                       g_utau(i) , iprocx, iprocy,iprocz
c           write(*,*)'blthickness :',i,g_theta(i),st_theta(i)
c     &                      , iprocx, iprocy,iprocz
c           i =1 
c           write(*,*)'blthickness :',i,g_utau(i),st_utau(i),
c     &   iprocx, iprocy,iproczc
c           enddo
c          endif
c
      return
      end
c
c----------------------------------------------------------------------
c     Diagnostics at each timestep
c----------------------------------------------------------------------    
      subroutine Diagnostics(ux,timestep,ndiag) 
      use mesh
      use mesh_common
      use time
      use properties
      use stats
      use InflowVelBC
      use virtualwallBC
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision:: ux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      integer:: timestep, ndiag
c
      double precision:: utmp(IXLO:IXHI,IZLO:IZHI,nvar)
      double precision:: uwall_t(1:nxlocal), guwall_t(1:nx)
      double precision:: g_del99t(1:nx),g_deltat(1:nx),g_thetat(1:nx)
      double precision:: g_utaut(1:nx)
c
      double precision:: g_u2(nz),g_u3(nz),g_u4(nz),g_Te(nz,6)
      double precision:: g_uv(nz),g_uw(nz),g_vv(nz),g_vw(nz),g_ww(nz)
c
      double precision:: Lundave
      integer:: i,k,l, ilocal
c
c     set average parameter

c     simple average for code B simulation
#ifdef CODEB
      Tave =  ttot -tstart
#endif

c     mean velocity
      call get_utmp(ux,utmp,uwall_t)
      do i=1, nxlocal+1, 1
        do k=izlo, izhi, 1
          st_umean(i,k,1)=Lundave(utmp(i,k,1),st_umean(i,k,1),dt,Tave)
          st_umean(i,k,2)=Lundave(utmp(i,k,2),st_umean(i,k,2),dt,Tave)
          st_umean(i,k,3)=Lundave(utmp(i,k,3),st_umean(i,k,3),dt,Tave)
        enddo
      enddo
      call get_guwall(uwall_t,guwall_t)
      do i=1, nx
         g_vwall(i) = Lundave(guwall_t(i),g_vwall(i),dt,Tave)
      enddo
      
c     fluctuation
      do ilocal=1, nxlocal, 1
        call get_ufluc(ilocal,ux,g_u2,g_u3,g_u4,
     &                 g_uv,g_uw,g_vv,g_vw,g_ww,g_Te)
        do k=1, nz
          g_fu2(ilocal,k)=Lundave(g_u2(k),g_fu2(ilocal,k),dt,Tave)
          g_fu3(ilocal,k)=Lundave(g_u3(k),g_fu3(ilocal,k),dt,Tave)
          g_fu4(ilocal,k)=Lundave(g_u4(k),g_fu4(ilocal,k),dt,Tave)
          g_fuv(ilocal,k)=Lundave(g_uv(k),g_fuv(ilocal,k),dt,Tave)
          g_fuw(ilocal,k)=Lundave(g_uw(k),g_fuw(ilocal,k),dt,Tave)
          g_fvv(ilocal,k)=Lundave(g_vv(k),g_fvv(ilocal,k),dt,Tave)
          g_fvw(ilocal,k)=Lundave(g_vw(k),g_fvw(ilocal,k),dt,Tave)
          g_fww(ilocal,k)=Lundave(g_ww(k),g_fww(ilocal,k),dt,Tave)
          do l=1,6
         g_Ten(ilocal,k,l)=Lundave(g_Te(k,l),g_Ten(ilocal,k,l),dt,Tave)
          enddo
        enddo
      enddo
c
c     boundary layer thickness 
      call get_blthickness(utmp,uwall_t,g_del99t,g_deltat,g_thetat)
      do i=1, nx
        g_del99(i)=Lundave(g_del99t(i),g_del99(i),dt,Tave)
        g_delta(i)=Lundave(g_deltat(i),g_delta(i),dt,Tave)
        g_theta(i)=Lundave(g_thetat(i),g_theta(i),dt,Tave)
      enddo
      if(iprocx.lt.xprocs) then
      st_del99(1:nxlocal+1)=
     &        g_del99((iprocx-1)*nxlocal+1:iprocx*nxlocal+1)
      st_delta(1:nxlocal+1)=
     &        g_delta((iprocx-1)*nxlocal+1:iprocx*nxlocal+1)
      st_theta(1:nxlocal+1)=
     &        g_theta((iprocx-1)*nxlocal+1:iprocx*nxlocal+1)
      else
      st_del99(1:nxlocal)=g_del99((iprocx-1)*nxlocal+1:iprocx*nxlocal)
      st_del99(nxlocal+1)=g_del99(nx)
      st_delta(1:nxlocal)=g_delta((iprocx-1)*nxlocal+1:iprocx*nxlocal)
      st_delta(nxlocal+1)=g_delta(nx)
      st_theta(1:nxlocal)=g_theta((iprocx-1)*nxlocal+1:iprocx*nxlocal)
      st_theta(nxlocal+1)=g_theta(nx)
      endif
c
c     utau  
c
c     utau
      call get_utau_mean(g_utaut)
      do i=1, nx, 1
        g_utau(i) = Lundave(g_utaut(i),g_utau(i),dt,Tave)
      enddo
      if(iprocx.lt.xprocs) then
      st_utau(1:nxlocal+1)=g_utau((iprocx-1)*nxlocal+1:iprocx*nxlocal+1)
      else
      st_utau(1:nxlocal)=g_utau((iprocx-1)*nxlocal+1:iprocx*nxlocal) 
      st_utau(nxlocal+1)=st_utau(nxlocal) ! the last point, give this value is ok   
      endif
c

       call get_wallpara

c      output part
       if(mod(timestep,ndiag).eq.1) then
c         call SpectralU(ux,timestep)
c         call SpectralU1Dy(ux,timestep)
c         call SpectralU1Dz(ux,timestep)
       endif
c          write(*,*) 'test'
       call output_xslice(ux,timestep,ndiag)
c
       return
       end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
      subroutine get_utmp(ux,utmp,uwall_t)
      use mesh
      use mesh_common
      use time
      use stats
      use properties
      use virtualwallBC
#ifdef PARALLEL
      use mpistuff
#endif
c
      implicit none
      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI, nvar)
      double precision::   utmp(IXLO:IXHI,IZLO:IZHI,nvar)        
      double precision::   uwall_t(1:nxlocal)
c
      double precision:: mpi_ut(IXLO:IXHI,IZLO:IZHI,nvar)
      double precision:: mpi_uwall(1:nxlocal)
      integer:: i, j, utsize
c
c     velocity for statistics
      utmp = 0.d0
      do j = 1, nylocal
        utmp(ixlo:ixhi,  izlo:izhi,1:3) = utmp(ixlo:ixhi,izlo:izhi,1:3)
     &   +ux(ixlo:ixhi,j,izlo:izhi,1:3)/nylocal    
      enddo 
      uwall_t = 0.d0
      if(iprocz.eq.1) then
        do i=1, nxlocal, 1
          do j=1, nylocal, 1
            uwall_t(i) = uwall_t(i)+ bc_uwall(i,j,1)/nylocal
          enddo
        enddo
      endif
#ifdef PARALLEL
c      if(yprocs.gt.1) then          
        utsize = (ixhi-ixlo+1)*(izhi-izlo+1)*nvar
        call MPI_Allreduce(utmp, mpi_ut, utsize,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm1d, ierr)
        call ErrorHandler(ierr,ERROR_ALLREDUCE)
c
        call MPI_Allreduce(uwall_t, mpi_uwall, nxlocal,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm1d, ierr)
        call ErrorHandler(ierr,ERROR_ALLREDUCE)
        utmp = mpi_ut/yprocs
        uwall_t = mpi_uwall /yprocs
c      endif
#endif
c
       return
       end
c----------------------------------------------------------------------
c      global uwall
c----------------------------------------------------------------------
       subroutine get_guwall(uwall_t,guwall_t)
      use mesh
      use mesh_common
      use virtualwallBC
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision:: uwall_t(1:nxlocal), guwall_t(1:nx), g_t(1:nx)
c      integer:: i
c
      g_t = 0.d0
c
      if(iprocy.eq.1 .and. iprocz.eq.1) 
     &   g_t((iprocx-1)*nxlocal+1:iprocx*nxlocal) = uwall_t(1:nxlocal)
#ifdef PARALLEL
      call MPI_Allreduce(g_t, guwall_t, nx,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
       call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
      guwall_t = g_t
#endif      
c
c      if(iprocx.eq.1 .and .iprocy.eq.1 .and. iprocz.eq.1) then
c        do i=1, nxlocal
c            write(*,*) i, guwall_t(i), uwall_t(i)
c         enddo
c      endif

       return
       end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
      subroutine get_wallpara
      use mesh
      use mesh_common
      use time
      use virtualwallBC
      use stats
      use GClaw
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision:: g_t(1:nx), g_n(1:nx), v_t(1:nxlocal+1)
      double precision:: Lundave
      integer:: i, j
c
      g_t = 0.d0
      v_t = 0.d0
      do j = 1, nylocal
        v_t(1:nxlocal) = v_t(1:nxlocal)+bc_kappa(1:nxlocal,j)/nylocal
      enddo 
      if(iprocy.eq.1 .and. iprocz.eq.1)
     &  g_t((iprocx-1)*nxlocal+1:iprocx*nxlocal)=v_t(1:nxlocal)
#ifdef PARALLEL
      g_n = 0.d0
      call MPI_Allreduce(g_t, g_n, nx,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
       call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
      g_n = g_t
#endif   
c      g_kappa(1:nx) = g_n(1:nx)
      do i=1, nx
        g_kappa(i)=Lundave(g_n(i),g_kappa(i),dt,Tave)
      enddo
c
      g_t = 0.d0
      v_t = 0.d0
      do j = 1, nylocal
        v_t(1:nxlocal) = v_t(1:nxlocal)+GC_co(1:nxlocal,j)/nylocal    
      enddo 
      if(iprocy.eq.1 .and. iprocz.eq.1)
     &  g_t((iprocx-1)*nxlocal+1:iprocx*nxlocal)=v_t(1:nxlocal)
#ifdef PARALLEL
      call MPI_Allreduce(g_t, g_n, nx,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
       call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
      g_n = g_t
#endif   
c      g_GCco(1:nx) = g_n(1:nx)
      do i=1, nx
        g_GCco(i)=Lundave(g_n(i),g_GCco(i),dt,Tave)
      enddo
c
      g_t = 0.d0
      v_t = 0.d0
      do j = 1, nylocal
        v_t(1:nxlocal) = v_t(1:nxlocal)+GC_ci(1:nxlocal,j)/nylocal    
      enddo 
      if(iprocy.eq.1 .and. iprocz.eq.1)
     &  g_t((iprocx-1)*nxlocal+1:iprocx*nxlocal)=v_t(1:nxlocal)
#ifdef PARALLEL
      call MPI_Allreduce(g_t, g_n, nx,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
       call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
      g_n = g_t
#endif   
c      g_GCci(1:nx) = g_n(1:nx)
      do i=1, nx
        g_GCci(i)=Lundave(g_n(i),g_GCci(i),dt,Tave)
      enddo      
c
      g_t = 0.d0
      v_t = 0.d0
      do j = 1, nylocal
        v_t(1:nxlocal) = v_t(1:nxlocal)+GC_gam(1:nxlocal,j)/nylocal    
      enddo 
      if(iprocy.eq.1 .and. iprocz.eq.1)
     &  g_t((iprocx-1)*nxlocal+1:iprocx*nxlocal)=v_t(1:nxlocal)
#ifdef PARALLEL
      call MPI_Allreduce(g_t, g_n, nx,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
       call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
      g_n = g_t
#endif   
c      g_GCgam(1:nx) = g_n(1:nx)
      do i=1, nx
        g_GCgam(i)=Lundave(g_n(i),g_GCgam(i),dt,Tave)
      enddo
c
      g_t = 0.d0
      v_t = 0.d0
      do j = 1, nylocal
        v_t(1:nxlocal) = v_t(1:nxlocal)+GC_pre(1:nxlocal,j)/nylocal    
      enddo 
      if(iprocy.eq.1 .and. iprocz.eq.1)
     &  g_t((iprocx-1)*nxlocal+1:iprocx*nxlocal)=v_t(1:nxlocal)
#ifdef PARALLEL
      call MPI_Allreduce(g_t, g_n, nx,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
       call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
      g_n = g_t
#endif   
c      g_GCpre(1:nx) = g_n(1:nx)
      do i=1, nx
        g_GCpre(i)=Lundave(g_n(i),g_GCpre(i),dt,Tave)
      enddo
c

   
c
c      if(iprocx.eq.1 .and .iprocy.eq.1 .and. iprocz.eq.1) then
c        do i=1, nxlocal
c            write(*,*) i, guwall_t(i), uwall_t(i)
c         enddo
c      endif

       return
       end       
c
c----------------------------------------------------------------------
c      fluctuation
c----------------------------------------------------------------------
       subroutine get_ufluc(ilocal,ux,g_u2,g_u3,g_u4,
     &                      g_uv,g_uw,g_vv,g_vw,g_ww,g_T)
       use mesh
       use mesh_common
       use time
       use stats
       use properties
       use virtualwallBC
#ifdef PARALLEL
      use mpistuff
#endif
       implicit none
c=======Declarations=========
      double precision:: ux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: g_u2(nz),g_u3(nz),g_u4(nz),g_T(nz,6)
      double precision:: g_uv(nz),g_uw(nz),g_vv(nz),g_vw(nz),g_ww(nz)
c
      double precision:: ufl_u1(nylocal,nzlocal,3)
      double precision:: ufl_u2(nylocal,nzlocal)
      double precision:: ufl_u3(nylocal,nzlocal)
      double precision:: ufl_u4(nylocal,nzlocal)
      double precision:: ufl_uv(nylocal,nzlocal)
      double precision:: ufl_uw(nylocal,nzlocal)
      double precision:: ufl_vv(nylocal,nzlocal)
      double precision:: ufl_vw(nylocal,nzlocal)
      double precision:: ufl_ww(nylocal,nzlocal)
      double precision:: Tijtmp(nylocal,nzlocal,6)
c
c      double precision:: gta(nz), gt(nz)
      double precision:: gt_Tt(nz,6)
c

      double precision:: g_uft(nz), gt_uft(nz)
      double precision:: utmp
      double precision:: aa,bb
c
      integer:: i,j,k,l, iplocal, ilocal
c
      aa = 9.d0/8.d0
      bb =-1.d0/8.d0
c     
c      iplocal = int(ioutput/nxlocal)+1
c      ilocal  = mod(ioutput,nxlocal)
c
      ufl_u1 = 0.d0
      ufl_u2 = 0.d0
      ufl_u3 = 0.d0
      ufl_u4 = 0.d0
      ufl_uv = 0.d0
      ufl_uw = 0.d0
      ufl_vv = 0.d0
      ufl_vw = 0.d0
      ufl_ww = 0.d0
c
c      if( iprocx.eq.iplocal ) then
        i = ilocal
        do j=1, nylocal, 1
          do k=1, nzlocal, 1
            ! velocity fluctuation at center
            ufl_u1(j,k,1)= (
     &            ux(i  ,j,k,1)-st_umean(i  ,k,1) +
     &            ux(i+1,j,k,1)-st_umean(i+1,k,1)  )/2.d0
            ufl_u1(j,k,2)= 
     &            aa*( ux(i,j  ,k,2)+ux(i,j+1,k,2) )/2.d0+
     &            bb*( ux(i,j-1,k,2)+ux(i,j+2,k,2) )/2.d0
            ufl_u1(j,k,3)= (
     &            ux(i,j,k  ,3)-st_umean(i,k  ,3) +
     &            ux(i,j,k+1,3)-st_umean(i,k+1,3)  )/2.d0
            ! normalize
            ufl_u1(j,k,1:3)=ufl_u1(j,k,1:3)
     &             / ( st_utau(i)/2.d0+st_utau(i+1)/2.d0 )
            Tijtmp(j,k,1:6)=st_Tijct(i,j,k,1:6)
     &             / ( st_utau(i)/2.d0+st_utau(i+1)/2.d0 )
     &             / ( st_utau(i)/2.d0+st_utau(i+1)/2.d0 )
            ! void small values
            if(dabs(ufl_u1(j,k,1)).lt.1.d-10) ufl_u1(j,k,1)=1.d-10
            if(dabs(ufl_u1(j,k,2)).lt.1.d-10) ufl_u1(j,k,2)=1.d-10
            if(dabs(ufl_u1(j,k,3)).lt.1.d-10) ufl_u1(j,k,3)=1.d-10
            ! Reynolds stress
            ufl_u2(j,k)=ufl_u1(j,k,1)*ufl_u1(j,k,1)
            ufl_uv(j,k)=ufl_u1(j,k,1)*ufl_u1(j,k,2)
            ufl_uw(j,k)=ufl_u1(j,k,1)*ufl_u1(j,k,3)
            ufl_vv(j,k)=ufl_u1(j,k,2)*ufl_u1(j,k,2)
            ufl_vw(j,k)=ufl_u1(j,k,2)*ufl_u1(j,k,3)
            ufl_ww(j,k)=ufl_u1(j,k,3)*ufl_u1(j,k,3)
            ! skew
            ufl_u3(j,k)=ufl_u1(j,k,1)**3.d0
            ! flat
            ufl_u4(j,k)=ufl_u1(j,k,1)**4.d0
          enddo
        enddo
c
c        
        g_uft = 0.d0
        gt_uft  = 0.d0 
#ifdef PARALLEL
c-----------------------------
c  uu
        do k= 1, nzlocal, 1
          utmp = 0.d0
          do j = 1, nylocal, 1
            utmp = utmp + ufl_u2(j,k)/nylocal
          enddo
          gt_uft((iprocz-1)*nzlocal+k) = utmp
        enddo
        call MPI_Allreduce(gt_uft, g_uft, nz,
     &         MPI_DOUBLE_PRECISION, MPI_SUM, comm2dyz, ierr)
        call ErrorHandler(ierr,ERROR_ALLREDUCE)
        g_u2 = g_uft / yprocs
c  uv
        do k= 1, nzlocal, 1
          utmp = 0.d0
          do j = 1, nylocal, 1
            utmp = utmp + ufl_uv(j,k)/nylocal
          enddo
          gt_uft((iprocz-1)*nzlocal+k) = utmp
        enddo
        call MPI_Allreduce(gt_uft, g_uft, nz,
     &         MPI_DOUBLE_PRECISION, MPI_SUM, comm2dyz, ierr)
        call ErrorHandler(ierr,ERROR_ALLREDUCE)
        g_uv = g_uft / yprocs
c  uw
        do k= 1, nzlocal, 1
          utmp = 0.d0
          do j = 1, nylocal, 1
            utmp = utmp + ufl_uw(j,k)/nylocal
          enddo
          gt_uft((iprocz-1)*nzlocal+k) = utmp
        enddo
        call MPI_Allreduce(gt_uft, g_uft, nz,
     &         MPI_DOUBLE_PRECISION, MPI_SUM, comm2dyz, ierr)
        call ErrorHandler(ierr,ERROR_ALLREDUCE)
        g_uw = g_uft / yprocs
c  vv
        do k= 1, nzlocal, 1
          utmp = 0.d0
          do j = 1, nylocal, 1
            utmp = utmp + ufl_vv(j,k)/nylocal
          enddo
          gt_uft((iprocz-1)*nzlocal+k) = utmp
        enddo
        call MPI_Allreduce(gt_uft, g_uft, nz,
     &         MPI_DOUBLE_PRECISION, MPI_SUM, comm2dyz, ierr)
        call ErrorHandler(ierr,ERROR_ALLREDUCE)
        g_vv = g_uft / yprocs
c  vw
        do k= 1, nzlocal, 1
          utmp = 0.d0
          do j = 1, nylocal, 1
            utmp = utmp + ufl_vw(j,k)/nylocal
          enddo
          gt_uft((iprocz-1)*nzlocal+k) = utmp
        enddo
        call MPI_Allreduce(gt_uft, g_uft, nz,
     &         MPI_DOUBLE_PRECISION, MPI_SUM, comm2dyz, ierr)
        call ErrorHandler(ierr,ERROR_ALLREDUCE)
           g_vw = g_uft / yprocs
c  ww
        do k= 1, nzlocal, 1
          utmp = 0.d0
          do j = 1, nylocal, 1
            utmp = utmp + ufl_ww(j,k)/nylocal
          enddo
          gt_uft((iprocz-1)*nzlocal+k) = utmp
        enddo
        call MPI_Allreduce(gt_uft, g_uft, nz,
     &         MPI_DOUBLE_PRECISION, MPI_SUM, comm2dyz, ierr)
        call ErrorHandler(ierr,ERROR_ALLREDUCE)
        g_ww = g_uft / yprocs
c  f3
        do k= 1, nzlocal, 1
          utmp = 0.d0
          do j = 1, nylocal, 1
            utmp = utmp + ufl_u3(j,k)/nylocal
          enddo
          gt_uft((iprocz-1)*nzlocal+k) = utmp
        enddo
        call MPI_Allreduce(gt_uft, g_uft, nz,
     &         MPI_DOUBLE_PRECISION, MPI_SUM, comm2dyz, ierr)
        call ErrorHandler(ierr,ERROR_ALLREDUCE)
        g_u3 = g_uft / yprocs
c  f4
        do k= 1, nzlocal, 1
          utmp = 0.d0
          do j = 1, nylocal, 1
            utmp = utmp + ufl_u4(j,k)/nylocal
          enddo
          gt_uft((iprocz-1)*nzlocal+k) = utmp
        enddo
        call MPI_Allreduce(gt_uft, g_uft, nz,
     &         MPI_DOUBLE_PRECISION, MPI_SUM, comm2dyz, ierr)
        call ErrorHandler(ierr,ERROR_ALLREDUCE)
        g_u4 = g_uft / yprocs
c  T
        gt_Tt = 0.d0
        do k= 1, nzlocal, 1
          do l=1, 6, 1
            utmp = 0.d0
            do j = 1, nylocal, 1
              utmp = utmp + Tijtmp(j,k,l)/nylocal
            enddo
            gt_Tt((iprocz-1)*nzlocal+k,l) = utmp 
          enddo
        enddo
        call MPI_Allreduce(gt_Tt, g_T, nz*6,
     &        MPI_DOUBLE_PRECISION, MPI_SUM, comm2dyz, ierr)
        call ErrorHandler(ierr,ERROR_ALLREDUCE)
        g_T = g_T / yprocs
#else
        g_u2 = 0.d0
        do k= 1, nzlocal, 1
          do j = 1, nylocal, 1
            g_u2(k) = g_u2(k) + ufl_u2(j,k)/nylocal
          enddo
        enddo
        g_u3 = 0.d0
        do k= 1, nzlocal, 1
          do j = 1, nylocal, 1
            g_u3(k) = g_u3(k) + ufl_u3(j,k)/nylocal
          enddo
        enddo
        g_u4 = 0.d0
        do k= 1, nzlocal, 1
          do j = 1, nylocal, 1
            g_u4(k) = g_u4(k) + ufl_u4(j,k)/nylocal
          enddo
        enddo
        g_uv = 0.d0
        do k= 1, nzlocal, 1
          do j = 1, nylocal, 1
            g_uv(k) = g_uv(k) + ufl_uv(j,k)/nylocal
          enddo
        enddo
        g_uw = 0.d0
        do k= 1, nzlocal, 1
          do j = 1, nylocal, 1
            g_uw(k) = g_uw(k) + ufl_uw(j,k)/nylocal
          enddo
        enddo
        g_vv = 0.d0
        do k= 1, nzlocal, 1
          do j = 1, nylocal, 1
            g_vv(k) = g_vv(k) + ufl_vv(j,k)/nylocal
          enddo
        enddo
        g_vw = 0.d0
        do k= 1, nzlocal, 1
          do j = 1, nylocal, 1
            g_vw(k) = g_vw(k) + ufl_vw(j,k)/nylocal
          enddo
        enddo      
        g_ww = 0.d0
        do k= 1, nzlocal, 1
          do j = 1, nylocal, 1
            g_ww(k) = g_ww(k) + ufl_ww(j,k)/nylocal
          enddo
        enddo
        g_T = 0.d0
        do k=1, nzlocal, 1
          do j=1, nylocal, 1
            g_T(k,1:6)=g_T(k,1:6)+Tijtmp(j,k,1:6)/nylocal 
          enddo
        enddo
#endif

c        do k=1, nz, 1
c          g_u2(k) = g_u2(k) + g_T(k,1)
c          g_uv(k) = g_uv(k) + g_T(k,2)
c          g_uw(k) = g_uw(k) + g_T(k,3) 
c          g_vv(k) = g_vv(k) + g_T(k,4)
c          g_vw(k) = g_vw(k) + g_T(k,5)
c          g_ww(k) = g_ww(k) + g_T(k,6) 
c        enddo
c      endif  ! iprocx=iplocal
c
       return
       end
c
c......................................--------------------------------
c      calculate boundary layer thickness
c......................................--------------------------------
#ifdef LES
       subroutine get_blthickness(utmp,uwall_t,
     &                g_del99t,g_deltat,g_thetat)
       use mesh
       use mesh_common
#ifdef PARALLEL
        use mpistuff
#endif
       use time
       use properties
       use stats
#ifdef LES
       use virtualwallBC
#endif
       implicit none
c=======Declarations=========
       double precision:: utmp(IXLO:IXHI,IZLO:IZHI,nvar)
       double precision:: delta_t(1:nxlocal)
       double precision:: del99_t(1:nxlocal)
       double precision:: theta_t(1:nxlocal)
       double precision:: g_del99t(1:nx),g_deltat(1:nx),g_thetat(1:nx)
       double precision:: g_t1d(1:nx)
       double precision:: uwall_t(1:nxlocal)
       double precision:: z_modify, ufree, vf, ut
       integer:: i, k, l  !,timestep
c       character*15 filename
c
c      local calculation of del99, delta, theta
       z_modify = bc_h0dz*dz
c    
       ufree = velRef 
       vf    = velRef
c
       delta_t = 0.d0
       del99_t = 0.d0
       theta_t = 0.d0
       do i=1, nxlocal, 1
c  del99
         l=0
         do k=1, nzlocal, 1
           if(utmp(i,k,1)/vf .lt. 0.99d0 ) then
             l= l+1
           endif
         enddo
         if(iprocz.eq.1) then
           del99_t(i) = (l-0.5d0 )*dz + z_modify
         else
           del99_t(i) = l * dz
         endif
         if( l.gt.0 .and. l.lt.nzlocal) del99_t(i) = del99_t(i) + 
     &       (0.99d0*vf-utmp(i,l,1))/(utmp(i,l+1,1)-utmp(i,l,1))*dz
c  delta
         if(iprocz.eq.1) then
          ut = uwall_t(i)
          delta_t(i) = delta_t(i) + (1.d0-ut/vf +1.d0)*z_modify/2.d0
         endif
         do k=1, nzlocal-1, 1
          ut = utmp(i,k,1) 
          delta_t(i) = delta_t(i) + (1.d0-ut/vf)*dz
         enddo
         ut = utmp(i,nzlocal,1)
         if(iprocz.eq.zprocs) then
          delta_t(i) = delta_t(i) + (1.d0-ut/vf)*dz*3.d0/4.d0
         else
          delta_t(i) = delta_t(i) + (1.d0-ut/vf)*dz
         endif
c   theta
         if(iprocz.eq.1) then
          ut = uwall_t(i)
          theta_t(i)=theta_t(i)+ut/vf*(1.d0-ut/vf+0.d0)*z_modify/2.d0
         endif
         do k=1, nzlocal-1, 1
          ut = utmp(i,k,1)
          theta_t(i) = theta_t(i) + ut/vf*(1.d0-ut/vf)*dz
         enddo
         ut = utmp(i,nzlocal,1)
         if(iprocz.eq.zprocs) then
          theta_t(i) = theta_t(i) + ut/vf*(1.d0-ut/vf)*dz*3.d0/4.d0
         else
          theta_t(i) = theta_t(i) + ut/vf*(1.d0-ut/vf)*dz 
         endif
c         write(*,*) i, del99_t(i), delta_t(i), theta_t(i)
       enddo

c      sum of local to get global
#ifdef PARALLEL
       g_t1d = 0.d0
       g_t1d((iprocx-1)*nxlocal+1:iprocx*nxlocal)=del99_t(1:nxlocal)
         call MPI_Allreduce(g_t1d, g_del99t, nx,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
         Call ErrorHandler(ierr,ERROR_ALLREDUCE)       
       g_del99t = g_del99t/yprocs
c  
       g_t1d = 0.d0
       g_t1d((iprocx-1)*nxlocal+1:iprocx*nxlocal)=delta_t(1:nxlocal)
         call MPI_Allreduce(g_t1d, g_deltat, nx,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
         Call ErrorHandler(ierr,ERROR_ALLREDUCE)    
       g_deltat = g_deltat/yprocs       
c
       g_t1d = 0.d0
       g_t1d((iprocx-1)*nxlocal+1:iprocx*nxlocal)=theta_t(1:nxlocal)
         call MPI_Allreduce(g_t1d, g_thetat, nx,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
         Call ErrorHandler(ierr,ERROR_ALLREDUCE)    
       g_thetat = g_thetat/yprocs
#else
       g_del99t = del99_t
       g_deltat = delta_t
       g_thetat = theta_t
#endif 
c
c       if(iprocx.eq.1 .and. iprocy.eq.1 .and. iprocz.eq.1) then
c         do i=1, nx
c          write(*,*) i, g_del99(i), g_delta(i), g_theta(i)
c         enddo 
c       endif
c
       return
       end
c
#endif
c----------------------------------------------------------------------
c      mean viscous velocity
c----------------------------------------------------------------------
#ifdef LES
       subroutine get_utau_mean(g_utaut)
       use mesh
       use mesh_common
#ifdef PARALLEL
        use mpistuff
#endif
       use time
       use properties
       use stats
#ifdef LES
       use virtualwallBC
#endif
       implicit none
c=======Declarations=========
       double precision:: utau_t(1:nxlocal)
       double precision:: g_utaut(1:nx), g_t(1:nx)
       integer::  jstat
c
       g_t = 0.d0
       utau_t = 0.d0
       if(iprocz.eq.1) then
       do jstat = 1, nylocal
       utau_t(1:nxlocal)= utau_t(1:nxlocal)+
     &      bc_utau(1:nxlocal,jstat)/nylocal 
       enddo
       endif
#ifdef PARALLEL
       g_t((iprocx-1)*nxlocal+1:iprocx*nxlocal)=utau_t(1:nxlocal)
         call MPI_Allreduce(g_t, g_utaut, nx,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
         Call ErrorHandler(ierr,ERROR_ALLREDUCE)    
       g_utaut = g_utaut/yprocs
#else
       g_utaut = utau_t
#endif       
c
c
       return
       end
c
#endif
c
c----------------------------------------------------------------------
c     1d average used by Lund
c----------------------------------------------------------------------
      double precision function Lundave(vnew, vold, dt, Tave)
      implicit none
      double precision vnew, vold, dt, Tave
c
      vnew = dt/Tave*vnew + (1.d0-dt/Tave)*vold  
      Lundave = vnew    
c
      return
      end






c----------------------------------------------------------------------
c   Part 3   
c        Some useful measurements
c----------------------------------------------------------------------
c----------------------------------------------------------
c      eigenvalues of Tijcn
c----------------------------------------------------------
       subroutine eigval_Tijcn(Tijcn,timestep)
       use mesh
       use mesh_common
       use properties    
       use diags_io        
       double precision:: Tijcn(0:nxlocal+2,0:nylocal+2,0:nzlocal+2,6)
       double precision:: Txx, Tyy, Tzz, Txy, Tyz, Tzx
       double precision:: eigval(3), maxeigval(3)
       integer::  coor(3,3)
c
       integer:: i,j,k, l,timestep
    
       maxeigval = 0.d0   
       do i= 0, nxlocal+2, 1
         do j= 0, nylocal+2, 1
           do k= 0, nzlocal+2, 1
             Txx = Tijcn(i,j,k,1)
             Txy = Tijcn(i,j,k,2)
             Tzx = Tijcn(i,j,k,3)
             Tyy = Tijcn(i,j,k,4)
             Tyz = Tijcn(i,j,k,5)
             Tzz = Tijcn(i,j,k,6)
             call Tij_eigenvalue_symm(Txx,Tyy,Tzz,Txy,Tyz,Tzx,eigval)
c
             l=3
               if(eigval(l).gt.maxeigval(l)) then
                maxeigval(l)=eigval(l)
                coor(l,1) = i
                coor(l,2) = j
                coor(l,3) = k
               endif
c
             l=1
               if(eigval(l).lt.maxeigval(l)) then
                maxeigval(l)=eigval(l)
                coor(l,1) = i
                coor(l,2) = j
                coor(l,3) = k
               endif

           enddo
         enddo
       enddo
       
       write(itijmax,*) timestep,maxeigval(3),
     &      coor(3,1),coor(3,2),coor(3,3)
       write(itijmin,*) timestep,maxeigval(1),
     &       coor(1,1),coor(1,2),coor(1,3)
c       write(*,*) coor(1,1), coor(1,2), coor(1,3)
c 
c           stop
       return
       end
c----------------------------------------------------------
c      Eigenvalues of Tij
c----------------------------------------------------------
       subroutine Tij_eigenvalue_symm
     &   (Sxx,Syy,Szz,Sxy,Syz,Szx,eigval)
c
       implicit none
       double precision:: eigval(0:2)
       double precision:: Sxx, Syy, Szz, Sxy, Syz, Szx
       double precision:: a, b, c, q, r, pai
       double precision:: costheta, theta
       double precision:: tmp, eps
       eps = 1.d-16
       pai = 4.d0*datan(1.d0)
c        x^3 = a* x^2 + b* x + c =0, where x is the eigenvalue
       a = -(Sxx+Syy+Szz)
       b = Sxx*Syy - Sxy*Sxy + Syy*Szz - Syz*Syz + Szz*Sxx - Szx*Szx
       c = -(   Sxx*(Syy*Szz-Syz*Syz)
     &         +Sxy*(Syz*Szx-Sxy*Szz)
     &         +Szx*(Sxy*Syz-Syy*Szx) )
       q = ( 3.d0*b - a*a )/ 9.d0
       r = ( 9.d0*a*b-27.d0*c-2.d0*a*a*a)/54.d0
c
       if(q .gt. 0.d0) then
          write(*,*) 'spiral eigenvalue symm q>= 0'
          stop
       endif 
c
       costheta = r/dsqrt(eps-q*q*q)
c      |costheta|>1 should not occur, except from round-off errors
       if( costheta .gt. 1.d0 ) then
             theta = 0.d0
       else if( costheta .lt. -1.d0 ) then
             theta = pai
       else
             theta = dacos(costheta) 
       endif
c
       eigval(0)=2.d0*dsqrt(eps-q)*dcos((theta         )/3.d0)-a/3.d0
       eigval(1)=2.d0*dsqrt(eps-q)*dcos((theta+2.d0*pai)/3.d0)-a/3.d0
       eigval(2)=2.d0*dsqrt(eps-q)*dcos((theta+4.d0*pai)/3.d0)-a/3.d0
c       write(*,*) eigval(0), eigval(1), eigval(2)

c      Sort eigenvalues: eigval(0)<eigval(1)<eigval(2)
       if ( eigval(0) .gt. eigval(1) )    then
         tmp = eigval(0)
         eigval(0) = eigval(1)
         eigval(1) = tmp
       endif
       if ( eigval(1) .gt. eigval(2) ) then
         tmp = eigval(1)
         eigval(1) = eigval(2)
         eigval(2) = tmp
       endif
       if ( eigval(0) .gt. eigval(1) ) then
         tmp = eigval(0)
         eigval(0) = eigval(1)
         eigval(1) = tmp
       endif
c
c      write(*,*) eigval(0),eigval(1),eigval(2),q,theta
       return
       end

c-0-------------------------------------------------------------------------------
! 1D
       subroutine SpectralU1Dy(ux,timestep,idata)
       use mesh
       use mesh_common
       use properties    
       use diags_io
#ifdef PARALLEL
      use mpistuff
#endif
c
       double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
       integer:: timestep, idata
c
       double precision:: uc(nxlocal,nylocal,nzlocal,nvar)
       double precision:: uyfft(ny,nvar),uy(nylocal,nvar),uytmp(ny,nvar)
       double precision:: ef(ny)
       integer:: i,j,k,l, nyfft, nsize, idir     
       integer:: ispectral, kspectral, iplocal, ilocal
       double complex::   utmp(ny), uftmp(ny), uf(ny,nvar)
c
       nyfft = ny
c
c       ek = 0.d0
c
       do i=1, nxlocal, 1
         do j=1, nylocal, 1
           do k=1, nzlocal, 1
           uc(i,j,k,1)=0.5d0*(ux(i,j,k,1)+ux(i+1,j,k,1))
           uc(i,j,k,2)=0.5d0*(ux(i,j,k,2)+ux(i,j+1,k,2))
           uc(i,j,k,3)=0.5d0*(ux(i,j,k,3)+ux(i,j,k+1,3))
           enddo
         enddo
       enddo
c
       ispectral = idata
       kspectral = 5  
       iplocal = int(idata/nxlocal)+1
       ilocal  = mod(idata,nxlocal)
c
      uy = 0.d0
      uytmp = 0.d0
      if(iprocx.eq. iplocal .and. iprocz.eq.1) then
       do j=1, nylocal, 1
         uy(j,1:3)= uc(ilocal,j,kspectral,1:3)
         uytmp((iprocy-1)*nylocal+j,1:3)=uy(j,1:3)
       enddo
      endif

       uyfft = 0.d0
       nsize = ny*3
#ifdef PARALLEL
       call MPI_Allreduce(uytmp, uyfft, nsize,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
              Call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
       uyfft = uytmp
#endif
c       
c      do j=1, ny, 1
c          write(*,*) j, uytmp(j,1),uytmp(j,2),uytmp(j,3),iprocx,iprocz
c          write(*,*) j, uyfft(j,1),uyfft(j,2),uyfft(j,3),iprocx,iprocz
c      enddo  
c
      if(iprocx.eq.1 .and. iprocy.eq.1 .and. iprocz.eq.1) then
       idir = 1
       do l=1, nvar, 1
         utmp(:) = uyfft(:,l)
c         call FFT1D(utmp,uftmp,ny,idir)   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         uf(:,l)=uftmp(:)
       enddo
       do j=1, ny, 1
         ef(j) = ( uf(j,1)*conjg(uf(j,1)) +
     &             uf(j,2)*conjg(uf(j,2)) +
     &             uf(j,3)*conjg(uf(j,3)) )/2.d0
       enddo
       write(ispy01,*) 'output at time :', timestep
       do j=2, ny/2, 1
         write(ispy01,*) j, ef(j)
       enddo
c       do j=1, ny, 1
c         write(ispy01,*) uyfft(j,1),uyfft(j,2),uyfft(j,3)
c       enddo
      endif ! iproc=1
c
c  check output
c       write(*,*) 'Max/Min ef:', maxval(ef),minval(ef)
c
       return
       end



c-------------------------------------------------------------------------------
c
c----------------------------------------------------------------------
c   Part 4   
c        Diagnostics for output results in y-z plane
c----------------------------------------------------------------------
c--------------------------------------------------------------------------
      subroutine output_xslice(ux,timestep,ndiag)
      use time
      use mesh
      use mesh_common
      use virtualwallBC
      use stats
      use diags_io
      use properties
      use GClaw
#ifdef PARALLEL
      use mpistuff
#endif
c
      integer:: timestep
      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision:: kappat, uept, hfact, utaut, Rethet
      integer:: idata, i,j,k, iplocal, ilocal, ndiag
      double precision:: nu, ut, vt, wt, et, g_et, spany, g_spany
c
      nu = 1.d0/Reynolds
c
      idata   = ioutput
      iplocal = int(idata/nxlocal)+1
      ilocal  = mod(idata,nxlocal)
c
      kappat = 0.d0
      uept = 0.d0
      hfact = 0.d0
      utaut = 0.d0
      Rethet = 0.d0
c
c-------------------------------------------
c  itestv :  total kinetic energy
c-------------------------------------------
      call kineticE(ux,timestep)
      call kineticEf(ux,timestep)
c
c-------------------------------------------
c
c--------------------------------------------------
      spany = 0.d0 
      do i=1, nxlocal, 1
        do j=1, nylocal, 1
          do k=1, nzlocal, 1
            spany = spany+ux(i,j,k,2) 
          enddo
        enddo
      enddo
#ifdef PARALLEL
       g_spany = 0.d0
       call MPI_Allreduce(spany, g_spany, 1,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
       call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
       g_spany = spany
#endif
      if(iprocx.eq.1 .and. iprocy.eq.1 .and. iprocz.eq.1) then
      write(itestout01,*) ttot, g_spany
      endif
c
c-------------------------------------------
      if(iprocx.eq.1 .and. iprocy.eq.1 .and. iprocz.eq.1) then
c
c-------------------------------------------
c  itestv1 : x - Re_\theta relation
c-------------------------------------------
       if(mod(timestep,ndiag).eq.1) then
       write(itestv1,*)'write data at', ttot
       do i=1, nx, 1
         Rethet = g_theta(i)/nu
         uept   = 1.d0/g_utau(i)
         write(itestv1,*) i, g_utau(i),g_vwall(i),
     &           g_del99(i),g_delta(i),g_theta(i)
       enddo
       endif

c-------------------------------------------
c
c-------------------------------------------
       if(mod(timestep,ndiag).eq.1) then
       write(itestv5,*)'write data at', ttot
       do i=1, nx, 1
         write(itestv5,*) i, g_kappa(i),g_GCci(i),
     &           g_GCco(i),g_GCgam(i),g_GCpre(i)
       enddo
       endif


c-------------------------------------------
c  itestv2 : utau at inlet
c-------------------------------------------
       utaut = g_utau(1)
       write(itestv2,*) ttot, utaut
c
c-------------------------------------------
c  itestv3 : Re_\theta at the selected plane
c------------------------------------------- 
       Rethet = g_theta(idata)/nu
       write(itestv3,*) ttot, g_kappa(idata)
c
c-------------------------------------------
c  itestv4 : Ue plus
c-------------------------------------------
       uept = 1.d0/g_utau(idata)
       write(itestv4,*) ttot, uept
c
c-------------------------------------------
c  itestv6: H factor
c-------------------------------------------
       hfact = g_delta(idata)/g_theta(idata)
       write(itestv6,*) ttot, hfact
c-------------------------------------------

      endif
c
c-----------------------------------------
c  itestv5, kappa
c------------------------------------------
c       kappat= 0.d0
c       if(iprocx.eq.iplocal .and. iprocz.eq.1) then
c       do j=1, nylocal, 1
c         kappat = kappat + bc_kappa(ilocal,j)/nylocal
c       enddo 
c       endif
c#ifdef PARALLEL
cc       if(yprocs.gt.1) then
c       g_kappa = 0.d0
c       call MPI_Allreduce(kappat, g_kappa, 1,
c     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
c              Call ErrorHandler(ierr,ERROR_ALLREDUCE)
c       g_kappa = g_kappa/yprocs
c       else
c       g_kappa = kappat
c       endif
c#else
c       g_kappa = kappat
c#endif   

c
c-------------------------------------------
c  out:  output velocity relevant.
c--------------------------------------------
      if(mod(timestep,ndiag).eq.1) then
        call outputxsliceuvw(timestep,idata)
        call SpectralU1Dy(ux,timestep,idata)
      endif
c
      return
      end
c
c----------------------------------------------------------------------
c     instantaneous kinetic energy
c----------------------------------------------------------------------
      subroutine kineticE(ux,timestep)
      use time
      use mesh
      use mesh_common
      use properties
      use diags_io
#ifdef PARALLEL
      use mpistuff
#endif
c
      integer:: timestep
      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      integer:: i,j,k, ilo, ihi, jlo, jhi, klo, khi
      double precision:: ut, vt, wt, ep, ef ! epoint, eflux
      double precision:: et, g_et, gl_et, gr_et, gt_et, gb_et,gd_et
      double precision:: dudx(3,3), St(3,3)
      double precision:: nu
      
      nu = 1.d0/Reynolds
c
c     total energy
       et = 0.d0
       do i=1, nxlocal, 1
         do j=1, nylocal, 1
           do k=1, nzlocal, 1
           ut = 0.5d0*(ux(i,j,k,1)+ux(i+1,j,k,1))
           vt = 0.5d0*(ux(i,j,k,2)+ux(i,j+1,k,2))
           wt = 0.5d0*(ux(i,j,k,3)+ux(i,j,k+1,3))
           et =et + 0.5d0*dsqrt(ut*ut+vt*vt+wt*wt)*dx*dy*dz
c            write(*,*) i,j,k,et
           enddo
         enddo
       enddo
#ifdef PARALLEL
       g_et = 0.d0
       call MPI_Allreduce(et, g_et, 1,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
              Call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
       g_et = et
#endif
c
c     viscous dissipation
      ef = 0.d0
      ilo = 1
      ihi = nxlocal
      klo = 1
      khi = nzlocal
      if(iprocx.eq.1) ilo = 2
      if(iprocx.eq.xprocs) ihi = nxlocal-1
      if(iprocz.eq.1) klo = 2
      if(iprocz.eq.zprocs) khi = nzlocal-1
      do i=ilo, ihi, 1
        do j=1, nylocal, 1
          do k=klo,khi, 1
            dudx(1,1)=(ux(i+1,j,k,1)-ux(i,j,k,1))*dxi
            dudx(2,2)=(ux(i,j+1,k,2)-ux(i,j,k,2))*dyi
            dudx(3,3)=(ux(i,j,k+1,3)-ux(i,j,k,3))*dzi
            dudx(1,2)=((ux(i+1,j+1,k,1)+ux(i,j+1,k,1))/2.d0-
     &                 (ux(i+1,j-1,k,1)+ux(i,j-1,k,1))/2.d0)/2.d0*dyi
            dudx(1,3)=((ux(i+1,j,k+1,1)+ux(i,j,k+1,1))/2.d0-
     &                 (ux(i+1,j,k-1,1)+ux(i,j,k-1,1))/2.d0)/2.d0*dzi
            dudx(2,1)=((ux(i+1,j+1,k,2)+ux(i+1,j,k,2))/2.d0-
     &                 (ux(i-1,j+1,k,2)+ux(i-1,j,k,2))/2.d0)/2.d0*dxi
            dudx(2,3)=((ux(i,j+1,k+1,2)+ux(i,j,k+1,2))/2.d0-
     &                 (ux(i,j+1,k-1,2)+ux(i,j,k-1,2))/2.d0)/2.d0*dzi
            dudx(3,1)=((ux(i+1,j,k+1,3)+ux(i+1,j,k,3))/2.d0-
     &                 (ux(i-1,j,k+1,3)+ux(i-1,j,k,3))/2.d0)/2.d0*dxi
            dudx(3,2)=((ux(i,j+1,k+1,3)+ux(i,j+1,k,3))/2.d0-
     &                 (ux(i,j-1,k+1,3)+ux(i,j-1,k,3))/2.d0)/2.d0*dyi 
            St(1,1) = 0.5d0* ( dudx(1,1)+dudx(1,1) )
            St(2,2) = 0.5d0* ( dudx(2,2)+dudx(2,2) )
            St(3,3) = 0.5d0* ( dudx(3,3)+dudx(3,3) )
            St(1,2) = 0.5d0* ( dudx(1,2)+dudx(2,1) )
            St(2,3) = 0.5d0* ( dudx(2,3)+dudx(3,2) )
            St(1,3) = 0.5d0* ( dudx(1,3)+dudx(3,1) )
            St(2,1) = St(1,2)
            St(3,1) = St(1,3)
            St(3,2) = St(2,3) 
            ef = ef + 2.d0*nu*dx*dy*dz*( 
     &          St(1,1)**2.d0+St(1,2)**2.d0+St(1,3)**2.d0+
     &          St(2,1)**2.d0+St(2,2)**2.d0+St(2,3)**2.d0+
     &          St(3,1)**2.d0+St(3,2)**2.d0+St(3,3)**2.d0 )
          enddo
        enddo
      enddo
#ifdef PARALLEL
       gd_et = 0.d0
       call MPI_Allreduce(ef, gd_et, 1,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
              Call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
       gd_et = ef
#endif

c
c     kinetic energy flux at inflow
      ef = 0.d0
      if(iprocx.eq.1) then
        i=2
        jlo = 1
        jhi = nylocal
        klo = 1
        khi = nzlocal
        if(iprocz.eq.1) klo = 2
        if(iprocz.eq.zprocs) khi = nzlocal-1 
        do j=jlo, jhi, 1
          do k=klo, khi, 1
            ut = ux(i,j,k,1)
            vt = 0.25d0*(ux(i,j  ,k,2)+ux(i-1,j  ,k,2)
     &                  +ux(i,j+1,k,2)+ux(i-1,j+1,k,2))
            wt = 0.25d0*(ux(i,j,k  ,3)+ux(i-1,j,k  ,3)
     &                  +ux(i,j,k+1,3)+ux(i-1,j,k+1,3))
            ep = 0.5d0* dsqrt(ut*ut+vt*vt+wt*wt)*dy*dz
            ef = ef + ut*ep
          enddo
        enddo
      endif
#ifdef PARALLEL
       gl_et = 0.d0
       call MPI_Allreduce(ef, gl_et, 1,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
              Call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
       gl_et = ef
#endif
c
c     kinetic energy flux at outflow
      ef = 0.d0
      if(iprocx.eq.xprocs) then
        i = nxlocal
        jlo = 1
        jhi = nylocal
        klo = 1
        khi = nzlocal
        if(iprocz.eq.1) klo = 2
        if(iprocz.eq.zprocs) khi = nzlocal-1
        do j=jlo, jhi, 1
          do k=klo, khi, 1
            ut = ux(i,j,k,1)
            vt = 0.25d0*(ux(i,j  ,k,2)+ux(i-1,j  ,k,2)
     &                  +ux(i,j+1,k,2)+ux(i-1,j+1,k,2))
            wt = 0.25d0*(ux(i,j,k  ,3)+ux(i-1,j,k  ,3)
     &                  +ux(i,j,k+1,3)+ux(i-1,j,k+1,3))
            ep = 0.5d0* dsqrt(ut*ut+vt*vt+wt*wt)*dy*dz
            ef = ef + ut*ep          
          enddo
        enddo
      endif
#ifdef PARALLEL
       gr_et = 0.d0
       call MPI_Allreduce(ef, gr_et, 1,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
              Call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
       gr_et = ef
#endif
c
c     kinetic energy flux at bottom
      ef = 0.d0
      if(iprocz.eq.1) then
        k = 2
        jlo = 1
        jhi = nylocal
        ilo = 1
        ihi = nxlocal
        if(iprocx.eq.1) ilo = 2
        if(iprocx.eq.xprocs) ihi = nxlocal-1
        do j=jlo, jhi, 1
          do i=ilo, ihi, 1
            ut = 0.25d0*(ux(i  ,j,k,1)+ux(i  ,j,k-1,1)
     &                  +ux(i-1,j,k,1)+ux(i-1,j,k-1,1))
            vt = 0.25d0*(ux(i,j  ,k,2)+ux(i,j  ,k-1,2)
     &                  +ux(i,j+1,k,2)+ux(i,j+1,k-1,2))
            wt = ux(i,j,k,3)
            ep = 0.5d0* dsqrt(ut*ut+vt*vt+wt*wt)*dy*dx
            ef = ef + wt*ep          
          enddo
        enddo
      endif
#ifdef PARALLEL
       gb_et = 0.d0
       call MPI_Allreduce(ef, gb_et, 1,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
              Call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
       gb_et = ef
#endif      
c
c      kinetic energy flux at top
      ef = 0.d0
      if(iprocz.eq.zprocs) then
        k = nzlocal
        jlo = 1
        jhi = nylocal
        ilo = 1
        ihi = nxlocal
        if(iprocx.eq.1) ilo = 2
        if(iprocx.eq.xprocs) ihi = nxlocal-1
        do j=jlo, jhi, 1
          do i=ilo, ihi, 1
            ut = 0.25d0*(ux(i  ,j,k,1)+ux(i  ,j,k-1,1)
     &                  +ux(i-1,j,k,1)+ux(i-1,j,k-1,1))
            vt = 0.25d0*(ux(i,j  ,k,2)+ux(i,j  ,k-1,2)
     &                  +ux(i,j+1,k,2)+ux(i,j+1,k-1,2))
            wt = ux(i,j,k,3)
            ep = 0.5d0* dsqrt(ut*ut+vt*vt+wt*wt)*dy*dx
            ef = ef + wt*ep          
          enddo
        enddo
      endif
#ifdef PARALLEL
       gt_et = 0.d0
       call MPI_Allreduce(ef, gt_et, 1,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
              Call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
       gt_et = ef
#endif          


c
c     output
      if(iprocx.eq.1 .and. iprocy.eq.1 .and. iprocz.eq.1) then
       write(itestv,*)ttot,g_et,gl_et,-gr_et,gb_et,-gt_et, gd_et
      endif
c
      return
      end
c
c----------------------------------------------------------------------
c     turbulent kinetic energy
c----------------------------------------------------------------------
      subroutine kineticEf(ux,timestep)
      use time
      use mesh
      use mesh_common
      use properties
      use diags_io
      use stats
#ifdef PARALLEL
      use mpistuff
#endif
c
      integer:: timestep
      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision:: ufl_u1(nxlocal,nylocal,nzlocal,3)
      integer:: i,j,k, ilo, ihi, jlo, jhi, klo, khi
      double precision:: ut, vt, wt, ep, ef ! epoint, eflux
      double precision:: et, g_et, gl_et, gr_et, gt_et, gb_et,gd_et
      double precision:: dudx(3,3), St(3,3)
      double precision:: nu
      
      nu = 1.d0/Reynolds
c
      ufl_u1 = 0.d0
      et = 0.d0
      do i=1, nxlocal, 1
        do j=1, nylocal, 1
          do k=1, nzlocal, 1
            ! velocity fluctuation at center
            ufl_u1(i,j,k,1)= (
     &            ux(i  ,j,k,1)-st_umean(i  ,k,1) +
     &            ux(i+1,j,k,1)-st_umean(i+1,k,1)  )/2.d0
            ufl_u1(i,j,k,2)= (
     &            ux(i  ,j,k,2)-st_umean(i  ,k,2) +
     &            ux(i,j+1,k,2)-st_umean(i  ,k,2)  )/2.d0
            ufl_u1(i,j,k,3)= (
     &            ux(i,j,k  ,3)-st_umean(i,k  ,3) +
     &            ux(i,j,k+1,3)-st_umean(i,k+1,3)  )/2.d0
           et =et + 0.5d0*dsqrt(
     &         ufl_u1(i,j,k,1)*ufl_u1(i,j,k,1)+
     &         ufl_u1(i,j,k,2)*ufl_u1(i,j,k,2)+
     &         ufl_u1(i,j,k,3)*ufl_u1(i,j,k,3)
     &                          )*dx*dy*dz
          enddo
        enddo
      enddo
#ifdef PARALLEL
       g_et = 0.d0
       call MPI_Allreduce(et, g_et, 1,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
              Call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
       g_et = et
#endif
c
c     turbulent viscous dissipation
      ef = 0.d0
      ilo = 2
      ihi = nxlocal-1
      klo = 2
      khi = nzlocal-1
      if(iprocx.eq.1) ilo = 2
      if(iprocx.eq.xprocs) ihi = nxlocal-1
      if(iprocz.eq.1) klo = 2
      if(iprocz.eq.zprocs) khi = nzlocal-1
      do i=ilo, ihi, 1
        do j=2, nylocal-1, 1
          do k=klo,khi, 1
            dudx(1,1)=( (ux(i+1,j,k,1)-st_umean(i+1,k,1))-
     &                  (ux(i  ,j,k,1)-st_umean(i  ,k,1)) )*dxi
            dudx(2,2)=( (ux(i,j+1,k,2)-st_umean(i  ,k,2))-
     &                  (ux(i,j  ,k,2)-st_umean(i  ,k,2)) )*dyi
            dudx(3,3)=( (ux(i,j,k+1,3)-st_umean(i,k+1,3))-
     &                  (ux(i,j,k  ,3)-st_umean(i,k  ,3)) )*dzi
c            dudx(1,2)=(ufl_u1(i,j+1,k,1)-ufl_u1(i,j-1,k,1))/2.d0*dyi
c            dudx(1,3)=(ufl_u1(i,j,k+1,1)-ufl_u1(i,j,k-1,1))/2.d0*dzi
c            dudx(2,1)=(ufl_u1(i+1,j,k,2)-ufl_u1(i-1,j,k,2))/2.d0*dxi
c            dudx(2,3)=(ufl_u1(i,j,k+1,2)-ufl_u1(i,j,k-1,2))/2.d0*dzi
c            dudx(3,1)=(ufl_u1(i+1,j,k,3)-ufl_u1(i-1,j,k,3))/2.d0*dxi
c            dudx(3,2)=(ufl_u1(i,j+1,k,3)-ufl_u1(i,j-1,k,3))/2.d0*dyi 

            dudx(1,2)=( (ux(i+1,j+1,k,1)-st_umean(i+1,k,1)+
     &                   ux(i  ,j+1,k,1)-st_umean(i  ,k,1))/2.d0-
     &                  (ux(i+1,j-1,k,1)-st_umean(i+1,k,1)+
     &                   ux(i  ,j-1,k,1)-st_umean(i  ,k,1))/2.d0
     &                )/2.d0*dyi
            dudx(1,3)=( (ux(i+1,j,k+1,1)-st_umean(i+1,k+1,1)+
     &                   ux(i  ,j,k+1,1)-st_umean(i  ,k+1,1))/2.d0-
     &                  (ux(i+1,j,k-1,1)-st_umean(i+1,k-1,1)+
     &                   ux(i  ,j,k-1,1)-st_umean(i  ,k-1,1))/2.d0 
     &                )/2.d0*dzi
            dudx(2,1)=( (ux(i+1,j+1,k,2)-st_umean(i+1,k,2)+
     &                   ux(i+1,j  ,k,2)-st_umean(i+1,k,2))/2.d0-
     &                  (ux(i-1,j+1,k,2)-st_umean(i-1,k,2)+
     &                   ux(i-1,j  ,k,2)-st_umean(i-1,k,2))/2.d0
     &                )/2.d0*dxi
            dudx(2,3)=( (ux(i,j+1,k+1,2)-st_umean(i,k+1,2)+
     &                   ux(i,j  ,k+1,2)-st_umean(i,k+1,2))/2.d0-
     &                  (ux(i,j+1,k-1,2)-st_umean(i,k-1,2)+
     &                   ux(i,j  ,k-1,2)-st_umean(i,k-1,2))/2.d0
     &                )/2.d0*dzi
            dudx(3,1)=( (ux(i+1,j,k+1,3)-st_umean(i+1,k+1,3)+
     &                   ux(i+1,j,k  ,3)-st_umean(i+1,k  ,3))/2.d0-
     &                  (ux(i-1,j,k+1,3)-st_umean(i-1,k+1,3)+
     &                   ux(i-1,j,k  ,3)-st_umean(i-1,k  ,3))/2.d0
     &                )/2.d0*dxi
            dudx(3,2)=( (ux(i,j+1,k+1,3)-st_umean(i,k+1,3)+
     &                   ux(i,j+1,k  ,3)-st_umean(i,k  ,3))/2.d0-
     &                  (ux(i,j-1,k+1,3)-st_umean(i,k+1,3)+
     &                   ux(i,j-1,k  ,3)-st_umean(i,k  ,3))/2.d0
     &                )/2.d0*dyi 

            St(1,1) = 0.5d0* ( dudx(1,1)+dudx(1,1) )
            St(2,2) = 0.5d0* ( dudx(2,2)+dudx(2,2) )
            St(3,3) = 0.5d0* ( dudx(3,3)+dudx(3,3) )
            St(1,2) = 0.5d0* ( dudx(1,2)+dudx(2,1) )
            St(2,3) = 0.5d0* ( dudx(2,3)+dudx(3,2) )
            St(1,3) = 0.5d0* ( dudx(1,3)+dudx(3,1) )
            St(2,1) = St(1,2)
            St(3,1) = St(1,3)
            St(3,2) = St(2,3) 
            ef = ef + 2.d0*nu*dx*dy*dz*( 
     &          St(1,1)**2.d0+St(1,2)**2.d0+St(1,3)**2.d0+
     &          St(2,1)**2.d0+St(2,2)**2.d0+St(2,3)**2.d0+
     &          St(3,1)**2.d0+St(3,2)**2.d0+St(3,3)**2.d0 )
          enddo
        enddo
      enddo
#ifdef PARALLEL
       gd_et = 0.d0
       call MPI_Allreduce(ef, gd_et, 1,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
              Call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
       gd_et = ef
#endif
c
c     output
      if(iprocx.eq.1 .and. iprocy.eq.1 .and. iprocz.eq.1) then
       write(itestv,*)ttot,g_et, gd_et
      endif
c
      return
      end
c----------------------------------------------------------------------
c     output mean velocity and fluctuation related values
c----------------------------------------------------------------------
      subroutine outputxsliceuvw(timestep,idata)
      use time
      use mesh
      use mesh_common
      use virtualwallBC
      use stats
      use diags_io
      use properties
#ifdef PARALLEL
      use mpistuff
#endif
c
      integer:: timestep
c      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      integer:: idata, i,j,k,l, iplocal, ilocal
c
      double precision:: g_umean(nx,nz,3), gt_umean(nx,nz,3)
      double precision:: g_zc(nz), gt_zc(nz)
      double precision:: zp, up, vp, wp
      double precision:: g_bc(nx,3), gt_bc(nx,3)
c
      double precision:: g_u2(nx,nz),g_u3(nx,nz),g_u4(nx,nz)
      double precision:: g_uv(nx,nz),g_uw(nx,nz),g_vv(nx,nz),g_vw(nx,nz)
      double precision:: g_ww(nx,nz), gta(nx,nz), gt(nx,nz),g_T(nx,nz,6)
c      double precision:: l_utauc(nzlocal),g_utauc(nz)
c
      character*50 DumpFile

c
c    for others
       double precision:: nu
c
      nu = 1.d0/Reynolds
c
      iplocal = int(idata/nxlocal)+1
      ilocal  = mod(idata,nxlocal)

c---------------------
c   get z coordinates
c--------------------------------------------
        gt_zc = 0.d0
        if(iprocx.eq.iplocal .and. iprocy.eq.1) then
          gt_zc((iprocz-1)*nzlocal+1:iprocz*nzlocal) = zc(1:nzlocal) 
        endif
#ifdef PARALLEL
        g_zc = 0.d0
        call MPI_Allreduce(gt_zc, g_zc, nz,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
        Call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
        g_zc = gt_zc
#endif  
c---------------------------
c    get velocity fluctuation
c-----------------------------------------

c----------------------------


#ifdef PARALLEL
       gt = 0.d0
       if( iprocy.eq.1.and.iprocz.eq.1 ) then
         gt((iprocx-1)*nxlocal+1:iprocx*nxlocal,1:nz)=
     &                          g_fu2(1:nxlocal,1:nz)
       endif
       call MPI_Allreduce(gt, gta, nx*nz,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
              Call ErrorHandler(ierr,ERROR_ALLREDUCE)
        g_u2 = gta
c
       gt = 0.d0
       if( iprocy.eq.1.and.iprocz.eq.1 ) then
         gt((iprocx-1)*nxlocal+1:iprocx*nxlocal,1:nz)=
     &                          g_fu3(1:nxlocal,1:nz)
       endif
       call MPI_Allreduce(gt, gta, nx*nz,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
              Call ErrorHandler(ierr,ERROR_ALLREDUCE)
        g_u3 = gta
c
       gt = 0.d0
       if( iprocy.eq.1.and.iprocz.eq.1 ) then 
         gt((iprocx-1)*nxlocal+1:iprocx*nxlocal,1:nz)=
     &                          g_fu4(1:nxlocal,1:nz) 
       endif
       call MPI_Allreduce(gt, gta, nx*nz,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
              Call ErrorHandler(ierr,ERROR_ALLREDUCE)
        g_u4 = gta
c
       gt = 0.d0
       if( iprocy.eq.1.and.iprocz.eq.1 ) then
         gt((iprocx-1)*nxlocal+1:iprocx*nxlocal,1:nz)=
     &                          g_fuv(1:nxlocal,1:nz) 
       endif
       call MPI_Allreduce(gt, gta, nx*nz,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
              Call ErrorHandler(ierr,ERROR_ALLREDUCE)
        g_uv = gta
c
       gt = 0.d0
       if( iprocy.eq.1.and.iprocz.eq.1 ) then
         gt((iprocx-1)*nxlocal+1:iprocx*nxlocal,1:nz)=
     &                          g_fuw(1:nxlocal,1:nz) 
       endif
       call MPI_Allreduce(gt, gta, nx*nz,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
              Call ErrorHandler(ierr,ERROR_ALLREDUCE)
        g_uw = gta
c
       gt = 0.d0
       if( iprocy.eq.1.and.iprocz.eq.1 ) then
         gt((iprocx-1)*nxlocal+1:iprocx*nxlocal,1:nz)=
     &                          g_fvv(1:nxlocal,1:nz) 
       endif
       call MPI_Allreduce(gt, gta, nx*nz,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
              Call ErrorHandler(ierr,ERROR_ALLREDUCE)
        g_vv = gta
c
       gt = 0.d0
       if( iprocy.eq.1.and.iprocz.eq.1 ) then
         gt((iprocx-1)*nxlocal+1:iprocx*nxlocal,1:nz)=
     &                          g_fvw(1:nxlocal,1:nz) 
       endif
       call MPI_Allreduce(gt, gta, nx*nz,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
              Call ErrorHandler(ierr,ERROR_ALLREDUCE)
        g_vw = gta
c
       gt = 0.d0
       if( iprocy.eq.1.and.iprocz.eq.1) then
         gt((iprocx-1)*nxlocal+1:iprocx*nxlocal,1:nz)=
     &                          g_fww(1:nxlocal,1:nz) 
       endif
       call MPI_Allreduce(gt, gta, nx*nz,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
              Call ErrorHandler(ierr,ERROR_ALLREDUCE)
        g_ww = gta
c
        do l=1, 6
        gt = 0.d0
        if( iprocy.eq.1.and.iprocz.eq.1)  then
          gt((iprocx-1)*nxlocal+1:iprocx*nxlocal,1:nz)=
     &                        g_Ten(1:nxlocal,1:nz,l)
        endif
        call MPI_Allreduce(gt, gta, nx*nz,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
               Call ErrorHandler(ierr,ERROR_ALLREDUCE)
        g_T(:,:,l) = gta
        enddo

#else
        g_u2 = g_fu2
        g_u3 = g_fu3
        g_u4 = g_fu4
        g_uv = g_fuv
        g_uw = g_fuw
        g_vv = g_fvv
        g_vw = g_fvw
        g_ww = g_fww
        g_T  = g_Ten
c
#endif


c--------------------------------------------
c  get  mean u velocity
c--------------------------------------------
c output umean
      g_umean = 0.d0
      gt_umean = 0.d0
        if( iprocy.eq.1) then
        gt_umean((iprocx-1)*nxlocal+1:iprocx*nxlocal,
     &           (iprocz-1)*nzlocal+1:iprocz*nzlocal,1:3) = 
     &     st_umean(1:nxlocal,1:nzlocal,1:3)
        endif      
#ifdef PARALLEL
        call MPI_Allreduce(gt_umean, g_umean, nx*nz*3,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
        Call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
        g_umean = gt_umean
#endif 
!  virtuwal wall point, for  output   
      g_bc = 0.d0
      gt_bc= 0.d0
      if(  iprocy.eq.1 .and. iprocz.eq.1) then
        do i=1, nxlocal, 1
        do j=1, nylocal, 1
          gt_bc((iprocx-1)*nxlocal+i,1:3)=
     &    gt_bc((iprocx-1)*nxlocal+i,1:3) + bc_uwall(i,j,1:3)/nylocal
        enddo
        enddo
      endif
#ifdef PARALLEL
        call MPI_Allreduce(gt_bc, g_bc, nx*3,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
        Call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
        g_bc = gt_bc
#endif   
       
c----------------------------------------------------------------------
c         output part
c----------------------------------------------------------------------

c--------------------------------------------
c iout1: mean u velocity
c--------------------------------------------
      if(iprocx.eq.1 .and. iprocy.eq.1 .and. iprocz.eq.1) then
      write(DumpFile, 1101) float(timestep) / 1000000
1101  format('sum', f7.6)
      open(803,file=DumpFile)
        write(803,*) 'output at time :', ttot, timestep
        do i=1, nx, 1
          write(803,*) 'x position :', i
          zp = bc_h0dz*dz*g_utau(i)*Reynolds 
          up = g_bc(i,1)/g_utau(i)
          vp = g_bc(i,2)/g_utau(i)
          wp = g_bc(i,3)/g_utau(i)        
          write(803,*) zp,up,vp,wp
          do k=1, nz, 1
            zp = (g_zc(k)+bc_h0dz*dz)*g_utau(i)*Reynolds
            up = g_umean(i,k,1)/g_utau(i)
            vp = g_umean(i,k,2)/g_utau(i)
            wp = g_umean(i,k,3)/g_utau(i)
            write(803,*) zp, up, vp, wp
c          write(*,*) g_zc(k), up 
          enddo
        enddo
      close(803)  
      endif
cc
c---------------------------------------
c     fluctuation
c-----------------------------------------
      if(iprocx.eq.1 .and. iprocy.eq.1 .and. iprocz.eq.1) then
      write(DumpFile, 1102) float(timestep) / 1000000
1102  format('rms', f7.6)
      open(804,file=DumpFile)
        write(804,*) 'output uu at time :',  timestep
        do i=1, nx, 1
           write(804,*) 'x position :', i
        do k=1, nz, 1
          write(804,'(15f18.8)') g_zc(k),g_u2(i,k),g_T(i,k,1),
     &                   g_uv(i,k),g_T(i,k,2),g_uw(i,k),g_T(i,k,3),
     &                   g_vv(i,k),g_T(i,k,4),g_vw(i,k),g_T(i,k,5),
     &                   g_ww(i,k),g_T(i,k,6),g_u3(i,k),g_u4(i,k)  
        enddo 
        enddo
      close(804)
      endif  

c
      return
      end


