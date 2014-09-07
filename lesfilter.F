c----------------------------------------------------------------------
c     Wan Cheng 
c     August 2012
c----------------------------------------------------------------------
c
c
c----------------------------------------------------------------------
      subroutine fftnonlinear(Adu)
      use mesh
      use mesh_common
c
      implicit none
      double precision:: Adu(nxlocal+1,nylocal+1,nzlocal+1,nvar)
      double precision:: adu1d(1:nylocal)
      double precision:: g_adu1d(1:ny)
      integer:: i,j,k,l
c
      adu1d   = 0.d0
      g_adu1d = 0.d0

c    u part
      do i=1, nxlocal+1, 1
        do k=1, nzlocal, 1
          do j=1, nylocal, 1
            adu1d(j)= Adu(i,j,k,1)
          enddo
          call mpiysum(adu1d,g_adu1d)
          call trancate(g_adu1d)
          do j=1, nylocal, 1
            Adu(i,j,k,1) = g_adu1d((iprocy-1)*nylocal+j)
          enddo
        enddo
      enddo
c   v part
      do i=1, nxlocal, 1
        do k=1, nzlocal, 1
          do j=1, nylocal, 1
            adu1d(j)= Adu(i,j,k,2)
          enddo
          call mpiysum(adu1d,g_adu1d)
          call trancate(g_adu1d)
          if(iprocy.lt.yprocs) then
            do j=1, nylocal+1, 1
              Adu(i,j,k,2) = g_adu1d((iprocy-1)*nylocal+j)
            enddo
          else
            do j=1, nylocal, 1
              Adu(i,j,k,2) = g_adu1d((iprocy-1)*nylocal+j)
            enddo
            Adu(i,nylocal+1,k,2) = g_adu1d(1)
          endif
        enddo
      enddo
c   w part
      do i=1, nxlocal, 1
        do k=1, nzlocal+1, 1
          do j=1, nylocal, 1
            adu1d(j)= Adu(i,j,k,3)
          enddo
          call mpiysum(adu1d,g_adu1d)
          call trancate(g_adu1d)
          do j=1, nylocal, 1
            Adu(i,j,k,3) = g_adu1d((iprocy-1)*nylocal+j)
          enddo
        enddo
      enddo
c
      return
      end
c--------------------
c
c---------------------
      subroutine ffttensor(Tijn)
      use mesh
      use mesh_common
c
      implicit none 
      double precision:: Tijn(0:nxlocal+2,0:nylocal+2,0:nzlocal+2,6)
      double precision:: Tijn1d(1:nylocal)
      double precision:: g_tijn1d(1:ny), g_t(1:ny)
      integer:: i,j,k,l, jlo, jhi
c     
      Tijn1d   = 0.d0
      g_tijn1d = 0.d0
      g_t      = 0.d0
c     
      do l=1, 6, 1
        do i=0, nxlocal+2, 1
          do k=0, nzlocal+2, 1
            do j=1, nylocal, 1
              Tijn1d(j)=Tijn(i,j,k,l)
            enddo
            call mpiysum(Tijn1d,g_tijn1d)
            call trancate(g_tijn1d)
            jlo = 0
            jhi = nylocal+2
            if(iprocy.eq.1) jlo = 1
            if(iprocy.eq.yprocs) jhi = nylocal
            do j=jlo, jhi, 1
              Tijn(i,j,k,l)=g_tijn1d((iprocy-1)*nylocal+j)
            enddo
            if(iprocy.eq.1) Tijn(i,0,k,l) = g_tijn1d(ny)
            if(iprocy.eq.yprocs) then
              Tijn(i,nylocal+1,k,l) = g_tijn1d(1)
              Tijn(i,nylocal+2,k,l) = g_tijn1d(2)
            endif
          enddo
        enddo
      enddo
c     
      return
      end        
c--------------------
c
c--------------------------
      subroutine fftbcnonlinear
      use mesh
      use mesh_common
      use virtualwallBC
c
      implicit none
      double precision:: adv1d(nylocal)
      double precision:: g_adv1d(ny), g_t(1:ny)
      integer:: i, j, k
c     
      adv1d   = 0.d0
      g_adv1d = 0.d0
      g_t     = 0.d0
c     
      do i=1, nxlocal+1, 1
        do j=1, nylocal, 1
          adv1d(j)=bc_convec(i,j)
        enddo
        call mpiysum(adv1d,g_adv1d)
        call trancate(g_adv1d)
        do j=1, nylocal, 1
          bc_convec(i,j)=g_adv1d((iprocy-1)*nylocal+j)
        enddo
      enddo
c
      return
      end
c--------------------
c
c--------------------------
      subroutine fftbcuw
      use mesh
      use mesh_common
      use virtualwallBC
c
      implicit none
      double precision:: uw1d(nylocal)
      double precision:: g_uw1d(ny), g_t(1:ny)
      integer:: i, j, k
c
      uw1d   = 0.d0
      g_uw1d = 0.d0
      g_t    = 0.d0
c
      do i=1, nxlocal+1, 1
        do j=1, nylocal, 1
          uw1d(j)=bc_uw(i,j)
        enddo
        call mpiysum(uw1d,g_uw1d)
        call trancate(g_uw1d)
        do j=1, nylocal, 1
          bc_uw(i,j)=g_uw1d((iprocy-1)*nylocal+j)
        enddo
      enddo
c
      return
      end
c--------------------
c
c--------------------------
      subroutine fftbcdudzh
      use mesh
      use mesh_common
      use virtualwallBC
c
      implicit none
      double precision:: dudzh1d(nylocal)
      double precision:: g_dudzh1d(ny), g_t(1:ny)
      integer:: i, j, k
c
      dudzh1d   = 0.d0
      g_dudzh1d = 0.d0
      g_t       = 0.d0
c
      do i=1, nxlocal+1, 1
        do j=1, nylocal, 1
          dudzh1d(j)=bc_dudzh(i,j)
        enddo
        call mpiysum(dudzh1d,g_dudzh1d)
        call trancate(g_dudzh1d)
        do j=1, nylocal, 1
          bc_dudzh(i,j)=g_dudzh1d((iprocy-1)*nylocal+j)
        enddo
      enddo
c
      return
      end
c--------------------
c
c--------------------------
      subroutine fftbcdpdx
      use mesh
      use mesh_common
      use virtualwallBC
c
      implicit none
      double precision:: dpdx1d(nylocal)
      double precision:: g_dpdx1d(ny), g_t(1:ny)
      integer:: i, j, k
c
      dpdx1d   = 0.d0
      g_dpdx1d = 0.d0
      g_t      = 0.d0
c
      do i=1, nxlocal+1, 1
        do j=1, nylocal, 1
          dpdx1d(j)=bc_dpdx(i,j)
        enddo
        call mpiysum(dpdx1d,g_dpdx1d)
        call trancate(g_dpdx1d)
        do j=1, nylocal, 1
          bc_dpdx(i,j)=g_dpdx1d((iprocy-1)*nylocal+j)
        enddo
      enddo
c
      return
      end

c--------------------
c
c--------------------------
      subroutine fftbcrhs
      use mesh
      use mesh_common
      use virtualwallBC
c
      implicit none
      double precision:: rhs1d(nylocal)
      double precision:: g_rhs1d(ny), g_t(1:ny)
      integer:: i, j, k
c
      rhs1d   = 0.d0
      g_rhs1d = 0.d0
      g_t     = 0.d0
c
      do i=1, nxlocal+1, 1
        do j=1, nylocal, 1
          rhs1d(j)=bc_rhs(i,j)
        enddo
        call mpiysum(rhs1d,g_rhs1d)
        call trancate(g_rhs1d)
        do j=1, nylocal, 1
          bc_rhs(i,j)=g_rhs1d((iprocy-1)*nylocal+j)
        enddo
      enddo
c
      return
      end
c--------------------
c
c--------------------------
      subroutine fftbcdudz
      use mesh
      use mesh_common
      use virtualwallBC
c
      implicit none
      double precision:: dudz1d(nylocal)
      double precision:: g_dudz1d(ny), g_t(1:ny)
      integer:: i, j, k
c
      dudz1d   = 0.d0
      g_dudz1d = 0.d0
      g_t      = 0.d0

c
      do i=1, nxlocal+1, 1
        do j=1, nylocal, 1
          dudz1d(j)=bc_dudz(i,j)
        enddo
        call mpiysum(dudz1d,g_dudz1d)
        call trancate(g_dudz1d)
        do j=1, nylocal, 1
          bc_dudz(i,j)=g_dudz1d((iprocy-1)*nylocal+j)
        enddo
      enddo
c
      return
      end
c--------------------
c
c--------------------------
      subroutine fftbcuwall
      use mesh
      use mesh_common
      use virtualwallBC
c
      implicit none
      double precision:: uwall1d(nylocal)
      double precision:: g_uwall1d(ny), g_t(1:ny)
      integer:: i, j, k, l
c
      uwall1d   = 0.d0
      g_uwall1d = 0.d0
      g_t       = 0.d0
c
      do l = 1, 3 
c         do i=1, nxlocal+1, 1
        do i=ixlo, ixhi, 1
          do j=1, nylocal, 1
            uwall1d(j)=bc_uwall(i,j,l)
          enddo
          call mpiysum(uwall1d,g_uwall1d)
          call trancate(g_uwall1d)
          do j=1, nylocal, 1
            bc_uwall(i,j,l)=g_uwall1d((iprocy-1)*nylocal+j)
          enddo
        enddo
      enddo 
c
      return
      end
c--------------
c 
c----------------
      subroutine mpiysum(vy,g_vy)
      use mesh
      use mesh_common
#ifdef PARALLEL
      use mpistuff
#endif
c
      implicit none
      double precision:: vy(1:nylocal)
      double precision:: g_vy(1:ny), g_t(1:ny)
c
      g_t = 0.d0
c
#ifdef PARALLEL
      g_t((iprocy-1)*nylocal+1:iprocy*nylocal)=vy(1:nylocal)
      if(yprocs.gt.1) then
        call MPI_Allreduce(g_t, g_vy, ny,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm1d, ierr)
        call ErrorHandler(ierr,ERROR_ALLREDUCE)
      else
        g_vy = vy 
      endif
#else
      g_vy = vy
#endif   
c
      return
      end
c-------------------------------
c
c----------------------------------------
      subroutine trancate(vp)
      use mesh
      use mesh_common
c
      implicit none
      double precision:: vp(ny)
      double complex::   v_p(ny) ! physical
      double complex::   v_w(ny) ! wavenumber
      double complex::   v_wn(ny)
      double complex:: v_pn(ny)
      integer:: j, idir
c
      v_p  = 0.d0
      v_w  = 0.d0
      v_wn = 0.d0
      v_pn = 0.d0
c
      v_p(:) = vp(:)
      idir = 1
      call FFT1D(v_p,v_w,ny,idir)     
c   trancate to make sure 2/3 rule    
      call trancatep1(v_w)         
      idir = -1
      call FFT1D(v_w,v_p,ny,idir)
      vp(:) = v_p(:)
c
      return
      end
c--------------------------------
c
c--------------------------------
      subroutine trancatep1(uf)
      use mesh
      use mesh_common
c
      implicit none
      double complex::   uf(ny) 
      double precision:: pai, xi, alpha, vfilter
      integer:: j, p1
c
      pai   = 0.d0
      xi    = 0.d0
      alpha = 0.d0
      p1   = 0
      j    = 0
c
      p1 = 36
      pai = 4.d0*datan(1.d0)
      alpha = 0.6931d0*1.5d0**p1
c
      do j=1, ny, 1
        if(j.eq.1) xi = 0.d0
        if(j.le.ny/2) then         
          xi = 2.d0*(j-1)/ny
        else
          xi = 2.d0*(ny+1-j)/ny
        endif
        vfilter = dexp(-alpha*xi**p1)
        if(j.eq.ny/2+1) vfilter = 0.d0
c  filter
        uf(j) = uf(j)*vfilter
      enddo
c
      return
      end

