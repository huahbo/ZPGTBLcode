c----------------------------------------------------------------------
c     Wan Cheng 
c     August 2012
c----------------------------------------------------------------------
c
c
c----------------------------------------------------------------------
c     set the virtual wall boundary condition
c----------------------------------------------------------------------
      subroutine setinflow(ux,timestep)
      use mesh
      use mesh_common
      use properties
      use InflowVelBC
      use stats
      use virtualwallBC
      use recycle        
c
      implicit none 
      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
c
      double precision:: uinletn(nghost+1,iylo:iyhi,izlo:izhi,3)
      double precision:: uinlett(iylo:iyhi,izlo:izhi,3)
      double precision:: uinlettm1(iylo:iyhi,izlo:izhi,3)
      double precision:: utau_inlet(1:nylocal)
      integer:: irecyclem1, irelocalm1, irecycle, irelocal, iprocxre
      integer:: idata , timestep
      integer:: i, j, k
      double precision:: nscalet
c  1.
      do i=0, nghost, 1
        irecycle = g_ire-nghost+i
        iprocxre = ip_re
        irelocal = i_re-nghost+i
c        if(timestep.lt.1000) then
c        nscalet = nscale
c        else
c        idata   = ioutput
c        nscalet = dlog(g_utau(idata)/g_utau(irecycle))/
c     &          dlog(g_del99(idata)/g_del99(irecycle))
c     &            dlog(g_theta(idata)/g_theta(irecycle))
c        endif
c
c  2. rescaling
        uinlett = 0.d0
        call rescale_inflow(nghost-i,irecycle,iprocxre,irelocal,
     &                      ux,uinlett)
c  3. rescale, recycle  utau
        call setutauin(i,uinlett,irecycle,iprocxre,irelocal)
        uinletn(nghost+1-i,:,:,:) = uinlett(:,:,:)
      enddo

c     consider the y z velo component for inflow location, 
c     this part is not used for multi-recycling
      irecyclem1 = g_ire -1
      irelocalm1 = i_re -1
      uinlettm1 = 0.d0
      call rescale_inflow(nghost,irecyclem1,iprocxre,irelocalm1,
     &                    ux,uinlettm1)
      uinlett(:,:,2:3)=0.5d0*(uinlett(:,:,2:3)+uinlettm1(:,:,2:3))
      uinletn(1,:,:,2:3) = uinlett(:,:,2:3)



c  3. recycling
#ifdef PARALLEL
c      call exchangeinflow5_mirror(iprocxre,uinletn)
      call exchangeinflow5(iprocxre,uinletn)
#else
      do j=iylo, iyhi, 1
        do k=izlo, izhi, 1
          Uin(1,j,k,1) =  uinlett(nylocal+1-j,k,1)
          Uin(1,j,k,2) = -uinlett(nylocal+1-j,k,2)
          Uin(1,j,k,3) =  uinlett(nylocal+1-j,k,3)
        enddo 
      enddo   

c      do j=iylo, iyhi, 1
c        do k=izlo, izhi, 1
c          Uin(1,j,k,1) =  uinlett(j,k,1)
c          Uin(1,j,k,2) =  uinlett(j,k,2)
c          Uin(1,j,k,3) =  uinlett(j,k,3)
c        enddo 
c      enddo         
#endif          
c
      if(iprocx.eq.1 .and. iprocz.eq.1) then
c       bc_uwall
        do i= 1, -nghost+1, -1
          bc_uwall(i,iylo:iyhi,1)=Uin(2-i,iylo:iyhi,0,1)
        enddo
c       bc_utau
        bc_utau(1,1:nylocal)= Uin(1,1:nylocal,-1,1)
        bc_dudz(1,1:nylocal)= bc_utau(1,1:nylocal)**2.d0*Reynolds
       endif
c        


      return
      end
c
c----------------------------------------------------------------------
c     rescale inflow 
c----------------------------------------------------------------------
      subroutine rescale_inflow(ighost,irecycle,iprocxre,irelocal,
     &                          ux,uinlett)
      use mesh
      use mesh_common
      use properties
      use recycle
      use stats
      use InflowVelBC
      use virtualwallBC
      use iounits
      use time
#ifdef PARALLEL
      use mpistuff
#endif
c
      implicit none
      integer:: irecycle,iprocxre, irelocal
      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision:: uinlett(iylo:iyhi,izlo:izhi,3)
c     original velocity in the recycling plane
c       mean
      double precision::  g_um(0:nz),  g_wm(0:nz),  g_zc(0:nz)
      double precision:: gt_um(0:nz), gt_wm(0:nz), gt_zc(0:nz)
      double precision:: uf_re(izlo:izhi,nvar)
      double precision:: gt_uf(0:nz,nvar)
      double precision:: gt_uf2d(iylo:iyhi,0:nz,nvar)
      double precision::  g_uf2d(iylo:iyhi,0:nz,nvar)
c     rescaled fluc velocity in the recycling plane
      double precision:: uf_out(iylo:iyhi,izlo:izhi,nvar) ! ux fluctuation
      double precision:: uf_inn(iylo:iyhi,izlo:izhi,nvar)
      double precision:: uft_out(izlo:izhi,nvar)
      double precision:: uft_inn(izlo:izhi,nvar)
c     rescaled mean velocity in the recycling plane
      double precision:: um_out(izlo:izhi)
      double precision:: um_inn(izlo:izhi)
      double precision:: wm_out(izlo:izhi)
      double precision:: wm_inn(izlo:izhi)
c       
      integer::  utsize
c
      double precision:: delta_re, del99_re, theta_re
      double precision:: delta_in, del99_in, theta_in
      double precision:: utau_re, utau_in
      double precision:: eta,heavyside,weight !,get_heavyside
      double precision:: get_weight_dis, get_weight, get_weight_theta
      integer:: i,j,k, ire, iin, iytlo,iztlo,iythi,izthi,ighost
c
      double precision:: zvxy, zvz, get_heavyside
c
      um_inn = 0.d0
      um_out = 0.d0
      wm_inn = 0.d0
      wm_out = 0.d0

c     g_zc
      gt_zc = 0.d0
      if(iprocx.eq.iprocxre .and. iprocy.eq.1) then
        gt_zc((iprocz-1)*nzlocal+1:iprocz*nzlocal) = zc(1:nzlocal) 
      endif
#ifdef PARALLEL
      g_zc = 0.d0
      call MPI_Allreduce(gt_zc, g_zc, nz+1,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
      call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
      g_zc = gt_zc
#endif  
c
c 0. get mean velocity
      gt_um = 0.d0
      gt_wm = 0.d0
      if(iprocx.eq.iprocxre .and. iprocy.eq.1) then
        if(iprocz.eq.1) gt_um(0) = g_vwall(irecycle)
        gt_um((iprocz-1)*nzlocal+1:iprocz*nzlocal) = 
     &     st_umean(irelocal,1:nzlocal,1)
        gt_wm((iprocz-1)*nzlocal+1:iprocz*nzlocal) = 
     &     st_umean(irelocal,1:nzlocal,3)
      endif
#ifdef PARALLEL
      g_um = 0.d0
      call MPI_Allreduce(gt_um, g_um, nz+1,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
      call ErrorHandler(ierr,ERROR_ALLREDUCE)
      g_wm = 0.d0
      call MPI_Allreduce(gt_wm, g_wm, nz+1,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
      call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
      g_um = gt_um
      g_wm = gt_wm
#endif    
c
c  01. get destination inflow parameters from global statistics
      if(iprocx.eq.iprocxre) then
        iin = 1
c        del99_in = g_del99(iin)
        del99_in = 1.d0-(g_del99(2)-g_del99(1))*ighost  !destination del99
        delta_in = 1.d0/8.d0 !g_delta(iin)
        theta_in = g_theta(1)- (g_theta(2)-g_theta(1))*ighost
c  02. calculate recycle information with inlet information
        ire = irelocal
        del99_re = g_del99(irecycle)
        delta_re = g_delta(irecycle)
        utau_re  = g_utau(irecycle)
        theta_re = g_theta(irecycle)
        utau_in  = utau_re *(theta_re/theta_in)**nscale
c  03. rescale mean velocities
        call rescale_mean(del99_in,del99_re,utau_in,utau_re,g_zc,
     &                   g_um,g_wm,um_inn,um_out,wm_inn,wm_out)
c        call rescale_mean(delta_in,delta_re,utau_in,utau_re,g_zc,
c     &                   g_um,g_wm,um_inn,um_out,wm_inn,wm_out)
c        call rescale_mean(theta_in,theta_re,utau_in,utau_re,g_zc,
c     &                   g_um,g_wm,um_inn,um_out,wm_inn,wm_out)

      endif
c 
c 1.  get flucutation velocity
      g_uf2d = 0.d0
      gt_uf2d(iylo:iyhi,0:nz,1:3) = 0.d0   
      if(iprocx.eq.iprocxre)  then
        do j = iylo, iyhi, 1
          do k= 1, nzlocal, 1
            uf_re(k,1) = ux(irelocal,j,k,1) - st_umean(irelocal,k,1)
            uf_re(k,2) = ux(irelocal,j,k,2) !- st_umean(irelocal,j,2)
            uf_re(k,3) = ux(irelocal,j,k,3) - st_umean(irelocal,k,3)
            gt_uf2d(j,(iprocz-1)*nzlocal+k,1:3) = uf_re(k,1:3)
          enddo
          if(iprocz.eq.1) 
     &      gt_uf2d(j,0,1) = bc_uwall(irelocal,j,1)-g_vwall(irecycle)
        enddo
      endif
#ifdef PARALLEL     
        call MPI_BARRIER(comm3d,ierr)   
        if(iprocx.eq.iprocxre) then 
          utsize = (iyhi-iylo+1)*(nz+1)*nvar
          call MPI_Allreduce(gt_uf2d, g_uf2d, utsize,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm1dre, ierr)
          call ErrorHandler(ierr,ERROR_ALLREDUCE)
        endif
#else
      g_uf2d= gt_uf2d
#endif

c  11.rescale velocity fluctuation.
      if(iprocx.eq.iprocxre) then
        do j=iylo,iyhi, 1
          gt_uf(0:nz,1:3) = g_uf2d(j,0:nz,1:3) 
          call rescale_fluc(del99_in,del99_re,
     &                    utau_in,utau_re,g_zc,gt_uf,uft_inn,uft_out)
c          call rescale_fluc(delta_in,delta_re,
c     &                    utau_in,utau_re,g_zc,gt_uf,uft_inn,uft_out)
c          call rescale_fluc(theta_in,theta_re,
c     &                    utau_in,utau_re,g_zc,gt_uf,uft_inn,uft_out)

          uf_inn(j,izlo:izhi,1:3)= uft_inn(izlo:izhi,1:3)
          uf_out(j,izlo:izhi,1:3)= uft_out(izlo:izhi,1:3)  
        enddo  
      endif      
c 
c     points near upper wall
      if(iprocx.eq.iprocxre .and. iprocz.eq.zprocs) then
        wm_inn(nzlocal+1) = wm_inn(nzlocal)
        wm_out(nzlocal+1) = wm_out(nzlocal)
        uf_inn(iylo:iyhi,nzlocal+1,1:3)=uf_inn(iylo:iyhi,nzlocal,1:3)
        uf_out(iylo:iyhi,nzlocal+1,1:3)=uf_out(iylo:iyhi,nzlocal,1:3)
      endif

c
c  2. get the recycling uin.
      uinlett = 0.d0
      if(iprocx.eq.iprocxre) then
        do j=iylo, iyhi, 1
          do k=izlo, izhi, 1
            if(zc(k).gt. bcupper ) then
              um_inn(k)=1.d0
              um_out(k)=1.d0
              uf_inn(j,k,1) = 0.d0
              uf_out(j,k,1) = 0.d0
              uf_inn(j,k,2) = 0.d0
              uf_out(j,k,2) = 0.d0
              uf_inn(j,k,3) = 0.d0
              uf_out(j,k,3) = 0.d0
            endif
              eta = (zc(k)+bc_h0dz*dz ) / del99_in
              weight = get_weight(eta)
c             eta = (zc(k)+bc_h0dz*dz ) / delta_in
c             weight = get_weight_dis(eta)
c              eta = (zc(k)+bc_h0dz*dz ) / theta_in
c              weight = get_weight_theta(eta)

              zvxy = zc(k)+ bc_h0dz*dz
              heavyside = get_heavyside(zvxy,1.2d0,1.3d0)
            uinlett(j,k,1)=
     &              (um_inn(k)+uf_inn(j,k,1)*heavyside)*(1.d0-weight)
     &             +(um_out(k)+uf_out(j,k,1)*heavyside)*weight   

            uinlett(j,k,2)=
     &              (         +uf_inn(j,k,2)*heavyside)*(1.d0-weight)
     &             +(         +uf_out(j,k,2)*heavyside)*weight

              eta = (zc(k)+bc_h0dz*dz -0.5d0*dz ) / del99_in
              weight = get_weight(eta)
c              eta = (zc(k)+bc_h0dz*dz -0.5d0*dz ) / delta_in
c              weight = get_weight_dis(eta)
c              eta = (zc(k)+bc_h0dz*dz -0.5d0*dz ) / theta_in
c              weight = get_weight_theta(eta)

              zvz = zc(k)+ bc_h0dz*dz-0.5d0*dz
              heavyside = get_heavyside(zvz,1.2d0,1.3d0)
            uinlett(j,k,3)= 
     &              (wm_inn(k)+uf_inn(j,k,3)*heavyside)*(1.d0-weight)
     &             +(wm_out(k)+uf_out(j,k,3)*heavyside)*weight
           enddo
         enddo
       endif
c
       return
c
c--- check part
c       write(*,*) 'MAX/MIN rescale u',
c     &              maxval(uinlett(1:nylocal,1:nzlocal,1)),
c     &              minval(uinlett(1:nylocal,1:nzlocal,1))
c
       end
c
c----------------------------------------------------------------------
c     rescale the mean velocity 
c----------------------------------------------------------------------
      subroutine rescale_mean(del_in,del_re,utau_in,utau_re,g_zc,
     &            g_um,g_wm,um_inn,um_out,wm_inn,wm_out)
      use mesh
      use time
      use mesh_common
      use properties
      use stats
      use virtualwallBC
      implicit none
      double precision:: g_um(0:nz),g_wm(0:nz),g_zc(0:nz)
      double precision:: um_out(izlo:izhi)
      double precision:: um_inn(izlo:izhi)
      double precision:: wm_inn(izlo:izhi)
      double precision:: wm_out(izlo:izhi)
      double precision:: del_in, del_re, utau_in, utau_re
c
      double precision:: gt_in(0:nz)
      double precision:: nu, normal_in, normal_re
      double precision:: z_modify, xi 
      integer::  k
c
      gt_in = 0.d0
c
      nu = 1.d0/Reynolds    
      xi = utau_in/utau_re
c  u mean
c    inn
      z_modify  = bc_h0dz*dz ! z in virtual wall
      normal_in = nu/utau_in
      normal_re = nu/utau_re  ! normal for rescale
      call nz_scale_u(normal_in,normal_re,g_zc,g_um,gt_in,z_modify)
      gt_in(:)  = gt_in(:)  * xi
      call set_uin_local(gt_in, um_inn)
c    outer 
      z_modify  = bc_h0dz*dz
      normal_in = del_in
      normal_re = del_re
      call nz_scale_u(normal_in,normal_re,g_zc,g_um,gt_in,z_modify)
      gt_in(:) = gt_in(:) * xi + (1.d0-xi)* velRef
      call set_uin_local(gt_in, um_out)
c
      do k=izlo, izhi, 1
        if(zc(k).gt. bcupper )  then
          um_inn(k)=1.d0
          um_out(k)=1.d0
        endif
      enddo
c  
c  w mean
c    inner
      z_modify = -0.5d0*dz + bc_h0dz*dz
      normal_in= nu/utau_in
      normal_re= nu/utau_re  
      call nz_scale_vw(normal_in,normal_re,g_zc,g_wm,gt_in,z_modify)
      call set_uin_local(gt_in,wm_inn)
c    outer 
      z_modify = -0.5d0*dz + bc_h0dz*dz
      normal_in= del_in
      normal_re= del_re
      call nz_scale_vw(normal_in,normal_re,g_zc,g_wm,gt_in,z_modify)
      call set_uin_local(gt_in,wm_out) 
c
      return
c
c ----check mean
c     do k=1, nzlocal,1
c       write(*,*)'mean :',k,um_inn(k), um_out(k), wm_inn(k),wm_out(k),
c     &           iprocx, iprocy, iprocz
c     enddo
c
      end
c
c----------------------------------------------------------------------
c      rescale fluctuation
c----------------------------------------------------------------------
      subroutine rescale_fluc(del_in,del_re,utau_in,utau_re,
     &            g_zc,uf,uf_inn,uf_out)
      use mesh
      use mesh_common
      use properties
      use stats
      use virtualwallBC
c
      implicit none
      double precision:: uf(0:nz,nvar), g_zc(0:nz)
      double precision:: uf_out(izlo:izhi,nvar)
      double precision:: uf_inn(izlo:izhi,nvar)
      double precision:: del_re, del_in, utau_in, utau_re
c
      double precision:: uf_t(0:nz), uf_inl(0:nz)
      double precision:: nu, normal
      double precision:: z_in, z_modify, z_max, z_min
      double precision:: normal_in, normal_re
      double precision:: xi,  vel_inn,vel_out
      integer:: k 
c
      nu = 1.d0/Reynolds
      xi = utau_in/utau_re
c u fluc      
c    inn  
      z_modify = bc_h0dz*dz
      normal_in = nu/utau_in
      normal_re = nu/utau_re
      call nz_scale_u(normal_in,normal_re,g_zc,uf(:,1),uf_inl,z_modify)
      call set_uin_local(uf_inl,uf_inn(:,1))
      uf_inn(:,1) = uf_inn(:,1) * xi
c    out
      normal_in = del_in
      normal_re = del_re
      call nz_scale_u(normal_in,normal_re,g_zc,uf(:,1),uf_inl,z_modify)
      call set_uin_local(uf_inl,uf_out(:,1))
      uf_out(:,1) = uf_out(:,1) * xi  
c--------
c v fluc
c   inn
      z_modify = bc_h0dz*dz
      normal_in = nu/utau_in
      normal_re = nu/utau_re       
      call nz_scale_vw(normal_in,normal_re,g_zc,uf(:,2),uf_inl,z_modify)
      call set_uin_local(uf_inl,uf_inn(:,2))
      uf_inn(:,2) = uf_inn(:,2) * xi
c   outer 
      normal_in = del_in
      normal_re = del_re
      call nz_scale_vw(normal_in,normal_re,g_zc,uf(:,2),uf_inl,z_modify)
      call set_uin_local(uf_inl,uf_out(:,2))
      uf_out(:,2) = uf_out(:,2) * xi  
c---------
c w fluc
c   inn
      z_modify = -0.5d0*dz + bc_h0dz*dz
      normal_in = nu/utau_in
      normal_re = nu/utau_re  
      call nz_scale_vw(normal_in,normal_re,g_zc,uf(:,3),uf_inl,z_modify)
      call set_uin_local(uf_inl,uf_inn(:,3))
      uf_inn(:,3) = uf_inn(:,3) * xi
c   outer
      normal_in = del_in
      normal_re = del_re
      call nz_scale_vw(normal_in,normal_re,g_zc,uf(:,3),uf_inl,z_modify)
      call set_uin_local(uf_inl,uf_out(:,3))
      uf_out(:,3) = uf_out(:,3) * xi  
c
      return
c--------
c      check fluc
c      do k = 1, nzlocal, 1
c       write(*,*)'fluc :', k, uf_out(k,1),uf_inn(k,1)
c      enddo
c
      end
c
c----------------------------------------------------------------------
c     scale module
c----------------------------------------------------------------------
      subroutine nz_scale_u(nor_in,normal,g_zc,uf_t,gt_in,z_modify_in)
      use mesh
      use mesh_common
c
      implicit none
      double precision:: g_zc(0:nz), uf_t(0:nz), gt_in(0:nz)
      double precision:: nor_in, normal
c
      double precision:: z_in, z_max, z_min, z_modify, z_modify_in
      double precision:: z_max_near, z_min_near
      integer:: k
c
      do k=0, nz, 1
c        if(k.eq.0) then
          z_modify = z_modify_in
c        else 
c          z_modify = 0.d0
c        endif
        z_in =   (g_zc(k)+z_modify)/nor_in
        z_max =  (g_zc(nz)+z_modify)/normal
        z_min =  (g_zc(0)+z_modify)/normal
        if( z_in .lt. z_max .and. z_in .gt. z_min ) then
          call linear_interpolation
     &           (g_zc,uf_t,z_in,gt_in(k),normal,z_modify)
        endif
        if(z_in .ge. z_max) then
          z_max_near = (g_zc(nz-1)+z_modify)/normal
          call linear_inter_bc
     &           (uf_t(nz-1),uf_t(nz),z_max_near,z_max,z_in,gt_in(k))
        endif
        if (z_in .le. z_min ) then
          z_min_near = (g_zc(1)+z_modify)/normal 
          call linear_inter_bc
     &           (uf_t(1),uf_t(0),z_min_near,z_min,z_in,gt_in(k))
        endif
      enddo
c
      return
      end
c----------------------------------------------------------------------
c     scale  for w
c----------------------------------------------------------------------
      subroutine nz_scale_vw(nor_in,normal,g_zc,uf_t,gt_in,z_modify_in)
      use mesh
      use mesh_common
c
      implicit none
      double precision:: g_zc(0:nz), uf_t(0:nz), gt_in(0:nz)
      double precision:: nor_in, normal
c
      double precision:: z_in, z_max, z_min, z_modify, z_modify_in
      double precision:: z_max_near, z_min_near
      integer:: k
c
      do k=1, nz, 1
c        if(k.eq.1) then
          z_modify = z_modify_in
c        else
c          z_modify =  -0.5d0*dz
c        endif
        z_in = (g_zc(k)+z_modify)/nor_in
        z_max =  (g_zc(nz)+z_modify)/normal
        z_min =  (g_zc(1)+z_modify)/normal
        if( z_in .lt. z_max .and. z_in .gt. z_min ) then
          call linear_interpolation
     &           (g_zc,uf_t,z_in,gt_in(k),normal,z_modify)
        endif
        if(z_in .ge. z_max) then
          z_max_near = (g_zc(nz-1)+z_modify)/normal
          call linear_inter_bc
     &           (uf_t(nz-1),uf_t(nz),z_max_near,z_max,z_in,gt_in(k))
        endif
        if (z_in .le. z_min ) then
          z_min_near = (g_zc(2)+z_modify)/normal 
          call linear_inter_bc
     &           (uf_t(2),uf_t(1),z_min_near,z_min,z_in,gt_in(k))
        endif
      enddo
c
      return
      end
c


c
c----------------------------------------------------------------------
c     linear interpolation  part
c----------------------------------------------------------------------
      subroutine linear_interpolation
     &            (g_zc,uf,z_aim,vel_aim, normal,z_modify)
      use mesh
      use mesh_common
      use properties
      use stats
      implicit none       
      double precision:: uf(0:nz)
      double precision:: g_zc(0:nz)
      double precision:: z_aim, vel_aim, normal, z_re, z0, z1, alpha
      double precision:: z_modify
      integer:: k, inde
c
      inde = 0
      do k=0, nz, 1
        z_re = (g_zc(k)+z_modify) / normal
        if(z_re .lt. z_aim) inde = inde+1
      enddo
      z0 = (g_zc(inde-1)+z_modify)/normal
      z1 = (g_zc(inde)+z_modify)/normal
      alpha = (z_aim-z0)/(z1-z0)
      vel_aim = (1.d0-alpha)*uf(inde-1)+alpha*uf(inde)
c
      return
      end
c
c----------------------------------------------------------------------
c     linear interpolation  part boundary
c----------------------------------------------------------------------
      subroutine linear_inter_bc
     &            (u1,u2,z1,z2,z_aim,v_aim)
      implicit none       
      double precision:: u1, u2, z1, z2, z_aim, v_aim
      double precision:: alpha
c     interpolation direction: from 1 to 2
c     aim > 2 > 1, or, 1 > 2 > aim
      alpha = (z_aim-z1)/(z2-z1)
      v_aim = alpha*(u2-u1)+u1
c    
      return
      end
c
c----------------------------------------------------------------------
c     from global to local
c----------------------------------------------------------------------
      subroutine set_uin_local(gt, v_aim)
      use mesh
      use mesh_common
c
      implicit none
      double precision:: gt(0:nz)
      double precision:: v_aim(izlo:izhi)
c
      if(zprocs.eq.1) then
        v_aim(0:nzlocal) = gt(0:nzlocal)
      else
        if(iprocz.eq.1) v_aim(0:izhi) = gt(0:izhi)
        if(iprocz.eq.zprocs) v_aim(izlo:nzlocal) = 
     &            gt( (iprocz-1)*nzlocal+izlo:nz)
        if(iprocz.gt.1 .and. iprocz.lt. zprocs) v_aim(izlo:izhi) =
     &       gt( (iprocz-1)*nzlocal+izlo:(iprocz-1)*nzlocal+izhi)  
      endif       
c
      return
      end
c
c----------------------------------------------------------------------
c     weight function
c----------------------------------------------------------------------
      double precision function get_weight(eta)
      implicit none
      double precision eta, a, b
      a = 4.d0
      b = 0.2d0
      get_weight = 1.d0/2.d0
     &    *( 1.d0+dtanh( a*(eta-b)/( (1.d0-2.d0*b)*eta+b ) )/dtanh(a))
      if(eta .gt. 1.d0) get_weight = 1.d0
c
      return
      end
c
c----------------------------------------------------------------------
c     weight function for displacement thickness
c----------------------------------------------------------------------
      double precision function get_weight_dis(eta0)
      implicit none
      double precision eta0, eta, a, b
      a = 4.d0
      b = 0.2d0
      eta = eta0/8.d0 ! approximate relation between measures of BL
      get_weight_dis = 1.d0/2.d0
     &    *( 1.d0+dtanh( a*(eta-b)/( (1.d0-2.d0*b)*eta+b ) )/dtanh(a))
      if(eta .gt. 1.d0) get_weight_dis = 1.d0
c
      return
      end
c
c----------------------------------------------------------------------
c     weight function for momentum thickness
c----------------------------------------------------------------------
      double precision function get_weight_theta(eta0)
      implicit none
      double precision eta0, eta, a, b
      a = 4.d0
      b = 0.2d0
      eta = eta0/10.d0 ! approximate relation between measures of BL
      get_weight_theta = 1.d0/2.d0
     &    *( 1.d0+dtanh( a*(eta-b)/( (1.d0-2.d0*b)*eta+b ) )/dtanh(a))
      if(eta .gt. 1.d0) get_weight_theta = 1.d0
c
      return
      end
c
c----------------------------------
c     heavyside function for combining the inner region and outer region
c     used by Inoue, now not used in our code
c-----------------------------------------
      double precision function get_heavyside(z, z100, z0)
      implicit none
      double precision delta99_in, e, phi, H, z100, z0, z, pai
c
      pai = 4.d0*datan(1.d0)
      delta99_in = 1.d0
      e = (z0-z100)*delta99_in/2.d0
      phi = z - z100*delta99_in - e
      if(phi.lt. -e) then
        H = 1.d0
      else if(phi .gt. e) then
        H = 0.d0
      else
        H = 0.5d0*(1.d0-phi/e-1.d0/pai*dsin(pai*phi/e))
      endif
      get_heavyside = H
c
      return
      end
c
c----------------------------------------------------------------------
c     set utau at inlet                                         
c----------------------------------------------------------------------
      subroutine setutauin_old(irecycle,iprocxre,irelocal)
      use mesh
      use mesh_common
      use properties
      use stats
      use virtualwallBC
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      integer:: irecycle, iprocxre, irelocal
c
      double precision:: utaure(1:ny)
      double precision:: g_utauin(1:ny)
      double precision:: utaurelocal(1:nylocal)
      double precision:: theta_in, theta_re, utre
      integer:: j
c
#ifdef PARALLEL
      utaure = 0.d0
      g_utauin = 0.d0
      if(iprocx.eq.iprocxre .and. iprocz.eq.1) then
        utaure( (iprocy-1)*nylocal+1 : iprocy*nylocal ) = 
     &   bc_utau(irelocal,1:nylocal)
      endif
      call MPI_Allreduce(utaure, g_utauin, ny,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
      call ErrorHandler(ierr,ERROR_ALLREDUCE)
c      do j=1, ny
c       g_utauin(j) = g_utau(irecycle)
c      enddo
#else
      utaure(1:ny) = bc_utau(irelocal,1:nylocal)
      utre = g_utau(irecycle)  ! statistical
      g_utauin = utaure  ! utre  
#endif
      if(iprocx.eq.1) then
        theta_in =  0.1d0  ! g_theta(1)   cccccccccccccccccccccccccccccccccccccccccc
        theta_re = g_theta(irecycle)
        do j=1,ny
          g_utauin(j) = g_utauin(j) *(theta_re/theta_in)**
     &           (1.d0/2.d0/(7.d0-1.d0)) 
        enddo
#ifdef PARALLEL
        bc_utau(1,1:nylocal)=
     &  g_utauin( (iprocy-1)*nylocal+1 : iprocy*nylocal )
#else
        bc_utau(1,1:nylocal)=g_utauin(1:nylocal)
#endif
      endif
c
      do j=1, nylocal, 1
        bc_dudz(1,j) = bc_utau(1,j)*bc_utau(1,j)*Reynolds
      enddo

c
      return
      end
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c     set utau at inlet                                         
c----------------------------------------------------------------------
      subroutine setutauin(ighost,uinlett,irecycle,iprocxre,irelocal)
      use mesh
      use mesh_common
      use properties
      use recycle
      use stats
      use virtualwallBC
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      integer:: irecycle, iprocxre, irelocal
c
      double precision:: uinlett(iylo:iyhi,izlo:izhi,3)
      double precision:: utaure(1:ny)
      double precision:: g_utauin(1:ny)
      double precision:: utaurelocal(1:nylocal)
      double precision:: theta_in, theta_re, utre
      integer:: j, k, ighost
c
      if(iprocz.eq.1 .and. iprocx.eq. iprocxre) then
c        theta_in = 0.1d0 - (g_theta(2)-g_theta(1))*ighost !!!!!!!!
        theta_in = g_theta(1)-(g_theta(2)-g_theta(1))*ighost
        theta_re = g_theta(irecycle)
        k = -1
        do j=1, nylocal, 1
          uinlett(j,k,1) = bc_utau(irelocal,j)*(theta_re/theta_in)**
     &           nscale
        enddo
      endif
c
      return
      end
c
c----------------------------------------------------------------------

