c----------------------------------------------------------------------
c     Wan Cheng 
c     August 2012
c----------------------------------------------------------------------
c
c
c----------------------------------------------------------------------
c     calculate the sgs tensor Tij
c----------------------------------------------------------------------
      subroutine sgsmodel(uu,Tij,Tij_cn,timestep)
      use mesh
      use mesh_common
      use properties
      use virtualwallBC
      use stats
#ifdef PARALLEL
      use mpistuff
#endif
c
      implicit none
      double precision:: uu(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar) 
      double precision:: Tij(-1:nxlocal+2,-1:nylocal+2,-1:nzlocal+2,6)
      double precision:: Tij_cn(0:nxlocal+2,0:nylocal+2,0:nzlocal+2,6)
c
      double precision:: ucn(-1:nxlocal+3,-1:nylocal+3,-1:nzlocal+3,3)
      double precision:: dudx(0:nxlocal+2,0:nylocal+2,0:nzlocal+2,3,3)
      double precision:: K_cn(0:nxlocal+2,0:nylocal+2,0:nzlocal+2)
c
      double precision:: u_near(27),v_near(27),w_near(27)
      double precision:: x_near(27),y_near(27),z_near(27)
      double precision:: dudxtmp(3,3)
      double precision:: e(3)
      double precision:: Txx, Tyy, Tzz, Txy, Tyz, Tzx
      double precision:: K_tmp, del, les_del
c     
      integer:: yes_near(27)
      integer:: ix_near, size_near, nx_near, ny_near, nz_near, ix_this
      integer:: ip, iq, ir, i_near, j_near, k_near
      integer:: i, j, k, timestep, in1 , klo , l !, itmp
      integer:: iplocal, ilocal
c
      double precision:: aa, bb
c
      aa= 9.D0/8.D0
      bb=-1.D0/8.D0
c
      ucn = 1.d7
      dudx= 1.d7
      Tij_cn = 1.d7
      Tij    = 1.d7

c      interpolation to corner
      if(iprocz.ne.1) then
        call interpolate_ududx_cn(uu,dudx,ucn)
      else
        call interpolate_ududx_cn_wall(uu,dudx,ucn)
      endif
c     cutoff length , del
      del = dz !  (dx*dy*dz)**(1.d0/3.d0)
      les_del =  (dx*dy*dz)**(1.d0/3.d0)
c
      size_near = 27
      nx_near = 1
      ny_near = 1
      nz_near = 1
      ix_this = 100
      K_tmp = 0.d0
c
c  part 1  getting Tij_cn: tensor at corner
c     Getting Tij from spiral SGS model
      klo = 0
      if(iprocz.eq.1) klo = 1
      do i=0,nxlocal+2,1
        do j=0,nylocal+2,1
          do k=klo,nzlocal+2,1
          ix_near = 1
          if(zc(k).gt. bcupper ) then
            Tij_cn(i,j,k,1) = 0.d0
            Tij_cn(i,j,k,2) = 0.d0
            Tij_cn(i,j,k,3) = 0.d0
            Tij_cn(i,j,k,4) = 0.d0
            Tij_cn(i,j,k,5) = 0.d0
            Tij_cn(i,j,k,6) = 0.d0
            K_cn(i,j,k) = 0.d0         
          else   !    sgs model
            do ip = - nx_near, nx_near, 1
              do iq = -ny_near, ny_near, 1
                do ir = -nz_near, nz_near, 1 
                i_near = i+ ip
                j_near = j+ iq
                k_near = k+ ir
                if( (iprocx.eq.1     .and.i.eq.1.and.i_near.lt.1).or.
     &              (iprocx.eq.xprocs.and.i.eq.nxlocal+1 
     &         .and. i_near.gt.nxlocal+1))  i_near = i-ip
                if( (iprocz.eq.1     .and.k.eq.1.and.k_near.lt.1).or. 
     &              (iprocz.eq.zprocs.and.k.eq.nzlocal+1 
     &          .and. k_near.gt.nzlocal+1))  k_near = k-ir
                x_near(ix_near)= xc(i_near)-dx/2.d0
                y_near(ix_near)= yc(j_near)-dy/2.d0  !!!!!!!!!!!!!!!!!!!!!!!!
                z_near(ix_near)= zc(k_near)-dz/2.d0
                u_near(ix_near)= ucn(i_near,j_near,k_near,1)
                v_near(ix_near)= ucn(i_near,j_near,k_near,2)
                w_near(ix_near)= ucn(i_near,j_near,k_near,3)
                if(ip.eq.0 .and. iq.eq.0 .and. ir.eq.0)ix_this= ix_near
                yes_near(ix_near) = 1
                if(ix_near.eq.ix_this) yes_near(ix_near)=0
                ix_near = ix_near +1
                enddo
              enddo
            enddo
            dudxtmp(1:3,1:3) = dudx(i,j,k,1:3,1:3)
c       
            e(1) = 0.d0
            e(2) = 0.d0
            e(3) = 0.d0
            call spiral_sgs_stress_v3
     &       (u_near, v_near, w_near, x_near, y_near, z_near, 
     &       yes_near, size_near, ix_this, dudxtmp, e, les_del, K_tmp,
     &       Txx, Tyy, Tzz, Txy, Tyz, Tzx)
            Tij_cn(i,j,k,1) = Txx     !                -- 1  2  3 --   
            Tij_cn(i,j,k,2) = Txy     !         Tij =  |  2  4  5  |
            Tij_cn(i,j,k,3) = Tzx     !                -- 3  5  6 --
            Tij_cn(i,j,k,4) = Tyy
            Tij_cn(i,j,k,5) = Tyz
            Tij_cn(i,j,k,6) = Tzz
            K_cn(i,j,k) = K_tmp
c
            if(iprocz.eq.1) then  ! virtual wall
              if(k.eq.2) then
c              e(1) = 1.d0/dsqrt(2.d0)
c              e(2) = 0.d0
c              e(3) = 1.d0/dsqrt(2.d0)
c              call spiral_sgs_stress_v4
c     &         (u_near, v_near, w_near, x_near, y_near, z_near, 
c     &         yes_near, size_near, ix_this, dudxtmp, e, del, K_tmp,
c     &         Txx, Tyy, Tzz, Txy, Tyz, Tzx) 
              e(1) = 1.d0
              e(2) = 0.d0
              e(3) = 0.d0
              call spiral_sgs_stress_v4
     &         (u_near, v_near, w_near, x_near, y_near, z_near, 
     &         yes_near, size_near, ix_this, dudxtmp, e, del, K_tmp,
     &         Txx, Txx, Txx, Txx, Txx, Txx) 
              if(i.ge.1.and.i.le.nxlocal+1.and.j.ge.1.and.j.le.nylocal)
     &        bc_Kappa(i,j)= 0.45d0/2.d0*dsqrt(K_tmp)/dsqrt(dabs(Tzx))
              endif
            endif    ! virtual wall
          endif          
c
          enddo
        enddo
      enddo
c                   
c Part 2 Get Tij at center
      call interpolateT(Tij_cn, Tij)
c
c Part 3 values for virtualwall
c     Txx at h is the first cell, at the center
      bc_Txx(-1:nxlocal+2,1:nylocal) 
     & = Tij(-1:nxlocal+2,1:nylocal,1,1)
c     Txy in the middle of two corners in z-direction
      k = 1
      do i= 1, nxlocal+1, 1
        do j= 0, nylocal+2, 1
          bc_Txy(i,j)=aa*(Tij_cn(i,j,k  ,2)+Tij_cn(i,j,k+1,2))/2.d0+
     &                bb*(Tij_cn(i,j,k-1,2)+Tij_cn(i,j,k+2,2))/2.d0
        enddo
      enddo
c     Txz in the center of the the u-edge
      k = 1
      do i= 1, nxlocal+1, 1
        do j= 1, nylocal, 1
          bc_Txz(i,j) = 
     &       ( aa*(Tij_cn(i,j  ,k,3)+Tij_cn(i,j+1,k,3))/2.d0+
     &         bb*(Tij_cn(i,j-1,k,3)+Tij_cn(i,j+2,k,3))/2.d0
     &       + aa*(Tij_cn(i,j  ,k+1,3)+Tij_cn(i,j+1,k+1,3))/2.d0+
     &         bb*(Tij_cn(i,j-1,k+1,3)+Tij_cn(i,j+2,k+1,3))/2.d0 )/2.d0
        enddo
      enddo
c
      iplocal = int(ioutput/nxlocal)+1
      ilocal  = mod(ioutput,nxlocal)
      st_Tijct = 0.d0
      do i=1, nxlocal, 1      
        do j=1, nylocal, 1
          do k=1, nzlocal, 1
          st_Tijct(i,j,k,1:6)=Tij(i,j,k,1:6)
          enddo
        enddo
      enddo

c
c      check the min/max
c       write(*,*)maxval(Tij_cn), minval(Tij_cn), maxval(Tij),minval(Tij)
c       stop
c       write(*,*) 'Min/max ',maxval(bc_Txy),minval(bc_Txy),
c     &   maxval(bc_Txz),minval(bc_Txz),maxval(bc_Txx),minval(bc_Txx)
c      
      return
      end
c
c----------------------------------------------------------------------
c     getting u and dudx at CorNer
c----------------------------------------------------------------------
      subroutine interpolate_ududx_cn(uu,dudx,ucn)
      use mesh
      use mesh_common
      use properties
      use virtualwallBC
#ifdef PARALLEL
      use mpistuff
#endif
c
      implicit none
      double precision:: uu(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: ucn(-1:nxlocal+3,-1:nylocal+3,-1:nzlocal+3,3)
      double precision:: dudx(0:nxlocal+2,0:nylocal+2,0:nzlocal+2,3,3)
c
      double precision:: uca(-2:nxlocal+4, -1:nylocal+3, -1:nzlocal+3)
      double precision:: vca(-1:nxlocal+3, -2:nylocal+4, -1:nzlocal+3)
      double precision:: wca(-1:nxlocal+3, -1:nylocal+3, -2:nzlocal+4)
      double precision:: tmpx(-1:nylocal+3,-3:nzlocal+4)
      double precision:: tmpy(-1:nxlocal+3,-3:nzlocal+4)
      double precision:: tmpz(-1:nxlocal+3,-3:nylocal+4)
      double precision:: ux_y( 0:nylocal+2,-2:nzlocal+3)
      double precision:: ux_z(-2:nylocal+3, 0:nzlocal+2)
      double precision:: uy_x( 0:nxlocal+2,-2:nzlocal+3)
      double precision:: uy_z(-2:nxlocal+3, 0:nzlocal+2)
      double precision:: uz_x( 0:nxlocal+2,-2:nylocal+3)
      double precision:: uz_y(-2:nxlocal+3, 0:nylocal+2)
c
       double precision:: bctmp(1:nxlocal+1,1:nylocal,1:2,2)
       double precision:: bcduzt(1:nxlocal+1,1:nylocal,2)
c
      double precision:: aa,bb,cc,dd
      integer:: i,j,k
c
      aa= 9.D0/8.D0
      bb=-1.D0/8.D0
      cc= 8.d0/12.d0
      dd=-1.d0/12.d0
c      velocity gradients
c      (u,x)
       do i=-2,nxlocal+4,1
         do j=-1, nylocal+3, 1
           do k= -3, nzlocal+4, 1
           tmpx(j,k)= aa*(uu(i,j-1,k,1)+uu(i,j  ,k,1))/2.d0+
     &                bb*(uu(i,j-2,k,1)+uu(i,j+1,k,1))/2.d0
           enddo
         enddo
         do j=-1, nylocal+3, 1
           do k= -1, nzlocal+3, 1
           uca(i,j,k)=aa*(tmpx(j,k-1)+tmpx(j,k  ))/2.d0+
     &                bb*(tmpx(j,k-2)+tmpx(j,k+1))/2.d0
           enddo 
         enddo
       enddo

       do i=0, nxlocal+2, 1 
         do j=0, nylocal+2, 1
           do k=0, nzlocal+2, 1
           dudx(i,j,k,1,1)=cc*(uca(i+1,j,k)-uca(i-1,j,k))*dxi+
     &                     dd*(uca(i+2,j,k)-uca(i-2,j,k))*dxi
           enddo
         enddo
       enddo
c       (u,y)
       do i=0,nxlocal+2,1
         do j= -2, nylocal+3, 1
           do k= 0, nzlocal+2, 1
           ux_z(j,k)= aa*(uu(i,j,k-1,1)+uu(i,j,k  ,1))/2.d0+
     &                bb*(uu(i,j,k-2,1)+uu(i,j,k+1,1))/2.d0
           enddo
         enddo
         do j=0, nylocal+2, 1
           do k= 0, nzlocal+2, 1
           dudx(i,j,k,1,2)=aa*(ux_z(j  ,k)-ux_z(j-1,k))*dyi+
     &                     bb*(ux_z(j+1,k)-ux_z(j-2,k))*dyi3
           enddo 
         enddo
       enddo
c       (u,z)
       do i=0,nxlocal+2,1
         do j=0, nylocal+2, 1
           do k= -2, nzlocal+3, 1
           ux_y(j,k)= aa*(uu(i,j-1,k,1)+uu(i,j  ,k,1))/2.d0+
     &                bb*(uu(i,j-2,k,1)+uu(i,j+1,k,1))/2.d0
           enddo
         enddo
         do j=0, nylocal+2, 1
           do k= 0, nzlocal+2, 1
           dudx(i,j,k,1,3)=aa*(ux_y(j,k  )-ux_y(j,k-1))*dzi+
     &                     bb*(ux_y(j,k+1)-ux_y(j,k-2))*dzi3
           enddo 
         enddo
       enddo
c
c      (v,x)
       do j=0, nylocal+2, 1
         do i=-2, nxlocal+3, 1
           do k= 0,nzlocal+2, 1
           uy_z(i,k)= aa*(uu(i,j,k-1,2)+uu(i  ,j,k,2))/2.d0+
     &                bb*(uu(i,j,k-2,2)+uu(i+1,j,k,2))/2.d0
           enddo
         enddo
         do i=0, nxlocal+2, 1
           do k=0, nzlocal+2, 1
           dudx(i,j,k,2,1)=aa*(uy_z(i  ,k)-uy_z(i-1,k))*dxi+
     &                     bb*(uy_z(i+1,k)-uy_z(i-2,k))*dxi3
           enddo
         enddo
       enddo
c      (v,y)
       do j=-2, nylocal+4, 1
         do i=-1, nxlocal+3,1
           do k= -3, nzlocal+4, 1
           tmpy(i,k)= aa*(uu(i-1,j,k,2)+uu(i  ,j,k,2))/2.d0+
     &                bb*(uu(i-2,j,k,2)+uu(i+1,j,k,2))/2.d0
           enddo
         enddo
         do i=-1, nxlocal+3, 1
           do k=-1, nzlocal+3, 1
           vca(i,j,k)=aa*(tmpy(i,k-1)+tmpy(i,k  ))/2.d0+
     &                  bb*(tmpy(i,k-2)+tmpy(i,k+1))/2.d0
           enddo
         enddo
       enddo 
       do i=0, nxlocal+2, 1
         do j=0, nylocal+2, 1
           do k=0, nzlocal+2, 1
           dudx(i,j,k,2,2)=cc*(vca(i,j+1,k)-vca(i,j-1,k))*dyi+
     &                     dd*(vca(i,j+2,k)-vca(i,j-2,k))*dyi
           enddo
         enddo
       enddo
c      (v,z)
       do j=0, nylocal+2, 1
         do i=0, nxlocal+2, 1
           do k=-2, nzlocal+3, 1
           uy_x(i,k)= aa*(uu(i-1,j,k,2)+uu(i  ,j,k,2))/2.d0+
     &                bb*(uu(i-2,j,k,2)+uu(i+1,j,k,2))/2.d0
           enddo
         enddo
         do i=0, nxlocal+2, 1
           do k=0, nzlocal+2, 1
           dudx(i,j,k,2,3)=aa*(uy_x(i,k  )-uy_x(i,k-1))*dzi+
     &                     bb*(uy_x(i,k+1)-uy_x(i,k-2))*dzi3
           enddo
         enddo
       enddo

c      (w,x)
       do k=0, nzlocal+2, 1
         do i=-2, nxlocal+3, 1
           do j=0, nylocal+2, 1
           uz_y(i,j)= aa*(uu(i,j-1,k,3)+uu(i,j  ,k,3))/2.d0+
     &                bb*(uu(i,j-2,k,3)+uu(i,j+1,k,3))/2.d0
           enddo
         enddo 
         do i=0, nxlocal+2, 1
           do j=0, nylocal+2,1 
           dudx(i,j,k,3,1)=aa*(uz_y(i  ,j)-uz_y(i-1,j))*dxi+
     &                     bb*(uz_y(i+1,j)-uz_y(i-2,j))*dxi3
           enddo
         enddo
       enddo 
c      (w,y)
       do k=0, nzlocal+2, 1
         do i=0, nxlocal+2, 1
           do j= -2, nylocal+3,1
           uz_x(i,j)= aa*(uu(i-1,j,k,3)+uu(i  ,j,k,3))/2.d0+
     &                bb*(uu(i-2,j,k,3)+uu(i+1,j,k,3))/2.d0
           enddo
         enddo
         do i=0, nxlocal+2, 1
           do j=0, nylocal+2, 1
           dudx(i,j,k,3,2)=aa*(uz_x(i,j  )-uz_x(i,j-1))*dyi+
     &                     bb*(uz_x(i,j+1)-uz_x(i,j-2))*dyi3
           enddo
         enddo
       enddo
c      (w,z)
       do k=-2, nzlocal+4, 1
         do i=-1, nxlocal+3, 1
           do j=-3, nylocal+4, 1
           tmpz(i,j)= aa*(uu(i-1,j,k,3)+uu(i  ,j,k,3))/2.d0+
     &                bb*(uu(i-2,j,k,3)+uu(i+1,j,k,3))/2.d0
           enddo
         enddo
         do i=-1, nxlocal+3, 1
           do j=-1, nylocal+3, 1
           wca(i,j,k)=aa*(tmpz(i,j-1)+tmpz(i,j  ))/2.d0+
     &                  bb*(tmpz(i,j-2)+tmpz(i,j+1))/2.d0
           enddo
         enddo
       enddo
       do i=0, nxlocal+2, 1
         do j=0, nylocal+2, 1
           do k=0, nzlocal+2, 1
           dudx(i,j,k,3,3)=cc*(wca(i,j,k+1)-wca(i,j,k-1))*dzi+
     &                     dd*(wca(i,j,k+2)-wca(i,j,k-2))*dzi
           enddo
         enddo
       enddo
c
c      velocites at corner
       ucn(-1:nxlocal+3,-1:nylocal+3,-1:nzlocal+3,1) =
     & uca(-1:nxlocal+3,-1:nylocal+3,-1:nzlocal+3)
       ucn(-1:nxlocal+3,-1:nylocal+3,-1:nzlocal+3,2) =
     & vca(-1:nxlocal+3,-1:nylocal+3,-1:nzlocal+3)
       ucn(-1:nxlocal+3,-1:nylocal+3,-1:nzlocal+3,3) = 
     & wca(-1:nxlocal+3,-1:nylocal+3,-1:nzlocal+3)
c


c      for wall boundary condition
       do k= 1, 2, 1
         do i=1, nxlocal+1, 1
           do j=1, nylocal
           bctmp(i,j,k,1)=aa*(ucn(i,j  ,k,1)+ucn(i,j+1,k,1))/2.d0+
     &                    bb*(ucn(i,j-1,k,1)+ucn(i,j+2,k,1))/2.d0
           bctmp(i,j,k,2)=aa*(ucn(i,j  ,k,3)+ucn(i,j+1,k,3))/2.d0+
     &                    bb*(ucn(i,j-1,k,3)+ucn(i,j+2,k,3))/2.d0
           enddo
         enddo
       enddo
       k =1
       do i=1, nxlocal+1, 1
         do j=1, nylocal,1
         bc_uw(i,j)=   ( bctmp(i,j,k  ,1)*bctmp(i,j,k+1,2) +
     &                   bctmp(i,j,k+1,1)*bctmp(i,j,k  ,2) )/2.d0 
         enddo
       enddo 
c
       do k= 1, 2, 1
         do i=1, nxlocal+1, 1
           do j=1, nylocal
           bcduzt(i,j,k)=
     &           aa*(dudx(i,j  ,k,1,3)+dudx(i,j+1,k,1,3))/2.d0+
     &           bb*(dudx(i,j-1,k,1,3)+dudx(i,j+2,k,1,3))/2.d0
           enddo
         enddo
       enddo
       do i=1, nxlocal+1, 1
         do j=1, nylocal,1
         bc_dudzh(i,j)= (bcduzt(i,j,1)+bcduzt(i,j,2))/2.d0
         enddo
       enddo 


c
       return
       end
c
c----------------------------------------------------------------------
c      for wall processor
c----------------------------------------------------------------------
       subroutine interpolate_ududx_cn_wall(uu,dudx,ucn)
       use mesh
       use mesh_common
       use properties
       use virtualwallBC
       use InflowVelBC
#ifdef PARALLEL
       use mpistuff
#endif
c
       implicit none
c=======Declarations=========
       double precision:: uu(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
       double precision:: ucn(-1:nxlocal+3,-1:nylocal+3,-1:nzlocal+3,3)
       double precision:: dudx(0:nxlocal+2,0:nylocal+2,0:nzlocal+2,3,3)
       double precision:: uca(-2:nxlocal+4, -1:nylocal+3, -1:nzlocal+3)
       double precision:: vca(-1:nxlocal+3, -2:nylocal+4, -1:nzlocal+3)
       double precision:: wca(-1:nxlocal+3, -1:nylocal+3, -2:nzlocal+4)
       double precision:: tmpx(-1:nylocal+3,-3:nzlocal+4)
       double precision:: tmpy(-1:nxlocal+3,-3:nzlocal+4)
       double precision:: tmpz(-1:nxlocal+3,-3:nylocal+4)
       double precision:: ux_y( 0:nylocal+2,-2:nzlocal+3)
       double precision:: ux_z(-2:nylocal+3, 0:nzlocal+2)
       double precision:: uy_x( 0:nxlocal+2,-2:nzlocal+3)
       double precision:: uy_z(-2:nxlocal+3, 0:nzlocal+2)
       double precision:: uz_x( 0:nxlocal+2,-2:nylocal+3)
       double precision:: uz_y(-2:nxlocal+3, 0:nylocal+2)
       double precision:: bctmp(1:nxlocal+1,1:nylocal,1:2,2)
       double precision:: bcduzt(1:nxlocal+1,1:nylocal,2)
       double precision:: aa,bb,cc,dd
       integer:: i,j,k, ilo, klo, ihi, khi
c
       aa= 9.D0/8.D0
       bb=-1.D0/8.D0
       cc= 8.d0/12.d0
       dd=-1.d0/12.d0
c---------------------------------------------------------------------
c
       do i=-2,nxlocal+4,1
c   uca : k=1
         k = 1
         do j= -1, nylocal+3, 1
           uca(i,j,k) = 
     &        aa*(bc_uwall(i,j  ,1)+bc_uwall(i,j-1,1))/2.d0 +
     &        bb*(bc_uwall(i,j+1,1)+bc_uwall(i,j-2,1))/2.d0
         enddo
c   uca : k=3, ...
         do j=-1, nylocal+3, 1
           do k= 0, nzlocal+4, 1
           tmpx(j,k)= aa*(uu(i,j-1,k,1)+uu(i,j  ,k,1))/2.d0+
     &                bb*(uu(i,j-2,k,1)+uu(i,j+1,k,1))/2.d0
           enddo
         enddo
         do j=-1, nylocal+3, 1
           do k= 2, nzlocal+3, 1
           uca(i,j,k)=aa*(tmpx(j,k-1)+tmpx(j,k  ))/2.d0+
     &                bb*(tmpx(j,k-2)+tmpx(j,k+1))/2.d0
           enddo 
         enddo
c   uca : k=2
        k = 2
         do j= -1, nylocal+3, 1
           uca(i,j,k) = 1.d0/20.d0*( 
     &         - 4.d0 * uca(i,j,1)
     &         +15.d0 * tmpx(j,1)
     &         +10.d0 * tmpx(j,2)
     &         - 1.d0 * tmpx(j,3) )
         enddo 
       enddo
c   dudx(u,x)
       do i=0, nxlocal+2, 1 
         do j=0, nylocal+2, 1
           do k=1, nzlocal+2, 1
           dudx(i,j,k,1,1)=cc*(uca(i+1,j,k)-uca(i-1,j,k))*dxi+
     &                     dd*(uca(i+2,j,k)-uca(i-2,j,k))*dxi
           enddo
         enddo
       enddo
c    dudx(u,y)
       do i=0,nxlocal+2,1
         do j= -2, nylocal+3, 1
           do k= 2, nzlocal+2, 1
           ux_z(j,k)= aa*(uu(i,j,k-1,1)+uu(i,j,k  ,1))/2.d0+
     &                bb*(uu(i,j,k-2,1)+uu(i,j,k+1,1))/2.d0
           enddo
         enddo
         k = 2
         do j=-2, nylocal+3, 1
           ux_z(j,k)= 1.d0/20.d0*( 
     &           - 4.d0 * bc_uwall(i,j,1)
     &           +15.d0 * uu(i,j,1,1)
     &           +10.d0 * uu(i,j,2,1)
     &           - 1.d0 * uu(i,j,3,1) )
         enddo
         do j=0, nylocal+2, 1
           k = 1
            dudx(i,j,k,1,2) = 
     &         aa*(bc_uwall(i,j  ,1)-bc_uwall(i,j-1,1))*dyi +
     &         bb*(bc_uwall(i,j+1,1)-bc_uwall(i,j-2,1))*dyi3
           do k= 2, nzlocal+2, 1
           dudx(i,j,k,1,2)=aa*(ux_z(j  ,k)-ux_z(j-1,k))*dyi+
     &                     bb*(ux_z(j+1,k)-ux_z(j-2,k))*dyi3
           enddo 
         enddo
       enddo

c       (u,z)
       do i=0,nxlocal+2,1
         do j=0, nylocal+2, 1
           do k= -1, nzlocal+3, 1
           ux_y(j,k)= aa*(uu(i,j-1,k,1)+uu(i,j  ,k,1))/2.d0+
     &                bb*(uu(i,j-2,k,1)+uu(i,j+1,k,1))/2.d0
           enddo
         enddo
         do j=0, nylocal+2, 1
           do k= 1, nzlocal+2, 1
           dudx(i,j,k,1,3)=aa*(ux_y(j,k  )-ux_y(j,k-1))*dzi+
     &                     bb*(ux_y(j,k+1)-ux_y(j,k-2))*dzi3
           enddo 
         enddo
         k=1
         do j=0, nylocal+2, 1
            dudx(i,j,k,1,3) = 1.d0/6.d0/dz*(
     &        -11.d0 * uca(i,j,k  ) + 18.d0 * uca(i,j,k+1)
     &        - 9.d0 * uca(i,j,k+2) +  2.d0 * uca(i,j,k+3) )
            dudx(i,j,k+1,1,3) = 1.d0/6.d0/dz*(
     &        - 2.d0 * uca(i,j,k  ) -  3.d0 * uca(i,j,k+1)
     &        + 6.d0 * uca(i,j,k+2) -  1.d0 * uca(i,j,k+3) )
         enddo
       enddo
c
c-------------------------------------------------------------
c  v
       do j=-2, nylocal+4, 1
         do i=-1, nxlocal+3,1
           do k= -3, nzlocal+4, 1
           tmpy(i,k)= aa*(uu(i-1,j,k,2)+uu(i  ,j,k,2))/2.d0+
     &                bb*(uu(i-2,j,k,2)+uu(i+1,j,k,2))/2.d0
           enddo
         enddo
         do i=-1, nxlocal+3, 1
           do k=-1, nzlocal+3, 1
           vca(i,j,k)=aa*(tmpy(i,k-1)+tmpy(i,k  ))/2.d0+
     &                  bb*(tmpy(i,k-2)+tmpy(i,k+1))/2.d0
           enddo
         enddo
       enddo 
       do i=0, nxlocal+2, 1
         do j=0, nylocal+2, 1
           do k=0, nzlocal+2, 1
           dudx(i,j,k,2,2)=cc*(vca(i,j+1,k)-vca(i,j-1,k))*dyi+
     &                     dd*(vca(i,j+2,k)-vca(i,j-2,k))*dyi
           enddo
         enddo
       enddo
c   dudx(v,y)
       do i=0, nxlocal+2, 1
         do j=0, nylocal+2, 1
           do k=1, nzlocal+2, 1
           dudx(i,j,k,2,2)=cc*(vca(i,j+1,k)-vca(i,j-1,k))*dyi+
     &                     dd*(vca(i,j+2,k)-vca(i,j-2,k))*dyi
           enddo
         enddo
       enddo
c    dudx(v,x)
       do j=0, nylocal+2, 1
         do i=-2, nxlocal+3, 1
           do k= 2,nzlocal+2, 1
           uy_z(i,k)= aa*(uu(i,j,k-1,2)+uu(i  ,j,k,2))/2.d0+
     &                bb*(uu(i,j,k-2,2)+uu(i+1,j,k,2))/2.d0
           enddo
         enddo
c         k = 2
c         do i=-2, nxlocal+3, 1
c           uy_z(j,k)= 1.d0/20.d0*( 
c     &           - 4.d0 *  0.d0 !bc_uwall(i,j,2)
c     &           +15.d0 * uu(i,j,1,2)
c     &           +10.d0 * uu(i,j,2,2)
c     &           - 1.d0 * uu(i,j,3,2) )
c         enddo
         do i=0, nxlocal+2, 1
           k = 1
           dudx(i,j,k,2,1) =  0.d0
c     &         aa*(bc_uwall(i  ,j,2)-bc_uwall(i-1,j,2))*dxi +
c     &         bb*(bc_uwall(i+1,j,2)-bc_uwall(i-2,j,2))*dxi3
           do k=2, nzlocal+2, 1
           dudx(i,j,k,2,1)=aa*(uy_z(i  ,k)-uy_z(i-1,k))*dxi+
     &                     bb*(uy_z(i+1,k)-uy_z(i-2,k))*dxi3
           enddo
         enddo
       enddo

c   dudx(v,z)
       do j=0, nylocal+2, 1
         do i=0, nxlocal+2, 1
           do k=-1, nzlocal+3, 1
           uy_x(i,k)= aa*(uu(i-1,j,k,2)+uu(i  ,j,k,2))/2.d0+
     &                bb*(uu(i-2,j,k,2)+uu(i+1,j,k,2))/2.d0
           enddo
         enddo
         do i=0, nxlocal+2, 1
           do k=1, nzlocal+2, 1
           dudx(i,j,k,2,3)=aa*(uy_x(i,k  )-uy_x(i,k-1))*dzi+
     &                     bb*(uy_x(i,k+1)-uy_x(i,k-2))*dzi3
           enddo
         enddo
         k=1
         do i=0, nxlocal+2, 1
            dudx(i,j,k,2,3) = 1.d0/6.d0/dz*(
     &        -11.d0 * vca(i,j,k  ) + 18.d0 * vca(i,j,k+1)
     &        - 9.d0 * vca(i,j,k+2) +  2.d0 * vca(i,j,k+3) )
            dudx(i,j,k+1,2,3) = 1.d0/6.d0/dz*(
     &        - 2.d0 * vca(i,j,k  ) -  3.d0 * vca(i,j,k+1)
     &        + 6.d0 * vca(i,j,k+2) -  1.d0 * vca(i,j,k+3) )
         enddo
       enddo
c
c-------------------------------------------------------------
c  w  
       do k=-1, nzlocal+4, 1
         do i=-1, nxlocal+3, 1
           do j=-3, nylocal+4, 1
           tmpz(i,j)= aa*(uu(i-1,j,k,3)+uu(i  ,j,k,3))/2.d0+
     &                bb*(uu(i-2,j,k,3)+uu(i+1,j,k,3))/2.d0
           enddo
         enddo
         do i=-1, nxlocal+3, 1
           do j=-1, nylocal+3, 1
           wca(i,j,k)=aa*(tmpz(i,j-1)+tmpz(i,j  ))/2.d0+
     &                  bb*(tmpz(i,j-2)+tmpz(i,j+1))/2.d0
           enddo
         enddo
       enddo
c   dudx(w,z)
       do i=0, nxlocal+2, 1
         do j=0, nylocal+2, 1
           do k=1, nzlocal+2, 1
           dudx(i,j,k,3,3)=cc*(wca(i,j,k+1)-wca(i,j,k-1))*dzi+
     &                     dd*(wca(i,j,k+2)-wca(i,j,k-2))*dzi
           enddo
           k=1
           dudx(i,j,k,3,3) = 1.d0/6.d0/dz*(
     &        -11.d0 * wca(i,j,k  ) + 18.d0 * wca(i,j,k+1)
     &        - 9.d0 * wca(i,j,k+2) +  2.d0 * wca(i,j,k+3) )
           dudx(i,j,k+1,3,3) = 1.d0/6.d0/dz*(
     &        - 2.d0 * wca(i,j,k  ) -  3.d0 * wca(i,j,k+1)
     &        + 6.d0 * wca(i,j,k+2) -  1.d0 * wca(i,j,k+3) )
         enddo
       enddo
c      (w,x)
       do k=1, nzlocal+2, 1
         do i=-2, nxlocal+3, 1
           do j=0, nylocal+2, 1
           uz_y(i,j)= aa*(uu(i,j-1,k,3)+uu(i,j  ,k,3))/2.d0+
     &                bb*(uu(i,j-2,k,3)+uu(i,j+1,k,3))/2.d0
           enddo
         enddo 
         do i=0, nxlocal+2, 1
           do j=0, nylocal+2,1 
           dudx(i,j,k,3,1)=aa*(uz_y(i  ,j)-uz_y(i-1,j))*dxi+
     &                     bb*(uz_y(i+1,j)-uz_y(i-2,j))*dxi3
           enddo
         enddo
       enddo 
c      (w,y)
       do k=1, nzlocal+2, 1
         do i=0, nxlocal+2, 1
           do j= -2, nylocal+3,1
           uz_x(i,j)= aa*(uu(i-1,j,k,3)+uu(i  ,j,k,3))/2.d0+
     &                bb*(uu(i-2,j,k,3)+uu(i+1,j,k,3))/2.d0
           enddo
         enddo
         do i=0, nxlocal+2, 1
           do j=0, nylocal+2, 1
           dudx(i,j,k,3,2)=aa*(uz_x(i,j  )-uz_x(i,j-1))*dyi+
     &                     bb*(uz_x(i,j+1)-uz_x(i,j-2))*dyi3
           enddo
         enddo
       enddo

c-----------------------------------------------------
c      velocites at corner
       ucn(-1:nxlocal+3,-1:nylocal+3, 1:nzlocal+3,1) =
     & uca(-1:nxlocal+3,-1:nylocal+3, 1:nzlocal+3)
       ucn(-1:nxlocal+3,-1:nylocal+3, 1:nzlocal+3,2) =
     & vca(-1:nxlocal+3,-1:nylocal+3, 1:nzlocal+3)
       ucn(-1:nxlocal+3,-1:nylocal+3, 1:nzlocal+3,3) = 
     & wca(-1:nxlocal+3,-1:nylocal+3, 1:nzlocal+3)
c-----------------------------------------------------
c
c only using inner points 
c       call setboundary_ududx_cn(uu,dudx,ucn)


c      for wall boundary condition
       do k= 1, 2, 1
         do i=1, nxlocal+1, 1
           do j=1, nylocal
           bctmp(i,j,k,1)=aa*(ucn(i,j  ,k,1)+ucn(i,j+1,k,1))/2.d0+
     &                    bb*(ucn(i,j-1,k,1)+ucn(i,j+2,k,1))/2.d0
           bctmp(i,j,k,2)=aa*(ucn(i,j  ,k,3)+ucn(i,j+1,k,3))/2.d0+
     &                    bb*(ucn(i,j-1,k,3)+ucn(i,j+2,k,3))/2.d0
           enddo
         enddo
       enddo
       k =1
       do i=1, nxlocal+1, 1
         do j=1, nylocal,1
         bc_uw(i,j)=   ( bctmp(i,j,k  ,1)*bctmp(i,j,k+1,2) +
     &                   bctmp(i,j,k+1,1)*bctmp(i,j,k  ,2) )/2.d0 
         enddo
       enddo 
c
       do k= 1, 2, 1
         do i=1, nxlocal+1, 1
           do j=1, nylocal
           bcduzt(i,j,k)=
     &           aa*(dudx(i,j  ,k,1,3)+dudx(i,j+1,k,1,3))/2.d0+
     &           bb*(dudx(i,j-1,k,1,3)+dudx(i,j+2,k,1,3))/2.d0
           enddo
         enddo
       enddo
       do i=1, nxlocal+1, 1
         do j=1, nylocal,1
         bc_dudzh(i,j)= (bcduzt(i,j,1)+bcduzt(i,j,2))/2.d0
         enddo
       enddo 
c
       return
       end
c
c
c----------------------------------------------------------------------
c     interpolation of tensor 
c----------------------------------------------------------------------
       subroutine interpolateT(Tijn, Tij)
       use mesh
       use mesh_common
       use properties
#ifdef PARALLEL
       use mpistuff
#endif       
c
       implicit none
       double precision:: Tij(-1:nxlocal+2,-1:nylocal+2,-1:nzlocal+2,6) 
       double precision:: Tijn(0:nxlocal+2,0:nylocal+2,0:nzlocal+2,6) 
       double precision:: tmpx(1:nxlocal,0:nylocal+2,0:nzlocal+2,6)
       double precision:: tmpy(1:nxlocal,1:nylocal,0:nzlocal+2,6)
       double precision:: Ttmp
       double precision:: Ti(1:nylocal,0:nzlocal+2,6)
c       double precision:: Tj(1:nxlocal,0:nzlocal+2,6)
       double precision:: Tk(1:nxlocal,0:nylocal+2,6)
       double precision aa, bb  !, cc, dd
       integer:: i, j, k, l, ig, jg, kg
c
       aa = 9.d0/8.d0 ! for  central diff. of 1st derivative at midpoint
       bb = -1.d0/8.d0 
c
c------------------------------------------------------------------------
      if(iprocx.eq.1) then
      i = 0
      do j=0, nylocal+2, 1
        do k=1, nzlocal+1, 1
          Tijn(i,j,k,1:6) = 
     &        4.d0*( Tijn(i+3,j,k,1:6)+Tijn(i+1,j,k,1:6) )
     &              -Tijn(i+4,j,k,1:6)
     &         -6.d0*Tijn(i+2,j,k,1:6)
        enddo
      enddo
      endif
      if(iprocx.eq.xprocs) then
      i=nxlocal+2
      do j=0, nylocal+2, 1
        do k=1, nzlocal+1, 1
          Tijn(i,j,k,1:6) =
     &        4.d0*( Tijn(i-3,j,k,1:6)+Tijn(i-1,j,k,1:6) )
     &              -Tijn(i-4,j,k,1:6)
     &         -6.d0*Tijn(i-2,j,k,1:6)
        enddo
      enddo
      endif
      if(iprocz.eq.1) then  
        k = 0
        do i=0,nxlocal+2,1
          do j=0,nylocal+2,1
          Tijn(i,j,k,1:6) = 
     &        4.d0*( Tijn(i,j,k+3,1:6)+Tijn(i,j,k+1,1:6) )
     &              -Tijn(i,j,k+4,1:6)
     &         -6.d0*Tijn(i,j,k+2,1:6) 
          enddo 
         enddo
       endif
c
c p1 filter
       call ffttensor(Tijn)
c--------------------------------------------------
c part 1 Interpolation from Tijn to Tij
       do i = 1, nxlocal, 1
         do j = 0, nylocal+2, 1
           do k = 0, nzlocal+2, 1
             do l = 1, 6, 1
             tmpx(i,j,k,l)=aa*(Tijn(i  ,j,k,l)+Tijn(i+1,j,k,l))/2.d0+
     &                     bb*(Tijn(i-1,j,k,l)+Tijn(i+2,j,k,l))/2.d0
             enddo  
           enddo
         enddo
       enddo      
c       
       do i = 1, nxlocal, 1
         do j = 1, nylocal, 1
           do k = 0, nzlocal+2, 1
             do l = 1, 6, 1
             tmpy(i,j,k,l)=aa*(tmpx(i,j  ,k,l)+tmpx(i,j+1,k,l))/2.d0+
     &                     bb*(tmpx(i,j-1,k,l)+tmpx(i,j+2,k,l))/2.d0  
             enddo
           enddo
         enddo
       enddo     
c
       do i = 1, nxlocal, 1
         do j = 1, nylocal, 1
           do k = 1, nzlocal, 1
             do l = 1, 6, 1
             Tij(i,j,k,l) =aa*(tmpy(i,j,k  ,l)+tmpy(i,j,k+1,l))/2.d0+
     &                     bb*(tmpy(i,j,k-1,l)+tmpy(i,j,k+2,l))/2.d0  
             enddo
           enddo
         enddo
       enddo        
c

c Part 2 parallel exchange Tij at center
#ifdef PARALLEL
       call exchangeTct(Tij,6,-1,nxlocal+2,-1,nylocal+2,-1,nzlocal+2)
#endif


c-----------------------------------------------------------------------
c part 3   B.C. :  Fill in Tij
c
c      i_proc = 1      
       ! i = 0
      if(iprocx .eq. 1) then
       ig = 0
       do j = 1, nylocal, 1
         do k = 0, nzlocal+2, 1
           do l= 1, 6, 1
           Ti(j,k,l)=aa*(Tijn(ig+2,j  ,k,l)+Tijn(ig+2,j+1,k,l))/2.d0+
     &               bb*(Tijn(ig+2,j-1,k,l)+Tijn(ig+2,j+2,k,l))/2.d0 
           enddo
         enddo
       enddo
       do j = 1, nylocal, 1
         do k = 1, nzlocal, 1
           do l = 1, 6, 1
           Ttmp=aa*(Ti(j,k,l)  +Ti(j,k+1,l))/2.d0+
     &          bb*(Ti(j,k-1,l)+Ti(j,k+2,l))/2.d0
           Tij(ig,j,k,l) = 
     &         - 16.d0*Ttmp
     &         +  9.d0*(Tij(ig+1,j,k,l)+Tij(ig+2,j,k,l))
     &         -  1.d0* Tij(ig+3,j,k,l) 
           enddo
         enddo
       enddo
       ! i=-1
       ig = -1
       do j = 1, nylocal, 1
         do k = 0, nzlocal+2, 1
           do l= 1, 6, 1
           Ti(j,k,l)=aa*(Tijn(ig+2,j  ,k,l)+Tijn(ig+2,j+1,k,l))/2.d0+
     &               bb*(Tijn(ig+2,j-1,k,l)+Tijn(ig+2,j+2,k,l))/2.d0 
           enddo
         enddo
       enddo
       do j = 1, nylocal, 1
         do k = 1, nzlocal, 1
           do l = 1, 6, 1
           Ttmp=aa*(Ti(j,k,l)  +Ti(j,k+1,l))/2.d0+
     &          bb*(Ti(j,k-1,l)+Ti(j,k+2,l))/2.d0
           Tij(ig,j,k,l) = 
     &         - 16.d0*Ttmp
     &         +  9.d0*(Tij(ig+1,j,k,l)+Tij(ig+2,j,k,l))
     &         -  1.d0* Tij(ig+3,j,k,l) 
           enddo
         enddo
       enddo
      endif
c
c      i_proc = i_max
       ! i = nxlocal +1
      if(iprocx.eq. xprocs) then
       ig = nxlocal+1
       do j=1, nylocal, 1
         do k = 0, nzlocal+2, 1
           do l = 1, 6, 1
           Ti(j,k,l)=aa*(Tijn(ig-1,j  ,k,l)+Tijn(ig-1,j+1,k,l))/2.d0+
     &               bb*(Tijn(ig-1,j-1,k,l)+Tijn(ig-1,j+2,k,l))/2.d0 
           enddo
         enddo
       enddo
       do j = 1, nylocal, 1
         do k = 1, nzlocal, 1
           do l = 1, 6, 1
           Ttmp=aa*(Ti(j,k,l)  +Ti(j,k+1,l))/2.d0+
     &          bb*(Ti(j,k-1,l)+Ti(j,k+2,l))/2.d0
           Tij(ig,j,k,l)= 
     &        - 16.d0*Ttmp
     &        +  9.d0*(Tij(ig-1,j,k,l)+Tij(ig-2,j,k,l))
     &        -  1.d0* Tij(ig-3,j,k,l)
           enddo
         enddo
       enddo
       ! i = nxlocal + 2
       ig = nxlocal+2
       do j=1, nylocal, 1
         do k=0, nzlocal+2, 1
           do l=1, 6, 1
           Ti(j,k,l)=aa*(Tijn(ig-1,j  ,k,l)+Tijn(ig-1,j+1,k,l))/2.d0+
     &               bb*(Tijn(ig-1,j-1,k,l)+Tijn(ig-1,j+2,k,l))/2.d0
           enddo
         enddo
       enddo
       do j=1, nylocal, 1
         do k=1, nzlocal, 1
           do l=1, 6, 1
           Ttmp=aa*(Ti(j,k  ,l)+Ti(j,k+1,l))/2.d0+
     &          bb*(Ti(j,k-1,l)+Ti(j,k+2,l))/2.d0
           Tij(ig,j,k,l) =
     &        - 16.d0*Ttmp
     &        +  9.d0*(Tij(ig-1,j,k,l)+Tij(ig-2,j,k,l))
     &        -  1.d0* Tij(ig-3,j,k,l)
           enddo
         enddo
       enddo
      endif
c
c     periodic in y direction
c     considering only one process
c      j = 0
      if(yprocs.eq.1) then  
       do i=1, nxlocal, 1
         do k=1, nzlocal, 1
           do l=1, 6, 1
           jg = 0
             Tij(i,jg,k,l)= Tij(i,nylocal+jg,k,l)
           jg = -1
             Tij(i,jg,k,l)= Tij(i,nylocal+jg,k,l)
           jg = nylocal+1
             Tij(i,jg,k,l)= Tij(i,jg-nylocal,k,l)
           jg = nylocal+2
             Tij(i,jg,k,l)= Tij(i,jg-nylocal,k,l)
           enddo
         enddo
       enddo
      endif
c
c      k_proc = 1
       ! k=0
      if(iprocz.eq.1)then
       kg = 0
       do i = 1, nxlocal, 1
         do j = 0, nylocal+2, 1
           do l= 1, 6, 1
           Tk(i,j,l)=aa*(Tijn(i  ,j,kg+2,l)+Tijn(i+1,j,kg+2,l))/2.d0+
     &               bb*(Tijn(i-1,j,kg+2,l)+Tijn(i+2,j,kg+2,l))/2.d0
           enddo
         enddo
       enddo
       do i = 1, nxlocal, 1
         do j= 1, nylocal, 1
           do l = 1, 6, 1
           Ttmp=aa*(Tk(i,j  ,l)+Tk(i,j+1,l))/2.d0+
     &          bb*(Tk(i,j-1,l)+Tk(i,j+2,l))/2.d0
           Tij(i,j,kg,l)= 
     &         - 16.d0*Ttmp
     &         +  9.d0*(Tij(i,j,kg+1,l)+Tij(i,j,kg+2,l))
     &         -  1.d0* Tij(i,j,kg+3,l)     
           enddo
         enddo
       enddo
       ! k=-1
       kg = -1
       do i = 1, nxlocal, 1
         do j = 0, nylocal+2, 1
           do l= 1, 6, 1
           Tk(i,j,l)=aa*(Tijn(i  ,j,kg+2,l)+Tijn(i+1,j,kg+2,l))/2.d0+
     &               bb*(Tijn(i-1,j,kg+2,l)+Tijn(i+2,j,kg+2,l))/2.d0
           enddo
         enddo
       enddo
       do i = 1, nxlocal, 1
         do j= 1, nylocal, 1
           do l = 1, 6, 1
           Ttmp=aa*(Tk(i,j  ,l)+Tk(i,j+1,l))/2.d0+
     &          bb*(Tk(i,j-1,l)+Tk(i,j+2,l))/2.d0
           Tij(i,j,kg,l)=  
     &         - 16.d0*Ttmp
     &         +  9.d0*(Tij(i,j,kg+1,l)+Tij(i,j,kg+2,l))
     &         -  1.d0* Tij(i,j,kg+3,l)     
           enddo
         enddo
       enddo
      endif
c      k_proc = k_max
       ! k=nzlocal+1
      if(iprocz.eq. zprocs) then
       kg = nzlocal+1
       do i = 1, nxlocal, 1
         do j = 0, nylocal+2, 1
           do l= 1, 6, 1
           Tk(i,j,l)=aa*(Tijn(i  ,j,kg-1,l)+Tijn(i+1,j,kg-1,l))/2.d0+
     &               bb*(Tijn(i-1,j,kg-1,l)+Tijn(i+2,j,kg-1,l))/2.d0
           enddo
         enddo
       enddo
       do i = 1, nxlocal, 1
         do j= 1, nylocal, 1
           do l = 1, 6, 1
           Ttmp=aa*(Tk(i,j  ,l)+Tk(i,j+1,l))/2.d0+
     &          bb*(Tk(i,j-1,l)+Tk(i,j+2,l))/2.d0
           Tij(i,j,kg,l)=
     &         - 16.d0*Ttmp
     &         +  9.d0*(Tij(i,j,kg-1,l)+Tij(i,j,kg-2,l))
     &         -  1.d0* Tij(i,j,kg-3,l)     
           enddo
         enddo
       enddo
       ! k = nzlocal +2
       kg = nzlocal +2
       do i=1, nxlocal, 1
         do j= 0, nylocal+2, 1
           do l= 1, 6, 1
           Tk(i,j,l)=aa*(Tijn(i  ,j,kg-1,l)+Tijn(i+1,j,kg-1,l))/2.d0+
     &               bb*(Tijn(i-1,j,kg-1,l)+Tijn(i+2,j,kg-1,l))/2.d0
           enddo
         enddo
       enddo
       do i = 1, nxlocal, 1
         do j= 1, nylocal, 1
           do l = 1, 6, 1
           Ttmp=aa*(Tk(i,j  ,l)+Tk(i,j+1,l))/2.d0+
     &          bb*(Tk(i,j-1,l)+Tk(i,j+2,l))/2.d0
           Tij(i,j,kg,l)=
     &         - 16.d0*Ttmp
     &         +  9.d0*(Tij(i,j,kg-1,l)+Tij(i,j,kg-2,l))
     &         -  1.d0* Tij(i,j,kg-3,l)     
           enddo
         enddo
       enddo
      endif
c
       return
       end
c
c----------------------------------------------------------------------
