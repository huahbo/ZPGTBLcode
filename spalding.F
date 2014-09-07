

        subroutine spalding(ux)
        use mesh
        use mesh_common
        use properties
        use virtualwallBC
        use InflowVelBC
c
        implicit none
        double precision:: ux(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
        double precision:: zp, up, vk, vb, va, zeta, delta
        double precision:: utau, nu, pai, Retheta, uep
        double precision:: zv, uv
        integer:: i,j,k
c
        double precision:: h0, hnup, usbc        
c        
        Retheta = Reynolds/11.0d0
        delta = 1.d0
c
        pai = 4.d0*datan(1.d0)
        vk = 0.378d0
        vb = 3.8d0
        va = 2.5d0
c
        h0 = bc_h0dz * dz
c
        nu = 1.d0/Reynolds
        uep = 1.d0/vk*dlog(Retheta)+vb
        utau = 1.d0/uep
c
        do i=ixlo, ixhi, 1
          do j=iylo, iyhi , 1
            do k=1, izhi, 1
              zv = zc(k) + h0
              zeta = zv/delta
              zp= zv*utau/nu
              up=1.d0/vk*dlog(zp)+vb+va*dsin(pai/2.d0*zeta)**2.d0
              uv = up*utau
              if(zc(k).gt.1.2d0 .or. uv.gt.1.d0) uv = 1.d0
              ux(i,j,k,1) = uv
            enddo
          enddo
        enddo
c  
        if(iprocz.eq.1) then
          do i=1, nxlocal+1, 1
              bc_uwall(i,nylocal+1,2) = 0.d0
          enddo
          hnup = bc_hnuplus 
          do i=1, nxlocal+1, 1
            do j=1, nylocal, 1
              zv= h0
              zp= zv*utau/nu
              if (zp .gt. hnup ) then
                usbc= utau*(1.d0/vk*dlog( zp/hnup)+hnup)
c                write(*,*) zp, zv, usbc, ux(i,j,1,1),ux(i,j,2,1)
              else    
                usbc= utau* zp
              endif 
              bc_uwall(i,j,1)= usbc
              bc_uwall(i,j,2)=0.d0
              bc_uwall(i,j,3)=0.d0
              bc_utau(i,j)=utau
              bc_dudz(i,j) = utau*utau*Reynolds
            enddo
          enddo
          do i= -nghost+1, 0, 1
            do j=1, nylocal, 1
              bc_uwall(i,j,1) = usbc
              bc_uwall(i,j,2:3) = 0.d0
            enddo
          enddo
        endif
c
        call setbcbcuwall
#ifdef PARALLEL
        call exchangebcuwall
#endif
c
        Uin(1,iylo:iyhi,izlo:izhi,1) = 
     &   ux(1,iylo:iyhi,izlo:izhi,1)
        Uin(1,iylo:iyhi,izlo:izhi,3) = 
     &   ux(1,iylo:iyhi,izlo:izhi,3)
        call SetBoundaryValuesST(ux)
        call exchange(ux,3)

c
        call addperturbation(ux)
c
        return
        end
c
c----------------------------------------------------------------------
