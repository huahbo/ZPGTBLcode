
c----------------------------------------------------------------------
c      this subroutine is for the derivative of the sgs term 
c                             in momentum equation
c----------------------------------------------------------------------
       subroutine sgsnonlinear(adv,uu,timestep)
       use mesh
       use mesh_common
       use properties
#ifdef PARALLEL
       use mpistuff
#endif
c
       implicit none
c=======Declarations=========
       double precision:: uu(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
       double precision:: Adv(nxlocal+1,nylocal+1,nzlocal+1,nvar)
       double precision:: Advt(nxlocal+1,nylocal+1,nzlocal+1,nvar)
       double precision:: Tij(-1:nxlocal+2,-1:nylocal+2,-1:nzlocal+2,6)
       double precision:: Tij_cn(0:nxlocal+2,0:nylocal+2,0:nzlocal+2,6)
       double precision:: tmp1(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
       double precision:: aa,bb, cc, dd
       integer:: i,j,k, l, timestep, izbo
c
       Advt = 0.d0

       aa = 9.D0/8.D0
       bb =-1.D0/8.D0
       cc = 8.d0/12.d0 
       dd =-1.d0/12.d0
c      cc and dd for central scheme:
c        dudx_0 = cc* (u_1-u_-1)/dx - dd*(u2-u_-2 )/dx
c            = 4/3 * dudx^2_0 - 1/3 *dudx^4_0
c         u_0 = 4/3 * u^2_0 - 1/3 u^4_0 
c
c      getting sgs tensor at center and corner.
c
       call sgsmodel(uu, Tij, Tij_cn,timestep)
c       write(*,*)'Max/Min Tij at corner:',maxval(Tij_cn),minval(Tij_cn)
c     &    ,iprocx, iprocy, iprocz
c
c       call eigval_Tijcn(Tij_cn,timestep)
c       call WriteHDF5FileTijcn(Tij_cn,timestep)
c
c      calculating  dTijdxj
c               -- 1  2  3 --   
c        Tij =  |  2  4  5  |
c               -- 3  5  6 --

       tmp1 = 1.d7

c      x-momentum equation        
       do i=1, nxlocal+1, 1
         do j=0, nylocal+2, 1
           do k=1, nzlocal, 1
           tmp1(i,j,k,2)=aa*(Tij_cn(i,j,k  ,2)+Tij_cn(i,j,k+1,2))/2.d0+ 
     &                   bb*(Tij_cn(i,j,k-1,2)+Tij_cn(i,j,k+2,2))/2.d0
           enddo
         enddo
         do j=1, nylocal, 1
           do k= 0, nzlocal+2, 1
           tmp1(i,j,k,3)=aa*(Tij_cn(i,j  ,k,3)+Tij_cn(i,j+1,k,3))/2.d0+ 
     &                   bb*(Tij_cn(i,j-1,k,3)+Tij_cn(i,j+2,k,3))/2.d0
           enddo
         enddo
       enddo
       do i=1, nxlocal+1, 1
         do j=1, nylocal, 1
           do k=1, nzlocal, 1
c           Adv(i,j,k,1) =  Adv(i,j,k,1) + 
            Advt(i,j,k,1) =
     &                aa*( Tij(i  ,j,k,1)- Tij(i-1,j,k,1))*dxi + 
     &                bb*( Tij(i+1,j,k,1)- Tij(i-2,j,k,1))*dxi3+
     &                aa*(tmp1(i,j+1,k,2)-tmp1(i,j  ,k,2))*dyi + 
     &                bb*(tmp1(i,j+2,k,2)-tmp1(i,j-1,k,2))*dyi3+
     &                aa*(tmp1(i,j,k+1,3)-tmp1(i,j,k  ,3))*dzi + 
     &                bb*(tmp1(i,j,k+2,3)-tmp1(i,j,k-1,3))*dzi3 
           enddo
         enddo
       enddo
c
c      y-momentum equation
       do j=1, nylocal+1, 1
         do i= 0,nxlocal+2, 1
           do k= 1, nzlocal, 1
           tmp1(i,j,k,1)=aa*(Tij_cn(i,j,k  ,2)+Tij_cn(i,j,k+1,2))/2.d0+ 
     &                   bb*(Tij_cn(i,j,k-1,2)+Tij_cn(i,j,k+2,2))/2.d0
           enddo
         enddo
         do i=1, nxlocal, 1
           do k= 0, nzlocal+2, 1
           tmp1(i,j,k,3)=aa*(Tij_cn(i  ,j,k,5)+Tij_cn(i+1,j,k,5))/2.d0+ 
     &                   bb*(Tij_cn(i-1,j,k,5)+Tij_cn(i+2,j,k,5))/2.d0
           enddo
         enddo
       enddo
       do i=1, nxlocal, 1
         do j=1, nylocal+1, 1
           do k=1, nzlocal, 1
c           Adv(i,j,k,2) = Adv(i,j,k,2) + 
            Advt(i,j,k,2) = 
     &                aa*(tmp1(i+1,j,k,1)-tmp1(i  ,j,k,1))*dxi + 
     &                bb*(tmp1(i+2,j,k,1)-tmp1(i-1,j,k,1))*dxi3+
     &                aa*( Tij(i,j  ,k,4)- Tij(i,j-1,k,4))*dyi + 
     &                bb*( Tij(i,j+1,k,4)- Tij(i,j-2,k,4))*dyi3+
     &                aa*(tmp1(i,j,k+1,3)-tmp1(i,j,k  ,3))*dzi + 
     &                bb*(tmp1(i,j,k+2,3)-tmp1(i,j,k-1,3))*dzi3 
           enddo
         enddo
       enddo
c
c      z-momentum equation
       do k= 1, nzlocal+1, 1
         do i= 0, nxlocal+2, 1
           do j= 1, nylocal, 1
           tmp1(i,j,k,1)=aa*(Tij_cn(i,j  ,k,3)+Tij_cn(i,j+1,k,3))/2.d0+ 
     &                   bb*(Tij_cn(i,j-1,k,3)+Tij_cn(i,j+2,k,3))/2.d0
           enddo
         enddo 
         do i=1, nxlocal, 1
           do j=0, nylocal+2, 1
           tmp1(i,j,k,2)=aa*(Tij_cn(i  ,j,k,5)+Tij_cn(i+1,j,k,5))/2.d0+ 
     &                   bb*(Tij_cn(i-1,j,k,5)+Tij_cn(i+2,j,k,5))/2.d0
           enddo
         enddo
       enddo
       izbo = 1
       if(iprocz.eq.1) izbo = 2
       do i=1, nxlocal, 1
         do j=1, nylocal, 1
           do k=izbo, nzlocal+1, 1
c           Adv(i,j,k,3) = Adv(i,j,k,3)  + 
            Advt(i,j,k,3) =
     &                aa*(tmp1(i+1,j,k,1)-tmp1(i  ,j,k,1))*dxi + 
     &                bb*(tmp1(i+2,j,k,1)-tmp1(i-1,j,k,1))*dxi3+
     &                aa*(tmp1(i,j+1,k,2)-tmp1(i,j  ,k,2))*dyi + 
     &                bb*(tmp1(i,j+2,k,2)-tmp1(i,j-1,k,2))*dyi3+
     &                aa*( Tij(i,j,k  ,6)- Tij(i,j,k-1,6))*dzi + 
     &                bb*( Tij(i,j,k+1,6)- Tij(i,j,k-2,6))*dzi3 
           enddo
         enddo
       enddo      
c      the z component does not include the bottom point
c
       call fftnonlinear(Advt)
       Adv = Adv + Advt




       return
       end
c
c----------------------------------------------------------------------
