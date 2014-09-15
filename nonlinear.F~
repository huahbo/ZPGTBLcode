c-----------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2011
c     King Abdullah University of Science and Technology
c     All Rights Reserved
c-----------------------------------------------------------------
c     $Log: nonlinear.F,v $
c     Revision 1.4  2011/10/19 15:14:49  samtanr
c     Added 2D option via TWO_D defines.
c
c     Revision 1.3  2011/10/03 08:35:03  chengw
c     remove commented lines.
c
c     Revision 1.2  2011/10/03 08:18:29  samtanr
c     Intermediate checkin.
c
c     Revision 1.1  2011/07/25 11:32:01  samtanr
c     Original source.
c     Wind code: Incompressible Navier-Stokes LES
c
c-----------------------------------------------------------------
       subroutine ComputeAdvectionTerms(Adv,ux)
c-----------------------------------------------------------------
       use mesh
       use properties
#ifdef PARALLEL
       use mpistuff
#endif

c
       implicit none
c=======Declarations=========
       double precision:: ux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c       double precision:: dt
c
#ifndef TWO_D
       double precision:: Adv(nxlocal+1,nylocal+1,nzlocal+1,nvar)
       double precision:: adv_t(nxlocal+1,nylocal+1,nzlocal+1)
#else
       double precision:: Adv(nxlocal+1,nylocal,nzlocal+1,nvar)
       double precision:: adv_t(nxlocal+1,nylocal,nzlocal+1)
#endif
c
c     
c     Compute terms using skew symmetric form and 
c     using finite differences in x,z and Fourier in y
c
c     Define velocity={u,v,w}
c     u is defined on i+1/2,j,k
c     v is defined on i, j+1/2,k
c     w is defined on i, j, k+1/2
c
c     cell center i,j is bounded by (i-1/2 j-1/2) = (i,j), 
c                                   (i+1/2,j-1/2) = (i+1,j)
c                                   (i+1/2,j+1/2) = (i+1,j+1),
c                                   (i-1/2,j+1/2) = (i,j+1) 
c

       adv_t(:,:,:) = Adv(:,:,:,1)
       call XMomentum(adv_t,ux)
       Adv(:,:,:,1) = adv_t(:,:,:)
#ifndef TWO_D
       adv_t(:,:,:) = Adv(:,:,:,2)
       call YMomentum(adv_t,ux)
       Adv(:,:,:,2) = adv_t(:,:,:)
#else
       adv(:,:,:,2)=0.D0
#endif       
       adv_t(:,:,:) = Adv(:,:,:,3)
       call ZMomentum(adv_t,ux)
       Adv(:,:,:,3) = adv_t(:,:,:)
c
       return
       end subroutine ComputeAdvectionTerms
c
c-----------------------------------------------------------------
       subroutine XMomentum(Adv,ux)
c-----------------------------------------------------------------
       use mesh
       use properties
#ifdef PARALLEL
       use mpistuff
#endif

c
       implicit none
c=======Declarations=========
       double precision:: ux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c       double precision:: dt
c
c       double precision:: Adv(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI)
#ifndef TWO_D
       double precision:: Adv(nxlocal+1,nylocal+1,nzlocal+1)
#else
       double precision:: Adv(nxlocal+1,nylocal,nzlocal+1)
#endif
c
c     
       double precision:: aa,bb
       double precision:: tmp,uc1,uc3, v1,v3,w1,w3,u1,u3
       double precision:: tmp1(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI)
       double precision:: tmp3(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI)
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
c     X - MOMENTUM
c     Need to solve x-momentum from i=1 to nx+1
c     div = d(uu)/dx + d(uv)/dy + d(uw)/dz
c     adv = u du/dx + vdu/dy + wdu/dz
c
c     consider d(uu)/dx term - these are first evaluated at
c     cell centers and then interpolated to i+1/2 faces
       do i=-1,nxlocal+2,1
          do j=1,nylocal,1
             do k=1,nzlocal,1
                uc1=0.5D0*(ux(i,j,k,1)+ux(i+1,j,k,1))
                uc3=0.5D0*(ux(i+2,j,k,1)+ux(i-1,j,k,1))
                tmp=aa*uc1+bb*uc3
                tmp1(i,j,k)=tmp*uc1
                tmp3(i,j,k)=tmp*uc3
             enddo
          enddo
       enddo
       do i=1,nxlocal+1,1
          do j=1,nylocal,1
             do k=1,nzlocal,1
                adv(i,j,k)=aa*(tmp1(i,j,k)-tmp1(i-1,j,k))*dxi+
     &               bb*(tmp3(i+1,j,k)-tmp3(i-2,j,k))*dxi3
             enddo
          enddo
       enddo
c     Advection term udu/dx
       do i=-1,nxlocal+2,1
          do j=1,nylocal,1
             do k=1,nzlocal,1
                uc1=0.5D0*(ux(i,j,k,1)+ux(i+1,j,k,1))
                uc3=0.5D0*(ux(i+2,j,k,1)+ux(i-1,j,k,1))
                tmp=aa*uc1+bb*uc3
                tmp1(i,j,k)=tmp*(ux(i+1,j,k,1)-ux(i,j,k,1))*dxi
                tmp3(i,j,k)=tmp*(ux(i+2,j,k,1)-ux(i-1,j,k,1))*dxi3
             enddo
          enddo
       enddo
       do i=1,nxlocal+1,1
          do j=1,nylocal,1
             do k=1,nzlocal,1
                adv(i,j,k)=0.5D0*adv(i,j,k)+
     &               0.25D0*(aa*(tmp1(i,j,k)+tmp1(i-1,j,k))
     &                      +bb*(tmp3(i+1,j,k)+tmp3(i-2,j,k)))
             enddo
          enddo
       enddo
c
#ifndef TWO_D
c     Second term d(uv)/dy: div part
c
       do i=1,nxlocal+1,1
          do j=0,nylocal+2,1
             do k=1,nzlocal,1
                u1=0.5D0*(ux(i,j-1,k,1)+ux(i,j,k,1))
                u3=0.5D0*(ux(i,j-2,k,1)+ux(i,j+1,k,1))
                v1=0.5D0*(ux(i-1,j,k,2)+ux(i,j,k,2))
                v3=0.5D0*(ux(i-2,j,k,2)+ux(i+1,j,k,2))
                tmp=aa*v1+bb*v3
                tmp1(i,j,k)=tmp*u1
                tmp3(i,j,k)=tmp*u3
             enddo
          enddo
       enddo
c
       do i=1,nxlocal+1,1
          do j=1,nylocal,1
             do k=1,nzlocal,1
                adv(i,j,k)=adv(i,j,k)+0.5D0*
     &               (aa*(tmp1(i,j+1,k)-tmp1(i,j,k))*dyi
     &               +bb*(tmp3(i,j+2,k)-tmp3(i,j-1,k))*dyi3)
             enddo
          enddo
       enddo
c
c     Second term vd(u)/dy: adv part
c
       do i=1,nxlocal+1,1
          do j=0,nylocal+2,1
             do k=1,nzlocal,1
                v1=0.5D0*(ux(i-1,j,k,2)+ux(i,j,k,2))
                v3=0.5D0*(ux(i-2,j,k,2)+ux(i+1,j,k,2))
                tmp=aa*v1+bb*v3
                tmp1(i,j,k)=tmp*(ux(i,j,k,1)-ux(i,j-1,k,1))*dyi
                tmp3(i,j,k)=tmp*(ux(i,j+1,k,1)-ux(i,j-2,k,1))*dyi3
             enddo
          enddo
       enddo
c
       do i=1,nxlocal+1,1
          do j=1,nylocal,1
             do k=1,nzlocal,1
                adv(i,j,k)=adv(i,j,k)+0.25D0*
     &               (aa*(tmp1(i,j+1,k)+tmp1(i,j,k))
     &               +bb*(tmp3(i,j+2,k)+tmp3(i,j-1,k)))
             enddo
          enddo
       enddo
c
#endif
c     Third term d(uw)/dz: div part
c
       do i=1,nxlocal+1,1
          do j=1,nylocal,1
             do k=0,nzlocal+2,1
                u1=0.5D0*(ux(i,j,k-1,1)+ux(i,j,k,1))
                u3=0.5D0*(ux(i,j,k-2,1)+ux(i,j,k+1,1))
                w1=0.5D0*(ux(i-1,j,k,3)+ux(i,j,k,3))
                w3=0.5D0*(ux(i-2,j,k,3)+ux(i+1,j,k,3))
                tmp=aa*w1+bb*w3
                tmp1(i,j,k)=tmp*u1
                tmp3(i,j,k)=tmp*u3
             enddo
          enddo
       enddo
c
       do i=1,nxlocal+1,1
          do j=1,nylocal,1
             do k=1,nzlocal,1
                adv(i,j,k)=adv(i,j,k)+0.5D0*
     &               (aa*(tmp1(i,j,k+1)-tmp1(i,j,k))*dzi
     &               +bb*(tmp3(i,j,k+2)-tmp3(i,j,k-1))*dzi3)
             enddo
          enddo
       enddo
c
c     Third term wd(u)/dz: adv part
c
       do i=1,nxlocal+1,1
          do j=1,nylocal,1
             do k=0,nzlocal+2,1
                w1=0.5D0*(ux(i-1,j,k,3)+ux(i,j,k,3))
                w3=0.5D0*(ux(i-2,j,k,3)+ux(i+1,j,k,3))
                tmp=aa*w1+bb*w3
                tmp1(i,j,k)=tmp*(ux(i,j,k,1)-ux(i,j,k-1,1))*dzi
                tmp3(i,j,k)=tmp*(ux(i,j,k+1,1)-ux(i,j,k-2,1))*dzi3
             enddo
          enddo
       enddo
c
       do i=1,nxlocal+1,1
          do j=1,nylocal,1
             do k=1,nzlocal,1
                adv(i,j,k)=adv(i,j,k)+0.25D0*
     &               (aa*(tmp1(i,j,k+1)+tmp1(i,j,k))
     &               +bb*(tmp3(i,j,k+2)+tmp3(i,j,k-1)))
             enddo
          enddo
       enddo
c
c
       return
       end
c     
c-----------------------------------------------------------------
       subroutine YMomentum(Adv,ux)
c-----------------------------------------------------------------
       use mesh
       use properties
#ifdef PARALLEL
       use mpistuff
#endif

c
       implicit none
c=======Declarations=========
       double precision:: ux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c       double precision:: dt
c
c       double precision:: Adv(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI)
       double precision:: Adv(nxlocal+1,nylocal+1,nzlocal+1)
c
c     
       double precision:: aa,bb
       double precision:: tmp,vc1,vc3, v1,v3,w1,w3,u1,u3
       double precision:: tmp1(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI)
       double precision:: tmp3(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI)
       integer:: i,j,k
c     
c     u is defined on i+1/2,j,k
c     v is defined on i, j+1/2,k
c     w is defined on i, j, k+1/2
c
#ifndef TWO_D
       aa=9.D0/8.D0
       bb=-1.D0/8.D0
c     Y - MOMENTUM
c     Need to solve y-momentum from j=1 to ny+1
c     div = d(uv)/dx + d(vv)/dy + d(vw)/dz
c     adv = u dv/dx + vdv/dy + wdv/dz
c
c     consider d(vv)/dy term - these are first evaluated at
c     cell centers and then interpolated to i+1/2 faces
c     Second term d(vv)/dy: div part
       do i=1,nxlocal,1
          do j=-1,nylocal+2,1
             do k=1,nzlocal,1
                vc1=0.5D0*(ux(i,j,k,2)+ux(i,j+1,k,2))
                vc3=0.5D0*(ux(i,j+2,k,2)+ux(i,j-1,k,2))
                tmp=aa*vc1+bb*vc3
                tmp1(i,j,k)=tmp*vc1
                tmp3(i,j,k)=tmp*vc3
             enddo
          enddo
       enddo
       do i=1,nxlocal,1
          do j=1,nylocal+1,1
             do k=1,nzlocal,1
                adv(i,j,k)=aa*(tmp1(i,j,k)-tmp1(i,j-1,k))*dyi+
     &               bb*(tmp3(i,j+1,k)-tmp3(i,j-2,k))*dyi3
             enddo
          enddo
       enddo
c     Second term vdv/dy: adv part
       do i=1,nxlocal,1
          do j=-1,nylocal+2,1
             do k=1,nzlocal,1
                vc1=0.5D0*(ux(i,j,k,2)+ux(i,j+1,k,2))
                vc3=0.5D0*(ux(i,j+2,k,2)+ux(i,j-1,k,2))
                tmp=aa*vc1+bb*vc3
                tmp1(i,j,k)=tmp*(ux(i,j+1,k,2)-ux(i,j,k,2))*dyi
                tmp3(i,j,k)=tmp*(ux(i,j+2,k,2)-ux(i,j-1,k,2))*dyi3
             enddo
          enddo
       enddo
       do i=1,nxlocal,1
          do j=1,nylocal+1,1
             do k=1,nzlocal,1
                adv(i,j,k)=0.5D0*adv(i,j,k)+
     &               0.25D0*(aa*(tmp1(i,j,k)+tmp1(i,j-1,k))
     &                      +bb*(tmp3(i,j+1,k)+tmp3(i,j-2,k)))
             enddo
          enddo
       enddo
c
c     First term d(uv)/dx: div part
c
       do i=0,nxlocal+2,1
          do j=1,nylocal+1,1
             do k=1,nzlocal,1
                u1=0.5D0*(ux(i,j-1,k,1)+ux(i,j,k,1))
                u3=0.5D0*(ux(i,j-2,k,1)+ux(i,j+1,k,1))
                v1=0.5D0*(ux(i-1,j,k,2)+ux(i,j,k,2))
                v3=0.5D0*(ux(i-2,j,k,2)+ux(i+1,j,k,2))
                tmp=aa*u1+bb*u3
                tmp1(i,j,k)=tmp*v1
                tmp3(i,j,k)=tmp*v3
             enddo
          enddo
       enddo
c
       do i=1,nxlocal,1
          do j=1,nylocal+1,1
             do k=1,nzlocal,1
                adv(i,j,k)=adv(i,j,k)+0.5D0*
     &               (aa*(tmp1(i+1,j,k)-tmp1(i,j,k))*dxi
     &               +bb*(tmp3(i+2,j,k)-tmp3(i-1,j,k))*dxi3)
             enddo
          enddo
       enddo
c
c     Second term vd(u)/dy: adv part
c
       do i=0,nxlocal+2,1
          do j=1,nylocal+1,1
             do k=1,nzlocal,1
                u1=0.5D0*(ux(i,j-1,k,1)+ux(i,j,k,1))
                u3=0.5D0*(ux(i,j-2,k,1)+ux(i,j+1,k,1))
                tmp=aa*u1+bb*u3
                tmp1(i,j,k)=tmp*(ux(i,j,k,2)-ux(i-1,j,k,2))*dxi
                tmp3(i,j,k)=tmp*(ux(i+1,j,k,2)-ux(i-2,j,k,2))*dxi3
             enddo
          enddo
       enddo
c
       do i=1,nxlocal,1
          do j=1,nylocal+1,1
             do k=1,nzlocal,1
                adv(i,j,k)=adv(i,j,k)+0.25D0*
     &               (aa*(tmp1(i+1,j,k)+tmp1(i,j,k))
     &               +bb*(tmp3(i+2,j,k)+tmp3(i-1,j,k)))
             enddo
          enddo
       enddo
c
c     Third term d(vw)/dz: div part
c
       do i=1,nxlocal,1
          do j=1,nylocal+1,1
             do k=0,nzlocal+2,1
                v1=0.5D0*(ux(i,j,k-1,2)+ux(i,j,k,2))
                v3=0.5D0*(ux(i,j,k-2,2)+ux(i,j,k+1,2))
                w1=0.5D0*(ux(i,j-1,k,3)+ux(i,j,k,3))
                w3=0.5D0*(ux(i,j-2,k,3)+ux(i,j+1,k,3))
                tmp=aa*w1+bb*w3
                tmp1(i,j,k)=tmp*v1
                tmp3(i,j,k)=tmp*v3
             enddo
          enddo
       enddo
c
       do i=1,nxlocal,1
          do j=1,nylocal+1,1
             do k=1,nzlocal,1
                adv(i,j,k)=adv(i,j,k)+0.5D0*
     &               (aa*(tmp1(i,j,k+1)-tmp1(i,j,k))*dzi
     &               +bb*(tmp3(i,j,k+2)-tmp3(i,j,k-1))*dzi3)
             enddo
          enddo
       enddo
c
c     Third term wd(v)/dz: adv part
c
       do i=1,nxlocal,1
          do j=1,nylocal+1,1
             do k=0,nzlocal+2,1
                w1=0.5D0*(ux(i,j-1,k,3)+ux(i,j,k,3))
                w3=0.5D0*(ux(i,j-2,k,3)+ux(i,j+1,k,3))
                tmp=aa*w1+bb*w3
                tmp1(i,j,k)=tmp*(ux(i,j,k,2)-ux(i,j,k-1,2))*dzi
                tmp3(i,j,k)=tmp*(ux(i,j,k+1,2)-ux(i,j,k-2,2))*dzi3
             enddo
          enddo
       enddo
c
       do i=1,nxlocal,1
          do j=1,nylocal+1,1
             do k=1,nzlocal,1
                adv(i,j,k)=adv(i,j,k)+0.25D0*
     &               (aa*(tmp1(i,j,k+1)+tmp1(i,j,k))
     &               +bb*(tmp3(i,j,k+2)+tmp3(i,j,k-1)))
             enddo
          enddo
       enddo
#endif
       return
       end subroutine YMomentum
c
cc-----------------------------------------------------------------
       subroutine ZMomentum(Adv,ux)
c-----------------------------------------------------------------
       use mesh
       use properties
#ifdef PARALLEL
       use mpistuff
#endif

c
       implicit none
c=======Declarations=========
       double precision:: ux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c       double precision:: dt
c
c       double precision:: Adv(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI)
#ifndef TWO_D
       double precision:: Adv(nxlocal+1,nylocal+1,nzlocal+1)
#else
       double precision:: Adv(nxlocal+1,nylocal,nzlocal+1)
#endif
c
c     
       double precision:: aa,bb
       double precision:: tmp,wc1,wc3, v1,v3,w1,w3,u1,u3
       double precision:: tmp1(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI)
       double precision:: tmp3(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI)
       integer:: i,j,k
c     
c     u is defined on i+1/2,j,k
c     v is defined on i, j+1/2,k
c     w is defined on i, j, k+1/2
c
       aa=9.D0/8.D0
       bb=-1.D0/8.D0
c
c     Z - MOMENTUM
c     Need to solve z-momentum from k=1 to nz+1
c     div = d(uw)/dx + d(wv)/dy + d(ww)/dz
c     adv = u dw/dx + vdw/dy + wdw/dz
c
c     consider d(ww)/dz term - these are first evaluated at
c     cell centers and then interpolated to i+1/2 faces
       do i=1,nxlocal,1
          do j=1,nylocal,1
             do k=-1,nzlocal+2,1
                wc1=0.5D0*(ux(i,j,k,3)+ux(i,j,k+1,3))
                wc3=0.5D0*(ux(i,j,k+2,3)+ux(i,j,k-1,3))
                tmp=aa*wc1+bb*wc3
                tmp1(i,j,k)=tmp*wc1
                tmp3(i,j,k)=tmp*wc3
             enddo
          enddo
       enddo
       do i=1,nxlocal,1
          do j=1,nylocal,1
             do k=1,nzlocal+1,1
                adv(i,j,k)=aa*(tmp1(i,j,k)-tmp1(i,j,k-1))*dzi+
     &               bb*(tmp3(i,j,k+1)-tmp3(i,j,k-2))*dzi3
             enddo
          enddo
       enddo
c     Advection term wdw/dz
       do i=1,nxlocal,1
          do j=1,nylocal,1
             do k=-1,nzlocal+2,1
                wc1=0.5D0*(ux(i,j,k,3)+ux(i,j,k+1,3))
                wc3=0.5D0*(ux(i,j,k+2,3)+ux(i,j,k-1,3))
                tmp=aa*wc1+bb*wc3
                tmp1(i,j,k)=tmp*(ux(i,j,k+1,3)-ux(i,j,k,3))*dzi
                tmp3(i,j,k)=tmp*(ux(i,j,k+2,3)-ux(i,j,k-1,3))*dzi3
             enddo
          enddo
       enddo
       do i=1,nxlocal,1
          do j=1,nylocal,1
             do k=1,nzlocal+1,1
                adv(i,j,k)=0.5D0*adv(i,j,k)+
     &               0.25D0*(aa*(tmp1(i,j,k)+tmp1(i,j,k-1))
     &                      +bb*(tmp3(i,j,k+1)+tmp3(i,j,k-2)))
             enddo
          enddo
       enddo
c
#ifndef TWO_D

c     Second term d(vw)/dy: div part
c
       do i=1,nxlocal,1
          do j=0,nylocal+2,1
             do k=1,nzlocal+1,1
                w1=0.5D0*(ux(i,j-1,k,3)+ux(i,j,k,3))
                w3=0.5D0*(ux(i,j-2,k,3)+ux(i,j+1,k,3))
                v1=0.5D0*(ux(i,j,k-1,2)+ux(i,j,k,2))
                v3=0.5D0*(ux(i,j,k-2,2)+ux(i,j,k+1,2))
                tmp=aa*v1+bb*v3
                tmp1(i,j,k)=tmp*w1
                tmp3(i,j,k)=tmp*w3
             enddo
          enddo
       enddo
c
       do i=1,nxlocal,1
          do j=1,nylocal,1
             do k=1,nzlocal+1,1
                adv(i,j,k)=adv(i,j,k)+0.5D0*
     &               (aa*(tmp1(i,j+1,k)-tmp1(i,j,k))*dyi
     &               +bb*(tmp3(i,j+2,k)-tmp3(i,j-1,k))*dyi3)
             enddo
          enddo
       enddo
c
c     Second term vd(w)/dy: adv part
c
       do i=1,nxlocal,1
          do j=0,nylocal+2,1
             do k=1,nzlocal+1,1
                v1=0.5D0*(ux(i,j,k-1,2)+ux(i,j,k,2))
                v3=0.5D0*(ux(i,j,k-2,2)+ux(i,j,k+1,2))
                tmp=aa*v1+bb*v3
                tmp1(i,j,k)=tmp*(ux(i,j,k,3)-ux(i,j-1,k,3))*dyi
                tmp3(i,j,k)=tmp*(ux(i,j+1,k,3)-ux(i,j-2,k,3))*dyi3
             enddo
          enddo
       enddo
c
       do i=1,nxlocal,1
          do j=1,nylocal,1
             do k=1,nzlocal+1,1
                adv(i,j,k)=adv(i,j,k)+0.25D0*
     &               (aa*(tmp1(i,j+1,k)+tmp1(i,j,k))
     &               +bb*(tmp3(i,j+2,k)+tmp3(i,j-1,k)))
             enddo
          enddo
       enddo
c
#endif
c     First term d(uw)/dx: div part
c
       do i=0,nxlocal+2,1
          do j=1,nylocal,1
             do k=1,nzlocal+1,1
                u1=0.5D0*(ux(i,j,k-1,1)+ux(i,j,k,1))
                u3=0.5D0*(ux(i,j,k-2,1)+ux(i,j,k+1,1))
                w1=0.5D0*(ux(i-1,j,k,3)+ux(i,j,k,3))
                w3=0.5D0*(ux(i-2,j,k,3)+ux(i+1,j,k,3))
                tmp=aa*u1+bb*u3
                tmp1(i,j,k)=tmp*w1
                tmp3(i,j,k)=tmp*w3
             enddo
          enddo
       enddo
c
       do i=1,nxlocal,1
          do j=1,nylocal,1
             do k=1,nzlocal+1,1
                adv(i,j,k)=adv(i,j,k)+0.5D0*
     &               (aa*(tmp1(i+1,j,k)-tmp1(i,j,k))*dxi
     &               +bb*(tmp3(i+2,j,k)-tmp3(i-1,j,k))*dxi3)
             enddo
          enddo
       enddo
c
c     First term ud(w)/dx: adv part
c
       do i=0,nxlocal+2,1
          do j=1,nylocal,1
             do k=1,nzlocal+1,1
                u1=0.5D0*(ux(i,j,k-1,1)+ux(i,j,k,1))
                u3=0.5D0*(ux(i,j,k-2,1)+ux(i,j,k+1,1))
                tmp=aa*u1+bb*u3
                tmp1(i,j,k)=tmp*(ux(i,j,k,3)-ux(i-1,j,k,3))*dxi
                tmp3(i,j,k)=tmp*(ux(i+1,j,k,3)-ux(i-2,j,k,3))*dxi3
             enddo
          enddo
       enddo
c
       do i=1,nxlocal,1
          do j=1,nylocal,1
             do k=1,nzlocal+1,1
                adv(i,j,k)=adv(i,j,k)+0.25D0*
     &               (aa*(tmp1(i+1,j,k)+tmp1(i,j,k))
     &               +bb*(tmp3(i+2,j,k)+tmp3(i-1,j,k)))
             enddo
          enddo
       enddo
c
       end subroutine ZMomentum
