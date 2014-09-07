      subroutine SetBoundaryValuesST(u)
      use mesh_parms
      use mesh_common
      use InflowVelBC
      use outflowBC
#ifdef LES
      use virtualwallBC
#endif
      double precision:: u(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision:: utmp(IXLO:IXHI,IYLO:IYHI,nvar)
c      double precision:: vt
      integer:: i,j,k,l
      integer:: ilo,ihi,jlo,jhi
c
c  In doing setboundary, due to the exist of virtual wall, 
c   the z direction is set before setting y periodic.
c
c     step 1:
c     x component: set iprocx =1 and iprocx = xprocs.
c
c     xl Boundary (reflecting)
      if (iprocx .eq. 1) then
         i=1
         do k=IZLO,IZHI,1
            do j=IYLO,IYHI,1
               u(i,j,k,1)=Uin(1,j,k,1)
               u(i-1,j,k,1)= Uin(2,j,k,1)
               u(i-2,j,k,1)= Uin(3,j,k,1)
               u(i-3,j,k,1)= Uin(4,j,k,1)
               u(i-4,j,k,1)= Uin(5,j,k,1)
c               u(i-1,j,k,1)=2.D0*Uin(1,j,k,1)-u(i+1,j,k,1)
c               u(i-2,j,k,1)=2.D0*Uin(1,j,k,1)-u(i+2,j,k,1)
c               u(i-3,j,k,1)=2.D0*Uin(1,j,k,1)-u(i+3,j,k,1)
c               u(i-4,j,k,1)=2.d0*Uin(1,j,k,1)-u(i+4,j,k,1)
            enddo
         enddo
c     
         i=0
         do k=IZLO,IZHI,1
            do j=IYLO,IYHI,1
               do l=0,nghost-1,1
                  u(i-l,j,k,2:3)= Uin(2+l,j,k,2:3)
c                  u(i-l,j,k,2:3)= 2.d0*Uin(1,j,k,2:3)-u(i+1+l,j,k,2:3)
               enddo
c     
            enddo
        enddo
      endif
#ifdef CODEB
      if (iprocx .eq. 1) then
         i=1
         do k=IZLO,IZHI,1
            do j=IYLO,IYHI,1
               u(i,j,k,1)  = Uin(1,j,k,1)
               u(i-1,j,k,1)= Uin(2,j,k,1)
               u(i-2,j,k,1)= Uin(3,j,k,1)
               u(i-3,j,k,1)= Uin(4,j,k,1)
               u(i-4,j,k,1)= Uin(5,j,k,1)
            enddo
         enddo
c     
         i=0
         do k=IZLO,IZHI,1
            do j=IYLO,IYHI,1
               do l=0,nghost-1,1
                  u(i-l,j,k,2:3)= Uin(2+l,j,k,2:3)
               enddo
c     
            enddo
        enddo
      endif
#endif
c     
c     xr Boundary: ZeroGradient
      if (iprocx .eq. xprocs) then
         i=nxlocal+1
         do k=IZLO,IZHI,1
            do j=IYLO,IYHI,1
#ifndef OUTFBC
               u(i+1,j,k,1)=2.D0*u(i,j,k,1)-u(i-1,j,k,1)
c               u(i+2,j,k,1)=2.D0*u(i+1,j,k,1)-u(i,j,k,1)
c               u(i+3,j,k,1)=2.D0*u(i+2,j,k,1)-u(i+1,j,k,1)
               u(i+2,j,k,1)=2.D0*u(i+1,j,k,1)-u(i-2,j,k,1)
               u(i+3,j,k,1)=2.D0*u(i+2,j,k,1)-u(i-3,j,k,1)

#else
               u(i+1,j,k,1)=bc_outflowu(1,j,k,1)
               u(i+2,j,k,1)=bc_outflowu(2,j,k,1)
               u(i+3,j,k,1)=bc_outflowu(3,j,k,1)
#endif
            enddo
         enddo   
c     
         i = nxlocal
         do k=IZLO,IZHI,1
            do j=IYLO,IYHI,1
#ifndef OUTFBC
c               u(i+1,j,k,2:3)=2.D0*u(i,j,k,2:3)-u(i-1,j,k,2:3)
c               u(i+2,j,k,2:3)=2.D0*u(i+1,j,k,2:3)-u(i,j,k,2:3)
c               u(i+3,j,k,2:3)=2.D0*u(i+2,j,k,2:3)-u(i+1,j,k,2:3)
c               u(i+4,j,k,2:3)=2.d0*u(i+3,j,k,2:3)-u(i+2,j,k,2:3)
               u(i+1,j,k,2:3)= u(i,j,k,2:3)
               u(i+2,j,k,2:3)= u(i-1,j,k,2:3)
               u(i+3,j,k,2:3)= u(i-2,j,k,2:3)
               u(i+4,j,k,2:3)= u(i-3,j,k,2:3)
#else
               u(i+1,j,k,2:3)=bc_outflowu(1,j,k,2:3)
               u(i+2,j,k,2:3)=bc_outflowu(2,j,k,2:3)
               u(i+3,j,k,2:3)=bc_outflowu(3,j,k,2:3)
               u(i+4,j,k,2:3)=bc_outflowu(4,j,k,2:3)
#endif
            enddo
         enddo
      endif
c
c     step 2:
c     z component: set iprocz =1 and iprocz = zprocs.
c
c     zl Boundary: Reflecting

c        write(*,*) maxval(utmp(:,:,1)),minval(utmp(:,:,1)),
c     &    maxval(utmp(:,:,2)),minval(utmp(:,:,2)),
c     &    maxval(utmp(:,:,3)),minval(utmp(:,:,3))

      if (iprocz .eq. 1) then
c
#ifdef LES
       utmp(ixlo:ixhi,iylo:iyhi,1:3)= !0.d0
     &  bc_uwall(ixlo:ixhi,iylo:iyhi,1:3)
#else
       utmp = 0.d0
#endif
c
         if(xprocs.eq.1) then 
           ilo = IXLO
           ihi = IXHI
         else
           if(iprocx.eq.1) then
             ilo = IXLO
           else
             ilo = 1
           endif
           if(iprocx.eq.xprocs) then
             ihi = IXHI
           else
             ihi = nxlocal+1
           endif 
         endif
c
         jlo = 1
         jhi = nylocal
c

         k=0
         do j=jlo, jhi,1
            do i=ilo, ihi,1
               do l=0,nghost-1,1
                  u(i,j,k-l,1)=2.d0*utmp(i,j,1)-u(i,j,k+1+l,1)
                  u(i,j,k-l,2)=2.d0*utmp(i,j,2)-u(i,j,k+1+l,2)
               enddo
c              u(i,j,k,1)=8.d0/3.d0*utmp(i,j,1)-2.d0*u(i,j,k+1,1)
c     &                  +1.d0/3.d0*u(i,j,k+2,1)
c              u(i,j,k-1,1) =  8.d0*utmp(i,j,1)-9.d0*u(i,j,k+1,1)
c     &                       +2.d0*u(i,j,k+2,1)
c              u(i,j,k-2,1) = 2.d0*utmp(i,j,1)-u(i,j,k+3,1)
c              u(i,j,k-3,1) = 2.d0*utmp(i,j,1)-u(i,j,k+4,1)
c              u(i,j,k,2)=8.d0/3.d0*utmp(i,j,2)-2.d0*u(i,j,k+1,2)
c     &                  +1.d0/3.d0*u(i,j,k+2,2)
c              u(i,j,k-1,2) =  8.d0*utmp(i,j,2)-9.d0*u(i,j,k+1,2)
c     &                       +2.d0*u(i,j,k+2,2)
c              u(i,j,k-2,2) = 2.d0*utmp(i,j,2)-u(i,j,k+3,2)
c              u(i,j,k-3,2) = 2.d0*utmp(i,j,2)-u(i,j,k+4,2)
            enddo
         enddo
         k=1
         do j=jlo, jhi,1
            do i=ilo, ihi,1
               u(i,j,k,3)= utmp(i,j,3)
               u(i,j,k-1,3)=2.d0*utmp(i,j,3)-u(i,j,k+1,3)
               u(i,j,k-2,3)=2.d0*utmp(i,j,3)-u(i,j,k+2,3)
               u(i,j,k-3,3)=2.d0*utmp(i,j,3)-u(i,j,k+3,3)
               u(i,j,k-4,3)=2.d0*utmp(i,j,3)-u(i,j,k+4,3)
            enddo
         enddo
#ifdef LESN
         k=0
         do j=jlo, jhi,1
            do i=ilo, ihi,1
               do l=0,nghost-1,1
                  u(i,j,k-l,1)= -bc_coef(l+1,1)*(u(i,j,k+1+l,1)-
     &                            utmp(i,j,1))-utmp(i,j,1)
                  u(i,j,k-l,2)= -bc_coef(l+1,2)*(u(i,j,k+1+l,2)-
     &                            utmp(i,j,2))-utmp(i,j,2)
               enddo
            enddo
         enddo
         k=0
         do j=jlo, jhi,1
            do i=ilo, ihi,1
               u(i,j,k+1,3)= utmp(i,j,3)
               do l=0, nghost-1, 1
                  u(i,j,k-l,3)= -bc_coef(l+1,3)*(u(i,j,k+2+l,3)-
     &                            utmp(i,j,3))-utmp(i,j,3)
               enddo
            enddo
         enddo
#endif
c
c
c         do j=IYLO, IYHI,1
c            do i=IXLO, IXHI,1
c               do l=0,nghost-1,1
c                  u(i,j,k-l,1:2)=2.d0*utmp(i,j,1:2)-u(i,j,k+1+l,1:2)
c               enddo
c     
c            enddo
c         enddo
c         k=1
c         do j=IYLO, IYHI,1
c            do i=IXLO, IXHI,1
c               u(i,j,k,3)= utmp(i,j,3)
c               u(i,j,k-1,3)=2.d0*utmp(i,j,3)-u(i,j,k+1,3)
c               u(i,j,k-2,3)=2.d0*utmp(i,j,3)-u(i,j,k+2,3)
c               u(i,j,k-3,3)=2.d0*utmp(i,j,3)-u(i,j,k+3,3)
c               u(i,j,k-4,3)=2.d0*utmp(i,j,3)-u(i,j,k+4,3)
cc     
c            enddo
c         enddo
c
      endif  !endif of iprocz = 1
c     
c     zr Boundary: Reflecting
      if (iprocz .eq. zprocs) then
         k=nzlocal+1
         do j=IYLO, IYHI,1
            do i=IXLO, IXHI,1
               do l=0,nghost-1,1
                  u(i,j,k+l,1)=u(i,j,k-1-l,1)
c                  u(i,j,k+l,1)= 2.D0-u(i,j,k-1-l,1)
                  u(i,j,k+l,2)=u(i,j,k-1-l,2)
               enddo
            enddo
         enddo
c     
         do j=IYLO, IYHI,1
            do i=IXLO, IXHI,1 
c               vt = 0.26121154d-3
c               u(i,j,k,3)=  vt !u(i,j,k-1,3) !0.D0
c               u(i,j,k+1,3)=2.d0*vt - u(i,j,k-1,3)
c               u(i,j,k+2,3)=2.d0*vt - u(i,j,k-2,3)
c               u(i,j,k+3,3)=2.d0*vt - u(i,j,k-3,3)
               u(i,j,k+1,3)= 2.d0*u(i,j,k,3)- u(i,j,k-1,3)
               u(i,j,k+2,3)= 2.d0*u(i,j,k,3)- u(i,j,k-2,3)
               u(i,j,k+3,3)= 2.d0*u(i,j,k,3)- u(i,j,k-3,3)
            enddo
         enddo
      endif
c     
c
c
c     step 3:
c     y component: set iprocy =1 when yprocs=1.
c
c     
#ifndef TWO_D
c     Periodic BC in Y
      if (yprocs .eq. 1) then
         j=0
         do k=IZLO,IZHI,1
            do i=IXLO,IXHI,1
               do l=0,nghost-1,1
                  u(i,j-l,k,1)=u(i,ny-l,k,1)
                  u(i,j-l,k,3)=u(i,ny-l,k,3)
               enddo
            enddo
         enddo
c         j=0
c     
         do k=IZLO,IZHI,1
            do i=IXLO,IXHI,1
               do l=0,nghost-1,1    
                  u(i,j-l,k,2)=u(i,ny-l,k,2)
               enddo
            enddo
         enddo
c
c
         j=ny+1
         do k=IZLO,IZHI,1
            do i=IXLO,IXHI,1
               do l=0,nghost-1,1
                  u(i,j+l,k,1)=u(i,l+1,k,1)
                  u(i,j+l,k,3)=u(i,l+1,k,3)
               enddo
            enddo
         enddo
c     
         do k=IZLO,IZHI,1
            do i=IXLO,IXHI,1
               do l=0,nghost-1,1
                  u(i,j+l,k,2)=u(i,l+1,k,2)
               enddo
            enddo
         enddo
      endif
c     
#endif
c
c
      return
      end subroutine
c
