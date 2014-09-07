c----------------------------------------------------------------------
c     Wan Cheng 
c     Feb 2013
c----------------------------------------------------------------------
c
c
c----------------------------------------------------------------------
c     set the virtual wall boundary condition
c----------------------------------------------------------------------
      subroutine les_virtualwall(ux,timestep,substep)
      use mesh
      use mesh_common
      use properties
      use time
      use virtualwallBC
      use stats
      use InflowVelBC
      use GClaw
#ifdef PARALLEL
      use mpistuff
#endif
c-----------------------------------------
      implicit none
      double precision:: ux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      integer:: timestep, substep
c
      double precision:: utaul,etal,usbc,vsbc,wsbc
      double precision:: k1, h0, hnup, h0p,gamma11,eps, nu, vtrans
      double precision:: GCprefac, del99
      double precision:: aa, bb 
      integer:: i,j, ilo  !, l
c
      if(iprocz.gt.1)  return
c
      nu = 1.d0/Reynolds
      h0 =  bc_h0dz * dz
      hnup = bc_hnuplus 
c
      eps = 1.0d-14
      aa =  9.d0/8.d0 
      bb = -1.d0/8.d0 
c       
c     get dudz 
      call les_dudz(ux,timestep,substep)
c
c     the virtual wall velocity: component u
      ilo = 1
      if(iprocx.eq.1) ilo = 2
      do i=ilo, nxlocal+1, 1
        do j=1, nylocal, 1
          utaul=   bc_utau(i,j)
          etal = bc_dudz(i,j)
          k1   = 0.384d0 ! bc_Kappa(i,j) 
          h0p  = h0*utaul/nu
          del99= st_del99(i)
          call getGCprefac(h0,utaul,nu,del99,GCprefac)
c          k1 at second point, but GCprefac is on h0.
          GCprefac = GCprefac/k1
          if (h0p .gt. GChnu ) then
c            usbc= utaul*(1.d0/(k1+eps)*dlog( h0p/hnup)+hnup)
            call getvtrans(GChnu,vtrans)
                 GC_ci(i,j)  = GCci
                 GC_co(i,j)  = GCco
                 GC_gam(i,j) = GCgam
                 GC_pre(i,j) = GCprefac             
c            usbc = utaul*(GCprefac*(h0p+GCap)**GCgam - 
c     &                    GCprefac*(GChnu+GCap)**GCgam +
c     &                    vtrans   )


            usbc = utaul*(GCci*(h0p+GCap)**GCgam -
     &                    GCci*(GChnu+GCap)**GCgam +
     &                    vtrans   )


c            note that in GC law GCap = - hnup
          else    
            call getvtrans(h0p,vtrans)
            usbc= utaul* vtrans
          endif
          bc_uwall(i,j,1)= usbc  
        enddo
      enddo 
c
c     the virtual wall velocity: component v
      do i=1, nxlocal, 1
        do j=1, nylocal+1, 1
          vsbc= 0.d0
          bc_uwall(i,j,2)= vsbc 
        enddo
      enddo

c     the virtual wall velocity: component w 
      do i=1, nxlocal, 1
        do j=1, nylocal, 1    
          wsbc =  -h0/2.d0*
     &         (bc_uwall(i,j,1)+bc_uwall(i+1,j,1))/2.d0
     &        *(bc_dudz(i+1,j) -bc_dudz(i,j))/dx
     &       /((bc_dudz(i+1,j) +bc_dudz(i,j))/2.d0)
c          wsbc =  -h0*
c     &         (bc_uwall(i,j,1)+bc_uwall(i+1,j,1))/2.d0
c     &        *(bc_utau(i+1,j) -bc_utau(i,j))/dx
c     &       /((bc_utau(i+1,j) +bc_utau(i,j))/2.d0)
c           wsbc = 0.d0
         bc_uwall(i,j,3)=  wsbc
        enddo
      enddo
c
c  bc uwall mpi
      call setbcbcuwall
c  fft bcuwall
c      call fftbcuwall
#ifdef PARALLEL
      call exchangebcuwall
#endif
c
c      return  ! code return here
c   check output part
c      write(*,*)'bc_Kappa MAX/MIN :',maxval(bc_kappa),minval(bc_kappa)
c     &            , iprocx, iprocy, iprocz
c      write(*,*)'bc_uwall1 MAX/MIN :', 
c     &          maxval(bc_uwall(1:nxlocal+1,1:nylocal,1)),
c     &          minval(bc_uwall(1:nxlocal+1,1:nylocal,1))
c      write(*,*) 'bc_uwall3 MAX/MIN :',
c     &          maxval(bc_uwall(1:nxlocal,1:nylocal,3)),
c     &          minval(bc_uwall(1:nxlocal,1:nylocal,3)) 
c
      return
      end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
      subroutine getvtrans(zp,vtrans)
      use GClaw
c
      implicit none
      double precision:: vtrans, zp
      double precision:: c4, c5, dpara
c
      c4 = -0.0003d0
      c5 = 1.35d-5
      dpara = 8.0d-8
c
       vtrans = (zp+c4*zp**4.d0+c5*zp**5.d0)*dexp(-dpara*zp**6.d0)+    
     &      GCci*zp**GCgam*
     &      (1.d0+GCgam*GCap/zp + 
     &       0.5d0*GCgam*(GCgam-1.d0)*(GCap/zp)**2.d0)
     &      *(1.d0-dexp(-dpara*zp**6.d0)  )

c        write(*,*) zp, 
c     & GCci,zp,GCgam, GCRe
c     &    , vtrans
c
      return 
      end
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
      subroutine getGCprefac(zh,utau,nu,del99,GCprefac)
      use GClaw
c
      implicit none
      double precision:: zh, zp, eta, utau, nu, del99, GCprefac
c
      zp    = zh*utau/nu
      GCRe  = del99*utau/nu
      GCgam = GCgaminfty + GCalpha*GCa/(dlog(GCRe)**(1.d0+GCalpha)) 
      GCprefac = (zp + GCap)**(1.d0-GCgam)/zp/GCgam
      GCco  = GCcoinfty*(1.d0+GCc1*dexp(GCc2*GCRe))
      GCci  = GCciinfty/GCcoinfty*GCco*
     &        dexp(-(1.d0+GCalpha)*GCa/(dlog(GCRe)**GCalpha))   
c
      return
      end
c
c----------------------------------------------------------------------
c      calculate dudz at wall 
c----------------------------------------------------------------------
      subroutine les_dudz(ux,timestep,substep)
      use mesh
      use mesh_common
      use properties
      use time
      use RKCoefficients
      use virtualwallBC
#ifdef PARALLEL
      use mpistuff
#endif
c------------------------------------
      implicit none
      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      integer:: timestep, substep
c
      double precision:: eta0(1:nxlocal+1,1:nylocal)
      double precision:: etap(1:nxlocal+1,1:nylocal)
      double precision:: nu
      integer:: i, j, ilo
c
      nu = 1.d0/Reynolds
c
c     get the right hand side of ode 
      call  paraeta(ux,timestep,substep)
c
      ilo = 1
      if(iprocx.eq.1)  ilo = 2
c
      do i = ilo, nxlocal+1, 1
        do j=1, nylocal, 1
          eta0(i,j) =  bc_dudz(i,j) 
c          if(substep.eq.1) then
c            etap(i,j) = eta0(i,j) + dt*  bc_rhs(i,j)
c            bc_rhs_old(i,j) = bc_rhs(i,j)
c          else
c            return
c          endif
          if(substep.eq.1) then
            etap(i,j) = eta0(i,j) + dt*  gamma(1)*bc_rhs(i,j)
            bc_rhs_old(i,j) = bc_rhs(i,j)
          endif
          if(substep.eq.2) then
            etap(i,j) = eta0(i,j) + dt* (gamma(2)*bc_rhs(i,j)
     &                                   + zeta(1)*bc_rhs_old(i,j) )
            bc_rhs_old(i,j) = bc_rhs(i,j)
          endif
          if(substep.eq.3) then
            etap(i,j) = eta0(i,j) + dt* (gamma(3)*bc_rhs(i,j)
     &                                   + zeta(2)*bc_rhs_old(i,j) )
            bc_rhs_old(i,j) = 1.d7
          endif         
c         update the viscous velocity
          bc_dudz(i,j) = etap(i,j)
        enddo
      enddo 

c  fft trancate
c      call fftbcdudz

      do i = 1, nxlocal+1, 1
        do j=1, nylocal, 1
          bc_utau(i,j) = dsqrt(dabs(bc_dudz(i,j))*nu)
        enddo
      enddo
c
      return   ! code return here 
c     check the values
      write(*,*) maxval(bc_utau),maxval(etap), minval(etap),
     &           maxval(eta0), iprocx, iprocy, iprocz
      end
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c      pressure gradient 
c----------------------------------------------------------------------
      subroutine les_dpdx(press)
      use mesh
      use properties
      use virtualwallBC
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision:: press(-2:nxlocal+3,-2:nylocal+3,-2:nzlocal+3)
c
      double precision:: aa,bb
      integer:: i,j,k
c
c     at the point of u
      aa=9.D0/8.D0
      bb=-1.D0/8.D0
      k = 1
      do i=1, nxlocal+1, 1
        do j =1, nylocal, 1
        bc_dpdx(i,j)= aa*(press(i  ,j,k)-press(i-1,j,k))*dxi
     &               +bb*(press(i+1,j,k)-press(i-2,j,k))*dxi3
        enddo
      enddo
c
      return
      end
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c      set the RHS term bc_rhs for solving \eta 
c----------------------------------------------------------------------
      subroutine paraeta(ux,timestep,substep)
      use mesh
      use mesh_common
      use properties
      use time
      use virtualwallBC
      use stats
#ifdef PARALLEL
      use mpistuff
#endif
c
      implicit none
      double precision:: ux(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      integer:: timestep, substep
c
      double precision::     eta(1:nxlocal+1,1:nylocal)
      double precision:: paralin(1:nxlocal+1,1:nylocal)
      double precision:: aa, bb, hh, uh, fuh, fGCRe, del99
      double precision:: linear, quadra, nu
      integer:: i,j  !,k
c
      nu = 1.d0/Reynolds
      aa = 9.d0/8.d0
      bb = -1.d0/8.d0
c
      call  les_bc_convec(ux)
c p1 filter
      call  fftbcnonlinear
      call  fftbcuw
c      call  fftbcdudzh
c      call  fftbcdpdx
c
      hh = bc_hdz*dz
      do i=1, nxlocal+1, 1
        do j=1, nylocal, 1
        bc_uw(i,j)     = bc_uw(i,j) + bc_Txz(i,j)
        bc_convec(i,j) = bc_convec(i,j) +
     &                   aa*(bc_Txx(i  ,j)-bc_Txx(i-1,j))*dxi+
     &                   bb*(bc_Txx(i+1,j)-bc_Txx(i-2,j))*dxi3+
     &                   aa*(bc_Txy(i,j+1)-bc_Txy(i,j  ))*dyi+
     &                   bb*(bc_Txy(i,j+2)-bc_Txy(i,j-1))*dyi3
        paralin(i,j)   =-bc_uw(i,j) + nu*bc_dudzh(i,j) -
     &                   hh*(bc_dpdx(i,j)+bc_convec(i,j))
        enddo
      enddo
      do i=1,nxlocal+1,1
        do j=1, nylocal, 1
        eta(i,j)   =  bc_dudz(i,j)
        linear     =  eta(i,j)*paralin(i,j)
        quadra     = -eta(i,j)*eta(i,j)*nu
        uh         =  ux(i,j,1,1)
c        fuh        = uh
        del99      =  st_del99(i)
        call getGCRe(hh,eta,nu,del99,fGCRe)
        fuh        = uh + uh*uh/2.d0*fGCRe
        bc_rhs(i,j)=  2.d0/(fuh*hh)*(linear+quadra)
        enddo
      enddo

c p1 filter
c      call fftbcrhs

      return  ! code return here
c
      end
c
c----------------------------------------------------------------------
c     get GCRe
c----------------------------------------------------------------------
      subroutine getGCRe(zh,eta,nu,del99,FRe)
      use GClaw
c
      implicit none
      double precision:: zh, zp, eta, utau, nu, del99
      double precision:: G1, G2, H1, FRe
c
      utau  = (eta*nu)**0.5d0
      zp    = zh*utau/nu
      GCRe  = del99*utau/nu
      GCgam = GCgaminfty + GCalpha*GCa/(dlog(GCRe)**(1.d0+GCalpha)) 
      GCco  = GCcoinfty*(1.d0+GCc1*dexp(GCc2*GCRe))
      GCci  = GCciinfty/GCcoinfty*GCco*
     &        dexp(-(1.d0+GCalpha)*GCa/(dlog(GCRe)**GCalpha))     
c
      G1 = GCc1*GCc2*GCRe/(dexp(-GCc2*GCRe)+GCc1)
      G2 = (1.d0+GCalpha)*GCa*GCalpha*(dlog(GCRe)-dlog(zp+GCap))
     &      /(dlog(GCRe)**(2.d0+GCalpha))
      H1 = GCgam
      FRe= (G1+G2)/H1 
c
      return
      end
c
c----------------------------------------------------------------------
c     get the convective term for RHS of ode
c----------------------------------------------------------------------
      subroutine les_bc_convec(ux)
      use mesh
      use properties
      use virtualwallBC
#ifdef PARALLEL
      use mpistuff
#endif
c
      implicit none
      double precision:: ux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c
      double precision:: tmp1(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI)
      double precision:: tmp3(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI)
      double precision:: adveta(1:nxlocal+1,1:nylocal)
      double precision:: tmp,uc1,uc3, v1,v3,u1,u3
      double precision:: aa,bb
      integer:: i,j,k
c
      aa=9.D0/8.D0
      bb=-1.D0/8.D0
c
      k =1
      do i=-1,nxlocal+2,1
        do j=1,nylocal,1
          uc1=0.5D0*(ux(i  ,j,k,1)+ux(i+1,j,k,1))
          uc3=0.5D0*(ux(i+2,j,k,1)+ux(i-1,j,k,1))
          tmp=aa*uc1+bb*uc3
          tmp1(i,j,k)=tmp*uc1
          tmp3(i,j,k)=tmp*uc3
         enddo
      enddo
      do i=1,nxlocal+1,1
        do j=1,nylocal,1
        adveta(i,j)=aa*(tmp1(i  ,j,k)-tmp1(i-1,j,k))*dxi+
     &               bb*(tmp3(i+1,j,k)-tmp3(i-2,j,k))*dxi3
        enddo
      enddo
c     Advection term udu/dx
      k = 1
      do i=-1,nxlocal+2,1
        do j=1,nylocal,1
          uc1=0.5D0*(ux(i,j,k,1)+ux(i+1,j,k,1))
          uc3=0.5D0*(ux(i+2,j,k,1)+ux(i-1,j,k,1))
          tmp=aa*uc1+bb*uc3
          tmp1(i,j,k)=tmp*(ux(i+1,j,k,1)-ux(i  ,j,k,1))*dxi
          tmp3(i,j,k)=tmp*(ux(i+2,j,k,1)-ux(i-1,j,k,1))*dxi3
        enddo
      enddo
      do i=1,nxlocal+1,1
        do j=1,nylocal,1
          adveta(i,j)=0.5D0*adveta(i,j)+
     &                0.25D0*(aa*(tmp1(i  ,j,k)+tmp1(i-1,j,k))
     &                       +bb*(tmp3(i+1,j,k)+tmp3(i-2,j,k)))
        enddo
      enddo
c     the second term
      k = 1
      do i=1,nxlocal+1,1
        do j=0,nylocal+2,1
          u1=0.5D0*(ux(i,j-1,k,1)+ux(i,j,k,1))
          u3=0.5D0*(ux(i,j-2,k,1)+ux(i,j+1,k,1))
          v1=0.5D0*(ux(i-1,j,k,2)+ux(i,j,k,2))
          v3=0.5D0*(ux(i-2,j,k,2)+ux(i+1,j,k,2))
          tmp=aa*v1+bb*v3
          tmp1(i,j,k)=tmp*u1
          tmp3(i,j,k)=tmp*u3
        enddo
      enddo
c
      k = 1
      do i=1,nxlocal+1,1
        do j=1,nylocal,1
          adveta(i,j)=adveta(i,j)+0.5D0*
     &               (aa*(tmp1(i,j+1,k)-tmp1(i,j  ,k))*dyi
     &               +bb*(tmp3(i,j+2,k)-tmp3(i,j-1,k))*dyi3)
        enddo
      enddo
c
c     Second term vd(u)/dy: adv part
c
      k = 1
      do i=1,nxlocal+1,1
        do j=0,nylocal+2,1
          v1=0.5D0*(ux(i-1,j,k,2)+ux(i,j,k,2))
          v3=0.5D0*(ux(i-2,j,k,2)+ux(i+1,j,k,2))
          tmp=aa*v1+bb*v3
          tmp1(i,j,k)=tmp*(ux(i,j  ,k,1)-ux(i,j-1,k,1))*dyi
          tmp3(i,j,k)=tmp*(ux(i,j+1,k,1)-ux(i,j-2,k,1))*dyi3
        enddo
      enddo
c
      k = 1
      do i=1,nxlocal+1,1
         do j=1,nylocal,1
           adveta(i,j)=adveta(i,j)+0.25D0*
     &               (aa*(tmp1(i,j+1,k)+tmp1(i,j  ,k))
     &               +bb*(tmp3(i,j+2,k)+tmp3(i,j-1,k)))
         enddo
      enddo
c
      bc_convec(1:nxlocal+1,1:nylocal) = adveta(1:nxlocal+1,1:nylocal) 
c
      return  ! code return here
c
      end
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c     set bc uwall for ghost points
c----------------------------------------------------------------------
      subroutine setbcbcuwall
      use mesh
      use mesh_common
      use virtualwallbc
      use InflowVelBC
c
      implicit none
      integer:: i,j,k

! no used interpolation, right now everything recycled
#ifdef CODEAN
      if (iprocz.eq.1 .and. iprocx .eq. 1) then
        i = 1
        do j=1,nylocal,1
          bc_uwall(i-1,j,1) = 2.d0*bc_uwall(i,j,1)-bc_uwall(i+1,j,1)
          bc_uwall(i-2,j,1) = 2.d0*bc_uwall(i,j,1)-bc_uwall(i+2,j,1)
          bc_uwall(i-3,j,1) = 2.d0*bc_uwall(i,j,1)-bc_uwall(i+3,j,1)
          bc_uwall(i-4,j,1) = 2.d0*bc_uwall(i,j,1)-bc_uwall(i+4,j,1)
c          bc_uwall(i-1,j,2) = 2.d0*bc_uwall(i,j,2)-bc_uwall(i+1,j,2)
c          bc_uwall(i-2,j,2) = 2.d0*bc_uwall(i,j,2)-bc_uwall(i+2,j,2)
c          bc_uwall(i-3,j,2) = 2.d0*bc_uwall(i,j,2)-bc_uwall(i+3,j,2)
c          bc_uwall(i-4,j,2) = 2.d0*bc_uwall(i,j,2)-bc_uwall(i+4,j,2)
           bc_uwall(i-4:i-1,j,2) = 0.d0
c         for w component, now consistent
          bc_uwall(i-1,j,3) = 2.d0*Uin(1,j,1,3)-bc_uwall(i,j,3)
          bc_uwall(i-2,j,3) = 2.d0*Uin(1,j,1,3)-bc_uwall(i+1,j,3)
          bc_uwall(i-3,j,3) = 2.d0*Uin(1,j,1,3)-bc_uwall(i+2,j,3)
          bc_uwall(i-4,j,3) = 2.d0*Uin(1,j,1,3)-bc_uwall(i+3,j,3)
c           bc_uwall(i-1,j,1) =  bc_uwall(i+1,j,1) ! uedge
c           bc_uwall(i-2,j,1) =  bc_uwall(i+2,j,1)
c           bc_uwall(i-3,j,1) =  bc_uwall(i+3,j,1)
c           bc_uwall(i-4,j,1) =  bc_uwall(i+4,j,1)
c           bc_uwall(i-1,j,2:3) = bc_uwall(i,j,2:3) ! w edge
c           bc_uwall(i-2,j,2:3) = bc_uwall(i+1,j,2:3)
c           bc_uwall(i-3,j,2:3) = bc_uwall(i+2,j,2:3)
c           bc_uwall(i-4,j,2:3) = bc_uwall(i+3,j,2:3)
        enddo
      endif
#endif
c
      if ( iprocz.eq.1 .and. iprocx .eq. xprocs) then
        i = nxlocal+1
        do j=1,nylocal,1
          bc_uwall(i+1,j,1) = 2.d0*bc_uwall(i,j,1)-bc_uwall(i-1,j,1)
          bc_uwall(i+2,j,1) = 2.d0*bc_uwall(i+1,j,1)-bc_uwall(i,j,1)
          bc_uwall(i+3,j,1) = 2.d0*bc_uwall(i+2,j,1)-bc_uwall(i+1,j,1)
          bc_uwall(i,j,2)   = 2.d0*bc_uwall(i-1,j,2)-bc_uwall(i-2,j,2)
          bc_uwall(i+1,j,2) = 2.d0*bc_uwall(i,j,2)-bc_uwall(i-1,j,2)
          bc_uwall(i+2,j,2) = 2.d0*bc_uwall(i+1,j,2)-bc_uwall(i,j,2)
          bc_uwall(i+3,j,2) = 2.d0*bc_uwall(i+2,j,2)-bc_uwall(i+1,j,2)
          bc_uwall(i,j,3)   = 2.d0*bc_uwall(i-1,j,3)-bc_uwall(i-2,j,3)
          bc_uwall(i+1,j,3) = 2.d0*bc_uwall(i,j,3)-bc_uwall(i-1,j,3)
          bc_uwall(i+2,j,3) = 2.d0*bc_uwall(i+1,j,3)-bc_uwall(i,j,3)
          bc_uwall(i+3,j,3) = 2.d0*bc_uwall(i+2,j,3)-bc_uwall(i+1,j,3)
        enddo
      endif
c spanwise only for yprocs=1
      if(yprocs.eq.1) then
      do i= ixlo, ixhi, 1
        do j= 1, nghost
          bc_uwall(i,nylocal+j,1:3) = bc_uwall(i,j,1:3)
          bc_uwall(i,1-j,1:3) = bc_uwall(i,nylocal+1-j,1:3)
        enddo
      enddo
      endif
c
      return  ! code return here
c
      end
c----------------------------------------------------------------------
