c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c       Original version 
c_______________________________________________________________________
c        Spiral SGS model by D. I. Pullin
c            transfer from C to Fortran
c_______________________________________________________________________
       subroutine spiral_sgs_stress
     &   (u_near,v_near,w_near,x_near,y_near,z_near,
     &    yes_near,size_near,ix_this, del,
     &    dudx,e, K,
     &    Txx, Tyy, Tzz, Txy, Tyz, Tzx)
       use mesh
       use properties
c
       implicit none
c=======Declarations=========
       double precision:: u_near(0:2), v_near(0:2), w_near(0:2)
       double precision:: x_near(0:2), y_near(0:2), z_near(0:2)
       double precision:: dudx(0:2,0:2)
       double precision:: e(0:2), e_comp(0:2)
       double precision:: eigval(0:2)
       double precision:: Txx, Txy, Tzz, Tyy, Tyz, Tzx
       double precision:: Sxx, Syy, Szz, Sxy, Syz, Szx
       double precision:: Kc, lv, K, prefac
       double precision:: length, invlength
       double precision:: a, nu, delx, dely, delz, del
       double precision:: F_copy, Qd_copy, F2, Qd
       double precision:: du, dv, dw, dx_t, dy_t, dz_t, dx2_t, dxe, d
       double precision:: spiral_ke_integral, spiral_sf_integral
       double precision:: eps, pi
c
       integer::  yes_near(0:26)
       integer:: size_near, ix_this
       integer:: i !, j
c      
c       initial
        F_copy = 0.d0
        Qd_copy = 0.d0
        F2 = 0.d0
        Qd = 0.d0
        du = 0.d0
        dv = 0.d0
        dw = 0.d0
        dx_t = 0.d0
        dy_t = 0.d0
        dz_t = 0.d0
        dx2_t = 0.d0
        dxe = 0.d0
        d = 0.d0   
c      some value    
       eps = 1.d-20
       nu = 1.d0/Reynolds
       delx = dx
       dely = dy
       delz = dz
       pi = 4.d0*datan(1.d0)
       e_comp=0.d0
c      strain-rate tensor
       Sxx = 0.5d0* ( dudx(0,0)+dudx(0,0) )
       Syy = 0.5d0* ( dudx(1,1)+dudx(1,1) )
       Szz = 0.5d0* ( dudx(2,2)+dudx(2,2) )
       Sxy = 0.5d0* ( dudx(0,1)+dudx(1,0) )
       Syz = 0.5d0* ( dudx(1,2)+dudx(2,1) )
       Szx = 0.5d0* ( dudx(0,2)+dudx(2,0) )
c      
       if( e(0).eq.0.d0 .and. e(1).eq.0.d0 .and. e(2).eq.0.d0 ) then
          call spiral_eigenvalue_symm
     &         ( Sxx, Syy, Szz, Sxy, Syz, Szx, eigval);
c         Default Alignment: most extensivee eigenvector       
          call spiral_eigenvector_symm
     &         ( Sxx, Syy, Szz, Sxy, Syz, Szx, eigval(2), e)    
c         most compressible eigenvector
          call spiral_eigenvector_symm
     &         ( Sxx, Syy, Szz, Sxy, Syz, Szx, eigval(0), e_comp)
        endif
c        write(*,*) 'sum_eigenvalue', sum(eigval)        
c
c       Make e(0:2) a unit vector
        length = sqrt( e(0)*e(0)+ e(1)*e(1)+ e(2)*e(2) )
        invlength = 1.d0/length
        e(0) = e(0) * invlength
        e(1) = e(1) * invlength
        e(1) = e(2) * invlength
        !
        length = sqrt( e_comp(0)**2.d0+e_comp(1)**2.d0+e_comp(2)**2.d0 )
        invlength = 1.0/length
        e_comp(0) = e_comp(0)*invlength
        e_comp(1) = e_comp(1)*invlength
        e_comp(2) = e_comp(2)*invlength        
c
c       Strain along vortex axis
        a=  e(0)*e(0)*Sxx+e(0)*e(1)*Sxy+e(0)*e(2)*Szx
     &     +e(1)*e(0)*Sxy+e(1)*e(1)*Syy+e(1)*e(2)*Syz
     &     +e(2)*e(0)*Szx+e(2)*e(1)*Syz+e(2)*e(2)*Szz
        lv = dsqrt(2.d0*nu/(3.d0*(dabs(a)+eps)))
c       write(*,*) 'lv', lv
c
        do i=0, size_near-1, 1
          if(yes_near(i) .ne. 0) then
           du = u_near(i)-u_near(ix_this)
           dv = v_near(i)-v_near(ix_this)
           dw = w_near(i)-w_near(ix_this)
           F_copy= F_copy + du*du+dv*dv+dw*dw
c           
           dx_t = x_near(i) - x_near(ix_this)
           dy_t = y_near(i) - y_near(ix_this)
           dz_t = z_near(i) - z_near(ix_this)
           dx2_t = dx_t*dx_t+dy_t*dy_t+dz_t*dz_t
           dxe = dx_t*e(0)+dy_t*e(1)+dz_t*e(2)
           d = dsqrt(dx2_t-dxe*dxe)/del
c
           Qd_copy=Qd_copy+spiral_sf_integral(d)
          endif
        enddo
         F2 = F_copy
         Qd = Qd_copy
         F2 = F2/(size_near-1)
         Qd = Qd/(size_near-1)
c
         prefac = F2/Qd ! \mathcal{K}_0 \epsilon^{2/3} k_c^{-2/3}
         kc = pi /del
c
         K= prefac* spiral_ke_integral(kc*lv)
c
c        T_{ij}=( \delta_{ij} -e_i^v e_j^v ) K
         Txx = ( 1.d0 - e(0)*e(0) ) * K
         Tyy = ( 1.d0 - e(1)*e(1) ) * K
         Tzz = ( 1.d0 - e(2)*e(2) ) * K  
         Txy = (      - e(0)*e(1) ) * K
         Tyz = (      - e(1)*e(2) ) * K
         Tzx = (      - e(2)*e(0) ) * K
c
       return
       end

c--------------------------------------------------------------
c Calculate the SGS scalar flux (*qx, *qy, *qz) given the resolved scalar
c gradient dsdx[3], vortex alignment e[3] (a unit vector), LES cutoff
c scale del and precalculated SGS kinetic energy K.
c----------------------------------------------------------------------
       subroutine spiral_sgs_flux
     &  ( dsdx, e, del, K, qx, qy, qz)
       use mesh
       use properties
c
       implicit none
c=======Declarations=========        
       double precision:: dsdx(0:2), e(0:2)
       double precision:: del, K, qx, qy, qz
       double precision:: gam, P
c
       gam = 1.d0
       P = -0.5d0*gam*del*dsqrt(K)
c
        qx = P*( (1.d0-e(0)*e(0))*dsdx(0)
     &        +  (    -e(0)*e(1))*dsdx(1)
     &        +  (    -e(0)*e(2))*dsdx(2) )
        qy = P*( (    -e(1)*e(0))*dsdx(0)
     &        +  (1.d0-e(1)*e(1))*dsdx(1)
     &        +  (    -e(1)*e(2))*dsdx(2) )
        qz = P*( (    -e(2)*e(0))*dsdx(0)
     &        +  (    -e(2)*e(1))*dsdx(1)
     &        +  (1.d0-e(2)*e(2))*dsdx(2) )
c
       return
       end    


c______________________________________________________________
c Approximation of (1/2) k^(2/3) Gamma[-1/3, k^2]
c    with maximum relative error of 0.17% at k=2.42806.
c_______________________________________________________________
       double precision function spiral_ke_integral(k)
c
       implicit none
       double precision:: k, k2, pade
c
       k2 = k*k
       if(k2 .lt. 2.42806d0) then
        pade=(3.0d0+2.5107d0*k2+0.330357d0*k2*k2+0.0295481d0*k2*k2*k2)
     &    /(1.d0+0.336901d0*k2+0.0416684d0*k2*k2+0.00187191d0*k2*k2*k2)
        spiral_ke_integral = 
     &   0.5d0* ( pade - 4.06235d0 * k2**(1.d0/3.d0) )
       else
        pade = (1.26429d0+0.835714d0*k2+0.0964286d0*k2*k2)
     &       / (1.d0+2.25d0*k2+0.964286d0*k2*k2+0.0964286d0*k2*k2*k2)
        spiral_ke_integral = 
     &   0.5d0*pade* dexp(-k2)
       endif
c
       return
       end
c

c______________________________________________________________
c Approximation of
c Integrate[4 x^(-5/3) (1 - BesselJ[0, x Pi d]), {x, 0, 1}]
c with maximum relative error of 2.71% at d=0.873469.
c______________________________________________________________
       double precision function spiral_sf_integral(d)
       implicit none
       double precision:: d, d2, pi
       pi = 4.d0* datan(1.d0)
c     Uncomment if spherical averaging and d=1.
c     if (d == 1.0) return 4.09047;
       d2 = d * d;
       if (d .lt. 0.873469d0) then
         spiral_sf_integral = 7.4022d0*d2-1.82642d0*d2*d2
       else
         spiral_sf_integral =
     &      12.2946d0 * d**(2.d0/3.d0) - 6.d0
     &   - 0.573159d0 * d**(-1.5d0) * dsin(3.1415926d0*d-0.785398d0)
       endif
c
       return
       end
c
c__________________________________________________________________________________
c Calculate the eigenvalues, eigval[0] < eigval[1] < eigval[2],
c of the 3 x 3 symmetric matrix,
c { { Sxx, Sxy, Szx }, { Sxy, Syy, Syz }, { Szx, Syz, Szz } },
c assuming distinct eigenvalues.
c_________________________________________________________________________
       subroutine spiral_eigenvalue_symm
     &   (Sxx,Syy,Szz,Sxy,Syz,Szx,eigval)
c
       implicit none
       double precision:: eigval(0:2)
       double precision:: Sxx, Syy, Szz, Sxy, Syz, Szx
       double precision:: a, b, c, q, r, pai
       double precision:: costheta, theta
       double precision:: tmp, eps
       eps = 1.d-20
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
c
c________________________________________________________________________
c Calculate the eigenvector (not normalized), eigvec[3],
c corresponding to the precalculated eigenvalue, eigval,
c of the 3 x 3 symmetric matrix,
c { { Sxx, Sxy, Szx }, { Sxy, Syy, Syz }, { Szx, Syz, Szz } },
c assuming distinct eigenvalues.
c________________________________________________________________________
        subroutine spiral_eigenvector_symm
     &  (Sxx, Syy, Szz, Sxy, Syz, Szx, eigval, eigvec)
c
        implicit none
        double precision:: eigvec(3), det(3), fabsdet(3)
        double precision:: Sxx, Syy, Szz, Sxy, Syz, Szx
        double precision:: norm, eigval, norm_num
        double precision:: eps_norm, eps
        eps_norm = 1.d-6
        eps = 1.d-20 
c
c       normalized by Frobenius norm
        norm= dsqrt ( Sxx*Sxx + Syy*Syy + Szz*Szz
     &                +Sxy*Sxy + Syz*Syz + Szx*Szx + eps )
        norm_num = 
     &       dabs( (Sxx-eigval)*( (Syy-eigval)*(Szz-eigval)-Syz*Syz )
     &            + Sxy*( Syz*Szx - Sxy*(Szz-eigval) )
     &            + Szx*( Sxy*Syz - (Syy-eigval)*Szx) )        
        if( norm_num /norm  .gt. eps_norm ) then
            write(*,*) 'norm :', norm_num/norm, norm
            write(*,*) 'norm :', Sxx, Syy, Szz, Sxy, Szx, Syz
            write(*,*) 'spiral eigenvector symm: invalid eigenvalue'
            stop
        endif
c
        det(1)=(Syy-eigval)*(Szz-eigval)-Syz*Syz
        det(2)=(Szz-eigval)*(Sxx-eigval)-Szx*Szx
        det(3)=(Sxx-eigval)*(Syy-eigval)-Sxy*Sxy
        fabsdet(1)=dabs(det(1))
        fabsdet(2)=dabs(det(2))
        fabsdet(3)=dabs(det(3))    
c          
          if( fabsdet(1).ge.fabsdet(2) .and. 
     &        fabsdet(1).ge.fabsdet(3))  then
                 eigvec(1) = 1.d0
                 eigvec(2) = (-Sxy*(Szz-eigval)+Szx*Syz)/det(1)
                 eigvec(3) = (-Szx*(Syy-eigval)+Sxy*Syz)/det(1)
          else if( fabsdet(2).ge.fabsdet(3) .and. 
     &        fabsdet(2).gt.fabsdet(1))  then
                 eigvec(1) = (-Sxy*(Szz-eigval)+Syz*Szx)/det(2)
                 eigvec(2) = 1.d0
                 eigvec(3) = (-Syz*(Sxx-eigval)+Sxy*Szx)/det(2)
          else if( fabsdet(3).ge.fabsdet(1) .and. 
     &        fabsdet(3).gt.fabsdet(2))  then
                 eigvec(1) = (-Szx*(Syy-eigval)+Syz*Sxy)/det(3)
                 eigvec(2) = (-Syz*(Sxx-eigval)+Szx*Sxy)/det(3)
                 eigvec(3) = 1.d0
          else
              write(*,*) fabsdet(1),fabsdet(2),fabsdet(3),eigval
              write(*,*) 'spiral eigenvector symm: error'
                stop
          endif
c
         return
         end
ccc
cc
c
c------------------------------------------------------------------
c   another version of spiral_sgs_model used by Inoue 
c                                                     ver3
c------------------------------------------------------------------
c______________________________________________________________________
c            transfer from C to Fortran
c              originally from Inoue's code
c_______________________________________________________________________
       subroutine spiral_sgs_stress_v3
     &    ( u_near, v_near, w_near, x_near, y_near, z_near, 
     &      yes_near, size_near, ix_this, dudx, e, del, K,
     &      Txx, Tyy, Tzz, Txy, Tyz, Tzx)
       use mesh
       use properties
       implicit none
c=======Declarations=========
       double precision:: u_near(27), v_near(27), w_near(27)
       double precision:: x_near(27), y_near(27), z_near(27)
       double precision:: dudx(3,3)
       double precision:: e(3), e_comp(3)
       double precision:: eigval(3)
       double precision:: Sxx, Syy, Szz, Sxy, Syz, Szx
       double precision:: Txx, Txy, Tzz, Tyy, Tyz, Tzx
       double precision:: Kc, lv, K, prefac
       double precision:: length, invlength
       double precision:: a, nu, delx, dely, delz, del
       double precision:: F_copy, Qd_copy, F2, Qd
       double precision:: du, dv, dw, dx_t, dy_t, dz_t, dx2_t, dxe, sf_d
       double precision:: spiral_ke_integral, spiral_sf_integral
       double precision:: eps, pi, eps_control
       integer:: yes_near(27)
       integer::  size_near, ix_this
       integer:: i
c      
        eps_control = 2d-3
        eps = 1d-20
        pi = 4.d0*datan(1.d0)
c       initial
        F_copy = 0.d0
        Qd_copy = 0.d0
        F2 = 0.d0
        Qd = 0.d0
        du = 0.d0
        dv = 0.d0
        dw = 0.d0
        dx_t = 0.d0
        dy_t = 0.d0
        dz_t = 0.d0
        dx2_t = 0.d0
        dxe = 0.d0
        sf_d = 0.d0   
c      some value
       nu = 1.d0/Reynolds
       delx = dx
       dely = dy
       delz = dz
       e_comp=0.d0
c      strain-rate tensor
       Sxx = 0.5d0* ( dudx(1,1)+dudx(1,1) )
       Syy = 0.5d0* ( dudx(2,2)+dudx(2,2) )
       Szz = 0.5d0* ( dudx(3,3)+dudx(3,3) )
       Sxy = 0.5d0* ( dudx(1,2)+dudx(2,1) )
       Syz = 0.5d0* ( dudx(2,3)+dudx(3,2) )
       Szx = 0.5d0* ( dudx(1,3)+dudx(3,1) )
c        write(*,*) Sxx, Syy, Szz, Sxy, Syz, Szx
c      
       if( e(1).eq.0.d0 .and. e(2).eq.0.d0 .and. e(3).eq.0.d0 ) then
          call spiral_eigenvalue_symm
     &         ( Sxx, Syy, Szz, Sxy, Syz, Szx, eigval);
c         Default Alignment: most extensivee eigenvector       
          call spiral_eigenvector_symm
     &         ( Sxx, Syy, Szz, Sxy, Syz, Szx, eigval(3), e)    
c         most compressible eigenvector
          call spiral_eigenvector_symm
     &         ( Sxx, Syy, Szz, Sxy, Syz, Szx, eigval(1), e_comp)
        endif
c        write(*,*) 'sum_eigenvalue', sum(eigval)        
c
c       Make e(0:2) a unit vector
        length = dsqrt( e(1)*e(1)+ e(2)*e(2)+ e(3)*e(3) )
        invlength = 1.d0/length
        e(1) = e(1) * invlength
        e(2) = e(2) * invlength
        e(3) = e(3) * invlength
        !
        length= dsqrt(e_comp(1)**2.d0+e_comp(2)**2.d0+e_comp(3)**2.d0)
        invlength = 1.d0/(length+eps)
        e_comp(1) = e_comp(1)*invlength
        e_comp(2) = e_comp(2)*invlength
        e_comp(3) = e_comp(3)*invlength        
c
c       Strain along vortex axis
        a=  e(1)*e(1)*Sxx+e(1)*e(2)*Sxy+e(1)*e(3)*Szx
     &     +e(2)*e(1)*Sxy+e(2)*e(2)*Syy+e(2)*e(3)*Syz
     &     +e(3)*e(1)*Szx+e(3)*e(2)*Syz+e(3)*e(3)*Szz
        lv = dsqrt(2.d0*nu/(3.d0*(dabs(a)+eps_control)))
           
c        write(*,*) 'lv', lv
c
        del = delx*dabs(e_comp(1))
     &      + dely*dabs(e_comp(2))
     &      + delz*dabs(e_comp(3))
        del = del*0.333d0
c
        do i=1, size_near, 1
        if(yes_near(i).ne.0) then
           du = u_near(i)-u_near(ix_this)
           dv = v_near(i)-v_near(ix_this)
           dw = w_near(i)-w_near(ix_this)
           F_copy= F_copy + du*du+dv*dv+dw*dw
c           
           dx_t = x_near(i) - x_near(ix_this)
           dy_t = y_near(i) - y_near(ix_this)
           dz_t = z_near(i) - z_near(ix_this)
           dx2_t = dx_t*dx_t+dy_t*dy_t+dz_t*dz_t
           dxe = dx_t*e(1)+dy_t*e(2)+dz_t*e(3)
           sf_d = dsqrt(dx2_t-dxe*dxe+eps)/del
c
           Qd_copy=Qd_copy+spiral_sf_integral(sf_d)
c      if also calculating yes_near(i)=0, error of spiral integral is 10^{-12}
         endif
        enddo
         F2 = F_copy
         Qd = Qd_copy
         F2 = F2/(size_near-1)
         Qd = Qd/(size_near-1)
c
         prefac = F2/Qd ! \mathcal{K}_0 \epsilon^{2/3} k_c^{-2/3}
         kc = pi /del
c
         K= prefac* spiral_ke_integral(kc*lv)
c      
c        T_{ij}=( \delta_{ij} -e_i^v e_j^v ) K
         Txx = ( 1.d0 - e(1)*e(1) ) * K
         Tyy = ( 1.d0 - e(2)*e(2) ) * K
         Tzz = ( 1.d0 - e(3)*e(3) ) * K  
         Txy = (      - e(1)*e(2) ) * K
         Tyz = (      - e(2)*e(3) ) * K
         Tzx = (      - e(3)*e(1) ) * K
c
       return
       end
cc

c---------------------------------------------------------------
c       sgs model for virtual wall

       subroutine spiral_sgs_stress_v3b
     &    ( u_near, v_near, w_near, x_near, y_near, z_near,
     &      yes_near, size_near, ix_this, dudx, e, del, K,
     &      Txx, Tyy, Tzz, Txy, Tyz, Tzx)
       use mesh
       use properties
       implicit none
c=======Declarations=========
       double precision:: u_near(27), v_near(27), w_near(27)
       double precision:: x_near(27), y_near(27), z_near(27)
       double precision:: dudx(3,3)
       double precision:: e(3), e_comp(3)
       double precision:: eigval(3)
       double precision:: Sxx, Syy, Szz, Sxy, Syz, Szx
       double precision:: Txx, Txy, Tzz, Tyy, Tyz, Tzx
       double precision:: Kc, lv, K, prefac
       double precision:: length, invlength
       double precision:: a, nu, delx, dely, delz, del
       double precision:: F_copy, Qd_copy, F2, Qd
       double precision:: du, dv, dw, dx_t, dy_t, dz_t, dx2_t, dxe, sf_d
       double precision:: spiral_ke_integral, spiral_sf_integral
       double precision:: eps, pi, eps_control
       integer:: yes_near(27)
       integer::  size_near, ix_this
       integer:: i
c      
        eps_control = 2d-3
        eps = 1d-20
        pi = 4.d0*datan(1.d0)
c       initial
        F_copy = 0.d0
        Qd_copy = 0.d0
        F2 = 0.d0
        Qd = 0.d0
        du = 0.d0
        dv = 0.d0
        dw = 0.d0
        dx_t = 0.d0
        dy_t = 0.d0
        dz_t = 0.d0
        dx2_t = 0.d0
        dxe = 0.d0
        sf_d = 0.d0
c      some value
       nu = 1.d0/Reynolds
       delx = dx
       dely = dy
       delz = dz
       e_comp=0.d0
c      strain-rate tensor
       Sxx = 0.5d0* ( dudx(1,1)+dudx(1,1) )
       Syy = 0.5d0* ( dudx(2,2)+dudx(2,2) )
       Szz = 0.5d0* ( dudx(3,3)+dudx(3,3) )
       Sxy = 0.5d0* ( dudx(1,2)+dudx(2,1) )
       Syz = 0.5d0* ( dudx(2,3)+dudx(3,2) )
       Szx = 0.5d0* ( dudx(1,3)+dudx(3,1) )
c        write(*,*) Sxx, Syy, Szz, Sxy, Syz, Szx
c      
       if( e(1).eq.0.d0 .and. e(2).eq.0.d0 .and. e(3).eq.0.d0 ) then
          call spiral_eigenvalue_symm
     &         ( Sxx, Syy, Szz, Sxy, Syz, Szx, eigval);
c         Default Alignment: most extensivee eigenvector       
          call spiral_eigenvector_symm
     &         ( Sxx, Syy, Szz, Sxy, Syz, Szx, eigval(3), e)
c         most compressible eigenvector
          call spiral_eigenvector_symm
     &         ( Sxx, Syy, Szz, Sxy, Syz, Szx, eigval(1), e_comp)
        endif
c        write(*,*) 'sum_eigenvalue', sum(eigval)        
c
c       Make e(0:2) a unit vector
        length = dsqrt( e(1)*e(1)+ e(2)*e(2)+ e(3)*e(3) )
        invlength = 1.d0/length
        e(1) = e(1) * invlength
        e(2) = e(2) * invlength
        e(3) = e(3) * invlength
        !
        length= dsqrt(e_comp(1)**2.d0+e_comp(2)**2.d0+e_comp(3)**2.d0)
        invlength = 1.d0/(length+eps)
        e_comp(1) = e_comp(1)*invlength
        e_comp(2) = e_comp(2)*invlength
        e_comp(3) = e_comp(3)*invlength
c
c       Strain along vortex axis
        a=  e(1)*e(1)*Sxx+e(1)*e(2)*Sxy+e(1)*e(3)*Szx
     &     +e(2)*e(1)*Sxy+e(2)*e(2)*Syy+e(2)*e(3)*Syz
     &     +e(3)*e(1)*Szx+e(3)*e(2)*Syz+e(3)*e(3)*Szz
        lv = dsqrt(2.d0*nu/(3.d0*(dabs(a)+eps_control)))

c        write(*,*) 'lv', lv
c
c        del = delx*dabs(e_comp(1))
c     &      + dely*dabs(e_comp(2))
c     &      + delz*dabs(e_comp(3))
c        del = del*0.333d0
c
        do i=1, size_near, 1
        if(yes_near(i).ne.0) then
           du = u_near(i)-u_near(ix_this)
           dv = v_near(i)-v_near(ix_this)
           dw = w_near(i)-w_near(ix_this)
           F_copy= F_copy + du*du+dv*dv+dw*dw
c           
           dx_t = x_near(i) - x_near(ix_this)
           dy_t = y_near(i) - y_near(ix_this)
           dz_t = z_near(i) - z_near(ix_this)
           dx2_t = dx_t*dx_t+dy_t*dy_t+dz_t*dz_t
           dxe = dx_t*e(1)+dy_t*e(2)+dz_t*e(3)
           sf_d = dsqrt(dx2_t-dxe*dxe+eps)/del
c
           Qd_copy=Qd_copy+spiral_sf_integral(sf_d)
c      if also calculating yes_near(i)=0, error of spiral integral is
c               10^{-12}
         endif
        enddo
         F2 = F_copy
         Qd = Qd_copy
         F2 = F2/(size_near-1)
         Qd = Qd/(size_near-1)
c
         prefac = F2/Qd ! \mathcal{K}_0 \epsilon^{2/3} k_c^{-2/3}
         kc = pi /del
c
         K= prefac* spiral_ke_integral(kc*lv)
c      
c        T_{ij}=( \delta_{ij} -e_i^v e_j^v ) K
         Txx = ( 1.d0 - e(1)*e(1) ) * K
         Tyy = ( 1.d0 - e(2)*e(2) ) * K
         Tzz = ( 1.d0 - e(3)*e(3) ) * K
         Txy = (      - e(1)*e(2) ) * K
         Tyz = (      - e(2)*e(3) ) * K
         Tzx = (      - e(3)*e(1) ) * K
c
       return
       end



c-----------------------------------------
       subroutine spiral_sgs_stress_v4
     &    ( u_near, v_near, w_near, x_near, y_near, z_near, 
     &      yes_near, size_near, ix_this, dudx, e, del, K,
     &      Txx, Tyy, Tzz, Txy, Tyz, Tzx)
       use mesh
       use properties
       implicit none
c=======Declarations=========
       double precision:: u_near(27), v_near(27), w_near(27)
       double precision:: x_near(27), y_near(27), z_near(27)
       double precision:: dudx(3,3)
       double precision:: e(3), e_comp(3)
       double precision:: eigval(3)
       double precision:: Sxx, Syy, Szz, Sxy, Syz, Szx
       double precision:: Txx, Txy, Tzz, Tyy, Tyz, Tzx
       double precision:: Kc, lv, K, prefac
       double precision:: length, invlength
       double precision:: a, nu, delx, dely, delz, del
       double precision:: F_copy, Qd_copy, F2, Qd
       double precision:: du, dv, dw, dx_t, dy_t, dz_t, dx2_t, dxe, sf_d
       double precision:: spiral_ke_integral, spiral_sf_integral
       double precision:: eps, pi, eps_control
       integer:: yes_near(27)
       integer::  size_near, ix_this
       integer:: i
c      
        eps_control = 2d-3
        eps = 1d-20
        pi = 4.d0*datan(1.d0)
c       initial
        F_copy = 0.d0
        Qd_copy = 0.d0
        F2 = 0.d0
        Qd = 0.d0
        du = 0.d0
        dv = 0.d0
        dw = 0.d0
        dx_t = 0.d0
        dy_t = 0.d0
        dz_t = 0.d0
        dx2_t = 0.d0
        dxe = 0.d0
        sf_d = 0.d0   
c      some value
       nu = 0.d0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
       delx = dx
       dely = dy
       delz = dz
       e_comp=0.d0
c      strain-rate tensor
       Sxx = 0.5d0* ( dudx(1,1)+dudx(1,1) )
       Syy = 0.5d0* ( dudx(2,2)+dudx(2,2) )
       Szz = 0.5d0* ( dudx(3,3)+dudx(3,3) )
       Sxy = 0.5d0* ( dudx(1,2)+dudx(2,1) )
       Syz = 0.5d0* ( dudx(2,3)+dudx(3,2) )
       Szx = 0.5d0* ( dudx(1,3)+dudx(3,1) )
c        write(*,*) Sxx, Syy, Szz, Sxy, Syz, Szx
c      
       if( e(1).eq.0.d0 .and. e(2).eq.0.d0 .and. e(3).eq.0.d0 ) then
          call spiral_eigenvalue_symm
     &         ( Sxx, Syy, Szz, Sxy, Syz, Szx, eigval);
c         Default Alignment: most extensivee eigenvector       
          call spiral_eigenvector_symm
     &         ( Sxx, Syy, Szz, Sxy, Syz, Szx, eigval(3), e)    
c         most compressible eigenvector
          call spiral_eigenvector_symm
     &         ( Sxx, Syy, Szz, Sxy, Syz, Szx, eigval(1), e_comp)
        endif
c        write(*,*) 'sum_eigenvalue', sum(eigval)        
c
c       Make e(0:2) a unit vector
        length = dsqrt( e(1)*e(1)+ e(2)*e(2)+ e(3)*e(3) )
        invlength = 1.d0/(length+eps)
        e(1) = e(1) * invlength
        e(2) = e(2) * invlength
        e(3) = e(3) * invlength
        !
        length= dsqrt(e_comp(1)**2.d0+e_comp(2)**2.d0+e_comp(3)**2.d0)
        invlength = 1.d0/(length+eps)
        e_comp(1) = e_comp(1)*invlength
        e_comp(2) = e_comp(2)*invlength
        e_comp(3) = e_comp(3)*invlength        
c
c       Strain along vortex axis
        a=  e(1)*e(1)*Sxx+e(1)*e(2)*Sxy+e(1)*e(3)*Szx
     &     +e(2)*e(1)*Sxy+e(2)*e(2)*Syy+e(2)*e(3)*Syz
     &     +e(3)*e(1)*Szx+e(3)*e(2)*Syz+e(3)*e(3)*Szz
c        lv = dsqrt(2.d0*nu/(3.d0*(dabs(a)+eps_control)))
         lv = 0.d0           
c        write(*,*) 'lv', lv
c
c        del = delx*dabs(e_comp(1))
c     &      + dely*dabs(e_comp(2))
c     &      + delz*dabs(e_comp(3))
c        del = del*0.333d0
c
        do i=1, size_near, 1
        if(yes_near(i).ne.0) then
           du = u_near(i)-u_near(ix_this)
           dv = v_near(i)-v_near(ix_this)
           dw = w_near(i)-w_near(ix_this)
           F_copy= F_copy + du*du+dv*dv+dw*dw
c           
           dx_t = x_near(i) - x_near(ix_this)
           dy_t = y_near(i) - y_near(ix_this)
           dz_t = z_near(i) - z_near(ix_this)
           dx2_t = dx_t*dx_t+dy_t*dy_t+dz_t*dz_t
           dxe = dx_t*e(1)+dy_t*e(2)+dz_t*e(3)
           sf_d = dsqrt(dx2_t-dxe*dxe+eps)/del
c
           Qd_copy=Qd_copy+spiral_sf_integral(sf_d)
c      if also calculating yes_near(i)=0, error of spiral integral is 10^{-12}
         endif
        enddo
         F2 = F_copy
         Qd = Qd_copy
         F2 = F2/(size_near-1)
         Qd = Qd/(size_near-1)
c
         prefac = F2/Qd ! \mathcal{K}_0 \epsilon^{2/3} k_c^{-2/3}
         kc = pi /del
c
         K= prefac* spiral_ke_integral(kc*lv)
c      
c        T_{ij}=( \delta_{ij} -e_i^v e_j^v ) K
         Txx = ( 1.d0 - e(1)*e(1) ) * K
         Tyy = ( 1.d0 - e(2)*e(2) ) * K
         Tzz = ( 1.d0 - e(3)*e(3) ) * K  
         Txy = (      - e(1)*e(2) ) * K
         Tyz = (      - e(2)*e(3) ) * K
         Tzx = (      - e(3)*e(1) ) * K
c
       return
       end
cc


c
       subroutine spiral_sgs_stress_v5
     &    ( u_near, v_near, w_near, x_near, y_near, z_near, 
     &      yes_near, size_near, ix_this, dudx, e, del, K,
     &      Txx, Tyy, Tzz, Txy, Tyz, Tzx)
       use mesh
       use properties
       implicit none
c=======Declarations=========
       double precision:: u_near(27), v_near(27), w_near(27)
       double precision:: x_near(27), y_near(27), z_near(27)
       double precision:: dudx(3,3)
       double precision:: e(3), e_comp(3)
       double precision:: eigval(3)
       double precision:: Sxx, Syy, Szz, Sxy, Syz, Szx
       double precision:: Txx, Txy, Tzz, Tyy, Tyz, Tzx
       double precision:: Kc, lv, K, prefac
       double precision:: length, invlength
       double precision:: a, nu, delx, dely, delz, del
       double precision:: F_copy, Qd_copy, F2, Qd
       double precision:: du, dv, dw, dx_t, dy_t, dz_t, dx2_t, dxe, sf_d
       double precision:: spiral_ke_integral, spiral_sf_integral
       double precision:: eps, pi
       integer:: yes_near(27)
       integer::  size_near, ix_this
       integer:: i
c      
       eps = 1d-20 !!!
       pi = 4.d0*datan(1.d0)
c       initial
        F_copy = 0.d0
        Qd_copy = 0.d0
        F2 = 0.d0
        Qd = 0.d0
        du = 0.d0
        dv = 0.d0
        dw = 0.d0
        dx_t = 0.d0
        dy_t = 0.d0
        dz_t = 0.d0
        dx2_t = 0.d0
        dxe = 0.d0
        sf_d = 0.d0   
c      some value
       nu = 1.d0/Reynolds
       delx = dx
       dely = dy
       delz = dz
c
        lv = 0.d0
c
        do i=1, size_near, 1
        if(yes_near(i).ne.0) then
           du = u_near(i)-u_near(ix_this)
           dv = v_near(i)-v_near(ix_this)
           dw = w_near(i)-w_near(ix_this)
           F_copy= F_copy + du*du+dv*dv+dw*dw
c           
           dx_t = x_near(i) - x_near(ix_this)
           dy_t = y_near(i) - y_near(ix_this)
           dz_t = z_near(i) - z_near(ix_this)
           dx2_t = dx_t*dx_t+dy_t*dy_t+dz_t*dz_t
           dxe = dx_t*e(1)+dy_t*e(2)+dz_t*e(3)
           sf_d = dsqrt(dx2_t-dxe*dxe+eps)/del
c
           Qd_copy=Qd_copy+spiral_sf_integral(sf_d)
         endif
        enddo
         F2 = F_copy
         Qd = Qd_copy
         F2 = F2/(size_near-1)
         Qd = Qd/(size_near-1)
c
         prefac = F2/Qd ! \mathcal{K}_0 \epsilon^{2/3} k_c^{-2/3}
         kc = pi /del
c
         K= prefac* spiral_ke_integral(kc*lv)
c
c        T_{ij}=( \delta_{ij} -e_i^v e_j^v ) K
         Txx = ( 1.d0 - e(1)*e(1) ) * K
         Tyy = ( 1.d0 - e(2)*e(2) ) * K
         Tzz = ( 1.d0 - e(3)*e(3) ) * K  
         Txy = (      - e(1)*e(2) ) * K
         Tyz = (      - e(2)*e(3) ) * K
         Tzx = (      - e(3)*e(1) ) * K
c
       return
       end
cc














