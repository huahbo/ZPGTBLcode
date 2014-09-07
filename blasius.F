c------------------------------------------------------
c      code to output blasius solution as 
c                     inflow for DNS or LES simulation 
c------------------------------------------------------
       subroutine blasiusinflow(ux)
       use mesh
       use mesh_common
       use properties
       use InflowVelBC
#ifdef PARALLEL
       use mpistuff
#endif
c
       implicit none
       double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)

       double precision, allocatable:: fv(:)
       double precision, allocatable:: fd(:) 
       double precision:: ut(nz,2)
       double precision:: etamax, etamin,  x_in, zdomain, etav, zbottom
       double precision:: pert, aa, bb, u_tau, pai
       integer:: etano, streamno, wallno
       integer:: i,j, k, l
c
       aa = 8.d0/12.d0
       bb = -1.d0/12.d0
       pai = 4.d0 *datan(1.d0)
c
c      delta = 5 sqrt(x/Re)=1 for Blasius solution
       x_in = 0.04*Reynolds 
c---DNS init
        zdomain = dz*nz
        etano = 2*nz*100
c
        allocate( fv(etano) )
        allocate( fd(etano) )
c
        etamax = zdomain*dsqrt(1.d0/2.d0*Reynolds/x_in)
c        write(*,*)'Blasius :', etamax, zdomain
c main solver
        call self_similar(etamax, etano, fv, fd)
c 
        do k=1, nz
          streamno = (k-1)*200+101 ! streamwise no. for staggered
          ut(k,1) = fd( streamno )
          wallno = (k-1)*200+1
          etav = etamax/nz*(k-1)
          ut(k,2) = dsqrt(1.d0/(x_in*Reynolds*2.d0))
     &       *(etav*fd(wallno)-fv(wallno))
        enddo

c set uin
        if(zprocs.eq.1)  then
          do j= iylo, iyhi
            Uin(1,j,izlo:0,1) = 0.d0
            Uin(1,j,1:nzlocal,1) = ut(1:nz,1)
            Uin(1,j,izlo:0,3) = 0.d0
            Uin(1,j,1:nzlocal,3) = ut(1:nz,2)
            Uin(1,j,nzlocal+1,3) = Uin(1,j,nzlocal,3)
          enddo
        else
          if(iprocz.eq.1) then
            do j=iylo, iyhi, 1
              Uin(1,j,izlo:0,1) = 0.d0
              Uin(1,j,izlo:0,3) = 0.d0
              Uin(1,j,1:izhi,1) = ut(1:izhi,1)
              Uin(1,j,1:izhi,3) = ut(1:izhi,2)
            enddo
          else
            if(iprocz.eq.zprocs) then
              do j=iylo,iyhi, 1
                Uin(1,j,izlo:nzlocal,1)=
     &          ut((iprocz-1)*nzlocal+izlo:iprocz*nzlocal,1)
                Uin(1,j,izlo:nzlocal,3)=
     &          ut((iprocz-1)*nzlocal+izlo:iprocz*nzlocal,2)
                Uin(1,j,nzlocal+1:izhi,1) = 1.d0
                Uin(1,j,nzlocal+1:izhi,3) = Uin(1,j,nzlocal,3) 
              enddo
            else
              do j=iylo,iyhi, 1
                Uin(1,j,izlo:izhi,1)=
     &          ut((iprocz-1)*nzlocal+izlo:(iprocz-1)*nzlocal+izhi,1)
                Uin(1,j,izlo:izhi,3)=
     &          ut((iprocz-1)*nzlocal+izlo:(iprocz-1)*nzlocal+izhi,2)
              enddo
            endif            
          endif  
        endif
        Uin(1,iylo:iyhi,izlo:izhi,2) = 0.d0
c--  init
        do i=ixlo, ixhi, 1
           do j=iylo, iyhi, 1
             do k=izlo, izhi, 1
             ux(i,j,k,1) = Uin(1,j,k,1) 
             ux(i,j,k,2) = 0.d0
             ux(i,j,k,3) = 0.d0
             enddo
           enddo
        enddo
c-- check 
c        j=1
c       do k =1, nzlocal
c         write(*,*) Uin(1,j,k,1), Uin(1,j,k,3)
c       enddo 
c       write(*,*) 'wall :', fd(19), Uin(1,1,1,3)
c
       return
       end subroutine  blasiusinflow
ccccccccccccccccccccccccc------------------------------------------------------------------------------

        subroutine self_similar(eta_max,j_max,f,g)
        implicit none
        integer::  j, j_max, n, n_max !, i
        double precision:: eta_min , eta_max , eta_cur, k , !eta_init,
     &     delta_j_init , delta_j , beta , p , q
        double precision:: eta(j_max), g(j_max),
     &       a(j_max),b(j_max),c(j_max),r(j_max),f(j_max) !,vv(j_max)
c
          delta_j_init = eta_max/j_max
          k = 1
          n_max = 1000
          beta = 0.d0
        eta_min = 0.0d0
        delta_j = delta_j_init
        eta_cur = eta_min+delta_j
c
        eta(1) = eta_min
        delta_j = delta_j_init
        g(1) = 0.0d0
        j = 2
        do while ( j .le. j_max )
           eta(j) = eta(j-1)+delta_j
           g(j) = 1.0d0
           delta_j = k* delta_j
           j = j+1
        enddo
C . . . Main loop for solving the self-similar equation
        n = 1
        do while (n .le. n_max )
          call integrate_g ( g , eta , f , j_max )
          j = 2
          do while(j.le. j_max -1)
             p = eta(j+1)-eta(j)
             q = eta(j)-eta(j-1)
           a(j) =2.0d0/((p+q)*q)-f(j)/(p+q)
           b(j) = -2.0d0/(p+q)*(1/p+1/q )-beta*2*g(j)
           c(j) = 2.0d0/((p+q)*p)+f(j)/(p+q)
           r(j) = -beta*((g(j)**2)+1)
           j = j+1
          enddo
          r(2) =r(2)-a(2)*g(1)
          r(j_max-1) = r(j_max-1) - c( j_max-1)*g(j_max )
          call thomas( a, b, c,r,g, j_max )
          n = n+1
      enddo

      end subroutine self_similar
C . . . Write a subroutine to out source the Thomas algorithm .
          subroutine thomas( a , b , c , r , g , j_max )
          double precision:: a(j_max) , b(j_max), c(j_max) , r(j_max) ,
     &          g(j_max) , q
          integer:: j, j_max
C . . . Gaussian elimination
          j = 3
          do while (j .le. j_max-1)
             q = a(j)/b(j-1)
             b(j) = b(j) - q*c(j-1)
             r(j) = r(j) - q*r(j-1)
             j = j+1
          enddo
C . . . Resubstituion
          g(j_max-1) = r(j_max-1)/b(j_max-1)
          j = j_max-2
          do while ( j .ge. 2 )
            g(j) = (r(j)-c(j)*g(j+1))/b(j)
            j= j-1
          enddo
       end

C . . . Write a subroutine to calculate f as the integral of g over eta .
        subroutine integrate_g(g,eta,f,j_max)
        double precision:: g(j_max) , eta(j_max) ,f(j_max )
        integer:: j, j_max
C . . . The initial condition for eta = 0 is f = 0
         f(1) = 0.0d0
          j = 2
         do while(j.le. j_max )
           f(j) = f(j-1) + (eta(j)-eta(j-1))*(g(j-1)+g(j))/2.0d0
            j = j+1
         enddo
         end

