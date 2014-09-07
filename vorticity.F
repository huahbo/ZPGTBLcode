c	Jaideep Ray
c	March 7, 1997

c-----------------------------------------------------------------------
c	Subroutine to calculate the cirulation in thex, y, z directions.
c	INPUTS : ux, xc [ux includes 2 updated ghost cells per cuboid face]
c	OUTPUTS : tvortx, tvorty, tvortz, tnvortx, tnvorty, tnvortz
c	          tvortxint, tvortyint, tvortzint, 
c	          tnvortxint, tnvortyint, tnvortzint
c------------------------------------------------------------------------

	subroutine VortSpaceTime(ux, iter, ndump)

	include 'mesh_uparms.h'
	include 'mesh_parms.h'
	include 'mesh_common.h'
	include 'mesh.h'
#ifdef PARALLEL
	include 'mpif.h'
	include 'mpistuff.h'
#endif
        real*8 ux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,6)
#ifdef TWO_D
	real*8 omega(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI)
	real*8 xi(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI)
	real*8 a(NXlsize, NYlsize, NZlsize)

	real*8 circ, circn
	real*8 circ_global, circn_global

	real*8 circint, circintn
	real*8 circint_global, circintn_global


	real*8 omega_max, omega_max_local
	real*8 omega_min, omega_min_local
	real*8 gammafun, gamma, u, v, w, p 	       
	real*8  wijp1k, wijm1k, vijkp1, vijkm1 
	real*8  wip1jk,wim1jk,uijkp1,uijkm1 
	real*8  vip1jk, vim1jk, uijp1k,uijm1k
	integer i, j, k, iter, ndump, index

	common/time/dt, ttot
	real*8 dt, ttot

	circ = 0.D0
	circn = 0.D0
	circint = 0.D0
	circintn = 0.D0

c	---- Determine xi (interface tracker)

	do k = 1, NZlocal
	   do j = 1, NYlocal
	      do i = 1, NXlocal
		 xi(i,j,k) = ux(i,j,k,6) / ux(i,j,k,1)
	      enddo
	   enddo
	enddo

c       ---- Determine vorticity

	do k = 1, NZlocal
	   do j = 1, NYlocal
	      do i = 1, NXlocal
		 
		 vip1jk = ux(i+1,j,k,3) / ux(i+1,j,k,1)
		 vim1jk = ux(i-1,j,k,3) / ux(i-1,j,k,1)
		 uijp1k = ux(i,j+1,k,2) / ux(i,j+1,k,1)
		 uijm1k = ux(i,j-1,k,2) / ux(i,j-1,k,1)
		 omega(i,j,k) = (vip1jk - vim1jk)/dx - 
     1	                         (uijp1k - uijm1k)/dy

	      enddo
	   enddo
	enddo
	
c       ----circulation
c
 	do k = 1, NZlocal
 	   do j = 1, NYlocal
 	      do i = 1, NXlocal
 		 circ = circ + omega(i,j,k) * dx*dy
		 if (omega(i,j,k) .lt. 0.D0) circn = circn + 
     1	            omega(i,j,k) * dx*dy
		 if ( (xi(i,j,k).gt.-0.999).and.(xi(i,j,k) .lt. 0.999)
     1	      		) then
 		    circint = circint + omega(i,j,k) * dx*dy
 		    if ( omega(i,j,k) .lt. 0.D0) circintn = circintn + 
     1	                                     omega(i,j,k) * dx*dy
 		 endif
 	      enddo
 	   enddo
 	enddo
	omega_max_local=-1.E16
	omega_min_local=1.E16
	do k = 1, NZlocal
	   do j = 1, NYlocal
	      do i = 1, NXlocal
		if(omega_max_local.lt.omega(i,j,k))
     1	           omega_max_local=omega(i,j,k)
		if(omega_min_local.gt.omega(i,j,k))
     1	           omega_min_local=omega(i,j,k)
	      enddo
	   enddo
	enddo
	
#ifdef PARALLEL
        call MPI_Allreduce(omega_max_local, omega_max, 1, 
     &       MPI_DOUBLE_PRECISION,
     &       MPI_MAX, comm3D, ierr)
        call MPI_Allreduce(omega_min_local, omega_min, 1, 
     &       MPI_DOUBLE_PRECISION,
     &       MPI_MIN, comm3D, ierr)
#else   
	omega_max=omega_max_local
	omega_min=omega_min_local
#endif  
	if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
	   write(61,*)ttot, omega_max,omega_min
	endif

c	---- Real circulations
c       ---- Global circulations.
c$$$
#ifdef PARALLEL
	circ_global = 0.D0
	circn_global = 0.D0
	circint_global = 0.D0
	circintn_global = 0.D0
	
	call MPI_Allreduce(circ, circ_global, 1, MPI_DOUBLE_PRECISION,
     &       MPI_SUM, comm3D, ierr)
        call MPI_Allreduce(circn, circn_global, 1,
     1       MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
	call MPI_Allreduce(circint, circint_global, 1,
     1       MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
	call MPI_Allreduce(circintn, circintn_global, 1,
     1       MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
#else
	circ_global = circ
	circn_global = circn
	circint_global = circint 
	circintn_global = circintn
#endif
        if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
	   write(11,111)ttot, circ_global, circn_global,
     1	         circint_global, circintn_global
 111	   format(1x, 5(E13.6, 1x))
        endif
#else
	real*8 omegax(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI)
	real*8 omegay(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI)
	real*8 omegaz(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI)
	real*8 omega(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI)
	real*8 xi(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI)
	real*8 a(NXlsize, NYlsize, NZlsize)

	real*8 circx, circnx, circy, circny, circz, circnz
	real*8 circx_global, circnx_global, circy_global, circny_global,
     1         circz_global, circnz_global

	real*8 circintx, circintnx, circinty, circintny, circintz,
     1         circintnz
	real*8 circintx_global, circintnx_global, circinty_global,
     1         circintny_global, circintz_global,circintnz_global

	real*8 circdomain, circint, circndomain, circnint
	real*8 circdomain_global, circint_global, circndomain_global,
     1         circnint_global

	real*8 omega_max, omega_max_local
	real*8 omegax_max, omegax_max_local
	real*8 omegay_max, omegay_max_local
	real*8 omegaz_max, omegaz_max_local
	real*8 omegax_min, omegax_min_local
	real*8 omegay_min, omegay_min_local
	real*8 omegaz_min, omegaz_min_local
	real*8 gammafun, gamma, u, v, w, p 	       
	real*8  wijp1k, wijm1k, vijkp1, vijkm1 
	real*8  wip1jk,wim1jk,uijkp1,uijkm1 
	real*8  vip1jk, vim1jk, uijp1k,uijm1k
	integer i, j, k, iter, ndump, index

	common/time/dt, ttot
	real*8 dt, ttot

	circx = 0.D0
	circnx = 0.D0
	circintx = 0.D0
	circintnx = 0.D0

	circy = 0.D0
	circny = 0.D0
	circinty = 0.D0
	circintny = 0.D0

	circz = 0.D0
	circnz = 0.D0
	circintz = 0.D0
	circintnz = 0.D0

c	---- Determine xi (interface tracker)

	do k = 1, NZlocal
	   do j = 1, NYlocal
	      do i = 1, NXlocal
		 xi(i,j,k) = ux(i,j,k,6) / ux(i,j,k,1)
	      enddo
	   enddo
	enddo

c       ---- Determine vorticity

	do k = 1, NZlocal
	   do j = 1, NYlocal
	      do i = 1, NXlocal
		 
		 wijp1k = ux(i,j+1,k,4) / ux(i,j+1,k,1)
		 wijm1k = ux(i,j-1,k,4) / ux(i,j-1,k,1)
		 vijkp1 = ux(i,j,k+1,3) / ux(i,j,k+1,1)
		 vijkm1 = ux(i,j,k-1,3) / ux(i,j,k-1,1)
		 omegax(i,j,k) = (wijp1k - wijm1k)/dy - 
     1	                         (vijkp1 - vijkm1)/dz

		 wip1jk = ux(i+1,j,k,4) / ux(i+1,j,k,1)
		 wim1jk = ux(i-1,j,k,4) / ux(i-1,j,k,1)
		 uijkp1 = ux(i,j,k+1,2) / ux(i,j,k+1,1)
		 uijkm1 = ux(i,j,k-1,2) / ux(i,j,k-1,1)
		 omegay(i,j,k) = - ( (wip1jk - wim1jk)/dx - 
     1	                             (uijkp1-uijkm1)/dz )

		 vip1jk = ux(i+1,j,k,3) / ux(i+1,j,k,1)
		 vim1jk = ux(i-1,j,k,3) / ux(i-1,j,k,1)
		 uijp1k = ux(i,j+1,k,2) / ux(i,j+1,k,1)
		 uijm1k = ux(i,j-1,k,2) / ux(i,j-1,k,1)
		 omegaz(i,j,k) = (vip1jk - vim1jk)/dx - 
     1	                         (uijp1k - uijm1k)/dy

	      enddo
	   enddo
	enddo
	
c       ---- X circulation
c$$$ 
c$$$ 	do k = 1, NZlocal
c$$$ 	   do j = 1, NYlocal
c$$$ 	      do i = 1, NXlocal
c$$$ 		 dy = y(i,j,k) - y(i,j-1,k)
c$$$ 		 dz = z(i,j,k) - z(i,j,k-1)
c$$$ 		 dA = dy * dz 
c$$$ 		 circx = circx + omegax(i,j,k) * dA
c$$$ 		 if (omegax(i,j,k) .lt. 0.D0) circnx = circnx + 
c$$$     1	                                     omegax(i,j,k) * dA
c$$$ 		 if ( (xi(i,j,k) .gt. 1e-3) .and. (xi(i,j,k) .lt. 0.999)
c$$$     1	      		) then
c$$$ 		    circintx = circintx + omegax(i,j,k) * dA
c$$$ 		    if ( omegax(i,j,k) .lt. 0.D0) circintnx = circintnx + 
c$$$     1	                                     omegax(i,j,k) * dA
c$$$ 		 endif
c$$$ 	      enddo
c$$$ 	   enddo
c$$$ 	enddo
c$$$ 
c$$$c       ---- Y circulation
c$$$ 
c$$$ 	do k = 1, NZlocal
c$$$ 	   do j = 1, NYlocal
c$$$ 	      do i = 1, NXlocal
c$$$ 		 dx = x(i,j,k) - x(i-1,j,k)
c$$$ 		 dz = z(i,j,k) - z(i,j,k-1)
c$$$ 		 dA = dx * dz
c$$$ 		 circy = circy + omegay(i,j,k) * dA
c$$$ 		 if (omegay(i,j,k) .lt. 0.D0) circny = circny + 
c$$$     1	                             omegay(i,j,k) * dA
c$$$ 		 if ( (xi(i,j,k) .gt. 1e-3) .and. (xi(i,j,k) .lt. 0.999)
c$$$     1	      	     ) then
c$$$ 		    circinty = circinty + omegay(i,j,k) * dA
c$$$ 		    if ( omegay(i,j,k) .lt. 0.D0) circintny = circintny + 
c$$$     1	                             omegay(i,j,k) * dA
c$$$ 		 endif
c$$$ 	      enddo
c$$$ 	   enddo
c$$$ 	enddo
c$$$
c$$$c       ---- Z circulation
c$$$
c$$$	do k = 1, NZlocal
c$$$	   do j = 1, NYlocal
c$$$	      do i = 1, NXlocal
c$$$		 dx = x(i,j,k) - x(i-1,j,k)
c$$$		 dy = y(i,j,k) - y(i,j-1,k)
c$$$		 dA = dx * dy
c$$$		 circz = circz + omegaz(i,j,k) * dA
c$$$		 if (omegaz(i,j,k) .lt. 0.D0) circnz = circnz + 
c$$$     1	                             omegaz(i,j,k) * dA
c$$$		 if ( (xi(i,j,k) .gt. 1e-3) .and. (xi(i,j,k) .lt. 0.999)
c$$$     1	      	    ) then
c$$$		    circintz = circintz + omegaz(i,j,k) * dA
c$$$		    if ( omegaz(i,j,k) .lt. 0.D0) circintnz = circintnz + 
c$$$     1	                             omegaz(i,j,k) * dA
c$$$		 endif
c$$$	      enddo
c$$$	   enddo
c$$$	enddo
c       Vorticity Magnitude
	do k = 1, NZlocal
	   do j = 1, NYlocal
	      do i = 1, NXlocal
		 omega(i,j,k)=dsqrt(omegax(i,j,k)**2.0+
     1	               omegay(i,j,k)**2.0+omegaz(i,j,k)**2.0)
	      enddo
	   enddo
	enddo
	omega_max_local=-1.E16
	omegax_max_local=-1.E16
	omegay_max_local=-1.E16
	omegaz_max_local=-1.E16
	omegax_min_local=1.E16
	omegay_min_local=1.E16
	omegaz_min_local=1.E16
	do k = 1, NZlocal
	   do j = 1, NYlocal
	      do i = 1, NXlocal
		if(omega_max_local.lt.omega(i,j,k))
     1	           omega_max_local=omega(i,j,k)
		if(omegax_max_local.lt.omegax(i,j,k))
     1	           omegax_max_local=omegax(i,j,k)
		if(omegay_max_local.lt.omegay(i,j,k))
     1	           omegay_max_local=omegay(i,j,k)
		if(omegaz_max_local.lt.omegaz(i,j,k))
     1	           omegaz_max_local=omegaz(i,j,k)
		if(omegax_min_local.gt.omegax(i,j,k))
     1	           omegax_min_local=omegax(i,j,k)
		if(omegay_min_local.gt.omegay(i,j,k))
     1	           omegay_min_local=omegay(i,j,k)
		if(omegaz_min_local.gt.omegaz(i,j,k))
     1	           omegaz_min_local=omegaz(i,j,k)
	      enddo
	   enddo
	enddo
	
#ifdef PARALLEL
        call MPI_Allreduce(omega_max_local, omega_max, 1, 
     &       MPI_DOUBLE_PRECISION,
     &       MPI_MAX, comm3D, ierr)
        call MPI_Allreduce(omegax_max_local, omegax_max, 1, 
     &       MPI_DOUBLE_PRECISION,
     &       MPI_MAX, comm3D, ierr)
        call MPI_Allreduce(omegay_max_local, omegay_max, 1, 
     &       MPI_DOUBLE_PRECISION,
     &       MPI_MAX, comm3D, ierr)
        call MPI_Allreduce(omegaz_max_local, omegaz_max, 1, 
     &       MPI_DOUBLE_PRECISION,
     &       MPI_MAX, comm3D, ierr)
        call MPI_Allreduce(omegax_min_local, omegax_min, 1, 
     &       MPI_DOUBLE_PRECISION,
     &       MPI_MIN, comm3D, ierr)
        call MPI_Allreduce(omegay_min_local, omegay_min, 1, 
     &       MPI_DOUBLE_PRECISION,
     &       MPI_MIN, comm3D, ierr)
        call MPI_Allreduce(omegaz_min_local, omegaz_min, 1, 
     &       MPI_DOUBLE_PRECISION,
     &       MPI_MIN, comm3D, ierr)
#else   
	omega_max=omega_max_local
	omegax_max=omegax_max_local
	omegay_max=omegay_max_local
	omegaz_max=omegaz_max_local
	omegax_min=omegax_min_local
	omegay_min=omegay_min_local
	omegaz_min=omegaz_min_local
#endif  
	if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
	write(61,*)ttot, omega_max,omegax_max,omegay_max,
     1	         omegaz_max
	write(61,*)ttot, omegax_min,omegay_min,
     1	         omegaz_min
	endif

c	---- Real circulations

c$$$	circdomain = circx + circy + circz
c$$$	circndomain = circnx + circny + circnz
c$$$	circint = circintx + circinty + circintz
c$$$	circnint = circintnx + circintny + circintnz
c$$$
c$$$c       ---- Global circulations.
c$$$
c$$$#ifdef PARALLEL
c$$$
c$$$c         ---- X circulation
c$$$
c$$$	circx_global = 0.D0
c$$$	circnx_global = 0.D0
c$$$	circintx_global = 0.D0
c$$$	circintnx_global = 0.D0
c$$$	
c$$$	call MPI_Allreduce(circx, circx_global, 1, MPI_DOUBLE_PRECISION,
c$$$     &       MPI_SUM, comm3D, ierr)
c$$$	call MPI_Allreduce(circnx, circnx_global, 1,
c$$$     1       MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
c$$$	call MPI_Allreduce(circintx, circintx_global, 1,
c$$$     1       MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
c$$$	call MPI_Allreduce(circintnx, circintnx_global, 1,
c$$$     1       MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
c$$$
c$$$c         ---- Y circulation
c$$$
c$$$	circy_global = 0.D0
c$$$	circny_global = 0.D0
c$$$	circinty_global = 0.D0
c$$$	circintny_global = 0.D0
c$$$	
c$$$	call MPI_Allreduce(circy, circy_global, 1, MPI_DOUBLE_PRECISION,
c$$$     &       MPI_SUM, comm3D, ierr)
c$$$	call MPI_Allreduce(circny, circny_global, 1,
c$$$     1       MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
c$$$	call MPI_Allreduce(circinty, circinty_global, 1,
c$$$     1       MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
c$$$	call MPI_Allreduce(circintny, circintny_global, 1,
c$$$     1       MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
c$$$
c$$$c         ---- Z circulation
c$$$
c$$$	circz_global = 0.D0
c$$$	circnz_global = 0.D0
c$$$	circintz_global = 0.D0
c$$$	circintnz_global = 0.D0
c$$$	
c$$$	call MPI_Allreduce(circz, circz_global, 1, MPI_DOUBLE_PRECISION,
c$$$     &       MPI_SUM, comm3D, ierr)
c$$$	call MPI_Allreduce(circnz, circnz_global, 1,
c$$$     1       MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
c$$$	call MPI_Allreduce(circintz, circintz_global, 1,
c$$$     1       MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
c$$$	call MPI_Allreduce(circintnz, circintnz_global, 1,
c$$$     1       MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
c$$$
c$$$c         ---- Total circulation; sum over all subdomains
c$$$
c$$$	call MPI_Allreduce(circdomain, circdomain_global, 1,
c$$$     1       MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
c$$$	call MPI_Allreduce(circndomain, circndomain_global, 1,
c$$$     1       MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
c$$$	call MPI_Allreduce(circint, circint_global, 1,
c$$$     1       MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
c$$$	call MPI_Allreduce(circnint, circnint_global, 1,
c$$$     1       MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
c$$$	
c$$$#else
c$$$
c$$$	circx_global = circx
c$$$	circnx_global = circnx
c$$$	circintx_global = circintx 
c$$$	circintnx_global = circintnx 
c$$$
c$$$	circy_global = circy
c$$$	circny_global = circny
c$$$	circinty_global = circinty
c$$$	circintny_global = circintny
c$$$
c$$$	circz_global = circz
c$$$	circnz_global = circnz
c$$$	circintz_global = circintz
c$$$	circintnz_global = circintnz
c$$$
c$$$	circdomain_global = circdomain
c$$$	circndomain_global = circndomain
c$$$	circint_global = circint
c$$$	circnint_global = circnint
c$$$			
c$$$#endif
c$$$
c$$$        if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
c$$$	   write(21,111)ttot, circdomain_global, circndomain_global,
c$$$     1	         circint_global, circnint_global
c$$$  111	   format(1x, 5(E13.6, 1x))
c$$$
c$$$	   write(*,112)ttot, circdomain_global, circint_global
c$$$  112	   format(1x,' Circulation -> ', 3(E13.6, 2x))
c$$$        endif
#endif
	return
	end






