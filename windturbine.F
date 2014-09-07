c----------------------------------------------------------------------
c     Wan Cheng 
c     November 2012
c----------------------------------------------------------------------
c
c
c
c
c----------------------------------------------------------------------
c     Subroutine to set the wind turbine model
c----------------------------------------------------------------------
c
      module wturbine
      save
       integer:: nturbine, npermesh
       integer,allocatable:: ipltxyz(:,:)
       double precision, allocatable:: wflag(:,:,:)
       double precision, allocatable:: wforce_old(:,:,:)
      end module wturbine
c----------------------------------------------------------------------
c     allocate values in module
c----------------------------------------------------------------------
      subroutine allocate_turbine
      use mesh
      use mesh_common
      use wturbine
c
      nturbine = 4
      npermesh = 10
c
      allocate(ipltxyz(6,nturbine))    
      allocate(wflag(nturbine,nylocal+1,nzlocal+1))
      allocate(wforce_old(nxlocal+1,nylocal+1,nzlocal+1))
c
      ipltxyz = 0
      wflag = 0.d0
      wforce_old = 0.d0
c
      return
      end    
c----------------------------------------------------------------------
c     locate turbine positions    
c----------------------------------------------------------------------
      subroutine locate_turbine
      use mesh
      use mesh_common
      use wturbine
#ifdef PARALLEL
      use mpistuff
#endif
c
      double precision:: xd(nturbine), yd(nturbine), zd(nturbine)
      double precision:: ydlocal,zdlocal, rcir
c
      double precision:: area(npermesh,npermesh)
      double precision:: yt(npermesh)
      double precision:: zt(npermesh)
c
      double precision:: xll,xlr,yll,ylr,zll,zlr
c
      integer:: g_i(nturbine), gt_i(nturbine)
      integer:: i,j,k,l,n,ly,lz,lno, lturbine, ixp, ix, ixgflag
      integer:: nxtmp, nytmp
      integer:: xtur, ytur
c
      ytur = 4
      xtur = nturbine/ytur
c
      rcir = 0.2d0
      zd = 1.0d0
      do l=1, nturbine
        xd(l) = 3.d0 + int((l-1)/4)* 3.d0
        yd(l) = 0.375d0 + mod(l-1,4)* 0.75d0
    
c        write(*,*) l, xd(l), yd(l)
      enddo
c
      xll = xc(1)-0.5d0*dx
      xlr = xc(nxlocal)+0.5d0*dx
      yll = yc(1)-0.5d0*dy
      ylr = yc(nylocal)+0.5d0*dy
      zll = zc(1)-0.5d0*dz
      zlr = zc(nzlocal)+0.5d0*dz
      do l=1, nturbine, 1
        if(xd(l).gt.xll .and. xd(l).lt.xlr) then
          ipltxyz(1,l)=iprocx
          do i=1, nxlocal+1, 1
            if(dabs(xd(l)-(xc(i)-0.5d0*dx)).le.0.5d0*dx) ipltxyz(4,l)=i
          enddo  
        endif
        if(yd(l).gt.yll .and. yd(l).lt.ylr) then 
          ipltxyz(2,l)=iprocy
          do j=1, nylocal+1, 1
            if(dabs(yd(l)-(yc(j)-0.5d0*dy)).le.0.5d0*dy) ipltxyz(5,l)=j
          enddo
        endif
        if(zd(l).gt.zll .and. zd(l).lt.zlr) then
          ipltxyz(3,l)=iprocz
          do k=1, nzlocal+1, 1
            if(dabs(zd(l)-(zc(k)-0.5d0*dz)).le.0.5d0*dz) ipltxyz(6,l)=k
          enddo
        endif
c         write(*,*) 'dddddddddddd', l, ipltxyz(1,l),
c     &        ipltxyz(2,l),ipltxyz(3,l),ipltxyz(4,l),ipltxyz(5,l),
c     &        ipltxyz(6,l),iprocx,iprocy,iprocz
      enddo
c
#ifdef PARALLEL
      do n=1, 6, 1
        gt_i(1:nturbine) = ipltxyz(n,1:nturbine)
        g_i = 0
        call MPI_Allreduce(gt_i, g_i, nturbine,
     &          MPI_INTEGER, MPI_SUM, comm3d, ierr)
        call ErrorHandler(ierr,ERROR_ALLREDUCE)
        ipltxyz(n,1:nturbine) = g_i(1:nturbine) 
      enddo
      ipltxyz(1,1:nturbine) = ipltxyz(1,1:nturbine)/yprocs/zprocs
      ipltxyz(2,1:nturbine) = ipltxyz(2,1:nturbine)/xprocs/zprocs
      ipltxyz(3,1:nturbine) = ipltxyz(3,1:nturbine)/xprocs/yprocs
      ipltxyz(4,1:nturbine) = ipltxyz(4,1:nturbine)/yprocs/zprocs
      ipltxyz(5,1:nturbine) = ipltxyz(5,1:nturbine)/xprocs/zprocs
      ipltxyz(6,1:nturbine) = ipltxyz(6,1:nturbine)/xprocs/yprocs      
#endif

      do n=1, 6, 1
        do l=1, nturbine, 1
          if(ipltxyz(n,l).le.0) then
            write(*,*) 'bad wind turbine location', l, n,ipltxyz(n,l)
            stop
          endif
        enddo
      enddo
c
ccccccccccccccccccccccccccccccccccccc
      do l=1, nturbine, 1
c        write(*,*) l, ipltxyz(1,l),ipltxyz(4,l), iprocx
      enddo
c      stop

cccccccccccccccccccccccccccccccccccccccccc     

c     get area ratio in y-z plane
      do l=1, nturbine, 1
        if( iprocx.eq.ipltxyz(1,l) ) then
        i=ipltxyz(4,l)
c        write(*,*) l, i
        do j=1, nylocal, 1
          do k=1, nzlocal, 1
            area = 0.d0
            do n=1, npermesh
              yt(n)= yc(j)-dy/2.d0+dy/npermesh*n-dy/npermesh/2.d0
              zt(n)= zc(k)-dz/2.d0+dz/npermesh*n-dz/npermesh/2.d0
c              write(*,*) j, k, n, yt(n), zt(n)
            enddo
            do ly=1,npermesh,1
              do lz=1,npermesh,1
              if( dsqrt((yt(ly)-yd(l))**2.d0+(zt(lz)-zd(l))**2.d0)
     &             .lt.rcir ) then
                area(ly,lz) = 1.d0
              else
                area(ly,lz) = 0.d0
              endif
c                 write(*,*) ly, lz, area(ly,lz)
              wflag(l,j,k)=wflag(l,j,k)+area(ly,lz)
              enddo
            enddo
            wflag(l,j,k)=wflag(l,j,k)/npermesh/npermesh
c            if(   wflag(l,j,k).gt. 0.d0) then
c                write(*,*) l, j,k, wflag(l,j,k)
c            endif
          enddo
        enddo
        endif
      enddo
c
      return
      end
c
c----------------------------------------------------------------------
c     get wind force term for streamwise momentum equation
c----------------------------------------------------------------------
      subroutine windforce(uxn,wforce,timestep,substep)
      use mesh
      use mesh_common
      use properties
      use wturbine
      use RKCoefficients
      use stats
      use iounits
#ifdef PARALLEL
      use mpistuff
#endif
c
      implicit none
      double precision:: uxn(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision:: wforce(nxlocal+1,nylocal+1,nzlocal+1,nvar)
      integer:: substep,timestep
c
      double precision:: wflagt(nxlocal,nylocal,nzlocal,nturbine)
      double precision:: ud(nturbine)
      double precision:: wt_t, wt_n, wt_o, w_area, ud_t, wtot1, wtot2
      double precision:: coeff
c
      integer:: i,j,k,l
      integer:: npoint(nturbine), gnpoint(nturbine)
c
      coeff = 4.d0/3.d0
c

c     get average velocity on the disk
      ud = 0.d0
      npoint = 0
      do l=1, nturbine, 1
        if( iprocx.eq.ipltxyz(1,l) ) then
        i=ipltxyz(4,l)
          do j=1, nylocal, 1
            do k=1, nzlocal, 1
              if( wflag(l,j,k).gt.0.d0 ) then
                ud(l) = ud(l) + wflag(l,j,k)*st_umean(i,k,1)
                npoint(l) = npoint(l) + 1
c                 if(l.eq.1) write(*,*) npoint(l), wflag(l,j,k)

              endif
            enddo
          enddo   
        endif
      enddo
#ifdef PARALLEL
      call MPI_Allreduce(npoint, gnpoint, nturbine,
     &          MPI_INTEGER, MPI_SUM, comm3d, ierr)
      call ErrorHandler(ierr,ERROR_ALLREDUCE)            
#else
      gnpoint = npoint
#endif
      do l=1, nturbine
        ud(l) = ud(l) / gnpoint(l)
c        write(*,*) l, ud(l), iprocx
      enddo
c
      wtot1 = 0.d0
c     calculate wforce
      do i=1, nxlocal, 1
        do j =1, nylocal, 1
          do k = 1, nzlocal, 1
            wt_n = 0.d0
            wt_t = 0.d0
            wt_o = 0.d0
c           calculate
            w_area = 0.d0
            ud_t = 0.d0
            do l=1, nturbine, 1
              if(iprocx.eq.ipltxyz(1,l) .and. i.eq.ipltxyz(4,l)) then
                w_area = w_area + wflag(l,j,k)
                if(wflag(l,j,k).gt.0.d0)  ud_t = ud(l)
              endif
            enddo
c
            wt_t=0.5d0*coeff*ud_t**2.d0*w_area/dx
c
            wt_o=wforce_old(i,j,k)
            if(substep.eq.1) then
              wt_n = gamma(1)*wt_t
              wt_o = wt_t
            endif
            if(substep.eq.2) then
              wt_n = gamma(2)*wt_t + zeta(1)*wt_o
              wt_o = wt_t
            endif
            if(substep.eq.3) then
              wt_n = gamma(3)*wt_t + zeta(2)*wt_o
              wt_o = 1.d10
            endif
            wforce(i,j,k,1) =  wt_n
            wforce_old(i,j,k) = wt_o
            wtot1=wtot1+wforce(i,j,k,1)*uxn(i,j,k,1)*dy*dz
          enddo         
        enddo
      enddo
c     force load---
#ifdef PARALLEL
      wtot2 = 0.d0
        call MPI_Allreduce(wtot1, wtot2, 1,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
        call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
      wtot2 = wtot1
#endif
            write(ihis,*) timestep, substep, wtot2
c
      return
      end
c----------------------------------------------------------------------







