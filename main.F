c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2010
c     King Abdullah University of Science and Technology
c     All Rights Reserved
c-----------------------------------------------------------------
c     $Log: main.F,v $
c     Revision 1.5  2011/10/19 15:14:49  samtanr
c     Added 2D option via TWO_D defines.
c
c     Revision 1.4  2011/10/18 13:01:39  samtanr
c     Added abort routine.
c
c     Revision 1.3  2011/10/08 13:37:51  samtanr
c     Changed no of processors.
c
c     Revision 1.2  2011/10/07 13:17:36  samtanr
c     Changed nx=32 to 128.
c
c     Revision 1.1  2011/07/25 11:32:01  samtanr
c     Original source.
c
c     Revision 1.1.1.1  2010/12/31 12:44:32  samtanr
c     Wind code: 2d, 3d and low Mach number
c
c-----------------------------------------------------------------
c     Code 
c-----------------------------------------------------------------
c
      program Wind
      use mesh_parms
      use mesh_common
      use iounits
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      logical:: error
      integer:: num_nodes
      integer:: ipx, ipy, ipz !, i
      integer,allocatable:: my_iproc_t(:)
c
      integer:: ndims
      logical:: periods(3), reorder
      integer:: NProcAxis(3)
      integer, external:: iprocnum
c
#ifdef DYNAMIC
      namelist/meshsize/nx,ny,nz
      namelist/meshprocs/xprocs,yprocs,zprocs
c
c      open(imesh,file='mesh.inp',form='formatted')
c      read(imesh,meshsize)
c      read(imesh,meshprocs)
c      close(imesh)
      nx=336
#ifndef TWO_D
      ny=32
#else
      ny=1
#endif
      nz=96
#ifdef PARALLEL
      xprocs=21
      yprocs=2
      zprocs=6
#else
      xprocs=1
      yprocs=1
      zprocs=1
#endif
c
#endif
c     Initialize the communicator, call it Comm3D
      ndims = 3                        
#ifdef DYNAMIC
      NPROCS=XPROCS*YPROCS*ZPROCS
      nxlocal=nx/xprocs
      nylocal=ny/yprocs
      nzlocal=nz/zprocs
      nxlsize=nxlocal
      nylsize=nylocal
      nzlsize=nzlocal
#else
      nxlocal=nxlsize
      nylocal=nylsize
      nzlocal=nzlsize
#endif
#ifdef PARALLEL
      NProcAxis(1) = XPROCS            
      NProcAxis(2) = YPROCS            
      NProcAxis(3) = ZPROCS
c      periods(1) = .true.             
      periods(1) = .false.             
#ifdef YPERIODIC
      periods(2) = .true.
#else
      periods(2) = .false.
#endif
c      periods(3) = .true.
      periods(3) = .false.
      reorder = .false.
      master = 0

      call MPI_init(ierr)          
      call MPI_cart_create( MPI_COMM_WORLD, ndims, NProcAxis, periods,
     &       reorder, comm3D, ierr)
      call MPI_comm_rank( comm3D, my_id, ierr)
      call MPI_comm_size( comm3D, num_nodes, ierr)
        
c      write(6,*) 'Started proc #', my_id, num_nodes

      if (num_nodes .lt. NPROCS) then
        write(istderr,*) 'SORRY: PROGRAM ONLY WORKS FOR ', NPROCS,
     &          ' PROCESSES.'
        call MPI_Finalize(ierr)
        stop
      endif
C     check for invalid combinations of dimensions and
C     process-grid dimensions
      error = .FALSE.
      if (NX - (XPROCS-1)*NXlsize .le. 0) then
        print*, 'invalid combination of NX, XPROCS'
        error = .TRUE.
      endif
      if (NY - (YPROCS-1)*NYlsize .le. 0) then
        print*, 'invalid combination of NY, YPROCS'
        error = .TRUE.
      endif
      if (NZ - (ZPROCS-1)*NZlsize .le. 0) then
        print*, 'invalid combination of NZ, ZPROCS'
        error = .TRUE.
      endif
      if (error) STOP
c
c main code
      if (my_id .lt. NPROCS) then
        call pxpypz (my_id+1, ipx, ipy, ipz)
c        write(6,*) 'PROGRAM START ...',iprocx,iprocy,iprocz
        iproc_idx = iprocnum(ipx, ipy, ipz)  ! process name in code
c       establish my_iproc array, used in recycling
        allocate( my_iproc(1:nprocs) )
        allocate( my_iproc_t(1:nprocs) )
        my_iproc = 0
        my_iproc_t=0
        my_iproc_t(iproc_idx)=my_id
        call MPI_Allreduce(my_iproc_t, my_iproc, nprocs,
     &          MPI_INTEGER, MPI_SUM, comm3D, ierr)
        Call ErrorHandler(ierr,ERROR_ALLREDUCE)
c---------------------------------------
        ! main code here
        call mesh_grid_proc(ipx,ipy,ipz)
c---------------------------------------
c        write(6,*) 'PROGRAM END',iprocx,iprocy,iprocz
        call MPI_Finalize(ierr)
      endif
#else
      ipx=1
      ipy=1
      ipz=1
c---------------------------------------
      call mesh_grid_proc(ipx,ipy,ipz)
c---------------------------------------
#endif
c
      stop
      end
c
c----------------------------------------------------------------------
c        main program for grid process
c----------------------------------------------------------------------
      subroutine mesh_grid_proc(a_iprocx, a_iprocy, a_iprocz)
      use mesh_parms
      use mesh_common
c      
      implicit none   
      integer:: a_iprocx, a_iprocy, a_iprocz
c
      iprocx = a_iprocx
      iprocy = a_iprocy
      iprocz = a_iprocz

      NXlocal = NXlsize
      IF (iprocx .EQ. XPROCS) NXlocal = NX-(XPROCS-1)*NXlsize
      NYlocal = NYlsize
      IF (iprocy .EQ. YPROCS) NYlocal = NY-(YPROCS-1)*NYlsize
      NZlocal = NZlsize
      IF (iprocz .EQ. ZPROCS) NZlocal = NZ-(ZPROCS-1)*NZlsize
c
      call WindMain
c
      return
      end

c
c----------------------------------------------------------------------
c     end of code
c----------------------------------------------------------------------
      subroutine Abortend
      use IOUnits
c      
      implicit none
      integer:: ierr
      write(ihis,*) 'ABORTING'
#ifdef PARALLEL
      call MPI_Finalize(ierr)
#endif
      return
      end
c
c-----------------------------------------------------------------------

