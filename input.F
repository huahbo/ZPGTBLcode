C -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2011
c     All Rights Reserved
c-----------------------------------------------------------------
c     $Log$
c     Read data from files at each time step
c-----------------------------------------------------------------

c
c-----------------------------------------------------------------------


#ifdef PARALLELHDF5
      subroutine ReadParallelHDF5File3D(ux,index)
      use mesh
      use mpistuff
      use mesh_common
      use iounits
c
      implicit none
      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI, nvar)
c      double precision:: ux(nxlsize, nylsize,nzlsize,nvar)
      double precision:: tmp(nxlsize, nylsize,nzlsize)
c     Local variables
      integer:: ndims,ierr1,dims(3),offset(3), count(3)
c     
      integer i,j,k,lloc,lg,ipv
      integer index,iunit, data_type
      integer nbytes, itmp,  nelements
      character*80 outfile,filename,dataname
      integer:: dbfile
c      integer,parameter::ifld=1
c     

      outfile = 'Init.h5'
      iunit=101
      filename=trim(outfile)//char(0)
 
c      write(6,*) 'Reading hdf5 file',filename
      call open_parhdf5_file_rdonly(iunit,comm3D,filename)
c
      ndims=3
      dims(1)=nx
      dims(2)=ny
      dims(3)=nz
      count(1)=nxlocal
      count(2)=nylocal
      count(3)=nzlocal
      offset(1)=(iprocx-1)*nxlocal
      offset(2)=(iprocy-1)*nylocal
      offset(3)=(iprocz-1)*nzlocal
c
      dataname='ux'//char(0)
      call read_parhdf5_f(iunit,tmp, data_type,ndims,dims,
     &     count,offset,comm3D,trim(dataname))
         ux(1:nxlocal,1:nylocal,1:nzlocal,1)
     & =tmp(1:nxlocal,1:nylocal,1:nzlocal)
c
      dataname='uy'//char(0)
      call read_parhdf5_f(iunit,tmp, data_type,ndims,dims,
     &     count,offset,comm3D,trim(dataname))
         ux(1:nxlocal,1:nylocal,1:nzlocal,2)
     & =tmp(1:nxlocal,1:nylocal,1:nzlocal)
c
      dataname='uz'//char(0)
      call read_parhdf5_f(iunit,tmp, data_type,ndims,dims,
     &     count,offset,comm3D,trim(dataname))
         ux(1:nxlocal,1:nylocal,1:nzlocal,3)
     & =tmp(1:nxlocal,1:nylocal,1:nzlocal)c        write(*,*) tmp(10,10,nzlocal-1), tmp(10,10,nzlocal)
c        write(*,*) ux(10,10,nzlocal-1,3),ux(10,10,nzlocal,3)
c
c
      call MPI_BARRIER(comm3D,ierr)
      call close_hdf5_file_rdonly(iunit)

      end
c
c-----------------------------------------------------------------------
#endif


