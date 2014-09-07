c       Ravi Samtaney
c       Copyright 2010
c-------------------------------------------------------------------------

c-----------------------------------------------------------------------
      subroutine WriteCheckPoint(ux,dt,ttot,index)
      use mesh
      use mesh_common

        double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI, nvar)
        double precision:: dt, ttot
        integer i,j,k
       integer index
        integer nbytes, offset,itmp, ierr, nelements
        integer:: FileNumber
       character*50 infile, outfile
c
        FileNumber = iprocx*10000 + iprocy*100 + iprocz
        write(inFile, 1000) float(FileNumber) / 1000000.0,
     &              float(index)/1000000.0
c1000    format(5Hcheck,f7.6,f7.6)
1000    format('check',f7.6,f7.6)
        open(21,file=infile,form='unformatted')
c
        write(21) dt,ttot
        write(21)ux
c        write(21) xs
c        write(21) ys
c        write(21) zs
c        write(21) us
c        write(21) vs
c        write(21) ws
        close(21)
      end
c----------------------------------------------------------------
      subroutine ReadCheckPoint(ux,dt,ttot,index)
      use mesh
      use mesh_common
        double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI, nvar)
        double precision:: dt, ttot
        integer i,j,k
       integer index
        integer nbytes, offset,itmp, ierr, nelements
        integer:: FileNumber
       character*50 infile, outfile
c
        FileNumber = iprocx*10000 + iprocy*100 + iprocz
        write(inFile, 1000) float(FileNumber) / 1000000.0,
     &              float(index)/1000000.0
c1000    format(5Hcheck,f7.6,f7.6)
1000    format('check',f7.6,f7.6)
        open(21,file=infile,form='unformatted')
c
        read(21) dt,ttot
        read(21)ux
c        read(21) xs
c        read(21) ys
c        read(21) zs
c        read(21) us
c        read(21) vs
c        read(21) ws
        close(21)
      end
c----------------------------------------------------------------
#ifdef HDF5
      subroutine WriteHDF5File(ux,index)
      use mesh
      use mesh_common
      use properties
      use iounits
c
      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI, nvar)
      double precision:: tmp(nxlsize, nylsize,nzlsize)
      double precision:: divu(nxlsize, nylsize,nzlsize)
      double precision:: vort(nxlsize, nylsize,nzlsize)
c     Local variables
      double precision:: ttot,dt
      real:: xn(nxlocal+1),yn(nylocal+1),zn(nzlocal+1)
      integer:: ndims,ierr1,dims(3)
c     
      integer i,j,k
      integer index,iunit
      integer nbytes, offset,itmp, ierr, nelements
      character*50 outfile,filename,dataname
      integer:: dbfile
c      integer,parameter::ifld=1
c     

c     File name for data 
      write(outfile,1000) 
     &     float(index)/10000000.0
 1000 format('output',f8.7,'.h5')
c
      iunit=101
      filename=trim(outfile)//char(0)
c      write(6,*) 'Writing hdf5 file',filename
      call open_hdf5_file(iunit,filename)
c
      ndims=3
      dims(1)=nxlocal
      dims(2)=nylocal
      dims(3)=nzlocal
c
c     Put the scalar vars
      tmp=ux(1:nxlocal,1:nylocal,1:nzlocal,1)
      dataname='ux'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      tmp=ux(1:nxlocal,1:nylocal,1:nzlocal,2)
      dataname='uy'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      tmp=ux(1:nxlocal,1:nylocal,1:nzlocal,3)
      dataname='uz'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      call close_hdf5_file(iunit)

      end
#endif
c-----------------------------------------------------------------------
#ifdef PARALLELHDF5
#ifdef TWO_D
      subroutine WriteParallelHDF5File(ux,index)
      use mesh
      use mpistuff
      use mesh_common
      use properties
      use iounits
       implicit none
c
      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI, nvar)
      double precision:: tmp(nxlsize, nzlsize)
c     Local variables
      double precision:: ttot,dt
      integer:: ndims,ierr1,dims(2),offset(2), count(2)
c     
      integer i,j,k
      integer index,iunit
      integer nbytes, itmp,  nelements
      character*50 outfile,filename,dataname
      integer:: dbfile
c      integer,parameter::ifld=1
c     

c     File name for data 
      write(outfile,1000) 
     &     float(index)/10000000.0
 1000 format('output',f8.7,'.h5')
c
      iunit=101
      filename=trim(outfile)//char(0)
c      write(6,*) 'Writing hdf5 file',filename
      call open_parhdf5_file(iunit,comm3D,filename)
c
      ndims=2
      dims(1)=nx
      dims(2)=nz
      count(1)=nxlocal
      count(2)=nzlocal
      offset(1)=(iprocx-1)*nxlocal
      offset(2)=(iprocz-1)*nzlocal
c
c
c     Put the scalar vars
      tmp=ux(1:nxlocal,1,1:nzlocal,1)
      dataname='ux'//char(0)
      call write_parhdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     count,offset,comm3d,trim(dataname))
c
      tmp=ux(1:nxlocal,1,1:nzlocal,2)
      dataname='uy'//char(0)
      call write_parhdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     count,offset,comm3d,trim(dataname))
c
      tmp=ux(1:nxlocal,1,1:nzlocal,3)
      dataname='uz'//char(0)
      call write_parhdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     count,offset,comm3d,trim(dataname))
c
      call MPI_BARRIER(comm3D,ierr)
      call close_hdf5_file(iunit)

      end
#else
c-----------------------------------------------------------------------
      subroutine WriteParallelHDF5File(ux,index)
      use mesh
      use mpistuff
      use mesh_common
      use properties
      use iounits
       implicit none
c
      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI, nvar)
      double precision:: tmp(nxlsize, nylsize,nzlsize)
c     Local variables
      double precision:: ttot,dt
      real:: xn(nxlocal+1),yn(nylocal+1),zn(nzlocal+1)
      integer:: ndims,ierr1,dims(3),offset(3), count(3)
c     
      integer i,j,k
      integer index,iunit
      integer nbytes, itmp,  nelements
      character*50 outfile,filename,dataname
      integer:: dbfile
c      integer,parameter::ifld=1
c     
c     File name for data 
      write(outfile,1000) 
     &     float(index)/10000000.0
 1000 format('output',f8.7,'.h5')
c
      iunit=101
      filename=trim(outfile)//char(0)
c      write(6,*) 'Writing hdf5 file',filename
      call open_parhdf5_file(iunit,comm3D,filename)
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
c      call ConservedToPrimitive(ux,vx)
c      call VelocityDivergence(ux,divu)
c      call Vorticity(ux,vort)
c
c     Put the scalar vars
c     Put the scalar vars
      do i=1,nxlocal,1
         tmp(i,:,:)=0.5D0*(ux(i,1:nylocal,1:nzlocal,1)+
     &        ux(i+1,1:nylocal,1:nzlocal,1))
      enddo
      dataname='ux'//char(0)
      call write_parhdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     count,offset,comm3d,trim(dataname))
c
      do j=1,nylocal,1
         tmp(:,j,:)=0.5D0*(ux(1:nxlocal,j,1:nzlocal,2)+
     &        ux(1:nxlocal,j+1,1:nzlocal,2))
      enddo
      dataname='uy'//char(0)
      call write_parhdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     count,offset,comm3d,trim(dataname))
c
      do k=1,nzlocal,1
         tmp(:,:,k)=0.5D0*(ux(1:nxlocal,1:nylocal,k,3)+
     &        ux(1:nxlocal,1:nylocal,k+1,3))
      enddo
      dataname='uz'//char(0)
      call write_parhdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     count,offset,comm3d,trim(dataname))
c
      call MPI_BARRIER(comm3D,ierr)
      call close_hdf5_file(iunit)

      end
#endif
#endif
c-----------------------------------------------------------------------
#ifdef SILO
#ifndef TWO_D
      subroutine WriteSiloFile(ux,index)
      use mesh
      use mesh_common
      use properties
      use iounits
      include "silo.inc"
c
      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI, nvar)
      double precision:: tmp(nxlsize, nylsize,nzlsize)
c      double precision:: divu(nxlsize, nylsize,nzlsize)
c      double precision:: vort(nxlsize, nylsize,nzlsize)
c      double precision:: Efld(nxlsize, nylsize,nzlsize)
c     Local variables
      double precision:: ttot,dt
      real:: xn(nxlocal+1),yn(nylocal+1), zn(nzlocal+1)
      integer:: ndims,ierr1,dims(3)
c     
      integer i,j,k
      integer index
      integer nbytes, offset,itmp, ierr, nelements
      character*50 outfile,fldfile, meshfile
      integer:: dbfile
c      integer,parameter::ifld=1
c     

c     File name for data 
      write(outfile,1000) 
     &     float(index)/10000000.0
c 1000 format('output',f7.6,'.silo')
 1000 format('output',f8.7,'.silo')
c
      write(6,*) 'Writing silo file',outfile
      ierr=dbcreate(outfile,19,DB_CLOBBER,DB_LOCAL,
     &     "Plasmoid data",22,DB_HDF5,dbfile)
      if(ierr.eq.-1) then
         write(6,*) 'Could not create silo file',outfile
         stop
      endif
c
c     Create nodal mesh
      do i=1,nxlocal+1,1
         xn(i)=xc(i)-0.5D0*dx         
      enddo
      do j=1,nylocal+1,1
         yn(j)=yc(j)-0.5D0*dy
      enddo
      do k=1,nzlocal+1,1
         zn(k)=zc(k)-0.5D0*dz
      enddo
      ndims=3
      dims(1)=nxlocal+1
      dims(2)=nylocal+1
      dims(3)=nzlocal+1
c
c     Write rectilinear mesh
      ierr1=dbputqm(dbfile,"quadmesh",8,"xc",2,"yc",2,
     &     "zc",2,xn,yn,zn,dims,ndims,
     &     DB_FLOAT,DB_COLLINEAR,DB_F77NULL,ierr)
c
c     Writing zone centered vars
      dims(1)=nxlocal
      dims(2)=nylocal
      dims(3)=nzlocal
c
c     Put the scalar vars
      do i=1,nxlocal,1
         tmp(i,:,:)=0.5D0*(ux(i,1:nylocal,1:nzlocal,1)+
     &        ux(i+1,1:nylocal,1:nzlocal,1))
      enddo
      ierr1=dbputqv1(dbfile,"ux",2,"quadmesh",8,
     &     tmp,dims,ndims,DB_F77NULL,0,DB_DOUBLE,
     &     DB_ZONECENT,DB_F77NULL,ierr)
c
      do j=1,nylocal,1
         tmp(:,j,:)=0.5D0*(ux(1:nxlocal,j,1:nzlocal,2)+
     &        ux(1:nxlocal,j+1,1:nzlocal,2))
      enddo
      ierr1=dbputqv1(dbfile,"uy",2,"quadmesh",8,
     &     tmp,dims,ndims,DB_F77NULL,0,DB_DOUBLE,
     &     DB_ZONECENT,DB_F77NULL,ierr)
c
      do k=1,nzlocal,1
         tmp(:,:,k)=0.5D0*(ux(1:nxlocal,1:nylocal,k,3)+
     &        ux(1:nxlocal,1:nylocal,k+1,3))
      enddo
      ierr1=dbputqv1(dbfile,"uz",2,"quadmesh",8,
     &     tmp,dims,ndims,DB_F77NULL,0,DB_DOUBLE,
     &     DB_ZONECENT,DB_F77NULL,ierr)
c
      ierr=dbclose(dbfile)

      end
c
c
c-----------------------------------------------------------------------
      subroutine WriteTmpSiloCC(ux,ncomp,index,ilo,ihi,jlo,jhi,
     &     klo,khi)
      use mesh
      use mesh_common
      use properties
      use iounits
      implicit none
      include "silo.inc"
c
c
      integer:: ncomp
c      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,ncomp)
      integer:: ilo,ihi,jlo,jhi,klo,khi
c      double precision::ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision::ux(ilo:ihi,jlo:jhi,klo:khi,ncomp)
c      double precision::ux(-2:nxlocal+3,-2:nylocal+3,-2:nzlocal+3,ncomp)
      double precision:: tmp(nxlocal,nylocal,nzlocal)
c     Local variables
      double precision:: ttot,dt
      real:: xn(nxlocal+1),yn(nylocal+1)
      real:: zn(nzlocal+1)
      integer:: ndims,ierr1,dims(3)
c     
      integer i,j,k
      integer index
      integer nbytes, offset,itmp, ierr, nelements
      character*50 outfile,fldfile, meshfile
      integer:: dbfile
c      integer,parameter::ifld=1
c     

c     File name for data 
      write(outfile,1000) 
     &     float(index)/10.0
 1000 format('tmp',f2.1,'.silo')
c
      write(6,*) 'Writing tmp silo file',outfile
      ierr=dbcreate(outfile,10,DB_CLOBBER,DB_LOCAL,
     &     "MHDAllSpeed data",16,DB_HDF5,dbfile)
      if(ierr.eq.-1) then
         write(6,*) 'Could not create silo file',outfile
         stop
      endif
c
c     Create nodal mesh
      do i=1,nxlocal+1,1
         xn(i)=xc(i)-0.5D0*dx         
      enddo
      do j=1,nylocal+1,1
         yn(j)=yc(j)-0.5D0*dy
      enddo
      do k=1,nzlocal+1,1
         zn(k)=zc(k)-0.5D0*dz
      enddo
      ndims=3
      dims(1)=nxlocal+1
      dims(2)=nylocal+1
      dims(3)=nzlocal+1
c
c     Write rectilinear mesh
      ierr1=dbputqm(dbfile,"quadmesh",8,"xc",2,"yc",2,
     &     "zc",2,xn,yn,zn,dims,ndims,
     &     DB_FLOAT,DB_COLLINEAR,DB_F77NULL,ierr)
c
c     Writing zone centered vars
      dims(1)=nxlocal
      dims(2)=nylocal
      dims(3)=nzlocal
c
c     Put the scalar vars
      tmp=ux(1:nxlocal,1:nylocal,1:nzlocal,1)
      ierr1=dbputqv1(dbfile,"comp1",5,"quadmesh",8,
     &     tmp,dims,ndims,DB_F77NULL,0,DB_DOUBLE,
     &     DB_ZONECENT,DB_F77NULL,ierr)
c
      if(ncomp.eq.2) then
      tmp=ux(1:nxlocal,1:nylocal,1:nzlocal,2)
      ierr1=dbputqv1(dbfile,"comp2",5,"quadmesh",8,
     &     tmp,dims,ndims,DB_F77NULL,0,DB_DOUBLE,
     &     DB_ZONECENT,DB_F77NULL,ierr)
      endif
c
      ierr=dbclose(dbfile)
      return
      end
c-----------------------------------------------------------------------
      subroutine WriteTmpSiloF(ux,index,indx,indy,indz,ilo,ihi,jlo,jhi,
     &     klo,khi)
      use mesh
      use mesh_common
      use properties
      use iounits
      implicit none
      include "silo.inc"
c
c
      integer:: ilo,ihi,jlo,jhi,klo,khi
      integer:: indx,indy,indz
c      double precision::ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision::ux(ilo:ihi,jlo:jhi,klo:khi)
      double precision:: tmp(nxlocal+indx,nylocal+indy,nzlocal+indz)
c     Local variables
      double precision:: ttot,dt
      real:: xn(nxlocal+indx),yn(nylocal+indy),zn(nzlocal+indz)
      integer:: ndims,ierr1,dims(3)
c     
      integer i,j,k
      integer index
      integer nbytes, offset,itmp, ierr, nelements
      character*50 outfile,fldfile, meshfile
      integer:: dbfile
c      integer,parameter::ifld=1
c     

c     File name for data 
      write(outfile,1000) 
     &     float(index)/10.0
 1000 format('tmpFX',f2.1,'.silo')
c
      write(6,*) 'Writing tmp silo file',outfile
      ierr=dbcreate(outfile,12,DB_CLOBBER,DB_LOCAL,
     &     "MHDAllSpeed data",16,DB_HDF5,dbfile)
      if(ierr.eq.-1) then
         write(6,*) 'Could not create silo file',outfile
         stop
      endif
c
c     Create nodal mesh
      do i=1,nxlocal+indx,1
         xn(i)=xc(i)-0.5D0*dx*max(0,indx)
      enddo
      do j=1,nylocal+indy,1
         yn(j)=yc(j)-0.5D0*dy*max(0,indy)
      enddo
      do k=1,nzlocal+indz,1
         zn(k)=zc(k)-0.5D0*dz*max(0,indz)
      enddo
      ndims=3
      dims(1)=nxlocal+indx
      dims(2)=nylocal+indy
      dims(3)=nzlocal+indz
c
c     Write rectilinear mesh
      ierr1=dbputqm(dbfile,"quadmesh",8,"xc",2,"yc",2,
     &     "zc",2,xn,yn,zn,dims,ndims,
     &     DB_FLOAT,DB_COLLINEAR,DB_F77NULL,ierr)
c
c     Put the scalar vars
      tmp=ux(1:nxlocal+indx,1:nylocal+indy,1:nzlocal+indz)
      ierr1=dbputqv1(dbfile,"comp1",5,"quadmesh",8,
     &     tmp,dims,ndims,DB_F77NULL,0,DB_DOUBLE,
     &     DB_NODECENT,DB_F77NULL,ierr)
c
      ierr=dbclose(dbfile)
      return
      end
c-----------------------------------------------------------------------
      subroutine WriteTmpSiloFY(ux,index)
      use mesh
      use mesh_common
      use properties
      use iounits
      implicit none
      include "silo.inc"
c
c
      double precision:: ux(IXLO:IXHI,IYLO:IYHI+1)
      double precision:: tmp(nxlocal,nylocal+1)
c     Local variables
      double precision:: ttot,dt
      real:: xn(nxlocal),yn(nylocal+1)
      integer:: ndims,ierr1,dims(2)
c     
      integer i,j,k
      integer index
      integer nbytes, offset,itmp, ierr, nelements
      character*50 outfile,fldfile, meshfile
      integer:: dbfile
c      integer,parameter::ifld=1
c     

c     File name for data 
      write(outfile,1000) 
     &     float(index)/10.0
 1000 format('tmpFY',f2.1,'.silo')
c
      write(6,*) 'Writing tmp silo file',outfile
      ierr=dbcreate(outfile,12,DB_CLOBBER,DB_LOCAL,
     &     "MHDAllSpeed data",16,DB_HDF5,dbfile)
      if(ierr.eq.-1) then
         write(6,*) 'Could not create silo file',outfile
         stop
      endif
c
c     Create nodal mesh
      do i=1,nxlocal,1
         xn(i)=xc(i)
      enddo
      do j=1,nylocal+1,1
         yn(j)=yc(j)-0.5D0*dy
      enddo
      ndims=2
      dims(1)=nxlocal
      dims(2)=nylocal+1
c
c     Write rectilinear mesh
      ierr1=dbputqm(dbfile,"quadmesh",8,"xc",2,"yc",2,
     &     "zc",2,xn,yn,DB_F77NULL,dims,ndims,
     &     DB_FLOAT,DB_COLLINEAR,DB_F77NULL,ierr)
c
c     Writing zone centered vars
      dims(1)=nxlocal
      dims(2)=nylocal+1
c
c     Put the scalar vars
      tmp=ux(1:nxlocal,1:nylocal+1)
      ierr1=dbputqv1(dbfile,"comp1",5,"quadmesh",8,
     &     tmp,dims,ndims,DB_F77NULL,0,DB_DOUBLE,
     &     DB_NODECENT,DB_F77NULL,ierr)
c
      ierr=dbclose(dbfile)
      return
      end
#else
c     two d
c-----------------------------------------------------------------------
      subroutine WriteSiloFile(ux,index)
      use mesh
      use mesh_common
      use properties
      use iounits
      include "silo.inc"
c
      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI, nvar)
      double precision:: tmp(nxlsize, nylsize,nzlsize)
c      double precision:: divu(nxlsize, nylsize,nzlsize)
c      double precision:: vort(nxlsize, nylsize,nzlsize)
c      double precision:: Efld(nxlsize, nylsize,nzlsize)
c     Local variables
      double precision:: ttot,dt
      real:: xn(nxlocal+1), zn(nzlocal+1)
      integer:: ndims,ierr1,dims(2)
c     
      integer i,j,k
      integer index
      integer nbytes, offset,itmp, ierr, nelements
      character*50 outfile,fldfile, meshfile
      integer:: dbfile
c      integer,parameter::ifld=1
c     

c     File name for data 
      write(outfile,1000) 
     &     float(index)/10000000.0
c 1000 format('output',f7.6,'.silo')
 1000 format('output',f8.7,'.silo')
c
      write(6,*) 'Writing silo file',outfile
      ierr=dbcreate(outfile,19,DB_CLOBBER,DB_LOCAL,
     &     "Plasmoid data",22,DB_HDF5,dbfile)
      if(ierr.eq.-1) then
         write(6,*) 'Could not create silo file',outfile
         stop
      endif
c
c     Create nodal mesh
      do i=1,nxlocal+1,1
         xn(i)=xc(i)-0.5D0*dx         
      enddo
      do k=1,nzlocal+1,1
         zn(k)=zc(k)-0.5D0*dz
      enddo
      ndims=2
      dims(1)=nxlocal+1
      dims(2)=nzlocal+1
c
c     Write rectilinear mesh
      ierr1=dbputqm(dbfile,"quadmesh",8,"xc",2,"yc",2,
     &     "zc",2,xn,zn,DB_F77NULL,dims,ndims,
     &     DB_FLOAT,DB_COLLINEAR,DB_F77NULL,ierr)
c
c     Writing zone centered vars
      dims(1)=nxlocal
      dims(2)=nzlocal
c
c     Put the scalar vars
      do i=1,nxlocal,1
         tmp(i,:,:)=0.5D0*(ux(i,1:nylocal,1:nzlocal,1)+
     &        ux(i+1,1:nylocal,1:nzlocal,1))
      enddo
      ierr1=dbputqv1(dbfile,"ux",2,"quadmesh",8,
     &     tmp,dims,ndims,DB_F77NULL,0,DB_DOUBLE,
     &     DB_ZONECENT,DB_F77NULL,ierr)
c
      do k=1,nzlocal,1
         tmp(:,:,k)=0.5D0*(ux(1:nxlocal,1:nylocal,k,3)+
     &        ux(1:nxlocal,1:nylocal,k+1,3))
      enddo
      ierr1=dbputqv1(dbfile,"uz",2,"quadmesh",8,
     &     tmp,dims,ndims,DB_F77NULL,0,DB_DOUBLE,
     &     DB_ZONECENT,DB_F77NULL,ierr)
c
      ierr=dbclose(dbfile)

      end
c
#endif
#endif
c-----------------------------------------------------------------------
