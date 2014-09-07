
c
      subroutine output_xslice2(ux)
      use time
      use mesh
      use mesh_common
      use virtualwallBC
      use stats
      use iounits
      use properties
#ifdef PARALLEL
      use mpistuff
#endif
c
      double precision:: ux(nxlocal+1,nylocal+1,nzlocal+1,nvar)
      double precision:: kappat, uept, hfact, utaut, Rethet
      integer:: idata, i,j,k
      double precision:: nu, ut, vt, wt, et, g_et
c
      nu = 1.d0/Reynolds

cccccccccccccccccccccccccccccccccccccc
       et = 0.d0
       if(iprocx.eq.1) et = ux(1,1,nzlocal-1,1)
       do i=1, nxlocal, 1
         do j=1, nylocal, 1
           do k=1, nzlocal, 1
           ut = 0.5d0*(ux(i,j,k,1)+ux(i+1,j,k,1))
           vt = 0.5d0*(ux(i,j,k,2)+ux(i,j+1,k,2))
           wt = 0.5d0*(ux(i,j,k,3)+ux(i,j,k+1,3))
c           et =et +  ut+vt+wt !0.5d0*(ut*ut+vt*vt+wt*wt)*dx*dy*dz
c            write(*,*) i,j,k,et
           enddo
         enddo
       enddo
#ifdef PARALLEL
       g_et = 0.d0
       call MPI_Allreduce(et, g_et, 1,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
              Call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
       g_et = et
#endif
       write(itestv,*)'adu', ttot,g_et
c
c
      return
      end


      subroutine output_bcuwall
      use time
      use mesh
      use mesh_common
      use virtualwallBC
      use stats
      use iounits
      use properties
#ifdef PARALLEL
      use mpistuff
#endif
c
      implicit none
      double precision:: uxau(nxlocal+1,nylocal)
      double precision:: kappat, uept, hfact, utaut, Rethet
      integer:: idata, i,j,k
      double precision:: nu, ut, vt, wt, et, g_et
c
      nu = 1.d0/Reynolds

        uxau(1:nxlocal+1,1:nylocal) = bc_utau(1:nxlocal+1,1:nylocal)
cccccccccccccccccccccccccccccccccccccc
       et = 0.d0
       if(iprocz.eq.1) then
       do i=1, nxlocal, 1
       do j=1, nylocal, 1
c         j=1
           ut = 0.5d0*(uxau(i,j)+uxau(i+1,j))
c           ut =  ux(i,j)
           et =et + ut !0.5d0*(ut*ut+vt*vt+wt*wt)*dx*dy*dz
c            write(*,*) i,j,k,et
         enddo
       enddo
       endif
       
#ifdef PARALLEL
       g_et = 0.d0
       call MPI_Allreduce(et, g_et, 1,
     &          MPI_DOUBLE_PRECISION, MPI_SUM, comm3D, ierr)
              Call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
       g_et = et
#endif
       write(itestv,*)'bcutau', ttot,g_et
c
c
      return
      end

#endif

c----------------------------------------------------------------------------------
c      some testing output: in HDF5 format
c---------------------------------------------------------------------------------
c
c--------------------------------------------------------------
c     gradient of pressure
c
#ifdef HDF5
      subroutine WriteHDF5Filep(gradp,index)
      use mesh
      use mesh_common
      use properties
      use iounits
c
      double precision::
     &    gradp(1:nxlocal+1,1:nylocal,1:nzlocal+1, nvar)
      double precision:: tmp(nxlsize, nylsize,nzlsize)
c      double precision:: vort(nxlsize, nylsize,nzlsize)
c     Local variables
c      double precision:: ttot ,dt
c      real:: xn(nxlocal+1),yn(nylocal+1),zn(nzlocal+1)
      integer:: ndims,dims(3)  ! ,ierr1
c     
c      integer i,j,k 
      integer index,iunit
c      integer nbytes, offset, nelements ,itmp, ierr
      character*50 outfile,filename,dataname
c      integer:: dbfile
c      integer,parameter::ifld=1
c     

c     File name for data 
      write(outfile,1000) 
     &     float(index)/10000000.0
 1000 format('optput',f8.7,'.h5')
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
      tmp=gradp(1:nxlocal,1:nylocal,1:nzlocal,1)
      dataname='ux'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      tmp=gradp(1:nxlocal,1:nylocal,1:nzlocal,2)
      dataname='uy'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      tmp=gradp(1:nxlocal,1:nylocal,1:nzlocal,3)
      dataname='uz'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      call close_hdf5_file(iunit)

      end
#endif
        
c-----------------------------------------------------------------
c     divergence of velocity 
c
#ifdef HDF5
      subroutine WriteHDF5Filedivu(divu,index)
      use mesh
      use mesh_common
      use properties
      use iounits
c
c      double precision::
c     &    divu(1:nxlocal,1:nylocal,1:nzlocal)
      double precision:: tmp(nxlsize, nylsize,nzlsize)
      double precision:: divu(nxlsize, nylsize,nzlsize)
c      double precision:: vort(nxlsize, nylsize,nzlsize)
c     Local variables
c      double precision:: dt,ttot
c      real:: xn(nxlocal+1),yn(nylocal+1),zn(nzlocal+1)
      integer:: ndims,dims(3)  ! ,ierr1
c     
c      integer i,j,k
      integer index,iunit
c      integer nbytes, offset,itmp, ierr, nelements
      character*50 outfile,filename,dataname
c      integer:: dbfile
c      integer,parameter::ifld=1
c     

c     File name for data 
      write(outfile,1000) 
     &     float(index)/10000000.0
 1000 format('divput',f8.7,'.h5')
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
      tmp=divu(1:nxlocal,1:nylocal,1:nzlocal)
      dataname='ux'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      tmp=divu(1:nxlocal,1:nylocal,1:nzlocal)
      dataname='uy'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      tmp=divu(1:nxlocal,1:nylocal,1:nzlocal)
      dataname='uz'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      call close_hdf5_file(iunit)

      end
#endif
        
c--------------------------------------------------------------
c     SGS tensor at CorNer
c
#ifdef HDF5
      subroutine WriteHDF5FileTijcn(Tijcn,index)
      use mesh
      use mesh_common
      use properties
      use iounits
c
c      double precision::
c     &    divu(1:nxlocal,1:nylocal,1:nzlocal)
      double precision:: tmp(nxlsize+3, nylsize+3,nzlsize+3)
      double precision:: Tijcn(0:nxlocal+2, 0:nylocal+2,0:nzlocal+2,6)
c      double precision:: vort(nxlsize, nylsize,nzlsize)
c     Local variables
c      double precision:: dt,ttot
c      real:: xn(nxlocal+1),yn(nylocal+1),zn(nzlocal+1)
      integer:: ndims,dims(3)  ! ,ierr1
c     
c      integer i,j,k
      integer index,iunit
c      integer nbytes, offset,itmp, ierr, nelements
      character*50 outfile,filename,dataname
c      integer:: dbfile
c      integer,parameter::ifld=1
c     

c     File name for data 
      write(outfile,1000) 
     &     float(index)/10000000.0
 1000 format('Tijcn',f8.7,'.h5')
c
      iunit=101
      filename=trim(outfile)//char(0)
c      write(6,*) 'Writing hdf5 file',filename
      call open_hdf5_file(iunit,filename)
c
      ndims=3
      dims(1)=nxlocal+3
      dims(2)=nylocal+3
      dims(3)=nzlocal+3
c
c     Put the scalar vars
      tmp=Tijcn(0:nxlocal+2, 0:nylocal+2,0:nzlocal+2,1)
      dataname='ux'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      tmp=Tijcn(0:nxlocal+2, 0:nylocal+2,0:nzlocal+2,2)
      dataname='uy'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      tmp=Tijcn(0:nxlocal+2, 0:nylocal+2,0:nzlocal+2,3)
      dataname='uz'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      call close_hdf5_file(iunit)

      end
#endif

        
c--------------------------------------------------------------
c     velocity at CorNer
c
#ifdef HDF5
      subroutine WriteHDF5FileUcn(ucn,index)
      use mesh
      use mesh_common
      use properties
      use iounits
c
c      double precision::
c     &    divu(1:nxlocal,1:nylocal,1:nzlocal)
      double precision:: tmp(nxlsize+5, nylsize+5,nzlsize+5)
      double precision:: Ucn(-1:nxlocal+3, -1:nylocal+3,-1:nzlocal+3,3)
c      double precision:: vort(nxlsize, nylsize,nzlsize)
c     Local variables
c      double precision:: dt,ttot
c      real:: xn(nxlocal+1),yn(nylocal+1),zn(nzlocal+1)
      integer:: ndims,dims(3)  ! ,ierr1
c     
c      integer i,j,k
      integer index,iunit
c      integer nbytes, offset,itmp, ierr, nelements
      character*50 outfile,filename,dataname
c      integer:: dbfile
c      integer,parameter::ifld=1
c     

c     File name for data 
      write(outfile,1000) 
     &     float(index)/10000000.0
 1000 format('Ucn',f8.7,'.h5')
c
      iunit=101
      filename=trim(outfile)//char(0)
c      write(6,*) 'Writing hdf5 file',filename
      call open_hdf5_file(iunit,filename)
c
      ndims=3
      dims(1)=nxlocal+5
      dims(2)=nylocal+5
      dims(3)=nzlocal+5
c
c     Put the scalar vars
      tmp=Ucn(-1:nxlocal+3, -1:nylocal+3,-1:nzlocal+3,1)
      dataname='ux'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      tmp=Ucn(-1:nxlocal+3, -1:nylocal+3,-1:nzlocal+3,2)
      dataname='uy'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      tmp=Ucn(-1:nxlocal+3, -1:nylocal+3,-1:nzlocal+3,3)
      dataname='uz'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      call close_hdf5_file(iunit)

      end
#endif

c-----------------------------------------------------------------
c     Kinetic energy in fourier phase space
c
#ifdef HDF5
      subroutine WriteHDF5Fileef(ef3,nxi,index)
      use mesh
      use mesh_common
      use properties
      use iounits
c
c      double precision::
c     &    divu(1:nxlocal,1:nylocal,1:nzlocal)
      double precision:: tmp(nxi, nylsize,nzlsize)
      double precision:: ef3(nxi, nylsize,nzlsize)
c      double precision:: vort(nxlsize, nylsize,nzlsize)
c     Local variables
c      double precision:: dt,ttot
c      real:: xn(nxlocal+1),yn(nylocal+1),zn(nzlocal+1)
      integer:: ndims,dims(3)  ! ,ierr1
c     
c      integer i,j,k
      integer index,iunit, nxi
c      integer nbytes, offset,itmp, ierr, nelements
      character*50 outfile,filename,dataname
c      integer:: dbfile
c      integer,parameter::ifld=1
c     

c     File name for data 
      write(outfile,1001) 
     &     float(index)/10000000.0
 1001 format('eftput',f8.7,'.h5')
c
      iunit=101
      filename=trim(outfile)//char(0)
c      write(6,*) 'Writing hdf5 file',filename



      call open_hdf5_file(iunit,filename)
c
      ndims=3
      dims(1)=nxi
      dims(2)=nylocal
      dims(3)=nzlocal
c
c     Put the scalar vars
      tmp=ef3(1:nxi,1:nylocal,1:nzlocal)
      dataname='ux'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      tmp=ef3(1:nxi,1:nylocal,1:nzlocal)
      dataname='uy'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      tmp=ef3(1:nxi,1:nylocal,1:nzlocal)
      dataname='uz'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      call close_hdf5_file(iunit)

      end
#endif

c-----------------------------------------------------------------------
c----------------------------------------------------------------
#ifdef HDF5
      subroutine WriteHDF5Fileuwall(uwall,indax,index)
      use mesh
      use mesh_common
      use properties
      use iounits
c
      double precision:: uwall(IXLO:IXHI,IYLO:IYHI, nvar)
      double precision:: uw3(IXLO:IXHI,IYLO:IYHI,1, nvar)
      double precision:: tmp(IXLO:IXHI,IYLO:IYHI,1)
c     Local variables
      double precision:: ttot,dt
      real:: xn(nxlocal+1),yn(nylocal+1),zn(nzlocal+1)
      integer:: ndims,ierr1,dims(3)
c     
      integer i,j,k
      integer index,iunit, indax
      integer nbytes, offset,itmp, ierr, nelements
      character*50 outfile,filename,dataname
      integer:: dbfile
c      integer,parameter::ifld=1
c     
      uw3(:,:,1,1)=uwall(:,:,1)
      uw3(:,:,1,2)=uwall(:,:,2)
      uw3(:,:,1,3)=uwall(:,:,3)

c     File name for data 
      write(outfile,1000)float(indax)/100000.0, 
     &     index
 1000 format('bcuwal',f7.6,i1,'.h5')
c
      iunit=101
      filename=trim(outfile)//char(0)
c      write(6,*) 'Writing hdf5 file',filename
      call open_hdf5_file(iunit,filename)
c
      ndims=3
      dims(1)=nxlocal+8
      dims(2)=nylocal+8
      dims(3)=1
c
c     Put the scalar vars
      tmp=uw3(ixlo:ixhi,iylo:iyhi,1:1,1)
      dataname='ux'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      tmp=uw3(ixlo:ixhi,iylo:iyhi,1:1,2)
      dataname='uy'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      tmp=uw3(ixlo:ixhi,iylo:iyhi,1:1,3)
      dataname='uz'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      call close_hdf5_file(iunit)

      end
#endif

c-----------------------------------------------------------------------
c----------------------------------------------------------------
#ifdef HDF5
      subroutine WriteHDF5Fileutau(utau,indax,index)
      use mesh
      use mesh_common
      use properties
      use iounits
c
      double precision:: utau(1:nxlocal+1,1:nylocal, nvar)
      double precision:: ut3(1:nxlocal+1,1:nylocal,1, nvar)
      double precision:: tmp(1:nxlocal+1,1:nylocal,1)
c     Local variables
      double precision:: ttot,dt
      real:: xn(nxlocal+1),yn(nylocal+1),zn(nzlocal+1)
      integer:: ndims,ierr1,dims(3)
c     
      integer i,j,k
      integer index,iunit, indax
      integer nbytes, offset,itmp, ierr, nelements
      character*50 outfile,filename,dataname
      integer:: dbfile
c      integer,parameter::ifld=1
c     
      ut3(:,:,1,1)=utau(:,:,1)
      ut3(:,:,1,2)=utau(:,:,2)
      ut3(:,:,1,3)=utau(:,:,3)

c     File name for data 
      write(outfile,1000)float(indax)/100000.0, 
     &     index
 1000 format('bcutau',f7.6,i1,'.h5')
c
      iunit=101
      filename=trim(outfile)//char(0)
c      write(6,*) 'Writing hdf5 file',filename
      call open_hdf5_file(iunit,filename)
c
      ndims=3
      dims(1)=nxlocal+1
      dims(2)=nylocal
      dims(3)=1
c
c     Put the scalar vars
      tmp=ut3(1:nxlocal+1,1:nylocal,1:1,1)
      dataname='ux'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      tmp=ut3(1:nxlocal+1,1:nylocal,1:1,2)
      dataname='uy'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      tmp=ut3(1:nxlocal+1,1:nylocal,1:1,3)
      dataname='uz'//char(0)
      call write_hdf5_f(iunit,tmp,'DOUBLE'//char(0),ndims,dims,
     &     trim(dataname))
c
      call close_hdf5_file(iunit)

      end
#endif












cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccc
       subroutine out_inlet(tstep)
       use mesh
       use properties
       use InflowVelBC
c 
       implicit none
c       double precision
       character*50 outfile,filename
       integer tstep , k

      write(outfile,1000) 
     &     float(tstep)/10000000.0
 1000 format('inlet',f8.7,'.curve')
      filename=trim(outfile)//char(0)
 
       open(12,file=filename)
       write(12,*) '#',filename
       do k =1 ,nzlocal, 1
         write(12,*) k, Uin(1,nylocal/2,k,1),Uin(1,nylocal/2,k,2),
     &  Uin(1,nylocal/2,k,3)
       enddo
       close(12)
       return
       end


