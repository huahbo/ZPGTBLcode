c----------------------------------------------------------------------
c     Wan Cheng 
c     August 2012
c----------------------------------------------------------------------
c
c
c
c----------------------------------------------------------------------
c     modules diagnostics
c----------------------------------------------------------------------
      module stats
      save
c  some parameters
      double precision:: Tave, bcupper
      integer:: nsteppre, nstats, numstat
      integer:: inflowstat, ioutput
c  local statistics
      double precision,allocatable:: st_umean(:,:,:)
      double precision,allocatable:: st_umeac(:,:,:)
      double precision,allocatable:: st_Tijct(:,:,:,:)
      double precision,allocatable:: st_delta(:)
      double precision,allocatable:: st_del99(:)
      double precision,allocatable:: st_theta(:)
      double precision,allocatable:: st_utau(:)
      double precision,allocatable:: g_fu2(:,:)  
      double precision,allocatable:: g_fu3(:,:)
      double precision,allocatable:: g_fu4(:,:)
      double precision,allocatable:: g_fuv(:,:)
      double precision,allocatable:: g_fuw(:,:)
      double precision,allocatable:: g_fvv(:,:)
      double precision,allocatable:: g_fvw(:,:)
      double precision,allocatable:: g_fww(:,:)                 
      double precision,allocatable:: g_Ten(:,:,:)
c  global statistics : only 1D allowed
      double precision,allocatable:: g_utau(:)
      double precision,allocatable:: g_delta(:)
      double precision,allocatable:: g_del99(:)
      double precision,allocatable:: g_theta(:)
      double precision,allocatable:: g_vwall(:)
      double precision,allocatable:: g_kappa(:)
      double precision,allocatable:: g_GCco(:)
      double precision,allocatable:: g_GCci(:)
      double precision,allocatable:: g_GCgam(:)
      double precision,allocatable:: g_GCpre(:) 
      end module stats
c
c-------------------------------------------
c      allocate stats
c-------------------------------------------
      subroutine allocate_stat
      use mesh
      use mesh_common
      use properties
      use stats
c
      implicit none
      ioutput = nx/2 - 5
      Tave = 10.d0
      bcupper = 2.5d0
      allocate(st_umean(1:nxlocal+1,izlo:izhi,3))
      allocate(st_umeac(1:nxlocal,1:nzlocal,3))
      allocate(st_Tijct(1:nxlocal,1:nylocal,1:nzlocal,6))
      allocate(st_delta(1:nxlocal+1))
      allocate(st_del99(1:nxlocal+1))
      allocate(st_theta(1:nxlocal+1))
      allocate(st_utau(1:nxlocal+1))
      allocate(g_fu2(1:nxlocal,1:nz))
      allocate(g_fu3(1:nxlocal,1:nz))
      allocate(g_fu4(1:nxlocal,1:nz))
      allocate(g_fuv(1:nxlocal,1:nz))
      allocate(g_fuw(1:nxlocal,1:nz))
      allocate(g_fvv(1:nxlocal,1:nz))
      allocate(g_fvw(1:nxlocal,1:nz))
      allocate(g_fww(1:nxlocal,1:nz))
      allocate(g_Ten(1:nxlocal,1:nz,6))
c
      allocate(g_utau(1:nx))
      allocate(g_delta(1:nx))
      allocate(g_del99(1:nx))
      allocate(g_theta(1:nx))
      allocate(g_vwall(1:nx))
c
      allocate(g_kappa(1:nx))
      allocate( g_GCco(1:nx))
      allocate( g_GCci(1:nx))
      allocate(g_GCgam(1:nx))
      allocate(g_GCpre(1:nx)) 
c

      st_umean = 0.d0
      st_delta = 0.d0
      st_del99 = 0.d0
      st_theta = 0.d0
      st_utau = 1.d7
c
      g_utau  = 0.d0
      g_delta = 0.d0
      g_del99 = 0.d0
      g_theta = 0.d0
      g_vwall = 0.d0
c
      g_kappa = 0.d0
      g_GCco  = 0.d0
      g_GCci  = 0.d0
      g_GCgam = 0.d0
      g_GCpre = 0.d0
c
      return
      end       
c
c----------------------------------------------------------------------
c     io  for diags
c----------------------------------------------------------------------
      module diags_io
c
c     eigenvalues of T
      integer, parameter:: itijmax=191
      integer, parameter:: itijmin=192   
c     x slice: y-z plane
      integer, parameter:: itestv = 1234
      integer, parameter:: itestv1 = 1240
      integer, parameter:: itestv2 = 1235
      integer, parameter:: itestv3 = 1236
      integer, parameter:: itestv4 = 1237
      integer, parameter:: itestv5 = 1238
      integer, parameter:: itestv6 = 1239
c     velocity in x slice
      integer, parameter:: iout01   = 123401
      integer, parameter:: iout02   = 123402
      integer, parameter:: iout03   = 123403
      integer, parameter:: iout04   = 123404
      integer, parameter:: iout05   = 123405
      integer, parameter:: iout06   = 123406
      integer, parameter:: iout07   = 123407
      integer, parameter:: iout08   = 123408
      integer, parameter:: iout09   = 123409   
c     spectral in y direction
      integer, parameter:: ispy01   = 123501
c
      integer, parameter:: itestout01 = 2345
      integer, parameter:: itestout02 = 2346
c        
      end module diags_io
c
c----------------------------------
c     init diags_io
c----------------------------------
      subroutine initdiagsio(new)
      use mesh_common
      use diags_io
#ifdef PARALLEL
      use mpistuff
#endif
      use OptionFlags
c   
      implicit none
      integer:: new
c      character*25 timestamp      
c      call fdate(timestamp)
c
      if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
        if(new.eq.-1) then
c
          open(itijmax,file='EigvalTij_max.txt')
          write(itijmax,*) '#MaxEigenvalueofTij'
c
          open(itijmin,file='Eigvaltij_min.txt')
          write(itijmin,*) '#MinEigenvalueofTij'         
c
          open(itestv,file='1kineticenergy.dat')
          write(itestv,*) '#kinetic energy'
c
          open(itestv1,file='7xretheta.dat')
          write(itestv1,*) '#xRethetarelations'
c
          open(itestv2,file='2utauin.dat')
          write(itestv2,*) '#utauatinlet'
c
          open(itestv3,file='3Retheta.dat')
          write(itestv3,*) '#reynoldsnumberinlet'
c
          open(itestv4,file='4Ueplus.dat')
          write(itestv4,*) '#ueplusatslice'
c
          open(itestv5,file='5kappa.dat')
          write(itestv5,*) '#Karmanconstantatslice'
c
          open(itestv6,file='6Hfactor.dat')
          write(itestv6,*) '#Hfactoratslice'
c------------------------------------------------
c
          open(iout01,  file='01meanu.dat')
          write(iout01,  *) '#meanu'
c------------------------------------------------
          open(ispy01,  file='11specty.dat')
          write(ispy01,  *) '#spectralinydirection'
c
          open(itestout01,  file='80spany.dat')
          write(itestout01,  *) '#spany'
          open(itestout02,  file='81umeanpoint.dat')
          write(itestout02,  *) '#umeanpoint'
c       
        else
          open(itestv,file='1kineticenergy.dat',position='Append')
          open(itestv2,file='2utauin.dat',      position='Append')
          open(itestv3,file='3Retheta.dat',     position='Append')
          open(itestv4,file='4Ueplus.dat',      position='Append')
          open(itestv5,file='5kappa.dat',       position='Append')
          open(itestv6,file='6Hfactor.dat',     position='Append')
          open(itestv1,file='7xretheta.dat',    position='Append')
          open(iout01, file='01meanu.dat',      position='Append')
        endif
      endif
c
      return
      end
c
c----------------------------------------------------------------------
c     module for recycling
c----------------------------------------------------------------------
      module recycle
      save
      integer:: g_ire, ip_re, i_re
      double precision:: nscale
      end module
c----------------------------------
c     init module recycle
c----------------------------------
      subroutine init_re
      use mesh
      use mesh_common
      use recycle

      g_ire = int(nx*0.8d0)
      ip_re = int(g_ire/nxlocal)+1
      i_re = mod(g_ire,nxlocal)
c
      nscale = 1.d0/12.d0 ! for 1/[2(n-1)] with n=7 
c
      return
      end
c
c----------------------------------------------------------------------





























