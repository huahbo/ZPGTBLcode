C -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2011
c     King Abdullah University of Science and Technology
c     All Rights Reserved
c-----------------------------------------------------------------
c     $Log: utils.F,v $
c     Revision 1.1  2011/07/25 11:32:01  samtanr
c     Original source.
c
c     Wind code: Incompressible Navier-Stokes LES
c
c-----------------------------------------------------------------
c     For FFT3D
c     Input uc
c     Output uf
c     idir=1 ==> forward
c     idir=-1==> backward fft
c-----------------------------------------------------------------------
      subroutine FFT3D(uc,uf,imax,jmax,kmax,idir)
      implicit none
      integer:: imax,jmax,kmax, idir
      double complex:: uc(imax,jmax,kmax)
      double complex:: uf(imax,jmax,kmax)
      double complex:: uctmp(imax)
      double complex:: uftmp(imax)
      double complex:: uctmp2(jmax)
      double complex:: uftmp2(jmax)
      double complex:: uctmp3(kmax)
      double complex:: uftmp3(kmax)
      integer:: i,j,k
c
c     3D FFT is basically calling 1D FFT in all directions
      do k=1,kmax,1
      do j=1,jmax,1
         uctmp(:)=uc(:,j,k)
         call FFT(uctmp,uftmp,imax,idir)
         uf(:,j,k)=uftmp(:)
      enddo
      enddo
c
      do k=1,kmax,1
      do i=1,imax,1
         uctmp2(:)=uf(i,:,k)
         call FFT(uctmp2,uftmp2,jmax,idir)
         uf(i,:,k)=uftmp2(:)
      enddo
      enddo
c
      do j=1,jmax,1
      do i=1,imax,1
         uctmp3(:)=uf(i,j,:)
         call FFT(uctmp3,uftmp3,kmax,idir)
         uf(i,j,:)=uftmp3(:)
      enddo
      enddo
c
      if(idir.eq.-1) uf=uf/imax/jmax/kmax
      return
      end

c--------------------------------------------------------
c     For FFT2D
c     Input uc
c     Output uf
c     idir=1 ==> forward
c     idir=-1==> backward fft
c-----------------------------------------------------------------------
      subroutine FFT2D(uc,uf,imax,jmax,idir)
      implicit none
      integer:: imax,jmax, idir
      double complex:: uc(imax,jmax)
      double complex:: uf(imax,jmax)
      double complex:: uctmp(imax)
      double complex:: uftmp(imax)
      double complex:: uctmp2(jmax)
      double complex:: uftmp2(jmax)
      integer:: i,j
c
c     2D FFT is basically calling 1D FFT in both directions
      do j=1,jmax,1
         uctmp(:)=uc(:,j)
         call FFT(uctmp,uftmp,imax,idir)
         uf(:,j)=uftmp(:)
      enddo
c
      do i=1,imax,1
         uctmp2(:)=uf(i,:)
         call FFT(uctmp2,uftmp2,jmax,idir)
         uf(i,:)=uftmp2(:)
      enddo
c
      if(idir.eq.-1) uf=uf/imax/jmax
      return
      end
c--------------------------------------------------------
c     For FFT1D
c     Input uc
c     Output uf
c     idir=1 ==> forward
c     idir=-1==> backward fft
c--------------------------------------------------------
c 
      subroutine FFT1D(uc,uf,jmax,idir)
      implicit none
      integer:: jmax, idir
      double complex:: uc(jmax)
      double complex:: uf(jmax)
      double complex:: uctmp(jmax)
      double complex:: uftmp(jmax)
c      integer:: i,j,k
c
         uctmp(:)=uc(:)
         call FFT(uctmp,uftmp,jmax,idir)
         uf(:)=uftmp(:)
c
      if(idir.eq.-1) uf=uf/jmax
      return
      end

c-----------------------------------------------------------------------
c     For FFT
c     Input uc
c     Output uf
c     idir=1 ==> forward
c     idir=-1==> backward fft
c-----------------------------------------------------------------------
      subroutine FFT(uc,uf,np,idir)
      integer:: np, idir
      double complex:: uc(np)
      double complex:: uf(np)
      integer*8 plan
      include 'fftw3.f'
c
      if(idir.eq.1) then
         call dfftw_plan_dft_1d(plan,np,uc,uf,
     &        FFTW_FORWARD,FFTW_ESTIMATE)
c         call dfftw_execute(plan)
         call dfftw_execute_dft(plan,uc,uf)
         call dfftw_destroy_plan(plan)
      else
         call dfftw_plan_dft_1d(plan,np,uc,uf,
     &        FFTW_BACKWARD,FFTW_ESTIMATE)
c         call dfftw_execute(plan)
         call dfftw_execute_dft(plan,uc,uf)
         call dfftw_destroy_plan(plan)
      endif
      return
      end

c--------------------------------------------------------
c     Does fft in j direction only
c      subroutine FFT1Dold(uc,uf,imax,jmax,kmax,idir)
c      implicit none
c      integer:: imax,jmax,kmax, idir
c      double complex:: uc(imax,jmax,kmax)
c      double complex:: uf(imax,jmax,kmax)
c      double complex:: uctmp(jmax)
c      double complex:: uftmp(jmax)
c      integer:: i,k
c
c      do k=1,kmax,1
c      do i=1,imax,1
c         uctmp(:)=uf(i,:,k)
c         call FFT(uctmp,uftmp,jmax,idir)
c         uf(i,:,k)=uftmp(:)
c      enddo
c      enddo
c
c      if(idir.eq.-1) uf=uf/jmax
c      return
c      end
c
