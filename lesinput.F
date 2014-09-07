c----------------------------------------------------------------------
c     Wan Cheng 
c     August 2012
c----------------------------------------------------------------------
c
c
c----------------------------------------------------------------------
c     init flow codeB
c----------------------------------------------------------------------
#ifdef CODEB
      subroutine readinflowB(timestep)
      use mesh
      use mesh_common
      use properties
      use InflowVelBC
      use virtualwallBC
c
      implicit none
      integer FileNumber
      character*50 DumpFile
      integer:: i,j,k,l, ilo, timestep
c

        FileNumber = 1*10000 + iprocy*100 + iprocz
        write(DumpFile, 1002) float(FileNumber) / 1000000,
     &                        float(timestep+1000)/1000000 !!!! 
c            1000 delete unphysical part
 1002   format('outcodeA1', f7.6,f7.6)
c       
        open(85,file=DumpFile)
        read(85,*) 
        do k=izlo, izhi, 1
          do j=iylo, iyhi, 1
            read(85,*)Uin(1,j,k,1),Uin(1,j,k,2),Uin(1,j,k,3)
            do i=2, nghost+1
            read(85,*)Uin(nghost+3-i,j,k,1),Uin(nghost+3-i,j,k,2),
     &                Uin(nghost+3-i,j,k,3)
            enddo
          enddo
        enddo     
        do j= iylo, iyhi, 1
            read(85,*) bc_uwall(1,j,1),bc_uwall(1,j,2),bc_uwall(1,j,3)
          do i=2, nghost+1
c            read(85,*) bc_uwall(2-i,j,1), bc_uwall(2-i,j,2),
c     &                                    bc_uwall(2-i,j,3)
            read(85,*) bc_uwall(-nghost-1+i,j,1),
     &        bc_uwall(-nghost-1+i,j,2),bc_uwall(-nghost-1+i,j,3)
          enddo
        enddo
        do j=1, nylocal, 1
          read(85,*) bc_utau(1,j)
          bc_dudz(1,j) = bc_utau(1,j)*bc_utau(1,j)*Reynolds
        enddo
        close(85)
c
        bc_uwall(1:ixhi,iylo:iyhi,2:3) = 0.d0  ! because we don't need this.
c
        return
        end
ccc--------------------------------------------------------------------
      subroutine initcodeB(ux)
      use mesh
      use mesh_common
      use properties
      use InflowVelBC
      use virtualwallBC
c
      implicit none
      double precision:: ux(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      integer:: i, j, k, ilo
c     init flow----------------------------------
      do i= ixlo, ixhi, 1
        do j= iylo, iyhi, 1
          do k= izlo, izhi, 1
            ux(i,j,k,1:3)=Uin(1,j,k,1:3)
          enddo
        enddo
      enddo
      ilo = ixlo
      if(iprocx.eq.1) ilo = 2
      do j=iylo, iyhi, 1
        bc_uwall(ilo:ixhi,j,1)=bc_uwall(1,j,1)
        bc_uwall(ilo:ixhi,j,2:3) = 0.d0 
      enddo
      do j=1, nylocal
        bc_utau(2:nxlocal+1,j)=bc_utau(1,j)
        bc_dudz(:,j) = bc_utau(:,j)*bc_utau(:,j)*Reynolds
      enddo
      
      call SetBoundaryValuesST(ux)
      call exchange(ux,3)
c------------------------------------------------
      return
      end
#endif
c
c----------------------------------------------------------------------
c     Adding ramdom perturbation to the flowfiled
c----------------------------------------------------------------------
      subroutine addperturbation(ux)
      use mesh
      use mesh_parms
      use mesh_common
c    
      implicit none
      double precision:: ux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: ut(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: pert, pai
      double precision:: ran(1000000)
      integer:: i,j,k,l, jin
c 
      call random_seed()
      call random_number(ran)
      pai = 4.d0 * datan(1.d0)
      pert = 0.05d0
      l = 1
      do i=ixlo,ixhi, 1
        do j=iylo, iyhi, 1
          do k= izlo, izhi, 1
            if( ux(i,j,k,1).lt.0.95d0 ) then
c            ux(i,j,k,1)=ux(i,j,k,1)*(1.d0 + pert* (ran(l)-0.5d0))
c            ux(i,nylocal+1-j,k,1) = ux(i,nylocal+1-j,k,1)
c     &                *(1.d0 - pert* (ran(l)-0.5d0))
c            ux(i,j,k,3)=ux(i,j,k,3)*(1.d0 + pert* (ran(l+1)-0.5d0))
           ux(i,j,k,1)=ux(i,j,k,1)*(1.d0+pert*dcos(2.d0*pai/3.d0*yc(j)))
            l = l+2
            endif
          enddo
        enddo
      enddo
c
      call SetBoundaryValuesST(ux)
      call exchange(ux,3)
c   
      ut = ux
      do l=1, 50, 1
        call SetBoundaryValuesST(ut)
        call exchange(ut,3)
        do i=1,nxlocal+1, 1
          do j=1,nylocal,1 
            do k=1, nzlocal, 1
             ux(i,j,k,1)=(2.d0*ut(i,j,k,1)+ut(i,j-1,k,1) 
     &                      +ut(i,j+1,k,1) )/4.d0
            enddo 
          enddo
        enddo
        ut = ux
      enddo
c       
      call SetBoundaryValuesST(ux)
      call exchange(ux,3)
c
      return
      end
c
c----------------------------------------------------------------------

