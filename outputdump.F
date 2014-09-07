c
c----------------------------------------------------------------------
c     output of codeA result per 100 timesteps
c          including 1 and 2, in case of writing failure
c----------------------------------------------------------------------
#ifdef CODEA
      subroutine outputcodeA(ux,timestep)
      use mesh_parms
      use mesh_common
      use virtualwallBC
      use InflowVelBC
      use time
c
      implicit none
      double precision:: ux(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      double precision:: uxdata(1:nghost+1,iylo:iyhi,izlo:izhi,nvar)
      double precision:: uwalld(1:nghost+1,iylo:iyhi,1:3)
c      
      integer:: icodeA, iplocal, ilocal
      integer:: ntout, timestep, idata
      integer:: i,j,k,l
      integer:: FileNumber
      character*50 DumpFile
c

c
      icodeA = int(nx*0.8d0) !nx/4*3 +10  !!!!!!!!!!!!!!!!! define here
      iplocal = int(icodeA/nxlocal)+1
      ilocal  = mod(icodeA,nxlocal)
c
c      ntout = mod(timestep,100)
      ntout = 1
      if(iprocx.eq.iplocal) then
      if(ntout.eq.1 .and. timestep .ge. 80000 ) then
c       data for output
        idata = ilocal
        uxdata(1,:,:,1)  = ux(idata,:,:,1)
        uxdata(1,:,:,2:3)=(ux(idata,:,:,2:3)+ux(idata-1,:,:,2:3))/2.d0
        uxdata(2:nghost+1,:,:,:)=ux(idata-nghost:idata-1,:,:,:)
        uwalld(1,:,1)  = bc_uwall(idata,:,1)
        uwalld(1,:,2:3)=
     &        (bc_uwall(idata,:,2:3)+bc_uwall(idata-1,:,2:3))/2.d0
        uwalld(2:nghost+1,:,:)=bc_uwall(idata-nghost:idata-1,:,:)

c       open file and write
        FileNumber = 1*10000 + iprocy*100 + iprocz
        write(DumpFile, 1102) float(FileNumber) / 1000000,
     &                        float(timestep) / 1000000
1102    format('outcodeA1', f7.6,f7.6)
c
        open(84,file=DumpFile)
        write(84,*) 'writing file for code B', timestep, ttot
c       ux
        do k=izlo, izhi, 1
          do j=iylo, iyhi, 1
            do i=1, nghost+1, 1
            write(84,*)
     &                 uxdata(i,j,k,1),uxdata(i,j,k,2),uxdata(i,j,k,3)
            enddo
          enddo
        enddo     
c       uwall
        do j= iylo, iyhi, 1
          do i=1, nghost+1
            write(84,*)uwalld(i,j,1),uwalld(i,j,2),uwalld(i,j,3)
          enddo
        enddo
c       utau
        do j=1, nylocal, 1
          write(84,*) bc_utau(ilocal,j)
        enddo         
c      
        close(84)
      endif
      endif
c 
      return
      end
#endif
c
c----------------------------------------------------------------------
c     Dumpout flow for restart case
c-------------------------------------------------------------------------
      subroutine DumpOutFluidB(ux,timestep,dflag)
      use mesh_parms
      use mesh_common
      use properties
      use time
#ifdef LES
      use virtualwallBC
#endif
c
      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      character*50 infile, outfile
      integer:: nbytes, offset,itmp
      integer:: FileNumber
      character*50 DumpFile
      integer:: i,j,k,l
      integer:: dflag, timestep
c       write(*,*) timestep, dflag, 'd'     
      if(dflag.eq.1) then
        FileNumber = iprocx*10000 + iprocy*100 + iprocz
        write(DumpFile, 1001) float(FileNumber)  / 1000000
1001    format('dump1', f7.6)
      else
        FileNumber = iprocx*10000 + iprocy*100 + iprocz
        write(DumpFile, 1002) float(FileNumber)  / 1000000
1002    format('dump2', f7.6)
      endif        

c        open(72, file=DumpFile, form='unformatted')
      open(72, file=DumpFile)
c      write(6,*) 'writing to file',iproc_idx,DumpFile
      write(72,*) timestep, ttot, Reynolds, velRef 
      do i = IXLO, IXHI
        do j = IYLO, IYHI
          do k = IZLO, IZHI
            do l = 1, 3
              write(72,*) ux(i, j, k, l)
            enddo
          enddo
        enddo
      enddo
#ifdef LES
      do i = IXLO, IXHI
        do j = IYLO, IYHI
          do l = 1, 3
            write(72,*) bc_uwall(i,j,l)
          enddo
        enddo
      enddo
      do i = 1, nxlocal+1
        do j = 1, nylocal
          write(72,*) bc_utau(i,j)
        enddo
      enddo
#endif

c        write(72)ux
      close(72)
      return
      end
c
c-------------------------------------------------------------------------
c     Dumpin flow for restart case    
c-------------------------------------------------------------------------
      subroutine DumpInFluidB(ux,dflag)
      use mesh_parms
      use mesh_common
      use properties
      use time
#ifdef LES
      use virtualwallBC
#endif
c
      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      character*50 infile, outfile
      integer:: nbytes, offset,itmp
      integer:: FileNumber
      character*50 DumpFile
      integer:: i,j,k,l
      integer:: dflag, timestep
      
      FileNumber = iprocx*10000 + iprocy*100 + iprocz
      write(DumpFile, 1002) float(FileNumber) / 1000000
1002   format('dump2', f7.6)

c        open(73, file=DumpFile, form='unformatted')
       open(73, file=DumpFile)
       write(6,*) 'reading from file',iproc_idx,DumpFile
       read(73,*) timestep, ttot, Reynolds, velRef 
       do i = IXLO, IXHI
         do j = IYLO, IYHI
           do k = IZLO, IZHI
             do l = 1, 3
               read(73,*) ux(i, j, k, l)
             enddo
           enddo
         enddo
       enddo
#ifdef LES
       do i = IXLO, IXHI
         do j = IYLO, IYHI
           do l = 1, 3
             read(73,*) bc_uwall(i,j,l)
c             write(*,*) bc_uwall(i,j,l)
           enddo
         enddo
       enddo
       do i = 1, nxlocal+1
         do j = 1, nylocal
           read(73,*) bc_utau(i,j)
           bc_dudz(i,j) = bc_utau(i,j)*bc_utau(i,j)*Reynolds
         enddo
       enddo
#endif       

       close(73)
c
       return
       end
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c     Dumpout flow for restart case
c-------------------------------------------------------------------------
      subroutine DumpOutFluidA(ux,timestep,dflag)
      use mesh_parms
      use mesh_common
      use properties
      use time
      use stats
      use recycle
#ifdef LES
      use virtualwallBC
#endif
c
      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      character*50 infile, outfile
      integer:: nbytes, offset,itmp
      integer:: FileNumber
      character*50 DumpFile
      integer:: i,j,k,l
      integer:: dflag, timestep
c       write(*,*) timestep, dflag, 'd'     
      if(dflag.eq.1) then
        FileNumber = iprocx*10000 + iprocy*100 + iprocz
        write(DumpFile, 1001) float(FileNumber)  / 1000000
1001    format('dump1', f7.6)
      else
        FileNumber = iprocx*10000 + iprocy*100 + iprocz
        write(DumpFile, 1002) float(FileNumber)  / 1000000
1002    format('dump2', f7.6)
      endif        

c        open(72, file=DumpFile, form='unformatted')
      open(72, file=DumpFile)
c      write(6,*) 'writing to file',iproc_idx,DumpFile
      write(72,*) timestep, ttot, Reynolds, velRef 
      do i = IXLO, IXHI
        do j = IYLO, IYHI
          do k = IZLO, IZHI
            do l = 1, 3
              write(72,*) ux(i, j, k, l)
            enddo
          enddo
        enddo
      enddo
#ifdef LES
      do i = IXLO, IXHI
        do j = IYLO, IYHI
          do l = 1, 3
            write(72,*) bc_uwall(i,j,l)
          enddo
        enddo
      enddo
      do i = 1, nxlocal+1
        do j = 1, nylocal
          write(72,*) bc_utau(i,j)
        enddo
      enddo
#endif
      do i= 1, nxlocal+1, 1
        do k= izlo, izhi, 1
          do l=1, 3, 1
            write(72,*) st_umean(i,k,l)
          enddo
        enddo
      enddo
      do i= 1, nx, 1
        write(72,*) g_del99(i), g_delta(i), g_theta(i), g_vwall(i),
     &              g_utau(i)
      enddo

c      write(*,*)'ddddddfffffffffffffff',g_ire, g_delta(g_ire)
c        write(72)ux
      close(72)
      return
      end
c
c-------------------------------------------------------------------------
c     Dumpin flow for restart case    
c-------------------------------------------------------------------------
      subroutine DumpInFluidA(ux,dflag)
      use mesh_parms
      use mesh_common
      use properties
      use time
      use stats
      use recycle
#ifdef LES
      use virtualwallBC
#endif
c
      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      character*50 infile, outfile
      integer:: nbytes, offset,itmp
      integer:: FileNumber
      character*50 DumpFile
      integer:: i,j,k,l
      integer:: dflag, timestep
      
      FileNumber = iprocx*10000 + iprocy*100 + iprocz
      write(DumpFile, 1002) float(FileNumber) / 1000000
1002   format('dump1', f7.6)

c        open(73, file=DumpFile, form='unformatted')
       open(73, file=DumpFile)
       write(6,*) 'reading from file',iproc_idx,DumpFile
       read(73,*) timestep, ttot, Reynolds, velRef 
       do i = IXLO, IXHI
         do j = IYLO, IYHI
           do k = IZLO, IZHI
             do l = 1, 3
               read(73,*) ux(i, j, k, l)
             enddo
           enddo
         enddo
       enddo
#ifdef LES
       do i = IXLO, IXHI
         do j = IYLO, IYHI
           do l = 1, 3
             read(73,*) bc_uwall(i,j,l)
c             write(*,*) bc_uwall(i,j,l)
           enddo
         enddo
       enddo
       do i = 1, nxlocal+1
         do j = 1, nylocal
           read(73,*) bc_utau(i,j)
           bc_dudz(i,j) = bc_utau(i,j)*bc_utau(i,j)*Reynolds
         enddo
       enddo
#endif       
      do i= 1, nxlocal+1, 1
        do k= izlo, izhi, 1
          do l=1, 3, 1
            read(73,*) st_umean(i,k,l)
          enddo
        enddo
      enddo
      do i= 1, nx, 1
        read(73,*) g_del99(i), g_delta(i), g_theta(i), g_vwall(i),
     &             g_utau(i)
      enddo
c      write(*,*)'ddddddfffffffffffffff',g_ire, g_delta(g_ire)
      close(73)
c-------------------------------------

c
       return
       end
c
