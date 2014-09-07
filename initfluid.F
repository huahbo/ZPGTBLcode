c---------------------------------------------------------------------------
c Subroutine to set up the initial conditions
c---------------------------------------------------------------------------
c
      subroutine InitFluid(ux)
      use mesh
      use mesh_common
      use properties
      implicit none
        double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
       integer:: i,j,k
c        double precision:: x,y,z
        double precision:: u0
c
       if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
       write(8,*) 'Mesh Size', nx, ny, nz
       write(8,*) 'Local Mesh Size', nxlsize, nylsize, nzlsize
       write(8,*) 'Initial conditions'
        endif
c
        Reynolds=1.0D5
        velRef = 1.d0
c
        u0=1.D0
      do 19 k=IZLO, IZHI,1
      do 20 j=IYLO, IYHI,1
      do 21 i=IXLO, IXHI,1
         ux(i,j,k,1)= u0 !xc(i) !u0
         ux(i,j,k,2)= 0.d0 !yc(j) !0.D0
         ux(i,j,k,3)= 0.d0 !zc(k) !0.D0
 21   continue
 20   continue
 19   continue
c
c----------------------------------
c        call blasiusinflow(ux)    
#ifdef CODEA
      call spalding(ux)
#endif
#ifdef CODEB
       call readinflowB(1)
       call initcodeB(ux,1)
#endif
c----------------------------------
c
c      write(6,*) 'EXCHANGE ...',iprocx,iprocy,iprocz
        call exchange(ux,3)
c      write(6,*) 'EXCHANGE ... done',iprocx,iprocy,iprocz
c        call WriteTmpSiloF(ux(ixlo,iylo,izlo,1),1,1,0,0)
      return
      end
c
c-----------------------------------------------------------------------






