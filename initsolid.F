c---------------------------------------------------------------------------
c Subroutine to set up the initial positions and velocities
c for the solid.
c---------------------------------------------------------------------------
      subroutine InitSolid(dts)
      use mesh
      use mesh_common
      use SignedDistanceModule
      use SignedDistanceModuleCoords
c
      implicit none
      double precision:: dts
      integer::iv
      double precision:: omega
c      return
      call read_mesh
      write(6,*) 'SOLID EXTENTS:'
      write(6,*) ' X:', minval(xs(1:nvert)),maxval(xs(1:nvert))
      write(6,*) ' Y:', minval(ys(1:nvert)),maxval(ys(1:nvert))
      write(6,*) ' Z:', minval(zs(1:nvert)),maxval(zs(1:nvert))
      call setup_sdf
c
c     Zero velocities
      us=0.D0
      vs=0.D0
      ws=0.D0
      omega=-10.D0
c      return
      do iv=1,nvert,1
c         vs(iv)=0.1D0
         vs(iv)=-zs(iv)*omega
         ws(iv)=ys(iv)*omega
      enddo
      write(6,*) 'INIT VEL',maxval(vs),maxval(ws)
        return
        end
c
c-----------------------------------------------------------------------


