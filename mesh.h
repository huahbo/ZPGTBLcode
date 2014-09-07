        module mesh
        use mesh_parms
        save
c
        double precision:: dx,dy,dz
        double precision:: dxi,dyi,dzi
        double precision:: dxi3,dyi3,dzi3
c
        double precision, allocatable,
     &      dimension(:):: xc,yc,zc
c
          integer:: hasxlo,hasylo,hasxhi,hasyhi
        double precision:: xbound(2), ybound(2),zbound(2)
        end module mesh
