C
        Module mesh_uparms
        save
        integer,parameter:: XPROCS=4
        integer,parameter:: YPROCS=1
        integer,parameter:: ZPROCS=2
C
        integer,parameter:: nx=64
        integer,parameter:: ny=64
        integer,parameter:: nz=64
C
C
        integer,parameter:: nghost=2
        integer,parameter:: nvar=6
        integer,parameter:: ndim=3
        end module mesh_uparms
