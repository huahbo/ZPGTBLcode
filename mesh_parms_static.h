        module mesh_parms
        use mesh_uparms
        save
C          number of processes
        integer,parameter:: NPROCS=XPROCS*YPROCS*ZPROCS

C          sizes for local sections without ghost boundaries
        integer,parameter:: NXlsize=NX/XPROCS
        integer,parameter:: NYlsize=NY/YPROCS
        integer,parameter:: NZlsize=NZ/ZPROCS
c
C          bounds for local sections with ghost boundaries
        integer,parameter:: IXLO =1-NGHOST
        integer,parameter:: IXHI=nxlsize+NGHOST
        integer,parameter:: IYLO =1-NGHOST
        integer,parameter:: IYHI=nylsize+NGHOST
c          integer,parameter:: IZLO =1
c          integer,parameter:: IZHI=nzlsize
        integer,parameter:: IZLO =1-NGHOST
        integer,parameter:: IZHI=nzlsize+NGHOST
c          integer,parameter:: INLO=IXLO
c          integer,parameter::  INHI=IXHI
c
        integer,parameter:: MSG_XCH_XLOW_TAG=1,MSG_XCH_XHI_TAG=2
        integer,parameter:: MSG_XCH_YLOW_TAG=3,MSG_XCH_YHI_TAG=4
        integer,parameter:: MSG_XCH_ZLOW_TAG=5,MSG_XCH_ZHI_TAG=6
        integer,parameter:: MSG_MAX_TAG=10

        integer:: NXlocal, NYlocal, NZlocal
c        integer:: ilo,ihi,jlo,jhi,klo,khi
c        integer:: hasxlo, hasxhi, hasylo,hasyhi
c
        end module mesh_parms
