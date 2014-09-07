        module mesh_parms
	use mesh_uparms
	save
C	number of processes
	integer:: NPROCS

C	sizes for local sections without ghost boundaries
	integer:: NXlsize
	integer:: NYlsize
        integer:: NZlsize
c
C	bounds for local sections with ghost boundaries
	integer:: IXLO
	integer:: IXHI
	integer:: IYLO
	integer:: IYHI
	integer:: IZLO
	integer:: IZHI
c
        integer,parameter:: MSG_XCH_XLOW_TAG=1,MSG_XCH_XHI_TAG=2
        integer,parameter:: MSG_XCH_YLOW_TAG=3,MSG_XCH_YHI_TAG=4
        integer,parameter:: MSG_XCH_ZLOW_TAG=5,MSG_XCH_ZHI_TAG=6
        integer,parameter:: MSG_MAX_TAG=10

	integer:: NXlocal, NYlocal, NZlocal
c
        end module mesh_parms
