c
c----------------------------------------------------------------------
c
c      variables used in virtual wall model
       module virtualwallBC
       save
        double precision:: bc_h0dz, bc_hdz, bc_hnuplus
        double precision,allocatable:: bc_utau(:,:)
        double precision,allocatable:: bc_dudz(:,:)
        double precision,allocatable:: bc_Txx(:,:)
        double precision,allocatable:: bc_Txy(:,:)
        double precision,allocatable:: bc_Txz(:,:)
        double precision,allocatable:: bc_Kske(:,:) 
        double precision,allocatable:: bc_Kappa(:,:)
        double precision,allocatable:: bc_dudzh(:,:)
        double precision,allocatable:: bc_uw(:,:)
        double precision,allocatable:: bc_convec(:,:)
        double precision,allocatable:: bc_dpdx(:,:)
        double precision,allocatable:: bc_rhs(:,:)
        double precision,allocatable:: bc_rhs_old(:,:)
        double precision,allocatable:: bc_uwall(:,:,:)
       end module virtualwallBC
c
c
c      allocate the variables in virtual wall model
       subroutine allocate_LES
       use mesh
       use mesh_common
       use properties
       use virtualwallBC
c
        allocate(  bc_utau(1:nxlocal+1,1:nylocal))
        allocate(  bc_dudz(1:nxlocal+1,1:nylocal))
        allocate(  bc_Txx(-1:nxlocal+2,1:nylocal))
        allocate(  bc_Txy( 1:nxlocal+1,0:nylocal+2))
        allocate(  bc_Txz( 1:nxlocal+1,1:nylocal))
        allocate( bc_Kske( 1:nxlocal+1,1:nylocal))
        allocate(bc_Kappa( 1:nxlocal+1,1:nylocal))
        allocate(bc_dudzh( 1:nxlocal+1,1:nylocal))
        allocate(   bc_uw( 1:nxlocal+1,1:nylocal))
        allocate(bc_convec(1:nxlocal+1,1:nylocal))
        allocate( bc_dpdx( 1:nxlocal+1,1:nylocal))
        allocate(  bc_rhs( 1:nxlocal+1,1:nylocal))
        allocate(bc_rhs_old(1:nxlocal+1,1:nylocal))
        allocate(bc_uwall( ixlo:ixhi,iylo:iyhi,3))
        bc_utau = 1.d7
        bc_dudz = 1.d7
        bc_Txx  = 1.d7
        bc_Txy  = 1.d7
        bc_Txz  = 1.d7
        bc_Kske = 1.d7
        bc_Kappa= 1.d7
        bc_dudzh= 1.d7
        bc_uw   = 1.d7        !--------------
        bc_convec= 1.d7
        bc_dpdx = 0.d0
        bc_rhs = 1.d7
        bc_rhs_old = 1.d7
c
        bc_uwall = 1.d7
c-------------------------
        bc_h0dz =  0.18d0
        bc_hdz = bc_h0dz + 0.5d0
        bc_hnuplus = 10.23d0 
c
        return
        end      
c
c
       module GClaw
       save
        double precision,allocatable::  GC_ci(:,:)
        double precision,allocatable::  GC_co(:,:)
        double precision,allocatable:: GC_gam(:,:)
        double precision,allocatable:: GC_pre(:,:)
        double precision:: GCap, GCgaminfty, GCa, GCalpha, GCcoinfty
        double precision:: GCciinfty, GCc1, GCc2, GChnu
        double precision:: GCgam, GCco, GCci, GCzp, GCRe 
       end module GClaw

       subroutine init_GClaw
       use mesh
       use mesh_common
       use GClaw
c
       implicit none
       allocate(  GC_ci(nxlocal+1,nylocal))
       allocate(  GC_co(nxlocal+1,nylocal))
       allocate( GC_gam(nxlocal+1,nylocal))
       allocate( GC_pre(nxlocal+1,nylocal))
c
       GC_ci  = 0.d0
       GC_co  = 0.d0
       GC_gam = 0.d0
       GC_pre = 0.d0
c
       GCap       = -16.d0
       GChnu      =  30.d0
c       GCgaminfty = 0.0362d0
       GCgaminfty = 0.045d0

       GCa        = 2.90d0
       GCalpha    = 0.46d0
       GCcoinfty   = 0.897d0
c       GCciinfty   = 55d0
       GCciinfty  = 48.d0
       GCc1       = 0.283d0
       GCc2       = -0.00598d0
c
       return
       end    
